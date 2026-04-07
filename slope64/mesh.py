"""Mesh generation for Slope64 trapezoidal slope geometry.

Generates node coordinates and 8-node quadrilateral element connectivity
for a two-zone mesh: embankment (trapezoid) + foundation (rectangle).

Coordinate system: x horizontal (left to right), y vertical (up positive).
Origin at bottom-left corner of the mesh.

Node numbering follows the Smith & Griffiths convention for 8-node quads:
    7---6---5
    |       |
    8       4
    |       |
    1---2---3

The mesh is built on a super-grid, but only nodes referenced by actual
elements are retained.
"""

from __future__ import annotations

import numpy as np

from slope64.parser import SlopeInput


def generate_mesh(inp: SlopeInput) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate the FE mesh for a Slope64 problem.

    Returns:
        coords: (n_nodes, 2) array of (x, y) coordinates
        connect: (n_elements, 8) array of node indices (0-based) per element
        etype_map: (n_elements,) array of material group index (0-based)
    """
    nx1, nx2 = inp.nx1, inp.nx2
    ny1, ny2 = inp.ny1, inp.ny2
    w1, s1, w2 = inp.w1, inp.s1, inp.w2
    h1, h2 = inp.h1, inp.h2

    nx_total = nx1 + nx2
    ny_total = ny1 + ny2

    # Super-grid dimensions
    n_sg_rows = 2 * ny_total + 1
    n_sg_cols = 2 * nx_total + 1

    # x-coordinates for super-grid columns
    x_emb = np.linspace(0.0, w1 + s1, 2 * nx1 + 1)
    if nx2 > 0 and w2 > 0:
        x_found = np.linspace(w1 + s1, w1 + s1 + w2, 2 * nx2 + 1)[1:]
        x_all = np.concatenate([x_emb, x_found])
    else:
        x_all = x_emb

    # y-coordinates for super-grid rows
    if ny2 > 0:
        y_found = np.linspace(0.0, h2, 2 * ny2 + 1)
        y_emb = np.linspace(h2, h2 + h1, 2 * ny1 + 1)[1:]
        y_all = np.concatenate([y_found, y_emb])
    else:
        y_all = np.linspace(0.0, h1, 2 * ny1 + 1)

    # Build super-grid coordinates with slope deformation
    sg_coords = np.zeros((n_sg_rows, n_sg_cols, 2))
    for j in range(n_sg_rows):
        y = y_all[j]
        for k in range(n_sg_cols):
            sg_coords[j, k, 1] = y
            if y > h2 + 1e-12 and k < 2 * nx1 + 1:
                # Embankment zone: compress x to fit within slope
                frac = (y - h2) / h1
                x_right = w1 + s1 * (1.0 - frac)
                sg_coords[j, k, 0] = x_all[k] * x_right / (w1 + s1)
            else:
                sg_coords[j, k, 0] = x_all[min(k, len(x_all) - 1)]

    # Build element connectivity on super-grid and collect used nodes
    elements_sg = []  # list of [8 super-grid (row, col) tuples]
    etype_list = []

    for ei in range(ny_total):
        # Mesh rows go bottom-to-top: ei=0 is bottom (foundation), ei=ny_total-1 is top
        nx_this_row = nx_total if ei < ny2 else nx1
        # etype rows in .dat go top-to-bottom: row 0 is top embankment
        etype_row_idx = ny_total - 1 - ei

        for ej in range(nx_this_row):
            r = 2 * ei
            c = 2 * ej
            sg_nodes = [
                (r, c),         # 1
                (r, c + 1),     # 2
                (r, c + 2),     # 3
                (r + 1, c + 2), # 4
                (r + 2, c + 2), # 5
                (r + 2, c + 1), # 6
                (r + 2, c),     # 7
                (r + 1, c),     # 8
            ]
            elements_sg.append(sg_nodes)

            if inp.np_types == 1:
                etype_list.append(0)
            else:
                etype_list.append(inp.etype[etype_row_idx][ej] - 1)

    # Collect unique super-grid positions used by elements
    used_sg = set()
    for sg_nodes in elements_sg:
        for rc in sg_nodes:
            used_sg.add(rc)

    # Sort for consistent numbering (row-major, bottom to top, left to right)
    used_sg_sorted = sorted(used_sg)
    sg_to_node = {rc: idx for idx, rc in enumerate(used_sg_sorted)}

    # Build coordinate and connectivity arrays
    n_nodes = len(used_sg_sorted)
    coords = np.zeros((n_nodes, 2))
    for idx, (r, c) in enumerate(used_sg_sorted):
        coords[idx] = sg_coords[r, c]

    connect = np.zeros((len(elements_sg), 8), dtype=np.int32)
    for el, sg_nodes in enumerate(elements_sg):
        for k, rc in enumerate(sg_nodes):
            connect[el, k] = sg_to_node[rc]

    etype_map = np.array(etype_list, dtype=np.int32)

    return coords, connect, etype_map


def get_boundary_conditions(inp: SlopeInput, coords: np.ndarray) -> np.ndarray:
    """Determine fixed DOFs for boundary conditions.

    Returns:
        fixed_dofs: array of DOF indices that are fixed (0-based).
            DOF numbering: node i has DOFs 2*i (x) and 2*i+1 (y).
    """
    tol = 1e-6
    fixed = []
    right_edge = inp.w1 + inp.s1 + inp.w2

    for i in range(coords.shape[0]):
        x, y = coords[i]
        # Bottom boundary: fully fixed
        if abs(y) < tol:
            fixed.append(2 * i)
            fixed.append(2 * i + 1)
        # Left boundary: roller (x fixed)
        elif abs(x) < tol:
            fixed.append(2 * i)
        # Right boundary of foundation: roller (x fixed)
        elif abs(x - right_edge) < tol:
            fixed.append(2 * i)

    return np.array(sorted(set(fixed)), dtype=np.int32)
