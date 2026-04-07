"""Global FEM assembly and visco-plastic iterative solver.

Uses the visco-plastic algorithm (Smith & Griffiths 1998):
1. Factor K once (LU decomposition).
2. Elastic solve with gravity loads.
3. Check Mohr-Coulomb at all Gauss points.
4. For yielded GPs, compute stress correction body loads.
5. Accumulate body loads, re-solve.
6. Repeat until convergence or iteration limit.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import splu

from slope64.parser import SlopeInput
from slope64.elements import (
    GAUSS_POINTS,
    GAUSS_WEIGHTS,
    shape_functions,
    shape_derivatives,
    jacobian,
    b_matrix,
    d_matrix_plane_strain,
    element_stiffness,
)
from slope64.material import factor_properties


def assemble_stiffness(
    coords: np.ndarray,
    connect: np.ndarray,
    etype_map: np.ndarray,
    properties: list[list[float]],
) -> csr_matrix:
    """Assemble global stiffness matrix using COO format."""
    n_dofs = 2 * coords.shape[0]
    n_elem = connect.shape[0]

    rows = np.zeros(n_elem * 256, dtype=np.int32)
    cols = np.zeros(n_elem * 256, dtype=np.int32)
    vals = np.zeros(n_elem * 256)

    for el in range(n_elem):
        mat = properties[etype_map[el]]
        E, v = mat[4], mat[5]
        elem_nodes = connect[el]
        elem_coords = coords[elem_nodes]
        ke = element_stiffness(elem_coords, E, v)

        dofs = np.empty(16, dtype=np.int32)
        dofs[0::2] = 2 * elem_nodes
        dofs[1::2] = 2 * elem_nodes + 1

        base = el * 256
        idx = 0
        for a in range(16):
            for b in range(16):
                rows[base + idx] = dofs[a]
                cols[base + idx] = dofs[b]
                vals[base + idx] = ke[a, b]
                idx += 1

    return coo_matrix((vals, (rows, cols)), shape=(n_dofs, n_dofs)).tocsr()


def assemble_gravity_loads(
    coords: np.ndarray,
    connect: np.ndarray,
    etype_map: np.ndarray,
    properties: list[list[float]],
    inp: SlopeInput,
    k_h: float = 0.0,
) -> np.ndarray:
    """Assemble gravity + seismic body-force load vector.

    Below the phreatic surface, uses buoyant weight (gamma - gam_w) for the
    vertical body force plus the horizontal seepage gradient force.  Above the
    phreatic surface (or with no water), total gamma is used.

    D*B*d then gives effective stresses directly; the yield check needs no
    pore-pressure correction.
    """
    from slope64.water import interpolate_surface, interpolate_surface_gradient

    n_dofs = 2 * coords.shape[0]
    F_global = np.zeros(n_dofs)
    has_water = inp.gam_w > 0 and inp.nosurf > 0
    crest_y = inp.h2 + inp.h1

    for el in range(connect.shape[0]):
        mat = properties[etype_map[el]]
        gamma = mat[3]
        elem_nodes = connect[el]
        elem_coords = coords[elem_nodes]

        fe = np.zeros(16)

        for gi in range(4):
            xi, eta = GAUSS_POINTS[gi]
            w = GAUSS_WEIGHTS[gi]
            N = shape_functions(xi, eta)
            dNdxi = shape_derivatives(xi, eta)
            _, det_j = jacobian(dNdxi, elem_coords)

            gamma_eff = gamma
            seepage_x = 0.0

            if has_water:
                x_gp = float(N @ elem_coords[:, 0])
                y_gp = float(N @ elem_coords[:, 1])
                y_surf_rel = interpolate_surface(x_gp, inp.surf)
                y_water = crest_y + y_surf_rel
                if y_gp < y_water - 1e-10:
                    gamma_eff = gamma - inp.gam_w
                    dy_water_dx = interpolate_surface_gradient(x_gp, inp.surf)
                    seepage_x = -inp.gam_w * dy_water_dx

            # Gravity (buoyant below water table, total above)
            for k in range(8):
                fe[2 * k + 1] += -gamma_eff * N[k] * w * det_j

            # Horizontal seepage force
            if abs(seepage_x) > 1e-12:
                for k in range(8):
                    fe[2 * k] += seepage_x * N[k] * w * det_j

            # Seismic horizontal inertia (always total gamma)
            if abs(k_h) > 1e-12:
                for k in range(8):
                    fe[2 * k] += gamma * k_h * N[k] * w * det_j

        for k in range(8):
            F_global[2 * elem_nodes[k]] += fe[2 * k]
            F_global[2 * elem_nodes[k] + 1] += fe[2 * k + 1]

    return F_global


def _element_horizontal_load(elem_coords: np.ndarray, force_per_vol: float) -> np.ndarray:
    fe = np.zeros(16)
    for gp in range(4):
        xi, eta = GAUSS_POINTS[gp]
        w = GAUSS_WEIGHTS[gp]
        N = shape_functions(xi, eta)
        dNdxi = shape_derivatives(xi, eta)
        _, det_j = jacobian(dNdxi, elem_coords)
        for k in range(8):
            fe[2 * k] += force_per_vol * N[k] * w * det_j
    return fe


def _precompute_element_data(coords, connect):
    """Pre-compute B matrices, det_j, weights, GP coords for all elements.

    Returns dict: B (n_elem,4,3,16), wdj (n_elem,4), xy (n_elem,4,2).
    """
    n_elem = connect.shape[0]
    B_all = np.zeros((n_elem, 4, 3, 16))
    wdj_all = np.zeros((n_elem, 4))
    xy_gp = np.zeros((n_elem, 4, 2))

    for el in range(n_elem):
        elem_nodes = connect[el]
        elem_coords = coords[elem_nodes]
        for gp in range(4):
            xi, eta = GAUSS_POINTS[gp]
            w = GAUSS_WEIGHTS[gp]
            N = shape_functions(xi, eta)
            dNdxi = shape_derivatives(xi, eta)
            jac_mat, det_j = jacobian(dNdxi, elem_coords)
            B_all[el, gp] = b_matrix(dNdxi, jac_mat, det_j)
            wdj_all[el, gp] = w * det_j
            xy_gp[el, gp, 0] = float(N @ elem_coords[:, 0])
            xy_gp[el, gp, 1] = float(N @ elem_coords[:, 1])

    return {"B": B_all, "wdj": wdj_all, "xy": xy_gp}


def solve_step(
    coords: np.ndarray,
    connect: np.ndarray,
    etype_map: np.ndarray,
    inp: SlopeInput,
    fos: float,
    limit: int,
) -> tuple[np.ndarray, float, int, bool]:
    """Perform one SRF step with visco-plastic iteration.

    Uses the Perzyna visco-plastic algorithm (Smith & Griffiths, Program 6.4):
    1. Elastic solve with gravity loads.
    2. At each GP: elastic strain = total strain - accumulated VP strain.
    3. Check Mohr-Coulomb yield on effective stress (D*B*d).
    4. For yielded GPs: accumulate VP strain, compute body loads.
    5. Solve increment: K * du = bdylds, u += du.
    6. Repeat until convergence or iteration limit.
    """
    from slope64.mesh import get_boundary_conditions

    n_dofs = 2 * coords.shape[0]
    n_elem = connect.shape[0]

    # Factor material properties (only phi and c)
    factored_props = []
    for mat in inp.properties:
        phi_deg, c, psi_deg, gamma, E, v = mat
        phi_f, c_f = factor_properties(phi_deg, c, fos)
        factored_props.append([phi_f, c_f, psi_deg, gamma, E, v])

    # Assemble K (E and v don't depend on FOS)
    K = assemble_stiffness(coords, connect, etype_map, inp.properties)
    F_ext = assemble_gravity_loads(coords, connect, etype_map, inp.properties, inp, inp.k_h)

    # Pre-compute element geometry data
    ed = _precompute_element_data(coords, connect)
    B_all = ed["B"]      # (n_elem, 4, 3, 16)
    wdj_all = ed["wdj"]  # (n_elem, 4)
    xy_gp = ed["xy"]     # (n_elem, 4, 2)

    # Element DOF arrays
    elem_dofs = np.zeros((n_elem, 16), dtype=np.int32)
    for el in range(n_elem):
        elem_nodes = connect[el]
        elem_dofs[el, 0::2] = 2 * elem_nodes
        elem_dofs[el, 1::2] = 2 * elem_nodes + 1

    # Boundary conditions
    fixed_dofs = get_boundary_conditions(inp, coords)
    free_dofs = np.setdiff1d(np.arange(n_dofs), fixed_dofs)

    # LU factorize K once
    K_ff = K[np.ix_(free_dofs, free_dofs)]
    lu = splu(K_ff.tocsc())

    # Initial elastic solve
    displ = np.zeros(n_dofs)
    displ[free_dofs] = lu.solve(F_ext[free_dofs])

    # Pre-compute D per element
    D_per_elem = np.zeros((n_elem, 3, 3))
    for idx, mat in enumerate(inp.properties):
        D_per_elem[etype_map == idx] = d_matrix_plane_strain(mat[4], mat[5])

    # Per-element factored material parameters
    sin_phi_arr = np.zeros(n_elem)
    cos_phi_arr = np.zeros(n_elem)
    sin_psi_arr = np.zeros(n_elem)
    c_f_arr = np.zeros(n_elem)
    for idx, mat in enumerate(factored_props):
        mask = etype_map == idx
        phi_rad = np.radians(mat[0])
        psi_rad = np.radians(mat[2])
        sin_phi_arr[mask] = np.sin(phi_rad)
        cos_phi_arr[mask] = np.cos(phi_rad)
        sin_psi_arr[mask] = np.sin(psi_rad)
        c_f_arr[mask] = mat[1]

    # Cormeau (1975) critical time step — uses factored phi per S&G p64
    dt = float("inf")
    for mat_f in factored_props:
        phi_f_deg = mat_f[0]
        E = mat_f[4]
        v = mat_f[5]
        snph = np.sin(np.radians(phi_f_deg))
        denom_dt = E * (1.0 - 2.0 * v + snph ** 2)
        if denom_dt < 1e-30:
            denom_dt = E * (1.0 - 2.0 * v)
        dt = min(dt, 4.0 * (1.0 + v) * (1.0 - 2.0 * v) / denom_dt)

    # Accumulated viscoplastic strain at each Gauss point
    evpt = np.zeros((n_elem, 4, 3))

    ext_norm = np.linalg.norm(F_ext[free_dofs])

    # Visco-plastic iteration (vectorized)
    for iteration in range(1, limit + 1):
        bdylds = np.zeros(n_dofs)

        # Gather element displacements: (n_elem, 16)
        u_elem = displ[elem_dofs]

        total_n_yielded = 0

        for gp in range(4):
            # Total strain at all elements for this GP: (n_elem, 3)
            strain = np.einsum("eij,ej->ei", B_all[:, gp], u_elem)

            # Elastic strain = total - accumulated VP
            elastic_strain = strain - evpt[:, gp]

            # Effective stress = D @ elastic_strain
            # (pore pressure handled in load assembly, so D*B*d = effective stress)
            stress = np.einsum("eij,ej->ei", D_per_elem, elastic_strain)

            sx = stress[:, 0]
            sy = stress[:, 1]
            txy = stress[:, 2]

            diff = sx - sy
            radius = np.sqrt(0.25 * diff ** 2 + txy ** 2)
            safe_r = np.maximum(radius, 1e-15)

            # Mohr-Coulomb yield on effective stress (pore pressure in loads → D*B*d = σ')
            F_val = radius + 0.5 * (sx + sy) * sin_phi_arr - c_f_arr * cos_phi_arr

            # Yielded mask
            yielded = F_val > 0.0
            n_y = np.count_nonzero(yielded)
            total_n_yielded += n_y

            if n_y == 0:
                continue

            # Plastic potential gradient dQ/dsigma (use psi)
            dQ = np.zeros((n_elem, 3))
            dQ[yielded, 0] = 0.25 * diff[yielded] / safe_r[yielded] + 0.5 * sin_psi_arr[yielded]
            dQ[yielded, 1] = -0.25 * diff[yielded] / safe_r[yielded] + 0.5 * sin_psi_arr[yielded]
            dQ[yielded, 2] = txy[yielded] / safe_r[yielded]

            # VP strain increment: evp = dt * F * dQ
            evp = np.zeros((n_elem, 3))
            evp[yielded] = dt * F_val[yielded, np.newaxis] * dQ[yielded]

            # Accumulate
            evpt[:, gp] += evp

            # Body loads: bdylds[dofs] += wdj * B^T @ (D @ evp)
            D_evp = np.einsum("eij,ej->ei", D_per_elem, evp)
            BtD = np.einsum("eij,ei->ej", B_all[:, gp], D_evp)
            contrib = wdj_all[:, gp, np.newaxis] * BtD  # (n_elem, 16)

            # Scatter to global: only for yielded elements
            y_idx = np.where(yielded)[0]
            for el in y_idx:
                bdylds[elem_dofs[el]] += contrib[el]

        # Converged if no Gauss points yielded
        if total_n_yielded == 0:
            return displ, _max_nodal_displ(displ), iteration, True

        # Solve for displacement increment
        du = np.zeros(n_dofs)
        du[free_dofs] = lu.solve(bdylds[free_dofs])
        displ += du

        # S&G checon convergence: max|du| / max|u| < tol
        du_max = np.max(np.abs(du))
        u_max = np.max(np.abs(displ))
        if iteration > 1 and u_max > 1e-30 and du_max / u_max < 0.0001:
            return displ, _max_nodal_displ(displ), iteration, True

        # Detect divergence early
        if u_max > 1e10:
            return displ, _max_nodal_displ(displ), iteration, False

    # Reached iteration limit — slope failure
    return displ, _max_nodal_displ(displ), limit, False


def _max_nodal_displ(displ: np.ndarray) -> float:
    ux = displ[0::2]
    uy = displ[1::2]
    return float(np.max(np.sqrt(ux ** 2 + uy ** 2)))
