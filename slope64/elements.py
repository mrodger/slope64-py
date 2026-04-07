"""8-node serendipity quadrilateral element formulations.

Shape functions, Jacobian, strain-displacement (B) matrix,
elastic stiffness (D) matrix, and 2x2 Gauss quadrature.

Local coordinate system: (xi, eta) in [-1, 1] x [-1, 1].

Node numbering (Smith & Griffiths convention):
    7---6---5
    |       |
    8       4
    |       |
    1---2---3

    Node 1: (-1, -1)  Node 2: (0, -1)  Node 3: (1, -1)
    Node 4: (1, 0)    Node 5: (1, 1)   Node 6: (0, 1)
    Node 7: (-1, 1)   Node 8: (-1, 0)
"""

from __future__ import annotations

import numpy as np

# 2x2 Gauss quadrature points and weights
_G = 1.0 / np.sqrt(3.0)
GAUSS_POINTS = np.array([
    [-_G, -_G],
    [_G, -_G],
    [_G, _G],
    [-_G, _G],
])
GAUSS_WEIGHTS = np.array([1.0, 1.0, 1.0, 1.0])

# Local coordinates of the 8 nodes
NODE_XI = np.array([-1.0, 0.0, 1.0, 1.0, 1.0, 0.0, -1.0, -1.0])
NODE_ETA = np.array([-1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 1.0, 0.0])


def shape_functions(xi: float, eta: float) -> np.ndarray:
    """Evaluate 8-node serendipity shape functions at (xi, eta).

    Returns:
        N: (8,) array of shape function values.
    """
    N = np.zeros(8)
    # Corner nodes (1, 3, 5, 7) — indices 0, 2, 4, 6
    for i in [0, 2, 4, 6]:
        xi_i, eta_i = NODE_XI[i], NODE_ETA[i]
        N[i] = 0.25 * (1.0 + xi_i * xi) * (1.0 + eta_i * eta) * (xi_i * xi + eta_i * eta - 1.0)
    # Midside nodes (2, 4, 6, 8) — indices 1, 3, 5, 7
    for i in [1, 5]:  # eta = ±1 side
        xi_i, eta_i = NODE_XI[i], NODE_ETA[i]
        N[i] = 0.5 * (1.0 - xi * xi) * (1.0 + eta_i * eta)
    for i in [3, 7]:  # xi = ±1 side
        xi_i, eta_i = NODE_XI[i], NODE_ETA[i]
        N[i] = 0.5 * (1.0 + xi_i * xi) * (1.0 - eta * eta)
    return N


def shape_derivatives(xi: float, eta: float) -> np.ndarray:
    """Evaluate derivatives of shape functions w.r.t. (xi, eta).

    Returns:
        dNdxi: (2, 8) array. dNdxi[0, :] = dN/dxi, dNdxi[1, :] = dN/deta.
    """
    dNdxi = np.zeros((2, 8))
    # Corner nodes
    for i in [0, 2, 4, 6]:
        xi_i, eta_i = NODE_XI[i], NODE_ETA[i]
        dNdxi[0, i] = 0.25 * xi_i * (1.0 + eta_i * eta) * (2.0 * xi_i * xi + eta_i * eta)
        dNdxi[1, i] = 0.25 * eta_i * (1.0 + xi_i * xi) * (xi_i * xi + 2.0 * eta_i * eta)
    # Midside: nodes on eta = const (indices 1, 5)
    for i in [1, 5]:
        eta_i = NODE_ETA[i]
        dNdxi[0, i] = -xi * (1.0 + eta_i * eta)
        dNdxi[1, i] = 0.5 * (1.0 - xi * xi) * eta_i
    # Midside: nodes on xi = const (indices 3, 7)
    for i in [3, 7]:
        xi_i = NODE_XI[i]
        dNdxi[0, i] = 0.5 * xi_i * (1.0 - eta * eta)
        dNdxi[1, i] = -(1.0 + xi_i * xi) * eta
    return dNdxi


def jacobian(dNdxi: np.ndarray, elem_coords: np.ndarray) -> tuple[np.ndarray, float]:
    """Compute the Jacobian matrix and its determinant.

    Args:
        dNdxi: (2, 8) shape function derivatives in local coords.
        elem_coords: (8, 2) nodal coordinates of the element.

    Returns:
        jac: (2, 2) Jacobian matrix.
        det_j: determinant of Jacobian.
    """
    jac = dNdxi @ elem_coords  # (2, 8) @ (8, 2) = (2, 2)
    det_j = jac[0, 0] * jac[1, 1] - jac[0, 1] * jac[1, 0]
    return jac, det_j


def b_matrix(dNdxi: np.ndarray, jac: np.ndarray, det_j: float) -> np.ndarray:
    """Compute the strain-displacement matrix B for plane strain.

    Args:
        dNdxi: (2, 8) shape function derivatives in local coords.
        jac: (2, 2) Jacobian matrix.
        det_j: Jacobian determinant.

    Returns:
        B: (3, 16) strain-displacement matrix.
            Rows: [eps_x, eps_y, gamma_xy]
            Cols: [u1, v1, u2, v2, ..., u8, v8]
    """
    # Inverse Jacobian
    inv_j = np.array([
        [jac[1, 1], -jac[0, 1]],
        [-jac[1, 0], jac[0, 0]],
    ]) / det_j

    # Shape function derivatives in global coords
    dNdx = inv_j @ dNdxi  # (2, 8): dN/dx and dN/dy

    B = np.zeros((3, 16))
    for k in range(8):
        B[0, 2 * k] = dNdx[0, k]       # eps_x = du/dx
        B[1, 2 * k + 1] = dNdx[1, k]   # eps_y = dv/dy
        B[2, 2 * k] = dNdx[1, k]        # gamma_xy = du/dy + dv/dx
        B[2, 2 * k + 1] = dNdx[0, k]
    return B


def d_matrix_plane_strain(E: float, v: float) -> np.ndarray:
    """Elastic stiffness matrix for plane strain.

    Returns:
        D: (3, 3) constitutive matrix relating [sigma_x, sigma_y, tau_xy]
           to [eps_x, eps_y, gamma_xy].
    """
    fac = E / ((1.0 + v) * (1.0 - 2.0 * v))
    D = fac * np.array([
        [1.0 - v, v, 0.0],
        [v, 1.0 - v, 0.0],
        [0.0, 0.0, 0.5 * (1.0 - 2.0 * v)],
    ])
    return D


def element_stiffness(elem_coords: np.ndarray, E: float, v: float) -> np.ndarray:
    """Compute the 16x16 element stiffness matrix using 2x2 Gauss quadrature.

    Args:
        elem_coords: (8, 2) nodal coordinates.
        E: Young's modulus.
        v: Poisson's ratio.

    Returns:
        ke: (16, 16) element stiffness matrix.
    """
    D = d_matrix_plane_strain(E, v)
    ke = np.zeros((16, 16))

    for gp in range(4):
        xi, eta = GAUSS_POINTS[gp]
        w = GAUSS_WEIGHTS[gp]

        dNdxi = shape_derivatives(xi, eta)
        jac_mat, det_j = jacobian(dNdxi, elem_coords)
        B = b_matrix(dNdxi, jac_mat, det_j)

        ke += w * det_j * (B.T @ D @ B)

    return ke


def element_gravity_load(elem_coords: np.ndarray, gamma: float) -> np.ndarray:
    """Compute the gravity load vector for an element.

    Gravity acts in the -y direction. Load = gamma * N * area.

    Returns:
        fe: (16,) element force vector.
    """
    fe = np.zeros(16)

    for gp in range(4):
        xi, eta = GAUSS_POINTS[gp]
        w = GAUSS_WEIGHTS[gp]

        N = shape_functions(xi, eta)
        dNdxi = shape_derivatives(xi, eta)
        _, det_j = jacobian(dNdxi, elem_coords)

        for k in range(8):
            # fy = -gamma (downward)
            fe[2 * k + 1] += -gamma * N[k] * w * det_j

    return fe


def stress_at_gauss_points(
    elem_coords: np.ndarray,
    elem_displ: np.ndarray,
    E: float,
    v: float,
) -> np.ndarray:
    """Compute stress at each Gauss point from nodal displacements.

    Returns:
        stresses: (4, 3) array. Each row is [sigma_x, sigma_y, tau_xy] at a GP.
    """
    D = d_matrix_plane_strain(E, v)
    stresses = np.zeros((4, 3))

    for gp in range(4):
        xi, eta = GAUSS_POINTS[gp]
        dNdxi = shape_derivatives(xi, eta)
        jac_mat, det_j = jacobian(dNdxi, elem_coords)
        B = b_matrix(dNdxi, jac_mat, det_j)
        strain = B @ elem_displ
        stresses[gp] = D @ strain

    return stresses
