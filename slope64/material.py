"""Mohr-Coulomb material model and visco-plastic stress redistribution.

Compression-negative convention (as in Griffiths & Lane 1999).

Mohr-Coulomb failure function (derived from Mohr circle geometry):
    F = R + p*sin(phi) - c*cos(phi)

where R = (sigma1 - sigma3)/2 (Mohr circle radius, always positive)
      p = (sigma1 + sigma3)/2 (mean stress, negative for compression)

F < 0: elastic, F = 0: yielding, F > 0: must redistribute.

In terms of stress components:
    R = sqrt(((sx-sy)/2)^2 + txy^2)
    p = (sx + sy) / 2
    F = R + p*sin(phi) - c*cos(phi)

The visco-plastic algorithm (Perzyna 1966):
    1. Compute elastic stresses from displacements.
    2. Check yield at each Gauss point.
    3. If yielded, compute excess stress and redistribute via body loads.
    4. Re-solve and iterate until convergence or iteration limit reached.
"""

from __future__ import annotations

import numpy as np


def principal_stresses(sigma: np.ndarray) -> tuple[float, float]:
    """Compute principal stresses from [sigma_x, sigma_y, tau_xy].

    Returns (sigma1, sigma3) where sigma1 >= sigma3 (algebraically).
    Compression-negative convention.
    """
    sx, sy, txy = sigma[0], sigma[1], sigma[2]
    mean = 0.5 * (sx + sy)
    radius = np.sqrt(0.25 * (sx - sy) ** 2 + txy ** 2)
    sigma1 = mean + radius  # Major (most tensile / least compressive)
    sigma3 = mean - radius  # Minor (most compressive)
    return sigma1, sigma3


def mohr_coulomb_F(sigma: np.ndarray, phi_rad: float, c: float) -> float:
    """Evaluate the Mohr-Coulomb failure function F.

    F = R + p*sin(phi) - c*cos(phi)
    where R = (sigma1 - sigma3)/2, p = (sigma1 + sigma3)/2

    F < 0: elastic (inside envelope)
    F = 0: yielding (on envelope)
    F > 0: stress state outside envelope (must redistribute)
    """
    s1, s3 = principal_stresses(sigma)
    sin_phi = np.sin(phi_rad)
    cos_phi = np.cos(phi_rad)
    R = 0.5 * (s1 - s3)
    p = 0.5 * (s1 + s3)
    F = R + p * sin_phi - c * cos_phi
    return F


def yield_correction(
    sigma: np.ndarray,
    phi_rad: float,
    c: float,
    psi_rad: float,
    E: float,
    v: float,
) -> np.ndarray:
    """Compute the visco-plastic stress correction for a yielded Gauss point.

    Uses the return mapping: delta_sigma = D @ dQ * F / (dF^T @ D @ dQ)
    where Q is the plastic potential (uses psi for non-associated flow).

    Returns:
        d_sigma: (3,) stress correction to subtract from current stress.
    """
    sx, sy, txy = sigma[0], sigma[1], sigma[2]
    sin_phi = np.sin(phi_rad)
    cos_phi = np.cos(phi_rad)
    sin_psi = np.sin(psi_rad)

    diff = sx - sy
    radius = np.sqrt(0.25 * diff ** 2 + txy ** 2)
    if radius < 1e-15:
        radius = 1e-15

    p = 0.5 * (sx + sy)
    F_val = radius + p * sin_phi - c * cos_phi
    if F_val <= 0.0:
        return np.zeros(3)

    # Derivatives of F w.r.t. stress components:
    # F = sqrt(((sx-sy)/2)^2 + txy^2) + (sx+sy)/2 * sin(phi) - c*cos(phi)
    #
    # dF/dsx = (sx-sy)/(4*R) + sin(phi)/2
    # dF/dsy = (sy-sx)/(4*R) + sin(phi)/2
    # dF/dtxy = txy/R

    dFdsx = 0.25 * diff / radius + 0.5 * sin_phi
    dFdsy = -0.25 * diff / radius + 0.5 * sin_phi
    dFdtxy = txy / radius
    dF = np.array([dFdsx, dFdsy, dFdtxy])

    # Plastic potential Q derivatives (replace phi with psi)
    dQdsx = 0.25 * diff / radius + 0.5 * sin_psi
    dQdsy = -0.25 * diff / radius + 0.5 * sin_psi
    dQdtxy = txy / radius
    dQ = np.array([dQdsx, dQdsy, dQdtxy])

    # D matrix for plane strain
    from slope64.elements import d_matrix_plane_strain
    D = d_matrix_plane_strain(E, v)

    # Stress correction: delta_sigma = D @ dQ * F / (dF^T @ D @ dQ)
    D_dQ = D @ dQ
    denom = dF @ D_dQ
    if abs(denom) < 1e-30:
        return np.zeros(3)

    d_sigma = D_dQ * F_val / denom
    return d_sigma


def factor_properties(
    phi_deg: float, c: float, fos: float
) -> tuple[float, float]:
    """Apply strength reduction factor to material properties.

    c_f = c / FOS
    phi_f = arctan(tan(phi) / FOS)

    Returns:
        (phi_f_deg, c_f): factored friction angle (degrees) and cohesion.
    """
    c_f = c / fos
    phi_rad = np.radians(phi_deg)
    phi_f_rad = np.arctan(np.tan(phi_rad) / fos)
    phi_f_deg = np.degrees(phi_f_rad)
    return phi_f_deg, c_f
