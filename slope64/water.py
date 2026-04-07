"""Pore pressure computation from a prescribed free surface.

The free surface is a polyline of (x, y) points defining the phreatic surface.
Pore pressure at any point below the surface is:
    u = gamma_w * (y_surface - y_point)
where y_surface is the interpolated surface elevation at x_point.

If the point is above the free surface, u = 0 (unsaturated).
"""

from __future__ import annotations

import numpy as np


def interpolate_surface(x: float, surf_points: list[tuple[float, float]]) -> float:
    """Linearly interpolate the free surface elevation at x.

    Args:
        x: horizontal coordinate to query.
        surf_points: list of (x, y) pairs defining the surface polyline.

    Returns:
        y_surface: interpolated surface elevation at x.
    """
    if len(surf_points) == 0:
        return 0.0

    # Sort by x (should already be sorted)
    pts = sorted(surf_points, key=lambda p: p[0])

    # Clamp to endpoints
    if x <= pts[0][0]:
        return pts[0][1]
    if x >= pts[-1][0]:
        return pts[-1][1]

    # Find segment
    for i in range(len(pts) - 1):
        x0, y0 = pts[i]
        x1, y1 = pts[i + 1]
        if x0 <= x <= x1:
            t = (x - x0) / (x1 - x0) if abs(x1 - x0) > 1e-15 else 0.0
            return y0 + t * (y1 - y0)

    return pts[-1][1]


def _soil_surface_elevation(
    x: float, w1: float, s1: float, w2: float, h1: float, h2: float,
) -> float:
    """Compute the soil surface elevation at horizontal coordinate x.

    The slope geometry:
        x in [0, w1]:         crest at y = h2 + h1
        x in [w1, w1+s1]:    slope face, linear from h2+h1 down to h2
        x in [w1+s1, w1+s1+w2]: foundation surface at y = h2
    """
    crest_y = h2 + h1
    toe_x = w1 + s1
    if x <= w1:
        return crest_y
    elif x <= toe_x:
        frac = (x - w1) / s1 if s1 > 1e-15 else 1.0
        return crest_y - frac * h1
    else:
        return h2


def interpolate_surface_gradient(
    x: float, surf_points: list[tuple[float, float]]
) -> float:
    """Return dy_surf_rel/dx (piecewise-constant slope of the surface polyline).

    Used to compute the horizontal seepage force = -gamma_w * dy_water/dx.
    Returns 0.0 for a horizontal water table or outside the defined range.
    """
    if len(surf_points) < 2:
        return 0.0

    pts = sorted(surf_points, key=lambda p: p[0])

    def _slope(p0: tuple[float, float], p1: tuple[float, float]) -> float:
        dx = p1[0] - p0[0]
        return (p1[1] - p0[1]) / dx if abs(dx) > 1e-15 else 0.0

    if x <= pts[0][0]:
        return _slope(pts[0], pts[1])
    if x >= pts[-1][0]:
        return _slope(pts[-2], pts[-1])

    for i in range(len(pts) - 1):
        if pts[i][0] <= x <= pts[i + 1][0]:
            return _slope(pts[i], pts[i + 1])

    return 0.0


def pore_pressure_at_point(
    x: float, y: float,
    surf_points: list[tuple[float, float]],
    gam_w: float,
    h2: float,
    h1: float,
    w1: float = 0.0,
    s1: float = 0.0,
    w2: float = 0.0,
) -> float:
    """Compute pore pressure at a point (x, y) for yield check.

    The free surface y-coordinates in the .dat file are relative to the
    crest of the embankment (y=0 at crest, negative = below crest).

    When the water surface is above the soil surface (submerged slope),
    the pore pressure is capped at the soil surface elevation. This is
    because gravity loads only include soil self-weight; the water weight
    above the soil surface creates equal total stress and pore pressure
    increments that cancel in the effective stress.

    Args:
        x, y: point coordinates in mesh system (y=0 at base).
        surf_points: (x, y_relative) pairs from the .dat file.
        gam_w: unit weight of water.
        h2: foundation thickness.
        h1: embankment height.
        w1, s1, w2: slope geometry for soil surface computation.

    Returns:
        u: pore pressure (positive = compression).
    """
    if gam_w <= 0 or len(surf_points) == 0:
        return 0.0

    # Convert surface y from crest-relative to mesh coordinates
    crest_y = h2 + h1
    y_surf_rel = interpolate_surface(x, surf_points)
    y_water = crest_y + y_surf_rel

    depth = y_water - y
    if depth <= 0:
        return 0.0

    return gam_w * depth
