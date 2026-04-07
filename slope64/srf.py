"""Strength Reduction Factor (SRF) bisection driver.

Finds the Factor of Safety by progressively increasing the SRF until
the visco-plastic solver fails to converge (slope failure).

The bisection starts with a coarse search (0.5, 1.0, 1.5, ...) to bracket
the failure point, then refines with bisection to the requested tolerance.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from slope64.parser import SlopeInput
from slope64.mesh import generate_mesh
from slope64.solver import solve_step


@dataclass
class SRFResult:
    """Result of the SRF analysis."""
    fos: float = 0.0
    converged: bool = True
    trials: list[dict] = field(default_factory=list)
    displacements: np.ndarray | None = None
    coords: np.ndarray | None = None
    connect: np.ndarray | None = None
    etype_map: np.ndarray | None = None


def run_srf(inp: SlopeInput, verbose: bool = False) -> SRFResult:
    """Run the SRF bisection to find the Factor of Safety.

    The algorithm:
    1. Coarse search: try FOS = 0.5, 1.0, 1.5, ... until non-convergence.
    2. Bracket the failure: last converging FOS is lower bound, first
       non-converging is upper bound.
    3. Bisect the bracket until (upper - lower) <= fos_tol.
    4. Return lower bound as the estimated FOS.
    """
    coords, connect, etype_map = generate_mesh(inp)
    result = SRFResult(coords=coords, connect=connect, etype_map=etype_map)

    if verbose:
        n_elem = connect.shape[0]
        n_nodes = coords.shape[0]
        print(f"Mesh: {n_elem} elements, {n_nodes} nodes, {2 * n_nodes} DOFs")
        print(f"{'trial factor':>14s} {'max displ':>14s} {'iterations':>12s}")

    def try_fos(fos: float) -> tuple[float, int, bool]:
        """Try a single FOS value. Returns (max_displ, iters, converged)."""
        displ, max_d, iters, conv = solve_step(
            coords, connect, etype_map, inp, fos, inp.limit
        )
        trial = {
            "fos": fos,
            "max_displ": max_d,
            "iterations": iters,
            "converged": conv,
        }
        result.trials.append(trial)
        if conv:
            result.displacements = displ
        if verbose:
            print(f"{fos:14.4f} {max_d:14.4e} {iters:12d}")
        return max_d, iters, conv

    # Phase 1: Coarse search with 0.5 increments
    lower = None
    upper = None

    fos_values = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]
    for fos in fos_values:
        _, _, conv = try_fos(fos)
        if conv:
            lower = fos
        else:
            upper = fos
            break

    if upper is None:
        result.fos = fos_values[-1]
        result.converged = True
        if verbose:
            print(f"\nSlope did not fail up to FOS = {fos_values[-1]}")
        return result

    if lower is None:
        result.fos = 0.5
        result.converged = False
        if verbose:
            print(f"\nSlope fails even at FOS = 0.5")
        return result

    # Phase 2: Bisection refinement
    while (upper - lower) > inp.fos_tol:
        mid = 0.5 * (upper + lower)
        _, _, conv = try_fos(mid)
        if conv:
            lower = mid
        else:
            upper = mid

    result.fos = round(0.5 * (lower + upper), 2)
    result.converged = True

    if verbose:
        print(f"\nEstimated Factor of Safety = {result.fos:.2f}")

    return result
