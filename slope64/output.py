"""Output generation: results text and matplotlib plots.

Replaces the PostScript output of the original Slope64 with modern
matplotlib-based visualizations.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from slope64.parser import SlopeInput
from slope64.srf import SRFResult
from slope64.mesh import generate_mesh


def write_results(inp: SlopeInput, result: SRFResult, outpath: Path) -> None:
    """Write a .res text file with the analysis results."""
    with open(outpath, "w") as f:
        f.write(f"Slope64-py Results\n")
        f.write(f"==================\n")
        f.write(f"Title: {inp.title}\n\n")
        f.write(f"Geometry:\n")
        f.write(f"  w1={inp.w1}  s1={inp.s1}  w2={inp.w2}\n")
        f.write(f"  h1={inp.h1}  h2={inp.h2}\n")
        f.write(f"  Mesh: {inp.nx1+inp.nx2} x {inp.ny1+inp.ny2} elements\n\n")
        f.write(f"Materials ({inp.np_types} groups):\n")
        for i, mat in enumerate(inp.properties):
            f.write(f"  Group {i+1}: phi={mat[0]}, c={mat[1]}, psi={mat[2]}, "
                    f"gamma={mat[3]}, E={mat[4]}, v={mat[5]}\n")
        f.write(f"\nAnalysis:\n")
        f.write(f"  k_h = {inp.k_h}\n")
        f.write(f"  gam_w = {inp.gam_w}\n")
        f.write(f"  Iteration limit = {inp.limit}\n")
        f.write(f"  FOS tolerance = {inp.fos_tol}\n\n")
        f.write(f"{'SRF':>10s} {'max displ':>12s} {'iterations':>12s}\n")
        f.write(f"{'-'*10} {'-'*12} {'-'*12}\n")
        for trial in result.trials:
            f.write(f"{trial['fos']:10.4f} {trial['max_displ']:12.4e} {trial['iterations']:12d}\n")
        f.write(f"\nEstimated Factor of Safety = {result.fos:.2f}\n")


def plot_mesh(inp: SlopeInput, outpath: Path) -> None:
    """Plot the finite element mesh."""
    coords, connect, _ = generate_mesh(inp)
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    for el in range(connect.shape[0]):
        nodes = connect[el]
        # Draw element outline: 1-2-3-4-5-6-7-8-1
        outline = [nodes[0], nodes[1], nodes[2], nodes[3],
                   nodes[4], nodes[5], nodes[6], nodes[7], nodes[0]]
        xs = coords[outline, 0]
        ys = coords[outline, 1]
        ax.plot(xs, ys, "b-", linewidth=0.3)

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"FE Mesh: {inp.title}")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_srf_curve(result: SRFResult, outpath: Path) -> None:
    """Plot SRF vs. max displacement."""
    fos_vals = [t["fos"] for t in result.trials]
    dmax_vals = [t["max_displ"] for t in result.trials]

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    ax.plot(fos_vals, dmax_vals, "bo-")
    ax.axvline(result.fos, color="r", linestyle="--", label=f"FOS = {result.fos:.2f}")
    ax.set_xlabel("Strength Reduction Factor")
    ax.set_ylabel("Max Displacement")
    ax.set_title("SRF vs. Maximum Displacement")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_deformed_mesh(result: SRFResult, outpath: Path, scale: float = 0.0) -> None:
    """Plot the deformed mesh at failure (like original .dis output)."""
    if result.displacements is None or result.coords is None:
        return
    coords = result.coords
    connect = result.connect
    displ = result.displacements

    # Auto-scale deformation for visibility
    ux = displ[0::2]
    uy = displ[1::2]
    max_d = np.max(np.sqrt(ux**2 + uy**2))
    if scale == 0.0 and max_d > 1e-15:
        mesh_size = np.max(coords[:, 0]) - np.min(coords[:, 0])
        scale = 0.05 * mesh_size / max_d

    deformed = coords.copy()
    deformed[:, 0] += scale * ux
    deformed[:, 1] += scale * uy

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    for el in range(connect.shape[0]):
        nodes = connect[el]
        outline = [nodes[0], nodes[1], nodes[2], nodes[3],
                   nodes[4], nodes[5], nodes[6], nodes[7], nodes[0]]
        ax.plot(deformed[outline, 0], deformed[outline, 1], "b-", linewidth=0.3)

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Deformed Mesh at Failure (scale: {scale:.1f}x)")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_displacement_vectors(result: SRFResult, outpath: Path) -> None:
    """Plot nodal displacement vectors at failure (like original .vec output)."""
    if result.displacements is None or result.coords is None:
        return
    coords = result.coords
    displ = result.displacements
    ux = displ[0::2]
    uy = displ[1::2]

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    mag = np.sqrt(ux**2 + uy**2)
    # Only plot vectors with non-trivial displacement
    threshold = 0.01 * np.max(mag) if np.max(mag) > 0 else 0
    mask = mag > threshold

    ax.quiver(
        coords[mask, 0], coords[mask, 1],
        ux[mask], uy[mask],
        mag[mask], cmap="hot_r", scale_units="xy",
        angles="xy",
    )
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Displacement Vectors at Failure")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
