"""Command-line interface for slope64-py."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from slope64.parser import parse_dat
from slope64.srf import run_srf
from slope64.output import (
    write_results, plot_mesh, plot_srf_curve,
    plot_deformed_mesh, plot_displacement_vectors,
)


def main():
    parser = argparse.ArgumentParser(description="Slope64-py: FEM slope stability analysis")
    parser.add_argument("datfile", help="Path to .dat input file or base name (e.g. ex1)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print iteration details")
    parser.add_argument("--no-plot", action="store_true", help="Skip plot generation")
    args = parser.parse_args()

    datpath = Path(args.datfile)
    if not datpath.exists():
        datpath = Path(args.datfile + ".dat")
    if not datpath.exists():
        datpath = Path("examples") / (args.datfile + ".dat")
    if not datpath.exists():
        print(f"Error: Cannot find '{args.datfile}'", file=sys.stderr)
        sys.exit(1)

    base = datpath.stem
    outdir = datpath.parent

    print(f"Slope64-py: {datpath.name}")
    inp = parse_dat(datpath)
    print(f"  Title: {inp.title}")
    print(f"  Mesh: {inp.nx_total} x {inp.ny_total} = {inp.n_elements} elements")
    print(f"  Materials: {inp.np_types} groups")
    print(f"  Running SRF analysis...")

    result = run_srf(inp, verbose=args.verbose)

    print(f"\n  Estimated Factor of Safety = {result.fos:.2f}")
    print(f"  Trials: {len(result.trials)}")

    # Write results
    res_path = outdir / f"{base}.res"
    write_results(inp, result, res_path)
    print(f"  Results: {res_path}")

    if not args.no_plot:
        mesh_path = outdir / f"{base}_mesh.png"
        srf_path = outdir / f"{base}_srf.png"
        dis_path = outdir / f"{base}_deformed.png"
        vec_path = outdir / f"{base}_vectors.png"

        plot_mesh(inp, mesh_path)
        plot_srf_curve(result, srf_path)
        plot_deformed_mesh(result, dis_path)
        plot_displacement_vectors(result, vec_path)

        print(f"  Plots: {mesh_path.name}, {srf_path.name}, {dis_path.name}, {vec_path.name}")


if __name__ == "__main__":
    main()
