# Slope64 — Python Rebuild

An attempt to convert D.V. Griffiths' Slope64 finite element slope stability program from Fortran to python, with the goal of making it more accessible to people who need it but can't run the original.
Still a work in progress, I have an alternative  for now.

---

## Background

Slope64 is a finite element slope stability program written in Fortran by **D.V. Griffiths** (Colorado School of Mines), based on the visco-plastic algorithm described in:

> Griffiths, D.V. & Lane, P.A. (1999). *Slope stability analysis by finite elements.* Géotechnique, 49(3), 387–403.

> Smith, I.M. & Griffiths, D.V. (1998). *Programming the Finite Element Method*, 3rd ed. Wiley.

The original program is distributed as Fortran 77/90 source and compiles cleanly on Linux with `gfortran`. It reads `.dat` input files and writes `.res` output with the computed Factor of Safety. The algorithm is elegant and well-validated — the problem is that most geotechnical practitioners working in Python (or without a Fortran toolchain) can't easily use it.

This is a conversion effort, not a reimagining. The goal is to preserve the algorithm and `.dat` file format exactly, so that existing input files and published benchmark results work without modification.

---

## What it does

Slope64 uses the **strength reduction method (SRM)** with a visco-plastic FEM solver to find the Factor of Safety (FoS) for a slope:

1. Reduce the material strength parameters (c′, φ′) by a trial factor
2. Apply gravity loading and solve the elastic system
3. Check Mohr-Coulomb yield at every Gauss point
4. For yielded points, compute stress corrections and iterate until convergence
5. Repeat with increasing reduction factors until the mesh can no longer converge (collapse)
6. The FoS is the last converged reduction factor

The mesh is generated internally from the `.dat` geometry — no external mesher required.

---

## Status

Working conversion. All 7 published benchmark examples from Griffiths & Lane (1999) are validated:

| Example | Expected FoS | Python result | Match |
|---------|-------------|---------------|-------|
| ex1 | 1.56 | 1.56 | ✓ |
| ex2 | 1.22 | 1.22 | ✓ |
| ex3 | 2.00 | 2.00 | ✓ |
| ex4 | 1.27 | 1.27 | ✓ |
| ex5 | 1.83 | 1.83 | ✓ |
| ex6 | 1.03 | 1.03 | ✓ |
| ex7 | 1.41 | 1.41 | ✓ |

Run the validation yourself:
```bash
python validate_all.py
```

---

## Install

```bash
git clone https://github.com/mrodger/slope64-py
cd slope64-py
pip install -e .
```

Requires Python 3.12+, numpy, scipy, matplotlib.

---

## Usage

```bash
# Run a single analysis
slope64 examples/ex1.dat

# With output plots
slope64 examples/ex1.dat --plot

# Python API
from slope64 import run
fos = run("examples/ex1.dat")
print(f"FoS = {fos:.2f}")
```

---

## Input format

The `.dat` format is unchanged from the original Fortran program — any existing Slope64 input files should work directly. Example:

```
'Example 1 - simple slope'
"Number of elements in x, y"
 10  5
"Element size (m)"
 1.0  1.0
"Slope geometry: toe_x crest_x height"
 3  7  5
"Material properties: E nu phi_deg c_kPa gamma"
 1e5  0.3  20.0  10.0  20.0
"Water table (0=none, 1=full)"
 0
"Convergence: max_iter tolerance"
 500  0.001
"SRF range and increment"
 1.0  3.0  0.01
```

See `examples/` for the 7 published benchmark cases.

---

## Project structure

```
slope64-py/
├── slope64/
│   ├── parser.py     # .dat file parser
│   ├── mesh.py       # mesh generation (8-node quads, Smith & Griffiths convention)
│   ├── elements.py   # shape functions, B-matrix, D-matrix, Gauss integration
│   ├── material.py   # Mohr-Coulomb yield, strength reduction
│   ├── water.py      # pore water pressure / water table
│   ├── solver.py     # visco-plastic iterative solver
│   ├── srf.py        # strength reduction factor loop
│   ├── output.py     # results and plots
│   └── cli.py        # command-line entry point
├── examples/         # 7 benchmark .dat files + expected .res outputs
├── tests/            # pytest unit tests
└── validate_all.py   # run all 7 benchmarks and compare FoS values
```

---

## Credit

All algorithmic credit belongs to **D.V. Griffiths** and co-authors. This repository is a Python translation of his work, intended to make the method available to a wider audience. The original Fortran program and associated publications remain the authoritative reference.

If you use this in research, cite the original:

```
Griffiths, D.V. & Lane, P.A. (1999). Slope stability analysis by finite elements.
Géotechnique, 49(3), 387–403. https://doi.org/10.1680/geot.1999.49.3.387
```

---

## License

MIT — but see the credit section above. The underlying algorithm and method are from Griffiths (1999).
