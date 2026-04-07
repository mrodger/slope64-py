#!/usr/bin/env python3
"""Validate all 7 examples against expected FOS values."""
import sys
import time
sys.path.insert(0, ".")

from slope64.parser import parse_dat
from slope64.srf import run_srf

EXPECTED = {
    "ex1": 1.56,
    "ex2": 1.22,
    "ex3": 2.00,
    "ex4": 1.27,
    "ex5": 1.83,
    "ex6": 1.03,
    "ex7": 1.41,
}

results = []
for name in sorted(EXPECTED):
    path = f"examples/{name}.dat"
    try:
        inp = parse_dat(path)
        t0 = time.time()
        res = run_srf(inp, verbose=True)
        elapsed = time.time() - t0
        diff = abs(res.fos - EXPECTED[name])
        status = "OK" if diff <= 0.02 else "MARGINAL" if diff <= 0.05 else "FAIL"
        print(f"FOS = {res.fos:.2f} (expected {EXPECTED[name]:.2f}, diff={diff:.2f}) [{status}] {elapsed:.1f}s\n")
        results.append((name, res.fos, EXPECTED[name], diff, status))
    except Exception as e:
        print(f"{name}: ERROR — {e}\n")
        results.append((name, None, EXPECTED[name], None, "ERROR"))

print("\n" + "="*60)
print(f"{'Example':>8s} {'Got':>6s} {'Expected':>8s} {'Diff':>6s} {'Status':>8s}")
print("-"*60)
for name, got, exp, diff, status in results:
    if got is not None:
        print(f"{name:>8s} {got:6.2f} {exp:8.2f} {diff:6.2f} {status:>8s}")
    else:
        print(f"{name:>8s}   ERR {exp:8.2f}    -   {status:>8s}")
print("="*60)
