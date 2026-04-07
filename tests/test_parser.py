"""Tests for the .dat file parser."""

from pathlib import Path

from slope64.parser import parse_dat

EXAMPLES = Path(__file__).parent.parent / "examples"


def test_parse_ex1():
    inp = parse_dat(EXAMPLES / "ex1.dat")
    assert inp.title == "Example 1: A homogeneous slope"
    assert inp.w1 == 12.0
    assert inp.s1 == 20.0
    assert inp.w2 == 12.0
    assert inp.h1 == 10.0
    assert inp.h2 == 10.0
    assert inp.nx1 == 32
    assert inp.nx2 == 12
    assert inp.ny1 == 10
    assert inp.ny2 == 10
    assert inp.np_types == 1
    assert len(inp.properties) == 1
    assert inp.properties[0] == [30.0, 5.0, 0.0, 20.0, 1e5, 0.3]
    assert inp.etype == []  # np_types=1, no etype data
    assert inp.k_h == 0.0
    assert inp.nosurf == 0
    assert inp.gam_w == 0.0
    assert inp.limit == 1000
    assert inp.fos_tol == 0.02


def test_parse_ex2_two_layer():
    inp = parse_dat(EXAMPLES / "ex2.dat")
    assert inp.np_types == 2
    assert len(inp.properties) == 2
    assert inp.properties[0][0] == 25.0  # phi group 1
    assert inp.properties[1][0] == 15.0  # phi group 2
    assert len(inp.etype) == 10  # ny1+ny2 = 5+5
    # First row should be all 1s (embankment)
    assert inp.etype[0] == [1, 1, 1, 1, 1]


def test_parse_ex4_free_surface():
    inp = parse_dat(EXAMPLES / "ex4.dat")
    assert inp.nosurf == 22
    assert len(inp.surf) == 22
    assert inp.surf[0] == (0.0, -2.40)
    assert inp.gam_w == 62.4  # Imperial units


def test_parse_ex5_submerged():
    inp = parse_dat(EXAMPLES / "ex5.dat")
    assert inp.nosurf == 2
    assert inp.surf[0] == (0.0, 2.0)  # Above crest
    assert inp.gam_w == 9.81


def test_parse_ex7_seismic():
    inp = parse_dat(EXAMPLES / "ex7.dat")
    assert inp.k_h == 0.25
    assert inp.limit == 2000
    assert inp.fos_tol == 0.02


def test_parse_all_examples():
    """Ensure all 7 examples parse without error."""
    for i in range(1, 8):
        inp = parse_dat(EXAMPLES / f"ex{i}.dat")
        assert inp.title != ""
        assert inp.nx1 > 0
        assert inp.ny1 > 0
        assert len(inp.properties) == inp.np_types
