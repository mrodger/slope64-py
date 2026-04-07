"""Parse Slope64 .dat input files.

The .dat format uses quoted string labels followed by values on the next line(s).
Labels are either single-quoted (title) or double-quoted (parameter descriptions).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class SlopeInput:
    """Parsed Slope64 input data."""

    title: str = ""

    # Geometry
    w1: float = 0.0  # Width of top of embankment
    s1: float = 0.0  # Width of sloping portion
    w2: float = 0.0  # Foundation extends right of toe
    h1: float = 0.0  # Height of embankment
    h2: float = 0.0  # Thickness of foundation layer

    # Mesh
    nx1: int = 0  # x-elements in embankment
    nx2: int = 0  # x-elements right of toe
    ny1: int = 0  # y-elements in embankment
    ny2: int = 0  # y-elements in foundation

    # Materials
    np_types: int = 1  # Number of property groups
    # Each row: [phi, c, psi, gamma, E, v]
    properties: list[list[float]] = field(default_factory=list)
    # Element type assignment grid (ny_total rows x nx_total cols), 1-indexed
    etype: list[list[int]] = field(default_factory=list)

    # Seismic
    k_h: float = 0.0  # Horizontal acceleration factor

    # Water
    nosurf: int = 0  # Number of free surface points
    surf: list[tuple[float, float]] = field(default_factory=list)  # (x, y) pairs
    gam_w: float = 0.0  # Unit weight of water

    # Solver
    limit: int = 500  # Iteration ceiling
    fos_tol: float = 0.02  # FOS accuracy tolerance

    @property
    def nx_total(self) -> int:
        return self.nx1 + self.nx2

    @property
    def ny_total(self) -> int:
        return self.ny1 + self.ny2

    @property
    def n_elements(self) -> int:
        return self.nx_total * self.ny_total


def parse_dat(filepath: str | Path) -> SlopeInput:
    """Parse a Slope64 .dat file and return a SlopeInput object."""
    filepath = Path(filepath)
    text = filepath.read_text(encoding="utf-8", errors="replace")
    lines = text.strip().splitlines()

    inp = SlopeInput()
    i = 0

    def skip_blank():
        nonlocal i
        while i < len(lines) and lines[i].strip() == "":
            i += 1

    def current_line() -> str:
        return lines[i].strip() if i < len(lines) else ""

    def read_value_line() -> str:
        nonlocal i
        skip_blank()
        val = current_line()
        i += 1
        return val

    def skip_label():
        """Skip a quoted label line."""
        nonlocal i
        skip_blank()
        line = current_line()
        if line.startswith('"') or line.startswith("'"):
            i += 1

    # Title (single-quoted)
    skip_blank()
    line = current_line()
    if line.startswith("'"):
        inp.title = line.strip("'").strip()
        i += 1
    elif line.startswith('"'):
        inp.title = line.strip('"').strip()
        i += 1

    # w1
    skip_label()
    inp.w1 = float(read_value_line())

    # s1
    skip_label()
    inp.s1 = float(read_value_line())

    # w2
    skip_label()
    inp.w2 = float(read_value_line())

    # h1
    skip_label()
    inp.h1 = float(read_value_line())

    # h2
    skip_label()
    inp.h2 = float(read_value_line())

    # nx1
    skip_label()
    inp.nx1 = int(read_value_line())

    # nx2
    skip_label()
    inp.nx2 = int(read_value_line())

    # ny1
    skip_label()
    inp.ny1 = int(read_value_line())

    # ny2
    skip_label()
    inp.ny2 = int(read_value_line())

    # np_types
    skip_label()
    inp.np_types = int(read_value_line())

    # Material properties
    skip_label()
    for _ in range(inp.np_types):
        vals = [float(x) for x in read_value_line().split()]
        inp.properties.append(vals)

    # Element type assignment
    skip_label()
    if inp.np_types > 1:
        # Embankment rows (ny1) have nx1 elements each.
        # Foundation rows (ny2) have (nx1+nx2) elements each.
        total_etype_count = inp.ny1 * inp.nx1 + inp.ny2 * (inp.nx1 + inp.nx2)
        all_etypes = []
        while len(all_etypes) < total_etype_count:
            skip_blank()
            line = current_line()
            if line.startswith('"') or line.startswith("'"):
                break
            vals = [int(x) for x in line.split()]
            all_etypes.extend(vals)
            i += 1
        # Reshape: first ny1 rows width nx1, then ny2 rows width (nx1+nx2)
        idx = 0
        for row_idx in range(inp.ny1 + inp.ny2):
            row_width = inp.nx1 if row_idx < inp.ny1 else (inp.nx1 + inp.nx2)
            row = all_etypes[idx : idx + row_width]
            inp.etype.append(row)
            idx += row_width

    # k_h
    skip_label()
    inp.k_h = float(read_value_line())

    # Free surface
    skip_label()
    inp.nosurf = int(read_value_line())
    for _ in range(inp.nosurf):
        vals = [float(x) for x in read_value_line().split()]
        inp.surf.append((vals[0], vals[1]))

    # gam_w
    skip_label()
    inp.gam_w = float(read_value_line())

    # Iteration ceiling
    skip_label()
    inp.limit = int(read_value_line())

    # FOS tolerance
    skip_label()
    inp.fos_tol = float(read_value_line())

    return inp
