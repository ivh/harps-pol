#!/usr/bin/env python3
"""Command-line demodulation, for use without pyesorex.

Takes the S2D files of a 2- or 4-exposure polarimetric sequence, both fibres.
The fibre is taken from the filename suffix and the sequence position from the
retarder angle in the header, so the order of the arguments does not matter.
"""

import sys

from pyespdr.demod import (
    TAGS_A,
    TAGS_B,
    demodulate_cycles,
    read_s2d,
    split_cycles,
    stokes_parameter,
    write_products,
)


def fiber_of(filename):
    for tag in TAGS_A:
        if filename.endswith(f"_{tag}.fits"):
            return "A"
    for tag in TAGS_B:
        if filename.endswith(f"_{tag}.fits"):
            return "B"
    raise SystemExit(f"Cannot tell which fibre {filename!r} is")


def main(filenames):
    specs = [read_s2d(f, fiber_of(f)) for f in filenames]
    cycles_a = split_cycles([s for s in specs if s.fiber == "A"])
    cycles_b = split_cycles([s for s in specs if s.fiber == "B"])
    if len(cycles_a) != len(cycles_b):
        raise SystemExit(f"{len(cycles_a)} fibre A but {len(cycles_b)} "
                         f"fibre B cycles")

    stokes = stokes_parameter(cycles_a[0])
    print(f"Stokes {stokes}, {len(cycles_a)} cycle(s) of "
          f"{len(cycles_a[0])} exposures")
    for i, cycle in enumerate(cycles_a, start=1):
        print(f"  cycle {i}: " + ", ".join(f"{s.angle:g}" for s in cycle))

    products = demodulate_cycles(cycles_a, cycles_b, null=len(cycles_a[0]) == 4)
    inputs = [s for cycle in zip(cycles_a, cycles_b)
              for pair in zip(*cycle) for s in pair]
    for filename, catg in write_products(cycles_b[0][0], products, inputs, stokes):
        print(f"  {catg:15s} {filename}")


if __name__ == "__main__":
    if len(sys.argv) < 5 or (len(sys.argv) - 1) % 4:
        raise SystemExit("Need both fibres of 2 or 4 exposures per cycle, "
                         "so a multiple of 4 files")
    main(sys.argv[1:])
