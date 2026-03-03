#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""Entry point for the avogadro-xtb-energy plugin.

Avogadro calls this as:
    avogadro-xtb-energy <identifier> [--lang <locale>] [--debug]

with the molecule bootstrap JSON on stdin (one compact JSON line).
"""

import argparse


def main():
    parser = argparse.ArgumentParser("avogadro-xtb-energy")
    parser.add_argument("feature")
    parser.add_argument("--lang", nargs="?", default="en")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    match args.feature:
        case "GFN1":
            from .gfn1 import run
            run()
        case "GFN2":
            from .gfn2 import run
            run()
        case "GFN-FF":
            from .gfnff import run
            run()
