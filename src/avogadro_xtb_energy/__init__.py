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
    parser.add_argument("--protocol", nargs="?", default="text-v1")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    binary = args.protocol.startswith("binary")

    match args.feature:
        case "GFN1":
            if binary:
                from .gfn1_binary import run
            else:
                from .gfn1 import run
            run()
        case "GFN2":
            if binary:
                from .gfn2_binary import run
            else:
                from .gfn2 import run
            run()
        case "GFN-FF":
            if binary:
                from .gfnff_binary import run
            else:
                from .gfnff import run
            run()
