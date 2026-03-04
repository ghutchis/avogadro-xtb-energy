#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""GFN-FF energy and gradient calculator.

GFN-FF ignores VERBOSITY_MUTED during Calculator construction and the first
singlepoint, so we redirect its C-level stdout to /dev/null for that phase to
prevent xtb's banner from polluting the AvogadroEnergy output stream.
"""

import json
import os
import sys
from contextlib import contextmanager

import numpy as np
from xtb.interface import Calculator, Param
from xtb.libxtb import VERBOSITY_MUTED

# Conversion factors
_ANG_TO_BOHR = 1.0 / 0.52917721067
_HARTREE_TO_KJ_MOL = 2625.4996394799
_HARTREE_BOHR_TO_KJ_MOL_ANG = _HARTREE_TO_KJ_MOL / _ANG_TO_BOHR

@contextmanager
def _suppress_c_stdout():
    """Redirect fd 1 to /dev/null for the duration of the block.

    GFN-FF writes initialization output via Fortran I/O which is fully-buffered
    when connected to a pipe. That buffer is only flushed at process exit (when
    Avogadro kills the process), so it never contaminates the energy/gradient
    output stream during a normal optimization run.
    """
    stdout_fd = sys.stdout.fileno()
    saved_fd = os.dup(stdout_fd)
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    try:
        sys.stdout.flush()
        os.dup2(devnull_fd, stdout_fd)
        yield
    finally:
        sys.stdout.flush()
        os.dup2(saved_fd, stdout_fd)
        os.close(saved_fd)
        os.close(devnull_fd)


def run():
    # Avogadro sends one compact JSON line on stdin, then keeps stdin open
    # for subsequent per-step coordinate updates (one line per atom: "x y z").
    bootstrap = json.loads(sys.stdin.readline())
    mol_cjson = bootstrap["cjson"]

    atoms = np.array(mol_cjson["atoms"]["elements"]["number"])
    coord_list = mol_cjson["atoms"]["coords"]["3d"]
    coordinates = np.array(coord_list, dtype=float).reshape(-1, 3)
    coordinates *= _ANG_TO_BOHR

    charge = bootstrap.get("charge", 0)
    # bootstrap["spin"] is spin multiplicity; xtb "uhf" is unpaired electrons
    uhf = bootstrap.get("spin", 1) - 1

    # Suppress C-level stdout during init: GFN-FF ignores VERBOSITY_MUTED here.
    with _suppress_c_stdout():
        calc = Calculator(Param.GFNFF, atoms, coordinates, charge=charge, uhf=uhf)
        calc.set_verbosity(VERBOSITY_MUTED)
        res = calc.singlepoint()

    # Loop forever — Avogadro kills the process when the optimization ends.
    while True:
        for i in range(len(atoms)):
            coordinates[i] = np.fromstring(input(), sep=" ", dtype=float)
        coordinates *= _ANG_TO_BOHR

        calc.update(coordinates)
        calc.singlepoint(res)

        print("AvogadroEnergy:", res.get_energy() * _HARTREE_TO_KJ_MOL)  
        print("AvogadroGradient:")
        grad = res.get_gradient() * _HARTREE_BOHR_TO_KJ_MOL_ANG
        output = np.array2string(grad)
        output = output.replace("[", "").replace("]", "")
        print(output)
        sys.stdout.flush()
