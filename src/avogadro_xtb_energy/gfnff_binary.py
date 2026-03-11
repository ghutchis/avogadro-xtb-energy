#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""GFN-FF energy and gradient calculator using the binary protocol.

GFN-FF ignores VERBOSITY_MUTED during Calculator construction and the first
singlepoint, so we redirect its C-level stdout to /dev/null for that phase to
prevent xtb's banner from polluting the output stream.
"""

import json
import os
import sys
from contextlib import contextmanager

import numpy as np
from xtb.interface import Calculator, Param
from xtb.libxtb import VERBOSITY_MUTED

from .energy import EnergyServer

# Conversion factors
_ANG_TO_BOHR = 1.0 / 0.52917721067
_HARTREE_TO_KJ_MOL = 2625.4996394799
_HARTREE_BOHR_TO_KJ_MOL_ANG = _HARTREE_TO_KJ_MOL / _ANG_TO_BOHR


@contextmanager
def _suppress_c_stdout():
    """Redirect fd 1 to /dev/null for the duration of the block."""
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
    bootstrap = json.loads(sys.stdin.buffer.readline())
    mol_cjson = bootstrap["cjson"]

    atoms = np.array(mol_cjson["atoms"]["elements"]["number"])
    coord_list = mol_cjson["atoms"]["coords"]["3d"]
    coordinates = np.array(coord_list, dtype=float).reshape(-1, 3)
    coordinates *= _ANG_TO_BOHR

    charge = bootstrap.get("charge", 0)
    uhf = bootstrap.get("spin", 1) - 1

    with _suppress_c_stdout():
        calc = Calculator(Param.GFNFF, atoms, coordinates, charge=charge, uhf=uhf)
        calc.set_verbosity(VERBOSITY_MUTED)
        res = calc.singlepoint()

    atom_count = len(atoms)

    with EnergyServer(sys.stdin.buffer, sys.stdout.buffer, atom_count) as server:
        for request in server.requests():
            coords_bohr = request.coords * _ANG_TO_BOHR

            calc.update(coords_bohr)
            calc.singlepoint(res)

            if request.wants_energy_and_gradient:
                energy = res.get_energy() * _HARTREE_TO_KJ_MOL
                grad = res.get_gradient() * _HARTREE_BOHR_TO_KJ_MOL_ANG
                request.send_energy_and_gradient(energy, grad)
            elif request.wants_gradient:
                grad = res.get_gradient() * _HARTREE_BOHR_TO_KJ_MOL_ANG
                request.send_gradient(grad)
            else:
                energy = res.get_energy() * _HARTREE_TO_KJ_MOL
                request.send_energy(energy)
