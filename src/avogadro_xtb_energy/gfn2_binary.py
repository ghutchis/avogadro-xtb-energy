#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""GFN2-xTB energy and gradient calculator using the binary protocol."""

import json
import sys

import numpy as np
from xtb.interface import Calculator, Param
from xtb.libxtb import VERBOSITY_MUTED

from .energy import EnergyServer

# Conversion factors
_ANG_TO_BOHR = 1.0 / 0.52917721067
_HARTREE_TO_KJ_MOL = 2625.4996394799
_HARTREE_BOHR_TO_KJ_MOL_ANG = _HARTREE_TO_KJ_MOL / _ANG_TO_BOHR


def run():
    # Avogadro sends one compact JSON line on stdin, then switches to binary
    bootstrap = json.loads(sys.stdin.readline())
    mol_cjson = bootstrap["cjson"]

    atoms = np.array(mol_cjson["atoms"]["elements"]["number"])
    coord_list = mol_cjson["atoms"]["coords"]["3d"]
    coordinates = np.array(coord_list, dtype=float).reshape(-1, 3)
    coordinates *= _ANG_TO_BOHR

    charge = bootstrap.get("charge", 0)
    # bootstrap["spin"] is spin multiplicity; xtb "uhf" is unpaired electrons
    uhf = bootstrap.get("spin", 1) - 1

    calc = Calculator(Param.GFN2xTB, atoms, coordinates, charge=charge, uhf=uhf)
    calc.set_verbosity(VERBOSITY_MUTED)
    res = calc.singlepoint()

    atom_count = len(atoms)

    with EnergyServer(sys.stdin.buffer, sys.stdout.buffer, atom_count) as server:
        for request in server.requests():
            coords_ang = request.coords  # shape (atom_count, 3), Angstrom
            coords_bohr = coords_ang * _ANG_TO_BOHR

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
