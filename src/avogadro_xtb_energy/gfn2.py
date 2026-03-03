#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""GFN2-xTB energy and gradient calculator."""

import json
import sys

import numpy as np
from xtb.interface import Calculator, Param
from xtb.libxtb import VERBOSITY_MUTED

# Angstrom → Bohr conversion factor
_ANG_TO_BOHR = 1.0 / 0.52917721067


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

    calc = Calculator(Param.GFN2xTB, atoms, coordinates, charge=charge, uhf=uhf)
    calc.set_verbosity(VERBOSITY_MUTED)
    res = calc.singlepoint()

    # Loop forever — Avogadro kills the process when the optimization ends.
    while True:
        for i in range(len(atoms)):
            coordinates[i] = np.fromstring(input(), sep=" ", dtype=float)
        coordinates *= _ANG_TO_BOHR

        calc.update(coordinates)
        calc.singlepoint(res)

        print("AvogadroEnergy:", res.get_energy())  # Hartree
        print("AvogadroGradient:")
        grad = res.get_gradient() * 4961.475  # Hartree/Bohr → kJ/mol/Å
        output = np.array2string(grad)
        output = output.replace("[", "").replace("]", "")
        print(output)
        sys.stdout.flush()
