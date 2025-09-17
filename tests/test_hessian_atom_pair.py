from pathlib import Path
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfhessiandriver import ScfHessianDriver


class TestHessianAtomPair:

    @pytest.mark.solvers
    def test_hessian_atom_pair(self):

        mol = Molecule.read_xyz_string("""15

        O       99.814000000   100.835000000   101.232000000
        H       99.329200000    99.976800000   101.063000000
        H       99.151600000   101.561000000   101.414000000
        O       98.804000000    98.512200000    97.758100000
        H       99.782100000    98.646900000    97.916700000
        H       98.421800000    99.326500000    97.321300000
        O       98.070300000    98.516900000   100.438000000
        H       97.172800000    98.878600000   100.690000000
        H       98.194000000    98.592200000    99.448100000
        O      102.360000000   101.551000000    99.964500000
        H      102.675000000   102.370000000   100.444000000
        H      101.556000000   101.180000000   100.430000000
        O      101.665000000    98.316100000    98.319400000
        H      101.904000000    99.233800000    98.002000000
        H      102.224000000    97.640900000    97.837700000""")

        basis = MolecularBasis.read(mol, 'sto-3g')

        atom_pairs = [(0, 3), (6, 9)]

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.compute(mol, basis)

        hess_drv = ScfHessianDriver(scf_drv)
        hess_drv.atom_pairs = atom_pairs
        hess_drv.compute(mol, basis)

        pair_hess = hess_drv.hessian

        if scf_drv.rank == mpi_master():
            cwd = Path(__file__).parent
            path = str(cwd / 'data' / 'water_hessian_pair_1-4_7-10.txt')
            reference_hess = np.loadtxt(path)
            np.testing.assert_allclose(pair_hess,
                                       reference_hess,
                                       rtol=1e-8,
                                       atol=1e-10)
