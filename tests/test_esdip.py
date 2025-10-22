import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.firstorderprop import FirstOrderProperties
from veloxchem.tddftgradientdriver import TddftGradientDriver


@pytest.mark.solvers
class TestExcitedStateDipole:

    def test_excited_state_dipole(self):

        xyzstr = """3
        xyz
        O      0.0000        0.0000        0.0000
        H      0.0000        0.7408        0.5821
        H      0.0000       -0.7408        0.5821
        """

        mol = Molecule.read_xyz_string(xyzstr)
        bas = MolecularBasis.read(mol, 'def2-svp', verbose=False)

        xcfun_label = 'pbe0'
        nstates = 5

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        tddft_drv = LinearResponseEigenSolver()
        tddft_drv.ostream.mute()
        tddft_drv.nstates = nstates
        tddft_drv.xcfun = xcfun_label
        tddft_results = tddft_drv.compute(mol, bas, scf_results)

        tddft_grad_drv = TddftGradientDriver(scf_drv)
        tddft_grad_drv.ostream.mute()
        tddft_grad_drv.state_deriv_index = [s + 1 for s in range(nstates)]
        tddft_grad_drv.xcfun = xcfun_label
        es_densities = tddft_grad_drv.compute_excited_state_densities(
            mol, bas, tddft_results)

        unrelaxed_es_dipoles = []
        relaxed_es_dipoles = []

        for s in range(nstates):

            unrelaxed_prop = FirstOrderProperties()
            relaxed_prop = FirstOrderProperties()

            if scf_drv.rank == mpi_master():
                unrelaxed_es_dens = (es_densities['unrelaxed_density'][s] +
                                     scf_results['D_alpha'] +
                                     scf_results['D_beta'])
                relaxed_es_dens = (es_densities['relaxed_density'][s] +
                                   scf_results['D_alpha'] +
                                   scf_results['D_beta'])
            else:
                unrelaxed_es_dens = None
                relaxed_es_dens = None

            unrelaxed_prop.compute(mol, bas, unrelaxed_es_dens)
            relaxed_prop.compute(mol, bas, relaxed_es_dens)

            unrelaxed_es_dipoles.append(
                unrelaxed_prop.properties['dipole_moment'])
            relaxed_es_dipoles.append(relaxed_prop.properties['dipole_moment'])

        unrelaxed_es_dipoles = np.array(unrelaxed_es_dipoles)
        relaxed_es_dipoles = np.array(relaxed_es_dipoles)

        ref_unrelaxed_es_dipoles = np.array([
            [0.0, 0.0, -0.623834],
            [0.0, 0.0, -0.427791],
            [0.0, 0.0, -0.837253],
            [0.0, 0.0, -0.646247],
            [0.0, 0.0, -0.265976],
        ])

        ref_relaxed_es_dipoles = np.array([
            [0.0, 0.0, -0.327991],
            [0.0, 0.0, -0.181910],
            [0.0, 0.0, -0.406894],
            [0.0, 0.0, -0.329932],
            [0.0, 0.0, -0.117703],
        ])

        assert np.max(np.abs(ref_unrelaxed_es_dipoles -
                             unrelaxed_es_dipoles)) < 1.0e-4
        assert np.max(np.abs(ref_relaxed_es_dipoles -
                             relaxed_es_dipoles)) < 1.0e-4
