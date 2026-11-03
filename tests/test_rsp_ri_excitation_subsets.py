from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.lreigensolverunrest import LinearResponseUnrestrictedEigenSolver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.linearsolver import LinearSolver
from veloxchem import rijkresponse
from veloxchem import rijresponse


@pytest.mark.solvers
class TestResponseRiExcitationSubsets:
    """Core excitation and a restricted subspace with RI: the RI-JK linear
    response builds take the factors of their densities, which give the
    densities the solver forms, and the RI-J builds take the densities."""

    @staticmethod
    def get_water(charge=0, multiplicity=1):

        molecule = Molecule.read_molecule_string(
            """
            O  0.0           0.0  0.0
            H   .7586020000  0.0  -.5042840000
            H   .7586020000  0.0   .5042840000
            """,
            units='angstrom')
        molecule.set_charge(charge)
        molecule.set_multiplicity(multiplicity)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        return molecule, basis

    @staticmethod
    def record(monkeypatch, name):
        """Records the densities and factors handed to a Fock build, and
        whether the RI-JK build was called."""

        calls = []
        real = getattr(LinearSolver, name)

        def recorded(self, dens, *args, dens_factors=None, **kwargs):
            if self.rank == mpi_master():
                calls.append((dens, dens_factors))
            return real(self, dens, *args, dens_factors=dens_factors,
                        **kwargs)

        monkeypatch.setattr(LinearSolver, name, recorded)

        ri_jk_calls = []
        for fname in ['fock_matrices', 'fock_matrices_unrestricted']:
            real_ri = getattr(rijkresponse, fname)

            def counted(*args, real_ri=real_ri, **kwargs):
                ri_jk_calls.append(1)
                return real_ri(*args, **kwargs)

            monkeypatch.setattr(rijkresponse, fname, counted)

        return calls, ri_jk_calls

    @staticmethod
    def rebuilt(factors):
        """The densities the factors give: left times the first right factor
        transposed plus the second right factor times left transposed."""

        left, rights, transposed_rights = factors

        return [
            np.matmul(left, r1.T) + np.matmul(r2, left.T)
            for r1, r2 in zip(rights, transposed_rights)
        ]

    @pytest.mark.parametrize('mode', ['core_excitation', 'restricted_subspace'])
    def test_restricted(self, monkeypatch, mode):

        molecule, basis = self.get_water()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = 'b3lyp'
        scf_results = scf_drv.compute(molecule, basis)

        calls, ri_jk_calls = self.record(monkeypatch, '_comp_lr_fock')

        rpa_drv = LinearResponseEigenSolver()
        rpa_drv.ostream.mute()
        rpa_drv.nstates = 3
        rpa_drv.num_core_orbitals = 1
        if mode == 'core_excitation':
            rpa_drv.core_excitation = True
        else:
            rpa_drv.restricted_subspace = True
            rpa_drv.num_valence_orbitals = 2
            rpa_drv.num_virtual_orbitals = 6
        rpa_drv.compute(molecule, basis, scf_results)

        assert rpa_drv._uses_ri_jk()

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert calls
            assert len(ri_jk_calls) == len(calls)

            # NOTE: the left factor has the occupied orbitals the trial
            # vectors excite from: one core orbital, or the core one and two
            # valence ones.
            nexc = 1 if mode == 'core_excitation' else 3

            for dens, factors in calls:
                assert factors[0].shape[1] == nexc
                for density, from_factors in zip(dens, self.rebuilt(factors)):
                    assert np.max(np.abs(density - from_factors)) < 1.0e-12

    def test_unrestricted_core_excitation(self, monkeypatch):

        molecule, basis = self.get_water(1, 2)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = 'b3lyp'
        scf_results = scf_drv.compute(molecule, basis)

        calls, ri_jk_calls = self.record(monkeypatch,
                                         '_comp_lr_fock_unrestricted')

        rpa_drv = LinearResponseUnrestrictedEigenSolver()
        rpa_drv.ostream.mute()
        rpa_drv.nstates = 3
        rpa_drv.core_excitation = True
        rpa_drv.num_core_orbitals = 1
        rpa_drv.compute(molecule, basis, scf_results)

        assert rpa_drv._uses_ri_jk()

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert calls
            assert len(ri_jk_calls) == len(calls)

            for (dens_a, dens_b), (factors_a, factors_b) in calls:
                for dens, factors in [(dens_a, factors_a),
                                      (dens_b, factors_b)]:
                    assert factors[0].shape[1] == 1
                    for density, from_factors in zip(dens,
                                                     self.rebuilt(factors)):
                        assert np.max(np.abs(density -
                                             from_factors)) < 1.0e-12

    @pytest.mark.parametrize('ri', [False, True])
    @pytest.mark.parametrize('solver_type',
                             [LinearResponseEigenSolver, TdaEigenSolver])
    @pytest.mark.parametrize('mode', ['core_excitation', 'restricted_subspace'])
    def test_pure_functional(self, monkeypatch, mode, solver_type, ri):
        """Core excitation and a restricted subspace with a pure functional,
        whose RI build is RI-J and takes the densities themselves."""

        ref_exc_enes = {
            ('core_excitation', 'LinearResponseEigenSolver', False):
                [18.79713995, 18.90052947, 19.28291539],
            ('core_excitation', 'LinearResponseEigenSolver', True):
                [18.79714871, 18.90028327, 19.28300715],
            ('core_excitation', 'TdaEigenSolver', False):
                [18.79719415, 18.90054945, 19.28292031],
            ('core_excitation', 'TdaEigenSolver', True):
                [18.79720245, 18.90030319, 19.28301205],
            ('restricted_subspace', 'LinearResponseEigenSolver', False):
                [0.26000615, 0.36295389, 0.41575479],
            ('restricted_subspace', 'LinearResponseEigenSolver', True):
                [0.25991965, 0.36256058, 0.41569659],
            ('restricted_subspace', 'TdaEigenSolver', False):
                [0.26041902, 0.36298762, 0.41846255],
            ('restricted_subspace', 'TdaEigenSolver', True):
                [0.26033415, 0.36259436, 0.41840296],
        }[(mode, solver_type.__name__, ri)]

        molecule, basis = self.get_water()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = 'blyp'
        scf_drv.ri = ri
        scf_drv.conv_thresh = 1.0e-8
        scf_results = scf_drv.compute(molecule, basis)

        ri_j_calls = []
        real = rijresponse.coulomb_matrices

        def counted(*args, **kwargs):
            ri_j_calls.append(1)
            return real(*args, **kwargs)

        monkeypatch.setattr(rijresponse, 'coulomb_matrices', counted)

        rsp_drv = solver_type()
        rsp_drv.ostream.mute()
        rsp_drv.nstates = 3
        rsp_drv.conv_thresh = 1.0e-6
        rsp_drv.num_core_orbitals = 1
        if mode == 'core_excitation':
            rsp_drv.core_excitation = True
        else:
            rsp_drv.restricted_subspace = True
            rsp_drv.num_valence_orbitals = 2
            rsp_drv.num_virtual_orbitals = 6
        rsp_results = rsp_drv.compute(molecule, basis, scf_results)

        assert rsp_drv._uses_ri_coulomb() == ri
        assert bool(ri_j_calls) == ri

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(
                np.abs(rsp_results['eigenvalues'] -
                       np.array(ref_exc_enes))) < 1.0e-6
