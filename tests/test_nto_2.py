from pathlib import Path
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecularorbitals import MolecularOrbitals
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.tdaeigensolverunrest import TdaUnrestrictedEigenSolver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.lreigensolverunrest import LinearResponseUnrestrictedEigenSolver
from veloxchem.resultsio import read_results
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis


@pytest.mark.solvers
class TestNTO2:

    def run_nto(self, mol, basis, flag, core_excitation, restricted_subspace):

        fname = '_vlx_test_nto_2_'

        # run SCF

        if mol.get_multiplicity() == 1:
            scf_drv = ScfRestrictedDriver()
        else:
            scf_drv = ScfUnrestrictedDriver()

        scf_drv.restart = False

        scf_drv.filename = fname
        scf_drv.ostream.mute()

        scf_results = scf_drv.compute(mol, basis)

        # run TDA/RPA

        if flag == 'tda':
            if mol.get_multiplicity() == 1:
                rsp_drv = TdaEigenSolver()
            else:
                rsp_drv = TdaUnrestrictedEigenSolver()
        elif flag == 'rpa':
            if mol.get_multiplicity() == 1:
                rsp_drv = LinearResponseEigenSolver()
            else:
                rsp_drv = LinearResponseUnrestrictedEigenSolver()

        rsp_drv.nstates = 3
        rsp_drv.nto = True

        if core_excitation:
            rsp_drv.core_excitation = True
            rsp_drv.num_core_orbitals = 1
        elif restricted_subspace:
            rsp_drv.restricted_subspace = True
            rsp_drv.num_core_orbitals = 1
            rsp_drv.num_valence_orbitals = 1
            rsp_drv.num_virtual_orbitals = 5
            rsp_drv.nstates = 8

        rsp_drv.filename = fname
        rsp_drv.ostream.mute()

        rsp_results = rsp_drv.compute(mol, basis, scf_results)

        if scf_drv.rank == mpi_master():
            rsp_results_2 = read_results(f'{fname}.h5', 'rsp')
        else:
            rsp_results_2 = {}

        for idx in range(rsp_drv.nstates):
            state_label = f'S{idx + 1}'

            nto = rsp_drv.get_nto(mol, basis, scf_results, rsp_results,
                                  state_label)
            nto_2 = rsp_drv.get_nto(mol, basis, scf_results, rsp_results_2,
                                    state_label)

            if scf_drv.rank == mpi_master():
                nto_3 = MolecularOrbitals.read_hdf5(
                    f'{fname}.h5', label=f'rsp/nto/NTO_{state_label}_')

                assert np.max(
                    np.abs(nto.alpha_to_numpy() -
                           nto_2.alpha_to_numpy())) < 1.0e-10
                assert np.max(
                    np.abs(nto.alpha_to_numpy() -
                           nto_3.alpha_to_numpy())) < 1.0e-10

                assert np.max(
                    np.abs(nto.beta_to_numpy() -
                           nto_2.beta_to_numpy())) < 1.0e-10
                assert np.max(
                    np.abs(nto.beta_to_numpy() -
                           nto_3.beta_to_numpy())) < 1.0e-10

                assert np.max(
                    np.abs(nto.occa_to_numpy() -
                           nto_2.occa_to_numpy())) < 1.0e-10
                assert np.max(
                    np.abs(nto.occa_to_numpy() -
                           nto_3.occa_to_numpy())) < 1.0e-10

                assert np.max(
                    np.abs(nto.occb_to_numpy() -
                           nto_2.occb_to_numpy())) < 1.0e-10
                assert np.max(
                    np.abs(nto.occb_to_numpy() -
                           nto_3.occb_to_numpy())) < 1.0e-10

                assert np.max(np.abs(nto.ea_to_numpy() -
                                     nto_2.ea_to_numpy())) < 1.0e-10
                assert np.max(np.abs(nto.ea_to_numpy() -
                                     nto_3.ea_to_numpy())) < 1.0e-10

                assert np.max(np.abs(nto.eb_to_numpy() -
                                     nto_2.eb_to_numpy())) < 1.0e-10
                assert np.max(np.abs(nto.eb_to_numpy() -
                                     nto_3.eb_to_numpy())) < 1.0e-10

        # clean up
        if scf_drv.rank == mpi_master():
            for fname in [f'{fname}.h5', f'{fname}_scf.h5', f'{fname}_rsp.h5']:
                fpath = Path(fname)
                if fpath.is_file():
                    fpath.unlink()
        scf_drv.comm.barrier()

    @pytest.mark.parametrize('multiplicity', [1, 3])
    @pytest.mark.parametrize('flag', ['rpa', 'tda'])
    @pytest.mark.parametrize('core_excitation_str', ['no_cvs', 'cvs'])
    @pytest.mark.parametrize('restricted_subspace_str', ['no_rsa', 'rsa'])
    def test_nto(self, multiplicity, flag, core_excitation_str,
                 restricted_subspace_str):

        core_excitation = (core_excitation_str == 'cvs')
        restricted_subspace = (restricted_subspace_str == 'rsa')

        if (multiplicity == 3) and restricted_subspace:
            pytest.skip('Restricted subspace not available for open-shell systems')

        if core_excitation and restricted_subspace:
            pytest.skip('Core excitation and restricted subspace are mutually exclusive')

        mol = Molecule.read_smiles('C=C')
        mol.set_multiplicity(multiplicity)

        basis = MolecularBasis.read(mol, 'def2-svp')

        self.run_nto(mol, basis, flag, core_excitation, restricted_subspace)
