from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.lrsolver import LinearResponseSolver


class TestH2Se:

    def test_h2se_scf(self):

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'h2se.inp'
        outfile = inpfile.with_suffix('.out')

        task = MpiTask([str(inpfile), str(outfile)])

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        e_scf = scf_drv.get_scf_energy()

        scf_drv.restart = True
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        e_scf_restart = scf_drv.get_scf_energy()

        if is_mpi_master(task.mpi_comm):
            scf_h5 = Path(task.input_dict['scf']['checkpoint_file'])
            if scf_h5.is_file():
                scf_h5.unlink()
            scf_final_h5 = scf_h5.with_suffix('.results.h5')
            if scf_final_h5.is_file():
                scf_final_h5.unlink()

        assert abs(-2400.70461320 - e_scf) < 1.0e-8
        assert abs(-2400.70461320 - e_scf_restart) < 1.0e-8

        # unrestricted scf

        task.molecule.set_charge(1)
        task.molecule.set_multiplicity(2)
        task.molecule.check_multiplicity()

        scf_unrest_drv = ScfUnrestrictedDriver(task.mpi_comm, task.ostream)
        scf_unrest_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        e_uhf = scf_unrest_drv.get_scf_energy()
        assert abs(-2400.38319890 - e_uhf) < 1.0e-8

        if is_mpi_master(task.mpi_comm):
            s2 = scf_unrest_drv.compute_s2(task.molecule,
                                           scf_unrest_drv.scf_tensors['S'],
                                           scf_unrest_drv.mol_orbs)
            assert abs(0.7619 - s2) < 1.0e-4

        task.finish()

    def test_h2se_rsp(self):

        # scf
        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'h2se.inp')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        scf_tensors = scf_drv.scf_tensors

        # TDA
        tda_exci = TDAExciDriver(task.mpi_comm, task.ostream)

        tda_exci.update_settings({
            'nstates': 3,
            'eri_thresh': 1.0e-12,
            'conv_thresh': 1.0e-4,
        })

        tda_results = tda_exci.compute(task.molecule, task.ao_basis,
                                       scf_tensors)

        if is_mpi_master(task.mpi_comm):

            reigs = tda_results['eigenvalues']
            osc_strs = tda_results['oscillator_strengths']
            trans_dipoles = tda_results['electric_transition_dipoles']

            ref_eigs = np.array([0.207436, 0.257474, 0.368358])
            ref_osc_strs = np.array([0.0000, 0.0003, 0.2797])
            ref_trans_dipoles = [
                np.array([0.0, 0.0, 0.0]),
                np.array([-0.043213, 0.0, 0.0]),
                np.array([0.0, -0.754589, 0.754589]),
            ]

            assert np.max(np.abs(reigs - ref_eigs)) < 1.0e-6
            assert np.max(np.abs(osc_strs - ref_osc_strs)) < 1.0e-4

            for td, ref_td in zip(trans_dipoles, ref_trans_dipoles):
                prefac = 1.0 if np.dot(td, ref_td) >= 0.0 else -1.0
                assert np.max(np.abs(td - ref_td * prefac)) < 2.0e-4

        # RPA
        lreig = LinearResponseEigenSolver(task.mpi_comm, task.ostream)

        lreig.update_settings({
            'nstates': 3,
            'eri_thresh': 1.0e-12,
            'conv_thresh': 1.0e-4,
        })

        rpa_results = lreig.compute(task.molecule, task.ao_basis, scf_tensors)

        if is_mpi_master(task.mpi_comm):

            reigs = rpa_results['eigenvalues']
            osc_strs = rpa_results['oscillator_strengths']
            trans_dipoles = rpa_results['electric_transition_dipoles']

            ref_eigs = np.array([0.20565979, 0.25474355, 0.36246841])
            ref_osc_strs = np.array([0.0000, 0.0012, 0.2477])
            ref_trans_dipoles = [
                np.array([0.0, 0.0, 0.0]),
                np.array([0.084626, 0.0, 0.0]),
                np.array([0.0, -0.715935, 0.715935])
            ]

            assert np.max(np.abs(reigs - ref_eigs)) < 1.0e-6
            assert np.max(np.abs(osc_strs - ref_osc_strs)) < 1.0e-4

            for td, ref_td in zip(trans_dipoles, ref_trans_dipoles):
                prefac = 1.0 if np.dot(td, ref_td) >= 0.0 else -1.0
                assert np.max(np.abs(td - ref_td * prefac)) < 2.0e-4

        # polarizability
        lr_solver = LinearResponseSolver(task.mpi_comm, task.ostream)

        lr_solver.update_settings({
            'frequencies': '0, 0.1',
            'eri_thresh': '1.0e-12',
            'conv_thresh': '1.0e-4',
        })

        lr_results = lr_solver.compute(task.molecule, task.ao_basis,
                                       scf_tensors)

        if is_mpi_master(task.mpi_comm):

            ref_polar = {
                ('x', 'x', 0): 15.26732,
                ('x', 'y', 0): 0.0,
                ('x', 'z', 0): 0.0,
                ('y', 'x', 0): 0.0,
                ('y', 'y', 0): 23.22477,
                ('y', 'z', 0): -0.759275,
                ('z', 'x', 0): 0.0,
                ('z', 'y', 0): -0.759275,
                ('z', 'z', 0): 23.22477,
                ('x', 'x', 0.1): 15.63176,
                ('x', 'y', 0.1): 0.0,
                ('x', 'z', 0.1): 0.0,
                ('y', 'x', 0.1): 0.0,
                ('y', 'y', 0.1): 24.27442,
                ('y', 'z', 0.1): -0.85701,
                ('z', 'x', 0.1): 0.0,
                ('z', 'y', 0.1): -0.85701,
                ('z', 'z', 0.1): 24.27442,
            }

            for key, val in lr_results['response_functions'].items():
                ref = -ref_polar[key]
                diff = abs(val - ref)
                assert diff < 1.0e-5

        task.finish()
