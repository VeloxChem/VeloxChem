from pathlib import Path
import sys
import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.rsppolarizability import Polarizability
from veloxchem.outputstream import OutputStream
from veloxchem.errorhandler import VeloxChemError


@pytest.mark.solvers
class TestPointCharges:

    @staticmethod
    def get_molecule_and_basis():

        xyz_string = """
        3
        xyz
        O    1.2361419   1.0137761  -0.0612424
        H    0.5104418   0.8944555   0.5514190
        H    1.9926927   1.1973129   0.4956931
        """
        mol = Molecule.read_xyz_string(xyz_string)
        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)
        return mol, bas

    def run_scf_with_point_charges(self, xcfun_label, ref_energy, ene_tol,
                                   ref_grad, grad_tol):

        mol, bas = self.get_molecule_and_basis()

        here = Path(__file__).parent
        potfile = str(here / 'data' / 'pe_water.pot')
        vdwfile = str(here / 'data' / 'pe_water.qm_vdw_params.txt')

        for scf_drv in [ScfRestrictedDriver(), ScfUnrestrictedDriver()]:

            scf_drv.xcfun = xcfun_label
            scf_drv.point_charges = potfile
            scf_drv.qm_vdw_params = vdwfile
            scf_drv.ostream.mute()

            scf_results = scf_drv.compute(mol, bas)

            if scf_drv.rank == mpi_master():
                assert abs(scf_results['scf_energy'] - ref_energy) < ene_tol

            grad_drv = ScfGradientDriver(scf_drv)
            grad_drv.compute(mol, bas, scf_results)
            grad = grad_drv.get_gradient()

            assert np.max(np.abs(grad - ref_grad)) < grad_tol

    def test_hf_with_point_charges(self):

        ref_energy = -75.9732334209
        ref_grad = np.array([
            [-0.012426949439, -0.009449823622, -0.013955922853],
            [-0.003881905016, -0.002485992070, 0.002890530027],
            [0.010380233206, 0.003716170298, 0.009359249589],
        ])

        self.run_scf_with_point_charges('hf', ref_energy, 1.0e-9, ref_grad,
                                        1.0e-6)

    def test_b3lyp_with_point_charges(self):

        ref_energy = -76.3690497939
        ref_grad = np.array([
            [-0.011154306762, -0.008028413518, 0.014921369370],
            [0.010672309330, -0.000331661909, -0.011756057376],
            [-0.006021550334, -0.000313431386, -0.004840032815],
        ])

        self.run_scf_with_point_charges('b3lyp', ref_energy, 1.0e-9, ref_grad,
                                        1.0e-4)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='requires standard single-process assertions')
    def test_invalid_point_charge_line_raises(self, tmp_path):

        mol, bas = self.get_molecule_and_basis()
        potfile = tmp_path / 'invalid_point_charges.pot'
        potfile.write_text('1\ninvalid\nQ 0.0 0.0 0.0\n')

        scf_drv = ScfRestrictedDriver()
        scf_drv.point_charges = str(potfile)
        scf_drv.ostream.mute()

        with pytest.raises(VeloxChemError,
                           match='potfile: Invalid data on point charge line 3'):
            scf_drv.compute(mol, bas)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='requires standard single-process assertions')
    def test_invalid_qm_vdw_line_raises(self, tmp_path):

        mol, bas = self.get_molecule_and_basis()
        potfile = tmp_path / 'valid_point_charges.pot'
        vdwfile = tmp_path / 'invalid_qm_vdw.txt'
        potfile.write_text('1\nvalid\nQ 0.0 0.0 0.0 0.1 1.0 2.0\n')
        vdwfile.write_text('1.0 2.0\n3.0\n')

        scf_drv = ScfRestrictedDriver()
        scf_drv.point_charges = str(potfile)
        scf_drv.qm_vdw_params = str(vdwfile)
        scf_drv.ostream.mute()

        with pytest.raises(VeloxChemError,
                           match='qm_vdw_params: Invalid data on line 2'):
            scf_drv.compute(mol, bas)

    def test_linear_response_with_point_charges(self):

        # point_charges only indirectly affect response

        mol, bas = self.get_molecule_and_basis()

        here = Path(__file__).parent
        potfile = str(here / 'data' / 'pe_water.pot')

        ref_energy_vac = -75.9610148052
        ref_energy_pc = -75.9798409316
        ref_alpha_zz_vac = 4.9696841693
        ref_alpha_zz_pc = 4.7783688210

        scf_settings = {'conv_thresh': 1.0e-8}
        rsp_settings = {'conv_thresh': 1.0e-5, 'frequencies': '0'}
        method_settings = {'xcfun': 'hf'}
        method_settings_pc = {'xcfun': 'hf', 'point_charges': potfile}

        scf_drv = ScfRestrictedDriver()
        scf_drv.update_settings(scf_settings, method_settings)
        scf_drv.ostream.mute()
        scf_results_vac = scf_drv.compute(mol, bas)

        scf_drv_pc = ScfRestrictedDriver()
        scf_drv_pc.update_settings(scf_settings, method_settings)
        scf_drv_pc.point_charges = potfile
        scf_drv_pc.ostream.mute()
        scf_results_pc = scf_drv_pc.compute(mol, bas)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert 'point_charges' in scf_results_pc
            assert abs(scf_results_vac['scf_energy'] -
                       ref_energy_vac) < 1.0e-8
            assert abs(scf_results_pc['scf_energy'] -
                       ref_energy_pc) < 1.0e-8

        lr_prop_vac = Polarizability(rsp_settings, method_settings)
        lr_prop_vac.init_driver(MPI.COMM_WORLD, OutputStream(None))
        lr_prop_vac.compute(mol, bas, scf_results_vac)

        # point_charges in the response method settings is ignored
        lr_prop_pc = Polarizability(rsp_settings, method_settings_pc)
        lr_prop_pc.init_driver(MPI.COMM_WORLD, OutputStream(None))
        lr_prop_pc.compute(mol, bas, scf_results_pc)

        lr_prop_pc_nokey = Polarizability(rsp_settings, method_settings)
        lr_prop_pc_nokey.init_driver(MPI.COMM_WORLD, OutputStream(None))
        lr_prop_pc_nokey.compute(mol, bas, scf_results_pc)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            alpha_zz_vac = -lr_prop_vac.get_property(
                'response_functions')[('z', 'z', 0)]
            alpha_zz_pc = -lr_prop_pc.get_property(
                'response_functions')[('z', 'z', 0)]
            alpha_zz_pc_nokey = -lr_prop_pc_nokey.get_property(
                'response_functions')[('z', 'z', 0)]

            assert abs(alpha_zz_vac - ref_alpha_zz_vac) < 1.0e-6
            assert abs(alpha_zz_pc - ref_alpha_zz_pc) < 1.0e-6
            assert abs(alpha_zz_pc_nokey - ref_alpha_zz_pc) < 1.0e-6

            # the response solver carries no point_charges attribute
            assert not hasattr(lr_prop_vac.rsp_driver, 'point_charges')
            assert not hasattr(lr_prop_pc.rsp_driver, 'point_charges')
            assert not hasattr(lr_prop_pc_nokey.rsp_driver, 'point_charges')

    def test_linear_response_with_electric_field(self):

        # electric_field is inherited from SCF results, and only indirectly
        # affect response

        mol, bas = self.get_molecule_and_basis()

        scf_settings = {'conv_thresh': 1.0e-8}
        rsp_settings = {'conv_thresh': 1.0e-5, 'frequencies': '0'}
        method_settings = {'xcfun': 'hf'}
        field = [0.0, 0.0, 0.001]

        scf_drv = ScfRestrictedDriver()
        scf_drv.update_settings(scf_settings, method_settings)
        scf_drv.electric_field = field
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert 'electric_field' in scf_results
            assert scf_results['electric_field'] == field

        lr_prop = Polarizability(rsp_settings, method_settings)
        lr_prop.init_driver(MPI.COMM_WORLD, OutputStream(None))
        lr_prop.compute(mol, bas, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert lr_prop.rsp_driver.electric_field == field

    def test_response_electric_field_conflict_warns(self, capsys):

        # a different electric field in response settings is overwritten
        # by the SCF value with a warning

        mol, bas = self.get_molecule_and_basis()

        scf_settings = {'conv_thresh': 1.0e-8}
        rsp_settings = {'conv_thresh': 1.0e-5, 'frequencies': '0'}
        method_settings = {'xcfun': 'hf'}
        field = [0.0, 0.0, 0.001]
        other_field = [0.0, 0.001, 0.0]

        scf_drv = ScfRestrictedDriver()
        scf_drv.update_settings(scf_settings, method_settings)
        scf_drv.electric_field = field
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        ostream = OutputStream(
            sys.stdout if MPI.COMM_WORLD.Get_rank() == mpi_master() else None)

        lr_prop = Polarizability(
            rsp_settings, {'xcfun': 'hf', 'electric_field': other_field})
        lr_prop.init_driver(MPI.COMM_WORLD, ostream)
        lr_prop.compute(mol, bas, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            # the SCF value wins
            assert lr_prop.rsp_driver.electric_field == field

        captured = capsys.readouterr()
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert 'electric_field' in captured.out
            assert 'overwritten' in captured.out

        # a matching field propagates silently
        lr_prop_match = Polarizability(
            rsp_settings, {'xcfun': 'hf', 'electric_field': field})
        lr_prop_match.init_driver(MPI.COMM_WORLD, ostream)
        lr_prop_match.compute(mol, bas, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert lr_prop_match.rsp_driver.electric_field == field

        captured = capsys.readouterr()
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert 'overwritten' not in captured.out
