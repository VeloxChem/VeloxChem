from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.rsppolarizability import Polarizability
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver


@pytest.mark.solvers
class TestQrfFD:

    def test_qrf_fd(self):

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)

        molecule_string = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """
        xcfun_label = 'bp86'
        basis_set_label = 'def2-svp'
        scf_conv_thresh = 1.0e-8
        rsp_conv_thresh = 1.0e-6

        molecule = Molecule.read_str(molecule_string, units='au')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=ostream)

        # QR driver

        scf_settings = {'conv_thresh': scf_conv_thresh}
        rsp_settings = {'conv_thresh': rsp_conv_thresh, 'damping': 0}
        method_settings = {'xcfun': xcfun_label}

        scfdrv = ScfRestrictedDriver(comm, ostream)
        scfdrv.update_settings(scf_settings, method_settings)
        scf_result = scfdrv.compute(molecule, basis)

        quad_solver = QuadraticResponseDriver(comm, ostream)
        quad_solver.update_settings(rsp_settings, method_settings)
        quad_result = quad_solver.compute(molecule, basis, scf_result)

        if is_mpi_master():
            beta_zzz = -quad_result[(0, 0)].real

        # Finite difference

        scf_settings = {'conv_thresh': scf_conv_thresh}
        rsp_settings = {'conv_thresh': rsp_conv_thresh}
        method_settings = {'xcfun': xcfun_label}

        delta_ef = 1.0e-4
        efield_plus = [0.0, 0.0, delta_ef]
        efield_minus = [0.0, 0.0, -delta_ef]

        method_dict_plus = dict(method_settings)
        method_dict_minus = dict(method_settings)
        method_dict_plus['electric_field'] = efield_plus
        method_dict_minus['electric_field'] = efield_minus

        scf_drv_plus = ScfRestrictedDriver(comm, ostream)
        scf_drv_plus.update_settings(scf_settings, method_dict_plus)
        scf_drv_plus.compute(molecule, basis)

        lr_prop_plus = Polarizability(rsp_settings, method_dict_plus)
        lr_prop_plus.init_driver(comm, ostream)
        lr_prop_plus.compute(molecule, basis, scf_drv_plus.scf_tensors)

        scf_drv_minus = ScfRestrictedDriver(comm, ostream)
        scf_drv_minus.update_settings(scf_settings, method_dict_minus)
        scf_drv_minus.compute(molecule, basis)

        lr_prop_minus = Polarizability(rsp_settings, method_dict_minus)
        lr_prop_minus.init_driver(comm, ostream)
        lr_prop_minus.compute(molecule, basis, scf_drv_minus.scf_tensors)

        if is_mpi_master():
            rsp_func_plus = lr_prop_plus.get_property('response_functions')
            rsp_func_minus = lr_prop_minus.get_property('response_functions')
            alpha_zz_plus = -rsp_func_plus[('z', 'z', 0)]
            alpha_zz_minus = -rsp_func_minus[('z', 'z', 0)]
            beta_zzz_fd = (alpha_zz_plus - alpha_zz_minus) / (2.0 * delta_ef)

            rel_diff = abs(beta_zzz - beta_zzz_fd) / abs(beta_zzz_fd)
            assert rel_diff < 1.0e-6
