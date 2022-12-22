from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.rsppolarizability import Polarizability
from veloxchem.firstorderprop import FirstOrderProperties



@pytest.mark.solvers
class TestLrfFD:

    def test_lrf_fd(self):

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)

        molecule_string = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """
        xcfun_label = 'SCAN'
        basis_set_label = 'def2-svp'
        scf_conv_thresh = 1.0e-8
        rsp_conv_thresh = 1.0e-8

        molecule = Molecule.read_str(molecule_string, units='au')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=ostream)

        # LR driver

        scf_settings = {'conv_thresh': scf_conv_thresh}
        method_settings = {'xcfun': xcfun_label,'grid_level': 4}

        scfdrv = ScfRestrictedDriver(comm, ostream)
        scfdrv.update_settings(scf_settings, method_settings)
        scfdrv.compute(molecule, basis)

        method_dict = dict(method_settings)
        lr_prop = Polarizability({'frequencies': '0'}, method_dict)
        lr_prop.init_driver(comm, ostream)
        lr_prop.compute(molecule, basis,scfdrv.scf_tensors)
        
        if is_mpi_master():
            alpha_lr = -lr_prop.get_property('response_functions')[('z', 'z', 0)]

        # Finite difference

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

        scf_prop_plus = FirstOrderProperties(comm, ostream)
        scf_prop_plus.compute_scf_prop(molecule, basis, scf_drv_plus.scf_tensors)

        scf_drv_minus = ScfRestrictedDriver(comm, ostream)
        scf_drv_minus.update_settings(scf_settings, method_dict_minus)
        scf_drv_minus.compute(molecule, basis)

        scf_prop_minus = FirstOrderProperties(comm, ostream)
        scf_prop_minus.compute_scf_prop(molecule, basis, scf_drv_minus.scf_tensors)

        if is_mpi_master():
            mu_plus = scf_prop_plus.get_property('dipole moment')
            mu_minus = scf_prop_minus.get_property('dipole moment')
            alpha_zz_fd = (mu_plus[2] - mu_minus[2]) / (2.0 * delta_ef)

            rel_diff = abs(alpha_lr - alpha_zz_fd) / abs(alpha_zz_fd)
            assert rel_diff < 1.0e-6
