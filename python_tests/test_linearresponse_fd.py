from mpi4py import MPI

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.outputstream import OutputStream
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.rsppolarizability import Polarizability
from veloxchem.firstorderprop import FirstOrderProperties


class TestLrfFD:

    def run_lrf_fd(self, xcfun_label):

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)

        molecule_string = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """
        basis_set_label = 'def2-svp'
        scf_conv_thresh = 1.0e-8
        rsp_conv_thresh = 1.0e-5

        molecule = Molecule.read_str(molecule_string, units='au')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        # LR driver

        scf_settings = {'conv_thresh': scf_conv_thresh}
        rsp_settings = {'conv_thresh': rsp_conv_thresh, 'frequencies': '0'}
        method_settings = {'xcfun': xcfun_label}

        scfdrv = ScfRestrictedDriver(comm, ostream)
        scfdrv.update_settings(scf_settings, method_settings)
        scf_results = scfdrv.compute(molecule, basis)

        lr_prop = Polarizability(rsp_settings, method_settings)
        lr_prop.init_driver(comm, ostream)
        lr_prop.compute(molecule, basis, scf_results)

        if is_mpi_master():
            alpha_zz = -lr_prop.get_property('response_functions')[('z', 'z',
                                                                    0)]

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
        scf_results_plus = scf_drv_plus.compute(molecule, basis)

        scf_prop_plus = FirstOrderProperties(comm, ostream)
        scf_prop_plus.compute_scf_prop(molecule, basis, scf_results_plus)

        scf_drv_minus = ScfRestrictedDriver(comm, ostream)
        scf_drv_minus.update_settings(scf_settings, method_dict_minus)
        scf_results_minus = scf_drv_minus.compute(molecule, basis)

        scf_prop_minus = FirstOrderProperties(comm, ostream)
        scf_prop_minus.compute_scf_prop(molecule, basis, scf_results_minus)

        if is_mpi_master():
            mu_plus = scf_prop_plus.get_property('dipole moment')
            mu_minus = scf_prop_minus.get_property('dipole moment')
            alpha_zz_fd = (mu_plus[2] - mu_minus[2]) / (2.0 * delta_ef)

            rel_diff = abs(alpha_zz - alpha_zz_fd) / abs(alpha_zz_fd)
            assert rel_diff < 1.0e-6

    def test_lda_lrf_fd(self):

        self.run_lrf_fd('slater')

    def test_gga_lrf_fd(self):

        self.run_lrf_fd('pbe0')

    def test_mgga_lrf_fd(self):

        self.run_lrf_fd('scan')
