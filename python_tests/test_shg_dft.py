import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.shgdriver import SHGDriver
from veloxchem.veloxchemlib import is_mpi_master


@pytest.mark.solvers
class TestSHG:

    def run_scf(self):

        molecule_string = """
        H      1.2001      0.0363      0.8431
        C      0.7031      0.0083     -0.1305
        H      0.9877      0.8943     -0.7114
        H      1.0155     -0.8918     -0.6742
        O     -0.6582     -0.0067      0.1730
        H     -1.1326     -0.0311     -0.6482
        """

        basis_set_label = '6-31G'

        scf_settings = {'conv_thresh': 1.0e-6}

        molecule = Molecule.read_str(molecule_string, units='au')
        molecule.set_charge(0)
        molecule.set_multiplicity(1)
        method_settings = {'xcfun': 'b3lyp', 'grid_level': 4}

        basis = MolecularBasis.read(molecule, basis_set_label)

        scf_drv = ScfRestrictedDriver()
        scf_drv.update_settings(scf_settings, method_settings)
        scf_drv.compute(molecule, basis)

        return scf_drv.scf_tensors, molecule, basis

    def run_shg(self, ref_result):

        method_settings = {'xcfun': 'b3lyp', 'grid_level': 4}

        scf_tensors, molecule, ao_basis = self.run_scf()

        rsp_settings = {
            'conv_thresh': 1.0e-6,
            'frequencies': [0.1],
            'damping': 0.1
        }

        shg_prop = SHGDriver()
        shg_prop.update_settings(rsp_settings, method_settings)

        shg_result = shg_prop.compute(molecule, ao_basis, scf_tensors)

        if is_mpi_master():

            # x-component

            assert abs(shg_result[0.1][0].real - ref_result['x'].real) < 1.0e-5
            assert abs(shg_result[0.1][0].imag - ref_result['x'].imag) < 1.0e-5

            # y-component

            assert abs(shg_result[0.1][1].real - ref_result['y'].real) < 1.0e-5
            assert abs(shg_result[0.1][1].imag - ref_result['y'].imag) < 1.0e-5

            # z-component

            assert abs(shg_result[0.1][2].real - ref_result['z'].real) < 1.0e-5
            assert abs(shg_result[0.1][2].imag - ref_result['z'].imag) < 1.0e-5

    def test_shg(self):

        ref_result = {
            'x': -45.07384219269649-39.17912107231103j,
            'y': -0.10548318636851473-0.022340142881491942j,
            'z': 28.863052804345973+28.176366259870544j,
        }

        self.run_shg(ref_result)
