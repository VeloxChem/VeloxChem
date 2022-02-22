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
            'damping': 0.1,
        }

        shg_prop = SHGDriver()
        shg_prop.update_settings(rsp_settings, method_settings)

        shg_result = shg_prop.compute(molecule, ao_basis, scf_tensors)

        if is_mpi_master():
            freq = list(shg_result['beta'].keys())[0]

            for ind, comp in enumerate('xyz'):
                assert abs(shg_result['beta'][freq][ind].real -
                           ref_result[comp].real) < 1.0e-5
                assert abs(shg_result['beta'][freq][ind].imag -
                           ref_result[comp].imag) < 1.0e-5

            assert abs(shg_result['beta_bar'][freq].real -
                       ref_result['beta_bar'].real) < 1.0e-4
            assert abs(shg_result['beta_bar'][freq].imag -
                       ref_result['beta_bar'].imag) < 1.0e-4

    def test_shg(self):

        ref_result = {
            'x': -45.07384219 - 39.17912107j,
            'y': -0.10548319 - 0.02234014j,
            'z': 28.86305280 + 28.17636626j,
            'beta_bar': -12.60290150 - 13.87516029j,
        }

        self.run_shg(ref_result)
