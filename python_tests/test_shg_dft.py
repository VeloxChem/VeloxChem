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
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """

        basis_set_label = 'cc-pVDZ'

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

            assert abs(shg_result[0.1][0].real - ref_result['x'].real) < 1.0e-6
            assert abs(shg_result[0.1][0].imag - ref_result['x'].imag) < 1.0e-6

            # y-component

            assert abs(shg_result[0.1][1].real - ref_result['y'].real) < 1.0e-6
            assert abs(shg_result[0.1][1].imag - ref_result['y'].imag) < 1.0e-6

            # z-component

            assert abs(shg_result[0.1][2].real - ref_result['z'].real) < 1.0e-6
            assert abs(shg_result[0.1][2].imag - ref_result['z'].imag) < 1.0e-6

    def test_shg(self):

        ref_result = {
            'x': 0 + 0j,
            'y': 0 + 0j,
            'z': 100.98592528 + 56.88987624j,
        }

        self.run_shg(ref_result)
