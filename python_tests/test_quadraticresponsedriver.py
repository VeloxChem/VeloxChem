import pytest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestQrf:

    def run_scf(self):

        molecule_string = """
            O  0.0           0.0  0.0
            H   .7586020000  0.0  -.5042840000
            H   .7586020000  0.0   .5042840000
        """

        basis_set_label = '6-31G'

        scf_settings = {'conv_thresh': 1.0e-6}

        molecule = Molecule.read_str(molecule_string, units='ang')
        molecule.set_charge(0)
        molecule.set_multiplicity(1)

        basis = MolecularBasis.read(molecule, basis_set_label)

        scf_drv = ScfRestrictedDriver()
        scf_drv.update_settings(scf_settings)
        scf_drv.compute(molecule, basis)

        return scf_drv.scf_tensors, molecule, basis

    def run_qrf(self, ref_result):

        scf_tensors, molecule, ao_basis = self.run_scf()

        rsp_settings = {
            'conv_thresh': 1.0e-4,
            'b_frequencies': [0.2],
            'c_frequencies': [0.2],
            'damping': 0.1,
            'a_components': 'x',
            'b_components': 'x',
            'c_components': 'x'
        }

        qrf_prop = QuadraticResponseDriver()
        qrf_prop.update_settings(rsp_settings)
        qrf_result_xxx = qrf_prop.compute(molecule, ao_basis, scf_tensors, method_settings)

        rsp_settings = {
            'conv_thresh': 1.0e-4,
            'b_frequencies': [0.2],
            'c_frequencies': [0.2],
            'damping': 0.1,
            'a_components': 'z',
            'b_components': 'z',
            'c_components': 'x'
        }

        qrf_prop.update_settings(rsp_settings)
        qrf_result_zzx = qrf_prop.compute(molecule, ao_basis, scf_tensors, method_settings)

        rsp_settings = {
            'conv_thresh': 1.0e-4,
            'b_frequencies': [0.2],
            'c_frequencies': [0.2],
            'damping': 0.1,
            'a_components': 'y',
            'b_components': 'y',
            'c_components': 'x'
        }

# master branch version
#        qrf_prop.update_settings(rsp_settings)
#        qrf_result_yyx = qrf_prop.compute(molecule, ao_basis, scf_tensors)
#
#        if is_mpi_master():
#
#            # x-component
#
#            assert abs(qrf_result_xxx[0.2].real -
#                       ref_result['xxx'].real) < 1.0e-4
#            assert abs(qrf_result_xxx[0.2].imag -
#                       ref_result['xxx'].imag) < 1.0e-4
#
#            # y-component
#
#            assert abs(qrf_result_yyx[0.2].real -
#                       ref_result['yyx'].real) < 1.0e-4
#            assert abs(qrf_result_yyx[0.2].imag -
#                       ref_result['yyx'].imag) < 1.0e-4
# qrf_dft branch version
        qrf_result_yyx = qrf_prop.compute(molecule, ao_basis, scf_tensors,method_settings)

        # x-component

        self.assertTrue(abs(qrf_result_xxx[(0.2,0.2)].real - ref_result['xxx'].real) < 1.0e-4)

        self.assertTrue(abs(qrf_result_xxx[(0.2,0.2)].imag - ref_result['xxx'].imag) < 1.0e-4)

        # y-component

        self.assertTrue(abs(qrf_result_yyx[(0.2,0.2)].real - ref_result['yyx'].real) < 1.0e-4)

        self.assertTrue(abs(qrf_result_yyx[(0.2,0.2)].imag - ref_result['yyx'].imag) < 1.0e-4)

        # z-component

        self.assertTrue(abs(qrf_result_zzx[(0.2,0.2)].real - ref_result['zzx'].real) < 1.0e-4)

        self.assertTrue(abs(qrf_result_zzx[(0.2,0.2)].imag - ref_result['zzx'].imag) < 1.0e-4)
#>>>>>>> qrf_dft

            # z-component

            assert abs(qrf_result_zzx[0.2].real -
                       ref_result['zzx'].real) < 1.0e-4
            assert abs(qrf_result_zzx[0.2].imag -
                       ref_result['zzx'].imag) < 1.0e-4

    def test_qrf(self):

        ref_result = {
            'xxx': 29.16175897 + 28.05788008j,
            'zzx': 27.37219617 + 32.23620966j,
            'yyx': -5.260931 + 2.081018j,
        }

        self.run_qrf(ref_result)


if __name__ == '__main__':
    unittest.main()
