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

        molecule = Molecule.read_molecule_string(molecule_string,
                                                 units='angstrom')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        return molecule, basis, scf_results

    def run_qrf(self, ref_result):

        molecule, basis, scf_results = self.run_scf()

        qrf_drv = QuadraticResponseDriver()
        qrf_drv.ostream.mute()

        rsp_settings = {
            'conv_thresh': 1.0e-4,
            'b_frequencies': [0.2],
            'c_frequencies': [0.2],
            'damping': 0.1,
            'a_components': 'x',
            'b_components': 'x',
            'c_components': 'x'
        }
        qrf_drv.update_settings(rsp_settings)
        qrf_result_xxx = qrf_drv.compute(molecule, basis, scf_results)

        rsp_settings = {
            'conv_thresh': 1.0e-4,
            'b_frequencies': [0.2],
            'c_frequencies': [0.2],
            'damping': 0.1,
            'a_components': 'z',
            'b_components': 'z',
            'c_components': 'x'
        }
        qrf_drv.update_settings(rsp_settings)
        qrf_result_zzx = qrf_drv.compute(molecule, basis, scf_results)

        rsp_settings = {
            'conv_thresh': 1.0e-4,
            'b_frequencies': [0.2],
            'c_frequencies': [0.2],
            'damping': 0.1,
            'a_components': 'y',
            'b_components': 'y',
            'c_components': 'x'
        }
        qrf_drv.update_settings(rsp_settings)
        qrf_result_yyx = qrf_drv.compute(molecule, basis, scf_results)

        if is_mpi_master():
            thresh = 1.0e-4

            # x-component
            assert abs(qrf_result_xxx[(0.2, 0.2)].real -
                       ref_result['xxx'].real) < thresh
            assert abs(qrf_result_xxx[(0.2, 0.2)].imag -
                       ref_result['xxx'].imag) < thresh

            # y-component
            assert abs(qrf_result_yyx[(0.2, 0.2)].real -
                       ref_result['yyx'].real) < thresh
            assert abs(qrf_result_yyx[(0.2, 0.2)].imag -
                       ref_result['yyx'].imag) < thresh

            # z-component
            assert abs(qrf_result_zzx[(0.2, 0.2)].real -
                       ref_result['zzx'].real) < thresh
            assert abs(qrf_result_zzx[(0.2, 0.2)].imag -
                       ref_result['zzx'].imag) < thresh

    def test_qrf(self):

        ref_result = {
            'xxx': 29.16175897 + 28.05788008j,
            'zzx': 27.37219617 + 32.23620966j,
            'yyx': -5.26093100 + 2.08101800j,
        }

        self.run_qrf(ref_result)
