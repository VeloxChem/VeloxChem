from mpi4py import MPI
import sys
import pytest

from veloxchem.veloxchemlib import is_single_node
from veloxchem.cubicresponsedriver import CubicResponseDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.skipif(not is_single_node(),
                    reason='This test only runs on single node')
class TestCrf:

    def run_scf(self):

        molecule_string = """
        O   0.0   0.0   0.0
        H   .7586020000 0.0  -.5042840000
        H   .7586020000  0.0   .5042840000"""

        basis_set_label = 'aug-cc-pVDZ'

        scf_settings = {'conv_thresh': 1.0e-8}

        molecule = Molecule.read_str(molecule_string, units='ang')
        molecule.set_charge(0)
        molecule.set_multiplicity(1)
        method_settings = {}

        basis = MolecularBasis.read(molecule, basis_set_label)

        comm = MPI.COMM_WORLD
        ostream = OutputStream(sys.stdout)

        scf_drv = ScfRestrictedDriver(comm, ostream)
        scf_drv.update_settings(scf_settings, method_settings)
        scf_drv.compute(molecule, basis)

        return scf_drv.scf_tensors, molecule, basis

    def run_crf(self, ref_result):

        comm = MPI.COMM_WORLD
        ostream = OutputStream(sys.stdout)

        scf_tensors, molecule, ao_basis = self.run_scf()

        method_settings = {}
        wb = -0.1
        wc = 0.3
        wd = -0.4
        rsp_settings = {
            'conv_thresh': 1.0e-8,
            'b_frequencies': [wb],
            'c_frequencies': [wc],
            'd_frequencies': [wd],
            'a_component': 'y',
            'b_component': 'y',
            'c_component': 'z',
            'd_component': 'z',
            'damping': 0.1
        }

        crf_prop = CubicResponseDriver(comm, ostream)

        crf_prop.update_settings(rsp_settings, method_settings)

        crf_result_yyzz = crf_prop.compute(molecule, ao_basis, scf_tensors)

        assert abs(crf_result_yyzz[('T4', wb, wc, wd)].real -
                   ref_result['T4'].real) < 1.0e-6

        assert abs(crf_result_yyzz[('T4', wb, wc, wd)].imag -
                   ref_result['T4'].imag) < 1.0e-6

        assert abs(crf_result_yyzz[('X3', wb, wc, wd)].real -
                   ref_result['X3'].real) < 1.0e-6

        assert abs(crf_result_yyzz[('X3', wb, wc, wd)].imag -
                   ref_result['X3'].imag) < 1.0e-6

        assert abs(crf_result_yyzz[('A3', wb, wc, wd)].real -
                   ref_result['A3'].real) < 1.0e-6

        assert abs(crf_result_yyzz[('A3', wb, wc, wd)].imag -
                   ref_result['A3'].imag) < 1.0e-6

    def test_crf(self):

        ref_result = {
            'E3': 0.33616818 - 24.75661969j,
            'T4': 8.58028589 - 21.58902361j,
            'X2': 211.57812207 + 56.41923653j,
            'X3': 35.38383531 - 40.21250147j,
            'A2': -262.74089736 + 405.99704670j,
            'A3': 24.59268588 - 10.94409800j,
        }

        self.run_crf(ref_result)
