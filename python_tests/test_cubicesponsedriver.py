from mpi4py import MPI
from pathlib import Path
import unittest
import pytest
import sys
import veloxchem as vlx

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.outputstream import OutputStream
from veloxchem.cubicresponsedriver import CubicResponseDriver
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestCrf(unittest.TestCase):
    """
    Tests the crf code 
    """

    def run_scf(self):

        molecule_string = """
        O   0.0   0.0   0.0
        H   .7586020000 0.0  -.5042840000
        H   .7586020000  0.0   .5042840000"""

        basis_set_label = 'aug-cc-pVDZ'

        scf_settings = {'conv_thresh': 1.0e-6}

        molecule = vlx.Molecule.read_str(molecule_string, units='ang')
        molecule.set_charge(0)
        molecule.set_multiplicity(1)
        method_settings = {}

        basis = vlx.MolecularBasis.read(molecule, basis_set_label)

        comm = MPI.COMM_WORLD
        ostream = vlx.OutputStream(sys.stdout)

        scf_drv = vlx.ScfRestrictedDriver(comm, ostream)
        scf_drv.update_settings(scf_settings, method_settings)
        scf_drv.compute(molecule, basis)

        return scf_drv.scf_tensors,molecule,basis

    def run_crf(self, w, ref_result):

        comm = MPI.COMM_WORLD
        ostream = vlx.OutputStream(sys.stdout)

        scf_tensors,molecule,ao_basis = self.run_scf()

        method_settings = {}

        rsp_settings = {'conv_thresh': 1.0e-8, 'b_frequencies': [0.2], 
                        'c_frequencies': [0.2], 'd_frequencies': [0.2],
                        'a_component': 'y','b_component': 'y','c_component': 'z',
                        'd_component': 'z', 'damping': 0.1}

        crf_prop = CubicResponseDriver(comm, ostream) 

        crf_prop.update_settings(rsp_settings, method_settings)

        crf_result_yyzz = crf_prop.compute(molecule, ao_basis, scf_tensors)
        
        # yyzz-component 

        self.assertTrue(abs(crf_result_yyzz[('T4',0.2)].real - ref_result['T4'].real) < 1.0e-4)

        self.assertTrue(abs(crf_result_yyzz[('T4',0.2)].imag - ref_result['T4'].imag) < 1.0e-4)
        
        self.assertTrue(abs(crf_result_yyzz[('X3',0.2)].real - ref_result['X3'].real) < 1.0e-4)

        self.assertTrue(abs(crf_result_yyzz[('X3',0.2)].imag - ref_result['X3'].imag) < 1.0e-4)

        self.assertTrue(abs(crf_result_yyzz[('A3',0.2)].real - ref_result['A3'].real) < 1.0e-4)

        self.assertTrue(abs(crf_result_yyzz[('A3',0.2)].imag - ref_result['A3'].imag) < 1.0e-4)

    def test_crf(self):

        w = 0.2

        # YYZZ component DALTON aug-cc-pVDZ
        ref_result = { 
            'T4': -2.65974591 + 11.26931601j,
             'X3': -3.53089306 + 17.45460417j,
             'A3': 10.88586429 + 3.80260981j
        }


        self.run_crf(w, ref_result)


if __name__ == '__main__':
    unittest.main()
