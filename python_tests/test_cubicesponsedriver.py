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

        basis_set_label = '6-31G'

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

        rsp_settings = {'conv_thresh': 1.0e-4, 'b_frequencies': [0.2],'c_frequencies': [0.2],'damping': 0.1,'a_component':'x','b_component':'x','c_component':'x'}

        crf_prop = CubicResponseDriver(comm, ostream) 

        crf_prop.update_settings(rsp_settings, method_settings)

        crf_result_xxx = crf_prop.compute(molecule, ao_basis, scf_tensors)

        rsp_settings = {'conv_thresh': 1.0e-4, 'b_frequencies': [0.2],'c_frequencies': [0.2],'damping': 0.1,'a_component':'z','b_component':'z','c_component':'x'}

        crf_prop.update_settings(rsp_settings, method_settings)

        crf_result_zzx = crf_prop.compute(molecule, ao_basis, scf_tensors)

        rsp_settings = {'conv_thresh': 1.0e-4, 'b_frequencies': [0.2],'c_frequencies': [0.2],'damping': 0.1,'a_component':'y','b_component':'y','c_component':'x'}

        crf_prop.update_settings(rsp_settings, method_settings)

        crf_result_yyx = crf_prop.compute(molecule, ao_basis, scf_tensors)
        
        # x-component 

        self.assertTrue(abs(crf_result_xxx[0.2].real - ref_result['xxx'].real) < 1.0e-4)

        self.assertTrue(abs(crf_result_xxx[0.2].imag - ref_result['xxx'].imag) < 1.0e-4)
        
        # y-component 
        
        self.assertTrue(abs(crf_result_yyx[0.2].real - ref_result['yyx'].real) < 1.0e-4)

        self.assertTrue(abs(crf_result_yyx[0.2].imag - ref_result['yyx'].imag) < 1.0e-4)

        # z-component 
        
        self.assertTrue(abs(crf_result_zzx[0.2].real - ref_result['zzx'].real) < 1.0e-4)

        self.assertTrue(abs(crf_result_zzx[0.2].imag - ref_result['zzx'].imag) < 1.0e-4)


    def test_crf(self):

        w = 0.2

        ref_result = { 
            'xxx': 29.16175897 + 28.05788008j,
            'zzx': 27.37219617 + 32.23620966j,
            'yyx': -5.260931+ 2.081018j,
        }

        self.run_crf(w, ref_result)


if __name__ == '__main__':
    unittest.main()
