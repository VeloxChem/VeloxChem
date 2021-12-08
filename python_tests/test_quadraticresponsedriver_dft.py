from mpi4py import MPI
from pathlib import Path
import unittest
import pytest
import sys
import veloxchem as vlx

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.outputstream import OutputStream
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestQrf(unittest.TestCase):
    """
    Tests the qrf code 
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
        method_settings = {'xcfun': 'SLATER', 'grid_level': 6}

        basis = vlx.MolecularBasis.read(molecule, basis_set_label)

        comm = MPI.COMM_WORLD
        ostream = vlx.OutputStream(sys.stdout)

        scf_drv = vlx.ScfRestrictedDriver(comm, ostream)
        scf_drv.update_settings(scf_settings, method_settings)
        scf_drv.compute(molecule, basis)

        return scf_drv.scf_tensors,molecule,basis

    def run_qrf(self, ref_result):

        comm = MPI.COMM_WORLD
        ostream = vlx.OutputStream(sys.stdout)

        scf_tensors,molecule,ao_basis = self.run_scf()

        method_settings = {'xcfun': 'SLATER', 'grid_level': 6}

        rsp_settings = {'conv_thresh': 1.0e-6, 'b_frequencies': [-0.1],'c_frequencies': [0.3],'damping': 0.1,'a_component':'x','b_component':'x','c_component':'x'}

        qrf_prop = QuadraticResponseDriver(comm, ostream) 

        qrf_prop.update_settings(rsp_settings, method_settings)

        qrf_result_xxx = qrf_prop.compute(molecule, ao_basis, scf_tensors,method_settings)

        rsp_settings = {'conv_thresh': 1.0e-6, 'b_frequencies': [-0.1],'c_frequencies': [0.3],'damping': 0.1,'a_component':'z','b_component':'z','c_component':'x'}

        qrf_prop.update_settings(rsp_settings, method_settings)

        qrf_result_zzx = qrf_prop.compute(molecule, ao_basis, scf_tensors,method_settings)

        rsp_settings = {'conv_thresh': 1.0e-6, 'b_frequencies': [-0.1],'c_frequencies': [0.3],'damping': 0.1,'a_component':'y','b_component':'y','c_component':'x'}

        qrf_prop.update_settings(rsp_settings, method_settings)

        qrf_result_yyx = qrf_prop.compute(molecule, ao_basis, scf_tensors,method_settings)
        
        # x-component 

        self.assertTrue(abs(qrf_result_xxx[(-0.1,0.3)].real - ref_result['xxx'].real) < 1.0e-5)

        self.assertTrue(abs(qrf_result_xxx[(-0.1,0.3)].imag - ref_result['xxx'].imag) < 1.0e-5)
        
        # y-component 
        
        self.assertTrue(abs(qrf_result_yyx[(-0.1,0.3)].real - ref_result['yyx'].real) < 1.0e-5)

        self.assertTrue(abs(qrf_result_yyx[(-0.1,0.3)].imag - ref_result['yyx'].imag) < 1.0e-5)

        # z-component 
        
        self.assertTrue(abs(qrf_result_zzx[(-0.1,0.3)].real - ref_result['zzx'].real) < 1.0e-5)

        self.assertTrue(abs(qrf_result_zzx[(-0.1,0.3)].imag - ref_result['zzx'].imag) < 1.0e-5)


    def test_qrf(self):

        ref_result = { 
            'xxx': 36.088147 + 29.283161j,
            'zzx': 32.539982 + 16.600985j,
            'yyx': 1.880527 + 5.129818j,
        }

        self.run_qrf(ref_result)


if __name__ == '__main__':
    unittest.main()
