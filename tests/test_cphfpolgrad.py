from pathlib import Path
import numpy as np
import pytest
import h5py
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
#from veloxchem.cphfsolver import CphfSolver
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.cppsolver import ComplexResponse
from veloxchem.polorbitalresponse import PolOrbitalResponse

try:
    import pyscf
except ImportError:
    pass

class TestCphfPolgrad:

    def run_cphfpolgrad_real(self, molecule, basis, xcfun=None, label=None):
        scf_drv = ScfRestrictedDriver()
        scf_dict = {}
        method_settings = {}
        if xcfun is not None:
            scf_drv._dft = True
            scf_drv.xcfun = xcfun
            method_settings = {'xcfun': xcfun}

        scf_drv.update_settings(scf_dict, method_settings)
        scf_drv.ostream.mute()
        scf_tensors = scf_drv.compute(molecule, basis)

        # linear response
        rsp_settings = {'conv_thresh': 1.0e-5, 'frequencies': (0.0, 0.4)}
        lr_drv = LinearResponseSolver()
        lr_drv.a_operator = "electric dipole"
        lr_drv.b_operator = "electric dipole"
        lr_drv.update_settings(rsp_settings, method_settings)
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(molecule, basis, scf_tensors)

        # test real analytical gradient
        cphfpolgrad_solver = PolOrbitalResponse()
        cphfpolgrad_settings = {'conv_thresh':2e-7, 'frequencies': (0.0, 0.4),
                                'use_subspace_solver': 'no'}
        cphfpolgrad_solver.update_settings(cphfpolgrad_settings, method_settings)
        cphfpolgrad_solver.ostream.mute()
        cphfpolgrad_solver.compute(molecule, basis, scf_tensors, lr_results)

        if scf_drv.rank == mpi_master():
            cphfpolgrad_results = cphfpolgrad_solver.cphf_results
            cphfpolgrad_coefficients = cphfpolgrad_results['cphf_ov']
            np.set_printoptions(suppress=True, precision=10)
            here = Path(__file__).parent
            hf_file_name = str(here /'data'/'cphfpolgrad_coefficients.h5')
            hf = h5py.File(hf_file_name, 'r')
            cphfpolgrad_reference = np.array(hf.get(label))
            hf.close()

            # Here we are comparing the CPHF coefficients in MO basis, so
            # there might be sign differences; we compare absolute values instead.
            assert np.max(np.abs(cphfpolgrad_coefficients) - np.abs(cphfpolgrad_reference)) < 1.0e-6

    def run_cphfpolgrad_complex(self, molecule, basis, xcfun=None, label=None):
        scf_drv = ScfRestrictedDriver()
        scf_dict = {}
        method_settings = {}
        if xcfun is not None:
            scf_drv._dft = True
            scf_drv.xcfun = xcfun
            method_settings = {'xcfun': xcfun}

        scf_drv.update_settings(scf_dict, method_settings)
        scf_drv.ostream.mute()
        scf_tensors = scf_drv.compute(molecule, basis)

        # linear response
        rsp_settings = {'conv_thresh': 1.0e-5, 'frequencies': (0.0, 0.4),
                        'damping': 0.5}
        lr_drv = ComplexResponse()
        lr_drv.a_operator = "electric dipole"
        lr_drv.b_operator = "electric dipole"
        lr_drv.update_settings(rsp_settings, method_settings)
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(molecule, basis, scf_tensors)

        # test complex analytical gradient
        cphfpolgrad_solver = PolOrbitalResponse()
        cphfpolgrad_settings = {'conv_thresh':2e-7, 'frequencies': (0.0, 0.4),
                                'is_complex': 'yes', 'damping': 0.5,
                                'use_subspace_solver': 'no'}
        cphfpolgrad_solver.update_settings(cphfpolgrad_settings, method_settings)
        cphfpolgrad_solver.ostream.mute()
        cphfpolgrad_solver.compute(molecule, basis, scf_tensors, lr_results)

        if scf_drv.rank == mpi_master():
            cphfpolgrad_results = cphfpolgrad_solver.cphf_results
            cphfpolgrad_coefficients = cphfpolgrad_results['cphf_ov']
            np.set_printoptions(suppress=True, precision=10)
            here = Path(__file__).parent
            hf_file_name = str(here /'data'/'cphfpolgrad_coefficients.h5')
            hf = h5py.File(hf_file_name, 'r')
            cphfpolgrad_reference = np.array(hf.get(label))
            hf.close()

            # Here we ar comparing the CPHF coefficients in MO basis, so
            # there might be sign differences; we compare absolute values instead.
            assert np.max(np.abs(cphfpolgrad_coefficients) - np.abs(cphfpolgrad_reference)) < 1.0e-6

    @pytest.mark.skipif('pyscf' not in sys.modules,
                        reason='pyscf for integral derivatives not available')
    def test_cphfpolgrad_coefficients(self):
        h2o_xyz = """3

        O     0.000000    0.000000    0.000000
        H     0.000000    0.504284    0.758602
        H     0.000000   -0.504284    0.758602 
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(h2o_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        #self.run_cphfpolgrad_real(molecule, basis, None, "cphfpolgrad_coefficients_real")
        self.run_cphfpolgrad_real(molecule, basis, None, "cphfpolgrad_coefficients_red_real")
        #self.run_cphfpolgrad_complex(molecule, basis, None, "cphfpolgrad_coefficients_complex")
        self.run_cphfpolgrad_complex(molecule, basis, None, "cphfpolgrad_coefficients_red_complex")

    @pytest.mark.skipif('pyscf' not in sys.modules,
                        reason='pyscf for integral derivatives not available')
    def test_cpkspolgrad_coefficients(self):
        h2o_xyz = """3

        O     0.000000    0.000000    0.000000
        H     0.000000    0.504284    0.758602
        H     0.000000   -0.504284    0.758602 
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(h2o_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        #self.run_cphfpolgrad_real(molecule, basis, "b3lyp", "cpkspolgrad_coefficients_b3lyp_real")
        self.run_cphfpolgrad_real(molecule, basis, "b3lyp", "cpkspolgrad_coefficients_red_b3lyp_real")
        #self.run_cphfpolgrad_complex(molecule, basis, "b3lyp", "cpkspolgrad_coefficients_b3lyp_complex")
        self.run_cphfpolgrad_complex(molecule, basis, "b3lyp", "cpkspolgrad_coefficients_red_b3lyp_complex")
