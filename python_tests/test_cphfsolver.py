import numpy as np
import unittest
import pytest
import sys
import h5py
from pathlib import Path

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.cphfsolver import CphfSolver

try:
    import pyscf
except ImportError:
    pass

class TestCphfSolver(unittest.TestCase):

    def run_cphfsolver(self, molecule, basis, xcfun, label):
        scf_drv = ScfRestrictedDriver()
        method_settings = {}
        if xcfun is not None:
            scf_drv.dft = True
            scf_drv.xcfun = xcfun
            method_settings = {'xcfun':xcfun}

        scf_drv.ostream.mute()
        scf_tensors = scf_drv.compute(molecule, basis)

        cphf_solver = CphfSolver(scfdrv=scf_drv)
        cphf_settings = {'conv_thresh':2e-7}
        cphf_solver.update_settings(cphf_settings, method_settings)
        cphf_solver.compute(molecule, basis, scf_tensors)

        if scf_drv.rank == mpi_master():
            cphf_coefficients = cphf_solver.cphf_results['cphf_ov']
        
            here = Path(__file__).parent
            hf_file_name = str(here /'inputs'/'cphf_coefficients.h5')
            hf = h5py.File(hf_file_name, 'r')
            print(hf.keys())
            cphf_reference = np.array(hf.get(label))
            hf.close()
             
            print("\ncphf_reference\n", cphf_reference)

            self.assertTrue(np.max(np.abs(cphf_coefficients - cphf_reference))
                             < 1.0e-6)
        
    @pytest.mark.skipif('pyscf' not in sys.modules,
                        reason='pyscf for integral derivatives not available')
    def test_cphf_coefficients(self):
        nh3_xyz = """4
        
        N     0.000000000     0.000000000     0.000000000
        H    -0.653401663     0.309213352     0.817609879
        H     0.695693936     0.071283622    -0.702632331
        H     0.330018952    -0.826347607     0.052270023
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(nh3_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        self.run_cphfsolver(molecule, basis, None, "cphf_coefficients")

    @pytest.mark.skipif('pyscf' not in sys.modules,
                        reason='pyscf for integral derivatives not available')
    def test_cpks_coefficients(self):
        nh3_xyz = """4
        
        N     0.000000000     0.000000000     0.000000000
        H    -0.653401663     0.309213352     0.817609879
        H     0.695693936     0.071283622    -0.702632331
        H     0.330018952    -0.826347607     0.052270023
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(nh3_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)
        nh3_xyz = """4
        
        N     0.000000000     0.000000000     0.000000000
        H    -0.653401663     0.309213352     0.817609879
        H     0.695693936     0.071283622    -0.702632331
        H     0.330018952    -0.826347607     0.052270023
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(nh3_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        self.run_cphfsolver(molecule, basis, "b3lyp", "cpks_coefficients_b3lyp")
