from mpi4py import MPI
from pathlib import Path
import numpy as np
import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver

try:
    import pyframe
except ImportError:
    pass


@pytest.mark.solvers
class TestPolarizableEmbeddingGradient:

    @staticmethod
    def get_molecule_and_basis():

        xyz_string = """3
        xyz
        O    1.206898105612    1.227024188072   -0.113458607566
        H    0.490554304233    1.138628097468    0.508036765414
        H    1.992568601678    1.153431301285    0.406410355109
        """
        basis_label = 'def2-svp'
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, basis_label, ostream=None)
        return mol, bas

    @staticmethod
    def get_embedding_dict(options_file):

        return {
            'settings': {
                'embedding_method': 'PE',
                'induced_dipoles': {
                    'solver': 'jacobi',
                    'mic': False,
                    'threshold': 1e-8,
                    'max_iterations': 100,
                },
            },
            'inputs': {
                'json_file': options_file,
            },
        }

    def run_scf_grad_with_pe(self, flag):

        mol, bas = self.get_molecule_and_basis()

        here = Path(__file__).parent
        options_file = str(here / 'data' / 'pe_water.json')

        if flag == 'restricted':
            scf_drv = ScfRestrictedDriver()
        elif flag == 'unrestricted':
            mol.set_charge(1)
            mol.set_multiplicity(2)
            scf_drv = ScfUnrestrictedDriver()

        scf_drv.xcfun = 'b3lyp'
        scf_drv.embedding = self.get_embedding_dict(options_file)
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.compute(mol, bas, scf_results)

        return grad_drv.get_gradient()

    @pytest.mark.skipif('pyframe' not in sys.modules,
                        reason='pyframe not available')
    def test_rest_scf_grad_with_pe(self):

        grad = self.run_scf_grad_with_pe('restricted')

        ref_grad = np.array([
            [0.008611638727, -0.002315029768, 0.024688869615],
            [0.010049749264, 0.001770366431, -0.010871444189],
            [-0.017517055291, 0.001949040005, -0.013629960165],
        ])

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(np.abs(grad - ref_grad)) < 5.0e-4

    @pytest.mark.skipif('pyframe' not in sys.modules,
                        reason='pyframe not available')
    def test_unrest_scf_grad_with_pe(self):

        grad = self.run_scf_grad_with_pe('unrestricted')

        ref_grad = np.array([
            [0.010596876656, -0.008838191086, 0.061742850690],
            [0.054171891115, 0.010786228266, -0.032188724788],
            [-0.057122709691, 0.005932945868, -0.028860672998],
        ])

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(np.abs(grad - ref_grad)) < 5.0e-4
