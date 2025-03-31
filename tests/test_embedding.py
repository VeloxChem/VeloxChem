from mpi4py import MPI
from pathlib import Path
import numpy as np
import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lrsolver import LinearResponseSolver

try:
    import pyframe
except ImportError:
    pass


@pytest.mark.solvers
class TestPolarizableEmbedding:

    @staticmethod
    def get_molecule_and_basis(name):

        if name == 'acrolein':
            xyz_string = """8
                C3OH4
                C             -0.145335   -0.546770    0.000607
                C              1.274009   -0.912471   -0.000167
                C              1.630116   -2.207690   -0.000132
                O             -0.560104    0.608977    0.000534
                H             -0.871904   -1.386459    0.001253
                H              2.004448   -0.101417   -0.000710
                H              0.879028   -3.000685    0.000484
                H              2.675323   -2.516779   -0.000673
            """
            basis_label = 'sto-3g'
            mol = Molecule.read_xyz_string(xyz_string)
            bas = MolecularBasis.read(mol, basis_label, ostream=None)
            return mol, bas
        else:
            return None, None

    def run_scf_with_pe(self, name):

        mol, bas = self.get_molecule_and_basis(name)

        here = Path(__file__).parent
        options_file = str(here / 'data' / f'{name}.json')

        scf_drv = ScfRestrictedDriver()
        scf_drv.embedding = {
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
        scf_drv.conv_thresh = 1.0e-8

        scf_drv.ostream.mute()

        return scf_drv.compute(mol, bas)

    def run_lrs_with_pe(self, scf_results, name, freqs):

        mol, bas = self.get_molecule_and_basis(name)

        here = Path(__file__).parent
        options_file = str(here / 'data' / f'{name}.json')

        lrsolver = LinearResponseSolver()
        lrsolver.embedding = {
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
        lrsolver.frequencies = freqs
        lrsolver.conv_thresh = 1.0e-8

        lrsolver.ostream.mute()

        return lrsolver.compute(mol, bas, scf_results)

    @pytest.mark.skipif('pyframe' not in sys.modules,
                        reason='pyframe not available')
    def test_scf_with_pe(self):

        scf_results = self.run_scf_with_pe(name='acrolein')

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert scf_results['scf_energy'] == pytest.approx(-188.314434428374,
                                                              rel=1e-10)

    @pytest.mark.skipif('pyframe' not in sys.modules,
                        reason='pyframe not available')
    def test_lrs_with_pe(self):

        scf_results = self.run_scf_with_pe(name='acrolein')

        ref_freqs = [0.0, 0.04]
        lrs_results = self.run_lrs_with_pe(name='acrolein',
                                           scf_results=scf_results,
                                           freqs=ref_freqs)

        reference_data = """
            2.097957814289E+01 -8.652608839280E+00 -7.686399537652E-03
            -8.652608839280E+00 3.353258603117E+01 -1.147335054308E-03
            -7.686399537652E-03 -1.147335054308E-03 4.586549080720E+00
            2.111014439306E+01 -8.811667762577E+00 -7.720196012987E-03
            -8.811667762577E+00 3.385544639017E+01 -1.118865578820E-03
            -7.720196012987E-03 -1.118865578820E-03 4.601670440902E+00
            """
        ref_prop = np.array([float(x) for x in reference_data.split()])

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            prop = np.array([
                -1.0 * lrs_results['response_functions'][(a, b, w)]
                for w in ref_freqs
                for a in 'xyz'
                for b in 'xyz'
            ])

            assert np.max(np.abs((prop - ref_prop) / prop)) < 5.0e-6
