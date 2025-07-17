from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.thgdriver import ThgDriver
from veloxchem.thgreddriver import ThgRedDriver
from veloxchem.outputstream import OutputStream


@pytest.mark.timeconsuming
class Testhgreduced:

    def run_thg_red(self, xcfun_label):

        molecule_string = """
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
        """
        basis_set_label = 'def2-svp'

        molecule = Molecule.read_molecule_string(molecule_string, units='au')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_settings = {}
        method_settings = {'xcfun': xcfun_label, 'grid_level': 1}

        scfdrv = ScfRestrictedDriver()
        scfdrv.ostream.mute()
        scfdrv.update_settings(scf_settings, method_settings)
        scf_results = scfdrv.compute(molecule, basis)

        rsp_settings = {
            'frequencies': [0.05],
            'damping': 0,
        }

        thg = ThgDriver(MPI.COMM_WORLD,OutputStream())
        thg.ostream.mute()
        thg.update_settings(rsp_settings, method_settings)
        thg_results = thg.compute(molecule, basis, scf_results)


        thgred = ThgRedDriver(MPI.COMM_WORLD,OutputStream())
        thgred.ostream.mute()
        thgred.update_settings(rsp_settings, method_settings)
        thgred_results = thgred.compute(molecule, basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():

            tol = 1.0e-5
            ref_val = thg_results["THG"][(0.05,0.05,0.05)]
            calc_val = thgred_results["THG"][(0.05,0.05,0.05)]
            assert abs(abs(calc_val / ref_val) - 1.0) < tol

    def test_thg_red_lda(self):

        self.run_thg_red('slda')

    def test_thg_red_gga(self):

        self.run_thg_red('pbe0')

    def test_thg_red_mgga(self):

        self.run_thg_red('tpssh')
