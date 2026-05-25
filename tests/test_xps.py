from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.xpsdriver import XPSDriver


@pytest.mark.solvers
class TestXpsDriver:

    @staticmethod
    def get_water_and_basis():

        xyz_string = """3
        xyz
        O      0.000000    0.000000    0.000000
        H      0.758602    0.000000    0.504284
        H     -0.758602    0.000000    0.504284"""

        molecule = Molecule.read_xyz_string(xyz_string)
        basis = MolecularBasis.read(molecule, 'def2-tzvp', ostream=None)

        return molecule, basis

    def test_water_oxygen_core_ionization_energy(self):

        molecule, basis = self.get_water_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = 'pbe0'
        scf_drv.ostream.mute()
        scf_drv.compute(molecule, basis)

        xps_drv = XPSDriver()
        xps_drv.ostream.mute()
        xps_results = xps_drv.compute(molecule, basis, scf_drv, element='O')

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert list(xps_results.keys()) == ['O']
            assert len(xps_results['O']) == 1

            mo_index, atom_index, ionization_energy, contribution = xps_results['O'][0]

            assert mo_index == 0
            assert atom_index == 0
            assert contribution == pytest.approx(0.6004496392, abs=1.0e-6)
            assert ionization_energy == pytest.approx(539.249, abs=1.0e-2)

    @staticmethod
    def get_methanol_and_basis():

        xyz_string = """6
        xyz
        C      0.000000    0.000000    0.000000
        O      1.420000    0.000000    0.000000
        H     -0.540000    0.940000    0.000000
        H     -0.540000   -0.470000    0.814000
        H     -0.540000   -0.470000   -0.814000
        H      1.760000    0.000000    0.940000"""

        molecule = Molecule.read_xyz_string(xyz_string)
        basis = MolecularBasis.read(molecule, 'def2-sv(p)', ostream=None)

        return molecule, basis

    def test_methanol_carbon_oxygen_core_ionization_energy(self):

        molecule, basis = self.get_methanol_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.compute(molecule, basis)

        xps_drv = XPSDriver()
        xps_drv.ostream.mute()
        xps_results = xps_drv.compute(molecule, basis, scf_drv, elements=['O', 'C'])

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert list(xps_results.keys()) == ['O', 'C']
            assert len(xps_results['O']) == 1
            assert len(xps_results['C']) == 1

            mo_index, atom_index, ionization_energy, contribution = xps_results['O'][0]

            assert mo_index == 0
            assert atom_index == 1
            assert contribution == pytest.approx(0.9824145671, abs=1.0e-6)
            assert ionization_energy == pytest.approx(540.410, abs=1.0e-2)

            mo_index, atom_index, ionization_energy, contribution = xps_results['C'][0]

            assert mo_index == 1
            assert atom_index == 0
            assert contribution == pytest.approx(0.9837630584, abs=1.0e-6)
            assert ionization_energy == pytest.approx(293.989, abs=1.0e-2)
