from mpi4py import MPI
import pytest
import h5py

from veloxchem.veloxchemlib import mpi_master
from veloxchem.resultsio import read_results
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

            record = xps_results['O'][0]

            assert record['mo_index'] == 0
            assert record['atom_index'] == 0
            assert record['contribution'] == pytest.approx(0.6004496392,
                                                           abs=1.0e-6)
            assert record['ionization_energy_ev'] == pytest.approx(539.249,
                                                                   abs=1.0e-2)
            assert record['is_delocalized'] is True

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
        xps_results = xps_drv.compute(molecule, basis, scf_drv, element=['O', 'C'])

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert list(xps_results.keys()) == ['O', 'C']
            assert len(xps_results['O']) == 1
            assert len(xps_results['C']) == 1

            oxygen_record = xps_results['O'][0]

            assert oxygen_record['mo_index'] == 0
            assert oxygen_record['atom_index'] == 1
            assert oxygen_record['contribution'] == pytest.approx(0.9824145671,
                                                                  abs=1.0e-6)
            assert oxygen_record['ionization_energy_ev'] == pytest.approx(
                540.410, abs=1.0e-2)
            assert oxygen_record['is_delocalized'] is False

            carbon_record = xps_results['C'][0]

            assert carbon_record['mo_index'] == 1
            assert carbon_record['atom_index'] == 0
            assert carbon_record['contribution'] == pytest.approx(0.9837630584,
                                                                  abs=1.0e-6)
            assert carbon_record['ionization_energy_ev'] == pytest.approx(
                293.989, abs=1.0e-2)
            assert carbon_record['is_delocalized'] is False

    def test_xps_results_hdf5_roundtrip(self, tmp_path):

        if MPI.COMM_WORLD.Get_rank() != mpi_master():
            return

        h5file = tmp_path / 'xps_results.h5'
        with h5py.File(h5file, 'w'):
            pass

        xps_results = {
            'O': [{
                'mo_index': 0,
                'atom_index': 0,
                'ionization_energy_ev': 539.249,
                'contribution': 0.6004496392,
                'is_delocalized': True,
            }],
            'C': [{
                'mo_index': 1,
                'atom_index': 2,
                'ionization_energy_ev': 293.989,
                'contribution': 0.9837630584,
                'is_delocalized': False,
            }, {
                'mo_index': 3,
                'atom_index': 1,
                'ionization_energy_ev': 295.123,
                'contribution': 0.5123,
                'is_delocalized': True,
            }],
        }

        xps_drv = XPSDriver()
        xps_drv.ostream.mute()
        xps_drv._write_final_hdf5(str(h5file), xps_results)

        recovered = read_results(str(h5file), 'xps')

        assert recovered.keys() == xps_results.keys()
        for element in xps_results:
            assert recovered[element] == xps_results[element]

    def test_xps_results_hdf5_roundtrip_without_suffix(self, tmp_path):

        if MPI.COMM_WORLD.Get_rank() != mpi_master():
            return

        h5stem = tmp_path / 'xps_results'
        h5file = tmp_path / 'xps_results.h5'
        with h5py.File(h5file, 'w'):
            pass

        xps_results = {
            'O': [{
                'mo_index': 0,
                'atom_index': 0,
                'ionization_energy_ev': 539.249,
                'contribution': 0.6004496392,
                'is_delocalized': True,
            }],
        }

        xps_drv = XPSDriver()
        xps_drv.ostream.mute()
        xps_drv._write_final_hdf5(str(h5stem), xps_results)

        recovered = read_results(str(h5file), 'xps')

        assert recovered == xps_results

    def test_normalize_element_input(self):

        assert XPSDriver._normalize_element_input('O,C,F') == ['O', 'C', 'F']
        assert XPSDriver._normalize_element_input('O, C, F') == ['O', 'C', 'F']
        assert XPSDriver._normalize_element_input(['C', 'O', 'F']) == ['C', 'O', 'F']
