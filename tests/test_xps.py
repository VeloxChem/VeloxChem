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
            assert record['contribution'] == pytest.approx(0.9998668939,
                                                           abs=1.0e-6)
            assert record['ionization_energy_ev'] == pytest.approx(539.249,
                                                                   abs=1.0e-2)
            assert record['intensity'] == pytest.approx(0.7807865496,
                                                        abs=1.0e-4)
            assert record['is_delocalized'] is False

    def test_water_oxygen_shake_ups(self):

        molecule, basis = self.get_water_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = 'pbe0'
        scf_drv.ostream.mute()
        scf_drv.compute(molecule, basis)

        xps_drv = XPSDriver()
        xps_drv.ostream.mute()
        xps_drv.shake_ups = True
        xps_drv.nstates = 5
        xps_results = xps_drv.compute(molecule, basis, scf_drv, element='O')

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert list(xps_results.keys()) == ['O']
            # mainline followed by the shake-up states in increasing energy
            assert len(xps_results['O']) == 6

            mainline = xps_results['O'][0]

            assert mainline['mo_index'] == 0
            assert mainline['atom_index'] == 0
            assert mainline['ionization_energy_ev'] == pytest.approx(539.249,
                                                                     abs=1.0e-2)
            assert mainline['intensity'] == pytest.approx(0.7807865496,
                                                          abs=1.0e-4)

            ref_shakeups = [
                (553.618, 0.0),
                (554.897, 0.0),
                (554.956, 2.422058e-4),
                (555.426, 0.0),
                (556.307, 4.8574839e-3),
            ]

            for record, (ref_ie, ref_intensity) in zip(xps_results['O'][1:],
                                                       ref_shakeups):
                assert record['mo_index'] == 0
                assert record['atom_index'] == 0
                assert record['ionization_energy_ev'] == pytest.approx(
                    ref_ie, abs=1.0e-2)
                assert record['intensity'] == pytest.approx(ref_intensity,
                                                            abs=2.0e-5)
                assert record['is_delocalized'] is False

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
            assert oxygen_record['contribution'] == pytest.approx(0.9993365235,
                                                                  abs=1.0e-6)
            assert oxygen_record['ionization_energy_ev'] == pytest.approx(
                540.410, abs=1.0e-2)
            assert oxygen_record['intensity'] == pytest.approx(0.7396281620,
                                                               abs=1.0e-4)
            assert oxygen_record['is_delocalized'] is False

            carbon_record = xps_results['C'][0]

            assert carbon_record['mo_index'] == 1
            assert carbon_record['atom_index'] == 0
            assert carbon_record['contribution'] == pytest.approx(0.9989394794,
                                                                  abs=1.0e-6)
            assert carbon_record['ionization_energy_ev'] == pytest.approx(
                293.989, abs=1.0e-2)
            assert carbon_record['intensity'] == pytest.approx(0.7452093207,
                                                               abs=1.0e-4)
            assert carbon_record['is_delocalized'] is False

    def test_methanol_shake_ups(self):

        molecule, basis = self.get_methanol_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.compute(molecule, basis)

        xps_drv = XPSDriver()
        xps_drv.ostream.mute()
        xps_drv.shake_ups = True
        xps_drv.nstates = 3
        xps_results = xps_drv.compute(molecule, basis, scf_drv, element=['O', 'C'])

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert list(xps_results.keys()) == ['O', 'C']
            # mainline followed by the shake-up states in increasing energy
            assert len(xps_results['O']) == 4
            assert len(xps_results['C']) == 4

            oxygen_mainline = xps_results['O'][0]

            assert oxygen_mainline['mo_index'] == 0
            assert oxygen_mainline['atom_index'] == 1
            assert oxygen_mainline['ionization_energy_ev'] == pytest.approx(
                540.410, abs=1.0e-2)
            assert oxygen_mainline['intensity'] == pytest.approx(0.7396281620,
                                                                 abs=1.0e-4)

            ref_oxygen_shakeups = [
                (550.117, 3.8184e-6),
                (550.355, 1.410530e-4),
                (551.418, 8.7338e-6),
            ]

            for record, (ref_ie, ref_intensity) in zip(xps_results['O'][1:],
                                                       ref_oxygen_shakeups):
                assert record['mo_index'] == 0
                assert record['atom_index'] == 1
                assert record['ionization_energy_ev'] == pytest.approx(
                    ref_ie, abs=1.0e-2)
                assert record['intensity'] == pytest.approx(ref_intensity,
                                                            abs=2.0e-5)
                assert record['is_delocalized'] is False

            carbon_mainline = xps_results['C'][0]

            assert carbon_mainline['mo_index'] == 1
            assert carbon_mainline['atom_index'] == 0
            assert carbon_mainline['ionization_energy_ev'] == pytest.approx(
                293.989, abs=1.0e-2)
            assert carbon_mainline['intensity'] == pytest.approx(0.7452093207,
                                                                 abs=1.0e-4)

            ref_carbon_shakeups = [
                (301.063, 2.203e-7),
                (303.119, 7.615e-7),
                (303.318, 3.7645967e-3),
            ]

            for record, (ref_ie, ref_intensity) in zip(xps_results['C'][1:],
                                                       ref_carbon_shakeups):
                assert record['mo_index'] == 1
                assert record['atom_index'] == 0
                assert record['ionization_energy_ev'] == pytest.approx(
                    ref_ie, abs=1.0e-2)
                assert record['intensity'] == pytest.approx(ref_intensity,
                                                            abs=2.0e-5)
                assert record['is_delocalized'] is False

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
                'intensity': 0.7807865496,
                'contribution': 0.9998668939,
                'is_delocalized': False,
            }],
            'C': [{
                'mo_index': 1,
                'atom_index': 2,
                'ionization_energy_ev': 293.989,
                'intensity': 0.7452093207,
                'contribution': 0.9989394794,
                'is_delocalized': False,
            }, {
                'mo_index': 1,
                'atom_index': 2,
                'ionization_energy_ev': 303.318,
                'intensity': 0.0037645967,
                'contribution': 0.9989394794,
                'is_delocalized': False,
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
                'intensity': 0.7807865496,
                'contribution': 0.9998668939,
                'is_delocalized': False,
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
