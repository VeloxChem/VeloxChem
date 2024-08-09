import tempfile
import numpy as np
import math as mt
import copy as cp
import pickle

from mpi4py import MPI
from pathlib import Path

from veloxchem import Point
from veloxchem import Molecule
from veloxchem import bohr_in_angstrom


class TestMolecule:

    def nh3_elements(self):
        return [7, 1, 1, 1]

    def nh3_coords(self):

        return [
            Point([-3.710, 3.019, -0.037]),
            Point([-3.702, 4.942, 0.059]),
            Point([-4.704, 2.415, 1.497]),
            Point([-4.780, 2.569, -1.573]),
        ]

    def nh3_xyzstr(self):

        return """N  -3.710   3.019  -0.037
                  H  -3.702   4.942   0.059
                  H  -4.704   2.415   1.497
                  H  -4.780   2.569  -1.573"""

    def h2o_xyzstr(self):

        return """O   0.000   0.000  -1.000
                  H   0.000   1.400  -2.100
                  H   0.000  -1.400  -2.100"""

    def nh3_h2o_xyzstr(self):

        return """N  -3.710   3.019  -0.037
                  H  -3.702   4.942   0.059
                  H  -4.704   2.415   1.497
                  H  -4.780   2.569  -1.573
                  O   0.000   0.000  -1.000
                  H   0.000   1.400  -2.100
                  H   0.000  -1.400  -2.100"""

    def ch4_xyzstr(self):

        return """C   0.0000000000    0.0000000000    0.0000000000
                  H   0.6298891440    0.6298891440    0.6298891440
                  H  -0.6298891440   -0.6298891440    0.6298891440
                  H   0.6298891440   -0.6298891440   -0.6298891440
                  H  -0.6298891440    0.6298891440   -0.6298891440"""

    def co2_xyzstr(self):

        return """C   0.0000000000    0.0000000000    0.0000000000
                  O   0.0000000000    0.0000000000    1.2000000000
                  O   0.0000000000    0.0000000000   -1.2000000000"""

    def nh3_molecule(self):

        return Molecule(self.nh3_elements(), self.nh3_coords(), "au")

    def test_constructors(self):

        tol = 1.0e-12

        elements = self.nh3_elements()
        coords = self.nh3_coords()
        xyzstr = self.nh3_xyzstr()

        # test static read_str method
        mol_a = Molecule(elements, coords, 'au')
        mol_b = Molecule.read_str(xyzstr, 'au')
        assert mol_a == mol_b

        # test unit changes in constructors
        f = bohr_in_angstrom()
        for r in coords:
            r.scale(f)
        mol_b = Molecule(elements, coords, 'angstrom')
        assert mol_a == mol_b

        # test copy constructor
        mol_b = Molecule(mol_a)
        assert mol_a == mol_b

        # test empty constructor
        mol_b = Molecule()
        assert mol_b.get_multiplicity() == 1
        assert mt.isclose(mol_b.get_charge(), 0.0, rel_tol=tol, abs_tol=tol)
        assert mol_b.number_of_atoms() == 0

        # test composition constructor
        mol_b = Molecule.read_str(self.h2o_xyzstr(), 'au')
        mol_c = Molecule(mol_a, mol_b)
        mol_d = Molecule.read_str(self.nh3_h2o_xyzstr(), 'au')
        assert mol_c == mol_d

    def test_pickle(self):

        mol_a = self.nh3_molecule()
        bobj = pickle.dumps(mol_a)
        mol_b = pickle.loads(bobj)
        assert mol_a == mol_b

    def test_add_atom(self):

        mol_a = self.nh3_molecule()

        # add_atom
        mol_b = Molecule()
        mol_b.add_atom(7, Point([-3.710, 3.019, -0.037]), 'au')
        mol_b.add_atom(1, Point([-3.702, 4.942, 0.059]), 'au')
        mol_b.add_atom(1, Point([-4.704, 2.415, 1.497]), 'au')
        mol_b.add_atom(1, Point([-4.780, 2.569, -1.573]), 'au')
        assert mol_a == mol_b

        # test unit changes in add_atom
        f = bohr_in_angstrom()
        coords = self.nh3_coords()
        for r in coords:
            r.scale(f)
        mol_b = Molecule()
        mol_b.add_atom(7, coords[0], 'angstrom')
        mol_b.add_atom(1, coords[1], 'angstrom')
        mol_b.add_atom(1, coords[2], 'angstrom')
        mol_b.add_atom(1, coords[3], 'angstrom')
        assert mol_a == mol_b

    def test_slice(self):

        mol_a = Molecule.read_str(self.nh3_h2o_xyzstr(), 'au')

        # slice NH3 fragment
        mol_b = mol_a.slice([0, 1, 2, 3])
        mol_c = Molecule.read_str(self.nh3_xyzstr(), 'au')
        assert mol_b == mol_c

        # slice H2O fragment
        mol_b = mol_a.slice([4, 5, 6])
        mol_c = Molecule.read_str(self.h2o_xyzstr(), 'au')
        assert mol_b == mol_c

    def test_charge(self):

        tol = 1.0e-12

        mol = self.nh3_molecule()
        q = mol.get_charge()
        assert mt.isclose(q, 0.0, rel_tol=tol, abs_tol=tol)
        mol.set_charge(1.0)
        q = mol.get_charge()
        assert mt.isclose(q, 1.0, rel_tol=tol, abs_tol=tol)

    def test_multiplicity(self):

        mol = self.nh3_molecule()
        assert mol.get_multiplicity() == 1
        mol.set_multiplicity(2)
        assert mol.get_multiplicity() == 2

    def test_number_of_atoms(self):

        mol = self.nh3_molecule()
        assert mol.number_of_atoms() == 4
        assert mol.number_of_atoms(1) == 3
        assert mol.number_of_atoms(7) == 1
        assert mol.number_of_atoms(0, 1, 1) == 0
        assert mol.number_of_atoms(0, 2, 1) == 1

    def test_elemental_composition(self):

        mol = self.nh3_molecule()
        assert mol.get_elemental_composition() == set((7, 1))

    def test_number_of_electrons(self):

        tol = 1.0e-12

        mol = self.nh3_molecule()
        n_e = mol.number_of_electrons()
        assert mt.isclose(n_e, 10.0, rel_tol=tol, abs_tol=tol)

    def test_identifiers(self):

        mol = self.nh3_molecule()
        assert mol.get_identifiers() == [7, 1, 1, 1]

    def test_coordinates(self):

        mol = self.nh3_molecule()

        coords = self.nh3_coords()

        coords_a = mol.get_coordinates()
        for r_a, r_b in zip(coords, coords_a):
            assert r_a == r_b

        coords_a = mol.get_coordinates('au')
        for r_a, r_b in zip(coords, coords_a):
            assert r_a == r_b

        f = bohr_in_angstrom()
        for r in coords:
            r.scale(f)
        coords_a = mol.get_coordinates('angstrom')
        for r_a, r_b in zip(coords, coords_a):
            assert r_a == r_b

    def test_charges(self):

        tol = 1.0e-12

        mol = self.nh3_molecule()
        charges = np.array(mol.get_charges())
        assert np.allclose(charges, np.array([7.0, 1.0, 1.0, 1.0]), tol, tol,
                           False)

    def test_masses(self):

        tol = 1.0e-12

        mol = self.nh3_molecule()
        masses = np.array(mol.get_masses())
        assert np.allclose(masses,
                           np.array([14.003074, 1.007825, 1.007825, 1.007825]),
                           tol, tol, False)

    def test_labels(self):

        mol = self.nh3_molecule()
        assert mol.get_labels() == ['N', 'H', 'H', 'H']

    def test_label(self):

        mol = self.nh3_molecule()
        assert mol.get_label(0) == 'N'
        assert mol.get_label(1) == 'H'
        assert mol.get_label(2) == 'H'
        assert mol.get_label(3) == 'H'

    def test_atom_coordinates(self):

        tol = 1.0e-12

        mol = self.nh3_molecule()

        # check nitrogen atom coordinates in au
        coords = mol.get_atom_coordinates(0, 'au')
        assert coords == Point([-3.710, 3.019, -0.037])

        # check last hydrogen atom coordinates in au
        coords = mol.get_atom_coordinates(3, 'au')
        assert coords == Point([-4.780, 2.569, -1.573])

        # check nitrogen atom coordinates in angstrom
        f = bohr_in_angstrom()
        coords = mol.get_atom_coordinates(0, 'angstrom')
        assert coords == Point([-3.710 * f, 3.019 * f, -0.037 * f])

        # check last hydrogen atom coordinates in au
        coords = mol.get_atom_coordinates(3, 'angstrom')
        assert coords == Point([-4.780 * f, 2.569 * f, -1.573 * f])

    def test_atom_indices(self):

        mol = self.nh3_molecule()
        assert mol.atom_indices('B') == []
        assert mol.atom_indices('N') == [
            0,
        ]
        assert mol.atom_indices('H') == [
            1,
            2,
            3,
        ]

    def test_nuclear_repulsion_energy(self):

        tol = 1.0e-12

        mol = Molecule.read_str(self.h2o_xyzstr(), 'au')
        assert mt.isclose(mol.nuclear_repulsion_energy(),
                          9.34363815797054450919,
                          rel_tol=tol,
                          abs_tol=tol)

    def test_check_proximity(self):

        mol = self.nh3_molecule()
        assert mol.check_proximity(0.1)
        assert not mol.check_proximity(5.1)

    def test_get_string(self):

        mol = self.nh3_molecule()
        molstr = mol.get_string()
        lines = molstr.splitlines()
        assert lines[0] == 'Molecular Geometry (Angstroms)'
        assert lines[1] == '================================'
        assert lines[2] == ''
        assert lines[
            3] == '  Atom         Coordinate X          Coordinate Y          Coordinate Z  '
        assert lines[4] == ''
        assert lines[
            5] == '  N          -1.963247452450        1.597585999716       -0.019579556803'
        assert lines[
            6] == '  H          -1.959014034763        2.615193776283        0.031221455443'
        assert lines[
            7] == '  H          -2.489249600088        1.277962964331        0.792178284722'
        assert lines[
            8] == '  H          -2.529467068116        1.359456254810       -0.832395752750'
        assert lines[9] == ''

    def test_write_xyz(self):

        with tempfile.TemporaryDirectory() as temp_dir:
            if MPI.COMM_WORLD.Get_rank() == 0:
                fname = str(Path(temp_dir, 'mol.xyz'))
                mol_a = self.nh3_molecule()
                mol_a.write_xyz(fname)
                mol_b = Molecule.read_xyz(fname)
                assert mol_a == mol_b

    def test_read_dict(self):

        rxyz = [
            'N  -3.710   3.019  -0.037', 'H  -3.702   4.942   0.059',
            'H  -4.704   2.415   1.497', 'H  -4.780   2.569  -1.573'
        ]
        mdict = {"xyz": rxyz, "charge": 3.0, "multiplicity": 2, "units": "au"}

        mol_a = Molecule.from_dict(mdict)
        mol_b = self.nh3_molecule()
        mol_b.set_charge(3.0)
        mol_b.set_multiplicity(2)
        assert mol_a == mol_b

    def test_moments_of_inertia(self):

        tol = 1.0e-12

        mol = Molecule.read_str(self.ch4_xyzstr(), 'au')
        imoms = mol.moments_of_inertia()
        rmoms = np.array(
            [3.198919866723860, 3.198919866723860, 3.198919866723860])
        assert np.allclose(rmoms, imoms, tol, tol, False)

    def test_is_linear(self):

        mol = Molecule.read_str(self.ch4_xyzstr(), 'au')
        assert not mol.is_linear()
        mol = Molecule.read_str(self.co2_xyzstr(), 'au')
        assert mol.is_linear()

    def test_center_of_mass_bohr(self):

        tol = 1.0e-12

        mol = self.nh3_molecule()
        mol_com = mol.center_of_mass_in_bohr()
        ref_com = np.array(
            [-3.831697485497502, 3.0704373126932536, -0.0314360099042971])
        assert np.allclose(ref_com, mol_com, tol, tol, False)

    def test_center_of_mass_angstrom(self):

        tol = 1.0e-12

        mol = self.nh3_molecule()
        mol_com = mol.center_of_mass_in_angstrom()
        f = bohr_in_angstrom()
        ref_com = np.array([
            -3.831697485497502 * f, 3.0704373126932536 * f,
            -0.0314360099042971 * f
        ])
        assert np.allclose(ref_com, mol_com, tol, tol, False)

    def test_more_info(self):

        mol = self.nh3_molecule()
        mol.set_multiplicity(3)
        molstr = mol.more_info()
        lines = molstr.splitlines()
        assert lines[0] == ('Molecular charge            : 0' + 39 * ' ')
        assert lines[1] == ('Spin multiplicity           : 3' + 39 * ' ')
        assert lines[2] == ('Number of atoms             : 4' + 39 * ' ')
        assert lines[3] == ('Number of alpha electrons   : 6' + 39 * ' ')
        assert lines[4] == ('Number of beta  electrons   : 4' + 39 * ' ')

    def test_number_of_alpha_electrons(self):

        mol = self.nh3_molecule()
        assert mol.number_of_alpha_electrons() == 5
        mol.set_multiplicity(3)
        assert mol.number_of_alpha_electrons() == 6

    def test_number_of_beta_electrons(self):

        mol = self.nh3_molecule()
        assert mol.number_of_beta_electrons() == 5
        mol.set_multiplicity(3)
        assert mol.number_of_beta_electrons() == 4

    def test_check_multiplicity(self):

        # singlet, doublet, and triplet NH3
        mol = self.nh3_molecule()
        assert mol.check_multiplicity()
        mol.set_multiplicity(2)
        assert not mol.check_multiplicity()
        mol.set_multiplicity(3)
        assert mol.check_multiplicity()

        # singlet, doublet, and triplet NH3+
        mol.set_charge(1.0)
        mol.set_multiplicity(1)
        assert not mol.check_multiplicity()
        mol.set_multiplicity(2)
        assert mol.check_multiplicity()
        mol.set_multiplicity(3)
        assert not mol.check_multiplicity()

    def test_get_aufbau_alpha_occupation(self):

        tol = 1.0e-12

        mol = Molecule.read_str(self.nh3_xyzstr())
        mol.set_charge(1.0)
        mol.set_multiplicity(2)

        nocc = mol.get_aufbau_alpha_occupation(6)
        rocc = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 0.0])
        assert np.allclose(rocc, nocc, tol, tol, False)

    def test_get_aufbau_beta_occupation(self):

        tol = 1.0e-12

        mol = Molecule.read_str(self.nh3_xyzstr())
        mol.set_charge(1.0)
        mol.set_multiplicity(2)

        nocc = mol.get_aufbau_beta_occupation(6)
        rocc = np.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.0])
        assert np.allclose(rocc, nocc, tol, tol, False)

    def test_get_aufbau_occupation(self):

        tol = 1.0e-12

        mol = Molecule.read_str(self.nh3_xyzstr())
        mol.set_charge(1.0)
        mol.set_multiplicity(2)

        nocc = mol.get_aufbau_occupation(6)
        rocc = np.array([1.0, 1.0, 1.0, 1.0, 0.5, 0.0])
        assert np.allclose(rocc, nocc, tol, tol, False)

        nocca, noccb = mol.get_aufbau_occupation(6, 'unrestricted')
        rocca = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 0.0])
        roccb = np.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.0])
        assert np.allclose(rocca, nocca, tol, tol, False)
        assert np.allclose(roccb, noccb, tol, tol, False)

    def test_coordinates_in_bohr(self):

        tol = 1.0e-12

        mol = self.nh3_molecule()
        coords_a = mol.get_coordinates_in_bohr()
        coords_b = np.array([[-3.710, 3.019, -0.037], [-3.702, 4.942, 0.059],
                             [-4.704, 2.415, 1.497], [-4.780, 2.569, -1.573]])
        assert np.allclose(coords_a, coords_b, tol, tol, False)

    def test_coordinates_in_angstrom(self):

        tol = 1.0e-12

        mol = self.nh3_molecule()
        coords_a = mol.get_coordinates_in_angstrom()
        coords_b = np.array([[-3.710, 3.019, -0.037], [-3.702, 4.942, 0.059],
                             [-4.704, 2.415, 1.497], [-4.780, 2.569, -1.573]])
        coords_b = coords_b * bohr_in_angstrom()
        assert np.allclose(coords_a, coords_b, tol, tol, False)

    def test_mpi_bcast(self):

        comm = MPI.COMM_WORLD

        mol_a = None
        if comm.Get_rank() == 0:
            mol_a = self.nh3_molecule()
        mol_a = comm.bcast(mol_a)
        mol_b = self.nh3_molecule()
        assert mol_a == mol_b
