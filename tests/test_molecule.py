from mpi4py import MPI
from pathlib import Path
from textwrap import dedent
import re
import numpy as np
import pickle
import pytest
import math
import sys
from unittest.mock import Mock

from veloxchem.veloxchemlib import Point
from veloxchem.veloxchemlib import bohr_in_angstrom, mpi_master
from veloxchem.veloxchemlib import get_vdw_radii_data_in_bohr
from veloxchem.veloxchemlib import chemical_element_identifier
from veloxchem.mpitask import MpiTask
from veloxchem.molecule import Molecule
from veloxchem.dispersionmodel import DispersionModel
from veloxchem.optimizationdriver import OptimizationDriver
from veloxchem.inputparser import get_random_string_serial

try:
    import rdkit
except ImportError:
    pass


class FakeBasis:

    def __init__(self, core_electrons):
        self._core_electrons = core_electrons

    def get_number_of_ecp_core_electrons(self):
        return self._core_electrons


class TestMolecule:

    def nh3_elements(self):
        return [7, 1, 1, 1]

    def nh3_labels(self):

        return ['N', 'H', 'H', 'H']

    def nh3_coords(self):

        return [
            Point([-3.710, 3.019, -0.037]),
            Point([-3.702, 4.942, 0.059]),
            Point([-4.704, 2.415, 1.497]),
            Point([-4.780, 2.569, -1.573]),
        ]

    def nh3_coords_np(self):

        return np.array([
            [-3.710, 3.019, -0.037],
            [-3.702, 4.942, 0.059],
            [-4.704, 2.415, 1.497],
            [-4.780, 2.569, -1.573],
        ])

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

    def cyclopropane_xyzstr(self):

        return """C              0.401049000000        -0.048676000000         0.720786000000
                  C              0.132408000000         0.924865000000         1.852698000000
                  C             -1.017462000000         0.171567000000         1.211199000000
                  H              0.853020000000        -1.011283000000         1.041107000000
                  H              0.629526000000         0.432541000000        -0.253752000000
                  H              0.198379000000         1.994972000000         1.562844000000
                  H              0.421836000000         0.551154000000         2.857718000000
                  H             -1.647101000000         0.786020000000         0.533335000000
                  H             -1.423602000000        -0.657837000000         1.828158000000"""

    def cyclobutane_xyzstr(self):

        return """C              1.191613000000        -0.128056000000        -0.605997000000
                  C             -0.183673000000        -0.778595000000        -0.609576000000
                  C             -0.472210000000        -0.161128000000        -1.969762000000
                  C              0.903078000000         0.489408000000        -1.966185000000
                  H              2.016724000000        -0.872771000000        -0.641091000000
                  H              1.301225000000         0.635243000000         0.195251000000
                  H             -0.846689000000        -0.380768000000         0.189663000000
                  H             -0.131187000000        -1.888779000000        -0.646682000000
                  H             -1.297316000000         0.583592000000        -1.934662000000
                  H             -0.581828000000        -0.924424000000        -2.771011000000
                  H              1.566090000000         0.091578000000        -2.765425000000
                  H              0.850600000000         1.599591000000        -1.929080000000"""

    def nh3_molecule(self):

        return Molecule(self.nh3_elements(), self.nh3_coords(), "au")

    def nh3_molecule_with_ghost_atom(self):

        nh3_str = self.nh3_xyzstr()
        nh3_str = nh3_str.replace('N', 'Bq_N')

        return Molecule.read_molecule_string(nh3_str, 'au')

    def test_constructors(self):

        tol = 1.0e-12

        elements = self.nh3_elements()
        labels = self.nh3_labels()
        coords = self.nh3_coords()
        coords_np = self.nh3_coords_np()
        xyzstr = self.nh3_xyzstr()

        # test static read_str method
        mol_a = Molecule(elements, coords, 'au')
        mol_b = Molecule(labels, coords_np, units='au')
        mol_c = Molecule.read_str(xyzstr, 'au')
        assert mol_a == mol_b
        assert mol_a == mol_c

        # test unit changes in constructors
        for r in coords:
            r.scale(bohr_in_angstrom())
        mol_b = Molecule(elements, coords, 'angstrom')
        assert mol_a == mol_b

        # test copy constructor
        mol_b = Molecule(mol_a)
        assert mol_a == mol_b

        # test empty constructor
        mol_b = Molecule()
        assert mol_b.get_multiplicity() == 1
        assert math.isclose(mol_b.get_charge(), 0.0, rel_tol=tol, abs_tol=tol)
        assert mol_b.number_of_atoms() == 0

        # test composition constructor
        mol_b = Molecule.read_str(self.h2o_xyzstr(), 'au')
        mol_c = Molecule(mol_a, mol_b)
        mol_d = Molecule.read_str(self.nh3_h2o_xyzstr(), 'au')
        assert mol_c == mol_d

    def test_read_molecule_string_ghost_atoms_and_basis_labels(self):

        mol = Molecule.read_molecule_string("""Bq_N  0.0  0.0  0.0  def2-svp
                                               H_Bq  0.0  0.0  1.0  aug-cc-pvdz
                                               O     0.0  1.0  0.0""")

        assert mol.get_labels() == ['Bq', 'Bq', 'O']
        assert mol.get_identifiers() == [0, 0, 8]
        assert mol.get_atom_basis_labels() == [('DEF2-SVP', 'N'),
                                               ('AUG-CC-PVDZ', 'H'), ('', 'O')]
        assert 'Bq_N' in mol.get_xyz_string()
        assert 'Bq_H' in mol.get_xyz_string()

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_read_xyz_string_validation(self):

        with pytest.raises(AssertionError,
                           match='Invalid number of atoms in XYZ input'):
            Molecule.read_xyz_string("""not-an-int
                                        comment
                                        H 0.0 0.0 0.0""")

        with pytest.raises(AssertionError,
                           match='Inconsistent number of atoms in XYZ input'):
            Molecule.read_xyz_string("""2
                                        comment
                                        H 0.0 0.0 0.0""")

    def test_read_xyz_file_success(self, tmp_path):

        xyz_path = tmp_path / 'water.xyz'
        xyz_path.write_text(
            dedent("""
                3
                water
                O  0.0  0.0  -1.0
                H  0.0  1.4  -2.1
                H  0.0 -1.4  -2.1
                """).lstrip())

        mol = Molecule.read_xyz_file(xyz_path)
        ref_mol = Molecule.read_molecule_string(self.h2o_xyzstr(), 'angstrom')
        assert mol == ref_mol

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_read_xyz_file_missing_file(self, tmp_path):

        # Note: use re.escape to avoid "incomplete escape \U" error on Windows
        missing_path = tmp_path / 'missing.xyz'
        with pytest.raises(
                AssertionError,
                match=re.escape(f'xyzfile {missing_path} does not exist')):
            Molecule.read_xyz_file(missing_path)

    def test_read_gro_file_element_guessing(self, tmp_path):

        gro_path = tmp_path / 'test.gro'
        gro_lines = [
            'test system',
            '3',
            f'{1:5d}{"ALA":<5}{"CA":>5}{1:5d}{0.100:8.3f}{0.200:8.3f}{0.300:8.3f}',
            f'{1:5d}{"LIG":<5}{"CL1":>5}{2:5d}{0.400:8.3f}{0.500:8.3f}{0.600:8.3f}',
            f'{1:5d}{"NA":<5}{"NA":>5}{3:5d}{0.700:8.3f}{0.800:8.3f}{0.900:8.3f}',
            '   1.00000   1.00000   1.00000',
        ]
        gro_path.write_text('\n'.join(gro_lines) + '\n')

        mol = Molecule.read_gro_file(gro_path)

        assert mol.get_labels() == ['C', 'Cl', 'Na']
        assert np.allclose(
            mol.get_coordinates_in_angstrom(),
            np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]))

    def test_read_pdb_file_element_guessing_and_explicit_columns(
            self, tmp_path):

        pdb_path = tmp_path / 'test.pdb'

        def pdb_record(record,
                       serial,
                       atom_name,
                       residue_name,
                       x,
                       y,
                       z,
                       element=''):
            return (f'{record:<6s}{serial:5d} '
                    f'{atom_name:^4s}'
                    f' '
                    f'{residue_name:>3s} A'
                    f'{1:4d}'
                    f'    '
                    f'{x:8.3f}{y:8.3f}{z:8.3f}'
                    f'{1.00:6.2f}{0.00:6.2f}'
                    f'          '
                    f'{element:>2s}')

        pdb_lines = [
            pdb_record('ATOM', 1, 'NA', 'NA', 1.0, 2.0, 3.0),
            pdb_record('HETATM', 2, 'CL1', 'UNL', 4.0, 5.0, 6.0),
            pdb_record('HETATM', 3, 'O', 'HOH', 7.0, 8.0, 9.0, 'O'),
        ]
        pdb_path.write_text('\n'.join(pdb_lines) + '\n')

        mol = Molecule.read_pdb_file(pdb_path)

        assert mol.get_labels() == ['Na', 'Cl', 'O']
        assert np.allclose(
            mol.get_coordinates_in_angstrom(),
            np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]))

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
        assert math.isclose(q, 0.0, rel_tol=tol, abs_tol=tol)
        mol.set_charge(1.0)
        q = mol.get_charge()
        assert math.isclose(q, 1.0, rel_tol=tol, abs_tol=tol)

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
        assert math.isclose(n_e, 10.0, rel_tol=tol, abs_tol=tol)

    def test_identifiers(self):

        mol = self.nh3_molecule()
        assert mol.get_identifiers() == [7, 1, 1, 1]

    def test_coordinates(self):

        mol = self.nh3_molecule()

        coords = self.nh3_coords()

        coords_a = mol._get_coordinates()
        for r_a, r_b in zip(coords, coords_a):
            assert r_a == r_b

        coords_a = mol._get_coordinates('au')
        for r_a, r_b in zip(coords, coords_a):
            assert r_a == r_b

        f = bohr_in_angstrom()
        for r in coords:
            r.scale(f)
        coords_a = mol._get_coordinates('angstrom')
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

        mol = self.nh3_molecule()

        # check nitrogen atom coordinates in au
        coords = mol.get_atom_coordinates(0, 'au')
        assert coords == [-3.710, 3.019, -0.037]

        # check last hydrogen atom coordinates in au
        coords = mol.get_atom_coordinates(3, 'au')
        assert coords == [-4.780, 2.569, -1.573]

        # check nitrogen atom coordinates in angstrom
        f = bohr_in_angstrom()
        coords = mol.get_atom_coordinates(0, 'angstrom')
        assert coords == [-3.710 * f, 3.019 * f, -0.037 * f]

        # check last hydrogen atom coordinates in au
        coords = mol.get_atom_coordinates(3, 'angstrom')
        assert coords == [-4.780 * f, 2.569 * f, -1.573 * f]

    def test_atom_indices(self):

        mol = self.nh3_molecule()
        assert mol.atom_indices('B') == []
        assert mol.atom_indices('N') == [0]
        assert mol.atom_indices('H') == [1, 2, 3]

    def test_nuclear_repulsion_energy(self):

        tol = 1.0e-12

        mol = Molecule.read_str(self.h2o_xyzstr(), 'au')
        assert math.isclose(mol.nuclear_repulsion_energy(),
                            9.34363815797054450919,
                            rel_tol=tol,
                            abs_tol=tol)

    def test_nuclear_repulsion_energy_with_ecp_basis(self):

        tol = 1.0e-12

        mol = Molecule.read_str(self.h2o_xyzstr(), 'au')
        basis = FakeBasis([2, 0, 0])

        coords = mol.get_coordinates_in_bohr()
        charges = np.array(mol.get_element_ids()) - np.array([2, 0, 0])
        ref_energy = 0.0

        for i in range(len(charges)):
            for j in range(i + 1, len(charges)):
                distance = np.linalg.norm(coords[j] - coords[i])
                ref_energy += charges[i] * charges[j] / distance

        assert math.isclose(mol.nuclear_repulsion_energy(basis),
                            ref_energy,
                            rel_tol=tol,
                            abs_tol=tol)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_nuclear_repulsion_energy_rejects_ecp_length_mismatch(self):

        mol = Molecule.read_str(self.h2o_xyzstr(), 'au')

        with pytest.raises(
                AssertionError,
                match='ECP core electron list must match number of atoms'):
            mol.nuclear_repulsion_energy(FakeBasis([2, 0]))

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_nuclear_repulsion_energy_rejects_negative_ecp_electrons(self):

        mol = Molecule.read_str(self.h2o_xyzstr(), 'au')

        with pytest.raises(AssertionError,
                           match='ECP core electrons must be non-negative'):
            mol.nuclear_repulsion_energy(FakeBasis([2, -1, 0]))

    def test_check_proximity(self):

        mol = self.nh3_molecule()
        assert mol.check_proximity(0.1)
        assert not mol.check_proximity(5.1)

    def test_get_string(self):

        mol = self.nh3_molecule()
        molstr = mol.get_string()
        lines = molstr.splitlines()
        assert lines[0].strip() == 'Molecular Geometry (Angstroms)'
        assert lines[1].strip() == '================================'
        assert lines[2].strip() == ''
        assert lines[3].lstrip().startswith('Atom')
        assert lines[3].rstrip().endswith('Coordinate Z')
        assert lines[4].strip() == ''
        assert lines[5].lstrip().startswith('N')
        assert lines[5].rstrip().endswith('-0.019579556803')
        assert lines[6].lstrip().startswith('H')
        assert lines[6].rstrip().endswith('0.031221455443')
        assert lines[7].lstrip().startswith('H')
        assert lines[7].rstrip().endswith('0.792178284722')
        assert lines[8].lstrip().startswith('H')
        assert lines[8].rstrip().endswith('-0.832395752750')
        assert lines[9].strip() == ''

    def test_get_xyz_string_formatting_and_aliases(self, tmp_path):

        mol = self.nh3_molecule_with_ghost_atom()
        xyz = mol.get_xyz_string(precision=4,
                                 comment='   comment with spaces   \nignored')
        lines = xyz.splitlines()

        assert lines[0].strip() == '4'
        assert lines[1].strip() == 'comment with spaces'
        assert lines[2].split()[0] == 'Bq_N'
        assert lines[2].rstrip().endswith('-0.0196')
        assert lines[3].rstrip().endswith('0.0312')
        assert mol.get_xyz_string(comment='x' *
                                  120).splitlines()[1].strip() == ('x' * 80)

        xyz_path = tmp_path / 'ghost.xyz'
        mol.write_xyz(xyz_path)

        assert Molecule.read_xyz(xyz_path) == Molecule.read_xyz_file(xyz_path)
        assert Molecule.from_xyz_string(xyz) == Molecule.read_xyz_string(xyz)
        assert Molecule.read_str(self.nh3_xyzstr(),
                                 'au') == Molecule.read_molecule_string(
                                     self.nh3_xyzstr(), 'au')

    def test_get_string_formats_ghost_atoms(self):

        mol = self.nh3_molecule_with_ghost_atom()
        lines = mol.get_string(title='Ghost Geometry', sep='-').splitlines()

        assert lines[0].strip() == 'Ghost Geometry (Angstroms)'
        assert lines[1].strip() == '----------------------------'
        assert lines[5].lstrip().startswith('Bq_N')
        assert lines[5].rstrip().endswith('-0.019579556803')

    def test_read_dict(self):

        rxyz = [
            'N  -3.710   3.019  -0.037', 'H  -3.702   4.942   0.059',
            'H  -4.704   2.415   1.497', 'H  -4.780   2.569  -1.573'
        ]
        mdict = {"xyz": rxyz, "charge": 3.0, "multiplicity": 2, "units": "au"}

        mol_a = Molecule.from_input_dict(mdict)
        mol_b = self.nh3_molecule()
        mol_b.set_charge(3.0)
        mol_b.set_multiplicity(2)
        assert mol_a == mol_b

    def test_read_dict_xyzfile(self, tmp_path):

        xyz_path = tmp_path / 'ammonia.xyz'
        xyz_path.write_text(
            dedent("""
                4
                ammonia
                N  -3.710   3.019  -0.037
                H  -3.702   4.942   0.059
                H  -4.704   2.415   1.497
                H  -4.780   2.569  -1.573
                """).lstrip())

        mol_a = Molecule.from_input_dict({
            "xyzfile": str(xyz_path),
            "charge": 3.0,
            "multiplicity": 2
        })
        mol_b = Molecule.read_molecule_string(self.nh3_xyzstr(), 'angstrom')
        mol_b.set_charge(3.0)
        mol_b.set_multiplicity(2)
        assert mol_a == mol_b

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_read_dict_rejects_xyz_and_xyzfile(self, tmp_path):

        xyz_path = tmp_path / 'water.xyz'
        xyz_path.write_text(
            dedent("""
                3
                water
                O  0.0  0.0  -1.0
                H  0.0  1.4  -2.1
                H  0.0 -1.4  -2.1
                """).lstrip())

        with pytest.raises(AssertionError,
                           match='Cannot have both "xyz" and "xyzfile" input'):
            Molecule.from_input_dict({
                "xyz": ['H 0.0 0.0 0.0'],
                "xyzfile": str(xyz_path)
            })

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_read_dict_rejects_units_with_xyzfile(self, tmp_path):

        xyz_path = tmp_path / 'water.xyz'
        xyz_path.write_text(
            dedent("""
                3
                water
                O  0.0  0.0  -1.0
                H  0.0  1.4  -2.1
                H  0.0 -1.4  -2.1
                """).lstrip())

        with pytest.raises(
                AssertionError,
                match='Cannot have both "units" and "xyzfile" input'):
            Molecule.from_input_dict({"xyzfile": str(xyz_path), "units": "au"})

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_read_dict_rejects_invalid_multiplicity_parity(self):

        with pytest.raises(
                AssertionError,
                match='Incompatible multiplicity and number of electrons'):
            Molecule.from_input_dict({
                "xyz": [
                    'N  -3.710   3.019  -0.037', 'H  -3.702   4.942   0.059',
                    'H  -4.704   2.415   1.497', 'H  -4.780   2.569  -1.573'
                ],
                "multiplicity": 2
            })

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_read_dict_rejects_close_atoms(self):

        with pytest.raises(AssertionError, match='Atoms too close'):
            Molecule.from_input_dict(
                {"xyz": ['H  0.0  0.0  0.0', 'H  0.0  0.0  0.05']})

    def test_moments_of_inertia(self):

        tol = 1.0e-12

        mol = Molecule.read_str(self.ch4_xyzstr(), 'au')
        imoms = mol.moments_of_inertia()
        rmoms = np.array(
            [3.198919866723860, 3.198919866723860, 3.198919866723860])
        assert np.allclose(rmoms, imoms, tol, tol, False)

    def test_moments_of_inertia_principal_axes(self):

        tol = 1.0e-12
        mol = Molecule.read_molecule_string(
            """O  0.000  0.000  0.000
               H  0.000  0.757  0.586
               H  0.000 -0.757  0.586""", 'angstrom')

        imoms, axes = mol.moments_of_inertia(principal_axes=True)

        masses = np.array(mol.get_masses())
        coords = mol.get_coordinates_in_bohr()
        center_of_mass = np.array(mol.center_of_mass_in_bohr())
        coords_com = coords - center_of_mass[np.newaxis, :]
        inertia_terms = []
        for i in range(mol.number_of_atoms()):
            coord = coords_com[i]
            diagonal = np.eye(3) * np.dot(coord, coord)
            outer = np.outer(coord, coord)
            inertia_terms.append(masses[i] * (diagonal - outer))
        inertia = np.sum(inertia_terms, axis=0)

        assert np.allclose(np.matmul(axes, axes.T), np.eye(3), tol, tol, False)
        principal_inertia = np.matmul(np.matmul(axes, inertia), axes.T)
        assert np.allclose(principal_inertia, np.diag(imoms), tol, tol, False)

    def test_is_linear(self):

        mol = Molecule.read_str(self.ch4_xyzstr(), 'au')
        assert not mol.is_linear()
        mol = Molecule.read_str(self.co2_xyzstr(), 'au')
        assert mol.is_linear()

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_is_linear_rejects_less_than_two_atoms(self):

        mol = Molecule.read_molecule_string('He 0.0 0.0 0.0')
        with pytest.raises(AssertionError,
                           match='Molecule.is_linear: Need at least two atoms'):
            mol.is_linear()

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
        assert lines[0].strip() == 'Molecular charge            : 0'
        assert lines[1].strip() == 'Spin multiplicity           : 3'
        assert lines[2].strip() == 'Number of atoms             : 4'
        assert lines[3].strip() == 'Number of alpha electrons   : 6'
        assert lines[4].strip() == 'Number of beta  electrons   : 4'

    def test_more_info_fixed_width_output(self):

        mol = self.nh3_molecule()
        mol.set_charge(-1.0)
        mol.set_multiplicity(2)
        lines = mol.more_info().splitlines()

        assert all(len(line) == 70 for line in lines)
        assert lines[0].strip() == 'Molecular charge            : -1'
        assert lines[1].strip() == 'Spin multiplicity           : 2'

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

    def test_number_of_occupied_orbitals_with_ecp_basis(self):

        mol = self.nh3_molecule()
        basis = FakeBasis([2, 0, 0, 0])

        assert mol.number_of_alpha_occupied_orbitals(basis) == 4
        assert mol.number_of_beta_occupied_orbitals(basis) == 4

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

    def test_get_aufbau_alpha_occupation_with_ecp_basis(self):

        tol = 1.0e-12

        mol = self.nh3_molecule()
        basis = FakeBasis([2, 0, 0, 0])

        nocc = mol.get_aufbau_alpha_occupation(5, basis)
        rocc = np.array([1.0, 1.0, 1.0, 1.0, 0.0])
        assert np.allclose(rocc, nocc, tol, tol, False)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_get_aufbau_alpha_occupation_rejects_too_few_mos_with_basis(self):

        mol = self.nh3_molecule()
        basis = FakeBasis([2, 0, 0, 0])

        with pytest.raises(AssertionError,
                           match='Number of molecular orbitals is too small ' +
                           'for explicit alpha electrons'):
            mol.get_aufbau_alpha_occupation(3, basis)

    def test_get_aufbau_beta_occupation(self):

        tol = 1.0e-12

        mol = Molecule.read_str(self.nh3_xyzstr())
        mol.set_charge(1.0)
        mol.set_multiplicity(2)

        nocc = mol.get_aufbau_beta_occupation(6)
        rocc = np.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.0])
        assert np.allclose(rocc, nocc, tol, tol, False)

    def test_get_aufbau_beta_occupation_with_ecp_basis(self):

        tol = 1.0e-12

        mol = self.nh3_molecule()
        basis = FakeBasis([2, 0, 0, 0])

        nocc = mol.get_aufbau_beta_occupation(5, basis)
        rocc = np.array([1.0, 1.0, 1.0, 1.0, 0.0])
        assert np.allclose(rocc, nocc, tol, tol, False)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_get_aufbau_beta_occupation_rejects_too_few_mos_with_basis(self):

        mol = self.nh3_molecule()
        basis = FakeBasis([2, 0, 0, 0])

        with pytest.raises(AssertionError,
                           match='Number of molecular orbitals is too small ' +
                           'for explicit beta electrons'):
            mol.get_aufbau_beta_occupation(3, basis)

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

    def test_get_connectivity_matrix(self):

        mol = Molecule.read_str(self.nh3_h2o_xyzstr(), 'au')
        connectivity = mol.get_connectivity_matrix()
        ref_connectivity = np.array([
            [0, 1, 1, 1, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 1],
            [0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0, 0],
        ])
        assert np.array_equal(connectivity, ref_connectivity)

    def test_get_connectivity_matrix_uses_h2_special_factor(self):

        mol = Molecule.read_molecule_string('H 0.0 0.0 0.0\nH 0.0 0.0 0.7',
                                            'angstrom')

        assert np.array_equal(mol.get_connectivity_matrix(),
                              np.array([[0, 1], [1, 0]]))
        assert np.array_equal(mol.get_connectivity_matrix(H2_factor=1.3),
                              np.array([[0, 0], [0, 0]]))

    def testfind_connected_atoms(self):

        mol = Molecule.read_str(self.nh3_h2o_xyzstr(), 'au')
        connectivity = mol.get_connectivity_matrix()

        assert mol.find_connected_atoms(0, connectivity) == {0, 1, 2, 3}
        assert mol.find_connected_atoms(4, connectivity) == {4, 5, 6}

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

    def test_vdw_radii_and_elem_ids(self):

        ref_labels = [
            'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
            'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Ga', 'Ge', 'As', 'Se',
            'Br', 'Kr', 'Rb', 'Sr', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs',
            'Ba', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra'
        ]
        ref_radii = [
            1.1, 1.4, 1.81, 1.53, 1.92, 1.7, 1.55, 1.52, 1.47, 1.54, 2.27, 1.73,
            1.84, 2.1, 1.8, 1.8, 1.75, 1.88, 2.75, 2.31, 1.87, 2.11, 1.85, 1.9,
            1.83, 2.02, 3.03, 2.49, 1.93, 2.17, 2.06, 2.06, 1.98, 2.16, 3.43,
            2.68, 1.96, 2.02, 2.07, 1.97, 2.02, 2.2, 3.48, 2.83
        ]

        vdw_data = np.array(get_vdw_radii_data_in_bohr()) * bohr_in_angstrom()

        for elem, radius in zip(ref_labels, ref_radii):
            if elem in ['Fr', 'Ra']:
                continue
            elem_id = chemical_element_identifier(elem.upper())
            assert abs(radius - vdw_data[elem_id]) < 1.0e-10

    @pytest.mark.skipif(not DispersionModel.is_available(),
                        reason='dftd4-python not available')
    def test_dispersion_model(self):

        ref_energies = {
            'b3lyp': (-0.2421986012033e-02, 0.4584078631740e-06),
            'blyp': (-0.3504378057810e-02, 0.9461341515956e-06),
            'hf': (-0.6671313029358e-02, 0.5473122492978e-05),
        }

        ref_gradient = {}

        ref_gradient['b3lyp'] = np.array([
            [-0.1299360785248e-03, 0.2173451050590e-03, -0.3709704540840e-05],
            [0.3994969804870e-05, -0.4278600323727e-04, -0.3004853785695e-05],
            [0.2248231831000e-04, 0.4826552264307e-04, -0.4026908304668e-04],
            [0.2585427033048e-04, 0.3687322138623e-04, 0.3605346888461e-04],
            [0.3668558637179e-04, -0.1301671081015e-03, 0.5463254511935e-05],
            [-0.3229412701673e-05, 0.4922085484071e-05, 0.5884251321327e-05],
            [0.1936253825266e-04, -0.5468305617267e-04, 0.4693862097277e-05],
            [0.1839250629302e-04, -0.6545014186048e-04, 0.3934710919238e-05],
            [0.6393301863664e-05, -0.1431962520046e-04, -0.9045906361170e-05],
        ])

        ref_gradient['blyp'] = np.array([
            [-0.1920479181714e-03, 0.2729735760896e-03, -0.2511659875033e-05],
            [0.9234017559838e-05, -0.7983770503499e-04, -0.4956308627117e-05],
            [0.4283192296416e-04, 0.6348404201255e-04, -0.6717298399052e-04],
            [0.4829075493265e-04, 0.4755520001250e-04, 0.6235137755103e-04],
            [0.4290323248918e-04, -0.1543006045860e-03, 0.6640266553842e-05],
            [-0.4698918137746e-05, 0.8662903487663e-05, 0.7308871341866e-05],
            [0.2376967810450e-04, -0.6486040126786e-04, 0.5846333490862e-05],
            [0.2249172537731e-04, -0.7949176018520e-04, 0.4807049578003e-05],
            [0.7225504881460e-05, -0.1418525052826e-04, -0.1231294602293e-04],
        ])

        ref_gradient['hf'] = np.array([
            [-0.3686765721713e-03, 0.3397769090112e-03, 0.4534616906695e-05],
            [0.3281830111142e-04, -0.2014976035809e-03, -0.1052292444208e-04],
            [0.1161479269369e-03, 0.9219122684571e-04, -0.1528315355729e-03],
            [0.1280011820621e-03, 0.6531340670539e-04, 0.1479566767179e-03],
            [0.4036674113369e-04, -0.1594768392328e-03, 0.1040317973416e-04],
            [0.4515452238650e-05, -0.5312330663685e-05, -0.6266158009138e-05],
            [0.4492572018223e-05, -0.5927187807585e-04, -0.3729549690795e-05],
            [0.3911077473617e-04, -0.6105721802557e-04, 0.2269617645742e-05],
            [0.3223621934107e-05, -0.1066567298349e-04, 0.8186076710362e-05],
        ])

        here = Path(__file__).parent
        inpfile = here / 'data' / 'dimer.inp'

        task = MpiTask([str(inpfile), None])
        molecule = task.molecule

        disp = DispersionModel()

        for xc_label in ref_energies:
            disp.compute(molecule, xc_label)

            e_disp = disp.get_energy()
            e_ref = sum(ref_energies[xc_label])
            assert abs(e_disp - e_ref) < 1.0e-11

            g_disp = disp.get_gradient()
            g_ref = ref_gradient[xc_label]
            max_diff = np.max(np.abs(g_disp - g_ref))
            assert max_diff < 1.0e-10

    def test_distance_matrix(self):

        xyz_string = """
            3
            water
            O       0.944703950608    0.070615471110   -0.056477643546
            H       1.903732770975    0.056427544902   -0.044142825330
            H       0.642753278417   -0.616673016012    0.540640468876
        """
        mol = Molecule.read_xyz_string(xyz_string)
        distance_matrix = mol.get_distance_matrix_in_angstrom()

        ref_distrance_matrix = np.array([[0., 0.95921308, 0.95921307],
                                         [0.95921308, 0., 1.54437856],
                                         [0.95921307, 1.54437856, 0.]])

        assert np.max(np.abs(distance_matrix - ref_distrance_matrix)) < 1.0e-8

    def test_get_ic_rmsd(self):

        init_xyz_str = """
            C          -2.723479309037        0.388400197499       -0.032041680121
            C          -1.198591674506        0.292402159666       -0.001957955680
            H          -3.161468702971        0.002148459476        0.916810101463
            H          -3.113340035638       -0.249559972670       -0.853012496915
            O          -3.131199766507        1.709538730510       -0.265228909885
            H          -0.880032285304       -0.767190536848        0.124520689502
            O          -0.670239980035        1.084030100114        1.027238801838
            H          -0.785679988597        0.645580322006       -0.970638007358
            H          -0.821288323120        0.587492539563        1.873631291844
            H          -3.112122928053        2.174373284354        0.611781773545
        """

        final_xyz_str = """
            C          -2.760169377434        0.322719976286       -0.001925407564
            C          -1.228941816065        0.304680448557       -0.002845411826
            H          -3.128135202340       -0.132424338111        0.933568102246
            H          -3.149271541311       -0.251138690089       -0.845341019529
            O          -3.252173881237        1.631288380215       -0.127670008266
            H          -0.863067365090       -0.729220194991        0.037102867132
            O          -0.717621604087        1.071775860470        1.064783463761
            H          -0.854080755853        0.782453035355       -0.910452256526
            H          -0.974559560638        0.659440589631        1.897920584004
            H          -2.669421889714        2.207640216349        0.385962694800
        """

        init_mol = Molecule.read_molecule_string(init_xyz_str, units='angstrom')
        final_mol = Molecule.read_molecule_string(final_xyz_str,
                                                  units='angstrom')

        ic_rmsd = OptimizationDriver.get_ic_rmsd(final_mol, init_mol)
        ref_ic_rmsd = {
            'bonds': {
                'rms': 0.017,
                'max': 0.028,
                'unit': 'Angstrom'
            },
            'angles': {
                'rms': 1.463,
                'max': 3.234,
                'unit': 'degree'
            },
            'dihedrals': {
                'rms': 20.573,
                'max': 45.089,
                'unit': 'degree'
            }
        }

        for ic in ['bonds', 'angles', 'dihedrals']:
            for val in ['rms', 'max']:
                assert abs(ic_rmsd[ic][val] - ref_ic_rmsd[ic][val]) < 1.0e-3

    def test_write_xyz(self):

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            here = Path(__file__).parent
            random_string = get_random_string_serial()
            fpath = here / 'data' / f'vlx_molecule_{random_string}.xyz'
            fname = str(fpath)

            for mol in [
                    self.nh3_molecule(),
                    self.nh3_molecule_with_ghost_atom(),
            ]:
                mol.write_xyz(fname)

                ref_labels = mol.get_labels()
                ref_coords_au = mol.get_coordinates_in_bohr()
                ref_coords_ang = mol.get_coordinates_in_angstrom()

                mol_2 = Molecule.read_xyz_file(fname)
                assert ref_labels == mol_2.get_labels()
                assert np.max(
                    np.abs(ref_coords_au -
                           mol_2.get_coordinates_in_bohr())) < 1.0e-10
                assert np.max(
                    np.abs(ref_coords_ang -
                           mol_2.get_coordinates_in_angstrom())) < 1.0e-10

                with open(fname, 'r') as f_xyz:
                    lines = f_xyz.readlines()
                    mol_3 = Molecule.read_xyz_string(''.join(lines))
                    assert ref_labels == mol_3.get_labels()
                    assert np.max(
                        np.abs(ref_coords_au -
                               mol_3.get_coordinates_in_bohr())) < 1.0e-10
                    assert np.max(
                        np.abs(ref_coords_ang -
                               mol_3.get_coordinates_in_angstrom())) < 1.0e-10

                if fpath.is_file():
                    fpath.unlink()

    def test_write_xyz_file_serializes_ghost_atoms(self, tmp_path):

        mol = self.nh3_molecule_with_ghost_atom()
        xyz_path = tmp_path / 'ghost_written.xyz'
        mol.write_xyz_file(xyz_path)

        lines = xyz_path.read_text().splitlines()
        assert lines[0] == '4'
        assert lines[2].split()[0] == 'Bq_N'

        mol.write_xyz(xyz_path)
        assert Molecule.read_xyz_file(xyz_path) == mol

    def test_rotate_around_vector(self):

        mol = self.nh3_molecule()
        coords = np.array([[1.0, 0.0, 0.0], [2.0, 1.0, 0.0]])
        origin = np.array([0.0, 0.0, 0.0])
        vector = np.array([0.0, 0.0, 1.0])

        rotated = mol._rotate_around_vector(coords, origin, vector, 90.0,
                                            'degree')

        assert np.allclose(rotated,
                           np.array([
                               [0.0, 1.0, 0.0],
                               [-1.0, 2.0, 0.0],
                           ]), 1.0e-12, 1.0e-12, False)

    def test_get_input_keywords(self):

        assert Molecule._get_input_keywords() == {
            'molecule': {
                'charge': ('int', 'net charge'),
                'multiplicity': ('int', 'spin multiplicity'),
                'units':
                    ('str_lower', 'unit of coordinates, default is Angstrom'),
                'xyz': ('list', 'atom and Cartesian coordinates'),
                'xyzfile': ('str', 'XYZ file name (conflicts with units/xyz)'),
            },
        }

    def test_partition_atoms_comm_self(self):

        mol = Molecule.read_molecule_string("""H  0.0  0.0  0.0
                                               O  1.0  0.0  0.0
                                               Li 2.0  0.0  0.0""")
        assert mol.partition_atoms(MPI.COMM_SELF) == [1, 2, 0]

    def test_distance_angle_dihedral(self):

        xyz_string = """
            9
            xyz
            O    1.086900000000    0.113880000000   -0.060730000000
            C    2.455250000000    0.132120000000   -0.071390000000
            C    3.171673900000   -0.838788100000    0.496389800000
            C    2.491492369403   -1.966504408464    1.155862765438
            O    1.664691816845   -2.650313648401    0.565927537003
            H    0.786520000000   -0.686240000000    0.407170000000
            H    2.871553600000    0.995167700000   -0.576074500000
            H    4.254316047964   -0.842628605717    0.498660431748
            H    2.767678583706   -2.159148582998    2.205812810612
        """
        mol = Molecule.read_xyz_string(xyz_string)
        assert abs(mol.get_dihedral_in_degrees((2, 3, 4, 5)) - 55.0) < 1e-4

        mol.set_dihedral_in_degrees((2, 3, 4, 5), 270.0)
        assert abs(mol.get_dihedral((2, 3, 4, 5), 'radian') +
                   math.pi / 2.0) < 1e-4

        mol.set_dihedral((2, 3, 4, 5), math.pi / 2.0, 'radian')
        assert abs(mol.get_dihedral_in_degrees((2, 3, 4, 5)) - 90.0) < 1e-4

        assert abs(mol.get_angle_in_degrees((1, 2, 3)) - 121.6) < 1e-1
        mol.set_angle_in_degrees((1, 2, 3), 110.0)
        assert abs(mol.get_angle_in_degrees((1, 2, 3)) - 110.0) < 1e-4

        assert abs(mol.get_distance_in_angstroms((1, 2)) - 1.3685) < 1e-4
        mol.set_distance_in_angstroms((1, 2), 1.7)
        assert abs(mol.get_distance_in_angstroms((1, 2)) - 1.7) < 1e-4

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_distance_angle_dihedral_invalid_arguments(self):

        mol = Molecule.read_molecule_string(self.nh3_xyzstr(), 'angstrom')

        with pytest.raises(
                AssertionError,
                match='Molecule.get_distance: Expecting two atom indices'):
            mol.get_distance((1,), 'angstrom')

        with pytest.raises(
                AssertionError,
                match='Molecule.get_distance: Invalid distance unit'):
            mol.get_distance((1, 2), 'nanometer')

        with pytest.raises(
                AssertionError,
                match='Molecule.set_distance: Expecting two atom indices'):
            mol.set_distance((1,), 1.1, 'angstrom')

        with pytest.raises(
                AssertionError,
                match='Molecule.set_distance: Invalid distance unit'):
            mol.set_distance((1, 2), 1.1, 'nanometer')

        with pytest.raises(AssertionError,
                           match='Molecule.get_angle: Expecting three atom '
                           'indices'):
            mol.get_angle((1, 2), 'degree')

        with pytest.raises(AssertionError,
                           match='Molecule.get_angle: Invalid angle unit'):
            mol.get_angle((1, 2, 3), 'gradian')

        with pytest.raises(AssertionError,
                           match='Molecule.set_angle: Expecting three atom '
                           'indices'):
            mol.set_angle((1, 2), 100.0, 'degree')

        with pytest.raises(AssertionError,
                           match='Molecule.set_angle: Invalid angle unit'):
            mol.set_angle((1, 2, 3), 100.0, 'gradian')

        with pytest.raises(AssertionError,
                           match='Molecule.get_dihedral: Expecting four atom '
                           'indices'):
            mol.get_dihedral((1, 2, 3), 'degree')

        with pytest.raises(AssertionError,
                           match='Molecule.get_dihedral: Invalid angle unit'):
            mol.get_dihedral((1, 2, 3, 4), 'gradian')

        with pytest.raises(AssertionError,
                           match='Molecule.set_dihedral: Expecting four atom '
                           'indices'):
            mol.set_dihedral((1, 2, 3), 120.0, 'degree')

        with pytest.raises(AssertionError,
                           match='Molecule.set_dihedral: Invalid angle unit'):
            mol.set_dihedral((1, 2, 3, 4), 120.0, 'gradian')

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_geometry_setters_reject_ring_rotations(self):

        cyclopropane = Molecule.read_molecule_string(self.cyclopropane_xyzstr(),
                                                     'angstrom')

        with pytest.raises(AssertionError,
                           match='Molecule.set_distance: Cannot set '
                           'distance'):
            cyclopropane.set_distance((1, 2), 1.6, 'angstrom')

        with pytest.raises(AssertionError,
                           match='Molecule.set_angle: Cannot set angle'):
            cyclopropane.set_angle((1, 2, 3), 70.0, 'degree')

        cyclobutane = Molecule.read_molecule_string(self.cyclobutane_xyzstr(),
                                                    'angstrom')

        with pytest.raises(AssertionError,
                           match='Molecule.set_dihedral: Cannot set '
                           'dihedral'):
            cyclobutane.set_dihedral((1, 2, 3, 4), 30.0, 'degree')

    def test_angle_and_dihedral_target_wraparound(self):

        xyz_string = """
            9
            xyz
            O    1.086900000000    0.113880000000   -0.060730000000
            C    2.455250000000    0.132120000000   -0.071390000000
            C    3.171673900000   -0.838788100000    0.496389800000
            C    2.491492369403   -1.966504408464    1.155862765438
            O    1.664691816845   -2.650313648401    0.565927537003
            H    0.786520000000   -0.686240000000    0.407170000000
            H    2.871553600000    0.995167700000   -0.576074500000
            H    4.254316047964   -0.842628605717    0.498660431748
            H    2.767678583706   -2.159148582998    2.205812810612
        """

        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_angle_in_degrees((1, 2, 3), 470.0, verbose=False)
        assert abs(mol.get_angle_in_degrees((1, 2, 3)) - 110.0) < 1e-4

        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_dihedral_in_degrees((2, 3, 4, 5), 450.0, verbose=False)
        assert abs(mol.get_dihedral_in_degrees((2, 3, 4, 5)) - 90.0) < 1e-4

        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_dihedral((2, 3, 4, 5),
                         -3.0 * math.pi / 2.0,
                         'radian',
                         verbose=False)
        assert abs(mol.get_dihedral((2, 3, 4, 5), 'radian') -
                   math.pi / 2.0) < 1e-4

    @pytest.mark.skipif("rdkit" not in sys.modules,
                        reason="rdkit not available")
    def test_is_water_molecule(self):

        mol = Molecule.read_smiles('O')
        assert mol.is_water_molecule()

        mol = Molecule.read_smiles('CO')
        assert not mol.is_water_molecule()

        mol = Molecule.read_smiles('OO')
        assert not mol.is_water_molecule()

    @pytest.mark.skipif("rdkit" not in sys.modules,
                        reason="rdkit not available")
    def test_smiles_to_xyz_returns_xyz_for_both_hydrogen_options(self):

        xyz_with_hydrogen = Molecule.smiles_to_xyz('O', hydrogen=True)
        xyz_without_hydrogen = Molecule.smiles_to_xyz('O', hydrogen=False)

        assert isinstance(xyz_with_hydrogen, str)
        assert isinstance(xyz_without_hydrogen, str)
        assert int(xyz_with_hydrogen.splitlines()[0]) == 3
        assert int(xyz_without_hydrogen.splitlines()[0]) == 1

    @pytest.mark.skipif("rdkit" not in sys.modules,
                        reason="rdkit not available")
    def test_read_smiles_with_multi_components(self):

        mol = Molecule.read_smiles('CCO.C1CCCC1.c1ccccc1')
        assert mol.number_of_atoms() == 36

    @pytest.mark.skipif("rdkit" not in sys.modules,
                        reason="rdkit not available")
    def test_read_smiles_with_reorder_hydrogens(self):

        # with reorder_hydrogens, the first atom will be C and
        # the second atom will be H
        mol = Molecule.read_smiles('CCO', reorder_hydrogens=True)
        labels = mol.get_labels()
        assert labels[0] != 'H'
        assert labels[1] == 'H'
        assert labels[2] == 'H'

        # H2 molecule
        mol = Molecule.read_smiles('[H][H]', reorder_hydrogens=True)
        labels = mol.get_labels()
        assert labels[0] == 'H'
        assert labels[1] == 'H'

    @pytest.mark.skipif("rdkit" not in sys.modules,
                        reason="rdkit not available")
    def test_draw_2d_uses_display(self, monkeypatch):

        ipython_display = pytest.importorskip("IPython.display")
        display_mock = Mock()
        monkeypatch.setattr(ipython_display, "display", display_mock)

        Molecule.draw_2d('O', width=200, height=100)

        display_mock.assert_called_once()

    def test_is_water_molecule_requires_water_connectivity(self):

        mol = Molecule.read_molecule_string("""O  0.000  0.000  0.000
                                               H  0.000  0.757  0.586
                                               H  0.000 -0.757  0.586""")
        assert mol.is_water_molecule()

        mol = Molecule.read_molecule_string("""O  0.000  0.000  0.000
                                               H  0.000  0.000  5.000
                                               H  0.000  5.000  0.000""")
        assert not mol.is_water_molecule()

    def test_contains_water_molecule(self):

        mol = Molecule.read_molecule_string(self.nh3_h2o_xyzstr(), 'au')
        assert mol.contains_water_molecule()

        mol = Molecule.read_molecule_string(self.nh3_xyzstr(), 'au')
        assert not mol.contains_water_molecule()

        mol = Molecule.read_molecule_string(
            """N  -3.710   3.019  -0.037
                                               H  -3.702   4.942   0.059
                                               H  -4.704   2.415   1.497
                                               H  -4.780   2.569  -1.573
                                               O   0.000   0.000  -1.000
                                               H   0.000   0.000   4.000
                                               H   0.000   4.000  -1.000""",
            'au')
        assert not mol.contains_water_molecule()
