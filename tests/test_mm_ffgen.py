from mpi4py import MPI
from pathlib import Path
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator
from veloxchem.errorhandler import VeloxChemError

skip_multi_rank_raises = pytest.mark.skipif(
    MPI.COMM_WORLD.Get_size() > 1,
    reason='skip pytest.raises for multiple MPI processes')


class TestMMForceFieldGenerator:

    @staticmethod
    def _ffgen_test_molecule():

        xyzstr = """10
            xyz
            C        1.560000    -0.075662     2.503629
            C        1.255506     0.490597     1.343469
            O        1.434318    -0.214973     0.154595
            C        1.118147     0.362263    -1.104238
            H        1.422185     0.469784     3.428142
            H        1.948080    -1.085628     2.537131
            H        0.867900     1.502691     1.326870
            H        1.336865    -0.370909    -1.907478
            H        1.732304     1.272846    -1.269421
            H        0.040051     0.626884    -1.140597
        """
        return Molecule.read_xyz_string(xyzstr)

    @skip_multi_rank_raises
    def test_ffgen_rejects_isolated_water_without_explicit_water_model(self):

        pytest.importorskip('rdkit')

        water = Molecule.read_smiles('O')

        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()

        with pytest.raises(VeloxChemError,
                           match='explicit water_model is required'):
            ff_gen.create_topology(water)

    def test_ffgen_water(self):

        pytest.importorskip('rdkit')

        water = Molecule.read_smiles('O')

        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()
        ff_gen.create_topology(water, water_model='cspce')

        if MPI.COMM_WORLD.Get_rank() == mpi_master():

            fname = '_vlx_test_water_'
            ff_gen.write_openmm_files(fname)

            here = Path(__file__).parent
            xmlfile_ref = str(here / 'data' / 'mmffgen_water.xml')

            assert self.match_xml(f'{fname}.xml', xmlfile_ref, 1e-10)

            self.clean_up_xml_and_pdb(fname)

    def test_ffgen_transition_metal(self):

        xyzstr = """25
        xyz
        Cl   -2.1166   -0.1503   -4.1637
        Ni   -1.8895    0.2693   -2.0458
        N    -1.6050    0.6134   -0.1842
        C    -2.4664    1.0748    0.7434
        C    -2.0815    1.2718    2.0886
        C    -0.7782    0.9909    2.5055
        C     0.1541    0.5027    1.5504
        C    -0.3211    0.3348    0.2191
        C     0.5368   -0.1488   -0.8081
        C     1.8942   -0.4781   -0.5331
        C     2.6710   -0.9531   -1.6243
        C     2.0665   -1.0663   -2.8785
        C     0.7098   -0.7158   -3.0608
        N    -0.0508   -0.2627   -2.0451
        Cl   -3.9666    0.8936   -1.9495
        H    -3.4816    1.2810    0.3624
        H    -2.8315    1.6511    2.7988
        H     2.6351   -1.4301   -3.7475
        H     0.1854   -0.7864   -4.0295
        C     1.5316    0.1656    1.8165
        C     2.3641   -0.3035    0.8199
        H     3.4123   -0.5541    1.0477
        H     1.9135    0.2905    2.8421
        H    -0.4707    1.1419    3.5522
        H     3.7280   -1.2245   -1.4752
        """

        tm_complex = Molecule.read_xyz_string(xyzstr)

        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()
        ff_gen.create_topology(tm_complex, resp=False)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():

            fname = '_vlx_test_tm_complex_'
            ff_gen.write_openmm_files(fname)

            here = Path(__file__).parent
            xmlfile_ref = str(here / 'data' / 'mmffgen_tm_noresp.xml')

            assert self.match_xml(f'{fname}.xml', xmlfile_ref, 1e-10)

            self.clean_up_xml_and_pdb(fname)

    def clean_up_xml_and_pdb(self, fname):

        ffgen_xml_file = Path(f'{fname}.xml')
        if ffgen_xml_file.is_file():
            ffgen_xml_file.unlink()

        ffgen_pdb_file = Path(f'{fname}.pdb')
        if ffgen_pdb_file.is_file():
            ffgen_pdb_file.unlink()

    def match_xml(self, fname_data, fname_ref, tolerance):

        data = self.read_xml(fname_data)
        data_ref = self.read_xml(fname_ref)

        if sorted(list(data.keys())) != sorted(list(data_ref.keys())):
            return False

        for key in data:
            assert len(data[key]) == len(data_ref[key])
            for val, val_ref in zip(data[key], data_ref[key]):
                if isinstance(val, str):
                    if val != val_ref:
                        return False
                else:
                    if abs(val - val_ref) >= tolerance:
                        return False

        return True

    @staticmethod
    def read_xml(fname):

        # This is a simple function that reads an xml file and collects the
        # double-quoted text.

        with open(fname, "r") as fh:
            xml_lines = fh.readlines()

        data = {}

        for line in xml_lines:
            if 'xml version' in line:
                continue

            if '="' not in line:
                continue

            key = line.split()[0].strip().replace('<', '')
            assert len(key) > 0

            if key not in data:
                data[key] = []

            content = line.split('"')
            assert len(content) > 0

            for idx in range((len(content) - 1) // 2):
                valstr = content[idx * 2 + 1]
                try:
                    val = float(valstr)
                    data[key].append(val)
                except ValueError:
                    data[key].append(valstr)

        return data

    def test_ffgen_set_params_and_get_params(self):

        pytest.importorskip('rdkit')

        mol = self._ffgen_test_molecule()

        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()
        ff_gen.create_topology(mol, resp=False)

        assert ff_gen.get_bond_params((1, 2)) == ff_gen.bonds[(0, 1)]
        assert ff_gen.get_angle_params((3, 2, 7)) == ff_gen.angles[(2, 1, 6)]
        assert ff_gen.get_dihedral_params(
            (6, 1, 2, 7)) == ff_gen.dihedrals[(5, 0, 1, 6)]

        bond_params = {
            'type': 'harmonic',
            'force_constant': 5.0e+5,
            'equilibrium': 0.2,
        }
        angle_params = {
            'type': 'harmonic',
            'force_constant': 450.0,
            'equilibrium': 100.0,
        }
        dihedral_params = {
            'type': 'Fourier',
            'multiple': False,
            'barrier': 20.0,
            'phase': 0.0,
            'periodicity': 2,
        }

        assert ff_gen.get_bond_params((1, 2)) != bond_params
        assert ff_gen.get_angle_params((3, 2, 7)) != angle_params
        assert ff_gen.get_dihedral_params((6, 1, 2, 7)) != dihedral_params

        ff_gen.set_bond_params((1, 2), bond_params)
        ff_gen.set_angle_params((3, 2, 7), angle_params)
        ff_gen.set_dihedral_params((6, 1, 2, 7), dihedral_params)

        assert ff_gen.get_bond_params((1, 2)) == bond_params
        assert ff_gen.get_angle_params((3, 2, 7)) == angle_params
        assert ff_gen.get_dihedral_params((6, 1, 2, 7)) == dihedral_params

    def test_ffgen_param_accessors_accept_reversed_indices(self):

        pytest.importorskip('rdkit')

        mol = self._ffgen_test_molecule()

        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()
        ff_gen.create_topology(mol, resp=False)

        bond_forward = ff_gen.get_bond_params((1, 2))
        angle_forward = ff_gen.get_angle_params((3, 2, 7))
        dihedral_forward = ff_gen.get_dihedral_params((6, 1, 2, 7))

        assert ff_gen.get_bond_params((2, 1)) == bond_forward
        assert ff_gen.get_angle_params((7, 2, 3)) == angle_forward
        assert ff_gen.get_dihedral_params((7, 2, 1, 6)) == dihedral_forward

        reversed_bond = {
            'type': 'harmonic',
            'force_constant': 5.5e+5,
            'equilibrium': 0.18,
        }
        reversed_angle = {
            'type': 'harmonic',
            'force_constant': 420.0,
            'equilibrium': 105.0,
        }
        reversed_dihedral = {
            'type': 'Fourier',
            'multiple': False,
            'barrier': 15.0,
            'phase': 180.0,
            'periodicity': 3,
        }

        ff_gen.set_bond_params((2, 1), reversed_bond)
        ff_gen.set_angle_params((7, 2, 3), reversed_angle)
        ff_gen.set_dihedral_params((7, 2, 1, 6), reversed_dihedral)

        assert ff_gen.get_bond_params((1, 2)) == reversed_bond
        assert ff_gen.get_angle_params((3, 2, 7)) == reversed_angle
        assert ff_gen.get_dihedral_params((6, 1, 2, 7)) == reversed_dihedral

        returned_dihedral = ff_gen.get_dihedral_params((6, 1, 2, 7))
        returned_dihedral['barrier'] = -1.0

        assert ff_gen.get_dihedral_params((6, 1, 2, 7)) == reversed_dihedral

    def test_forcefield_json_roundtrip_preserves_topology_data(self, tmp_path):

        ff_gen = MMForceFieldGenerator()
        ff_gen.atoms = {
            0: {
                'type': 'c3',
                'name': 'C1',
                'mass': 12.011,
                'charge': -0.12,
            },
            1: {
                'type': 'c2',
                'name': 'C2',
                'mass': 12.011,
                'charge': 0.02,
            },
            2: {
                'type': 'hc',
                'name': 'H1',
                'mass': 1.008,
                'charge': 0.12,
            },
            3: {
                'type': 'hc',
                'name': 'H2',
                'mass': 1.008,
                'charge': -0.02,
            },
        }
        ff_gen.bonds = {
            (0, 1): {
                'type': 'harmonic',
                'force_constant': 2.5e5,
                'equilibrium': 0.152,
                'comment': 'C-C',
            },
        }
        ff_gen.angles = {
            (0, 1, 2): {
                'type': 'harmonic',
                'force_constant': 350.0,
                'equilibrium': 109.5,
                'comment': 'test angle',
            },
        }
        ff_gen.dihedrals = {
            (0, 1, 2, 3): {
                'type': 'Fourier',
                'multiple': True,
                'barrier': [1.5, 0.5],
                'phase': [0.0, 180.0],
                'periodicity': [1, -3],
                'comment': ['first term', 'second term'],
            },
        }
        ff_gen.impropers = {
            (3, 1, 0, 2): {
                'type': 'Fourier',
                'multiple': False,
                'barrier': 2.0,
                'phase': 180.0,
                'periodicity': 2,
                'comment': 'improper',
            },
        }

        json_string = MMForceFieldGenerator.get_forcefield_as_json(ff_gen)
        loaded_from_string = MMForceFieldGenerator.load_forcefield_from_json_string(
            json_string)

        assert loaded_from_string.atoms == ff_gen.atoms
        assert loaded_from_string.bonds == ff_gen.bonds
        assert loaded_from_string.angles == ff_gen.angles
        assert loaded_from_string.dihedrals == ff_gen.dihedrals
        assert loaded_from_string.impropers == ff_gen.impropers

        json_path = tmp_path / 'nested' / 'forcefield.json'
        MMForceFieldGenerator.save_forcefield_as_json(ff_gen, str(json_path))
        loaded_from_file = MMForceFieldGenerator.load_forcefield_from_json_file(
            str(json_path))

        assert loaded_from_file.atoms == ff_gen.atoms
        assert loaded_from_file.bonds == ff_gen.bonds
        assert loaded_from_file.angles == ff_gen.angles
        assert loaded_from_file.dihedrals == ff_gen.dihedrals
        assert loaded_from_file.impropers == ff_gen.impropers

    @skip_multi_rank_raises
    def test_ffgen_canonicalize_helpers_reject_invalid_indices(self):

        with pytest.raises(VeloxChemError,
                           match='zero-based atom indices must be non-negative'):
            MMForceFieldGenerator._canonicalize_zero_based_bond_key((-1, 2))

        with pytest.raises(VeloxChemError,
                           match='zero-based atom indices must be non-negative'):
            MMForceFieldGenerator._canonicalize_zero_based_angle_key(
                [0, -1, 2])

        with pytest.raises(VeloxChemError,
                           match='zero-based atom indices must be non-negative'):
            MMForceFieldGenerator._canonicalize_zero_based_dihedral_key(
                (0, 1, 2, -3))

        with pytest.raises(VeloxChemError,
                           match='atom indices in a bond must be unique'):
            MMForceFieldGenerator._canonicalize_zero_based_bond_key((2, 2))

        with pytest.raises(VeloxChemError,
                           match='atom indices in an angle must be unique'):
            MMForceFieldGenerator._canonicalize_zero_based_angle_key(
                (0, 1, 0))

        with pytest.raises(VeloxChemError,
                           match='atom indices in a dihedral must be unique'):
            MMForceFieldGenerator._canonicalize_zero_based_dihedral_key(
                (0, 1, 0, 2))

        with pytest.raises(VeloxChemError,
                           match='one-based atom indices must be greater than 0'):
            MMForceFieldGenerator._to_zero_based_indices((1, 0))

        ff_gen = MMForceFieldGenerator()

        with pytest.raises(VeloxChemError,
                           match='one-based atom indices must be greater than 0'):
            ff_gen.get_bond_params((0, 1))

        with pytest.raises(VeloxChemError,
                           match='one-based atom indices must be greater than 0'):
            ff_gen.set_angle_params((1, -2, 3), {})

    def test_ffgen_update_settings_resolves_local_gaff_paths(self, tmp_path):

        ff_gen = MMForceFieldGenerator()
        ff_gen.molecule_name = str(tmp_path / 'inputs' / 'ligand.inp')

        force_field_file = tmp_path / 'inputs' / 'gaff2.dat'
        force_field_file.parent.mkdir(parents=True)
        force_field_file.write_text('; mock GAFF data\n')

        ff_gen.update_settings({
            'molecule_name': ff_gen.molecule_name,
            'force_field_data': 'gaff2.dat',
        })

        assert ff_gen.force_field_data == str(force_field_file)
        assert ff_gen.force_field_data_extension == str(
            force_field_file.with_name('gaff2_extension.dat'))
        assert ff_gen.resp_dict == {}

        ff_gen.update_settings({'filename': 'ligand_from_filename'})

        assert ff_gen.molecule_name == 'ligand_from_filename'

    @skip_multi_rank_raises
    def test_ffgen_update_settings_rejects_non_gaff_force_field(self, tmp_path):

        ff_gen = MMForceFieldGenerator()
        ff_gen.molecule_name = str(tmp_path / 'inputs' / 'ligand.inp')

        force_field_file = tmp_path / 'inputs' / 'amber99sb.dat'
        force_field_file.parent.mkdir(parents=True)
        force_field_file.write_text('; unsupported FF data\n')

        with pytest.raises(VeloxChemError, match='Only GAFF is supported'):
            ff_gen.update_settings({'force_field_data': 'amber99sb.dat'})

    @skip_multi_rank_raises
    def test_ffgen_rejects_malformed_scan_files(self, tmp_path):

        xyzstr = """10
            xyz
            C        1.560000    -0.075662     2.503629
            C        1.255506     0.490597     1.343469
            O        1.434318    -0.214973     0.154595
            C        1.118147     0.362263    -1.104238
            H        1.422185     0.469784     3.428142
            H        1.948080    -1.085628     2.537131
            H        0.867900     1.502691     1.326870
            H        1.336865    -0.370909    -1.907478
            H        1.732304     1.272846    -1.269421
            H        0.040051     0.626884    -1.140597
        """
        mol = Molecule.read_xyz_string(xyzstr)

        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()
        ff_gen.create_topology(mol, resp=False)

        missing_scan_file = tmp_path / '1-2-3-4.xyz'
        missing_scan_file.write_text("""2
        comment
        H 0.0 0.0 0.0
        H 0.0 0.0 0.7
        """)

        with pytest.raises(VeloxChemError,
                           match='scan file does not contain any Scan records'):
            ff_gen.reparameterize_dihedrals((1, 2),
                                            scan_file=str(missing_scan_file))

        wrong_name_file = tmp_path / 'wrong-name.xyz'
        wrong_name_file.write_text(
            'Scan Cycle 1/1 ; Dihedral 6-1-2-7 = 0.00 ; Iteration 1 Energy -1.0\n'
        )

        with pytest.raises(VeloxChemError,
                           match='scan file name does not match dihedral indices'):
            ff_gen.reparameterize_dihedrals((1, 2),
                                            scan_file=str(wrong_name_file))

    @skip_multi_rank_raises
    def test_ffgen_read_qm_scan_xyz_files_rejects_missing_scan_records(
            self, tmp_path):

        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()

        xyz_file = tmp_path / 'scan.xyz'
        xyz_file.write_text("""2
        comment
        H 0.0 0.0 0.0
        H 0.0 0.0 0.7
        """)

        with pytest.raises(VeloxChemError,
                           match='does not contain any Scan records'):
            ff_gen.read_qm_scan_xyz_files([xyz_file.name], inp_dir=tmp_path)

    @skip_multi_rank_raises
    def test_ffgen_read_qm_scan_xyz_files_rejects_inconsistent_scan_lengths(
            self, tmp_path):

        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()

        xyz_file = tmp_path / '1-2-3-4.xyz'
        xyz_file.write_text("2\n"
                            "Scan Cycle 1/1 ; Dihedral 1-2-3-4 = 0.00 ; "
                            "Iteration 1 Energy -1.0\n"
                            "H 0.0 0.0 0.0\n"
                            "H 0.0 0.0 0.7\n"
                            "Scan Cycle 2/2 ; Dihedral 1-2-3-4 = 180.00 ; "
                            "Iteration 1 Energy -0.5\n")

        with pytest.raises(VeloxChemError,
                           match='inconsistent number of geometries, energies, and dihedral angles'):
            ff_gen.read_qm_scan_xyz_files([xyz_file.name], inp_dir=tmp_path)

    @skip_multi_rank_raises
    def test_ffgen_reparameterize_dihedrals_rejects_scan_results_for_other_bond(
            self):

        pytest.importorskip('scipy')

        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()
        ff_gen.dihedrals = {
            (0, 1, 2, 3): {
                'multiple': False,
                'comment': 'c3-c3-c3-c3',
                'barrier': 0.5,
                'phase': 0.0,
                'periodicity': 1,
            },
        }

        scan_results = {
            'scan_dih_angles': [[0.0, 180.0]],
            'scan_energies': [[0.0, 1.0]],
            'scan_geometries': [[self._ffgen_test_molecule(),
                                 self._ffgen_test_molecule()]],
            'target_dihedrals': [[0, 4, 5, 6]],
        }

        with pytest.raises(VeloxChemError,
                           match='rotatable bond does not match the scan data'):
            ff_gen.reparameterize_dihedrals((2, 3),
                                            scan_results=scan_results,
                                            initial_validation=False,
                                            verbose=False)

    def test_ffgen_reparameterize_dihedrals_uses_matching_scan_results_entry(
            self):

        mol = self._ffgen_test_molecule()
        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()
        ff_gen.create_topology(mol, resp=False)

        rotatable_bond = (2, 3)
        matching_dihedral = None
        for dihedral in ff_gen.dihedrals:
            if sorted([dihedral[1] + 1, dihedral[2] + 1]) == list(
                    rotatable_bond):
                matching_dihedral = list(dihedral)
                break

        assert matching_dihedral is not None

        scan_results = {
            'scan_dih_angles': [[10.0, 20.0], [30.0, 60.0]],
            'scan_energies': [[0.1, 0.2], [0.3, 0.6]],
            'scan_geometries': [[mol, mol], [mol, mol]],
            'target_dihedrals': [[0, 4, 5, 6], matching_dihedral],
        }

        ff_gen.scan_dih_angles = scan_results['scan_dih_angles']
        ff_gen.scan_energies = scan_results['scan_energies']
        ff_gen.scan_geometries = scan_results['scan_geometries']
        ff_gen.target_dihedrals = scan_results['target_dihedrals']

        target_scan_index = ff_gen._get_target_scan_index(rotatable_bond)

        assert target_scan_index == 1
        assert ff_gen.scan_dih_angles[target_scan_index] == [30.0, 60.0]
        assert ff_gen.scan_energies[target_scan_index] == [0.3, 0.6]
        assert (ff_gen.scan_geometries[target_scan_index][0].get_xyz_string() ==
                mol.get_xyz_string())

    def test_ffgen_reparameterize_dihedrals_runs_with_matching_scan_results(
            self):

        pytest.importorskip('scipy')

        mol = self._ffgen_test_molecule()
        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()
        ff_gen.partial_charges = mol.get_partial_charges(mol.get_charge())
        ff_gen.create_topology(mol)

        rotatable_bond = (2, 3)
        matching_dihedral = None
        for dihedral in ff_gen.dihedrals:
            if sorted([dihedral[1] + 1, dihedral[2] + 1]) == list(
                    rotatable_bond):
                matching_dihedral = list(dihedral)
                break

        assert matching_dihedral is not None

        # Keep the constrained targets reasonably close to the starting
        # geometry. More aggressive jumps from this nearly planar structure can
        # trigger geomeTRIC internal-coordinate rebuild divergence in MPI runs.
        scan_results = {
            'scan_dih_angles': [[10.0, 20.0], [150.0, 120.0]],
            'scan_energies': [[0.1, 0.2], [0.3, 0.6]],
            'scan_geometries': [[mol, mol], [mol, mol]],
            'target_dihedrals': [[0, 4, 5, 6], matching_dihedral],
        }

        results = ff_gen.reparameterize_dihedrals(rotatable_bond,
                                                  scan_results=scan_results,
                                                  initial_validation=False,
                                                  verbose=False)

        assert results['dihedral_angles'].tolist() == [150.0, 120.0]

    def test_ffgen_file_writers_use_default_names_and_copy_file_creates_parent(
            self, tmp_path):

        ff_gen = MMForceFieldGenerator()
        ff_gen.molecule_name = str(tmp_path / 'ethanol_input.xyz')
        ff_gen.molecule = self._ffgen_test_molecule()
        ff_gen.atoms = {
            0: {
                'type': 'c3',
                'name': 'C1',
                'mass': 12.011,
                'charge': -0.2,
                'sigma': 0.34,
                'epsilon': 0.11,
                'equivalent_atom': 1,
            },
            1: {
                'type': 'hc',
                'name': 'H1',
                'mass': 1.008,
                'charge': 0.2,
                'sigma': 0.25,
                'epsilon': 0.03,
                'equivalent_atom': 1,
            },
        }
        ff_gen.unique_atom_types = ['c3', 'hc']
        ff_gen.bonds = {
            (0, 1): {
                'type': 'harmonic',
                'force_constant': 3.1e5,
                'equilibrium': 0.109,
                'comment': 'C-H',
            },
        }
        ff_gen.angles = {}
        ff_gen.dihedrals = {}
        ff_gen.impropers = {}
        ff_gen.pairs = []

        itp_file = tmp_path / 'defaults.itp'
        gro_file = tmp_path / 'defaults.gro'
        pdb_file = tmp_path / 'defaults.pdb'

        ff_gen.write_itp(itp_file)
        ff_gen.write_gro(gro_file)
        ff_gen.write_pdb(pdb_file)

        assert 'MOL' in itp_file.read_text()
        assert 'GRO file of MOL' in gro_file.read_text()
        assert 'PDB file of ethanol_input' in pdb_file.read_text()

        src = tmp_path / 'copied-source.txt'
        src.write_text('copied content')
        dest = tmp_path / 'nested' / 'dir' / 'copied-dest.txt'

        MMForceFieldGenerator.copy_file(src, dest)

        assert dest.read_text() == 'copied content'

    def test_ffgen_add_dihedral_promotes_existing_term_to_multiple(self):

        pytest.importorskip('rdkit')

        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()
        ff_gen.create_topology(self._ffgen_test_molecule(), resp=False)

        original_dihedral = ff_gen.get_dihedral_params((6, 1, 2, 7))

        assert original_dihedral['multiple'] is False

        ff_gen.add_dihedral((7, 2, 1, 6), barrier=2.5, phase=180.0, periodicity=3)

        updated_dihedral = ff_gen.get_dihedral_params((6, 1, 2, 7))

        assert updated_dihedral['multiple'] is True
        assert updated_dihedral['barrier'][0] == original_dihedral['barrier']
        assert updated_dihedral['phase'][0] == original_dihedral['phase']
        assert updated_dihedral['periodicity'][0] == original_dihedral[
            'periodicity']
        assert updated_dihedral['comment'][0] == original_dihedral['comment']
        assert updated_dihedral['barrier'][-1] == 2.5
        assert updated_dihedral['phase'][-1] == 180.0
        assert updated_dihedral['periodicity'][-1] == -3
        assert updated_dihedral['comment'][-1] == 'Added 6-1-2-7 dihedral'

    def test_ffgen_get_atom_names_gets_unique_suffixes(self):

        ff_gen = MMForceFieldGenerator()
        ff_gen.molecule = Molecule.read_xyz_string("""5
        labels
        C 0.0 0.0 0.0
        H 0.0 0.0 1.0
        H 1.0 0.0 0.0
        O 0.0 1.0 0.0
        N 1.0 1.0 1.0
        """)

        assert ff_gen.get_atom_names() == ['C', 'H1', 'H2', 'O', 'N']

    @skip_multi_rank_raises
    def test_ffgen_get_included_file_requires_include_directive(self, tmp_path):

        ff_gen = MMForceFieldGenerator()

        top_file = tmp_path / 'molecule.top'
        top_file.write_text('#include "molecule.itp"\n[ system ]\nTest\n')

        assert ff_gen.get_included_file(str(top_file)) == str(tmp_path /
                                                               'molecule.itp')

        missing_include_top = tmp_path / 'missing.top'
        missing_include_top.write_text('[ system ]\nTest\n')

        with pytest.raises(VeloxChemError, match='could not find included file'):
            ff_gen.get_included_file(str(missing_include_top))
