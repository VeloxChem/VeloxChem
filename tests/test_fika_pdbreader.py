import numpy as np
import pytest

from veloxchem.molecularbasis import MolecularBasis
from veloxchem.fikapdbreader import FikaPdbReader
from veloxchem.veloxchemlib import (bohr_in_angstrom, chemical_element_mass,
                                    FikaClassicalSystem, FikaForceField,
                                    FikaInducedDipolesDriver,
                                    FikaQmmmEmbeddingDriver,
                                    fika_classical_charges,
                                    fika_polarizable_sites)

ACROLEIN = [
    ('C', -0.145335, -0.546770, 0.000607),
    ('C', 1.274009, -0.912471, -0.000167),
    ('C', 1.630116, -2.207690, -0.000132),
    ('O', -0.560104, 0.608977, 0.000534),
    ('H', -0.871904, -1.386459, 0.001253),
    ('H', 2.004448, -0.101417, -0.000710),
    ('H', 0.879028, -3.000685, 0.000484),
    ('H', 2.675323, -2.516779, -0.000673),
]


def pdb_line(serial, atom, residue, chain, number, xyz, element):

    return (f'HETATM{serial:5d} {atom:<4} {residue:>3} {chain}{number:4d}    ' +
            f'{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00' +
            f'          {element:>2}')


def get_waters(count, seed=7):
    """
    Water molecules (O, H, H; angstrom) at random around acrolein, at least
    2.6 angstrom from its atoms and 2.8 angstrom apart (oxygens).
    """

    rng = np.random.default_rng(seed)
    atoms = np.array([a[1:] for a in ACROLEIN])
    centre = atoms.mean(axis=0)
    bond, angle = 0.9572, np.radians(104.52)
    local = np.array([[0.0, 0.0, 0.0],
                      [bond * np.sin(angle / 2), 0.0, bond * np.cos(angle / 2)],
                      [-bond * np.sin(angle / 2), 0.0, bond * np.cos(angle / 2)]])
    waters = []
    while len(waters) < count:
        oxygen = centre + rng.uniform(-11.0, 11.0, 3)
        if np.min(np.linalg.norm(atoms - oxygen, axis=1)) < 2.6:
            continue
        if waters and np.min(
                np.linalg.norm(np.array([w[0] for w in waters]) - oxygen,
                               axis=1)) < 2.8:
            continue
        rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        waters.append(oxygen + local @ rotation.T)
    return np.array(waters)


def get_pdb(waters, solute=True):
    """
    PDB text of acrolein (MOL, chain A) and the waters (HOH, chain W),
    the solute after the first ten waters.
    """

    lines = ['REMARK   test droplet']
    serial = 1
    for w, water in enumerate(waters):
        if w == 10 and solute:
            for name, x, y, z in ACROLEIN:
                lines.append(
                    pdb_line(serial, name + '1', 'MOL', 'A', 1, (x, y, z), name))
                serial += 1
            lines.append('TER')
        for atom, element, xyz in zip(['O', 'H1', 'H2'], ['O', 'H', 'H'],
                                      water):
            lines.append(
                pdb_line(serial, atom, 'HOH', 'W', w + 1, xyz, element))
            serial += 1
    lines.append('END')
    lines.append(pdb_line(serial, 'O', 'HOH', 'W', 9999, (0, 0, 0), 'O'))
    return '\n'.join(lines) + '\n'


def rounded_bohr(values):
    # PDB coordinates have three decimals; converted as the reader does.
    return np.round(np.asarray(values), 3) * (1.0 / bohr_in_angstrom())


def distances(waters):
    """
    Distances (angstrom) of the water centres of mass from acrolein's, with
    the rounded PDB coordinates.
    """

    def centre(elements, xyz):
        masses = np.array([chemical_element_mass(z) for z in elements])
        return masses @ xyz / np.sum(masses)

    solute = centre([6, 6, 6, 8, 1, 1, 1, 1],
                    rounded_bohr([a[1:] for a in ACROLEIN]))
    return np.array([
        np.linalg.norm(centre([8, 1, 1], rounded_bohr(w)) - solute)
        for w in waters
    ]) * bohr_in_angstrom()


def assert_error(text, line, message):

    with pytest.raises(ValueError) as error:
        FikaPdbReader().parse(text, 'test.pdb')
    assert str(error.value).startswith(f'test.pdb: line {line}: ')
    assert message in str(error.value)


class TestFikaPdbReader:

    def test_parse(self):

        waters = get_waters(20)
        residues = FikaPdbReader().parse(get_pdb(waters))

        # The record after END is not read.
        assert len(residues) == 21
        assert [r.name for r in residues[9:12]] == ['HOH', 'MOL', 'HOH']
        solute = residues[10]
        assert (solute.chain, solute.number, solute.insertion_code) == ('A', 1,
                                                                        ' ')
        assert solute.atom_names == [a[0] + '1' for a in ACROLEIN]
        assert solute.elements == [6, 6, 6, 8, 1, 1, 1, 1]
        water = residues[11]
        assert water.key() == ('HOH', 'W', 11, ' ')
        assert water.atom_names == ['O', 'H1', 'H2']
        assert water.elements == [8, 1, 1]
        assert np.max(np.abs(water.coordinates - rounded_bohr(waters[10]))) < 1e-12

    def test_errors(self):

        good = pdb_line(1, 'O', 'HOH', 'W', 1, (0, 0, 0), 'O')
        assert_error(good[:76] + '  ', 1, 'missing element symbol')
        assert_error(good[:76] + 'Qq', 1, "unknown element symbol 'Qq'")
        assert_error(good[:16] + 'A' + good[17:], 1, "alternate location 'A'")
        assert_error('MODEL        1\n' + good, 1, 'MODEL records are not supported')
        assert_error(good[:30] + '   x.000' + good[38:], 1, "invalid coordinate 'x.000'")
        assert_error(good[:17] + '   ' + good[20:], 1, 'missing residue name')
        assert_error(good[:22] + '  ab' + good[26:], 1, "invalid residue number 'ab'")
        assert_error(
            '\n'.join([
                good,
                pdb_line(2, 'O', 'HOH', 'W', 2, (3, 0, 0), 'O'),
                pdb_line(3, 'H1', 'HOH', 'W', 1, (1, 0, 0), 'H'),
            ]), 3, 'is split by another residue')
        with pytest.raises(ValueError, match='no ATOM or HETATM records'):
            FikaPdbReader().parse('REMARK nothing\nEND\n')
        with pytest.raises(FileNotFoundError):
            FikaPdbReader().read('/no/such/file.pdb')

    def test_solute(self, tmp_path):

        path = tmp_path / 'droplet.pdb'
        path.write_text(get_pdb(get_waters(12)))
        reader = FikaPdbReader()
        residues = reader.read(path)

        molecule = reader.get_solute(residues)
        assert molecule.get_labels() == [a[0] for a in ACROLEIN]
        assert np.max(
            np.abs(molecule.get_coordinates_in_bohr() -
                   rounded_bohr([a[1:] for a in ACROLEIN]))) < 1e-12

        with pytest.raises(ValueError, match='expected one HOH residue'):
            reader.get_solute(residues, 'HOH')
        with pytest.raises(ValueError, match='expected one MOL residue, found 0'):
            reader.get_solute(reader.parse(get_pdb(get_waters(3), False)))

    def test_all_solvent(self):

        reader = FikaPdbReader()
        waters = get_waters(15)

        # Without a solute: an MM-only system.
        residues = reader.parse(get_pdb(waters, solute=False))
        system, identities = reader.get_classical_system(residues,
                                                         {'HOH': 'tip3p'})
        assert system.number_of_residues(False) == 15
        assert system.number_of_residues(True) == 0
        assert identities['nonpolarizable'] == [('HOH', 'W', w + 1, ' ')
                                                for w in range(15)]
        charges, coordinates, offsets = fika_classical_charges(system)
        assert charges == [-0.834, 0.417, 0.417] * 15
        assert np.max(np.abs(coordinates - rounded_bohr(waters.reshape(-1, 3)))) < 1e-12

        # The solute is skipped; polarizable force fields fill the
        # polarizable region.
        residues = reader.parse(get_pdb(waters))
        system, identities = reader.get_classical_system(residues,
                                                         {'hoh': 'Ahlstrom'})
        assert system.number_of_residues(True) == 15
        positions, alphas, owners = fika_polarizable_sites(system)
        assert np.max(np.abs(positions - rounded_bohr(waters[:, 0, :]))) < 1e-12
        assert alphas == [9.718] * 15

    def test_radius_and_shells(self):

        reader = FikaPdbReader()
        # Residues in shuffled order, so selection must not depend on it.
        rng = np.random.default_rng(1)
        waters = get_waters(60)[rng.permutation(60)]
        residues = reader.parse(get_pdb(waters))
        d = distances(waters)
        inner, outer = 6.0, 9.0
        assert 0 < np.sum(d <= inner) < np.sum(d <= outer) < 60

        system, identities = reader.get_classical_system(residues,
                                                         {'HOH': 'ahlstrom'},
                                                         radius=outer)
        expected = [('HOH', 'W', w + 1, ' ') for w in np.flatnonzero(d <= outer)]
        assert identities['polarizable'] == expected
        assert system.number_of_residues(True) == len(expected)

        # The same radius in bohr.
        system, identities = reader.get_classical_system(
            residues, {'HOH': 'ahlstrom'},
            radius=outer / bohr_in_angstrom(),
            unit='bohr')
        assert identities['polarizable'] == expected

        system, identities = reader.get_classical_system(
            residues, {
                'polarizable': {
                    'HOH': 'ahlstrom'
                },
                'nonpolarizable': {
                    'HOH': 'tip3p'
                }
            },
            shells=(inner, outer))
        assert identities['polarizable'] == [
            ('HOH', 'W', w + 1, ' ') for w in np.flatnonzero(d <= inner)
        ]
        assert identities['nonpolarizable'] == [
            ('HOH', 'W', w + 1, ' ')
            for w in np.flatnonzero((d > inner) & (d <= outer))
        ]
        charges, _, _ = fika_classical_charges(system)
        assert charges == ([-0.669, 0.3345, 0.3345] *
                           len(identities['polarizable']) +
                           [-0.834, 0.417, 0.417] *
                           len(identities['nonpolarizable']))

        # An empty mapping skips its region.
        system, identities = reader.get_classical_system(
            residues, {'nonpolarizable': {
                'HOH': 'tip3p'
            }},
            shells=(inner, outer))
        assert identities['polarizable'] == []
        assert system.number_of_residues(False) == np.sum((d > inner) &
                                                          (d <= outer))

    def test_system_errors(self, tmp_path, monkeypatch):

        reader = FikaPdbReader()
        residues = reader.parse(get_pdb(get_waters(12)))
        shells = {'polarizable': {'HOH': 'ahlstrom'}}

        def raises(message, *args, **kwargs):
            with pytest.raises(ValueError, match=message):
                reader.get_classical_system(residues, *args, **kwargs)

        raises('either a radius or shells', {'HOH': 'tip3p'},
               radius=5.0,
               shells=(1.0, 2.0))
        raises("unknown unit 'nm'", {'HOH': 'tip3p'}, radius=5.0, unit='nm')
        raises('the radius must be positive', {'HOH': 'tip3p'}, radius=-1.0)
        raises('inner radius exceeds the outer', shells, shells=(5.0, 4.0))
        raises('no force field for residue HOH 1', {})
        raises('listed twice', {'HOH': 'tip3p', 'hoh': 'ahlstrom'})
        raises('no solvent force fields', {}, shells=(5.0, 9.0))
        raises('assigned to the polarizable region',
               {'polarizable': {
                   'HOH': 'tip3p'
               }},
               shells=(5.0, 9.0))
        raises("given as {'polarizable'", {'HOH': 'tip3p'}, shells=(5.0, 9.0))

        # A force field of another atom count (the solute named otherwise):
        # a three-atom MOL from a user force-field directory.
        (tmp_path / 'three.fff').write_text(
            'fika-force-field 1 three nonpolarizable\n' +
            'residue MOL 3\ncharges\n1 0\n2 0\n3 0\nend\n')
        monkeypatch.setenv('VLXFORCEFIELDPATH', str(tmp_path))
        raises('residue MOL 1 of chain \'A\' has 8 atoms, force field three 3', {
            'HOH': 'tip3p',
            'MOL': 'three'
        },
               solute='XYZ')

        waters_only = reader.parse(get_pdb(get_waters(3), solute=False))
        with pytest.raises(ValueError, match='no MOL residue'):
            reader.get_classical_system(waters_only, {'HOH': 'tip3p'},
                                        radius=5.0)
        twice = reader.parse(
            get_pdb(get_waters(12)).replace('W  12', 'Z  12').replace(
                'HOH W  11', 'MOL W  11'))
        with pytest.raises(ValueError, match='more than one MOL residue'):
            reader.get_classical_system(twice, {'HOH': 'tip3p'})

    def test_embedding_from_pdb(self):

        reader = FikaPdbReader()
        waters = get_waters(30)
        residues = reader.parse(get_pdb(waters))
        system, _ = reader.get_classical_system(residues, {'HOH': 'ahlstrom'})

        # The same system built directly.
        direct = FikaClassicalSystem()
        direct.add_force_field(
            FikaForceField('ahlstrom', 'HOH', [-0.669, 0.3345, 0.3345],
                           {0: 9.718}))
        direct.add_residues('HOH', 'ahlstrom', [8, 1, 1], rounded_bohr(waters),
                            True)
        a = FikaInducedDipolesDriver().compute(system)
        b = FikaInducedDipolesDriver().compute(direct)
        assert a.converged
        assert np.array_equal(a.dipoles, b.dipoles)

        molecule = reader.get_solute(residues)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)
        n = basis.get_dimensions_of_basis()
        density = np.eye(n) * 0.1
        result = FikaQmmmEmbeddingDriver().compute(molecule, basis, density,
                                                   system)
        assert result.induced.converged
        assert np.isfinite(result.energy())
        assert result.fock().shape == (n, n)
