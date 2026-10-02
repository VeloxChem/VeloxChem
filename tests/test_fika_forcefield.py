import numpy as np
import pytest

from veloxchem.fikaforcefield import FikaForceFieldReader
from veloxchem.veloxchemlib import (FikaClassicalSystem, FikaForceField,
                                    FikaInducedDipolesDriver,
                                    fika_classical_charges)

WATER = """
fika-force-field 1 Model polarizable   # a comment
# a two-residue force field

residue HOH 3
polarizabilities isotropic 2
1 9.5
3 1.25
charges
1 -0.8
2  0.4
3  0.4
end

residue ION 1
charges
1 1.0
polarizabilities isotropic 1
1 0.5
end
"""


def assert_error(text, line, message):

    with pytest.raises(ValueError) as error:
        FikaForceFieldReader().parse(text, 'test.fff')
    assert str(error.value).startswith(f'test.fff: line {line}: ')
    assert message in str(error.value)


class TestFikaForceFieldReader:

    def test_parse(self):

        force_fields = FikaForceFieldReader().parse(WATER)

        assert [f.residue for f in force_fields] == ['hoh', 'ion']
        water, ion = force_fields
        assert water.label == 'model'
        assert water.atom_count == 3
        assert water.polarizable
        assert water.charges == [-0.8, 0.4, 0.4]
        assert water.polarizabilities == {0: 9.5, 2: 1.25}
        assert ion.charges == [1.0]
        assert ion.polarizabilities == {0: 0.5}

    def test_errors(self):

        assert_error('fika-force-field 2 x nonpolarizable\n', 1,
                     'unsupported format version')
        assert_error('fika-forcefield 1 x nonpolarizable\n', 1,
                     "expected 'fika-force-field' header")
        assert_error('fika-force-field 1 x maybe\n', 1,
                     "expected 'polarizable' or 'nonpolarizable'")

        header = 'fika-force-field 1 x nonpolarizable\n'
        assert_error(header + 'atom HOH 3\n', 2, "expected 'residue'")
        assert_error(header + 'residue HOH 0\n', 2, 'at least one atom')
        assert_error(header + 'residue HOH 2\ncharges\n2 0.5\n1 -0.5\nend\n',
                     4, 'charges must list atoms 1..2 in order')
        assert_error(header + 'residue HOH 2\ncharges\n1 0.5\n3 -0.5\nend\n',
                     5, 'atom index 3 outside 1..2')
        assert_error(header + 'residue HOH 1\ncharges\n1 abc\nend\n', 4,
                     "invalid number 'abc'")
        assert_error(header + 'residue HOH 1\ncharges\n1 0.5\ncharges\n', 5,
                     "repeated section 'charges'")
        assert_error(header + 'residue HOH 1\nend\n', 3, 'has no charges')
        assert_error(header + 'residue HOH 1\ncharges\n1 0.5\nmass\n', 5,
                     "unknown section 'mass'")
        assert_error(
            header + 'residue HOH 1\ncharges\n1 0.5\nend\n' +
            'residue HOH 1\ncharges\n1 0.5\nend\n', 6, 'repeated residue')

        # Polarizabilities must match the force-field flag.
        assert_error(
            header + 'residue HOH 1\ncharges\n1 0.5\n' +
            'polarizabilities isotropic 1\n1 2.0\nend\n', 2,
            'a nonpolarizable force field has no polarizabilities')
        assert_error(
            'fika-force-field 1 x polarizable\nresidue HOH 1\n' +
            'charges\n1 0.5\nend\n', 2,
            'a polarizable force field needs polarizabilities')
        polarizable = 'fika-force-field 1 x polarizable\nresidue HOH 2\n'
        assert_error(
            polarizable + 'polarizabilities isotropic 2\n2 1.0\n1 1.0\n', 5,
            'increasing order')
        assert_error(polarizable + 'polarizabilities isotropic 1\n1 -1.0\n' +
                     'charges\n1 0\n2 0\nend\n', 2, 'residue HOH')

        # Parts of the format that are not supported yet.
        assert_error(polarizable + 'polarizabilities anisotropic 1\n', 3,
                     'anisotropic polarizabilities are not supported yet')
        for section in ['dipoles 1', 'quadrupoles 1', 'geometry bohr']:
            assert_error(header + f'residue HOH 1\n{section}\n', 3,
                         'is not supported yet')

        with pytest.raises(ValueError, match='unexpected end of file'):
            FikaForceFieldReader().parse(header + 'residue HOH 1\ncharges\n')
        with pytest.raises(ValueError, match='has no residues'):
            FikaForceFieldReader().parse(header)

    def test_library(self):

        reader = FikaForceFieldReader()

        (tip3p,) = reader.load('TIP3P')
        assert tip3p.label == 'tip3p'
        assert tip3p.charges == [-0.834, 0.417, 0.417]
        assert not tip3p.polarizable

        ahlstrom = reader.load('ahlstrom', 'HOH')
        assert ahlstrom.charges == [-0.669, 0.3345, 0.3345]
        assert ahlstrom.polarizabilities == {0: 9.718}

        with pytest.raises(FileNotFoundError):
            reader.load('no-such-force-field')
        with pytest.raises(ValueError, match='has no residue'):
            reader.load('tip3p', 'SOL')

    def test_search_path(self, tmp_path, monkeypatch):

        # A user directory takes precedence over the library.
        (tmp_path / 'tip3p.fff').write_text(
            'fika-force-field 1 tip3p nonpolarizable\n' +
            'residue HOH 3\ncharges\n1 -0.8\n2 0.4\n3 0.4\nend\n')
        monkeypatch.setenv('VLXFORCEFIELDPATH', f'/no/such/directory:{tmp_path}')

        reader = FikaForceFieldReader()
        assert reader.get_search_path()[1] == tmp_path
        assert reader.load('tip3p', 'hoh').charges == [-0.8, 0.4, 0.4]
        assert reader.load('ahlstrom', 'hoh').charges[0] == -0.669

        with pytest.raises(FileNotFoundError):
            reader.read(tmp_path / 'missing.fff')

    def test_read_into_classical_system(self, tmp_path):

        # Force fields read from a file act as those built in Python.
        path = tmp_path / 'model.fff'
        path.write_text(WATER)
        read_fields = FikaForceFieldReader().read(path)
        built_fields = [
            FikaForceField('model', 'HOH', [-0.8, 0.4, 0.4], {
                0: 9.5,
                2: 1.25
            }),
            FikaForceField('model', 'ION', [1.0], {0: 0.5}),
        ]

        rng = np.random.default_rng(3)
        local = np.array([[0.0, 0.0, 0.0], [1.43, 0.0, 1.1], [-1.43, 0.0,
                                                              1.1]])
        waters = np.array([local + rng.uniform(-12.0, 12.0, 3) for _ in range(6)])
        ions = rng.uniform(15.0, 20.0, (2, 3))

        results = []
        for force_fields in [read_fields, built_fields]:
            system = FikaClassicalSystem()
            for force_field in force_fields:
                system.add_force_field(force_field)
            system.add_residues('HOH', 'model', [8, 1, 1], waters, True)
            system.add_residues('ION', 'model', [11], ions, True)
            results.append(FikaInducedDipolesDriver().compute(system, None))
            charges = fika_classical_charges(system)[0]
            assert charges == [-0.8, 0.4, 0.4] * 6 + [1.0, 1.0]

        assert results[0].converged
        assert np.array_equal(results[0].dipoles, results[1].dipoles)
