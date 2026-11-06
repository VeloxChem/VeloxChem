from types import SimpleNamespace
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import (MolecularBasis, get_basis_set_category,
                                      _known_basis_set_categories)
from veloxchem.environment import get_basis_path
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
import veloxchem.scfdriver as scfdriver


@pytest.mark.solvers
class TestRiJkAtomLimit:
    """RI-JK falls back to the four-center integrals from 128 atoms on in a
    basis set of at most double-zeta valence without diffuse functions, from
    256 atoms on in a triple-zeta one without diffuse functions, and from 512
    atoms on in any basis set; force_ri keeps RI-JK, and RI-J is not
    limited."""

    @staticmethod
    def is_orbital_basis(fname):

        return not (fname.endswith('-RIFIT') or fname in [
            'DEF2-UNIVERSAL-JFIT', 'DEF2-UNIVERSAL-JKFIT', 'DEF2-ECP',
            'AO-START-GUESS', 'AO-START-GUESS-FOR-ECP'
        ])

    def test_every_orbital_basis_set_has_a_category(self):

        categories = _known_basis_set_categories()

        fnames = [
            path.name for path in get_basis_path().iterdir() if path.is_file()
        ]
        orbital = sorted(f for f in fnames if self.is_orbital_basis(f))

        assert orbital
        assert sorted(categories) == orbital

    @pytest.mark.parametrize('name, category', [
        ('def2-svp', (2, False)),
        ('DEF2-SV(P)', (2, False)),
        ('6-31G*', (2, False)),
        ('6-31+G**', (2, True)),
        ('def2-svpd', (2, True)),
        ('sto-3g', (1, False)),
        ('def2-tzvp', (3, False)),
        ('6-311G**', (3, False)),
        ('def2-tzvpd', (3, True)),
        ('aug-cc-pvtz', (3, True)),
        ('def2-qzvp', (4, False)),
        ('my-own-basis', None),
    ])
    def test_category(self, name, category):

        assert get_basis_set_category(name) == category

    @pytest.mark.parametrize('label, limit', [
        ('DEF2-SVP', 128),
        ('STO-3G', 128),
        ('6-31G*', 128),
        ('DEF2-TZVP', 256),
        ('DEF2-SVPD', 512),
        ('DEF2-TZVPD', 512),
        ('DEF2-QZVP', 512),
        ('MY-OWN-BASIS', 512),
    ])
    def test_limit_at_its_boundary(self, label, limit):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()

        basis = SimpleNamespace(get_main_basis_label=lambda: label)

        assert scf_drv._get_ri_jk_atom_limit(basis)[0] == limit

        for natoms, falls_back in [(limit - 1, False), (limit, True)]:
            molecule = SimpleNamespace(number_of_atoms=lambda: natoms)
            scf_drv.force_ri = False
            assert (scf_drv._ri_jk_exceeds_atom_limit(molecule, basis) ==
                    falls_back)
            scf_drv.force_ri = True
            assert not scf_drv._ri_jk_exceeds_atom_limit(molecule, basis)

    @staticmethod
    def run_scf(basis_label, xcfun_label, force_ri=False, ri=True):

        molecule = Molecule.read_str(
            """
            O 0.0000000000    0.0000000000   -0.0254395383
            H 0.0000000000    0.7695699584    0.5948147012
            H 0.0000000000   -0.7695699584    0.5948147012
            """,
            units='angstrom')
        basis = MolecularBasis.read(molecule, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.xcfun = xcfun_label
        scf_drv.force_ri = force_ri
        scf_drv.ostream.mute()

        warnings = []
        scf_drv.ostream.print_warning = warnings.append

        scf_drv.compute(molecule, basis)

        return scf_drv, molecule, basis, warnings

    @pytest.mark.parametrize('basis_label, limit_name', [
        ('def2-svp', '_RI_JK_FALLBACK_ATOMS_DOUBLE_ZETA'),
        ('def2-tzvp', '_RI_JK_FALLBACK_ATOMS_TRIPLE_ZETA'),
        ('def2-svpd', '_RI_JK_FALLBACK_ATOMS'),
    ])
    def test_scf_and_what_follows_fall_back(self, monkeypatch, basis_label,
                                            limit_name):

        # NOTE: the limit lowered to the 3 atoms of water.
        monkeypatch.setattr(scfdriver, limit_name, 3)

        scf_drv, molecule, basis, warnings = self.run_scf(basis_label, 'b3lyp')

        assert scf_drv.ri_jk
        assert scf_drv._ri_fallback
        assert not scf_drv._uses_ri_jk()
        assert scf_drv._ri_drv is None
        assert any('The molecule has 3 atoms' in w for w in warnings)

        if scf_drv.rank == mpi_master():
            assert scf_drv.scf_results['ri_fallback']

        rsp_drv = TdaEigenSolver()
        rsp_drv.nstates = 3
        rsp_drv.ostream.mute()
        rsp_drv.compute(molecule, basis, scf_drv.scf_results)
        assert not rsp_drv._uses_ri_jk()

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.ostream.mute()
        grad_drv.compute(molecule, basis, scf_drv.scf_results)
        assert not grad_drv._ri_fallback

        # NOTE: the energy is the four-center one.
        four_center = self.run_scf(basis_label, 'b3lyp', ri=False)[0]
        assert not four_center.ri_jk

        if scf_drv.rank == mpi_master():
            assert abs(scf_drv.get_scf_energy() -
                       four_center.get_scf_energy()) < 1.0e-8

    def test_force_ri_keeps_ri_jk(self, monkeypatch):

        monkeypatch.setattr(scfdriver, '_RI_JK_FALLBACK_ATOMS_DOUBLE_ZETA', 3)

        scf_drv, molecule, basis, warnings = self.run_scf(
            'def2-svp', 'b3lyp', True)

        assert scf_drv._uses_ri_jk()
        assert any('force_ri is set' in w for w in warnings)

        rsp_drv = TdaEigenSolver()
        rsp_drv.nstates = 3
        rsp_drv.ostream.mute()
        rsp_drv.compute(molecule, basis, scf_drv.scf_results)
        assert rsp_drv._uses_ri_jk()

    def test_ri_j_is_not_limited(self, monkeypatch):

        for name in [
                '_RI_JK_FALLBACK_ATOMS_DOUBLE_ZETA',
                '_RI_JK_FALLBACK_ATOMS_TRIPLE_ZETA', '_RI_JK_FALLBACK_ATOMS'
        ]:
            monkeypatch.setattr(scfdriver, name, 3)

        scf_drv, _, _, warnings = self.run_scf('def2-svp', 'blyp')

        assert scf_drv._uses_ri_coulomb()
        assert not scf_drv._ri_fallback
        assert not any('The molecule has' in w for w in warnings)
