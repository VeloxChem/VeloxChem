from mpi4py import MPI
from pathlib import Path
import h5py
import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.fikaembedding import (FikaEmbedding,
                                     fika_embedding_sanity_check)
from veloxchem.oneeints import compute_electric_field_integrals
from veloxchem.oneeints import compute_nuclear_potential_integrals
from veloxchem.veloxchemlib import mpi_master, bohr_in_angstroms
from veloxchem.veloxchemlib import (FikaChargeSummation, FikaClassicalSystem,
                                    FikaDipolePotentialDriver, FikaForceField,
                                    FikaInducedDipoleFockDriver,
                                    FikaInducedDipolesDriver,
                                    FikaInnerSwitch,
                                    FikaNuclearAttractionDriver,
                                    FikaQmmmEmbeddingDriver,
                                    FikaQmmmEmbeddingOptions, FikaQmmmSources,
                                    FikaQuadrupolePotentialDriver,
                                    FikaResidue, FikaTholeDamping,
                                    fika_classical_charges,
                                    fika_polarizable_sites)

# Ahlstrom polarizable water: charges on O, H, H and an isotropic
# polarizability on O (a.u.).
AHLSTROM_CHARGES = [-0.669, 0.3345, 0.3345]
AHLSTROM_ALPHA = 9.718


def get_acrolein():

    xyz_string = """8
        C3OH4
        C             -0.145335   -0.546770    0.000607
        C              1.274009   -0.912471   -0.000167
        C              1.630116   -2.207690   -0.000132
        O             -0.560104    0.608977    0.000534
        H             -0.871904   -1.386459    0.001253
        H              2.004448   -0.101417   -0.000710
        H              0.879028   -3.000685    0.000484
        H              2.675323   -2.516779   -0.000673
    """
    return Molecule.read_xyz_string(xyz_string)


def get_water_shell(molecule, count, seed=7):
    """
    Water molecules (O, H, H; bohr) placed at random around the molecule,
    at least 5 bohr from its atoms and 5.3 bohr apart (oxygens), with random
    orientations; deterministic for a given seed.
    """

    rng = np.random.default_rng(seed)
    atoms = molecule.get_coordinates_in_bohr()
    centre = atoms.mean(axis=0)

    bond, angle = 1.8088, np.radians(104.52)
    local = np.array([[0.0, 0.0, 0.0],
                      [bond * np.sin(angle / 2), 0.0, bond * np.cos(angle / 2)],
                      [-bond * np.sin(angle / 2), 0.0, bond * np.cos(angle / 2)]])

    waters = []
    while len(waters) < count:
        oxygen = centre + rng.uniform(-14.0, 14.0, 3)
        if np.min(np.linalg.norm(atoms - oxygen, axis=1)) < 5.0:
            continue
        if waters and np.min(
                np.linalg.norm(np.array([w[0] for w in waters]) - oxygen,
                               axis=1)) < 5.3:
            continue
        rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        waters.append(oxygen + local @ rotation.T)

    return np.array(waters)


def get_classical_system(waters):

    system = FikaClassicalSystem()
    system.add_force_field(
        FikaForceField('ahlstrom', 'HOH', AHLSTROM_CHARGES,
                       {0: AHLSTROM_ALPHA}))
    system.add_residues('HOH', 'ahlstrom', [8, 1, 1], waters, True)
    return system


def get_hf_density(molecule, basis):

    scf_drv = ScfRestrictedDriver()
    scf_drv.conv_thresh = 1.0e-8
    scf_drv.ostream.mute()
    results = scf_drv.compute(molecule, basis)
    density = None
    if scf_drv.rank == mpi_master():
        density = results['D_alpha'] + results['D_beta']
    return scf_drv.comm.bcast(density, root=mpi_master())


class TestFikaIntegrals:

    @pytest.mark.parametrize('basis_label',
                             ['def2-svp', 'def2-tzvp', 'cc-pvdz'])
    def test_point_charges_and_dipoles(self, basis_label):

        # fika's potential integrals of point charges and dipoles (no
        # electron charge) are the negatives of VeloxChem's.

        molecule = get_acrolein()
        basis = MolecularBasis.read(molecule, basis_label, ostream=None)

        rng = np.random.default_rng(11)
        centre = molecule.get_coordinates_in_bohr().mean(axis=0)
        coordinates = np.vstack([
            centre + rng.uniform(-8.0, 8.0, (20, 3)),
            centre + rng.normal(size=(20, 3)) * 40.0,
        ])
        charges = rng.uniform(-1.0, 1.0, len(coordinates))
        dipoles = rng.uniform(-1.0, 1.0, (len(coordinates), 3))

        direct = FikaChargeSummation.direct

        ref = -compute_nuclear_potential_integrals(molecule, basis, charges,
                                                   coordinates)
        fika = FikaNuclearAttractionDriver(0, direct).compute(
            molecule, basis, charges, coordinates, 0.0)
        assert np.max(np.abs(fika - ref)) < 1.0e-12 * np.max(np.abs(ref))

        ref = -compute_electric_field_integrals(molecule, basis, coordinates,
                                                dipoles)
        fika = FikaDipolePotentialDriver(0, direct).compute(
            molecule, basis, dipoles, coordinates, 0.0)
        assert np.max(np.abs(fika - ref)) < 1.0e-12 * np.max(np.abs(ref))

    def test_quadrupoles(self):

        # A traceless quadrupole eps (mu b^T + b mu^T), mu perpendicular to
        # b, is half the difference of the dipoles mu at C + eps b and C - eps
        # b (Richardson extrapolation over eps removes the eps^2 term).

        molecule = get_acrolein()
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        rng = np.random.default_rng(5)
        centre = molecule.get_coordinates_in_bohr().mean(axis=0)
        coordinates = centre + rng.uniform(-8.0, 8.0, (8, 3))
        b = rng.normal(size=(8, 3))
        b /= np.linalg.norm(b, axis=1)[:, None]
        mu = rng.normal(size=(8, 3))
        mu -= np.sum(mu * b, axis=1)[:, None] * b

        direct = FikaChargeSummation.direct
        dipole_drv = FikaDipolePotentialDriver(0, direct)

        def from_dipoles(eps):
            plus = dipole_drv.compute(molecule, basis, mu, coordinates + eps * b,
                                      0.0)
            minus = dipole_drv.compute(molecule, basis, mu,
                                       coordinates - eps * b, 0.0)
            return 0.5 * (plus - minus) / eps

        ref = (4.0 * from_dipoles(5.0e-4) - from_dipoles(1.0e-3)) / 3.0

        pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
        quadrupoles = np.array([[m[i] * v[j] + m[j] * v[i]
                                 for i, j in pairs]
                                for m, v in zip(mu, b)])
        fika = FikaQuadrupolePotentialDriver(0, direct).compute(
            molecule, basis, quadrupoles, coordinates, 0.0)

        assert np.max(np.abs(fika - ref)) < 1.0e-9 * np.max(np.abs(ref))

    def test_multipole_summation(self):

        molecule = get_acrolein()
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        rng = np.random.default_rng(9)
        centre = molecule.get_coordinates_in_bohr().mean(axis=0)
        coordinates = centre + rng.normal(size=(2000, 3)) * 30.0
        charges = rng.uniform(-1.0, 1.0, len(coordinates))

        direct = FikaNuclearAttractionDriver(
            0, FikaChargeSummation.direct).compute(molecule, basis, charges,
                                                   coordinates)
        multipole = FikaNuclearAttractionDriver(
            0, FikaChargeSummation.multipole).compute(molecule, basis, charges,
                                                      coordinates)

        assert np.max(np.abs(multipole - direct)) < 1.0e-11

        # All-zero dipoles or quadrupoles give a zero matrix.
        summation = FikaChargeSummation.multipole
        zero = FikaDipolePotentialDriver(0, summation).compute(
            molecule, basis, np.zeros((len(coordinates), 3)), coordinates)
        assert np.max(np.abs(zero)) == 0.0
        zero = FikaQuadrupolePotentialDriver(0, summation).compute(
            molecule, basis, np.zeros((len(coordinates), 6)), coordinates)
        assert np.max(np.abs(zero)) == 0.0

    def test_invalid_input(self):

        molecule = get_acrolein()
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        with pytest.raises(ValueError):
            FikaNuclearAttractionDriver().compute(molecule, basis, [1.0, 2.0],
                                                  [[5.0, 0.0, 0.0]])

        ghost = Molecule.read_molecule_string(
            'O 0 0 0\nBq_H 0 0.76 0.59\nH 0 -0.76 0.59')
        ghost_basis = MolecularBasis.read(ghost, 'def2-svp', ostream=None)
        with pytest.raises(ValueError):
            FikaNuclearAttractionDriver().compute(ghost, ghost_basis, [1.0],
                                                  [[5.0, 0.0, 0.0]])


class TestFikaClassicalSystem:

    def test_system(self):

        waters = get_water_shell(get_acrolein(), 10)
        system = get_classical_system(waters)

        force_field = system.force_field('ahlstrom', 'hoh')
        assert force_field.atom_count == 3
        assert force_field.polarizable
        assert force_field.charges == AHLSTROM_CHARGES
        assert force_field.polarizabilities == {0: AHLSTROM_ALPHA}

        assert system.number_of_residues(True) == 10
        assert system.number_of_residues(False) == 0

        positions, alphas, owners = fika_polarizable_sites(system)
        assert np.array_equal(positions, waters[:, 0, :])
        assert alphas == [AHLSTROM_ALPHA] * 10
        assert owners == list(range(10))

        charges, coordinates, offsets = fika_classical_charges(system)
        assert charges == AHLSTROM_CHARGES * 10
        assert np.array_equal(coordinates, waters.reshape(-1, 3))
        assert offsets == list(range(0, 31, 3))

        residue = FikaResidue('HOH', 0, 'ahlstrom', [8, 1, 1], waters[0], True)
        assert residue.name == 'hoh'
        assert residue.elements == [8, 1, 1]
        assert np.array_equal(residue.coordinates, waters[0])

    def test_errors(self):

        # A residue without its force field in the system.
        system = FikaClassicalSystem()
        system.add_residue(
            FikaResidue('HOH', 0, 'ahlstrom', [8, 1, 1],
                        get_water_shell(get_acrolein(), 1)[0], True))
        with pytest.raises(RuntimeError):
            fika_classical_charges(system)

        # Two undamped sites 1 bohr apart (polarizability 10) in the field of
        # an off-centre charge: B = alpha^-1 - T is not positive definite.
        system = FikaClassicalSystem()
        system.add_force_field(FikaForceField('x', 'x', [0.5], {0: 10.0}))
        system.add_force_field(FikaForceField('c', 'c', [1.0]))
        system.add_residues('x', 'x', [8], np.array([[0, 0, 0], [0, 0, 1.0]]),
                            True)
        system.add_residues('c', 'c', [11], np.array([[0, 0, 6.0]]), False)
        with pytest.raises(RuntimeError):
            FikaInducedDipolesDriver().compute(system, None)

        # Thole damping keeps it positive definite.
        result = FikaInducedDipolesDriver().compute(system, FikaTholeDamping())
        assert result.converged


class TestFikaSwitching:

    @staticmethod
    def _switch():

        return FikaInnerSwitch(1.0 / bohr_in_angstroms(),
                               1.9 / bohr_in_angstroms())

    @staticmethod
    def _charge_system(r_ang, alpha=2.0, q=0.5):
        # polarizable site at the origin, a charge of another residue at +x

        r = r_ang / bohr_in_angstroms()
        system = FikaClassicalSystem()
        system.add_force_field(FikaForceField('p', 'p', [0.0], {0: alpha}))
        system.add_force_field(FikaForceField('c', 'c', [q]))
        system.add_residue(
            FikaResidue('p', 0, 'p', [6], np.array([[0.0, 0.0, 0.0]]), True))
        system.add_residue(
            FikaResidue('c', 0, 'c', [6], np.array([[r, 0.0, 0.0]]), False))
        return system

    @staticmethod
    def _two_site_system(alpha1, alpha2, r12_ang, q, rq_ang):

        r12 = r12_ang / bohr_in_angstroms()
        rq = rq_ang / bohr_in_angstroms()
        system = FikaClassicalSystem()
        system.add_force_field(FikaForceField('p', 'p', [0.0], {0: alpha1}))
        system.add_force_field(FikaForceField('q', 'q', [0.0], {0: alpha2}))
        system.add_force_field(FikaForceField('c', 'c', [q]))
        system.add_residue(
            FikaResidue('p', 0, 'p', [6], np.array([[0.0, 0.0, 0.0]]), True))
        system.add_residue(
            FikaResidue('q', 0, 'q', [6], np.array([[r12, 0.0, 0.0]]),
                        True))
        system.add_residue(
            FikaResidue('c', 0, 'c', [6], np.array([[0.0, 0.0, rq]]), False))
        return system

    def test_charge_dipole_switch(self):

        alpha, q = 2.0, 0.5
        switch = self._switch()

        for r_ang in (1.0, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0):
            r = r_ang / bohr_in_angstroms()
            if r_ang <= 1.0:
                s = 0.0
            elif r_ang >= 1.9:
                s = 1.0
            else:
                x = (r_ang - 1.0) / 0.9
                s = x**3 * (10.0 + x * (-15.0 + 6.0 * x))
            result = FikaInducedDipolesDriver().compute(
                self._charge_system(r_ang), None, switch)
            # the fixed charge field at the origin is -q/r^2 along x,
            # scaled by S(r); mu = alpha F for a single site
            assert result.dipoles[0, 0] == pytest.approx(
                -alpha * q * s / r**2, abs=1.0e-10)
            assert result.dipoles[0, 1] == 0.0
            assert result.dipoles[0, 2] == 0.0

    def test_switch_none_is_undamped_field(self):

        r = 1.0 / bohr_in_angstroms()
        result = FikaInducedDipolesDriver().compute(self._charge_system(1.0),
                                                    None)
        assert result.dipoles[0, 0] == pytest.approx(-2.0 * 0.5 / r**2,
                                                     abs=1.0e-10)

    def test_invalid_switch(self):

        for r_on, r_off in [(-1.0, 1.8), (1.8, 1.2), (1.2, 1.2)]:
            switch = FikaInnerSwitch(r_on / bohr_in_angstroms(),
                                     r_off / bohr_in_angstroms())
            with pytest.raises(ValueError):
                FikaInducedDipolesDriver().compute(
                    self._charge_system(2.0), None, switch)

    def test_overlapping_charge_auto_excluded(self):

        # a permanent charge of another residue at the site position is a
        # degenerate overlap; S(0) = 0 auto-excludes it, like the FMM near
        # field excludes a source at the target point itself
        system = FikaClassicalSystem()
        system.add_force_field(FikaForceField('p', 'p', [0.0], {0: 2.0}))
        system.add_force_field(FikaForceField('c', 'c', [0.5]))
        system.add_residue(
            FikaResidue('p', 0, 'p', [6], np.array([[0.0, 0.0, 0.0]]), True))
        system.add_residue(
            FikaResidue('c', 0, 'c', [6], np.array([[0.0, 0.0, 0.0]]), False))

        direct = FikaInducedDipolesDriver().compute(system, None,
                                                    self._switch())
        assert direct.converged
        assert np.array_equal(direct.dipoles, np.zeros((1, 3)))

        multipole_drv = FikaInducedDipolesDriver()
        multipole_drv.options.summation = FikaChargeSummation.multipole
        multipole = multipole_drv.compute(system, None, self._switch())
        assert multipole.converged
        assert np.array_equal(multipole.dipoles, np.zeros((1, 3)))

    def test_multipole_near_field_switch(self):

        # A spread system whose FMM near radius covers r_off: the switch is
        # applied in the near-field P2P sums themselves, so a different-
        # residue charge at ~1e-12 A gives no spurious field from cancellation
        # (the (S - 1) correction of the far-field path would leave
        # ~eps * q / r^2, i.e. ~1e9 a.u. here).
        system = FikaClassicalSystem()
        system.add_force_field(FikaForceField('p', 'p', [0.0], {0: 2.0}))
        system.add_force_field(FikaForceField('c', 'c', [0.7]))
        system.add_residue(
            FikaResidue('p', 0, 'p', [6],
                        np.array([[0.1 + 1.0e-12, 0.0, 0.0]]), True))
        for x, y, z in [(0.1, 0.0, 0.0), (100.0, 0.0, 0.0),
                        (-100.0, 0.0, 0.0), (0.0, 100.0, 0.0),
                        (0.0, -100.0, 0.0), (0.0, 0.0, 100.0),
                        (0.0, 0.0, -100.0)]:
            system.add_residue(
                FikaResidue('c', 0, 'c', [8],
                            np.array([[x, y, z]]), False))

        direct_drv = FikaInducedDipolesDriver()
        direct_drv.options.summation = FikaChargeSummation.direct
        direct = direct_drv.compute(system, None, self._switch())

        multipole_drv = FikaInducedDipolesDriver()
        multipole_drv.options.summation = FikaChargeSummation.multipole
        multipole = multipole_drv.compute(system, None, self._switch())

        assert direct.converged and multipole.converged
        assert np.max(np.abs(direct.dipoles - multipole.dipoles)) < 1.0e-6

    def test_multipole_near_field_switch_dipoles(self):

        # Spread polarizable sites whose FMM near radius covers r_off. Part 1 checks that a pair
        # at an extreme overlap is switched off by the weighted near field (the correction path
        # would leave ~eps / r^3 spurious coupling); part 2 checks the S (lambda - 1) damping
        # part of the weighted path against the direct sum.

        switch = self._switch()
        r_on = 1.0 / bohr_in_angstroms()
        r_off = 1.9 / bohr_in_angstroms()
        alphas = [2.0, 3.0, 1.0, 1.0, 1.0]
        far = [[10.0, 0.0, 0.0], [-10.0, 0.0, 0.0], [0.0, 10.0, 0.0]]
        charges = np.array([[0.0, 0.0, 8.0], [0.0, 0.0, -4.0]])
        q = 0.5

        def build(r12):

            positions = np.vstack(([[0.0, 0.0, 0.0], [r12, 0.0, 0.0]], far))
            system = FikaClassicalSystem()
            for i, alpha in enumerate(alphas):
                system.add_force_field(
                    FikaForceField(f'p{i}', f'p{i}', [0.0], {0: alpha}))
                system.add_residue(
                    FikaResidue(f'p{i}', 0, f'p{i}', [6],
                                positions[i].reshape(1, 3), True))
            system.add_force_field(FikaForceField('c', 'c', [q]))
            for pos in charges:
                system.add_residue(
                    FikaResidue('c', 0, 'c', [8], pos.reshape(1, 3), False))
            return system, positions

        # part 1: extreme overlap, no damping, switched reference
        system, positions = build(1.0e-6)
        n = len(alphas)
        field = np.zeros((n, 3))
        for pos in charges:
            for i in range(n):
                d = positions[i] - pos
                field[i] += q * d / np.linalg.norm(d)**3
        b = np.zeros((3 * n, 3 * n))
        for i in range(n):
            b[3 * i:3 * i + 3, 3 * i:3 * i + 3] = np.eye(3) / alphas[i]
        for i in range(n):
            for j in range(i + 1, n):
                d = positions[i] - positions[j]
                r = np.linalg.norm(d)
                if r <= r_on:
                    s = 0.0
                elif r >= r_off:
                    s = 1.0
                else:
                    x = (r - r_on) / (r_off - r_on)
                    s = x**3 * (10.0 + x * (-15.0 + 6.0 * x))
                t = s * (3.0 * np.outer(d, d) - r * r * np.eye(3)) / r**5
                b[3 * i:3 * i + 3, 3 * j:3 * j + 3] = -t
                b[3 * j:3 * j + 3, 3 * i:3 * i + 3] = -t
        ref = np.linalg.solve(b, field.reshape(-1)).reshape(-1, 3)

        multipole_drv = FikaInducedDipolesDriver()
        multipole_drv.options.summation = FikaChargeSummation.multipole
        multipole = multipole_drv.compute(system, None, switch)
        assert multipole.converged
        assert np.max(np.abs(multipole.dipoles - ref)) < 1.0e-7

        # part 2: ramp pair, Thole damping on, direct against multipole
        system, _ = build(2.8)
        direct_drv = FikaInducedDipolesDriver()
        direct_drv.options.summation = FikaChargeSummation.direct
        direct = direct_drv.compute(system, FikaTholeDamping(), switch)
        multipole = multipole_drv.compute(system, FikaTholeDamping(), switch)
        assert direct.converged and multipole.converged
        assert np.max(np.abs(direct.dipoles - multipole.dipoles)) < 1.0e-7

    def test_dipole_dipole_switch(self):

        # two polarizable sites 1.4 A apart (their T coupling is switched)
        # plus a distant charge; compare with a direct solve of
        # (alpha^-1 - S(r) T) mu = F

        alpha1, alpha2, q = 2.0, 3.0, 0.8
        r12_ang, rq_ang = 1.4, 4.0
        r12 = r12_ang / bohr_in_angstroms()
        rq = rq_ang / bohr_in_angstroms()
        system = self._two_site_system(alpha1, alpha2, r12_ang, q, rq_ang)
        dipole_driver = FikaInducedDipolesDriver()
        dipole_driver.options.tolerance = 1.0e-12
        result = dipole_driver.compute(system, None, self._switch())

        pos = np.array([[0.0, 0.0, 0.0], [r12, 0.0, 0.0]])
        field = np.zeros((2, 3))
        for i in range(2):
            d = pos[i] - np.array([0.0, 0.0, rq])
            field[i] = q * d / np.linalg.norm(d)**3
        b = np.zeros((6, 6))
        b[0:3, 0:3] = np.eye(3) / alpha1
        b[3:6, 3:6] = np.eye(3) / alpha2
        d = pos[0] - pos[1]
        r = np.linalg.norm(d)
        x = (r12_ang - 1.0) / 0.9
        s = x**3 * (10.0 + x * (-15.0 + 6.0 * x))
        t = s * (3.0 * np.outer(d, d) - r * r * np.eye(3)) / r**5
        b[0:3, 3:6] = -t
        b[3:6, 0:3] = -t
        ref = np.linalg.solve(b, field.reshape(-1)).reshape(-1, 3)

        assert np.max(np.abs(result.dipoles - ref)) < 1.0e-8
        # without the switch the coupling (and the dipoles) differ
        bare = FikaInducedDipolesDriver().compute(system, None)
        assert np.max(np.abs(bare.dipoles - ref)) > 1.0e-6

    def test_direct_against_multipole(self):

        # direct and multipole summations must agree with the switch (and
        # with Thole damping on top of it)

        rng = np.random.default_rng(11)
        system = FikaClassicalSystem()
        for i in range(4):
            system.add_force_field(
                FikaForceField(f'p{i}', f'p{i}', [0.0], {0: 2.0 + i}))
            system.add_residue(
                FikaResidue(f'p{i}', 0, f'p{i}', [6],
                            rng.normal(size=(1, 3)) * 1.0, True))
        for i in range(6):
            system.add_force_field(
                FikaForceField(f'c{i}', f'c{i}',
                               [(-1.0)**i * 0.4]))
            system.add_residue(
                FikaResidue(f'c{i}', 0, f'c{i}', [8],
                            rng.normal(size=(1, 3)) * 1.0, False))

        direct = FikaInducedDipolesDriver().compute(system, None,
                                                    self._switch())
        multipole_drv = FikaInducedDipolesDriver()
        multipole_drv.options.summation = FikaChargeSummation.multipole
        multipole = multipole_drv.compute(system, None, self._switch())
        assert np.max(np.abs(direct.dipoles - multipole.dipoles)) < 1.0e-7

        # the same with Thole damping on top of the switch
        direct = FikaInducedDipolesDriver().compute(system,
                                                    FikaTholeDamping(),
                                                    self._switch())
        multipole = multipole_drv.compute(system, FikaTholeDamping(),
                                          self._switch())
        assert np.max(np.abs(direct.dipoles - multipole.dipoles)) < 1.0e-7

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_embedding_switching_settings(self):

        system = FikaClassicalSystem()
        system.add_force_field(FikaForceField('p', 'p', [0.0], {0: 1.0}))
        system.add_residue(
            FikaResidue('p', 0, 'p', [6], np.array([[0.0, 0.0, 0.0]]), True))
        options = {
            'settings': {
                'embedding_method': 'fika',
                'damping': None,
            },
            'inputs': {
                'objects': {
                    'classical_system': system,
                },
            },
        }
        molecule = get_acrolein()
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

        # The inner switching function defaults to 1.0/1.9 angstrom and
        # damping defaults to None.
        fika_embedding_sanity_check(options)
        embedding = FikaEmbedding(molecule, basis, options)
        assert embedding.damping is None
        assert embedding.switching.r_on == pytest.approx(
            1.0 / bohr_in_angstroms())
        assert embedding.switching.r_off == pytest.approx(
            1.9 / bohr_in_angstroms())

        # An explicit None disables the switching function.
        options['settings']['switching'] = None
        fika_embedding_sanity_check(options)
        embedding = FikaEmbedding(molecule, basis, options)
        assert embedding.switching is None

        options['settings']['switching'] = {
            'r_on': 1.0,
            'r_off': 1.9,
        }
        fika_embedding_sanity_check(options)
        embedding = FikaEmbedding(molecule, basis, options)
        assert embedding.switching.r_on == pytest.approx(
            1.0 / bohr_in_angstroms())
        assert embedding.switching.r_off == pytest.approx(
            1.9 / bohr_in_angstroms())

        # radii in bohr are used as given
        options['settings']['switching'] = {
            'r_on': 2.5,
            'r_off': 3.2,
            'unit': 'bohr',
        }
        fika_embedding_sanity_check(options)
        embedding = FikaEmbedding(molecule, basis, options)
        assert embedding.switching.r_on == pytest.approx(2.5)
        assert embedding.switching.r_off == pytest.approx(3.2)

        for switching in [
                {'r_on': 1.8, 'r_off': 1.0},
                {'r_on': 1.0},
                {'r_off': 1.8},
                {'r_on': 1.0, 'r_off': 1.8, 'unit': 'nm'},
                {'r_on': 1.0, 'r_off': 1.8, 'unknown': 1},
        ]:
            options['settings']['switching'] = switching
            with pytest.raises(Exception):
                fika_embedding_sanity_check(options)

        # Thole damping and switching are mutually exclusive.
        options['settings']['damping'] = 'thole'
        options['settings']['switching'] = {
            'r_on': 1.0,
            'r_off': 1.9,
        }
        with pytest.raises(Exception, match='mutually exclusive'):
            fika_embedding_sanity_check(options)

        # Thole damping suppresses the default switching function.
        options['settings'].pop('switching')
        fika_embedding_sanity_check(options)
        embedding = FikaEmbedding(molecule, basis, options)
        assert embedding.damping is not None
        assert embedding.switching is None


class TestFikaEmbedding:

    @staticmethod
    def get_options(sources=FikaQmmmSources.all):

        options = FikaQmmmEmbeddingOptions()
        options.sources = sources
        options.induced.tolerance = 1.0e-10
        options.induced.field_accuracy = 1.0e-12
        options.fock_threshold = 1.0e-14
        return options

    @pytest.mark.solvers
    def test_ground_state(self):

        molecule = get_acrolein()
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        density = get_hf_density(molecule, basis)
        waters = get_water_shell(molecule, 40)

        ref_file = Path(__file__).parent / 'data' / 'fika_pe_reference.h5'
        with h5py.File(ref_file, 'r') as hf:
            ref_mu = np.array(hf.get('acrolein_mu'))
            ref_f_elec = np.array(hf.get('acrolein_f_elec'))
            ref_v_emb = np.array(hf.get('acrolein_v_emb'))

        result = FikaQmmmEmbeddingDriver(molecule, basis,
                                         get_classical_system(waters), None,
                                         self.get_options()).compute(density)

        assert result.induced.converged
        assert result.electron_permanent_energy == pytest.approx(
            -0.03903501058147324, abs=1.0e-10)
        assert result.nuclear_permanent_energy == pytest.approx(
            0.03590023271653367, abs=1.0e-10)
        assert result.polarization_energy == pytest.approx(
            -0.006375170012274301, abs=1.0e-9)
        assert result.energy() == pytest.approx(-0.009509947877213869,
                                                abs=1.0e-9)

        # The reference dipoles use the same sign convention as fika.
        assert np.max(np.abs(result.induced.dipoles - ref_mu)) < 1.0e-8
        assert np.max(np.abs(result.permanent_fock - ref_f_elec)) < 1.0e-12
        assert np.max(np.abs(result.fock() - ref_v_emb)) < 1.0e-9

    @pytest.mark.solvers
    def test_response(self):

        molecule = get_acrolein()
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        waters = get_water_shell(molecule, 40)

        ref_file = Path(__file__).parent / 'data' / 'fika_pe_reference.h5'
        with h5py.File(ref_file, 'r') as hf:
            ref = np.array(hf.get('acrolein_response'))

        # A non-symmetric perturbed density: only its symmetric part acts.
        rng = np.random.default_rng(4)
        n = basis.get_dimensions_of_basis()
        perturbed = 1.0e-2 * rng.normal(size=(n, n))
        perturbed += 0.3 * perturbed.T

        result = FikaQmmmEmbeddingDriver(
            molecule, basis, get_classical_system(waters), None,
            self.get_options()).compute(perturbed,
                                        FikaQmmmSources.electrons_only)

        assert result.energy() == 0.0
        assert result.permanent_fock.shape == (0, 0)
        assert np.max(np.abs(result.induced_fock - ref)) < 1.0e-9
        assert np.max(np.abs(result.fock() - ref)) < 1.0e-9

    def test_warm_start_and_fock_driver(self):

        molecule = get_acrolein()
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        density = get_hf_density(molecule, basis)
        system = get_classical_system(get_water_shell(molecule, 40))

        # Converged tightly, so the restart needs no rescaling (s = 1).
        tight = FikaQmmmEmbeddingOptions()
        tight.induced.tolerance = 1.0e-11
        first = FikaQmmmEmbeddingDriver(molecule, basis, system,
                                        options=tight).compute(density)

        driver = FikaQmmmEmbeddingDriver(molecule, basis, system)
        second = driver.compute(density, initial_guess=first.induced.dipoles)

        assert second.induced.guess_used
        assert second.induced.guess_scale == pytest.approx(1.0, abs=1.0e-8)
        assert second.induced.iterations == 0

        positions, _, _ = fika_polarizable_sites(system)
        fock = FikaInducedDipoleFockDriver().compute(molecule, basis, positions,
                                                     second.induced.dipoles)
        assert np.max(np.abs(fock - second.induced_fock)) < 1.0e-14

        # The driver reuses its geometry-only parts: a repeated computation
        # gives the same bits, and the sites are those of the system.
        again = driver.compute(density, initial_guess=first.induced.dipoles)
        assert np.array_equal(again.induced.dipoles, second.induced.dipoles)
        assert np.array_equal(again.fock(), second.fock())
        assert np.array_equal(driver.site_positions(), positions)
