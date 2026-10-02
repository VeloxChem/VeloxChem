import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.oneeints import compute_electric_field_integrals
from veloxchem.oneeints import compute_nuclear_potential_integrals
from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import (FikaChargeSummation, FikaClassicalSystem,
                                    FikaDipolePotentialDriver, FikaForceField,
                                    FikaInducedDipoleFockDriver,
                                    FikaInducedDipolesDriver,
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


def get_pyframe_embedding(molecule, waters, tmp_path):
    """
    PyFraME input of the same system: point charges and isotropic
    polarizabilities (zero on hydrogens, which PyFraME needs) without
    damping, each water excluded from its own sites.
    """

    import json

    fragments = []
    for w, water in enumerate(waters):
        first = 3 * w + 1
        atoms = []
        for k, element in enumerate(['O', 'H', 'H']):
            alpha = AHLSTROM_ALPHA if element == 'O' else 0.0
            atoms.append({
                'index': first + k,
                'element': element,
                'coordinate': list(water[k]),
                'multipoles': {
                    'elements': [AHLSTROM_CHARGES[k]]
                },
                'exclusions': [first, first + 1, first + 2],
                'polarizabilities': {
                    'elements':
                        [0.0, 0.0, 0.0, 0.0, alpha, 0.0, 0.0, alpha, 0.0, alpha],
                    'order': [1, 1]
                },
            })
        fragments.append({'index': w + 1, 'name': 'HOH', 'atoms': atoms})

    nuclei = [{
        'index': i + 1,
        'element': label,
        'charge': float(charge),
        'coordinate': list(xyz)
    } for i, (label, charge, xyz) in enumerate(
        zip(molecule.get_labels(), molecule.get_element_ids(),
            molecule.get_coordinates_in_bohr()))]

    json_file = tmp_path / 'fika_pe.json'
    json_file.write_text(
        json.dumps({
            'quantum_subsystems': [{
                'nuclei': nuclei
            }],
            'classical_subsystems': [{
                'classical_fragments': fragments
            }],
        }))

    return {
        'settings': {
            'embedding_method': 'PE',
            'induced_dipoles': {
                'solver': 'jidiis',
                'threshold': 1e-10,
                'max_iterations': 200,
            },
        },
        'inputs': {
            'json_file': str(json_file),
        },
    }


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
    def test_ground_state_against_pyframe(self, tmp_path):

        pytest.importorskip('pyframe')
        from veloxchem.embedding import PolarizableEmbeddingSCF

        molecule = get_acrolein()
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        density = get_hf_density(molecule, basis)
        waters = get_water_shell(molecule, 40)

        pe = PolarizableEmbeddingSCF(
            molecule, basis, get_pyframe_embedding(molecule, waters, tmp_path))
        e_emb, v_emb = pe.compute_pe_contributions(density)

        result = FikaQmmmEmbeddingDriver(molecule, basis,
                                         get_classical_system(waters), None,
                                         self.get_options()).compute(density)

        assert result.induced.converged
        assert result.electron_permanent_energy == pytest.approx(
            np.sum(pe._f_elec_es * density), abs=1.0e-10)
        assert result.nuclear_permanent_energy == pytest.approx(pe._e_nuc_es,
                                                                abs=1.0e-10)
        assert result.polarization_energy == pytest.approx(pe._e_induction,
                                                           abs=1.0e-9)
        assert result.energy() == pytest.approx(e_emb, abs=1.0e-9)

        # PyFraME stores the induced dipoles as -mu.
        mu = -np.asarray(
            pe.classical_subsystem.induced_dipoles.induced_dipoles).reshape(
                -1, 3)[0::3]
        assert np.max(np.abs(result.induced.dipoles - mu)) < 1.0e-8
        assert np.max(np.abs(result.permanent_fock - pe._f_elec_es)) < 1.0e-12
        assert np.max(np.abs(result.fock() - v_emb)) < 1.0e-9

    @pytest.mark.solvers
    def test_response_against_pyframe(self, tmp_path):

        pytest.importorskip('pyframe')
        from veloxchem.embedding import PolarizableEmbeddingLRS

        molecule = get_acrolein()
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        waters = get_water_shell(molecule, 40)

        # A non-symmetric perturbed density: only its symmetric part acts.
        rng = np.random.default_rng(4)
        n = basis.get_dimensions_of_basis()
        perturbed = 1.0e-2 * rng.normal(size=(n, n))
        perturbed += 0.3 * perturbed.T

        lrs = PolarizableEmbeddingLRS(
            molecule, basis, get_pyframe_embedding(molecule, waters, tmp_path))
        ref = lrs.compute_pe_contributions(0.5 * (perturbed + perturbed.T))

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
