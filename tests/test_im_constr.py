"""
Pre-merge test phase, step 3: the construction workflow with active QM checks.

One construction as a user runs it: IMForceFieldGenerator.set_up_the_system
and compute, for H2O2 with PBE/def2-SVP. The seed stage optimises the input
structure. The collector then runs NVE dynamics on the interpolated surface,
checks it with QM along the trajectory, adds datapoints, keeps reference
structures and fits the confidence radii. The tests look at what the user
gets: the database, the reference file of the QM checks, the output, the
dynamics record of the generator, and MSI dynamics on the result.

Sections, one test class each:
    1  TestConstruction  one construction (in setup_class, about 90 s with
                         OMP_NUM_THREADS=4): known answers (the seed is a
                         minimum of the PBE surface, every datapoint and
                         reference holds the PBE data of its structure),
                         every QM check of the output is recorded, the
                         radius fit ran, OpenMMIMDynamics on the result
                         gives the constructed surface and conserves the
                         energy, and the run equals the reference run within
                         the run-to-run spread of the DFT
    2  TestRestart       a restart on the reference run's files continues
                         the construction
    3  the opt-in speed-ups of the reference growth (default off):
       TestReferenceGrowthCaches           on the data of the reference run:
                         the gate pre-filter gives the decisions of the full
                         scan; the cached radius fit gives the default fit to
                         the tolerance of the optimiser, and falls back where
                         it does not apply
       TestConstructionWithGatePrefilter   the tests of TestConstruction on a
                         construction with gate_distance_prefilter; it is the
                         reference run
       TestConstructionWithFitCache        the tests of TestConstruction on a
                         construction with trust_radius_fit_cache, except the
                         reference comparison: the run may amplify the
                         round-off of the fit beyond the reference tolerances
"""

from contextlib import contextmanager
from pathlib import Path
import json
import os
import random
import re
import shutil
import sys
import tempfile
import time

import h5py
import numpy as np
import openmm as mm
import openmm.unit as unit

import veloxchem
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.scfhessiandriver import ScfHessianDriver
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator
from veloxchem.imforcefieldgenerator import IMForceFieldGenerator
from veloxchem.imdatabasepointcollecter import IMDatabasePointCollecter
from veloxchem.interpolationdatapoint import InterpolationDatapoint
from veloxchem.interpolationdriver import InterpolationDriver
from veloxchem.openmmimdynamics import OpenMMIMDynamics
from veloxchem.veloxchemlib import bohr_in_angstrom, hartree_in_kcalpermol, hartree_in_kjpermol

DATA = Path(__file__).parent / 'data'
DATABASE = 'im_database_0.h5'
REFERENCES = 'QM_ref_along_traj.h5'
QUALITY_CHECK = 'random_xyz_structures.xyz'

H2O2_XYZ = """4

O   0.000   0.700   0.000
O   0.000  -0.700   0.000
H   0.900   0.950   0.300
H  -0.600  -0.950   0.800
"""

BASIS = 'def2-svp'

# the construction as a user sets it up; the input structure is the only seed
SETTINGS = {
    'gs_basis_set_label': BASIS,
    'add_conformal_structures': False,
    'open_mm_platform': 'Reference',
    'ensemble': 'NVE',
    'timestep': 0.5,
    'nsteps': 3000,
    'snapshots': 30,
    'duration': 500,
    'converged_cycle': 3,
    'desired_point_density': 15,
    'energy_threshold': 0.05,
    'gradient_rmsd_thrsh': 1.5,
    'use_opt_confidence_radius': [True, 'multi_grad', 0.5, 0.3],
    'load_system': 'start_state.xml',
    'sampling_settings': {'enabled': False, 'e_thrsh_kcal_per_atom': 0.1},
}
START_TEMPERATURE = 600.0

# geomeTRIC's default convergence criterion of the largest gradient component, Eh/bohr
GEOMETRIC_GMAX = 4.5e-4

# ---------------------------------------------------------------------------
# QM: PBE/def2-SVP
# ---------------------------------------------------------------------------


def pbe_driver():
    """The QM of the construction: a ScfRestrictedDriver for PBE."""

    scf_driver = ScfRestrictedDriver(ostream=OutputStream(None))
    scf_driver.xcfun = 'pbe'
    return scf_driver


def pbe_data(molecule, hessian=True):
    """
    E (Hartree), G ((natoms, 3), Hartree/bohr) and, with hessian, H ((3N,
    3N); else None) of molecule, PBE/def2-SVP with fresh drivers.
    """

    basis = MolecularBasis.read(molecule, BASIS, ostream=None)
    scf_driver = pbe_driver()
    scf_results = scf_driver.compute(molecule, basis)
    gradient_driver = ScfGradientDriver(scf_driver)
    gradient_driver.compute(molecule, basis, scf_results)
    if not hessian:
        return scf_driver.get_scf_energy(), gradient_driver.get_gradient(), None
    hessian_driver = ScfHessianDriver(scf_driver)
    hessian_driver.compute(molecule, basis)
    return scf_driver.get_scf_energy(), gradient_driver.get_gradient(), hessian_driver.hessian


# ---------------------------------------------------------------------------
# The construction
# ---------------------------------------------------------------------------


@contextmanager
def working_directory(path):
    """Runs the block in path and returns to the previous directory."""

    previous = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def h2o2():
    """The input structure: H2O2 at a dihedral of 150 degrees."""

    molecule = Molecule.read_xyz_string(H2O2_XYZ)
    molecule.set_dihedral_in_degrees([3, 1, 2, 4], 150.0, verbose=False)
    return molecule


def write_start_state(molecule, temperature, filename):
    """
    An OpenMM state of molecule with velocities drawn at temperature (seed
    1234), without centre-of-mass motion: the collector and OpenMMIMDynamics
    remove it, so a net momentum would be lost in the first step.
    """

    masses = np.array(molecule.get_masses())
    system = mm.System()
    for mass in masses:
        system.addParticle(float(mass))
    context = mm.Context(system, mm.VerletIntegrator(0.0005), mm.Platform.getPlatformByName('Reference'))
    context.setPositions(molecule.get_coordinates_in_angstrom() * 0.1)
    context.setVelocitiesToTemperature(temperature * unit.kelvin, 1234)
    velocities = context.getState(getVelocities=True).getVelocities(asNumpy=True).value_in_unit(
        unit.nanometer / unit.picosecond)
    velocities -= masses @ velocities / masses.sum()
    context.setVelocities(velocities)
    with open(filename, 'w') as fh:
        fh.write(mm.XmlSerializer.serialize(context.getState(getPositions=True, getVelocities=True)))


def construct(workdir, **settings):
    """
    Runs the construction in workdir as a user runs it, with SETTINGS and the
    given changes: set_up_the_system, then compute. The output is written to
    construction.out.

    :return:
        The generator.
    """

    with working_directory(workdir):
        np.random.seed(0)
        random.seed(0)
        molecule = h2o2()
        write_start_state(molecule, START_TEMPERATURE, 'start_state.xml')

        ostream = OutputStream('construction.out')
        generator = IMForceFieldGenerator(pbe_driver(), ostream=ostream)
        for name, value in {**SETTINGS, **settings}.items():
            setattr(generator, name, value)

        generator.set_up_the_system(molecule)
        generator.compute(molecule)
        ostream.close()

    return generator


def dynamics_record(generator):
    """
    The structure, the interpolated energy (kJ/mol), the kinetic energy
    (kJ/mol) and the temperature (K) at every MD step of the collector. The
    temperature is 2 E_kin / (dof R) with dof = 3N - 3, as the collector
    counts it (the centre-of-mass motion is removed).
    """

    molecules = generator.state_specific_molecules[0]
    kinetic = np.array(generator.kinetic_energies[0])
    dof = 3 * molecules[0].number_of_atoms() - 3
    gas_constant = unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoule_per_mole / unit.kelvin)
    return {
        'coordinates': np.array([m.get_coordinates_in_bohr() for m in molecules]),
        'energy': np.array(generator.qm_energies[0]),
        'kinetic_energy': kinetic,
        'temperature': 2.0 * kinetic / (dof * gas_constant),
    }


def write_reference():
    """Runs the construction and stores it in data/ as the reference run."""

    DATA.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as workdir:
        generator = construct(workdir)
        shutil.copy(Path(workdir, DATABASE), DATA / 'im_h2o2_pbe_database.h5')
        shutil.copy(Path(workdir, REFERENCES), DATA / 'im_h2o2_pbe_reference.h5')
        with h5py.File(DATA / 'im_h2o2_pbe_dynamics.h5', 'w') as h5f:
            for name, value in dynamics_record(generator).items():
                h5f[name] = value
            h5f.attrs['veloxchem_version'] = veloxchem.__version__
            h5f.attrs['written'] = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    print(f'reference run written to {DATA}')


# ---------------------------------------------------------------------------
# Reading the files
# ---------------------------------------------------------------------------


def labels(filename):
    """The datapoint labels of a database, in number order."""

    with h5py.File(filename, 'r') as h5f:
        found = {'_'.join(key.split('_')[:2]) for key in h5f.keys() if key.startswith('point_')}
    return sorted(found, key=lambda label: int(label.split('_')[1]))


def datapoints(filename):
    """{label: {'coordinates', 'energy', 'gradient', 'hessian', 'radius'}} of a database (inverse bond mode)."""

    points = {}
    with h5py.File(filename, 'r') as h5f:
        for label in labels(filename):
            key = f'{label}_rinv_dihedral'
            points[label] = {
                'coordinates': h5f[f'{key}_cartesian_coordinates'][()],
                'energy': float(np.ravel(h5f[f'{key}_energy'][()])[0]),
                'gradient': h5f[f'{key}_cartesian_gradient'][()],
                'hessian': h5f[f'{key}_cartesian_hessian'][()],
                'radius': float(h5f[f'{key}_confidence_radius'][()]),
            }
    return points


def references(filename):
    """The QM checks of a reference file: qm_energy, im_energy (Hartree), coordinates (bohr), qm_gradient."""

    with h5py.File(filename, 'r') as h5f:
        group = h5f['states/0']
        return {
            'qm_energy': group['qm_energy'][()],
            'im_energy': group['im_energy'][()],
            'coordinates': np.array([np.asarray(c).reshape(-1, 3) / bohr_in_angstrom()
                                     for c in group['coords_flat'][()]]),
            'qm_gradient': np.array([np.asarray(g).reshape(-1, 3) for g in group['qm_grad_flat'][()]]),
        }


def history(filename):

    with h5py.File(filename, 'r') as h5f:
        return json.loads(h5f.attrs['msi_construction_history'])


def molecule_at(coordinates):
    """An H2O2 molecule at coordinates (bohr)."""

    return Molecule(h2o2().get_labels(), coordinates, 'bohr')


def assert_close(new, reference, tol, name):

    new, reference = np.asarray(new, dtype=float), np.asarray(reference, dtype=float)
    assert new.shape == reference.shape, name
    assert np.max(np.abs(new - reference), initial=0.0) < tol, name


# The tolerances of the comparison with the reference run (atomic units; the
# energies of the dynamics record in kJ/mol, its temperature in K). Threaded
# DFT sums in a varying order, and the radius fit amplifies that round-off. In
# 8 runs with 1 to 8 threads the largest spread, after a fit that landed 4e-4
# bohr away in one radius, was 2e-5 bohr in a structure, 2e-7 Eh in an
# energy, 4e-6 Eh/bohr in a gradient, 2e-5 in a Hessian and 4e-4 bohr in a
# radius, with the same decision at every QM check; the tolerances are 25 to
# 65 times that.
REFERENCE_TOLERANCES = {
    'coordinates': 1.0e-3,
    'energy': 1.0e-5,
    'gradient': 1.0e-4,
    'hessian': 1.0e-3,
    'confidence_radius': 1.0e-2,
    'dynamics_energy': 1.0e-2,
    'temperature': 0.5,
}


def reference_tolerance(key):
    """The tolerance of a dataset of the database by its name; z-matrix and bond lengths: 1e-8."""

    for quantity in ('coordinates', 'energy', 'gradient', 'hessian', 'confidence_radius'):
        if key.endswith(quantity):
            return REFERENCE_TOLERANCES[quantity]
    return 1.0e-8


# ---------------------------------------------------------------------------
# 1. One construction
# ---------------------------------------------------------------------------


class TestConstruction:
    """
    One construction with SETTINGS: the optimised input structure is the
    seed, then 3000 NVE steps at about 700 K with QM checks.
    """

    # changes to SETTINGS (the subclasses of section 3 run the same tests on other settings)
    settings = {}

    @classmethod
    def setup_class(cls):

        cls.workdir = tempfile.TemporaryDirectory()
        cls.path = Path(cls.workdir.name)
        cls.generator = construct(cls.path, **cls.settings)
        cls.output = (cls.path / 'construction.out').read_text()

    @classmethod
    def teardown_class(cls):

        cls.workdir.cleanup()

    def test_seed_is_the_optimised_minimum(self):
        """
        The seed stage optimises the input structure (dihedral 150 degrees,
        largest gradient component 0.035 Eh/bohr) with geomeTRIC. point_1 is a
        minimum of the PBE surface: its gradient meets geomeTRIC's criterion,
        and its Hessian in the internal coordinates of the database (3N - 6 =
        6 of them) is positive definite. That these are the PBE data of the
        structure is test_datapoints_hold_the_qm_data.
        """

        seed = datapoints(self.path / DATABASE)['point_1']
        sqrt_masses = np.sqrt(np.repeat(h2o2().get_masses(), 3))
        with h5py.File(self.path / DATABASE, 'r') as h5f:
            internal_hessian = h5f['point_1_rinv_dihedral_hessian'][()]

        # the stored gradient is mass-weighted
        assert np.max(np.abs(seed['gradient'].reshape(-1) * sqrt_masses)) < GEOMETRIC_GMAX
        assert internal_hessian.shape == (6, 6)
        assert np.min(np.linalg.eigvalsh(internal_hessian)) > 0.0
        # the _org copy holds the seed only
        assert labels(self.path / 'im_database_0_org.h5') == ['point_1']
        assert [entry['event'] for entry in history(self.path / DATABASE)] == ['seed stage finished']

    def test_datapoints_hold_the_qm_data(self):
        """
        Every datapoint stores the PBE E, G and H of its geometry (G and H
        mass-weighted), recomputed here with fresh drivers. The SCF converges
        to 1e-6 (the default), so recomputed gradients and Hessians differ by
        up to about 5e-8.
        """

        stored = datapoints(self.path / DATABASE)
        announced = re.findall(r'Added point (point_\d+) to the database', self.output)

        assert len(stored) >= 5
        assert list(stored) == [f'point_{k}' for k in range(1, len(stored) + 1)]
        # the seed and every addition, each once and in this order
        assert announced == list(stored)
        assert self.generator.im_results['n_datapoints'] == {0: len(stored)}

        inv_sqrt_masses = 1.0 / np.sqrt(np.repeat(h2o2().get_masses(), 3))
        for label, point in stored.items():
            energy, gradient, hessian = pbe_data(molecule_at(point['coordinates']))
            assert abs(point['energy'] - energy) < 1.0e-9, label
            assert np.max(np.abs(point['gradient'].reshape(-1) - inv_sqrt_masses * gradient.reshape(-1))) < 1.0e-6, label
            assert np.max(np.abs(point['hessian'] - np.outer(inv_sqrt_masses, inv_sqrt_masses) * hessian)) < 1.0e-6, label

    def test_every_qm_check_is_recorded(self):
        """
        Every QM check of the output follows the thresholds and is recorded.
        A check that exceeds a threshold adds a datapoint ('yes') or, before
        20 references exist, is kept as a reference for the radius fit. A check
        within both thresholds is stored as a reference ('no'). So is the first
        check after an addition, whatever its result: it triggers the refit of
        the radii, and no datapoint is added before that.
        """

        checks = re.findall(r'Energy difference/atom:\s+(\S+) kcal/mol\s+Gradient RMSD:\s+(\S+) kcal/mol/angstrom\s+'
                            r'Add datapoint:\s+(yes|no \(threshold exceeded|no)', self.output)
        stored = references(self.path / REFERENCES)
        natoms = h2o2().number_of_atoms()

        assert len(checks) == self.output.count('Performing point correlation check')
        decisions = [decision for _, _, decision in checks]
        assert decisions.count('yes') == len(labels(self.path / DATABASE)) - 1

        not_added = []
        for k, (energy_error, gradient_error, decision) in enumerate(checks):
            exceeded = (float(energy_error) > SETTINGS['energy_threshold'] or
                        float(gradient_error) > SETTINGS['gradient_rmsd_thrsh'])
            if decision == 'no':
                assert not exceeded or decisions[k - 1] == 'yes', k
            else:
                assert exceeded, k
            if decision != 'yes':
                not_added.append(float(energy_error))

        # the references of the file, in the order of the checks
        assert len(stored['qm_energy']) == len(not_added) >= 20
        for k, coordinates in enumerate(stored['coordinates']):
            # the printed energy error is the one of the stored energies
            error = abs(stored['qm_energy'][k] - stored['im_energy'][k]) / natoms * hartree_in_kcalpermol()
            assert abs(error - not_added[k]) < 1.0e-6, k
            # the PBE data of the structure (recomputed gradients differ by up to about 1e-7)
            energy, gradient, _ = pbe_data(molecule_at(coordinates), hessian=False)
            assert abs(stored['qm_energy'][k] - energy) < 1.0e-9, k
            assert np.max(np.abs(stored['qm_gradient'][k] - gradient)) < 1.0e-6, k
            # no structure is stored twice
            for other in stored['coordinates'][:k]:
                assert np.max(np.abs(coordinates - other)) > 1.0e-6, k

    def test_radius_fit_ran(self):

        radii = [point['radius'] for point in datapoints(self.path / DATABASE).values()]

        assert self.output.count('Trust-Radius Optimization Result') >= 1
        # within the bounds of the fit, and not all at the initial radius
        assert all(1.0e-3 <= radius <= 1.0 for radius in radii)
        assert any(radius != SETTINGS['use_opt_confidence_radius'][2] for radius in radii)

    def frozen_dynamics(self, timestep, nsteps):
        """
        OpenMMIMDynamics on the constructed database, with the interpolation
        settings of its history, NVE from a start state at 300 K.
        """

        settings = dict(history(self.path / DATABASE)[-1]['interpolation_settings'],
                        imforcefield_file=DATABASE)

        with working_directory(self.path):
            molecule = h2o2()
            write_start_state(molecule, 300.0, 'frozen_start.xml')
            ff_gen = MMForceFieldGenerator(ostream=OutputStream(None))
            ff_gen.partial_charges = molecule.get_partial_charges(molecule.get_charge())
            ff_gen.create_topology(molecule)

            dynamics = OpenMMIMDynamics(ostream=OutputStream(None))
            dynamics.platform = 'Reference'
            dynamics.use_symmetry = False
            dynamics.create_system_from_molecule(molecule, {0: settings}, ff_gen, solvent='gas', qm_atoms='all')
            dynamics.run_immm(restart_file='frozen_start.xml', ensemble='NVE', timestep=timestep,
                              nsteps=nsteps, snapshots=10, traj_file='frozen.pdb',
                              state_file='frozen.xml', output_file='frozen')
        return dynamics

    def test_frozen_dynamics_gives_the_constructed_surface(self):
        """
        At the structures of the final quality check, the dynamics' energy is
        the interpolated energy the construction reported (random_xyz_structures.xyz).
        """

        dynamics = self.frozen_dynamics(0.5, 100)

        blocks = (self.path / QUALITY_CHECK).read_text().strip().split('\n4\n')
        assert len(blocks) == self.generator.nstruc_to_confirm_database_quality
        for block in blocks:
            lines = block.strip().splitlines()
            if lines[0].strip() == '4':
                lines = lines[1:]
            reported = float(re.search(r'IM: (\S+)', lines[0]).group(1)) * hartree_in_kjpermol()
            angstrom = np.array([[float(v) for v in line.split()[1:4]] for line in lines[1:5]])
            assert abs(dynamics.update_potential_energy(angstrom * 0.1) - reported) < 1.0e-6

    def test_frozen_dynamics_conserves_the_energy(self):
        """
        Velocity Verlet on a smooth surface: the fluctuation of the total
        energy scales as dt^2, so halving dt over the same 200 fs reduces it
        by 4. A force inconsistent with the energy, or a discontinuous surface,
        breaks the scaling.
        """

        coarse = self.frozen_dynamics(0.5, 400)
        fine = self.frozen_dynamics(0.25, 800)

        ratio = np.std(coarse.total_energies) / np.std(fine.total_energies)
        assert 3.0 < ratio < 5.0
        assert np.std(fine.total_energies) < 0.01 * np.mean(fine.kinetic_energies)

    def test_same_as_reference_run(self):
        """
        The same database, the same QM checks and the same dynamics (structure,
        energy and temperature at every step) as the reference run in data/,
        within the run-to-run spread of the threaded DFT (REFERENCE_TOLERANCES).
        Should this test alone fail on other hardware, with the known answers
        passing, check first whether a radius fit landed in another minimum:
        the fit can amplify round-off beyond these tolerances.
        """

        # the database: every dataset, its attributes and the recorded settings
        with h5py.File(self.path / DATABASE, 'r') as new, h5py.File(DATA / 'im_h2o2_pbe_database.h5', 'r') as ref:
            assert sorted(new.keys()) == sorted(ref.keys())
            for key in ref.keys():
                assert_close(new[key][()], ref[key][()], reference_tolerance(key), key)
                assert dict(new[key].attrs) == dict(ref[key].attrs), key
        for new_entry, ref_entry in zip(history(self.path / DATABASE), history(DATA / 'im_h2o2_pbe_database.h5'), strict=True):
            for name in ('event', 'bond_mode', 'interpolation_settings'):
                assert new_entry[name] == ref_entry[name], name

        # the QM checks
        new_checks, ref_checks = references(self.path / REFERENCES), references(DATA / 'im_h2o2_pbe_reference.h5')
        tol = REFERENCE_TOLERANCES
        assert_close(new_checks['qm_energy'], ref_checks['qm_energy'], tol['energy'], 'qm_energy')
        assert_close(new_checks['im_energy'], ref_checks['im_energy'], tol['energy'], 'im_energy')
        assert_close(new_checks['coordinates'], ref_checks['coordinates'], tol['coordinates'], 'coordinates')
        assert_close(new_checks['qm_gradient'], ref_checks['qm_gradient'], tol['gradient'], 'qm_gradient')

        # the dynamics, step by step
        record = dynamics_record(self.generator)
        with h5py.File(DATA / 'im_h2o2_pbe_dynamics.h5', 'r') as h5f:
            reference = {name: h5f[name][()] for name in record}
        assert_close(record['coordinates'], reference['coordinates'], tol['coordinates'], 'coordinates')
        assert_close(record['energy'], reference['energy'], tol['dynamics_energy'], 'energy')
        assert_close(record['kinetic_energy'], reference['kinetic_energy'], tol['dynamics_energy'], 'kinetic energy')
        assert_close(record['temperature'], reference['temperature'], tol['temperature'], 'temperature')

        # the temperatures the output prints
        printed = re.findall(r'Step: (\d+) / \d+.*?Temperature: (\S+) K', self.output, flags=re.S)
        assert len(printed) == SETTINGS['snapshots']
        for step, temperature in printed:
            assert abs(float(temperature) - reference['temperature'][int(step)]) < tol['temperature']


# ---------------------------------------------------------------------------
# 2. Restart
# ---------------------------------------------------------------------------


class TestRestart:
    """A second compute on the reference run's files continues the construction."""

    def test_restart_continues_the_construction(self, tmp_path):

        shutil.copy(DATA / 'im_h2o2_pbe_database.h5', tmp_path / DATABASE)
        shutil.copy(DATA / 'im_h2o2_pbe_reference.h5', tmp_path / REFERENCES)
        before = datapoints(tmp_path / DATABASE)
        checks_before = references(tmp_path / REFERENCES)

        generator = construct(tmp_path, desired_point_density=len(before) + 1)
        output = (tmp_path / 'construction.out').read_text()

        # new points continue the numbering; the old ones keep their data (the fit may change radii)
        after = datapoints(tmp_path / DATABASE)
        assert list(after) == [f'point_{k}' for k in range(1, len(before) + 2)]
        assert generator.im_results['n_datapoints'] == {0: len(before) + 1}
        for label, point in before.items():
            for name in ('coordinates', 'energy', 'gradient', 'hessian'):
                assert np.array_equal(after[label][name], point[name]), (label, name)

        # the reference structures are preloaded and new ones appended
        n_before = len(checks_before['qm_energy'])
        assert f'Preloaded {n_before} QM reference structures' in output
        checks_after = references(tmp_path / REFERENCES)
        assert len(checks_after['qm_energy']) > n_before
        for name, values in checks_before.items():
            assert np.array_equal(checks_after[name][:n_before], values), name

        assert ([entry['event'] for entry in history(tmp_path / DATABASE)] ==
                ['seed stage finished', 'construction restarted'])


# ---------------------------------------------------------------------------
# 3. The opt-in speed-ups of the reference growth
# ---------------------------------------------------------------------------


def gate_scan_sequence(collector, molecule_list, queries, threshold, bounded, K=5):
    """
    The proximity gate of IMDatabasePointCollecter.update_forces at a sequence
    of structures: first the K cached candidates of the previous structure,
    then, without a match, the full scan in list order (the oracle) or the
    bounded scan of gate_distance_prefilter. Returns per structure whether a
    reference lies within threshold and the K closest (index, RMSD) entries
    after the stable sort of update_forces.
    """

    n_atoms = queries[0].shape[0]

    def rmsd(idx, coords):
        distance = collector.cartesian_just_distance(molecule_list[idx].get_coordinates_in_bohr(), coords, [])
        return np.linalg.norm(distance) / np.sqrt(n_atoms) * bohr_in_angstrom()

    previous, results = [], []
    for coords in queries:
        closest, added, found = [], [], False
        for idx in previous:
            d = rmsd(idx, coords)
            if idx not in added:
                closest.append((idx, d))
                added.append(idx)
            if d <= threshold:
                found = True
                break
        if not found and bounded:
            found = collector._bounded_reference_scan(molecule_list, coords, n_atoms, threshold, K, closest)
        elif not found:
            for idx in range(len(molecule_list)):
                d = rmsd(idx, coords)
                closest.append((idx, d))
                if d <= threshold:
                    found = True
                    break
        closest.sort(key=lambda entry: entry[1])
        previous = [idx for idx, _ in closest[:K]]
        results.append((found, closest[:K]))
    return results


class _CountingCollector(IMDatabasePointCollecter):
    """IMDatabasePointCollecter that counts the alignments of the proximity gate."""

    def __init__(self):
        super().__init__(ostream=OutputStream(None))
        self.alignments = 0

    def cartesian_just_distance(self, coordinate_1, coordinate_2, non_core_atoms):
        self.alignments += 1
        return super().cartesian_just_distance(coordinate_1, coordinate_2, non_core_atoms)


def symmetry_information():
    """The symmetry information IMForceFieldGenerator builds for H2O2 (the 'gs' entry)."""

    with tempfile.TemporaryDirectory() as workdir, working_directory(workdir):
        generator = IMForceFieldGenerator(pbe_driver(), ostream=OutputStream(None))
        generator.add_conformal_structures = False
        generator.set_up_the_system(h2o2())
        return generator.symmetry_information['gs']


def replay_fit(output_file, cached, shift=0.0, **setting_changes):
    """
    IMDatabasePointCollecter.determine_trust_radius_gradient on the reference
    run in data/ as the collector calls it after its last addition: the
    datapoints (fresh objects from im_h2o2_pbe_database.h5), all references but the
    newest (im_h2o2_pbe_reference.h5), the interpolation settings of the history,
    the numpy seed 0. Returns the fitted radii and the output.

    :param shift:
        Moves the first atom of the first reference by shift bohr in x (a
        round-off control of the fit).
    """

    database = str(DATA / 'im_h2o2_pbe_database.h5')
    settings = dict(history(database)[-1]['interpolation_settings'], imforcefield_file=database, **setting_changes)
    driver = InterpolationDriver(ostream=OutputStream(None))
    driver.update_settings(settings)
    labels_of_db, z_matrix = driver.read_labels()
    inv_sqrt_masses = 1.0 / np.sqrt(np.repeat(h2o2().get_masses(), 3))
    points = []
    for label in sorted(labels_of_db, key=lambda label: int(label.split('_')[1])):
        point = InterpolationDatapoint(z_matrix, ostream=OutputStream(None))
        point.update_settings(settings)
        point.read_hdf5(database, label)
        point.inv_sqrt_masses = inv_sqrt_masses
        points.append(point)

    stored = references(DATA / 'im_h2o2_pbe_reference.h5')
    n = len(stored['qm_energy']) - 1
    stored['coordinates'][0][0, 0] += shift
    molecules = [molecule_at(stored['coordinates'][k]) for k in range(n)]

    ostream = OutputStream(output_file)
    collector = IMDatabasePointCollecter(ostream=ostream)
    collector.use_opt_confidence_radius = list(SETTINGS['use_opt_confidence_radius'])
    collector.trust_radius_fit_cache = cached
    np.random.seed(0)
    radii = collector.determine_trust_radius_gradient(
        molecules, list(stored['qm_energy'][:n]), list(stored['qm_gradient'][:n]), list(stored['im_energy'][:n]),
        points, settings, symmetry_information(), z_matrix,
        exponent_p_q=(int(settings['exponent_p']), int(settings['exponent_q'])))
    ostream.close()
    return np.asarray(radii, dtype=float), Path(output_file).read_text()


def printed_block(output, start, stops):
    """The stripped, non-empty lines of output from the line containing start up to a line containing one of stops."""

    block, inside = [], False
    for line in output.splitlines():
        if start in line:
            inside = True
        elif inside and any(stop in line for stop in stops):
            break
        if inside and line.strip() and not line.strip().startswith('* Info *'):
            block.append(line.strip())
    return block


class TestReferenceGrowthCaches:
    """
    The two opt-in speed-ups of the reference growth (default off), on the
    data of the reference run: gate_distance_prefilter must give the decisions
    and candidates of the full scan; trust_radius_fit_cache the fit of the
    default path, and the default path itself where it does not apply.
    """

    def test_gate_prefilter_matches_the_full_scan(self):
        """
        The bounded scan against the full scan, with the cached candidates of
        update_forces, at the trajectory of the reference run: (1) the
        reference structures against every 20th MD step, (2) every 10th MD step
        as references against the steps between them, at the gate's 0.05
        Angstrom and at 0.005 Angstrom (fewer matches, so the K closest are
        searched). The bound must be a lower bound of every RMSD.
        """

        with h5py.File(DATA / 'im_h2o2_pbe_dynamics.h5', 'r') as h5f:
            frames = h5f['coordinates'][()]
        reference_sets = [
            ([molecule_at(c) for c in references(DATA / 'im_h2o2_pbe_reference.h5')['coordinates']], frames[::20]),
            ([molecule_at(c) for c in frames[0::10]], frames[5::10]),
        ]
        for molecule_list, queries in reference_sets:
            for threshold in (0.05, 0.005):
                full, bounded = _CountingCollector(), _CountingCollector()
                expected = gate_scan_sequence(full, molecule_list, queries, threshold, bounded=False)
                result = gate_scan_sequence(bounded, molecule_list, queries, threshold, bounded=True)
                assert result == expected
                assert bounded.alignments <= full.alignments

            collector = _CountingCollector()
            for coords in queries[::10]:
                lower = collector._reference_rmsd_lower_bounds(molecule_list, coords)
                rmsd = np.array([np.linalg.norm(collector.cartesian_just_distance(m.get_coordinates_in_bohr(), coords, []))
                                 for m in molecule_list]) / np.sqrt(len(coords)) * bohr_in_angstrom()
                # round-off only where both are zero (identical structures)
                assert np.all(lower <= rmsd + 1.0e-14)

    def test_cached_fit_matches_the_default_fit(self, tmp_path):
        """
        The same fit both ways: the same groups, the same error tables before
        the fit, the same objective to the printed digits, and the same radii
        to the tolerance of the optimiser. On the reference run the two paths
        differ by 5e-9 bohr, and the default path itself moves its radii by
        3e-9 bohr when one reference is moved by 1e-12 bohr: round-off of the
        fit, not an error. The tolerance, 1e-4 bohr, allows for fits that
        amplify round-off more (on the earlier model data: 3e-6 bohr).
        """

        default, default_out = replay_fit(tmp_path / 'default.out', cached=False)
        cached, cached_out = replay_fit(tmp_path / 'cached.out', cached=True)

        assert 'Cached radius fit:' in cached_out
        assert 'trust_radius_fit_cache is not used' not in cached_out
        stops = ('Trust-Radius Optimization Progress', 'Cached radius fit')
        assert (printed_block(cached_out, 'Trust-Radius Optimization Groups', stops) ==
                printed_block(default_out, 'Trust-Radius Optimization Groups', stops))
        for scenario in ('without_last_datapoint', 'with_new_datapoint_pre_opt'):
            assert ([line for line in cached_out.splitlines() if line.strip().startswith(scenario)] ==
                    [line for line in default_out.splitlines() if line.strip().startswith(scenario)])
        assert (re.findall(r'Final objective\s+(\S+)', cached_out) ==
                re.findall(r'Final objective\s+(\S+)', default_out) != [])
        assert np.max(np.abs(cached - default)) < 1.0e-4
        assert np.any(default != np.array([p['radius'] for p in datapoints(DATA / 'im_h2o2_pbe_database.h5').values()]))

    def test_cached_fit_falls_back_on_unsupported_settings(self, tmp_path):
        """
        With weights in internal coordinates the closed form does not apply:
        the fit warns, names the reason, and runs the default path, with
        bitwise the same radii.
        """

        default, default_out = replay_fit(tmp_path / 'default.out', cached=False, weightfunction_type='internal')
        cached, cached_out = replay_fit(tmp_path / 'cached.out', cached=True, weightfunction_type='internal')

        assert ('trust_radius_fit_cache is not used: it needs Shepard interpolation with Cartesian weights'
                in cached_out)
        assert 'Cached radius fit:' not in cached_out
        assert 'Trust-Radius Optimization Result' in default_out
        assert np.array_equal(cached, default)


class TestConstructionWithGatePrefilter(TestConstruction):
    """
    The construction of TestConstruction with gate_distance_prefilter. Its
    decisions must be those of the full scan, so every test of
    TestConstruction holds, the comparison with the reference run included.
    """

    settings = {'gate_distance_prefilter': True}

    def test_settings_are_recorded(self):

        assert history(self.path / DATABASE)[-1]['gate_distance_prefilter'] is True


class TestConstructionWithFitCache(TestConstruction):
    """
    The construction of TestConstruction with trust_radius_fit_cache. The
    cached fit equals the default fit to the tolerance of the optimiser
    (TestReferenceGrowthCaches), and the dynamics and the later fits can
    amplify that difference. With PBE the run stays within
    REFERENCE_TOLERANCES of the reference run, but a fit that lands in another
    minimum would part the two runs (as it did with a model QM), so only the
    known answers of TestConstruction are required to hold.
    """

    settings = {'trust_radius_fit_cache': True}

    # the run is not required to be the reference run (see the class docstring)
    test_same_as_reference_run = None

    def test_every_fit_used_the_cache(self):

        fits = self.output.count('Trust-Radius Optimization Groups')
        assert fits >= 1
        assert self.output.count('Cached radius fit:') == fits
        assert 'trust_radius_fit_cache is not used' not in self.output
        assert history(self.path / DATABASE)[-1]['trust_radius_fit_cache'] is True


if __name__ == '__main__':

    if sys.argv[1:] == ['--write-reference']:
        write_reference()
    else:
        print('usage: python test_step3_construction.py --write-reference')
