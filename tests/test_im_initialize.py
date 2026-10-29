"""
Pre-merge test phase, step 2: state and persistence integrity, as a user
meets it.

A user builds a database with IMForceFieldGenerator (set_up_the_system,
compute) and restarts it the same way. What the user gets back are files: the
HDF5 database with its construction history. The tests go through these
entry points and check those files.

Sections, one test class each:
    1  TestEntryChecks     settings that cannot work are refused by
                           set_up_the_system, before any QM
    2  TestDatapointFile   a datapoint is read back from the database as it
                           was written
    3  TestFrozenSettings  settings changed after set_up_the_system are
                           reported and not used
    4  TestRestart         a second compute on the database continues it
    5  TestLabelOrder      with 12 points (point_10 sorts before point_2) the
                           collector keeps every datapoint and radius with its
                           label; the one test that looks inside the collector

The QM checks of the collector, the radius fit and the reference file are
checked in step 3, on the outputs of an active construction.

QM is PBE/def2-SVP: a ScfRestrictedDriver that counts its SCF calculations
(_CountingScfDriver); IMForceFieldGenerator builds the gradient and Hessian
drivers from it. The xTB pre-screen of the QM checks (sampling_settings) is
off, so no xTB runs and no xTB mirror of the database is written. By default
compute writes the input structure as point_1 and stops there
(desired_point_density = 1), so no dynamics runs. Every test works in its own
tmp_path, because the generator also writes into the working directory.
About 1 min with OMP_NUM_THREADS=4, most of it the 12 Hessians of
TestLabelOrder.
"""

from contextlib import contextmanager
import json
import os

import h5py
import numpy as np
import pytest

from veloxchem.errorhandler import VeloxChemError
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.scfhessiandriver import ScfHessianDriver
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator
from veloxchem.interpolationdatapoint import InterpolationDatapoint
from veloxchem.imforcefieldgenerator import IMForceFieldGenerator
from veloxchem.imdatabasepointcollecter import IMDatabasePointCollecter

H2O2_XYZ = """4

O    0.000000   0.000000   0.000000
O    0.000000   0.000000   1.455237
H    0.952886   0.000000  -0.172612
H   -0.402707   0.863608   1.627850
"""

# the z-matrix IMForceFieldGenerator builds for H2O2
H2O2_Z_MATRIX = {
    'bonds': [(0, 1), (0, 2), (1, 3)],
    'angles': [(1, 0, 2), (0, 1, 3)],
    'dihedrals': [(2, 0, 1, 3)],
    'impropers': [],
}

DATABASE = 'im_database_0.h5'

# ---------------------------------------------------------------------------
# QM: PBE/def2-SVP
# ---------------------------------------------------------------------------

BASIS = 'def2-svp'


class _CountingScfDriver(ScfRestrictedDriver):
    """ScfRestrictedDriver for PBE that counts its SCF calculations."""

    def __init__(self):
        super().__init__(ostream=silent())
        self.xcfun = 'pbe'
        self.calls = 0

    def compute(self, molecule, basis, *args, **kwargs):
        self.calls += 1
        return super().compute(molecule, basis, *args, **kwargs)


def pbe_data(molecule):
    """E (Hartree), G ((natoms, 3), Hartree/bohr) and H ((3N, 3N)) of molecule, PBE/def2-SVP."""

    basis = MolecularBasis.read(molecule, BASIS, ostream=None)
    scf_driver = _CountingScfDriver()
    scf_results = scf_driver.compute(molecule, basis)
    gradient_driver = ScfGradientDriver(scf_driver)
    gradient_driver.compute(molecule, basis, scf_results)
    hessian_driver = ScfHessianDriver(scf_driver)
    hessian_driver.compute(molecule, basis)
    return scf_driver.get_scf_energy(), gradient_driver.get_gradient(), hessian_driver.hessian


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def silent():
    """An output stream without output (it never closes pytest's stdout)."""

    return OutputStream(None)


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

    return Molecule.read_xyz_string(H2O2_XYZ)


def pbe_generator(ostream=None, **settings):
    """
    An IMForceFieldGenerator for H2O2 with PBE/def2-SVP (_CountingScfDriver)
    and without the xTB pre-screen. The defaults make compute write the input
    structure as the only seed point (no conformers, no optimisation) and stop
    there. Settings are assigned as attributes.
    """

    generator = IMForceFieldGenerator(_CountingScfDriver(), ostream=ostream or silent())
    settings = {'gs_basis_set_label': BASIS, 'open_mm_platform': 'Reference',
                'sampling_settings': {'enabled': False, 'e_thrsh_kcal_per_atom': 0.1},
                'add_conformal_structures': False, 'use_minimized_structures': [False, [], []],
                'desired_point_density': 1, **settings}
    for name, value in settings.items():
        setattr(generator, name, value)
    return generator


def construct(**settings):
    """A construction as a user runs it: set_up_the_system, then compute."""

    generator = pbe_generator(**settings)
    generator.set_up_the_system(h2o2())
    generator.compute(h2o2())
    return generator


def file_contents(filename):

    with h5py.File(filename, 'r') as h5f:
        return {key: np.array(h5f[key]) for key in h5f.keys()}


def history(filename):
    """The construction history of a database (a list of dicts)."""

    with h5py.File(filename, 'r') as h5f:
        return json.loads(h5f.attrs.get('msi_construction_history', '[]'))


# ---------------------------------------------------------------------------
# 1. Entry checks
# ---------------------------------------------------------------------------


class TestEntryChecks:
    """
    set_up_the_system refuses settings that cannot work, before any QM runs.
    The error then names the wrong setting, instead of a failure later in a
    long construction.
    """

    def test_two_bond_modes_are_refused(self, tmp_path):

        with working_directory(tmp_path):
            generator = pbe_generator(use_eq_bond_length=True)
            with pytest.raises(VeloxChemError):
                generator.set_up_the_system(h2o2())

        assert generator.drivers['gs'][0].calls == 0
        assert not (tmp_path / DATABASE).exists()

    def test_restart_in_another_bond_mode_is_refused(self, tmp_path):
        """
        The database holds points of the inverse bond mode (the default). A
        restart in the eq mode finds none of them; it must stop instead of
        continuing with an empty database.
        """

        with working_directory(tmp_path):
            construct()
            restart = pbe_generator(use_inverse_bond_length=False, use_eq_bond_length=True)
            with pytest.raises(VeloxChemError):
                restart.set_up_the_system(h2o2())

        assert restart.drivers['gs'][0].calls == 0
        assert [entry['event'] for entry in history(tmp_path / DATABASE)] == ['seed stage finished']


# ---------------------------------------------------------------------------
# 2. The datapoint file
# ---------------------------------------------------------------------------


class TestDatapointFile:
    """
    InterpolationDatapoint.write_hdf5 and read_hdf5: the database format, with
    the PBE/def2-SVP data of the input structure.
    """

    @classmethod
    def setup_class(cls):

        # one calculation for all cases
        cls.qm_data = pbe_data(h2o2())

    def datapoint(self, molecule, bond_mode, mass_weighted):

        settings = {'use_inverse_bond_length': bond_mode == 'inverse',
                    'use_eq_bond_length': bond_mode == 'eq'}
        energy, gradient, hessian = self.qm_data[0], self.qm_data[1].copy(), self.qm_data[2].copy()

        datapoint = InterpolationDatapoint(H2O2_Z_MATRIX, ostream=silent())
        datapoint.update_settings(settings)
        if bond_mode == 'eq':
            datapoint.eq_bond_lengths = np.array(
                [molecule.get_distance([i + 1, j + 1], 'bohr') for i, j in H2O2_Z_MATRIX['bonds']])
        if mass_weighted:
            inv_sqrt_masses = 1.0 / np.sqrt(np.repeat(molecule.get_masses(), 3))
            datapoint.inv_sqrt_masses = inv_sqrt_masses
            gradient = (inv_sqrt_masses * gradient.reshape(-1)).reshape(-1, 3)
            hessian = inv_sqrt_masses[:, None] * hessian * inv_sqrt_masses[None, :]
        datapoint.cartesian_coordinates = molecule.get_coordinates_in_bohr()
        datapoint.energy, datapoint.gradient, datapoint.hessian = energy, gradient, hessian
        datapoint.confidence_radius = 0.7
        datapoint.imp_int_coordinates = {'bonds': [(0, 1)], 'angles': [],
                                         'dihedrals': [(2, 0, 1, 3)], 'impropers': []}
        datapoint.transform_gradient_and_hessian()
        return datapoint, settings

    @pytest.mark.parametrize('bond_mode', ['inverse', 'plain', 'eq'])
    @pytest.mark.parametrize('mass_weighted', [False, True])
    def test_write_and_read_back(self, tmp_path, bond_mode, mass_weighted):

        written, settings = self.datapoint(h2o2(), bond_mode, mass_weighted)
        written.write_hdf5(tmp_path / DATABASE, 'point_1')

        read = InterpolationDatapoint(H2O2_Z_MATRIX, ostream=silent())
        read.update_settings(settings)
        read.read_hdf5(tmp_path / DATABASE, 'point_1')

        for name in ('energy', 'internal_gradient', 'internal_hessian', 'gradient', 'hessian',
                     'internal_coordinates_values', 'cartesian_coordinates', 'confidence_radius',
                     'eq_bond_lengths'):
            assert np.array_equal(np.asarray(getattr(read, name)), np.asarray(getattr(written, name))), name
        assert read.z_matrix_dict == written.z_matrix_dict
        assert read.imp_int_coordinates == written.imp_int_coordinates

        # the metric and the bond mode of the stored Cartesian data
        suffix = {'inverse': 'rinv', 'plain': 'r', 'eq': 'eq'}[bond_mode]
        with h5py.File(tmp_path / DATABASE, 'r') as h5f:
            attrs = h5f[f'point_1_{suffix}_dihedral_cartesian_hessian'].attrs
            assert bool(attrs['mass_weighted']) == mass_weighted
            assert attrs['bond_mode'] == suffix


# ---------------------------------------------------------------------------
# 3. Frozen settings
# ---------------------------------------------------------------------------


class TestFrozenSettings:
    """
    set_up_the_system freezes the interpolation settings. Attributes changed
    afterwards are not used (known issue K24, fixed), and compute warns about
    each of them.
    """

    changes = {'use_mass_weight': False, 'exponent_p': 4, 'use_tc_weights': False,
               'use_inverse_bond_length': False}

    def test_changes_after_set_up_are_reported_and_not_used(self, tmp_path):

        with working_directory(tmp_path):
            ostream = OutputStream(tmp_path / 'construction.out')
            generator = pbe_generator(ostream=ostream)
            generator.set_up_the_system(h2o2())
            frozen = dict(generator.states_interpolation_settings[0])

            for name, value in self.changes.items():
                setattr(generator, name, value)
            generator.compute(h2o2())
            ostream.close()

        output = (tmp_path / 'construction.out').read_text()
        for name in self.changes:
            assert f'"{name}" was changed after set_up_the_system' in output

        # the seed point has the bond mode and the metric of the set-up
        with h5py.File(tmp_path / DATABASE, 'r') as h5f:
            attrs = h5f['point_1_rinv_dihedral_cartesian_hessian'].attrs
            assert attrs['bond_mode'] == 'rinv'
            assert bool(attrs['mass_weighted'])
        frozen.pop('imforcefield_file')
        assert history(tmp_path / DATABASE)[-1]['interpolation_settings'] == frozen


# ---------------------------------------------------------------------------
# 4. Restart
# ---------------------------------------------------------------------------


class TestRestart:
    """A second construction on the same database continues it."""

    def test_restart_continues_the_database(self, tmp_path):
        """
        With desired_point_density reached, the restart records itself in the
        history and stops: no QM, and the database keeps its points.
        """

        with working_directory(tmp_path):
            construct()
            before = file_contents(DATABASE)
            restart = construct()

        assert restart.drivers['gs'][0].calls == 0
        assert restart.im_results['n_datapoints'] == {0: 1}
        assert ([entry['event'] for entry in history(tmp_path / DATABASE)] ==
                ['seed stage finished', 'construction restarted'])
        after = file_contents(tmp_path / DATABASE)
        assert sorted(after) == sorted(before)
        for key, value in before.items():
            assert np.array_equal(after[key], value), key


# ---------------------------------------------------------------------------
# 5. Label order with 10 or more points
# ---------------------------------------------------------------------------


def loaded_collector(generator, molecule):
    """
    The collector of generator.compute after it has loaded the database and
    the mirror: set up as compute sets it up, with the dynamics settings of
    the last compute, then one MD step of run_qmmm with the QM checks held
    off (start_collect).
    """

    ff_gen = MMForceFieldGenerator(ostream=silent())
    ff_gen.partial_charges = molecule.get_partial_charges(molecule.get_charge())
    ff_gen.create_topology(molecule)

    collector = IMDatabasePointCollecter(ostream=silent())
    collector.non_core_symmetry_groups = generator.symmetry_information
    collector.platform = generator.open_mm_platform
    collector.consider_locality = generator.consider_locality
    collector.identify_relevant_int_coordinates = (generator.identify_relevant_int_coordinates,
                                                   generator.use_minimized_structures[1])
    collector.use_opt_confidence_radius = generator.use_opt_confidence_radius
    collector.system_from_molecule(molecule, generator.roots_z_matrix, ff_gen,
                                   solvent=generator.solvent, qm_atoms='all')
    collector.density_around_data_point = generator.determine_datapoint_density(
        generator.states_interpolation_settings)

    settings = dict(generator.dynamics_settings, nsteps=1, snapshots=1, start_collect=10**9)
    collector.update_settings(settings, generator.states_interpolation_settings,
                              generator.sampling_states_interpolation_settings)
    collector.run_qmmm()
    return collector


class TestLabelOrder:
    """
    HDF5 lists the keys in byte order, so point_10 comes before point_2. The
    radius fit gives radius k to datapoint k, so the collector must load the
    database in number order (known issue K11, fixed). In the files this
    shows only after a fit, deep in a long construction; this test therefore
    looks at the collector's datapoints directly.
    """

    def test_collector_keeps_data_and_radius_with_the_label(self, tmp_path):

        labels = [f'point_{k}' for k in range(1, 13)]
        # distinct radii and geometries, so that data under another label is detected
        radii = {label: 0.30 + 0.03 * k for k, label in enumerate(labels, start=1)}
        geometries = {}

        with working_directory(tmp_path):
            molecule = h2o2()
            generator = pbe_generator(desired_point_density=12)
            generator.set_up_the_system(molecule)

            # a database of 12 points at displaced geometries
            basis = MolecularBasis.read(molecule, generator.gs_basis_set_label, ostream=None)
            rng = np.random.default_rng(1)
            for label in labels:
                geometries[label] = molecule.get_coordinates_in_bohr() + rng.normal(scale=0.05, size=(4, 3))
                point = Molecule(molecule.get_labels(), geometries[label], 'bohr')
                generator.add_point([(point, basis, [0], [])], generator.states_interpolation_settings)
                InterpolationDatapoint(H2O2_Z_MATRIX, ostream=silent()).update_confidence_radius(
                    DATABASE, label, radii[label])

            # the restart stops (12 points reached) and leaves the dynamics settings for the collector
            generator.compute(molecule)
            collector = loaded_collector(generator, molecule)

        assert collector.sorted_state_spec_im_labels[0] == labels
        for label, datapoint in zip(labels, collector.qm_data_point_dict[0], strict=True):
            assert float(datapoint.confidence_radius) == radii[label]
            assert np.array_equal(datapoint.cartesian_coordinates, geometries[label])
