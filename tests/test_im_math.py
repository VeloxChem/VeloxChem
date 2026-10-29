"""
Pre-merge test phase, step 1: the mathematics of the MSI model.

Each test checks one equation of the modified Shepard interpolation
(InterpolationDatapoint, InterpolationDriver), or the force units of the
dynamics, against an oracle that is not the code under test:
    - an analytic reference potential with exact derivatives (section 1),
    - the documented formula, written again in the test. It uses the
      geometry functions of Molecule for the internal coordinates
      (get_distance, get_angle, get_dihedral) and geomeTRIC for the Kabsch
      rotation, as the code does,
    - a finite difference.

Geometries are changed with the Molecule functions set_distance, set_angle
and set_dihedral.

The equations (from the docstrings of InterpolationDatapoint and
InterpolationDriver):
    (1) internal coordinates q: bonds 1/r (inverse mode), r (plain mode) or
        r_eq - r_eq^2/r (eq mode); angles and dihedrals in radians.
    (2) B = dq/dx and B2 = d2q/dx2; with mass weighting dq/dQ, Q = sqrt(m) x.
    (3) g_q = G^- B g and H_q = G^- B (H - sum_k g_q,k B2_k) B^T G^-,
        with G = B B^T.
    (4) T_i = E_i + g_i.s + 1/2 s.H_i.s, with s = Delta q for bonds and
        angles, 2 sin(Delta/2) for dihedrals and 2 tan(Delta/2) for
        impropers; Delta is wrapped to [-pi, pi).
    (5) grad T_i = B^T [c * (g_i + H_i s)], with c = ds/dDelta: 1,
        cos(Delta/2) or 1/cos^2(Delta/2).
    (6) w_i = 1/(u^p + u^q), with u = d_i^2/R_i^2 and d_i the
        Kabsch-aligned Cartesian distance. The target-customised weights
        multiply w_i by exp(-lambda A_i). The radius fit uses dw_i/dR_i and
        d(grad w_i)/dR_i.
    (7) E = sum_i W_i T_i with W_i = w_i / sum_j w_j, and grad E is the exact
        gradient of E.
    (8) F [kJ/mol/nm] = -grad E [Eh/bohr] * 2625.4996394798 / 0.0529177210903.

Sections, one test class each:
    1  TestReferencePotential   the oracle: exact derivatives, invariance
    2  TestInternalCoordinates  equations (1) and (2)
    3  TestTransformation       equation (3)
    4  TestTaylorModel          equations (4) and (5), at large torsion
                                displacements
    5  TestShepardWeights       equation (6)
    6  TestInterpolation        equation (7): the gradient is the derivative
                                of the energy; invariance
    7  TestFidelity             the interpolant reproduces the reference
                                potential
    8  TestForceUnits           equation (8) in OpenMMIMDynamics and
                                IMDatabasePointCollecter

Tolerances (none was tuned to make a test pass):
    - 1e-10 or tighter for exact identities. The quantities are of order
      0.01 to 1 in Hartree, bohr and radian, and round-off is about 1e-15.
    - 1e-8 for comparisons with a Richardson-extrapolated central
      difference, whose error is O(h^4); observed 1e-14 to 1e-10.
    - relative tolerances for large quantities (weights, forces).
    - the Taylor order: within 0.2 of the theoretical order.

The tests use no QM: they create no SCF or xTB driver, and the code they
call (InterpolationDatapoint, InterpolationDriver, and the update_forces
methods with a fixed gradient) does not run one.
"""

import numpy as np
import pytest

from veloxchem.veloxchemlib import bohr_in_angstrom
from veloxchem.molecule import Molecule
from veloxchem.outputstream import OutputStream
from veloxchem.interpolationdatapoint import InterpolationDatapoint
from veloxchem.interpolationdriver import InterpolationDriver
from veloxchem.openmmimdynamics import OpenMMIMDynamics
from veloxchem.imdatabasepointcollecter import IMDatabasePointCollecter

import geometric

# ---------------------------------------------------------------------------
# Molecules and z-matrices
# ---------------------------------------------------------------------------

WATER_XYZ = """3

O   0.000   0.000   0.117
H   0.000   0.757  -0.467
H   0.000  -0.757  -0.467
"""

# 2 bonds + 1 angle = 3N - 6: the internal coordinates are not redundant
WATER_Z_MATRIX = {
    'bonds': [(0, 1), (0, 2)],
    'angles': [(1, 0, 2)],
    'dihedrals': [],
    'impropers': [],
}

ETHANOL_XYZ = """9

C   -0.661535066309   1.128097375722   1.298956571720
C    0.160830459103  -0.078956848612   0.840928833654
O   -0.545706317697  -1.293534712398   0.947688936392
H   -1.201152069039   1.577440238953   0.467697680596
H   -0.006154222315   1.890329398828   1.712599756591
H   -1.363731740381   0.823716932687   2.070509517949
H    1.028694858953  -0.210814780197   1.491718794709
H    0.521388572669   0.070466580703  -0.186332296889
H   -1.346953592312  -1.233309727019   0.415095469234
"""

# 8 bonds + 13 angles + 12 dihedrals = 33 > 3N - 6 = 21: redundant
ETHANOL_Z_MATRIX = {
    'bonds': [(0, 1), (0, 3), (0, 4), (0, 5), (1, 2), (1, 6), (1, 7), (2, 8)],
    'angles': [(1, 0, 3), (1, 0, 4), (1, 0, 5), (3, 0, 4), (3, 0, 5),
               (4, 0, 5), (0, 1, 2), (0, 1, 6), (0, 1, 7), (2, 1, 6),
               (2, 1, 7), (6, 1, 7), (1, 2, 8)],
    'dihedrals': [(3, 0, 1, 2), (3, 0, 1, 6), (3, 0, 1, 7), (4, 0, 1, 2),
                  (4, 0, 1, 6), (4, 0, 1, 7), (5, 0, 1, 2), (5, 0, 1, 6),
                  (5, 0, 1, 7), (0, 1, 2, 8), (6, 1, 2, 8), (7, 1, 2, 8)],
    'impropers': [],
}

# important coordinates for the target-customised weights
WATER_IMPORTANT = {'bonds': [(0, 1)], 'angles': [], 'dihedrals': [],
                   'impropers': []}
ETHANOL_IMPORTANT = {'bonds': [(1, 2)], 'angles': [],
                     'dihedrals': [(0, 1, 2, 8)], 'impropers': []}
NO_IMPORTANT = {'bonds': [], 'angles': [], 'dihedrals': [], 'impropers': []}

H2O2_XYZ = """4

O    0.000000   0.000000   0.000000
O    0.000000   0.000000   1.455237
H    0.952886   0.000000  -0.172612
H   -0.402707   0.863608   1.627850
"""

H2O2_IMPORTANT = {'bonds': [(0, 1)], 'angles': [], 'dihedrals': [(2, 0, 1, 3)],
                  'impropers': []}


def h2o2_z_matrix(torsion_section='dihedrals'):
    """
    The z-matrix of H2O2 (3 bonds, 2 angles, 1 torsion), with the torsion
    H3-O1-O2-H4 listed as a proper dihedral or as an improper.
    """

    z_matrix = {
        'bonds': [(0, 1), (0, 2), (1, 3)],
        'angles': [(1, 0, 2), (0, 1, 3)],
        'dihedrals': [],
        'impropers': [],
    }
    z_matrix[torsion_section] = [(2, 0, 1, 3)]
    return z_matrix


def h2o2_molecule(phi_degrees, r_oo=None, r_oh=None, angle_degrees=None):
    """
    H2O2 with the torsion H3-O1-O2-H4 set to phi_degrees and, if given, the
    O1-O2 and O2-H4 bond lengths (bohr) and the O2-O1-H3 angle (degrees).
    The degrees of freedom are set with the Molecule functions (1-based atom
    indices); the torsion is set last.
    """

    molecule = Molecule.read_xyz_string(H2O2_XYZ)
    if r_oo is not None:
        molecule.set_distance([1, 2], r_oo, 'bohr')
    if r_oh is not None:
        molecule.set_distance([2, 4], r_oh, 'bohr')
    if angle_degrees is not None:
        molecule.set_angle_in_degrees([2, 1, 3], angle_degrees)
    molecule.set_dihedral_in_degrees([3, 1, 2, 4], phi_degrees)
    return molecule


def moved(molecule, coords):
    """The molecule with other coordinates (bohr)."""

    return Molecule(molecule.get_labels(), coords, 'bohr')


# ---------------------------------------------------------------------------
# Geometry of the internal coordinates (the oracle for equation (1))
# ---------------------------------------------------------------------------


def one_based(indices):

    return [i + 1 for i in indices]


def coordinate_kinds(z_matrix):

    return np.array(['bond'] * len(z_matrix['bonds']) +
                    ['angle'] * len(z_matrix['angles']) +
                    ['dihedral'] * len(z_matrix['dihedrals']) +
                    ['improper'] * len(z_matrix['impropers']))


def bond_lengths(molecule, z_matrix):
    """The bond lengths of the z-matrix in bohr (Molecule.get_distance)."""

    return np.array([molecule.get_distance(one_based(b), 'bohr')
                     for b in z_matrix['bonds']])


def coordinate_values(molecule, z_matrix, bond_mode='plain',
                      eq_bond_lengths=None):
    """
    Equation (1): the internal coordinates of molecule in the flat order,
    from Molecule.get_distance, get_angle and get_dihedral (bohr, radian).
    These are VeloxChem's own geometry functions, independent of the geomeTRIC
    primitives that InterpolationDatapoint uses.
    """

    r = bond_lengths(molecule, z_matrix)
    if bond_mode == 'inverse':
        bonds = 1.0 / r
    elif bond_mode == 'eq':
        r_eq = np.asarray(eq_bond_lengths)
        bonds = r_eq - r_eq**2 / r
    else:
        bonds = r
    angles = [molecule.get_angle(one_based(a), 'radian')
              for a in z_matrix['angles']]
    torsions = [molecule.get_dihedral(one_based(d), 'radian')
                for d in z_matrix['dihedrals'] + z_matrix['impropers']]
    return np.concatenate([bonds, angles, torsions])


def wrap_angle(delta):
    """Wraps an angle difference to [-pi, pi)."""

    return (delta + np.pi) % (2.0 * np.pi) - np.pi


def kabsch_distance(x, y):
    """
    |x_c - R y_c|: x_c and y_c are x and y centred on their geometric centre,
    and R is the rotation of y_c onto x_c from geomeTRIC, called as in
    InterpolationDriver.cartesian_distance.
    """

    xc, yc = x - x.mean(axis=0), y - y.mean(axis=0)
    rotation = geometric.rotate.get_rot(yc, xc)
    return np.linalg.norm(xc - yc @ rotation.T)


def random_rotation(rng):

    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0.0:
        q[:, 0] *= -1.0
    return q


def derivative(func, coords, h=1.0e-3):
    """
    d func / d coords by Richardson-extrapolated central differences,
    (4 D(h/2) - D(h)) / 3, with error O(h^4). func maps an (natoms, 3) array
    to a scalar or an array; the result has the shape of func with one more
    axis of length 3N.
    """

    x = np.array(coords, dtype=float)

    def central(step):
        columns = []
        for k in range(x.size):
            xp, xm = x.copy().reshape(-1), x.copy().reshape(-1)
            xp[k] += step
            xm[k] -= step
            fp = np.asarray(func(xp.reshape(x.shape)))
            fm = np.asarray(func(xm.reshape(x.shape)))
            columns.append((fp - fm) / (2.0 * step))
        return np.stack(columns, axis=-1)

    return (4.0 * central(0.5 * h) - central(h)) / 3.0


# ---------------------------------------------------------------------------
# The reference potential (the 'QM' of these tests)
# ---------------------------------------------------------------------------


class ReferencePotential:
    """
    An analytic potential that depends only on bond lengths, angles and
    dihedrals. It is therefore invariant under rotations and translations,
    and its energy, gradient and Hessian are exact:

        V = sum_bonds      1/2 kb (f(r) - f0)^2,   f(r) = r or 1/r
          + sum_angles     1/2 ka (theta - theta0)^2
          + sum_dihedrals  kd (1 + cos(3 phi + delta))

    The minimum lies near the geometry passed to the constructor. Units:
    Hartree, bohr and radian.
    """

    def __init__(self, z_matrix, coords, bond_function='r', torsions=True,
                 seed=1):

        rng = np.random.default_rng(seed)
        x = np.asarray(coords, dtype=float).reshape(-1)

        self.bond_function = bond_function
        self.bonds = [geometric.internal.Distance(*b)
                      for b in z_matrix['bonds']]
        self.angles = [geometric.internal.Angle(*a)
                       for a in z_matrix['angles']]
        self.dihedrals = [geometric.internal.Dihedral(*d)
                          for d in z_matrix['dihedrals']]

        self.kb = rng.uniform(0.2, 0.5, len(self.bonds))
        self.ka = rng.uniform(0.05, 0.15, len(self.angles))
        self.kd = (rng.uniform(0.001, 0.003, len(self.dihedrals))
                   if torsions else np.zeros(len(self.dihedrals)))

        self.f0 = np.array([self.f(b.value(x)) for b in self.bonds])
        self.f0 *= rng.uniform(0.97, 1.03, len(self.bonds))
        self.theta0 = np.array([a.value(x) for a in self.angles])
        self.theta0 *= rng.uniform(0.97, 1.03, len(self.angles))
        self.delta = rng.uniform(-np.pi, np.pi, len(self.dihedrals))

    def f(self, r):

        return r if self.bond_function == 'r' else 1.0 / r

    def compute(self, coords):
        """
        Returns E (Hartree), G ((natoms, 3), Hartree/bohr) and
        H ((3N, 3N), Hartree/bohr^2).
        """

        x = np.asarray(coords, dtype=float).reshape(-1)
        n = x.size

        # every term is e(v(x)) of a primitive v:
        # (dv/dx, d2v/dx2, e, de/dv, d2e/dv2)
        terms = []

        for bond, k, f0 in zip(self.bonds, self.kb, self.f0):
            r = bond.value(x)
            dr = bond.derivative(x).reshape(-1)
            d2r = bond.second_derivative(x).reshape(n, n)
            if self.bond_function == 'r':
                f, df, d2f = r, dr, d2r
            else:
                f = 1.0 / r
                df = -dr / r**2
                d2f = -d2r / r**2 + 2.0 * np.outer(dr, dr) / r**3
            terms.append((df, d2f, 0.5 * k * (f - f0)**2, k * (f - f0), k))

        for angle, k, theta0 in zip(self.angles, self.ka, self.theta0):
            t = angle.value(x)
            dt = angle.derivative(x).reshape(-1)
            d2t = angle.second_derivative(x).reshape(n, n)
            terms.append((dt, d2t, 0.5 * k * (t - theta0)**2,
                          k * (t - theta0), k))

        for dihedral, k, delta in zip(self.dihedrals, self.kd, self.delta):
            phi = dihedral.value(x)
            dphi = dihedral.derivative(x).reshape(-1)
            d2phi = dihedral.second_derivative(x).reshape(n, n)
            terms.append((dphi, d2phi, k * (1.0 + np.cos(3.0 * phi + delta)),
                          -3.0 * k * np.sin(3.0 * phi + delta),
                          -9.0 * k * np.cos(3.0 * phi + delta)))

        energy, gradient, hessian = 0.0, np.zeros(n), np.zeros((n, n))
        for dv, d2v, e, de, d2e in terms:
            energy += e
            gradient += de * dv
            hessian += d2e * np.outer(dv, dv) + de * d2v

        return energy, gradient.reshape(-1, 3), hessian


# ---------------------------------------------------------------------------
# Building datapoints and drivers as the construction does
# ---------------------------------------------------------------------------

# (bond_mode, mass_weighted, exponent_p, exponent_q, tc_weights); the second
# row of each mode is the production default of IMForceFieldGenerator
SETTINGS_CASES = [
    ('inverse', False, 2, 2, False),
    ('inverse', True, 6, 2, True),
    ('plain', False, 2, 2, False),
    ('plain', True, 6, 2, True),
    ('eq', False, 2, 2, False),
    ('eq', True, 6, 2, True),
]
SETTINGS_IDS = ['-'.join(str(v) for v in case) for case in SETTINGS_CASES]


def silent():
    """
    An output stream without output. It also keeps the test objects from
    closing pytest's captured stdout (known issue K45).
    """

    return OutputStream(None)


def im_settings(bond_mode='inverse', p=2, q=2, tc=False):

    return {
        'interpolation_type': 'shepard',
        'weightfunction_type': 'cartesian',
        'exponent_p': p,
        'exponent_q': q,
        'use_inverse_bond_length': bond_mode == 'inverse',
        'use_eq_bond_length': bond_mode == 'eq',
        'use_tc_weights': tc,
    }


def inverse_sqrt_masses(molecule):
    """1/sqrt(m) per Cartesian component, as stored with a mass-weighted database."""

    return 1.0 / np.sqrt(np.repeat(molecule.get_masses(), 3))


def make_datapoint(z_matrix, settings, coords, energy, gradient, hessian,
                   radius=1.0, inv_sqrt_masses=None, eq_bond_lengths=None,
                   important=None):
    """
    A datapoint built as IMForceFieldGenerator builds it: from the Cartesian
    E, G and H (mass-weighted with m^-1/2 if inv_sqrt_masses is given) to
    g_q and H_q.
    """

    datapoint = InterpolationDatapoint(z_matrix, ostream=silent())
    datapoint.update_settings(settings)
    datapoint.eq_bond_lengths = eq_bond_lengths
    datapoint.inv_sqrt_masses = inv_sqrt_masses

    if inv_sqrt_masses is not None:
        gradient = (inv_sqrt_masses * gradient.reshape(-1)).reshape(-1, 3)
        hessian = (inv_sqrt_masses[:, None] * hessian *
                   inv_sqrt_masses[None, :])

    datapoint.cartesian_coordinates = np.array(coords, dtype=float)
    datapoint.energy = energy
    datapoint.gradient = gradient
    datapoint.hessian = hessian
    datapoint.confidence_radius = radius
    datapoint.imp_int_coordinates = (important if important is not None else
                                     NO_IMPORTANT)
    datapoint.transform_gradient_and_hessian()
    return datapoint


def make_driver(z_matrix, settings, datapoints, natoms, inv_sqrt_masses=None,
                eq_bond_lengths=None):

    driver = InterpolationDriver(z_matrix, ostream=silent())
    driver.update_settings(settings)
    # the driver reads entry 3 (the atoms in the weights) and entry 4 (the
    # atoms excluded from the weights)
    atoms = list(range(natoms))
    driver.symmetry_information = [atoms, [], [], atoms, [], [], [], {}, {}, []]
    driver.impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
    driver.impes_coordinate.eq_bond_lengths = eq_bond_lengths
    driver.qm_data_points = list(datapoints)
    return driver


def interpolate(driver, labels, coords):
    """E (Hartree) and grad E ((natoms, 3), Hartree/bohr) of the interpolant."""

    driver.compute(Molecule(labels, coords, 'bohr'))
    return driver.get_energy(), np.array(driver.get_gradient())


def reference_database(xyz, z_matrix, case, n_points=3, scale=0.08, seed=7):
    """
    A small database of n_points datapoints of the reference potential near
    the geometry of xyz, for one of SETTINGS_CASES.

    :return:
        (labels, x0, potential, datapoints, driver), x0 in bohr.
    """

    bond_mode, mass_weighted, p, q, tc = case
    rng = np.random.default_rng(seed)
    molecule = Molecule.read_xyz_string(xyz)
    labels = molecule.get_labels()
    x0 = molecule.get_coordinates_in_bohr()

    settings = im_settings(bond_mode, p, q, tc)
    # the reference potential is quadratic in the bond coordinate of the mode
    potential = ReferencePotential(
        z_matrix, x0, bond_function='r' if bond_mode == 'plain' else '1/r',
        torsions=bool(z_matrix['dihedrals']))
    masses = inverse_sqrt_masses(molecule) if mass_weighted else None
    eq_lengths = (bond_lengths(molecule, z_matrix)
                  if bond_mode == 'eq' else None)
    important = ((ETHANOL_IMPORTANT if z_matrix['dihedrals'] else
                  WATER_IMPORTANT) if tc else None)

    datapoints = []
    for _ in range(n_points):
        x = x0 + rng.normal(scale=scale, size=x0.shape)
        datapoints.append(
            make_datapoint(z_matrix, settings, x, *potential.compute(x),
                           radius=1.0, inv_sqrt_masses=masses,
                           eq_bond_lengths=eq_lengths, important=important))

    driver = make_driver(z_matrix, settings, datapoints, len(labels),
                         inv_sqrt_masses=masses, eq_bond_lengths=eq_lengths)
    return labels, x0, potential, datapoints, driver


# ---------------------------------------------------------------------------
# 1. The reference potential
# ---------------------------------------------------------------------------


class TestReferencePotential:
    """The oracle of sections 3, 6 and 7 must itself be exact."""

    def test_gradient_and_hessian_are_exact(self):

        molecule = Molecule.read_xyz_string(ETHANOL_XYZ)
        x0 = molecule.get_coordinates_in_bohr()
        potential = ReferencePotential(ETHANOL_Z_MATRIX, x0, '1/r')
        x = x0 + np.random.default_rng(1).normal(scale=0.05, size=x0.shape)

        _, gradient, hessian = potential.compute(x)
        fd_gradient = derivative(lambda y: potential.compute(y)[0], x)
        fd_hessian = derivative(lambda y: potential.compute(y)[1], x)

        assert np.max(np.abs(gradient.reshape(-1) - fd_gradient)) < 1.0e-8
        assert np.max(np.abs(hessian - fd_hessian.reshape(hessian.shape))) < 1.0e-8

    def test_invariance_under_rigid_motion(self):

        molecule = Molecule.read_xyz_string(ETHANOL_XYZ)
        x0 = molecule.get_coordinates_in_bohr()
        potential = ReferencePotential(ETHANOL_Z_MATRIX, x0, 'r')
        rng = np.random.default_rng(2)
        x = x0 + rng.normal(scale=0.05, size=x0.shape)
        rotation = random_rotation(rng)

        energy, gradient, _ = potential.compute(x)
        moved_energy, moved_gradient, _ = potential.compute(
            x @ rotation.T + np.array([0.3, -1.2, 2.0]))

        assert abs(moved_energy - energy) < 1.0e-12
        assert np.max(np.abs(moved_gradient - gradient @ rotation.T)) < 1.0e-10
        # no net force and no net torque
        assert np.max(np.abs(gradient.sum(axis=0))) < 1.0e-12
        assert np.max(np.abs(np.cross(x, gradient).sum(axis=0))) < 1.0e-12


# ---------------------------------------------------------------------------
# 2. Internal coordinates and the B and B2 matrices
# ---------------------------------------------------------------------------


class TestInternalCoordinates:
    """Equations (1) and (2) of InterpolationDatapoint."""

    def setup_datapoint(self, bond_mode, mass_weighted):
        """
        A datapoint of ethanol at a displaced geometry x. Returns the
        datapoint, the molecule at x and the r_eq of the eq mode.
        """

        molecule = Molecule.read_xyz_string(ETHANOL_XYZ)
        x0 = molecule.get_coordinates_in_bohr()
        x = x0 + np.random.default_rng(3).normal(scale=0.05, size=x0.shape)
        eq_lengths = (bond_lengths(molecule, ETHANOL_Z_MATRIX)
                      if bond_mode == 'eq' else None)

        datapoint = InterpolationDatapoint(ETHANOL_Z_MATRIX, ostream=silent())
        datapoint.update_settings(im_settings(bond_mode))
        datapoint.eq_bond_lengths = eq_lengths
        datapoint.inv_sqrt_masses = (inverse_sqrt_masses(molecule)
                                     if mass_weighted else None)
        datapoint.reset_coordinates(x)
        return datapoint, moved(molecule, x), eq_lengths

    @pytest.mark.parametrize('bond_mode', ['inverse', 'plain', 'eq'])
    def test_coordinate_values(self, bond_mode):

        datapoint, molecule, eq_lengths = self.setup_datapoint(bond_mode, False)

        expected = coordinate_values(molecule, ETHANOL_Z_MATRIX, bond_mode,
                                     eq_lengths)

        assert np.max(np.abs(datapoint.internal_coordinates_values -
                             expected)) < 1.0e-12

    @pytest.mark.parametrize('bond_mode', ['inverse', 'plain', 'eq'])
    @pytest.mark.parametrize('mass_weighted', [False, True])
    def test_b_matrix_is_derivative_of_coordinates(self, bond_mode,
                                                   mass_weighted):

        datapoint, molecule, eq_lengths = self.setup_datapoint(bond_mode,
                                                               mass_weighted)
        x = molecule.get_coordinates_in_bohr()
        q0 = coordinate_values(molecule, ETHANOL_Z_MATRIX, bond_mode,
                               eq_lengths)

        # wrap the dihedral differences, so that a dihedral near +-pi is
        # differentiated correctly
        def q(y):
            return q0 + wrap_angle(
                coordinate_values(moved(molecule, y), ETHANOL_Z_MATRIX,
                                  bond_mode, eq_lengths) - q0)

        expected = derivative(q, x)
        if mass_weighted:
            # dq/dQ = dq/dx m^-1/2
            expected = expected * datapoint.inv_sqrt_masses[None, :]

        assert np.max(np.abs(datapoint.b_matrix - expected)) < 1.0e-8

    @pytest.mark.parametrize('bond_mode', ['inverse', 'plain', 'eq'])
    @pytest.mark.parametrize('mass_weighted', [False, True])
    def test_b2_matrix_is_derivative_of_b_matrix(self, bond_mode,
                                                 mass_weighted):

        datapoint, molecule, _ = self.setup_datapoint(bond_mode, mass_weighted)
        x = molecule.get_coordinates_in_bohr()
        b2_matrix = datapoint.b2_matrix.copy()

        def b_matrix(y):
            datapoint.reset_coordinates(y)
            return datapoint.b_matrix.copy()

        expected = derivative(b_matrix, x)
        if mass_weighted:
            # d2q/dQ2 = m^-1/2 d2q/dx2 m^-1/2; b_matrix already holds one factor
            expected = expected * datapoint.inv_sqrt_masses[None, None, :]

        assert np.max(np.abs(b2_matrix - expected)) < 1.0e-8


# ---------------------------------------------------------------------------
# 3. Transformation of the gradient and Hessian to internal coordinates
# ---------------------------------------------------------------------------


class TestTransformation:
    """
    Equation (3). For water, the internal coordinates are a complete,
    non-redundant set, so g_q and H_q of the reference potential are known
    exactly: g_q = dV/dq and H_q = d2V/dq2. A missing B2 term in the Hessian
    transformation (mutant M1) shows up here. The geometry is not a
    stationary point, so the B2 term is not zero.
    """

    @pytest.mark.parametrize('bond_mode', ['inverse', 'plain', 'eq'])
    @pytest.mark.parametrize('mass_weighted', [False, True])
    def test_internal_gradient_and_hessian(self, bond_mode, mass_weighted):

        molecule = Molecule.read_xyz_string(WATER_XYZ)
        x0 = molecule.get_coordinates_in_bohr()
        x = x0 + np.random.default_rng(4).normal(scale=0.1, size=x0.shape)
        potential = ReferencePotential(
            WATER_Z_MATRIX, x0, 'r' if bond_mode == 'plain' else '1/r')
        eq_lengths = (bond_lengths(molecule, WATER_Z_MATRIX)
                      if bond_mode == 'eq' else None)

        datapoint = make_datapoint(
            WATER_Z_MATRIX, im_settings(bond_mode), x, *potential.compute(x),
            inv_sqrt_masses=(inverse_sqrt_masses(molecule)
                             if mass_weighted else None),
            eq_bond_lengths=eq_lengths)

        # dV/dq and d2V/dq2 in the coordinate of the mode. V is quadratic in
        # f = r (plain) or f = 1/r (inverse, eq); in eq mode
        # q = r_eq - r_eq^2/r, so df/dq = -1/r_eq^2.
        r = bond_lengths(moved(molecule, x), WATER_Z_MATRIX)
        df_dq = -1.0 / eq_lengths**2 if bond_mode == 'eq' else np.ones(2)
        theta = moved(molecule, x).get_angle([2, 1, 3], 'radian')
        expected_gradient = np.concatenate([
            potential.kb * (potential.f(r) - potential.f0) * df_dq,
            potential.ka * (theta - potential.theta0)])
        expected_hessian = np.diag(np.concatenate([potential.kb * df_dq**2,
                                                   potential.ka]))

        assert np.max(np.abs(datapoint.internal_gradient -
                             expected_gradient)) < 1.0e-10
        assert np.max(np.abs(datapoint.internal_hessian -
                             expected_hessian)) < 1.0e-10


# ---------------------------------------------------------------------------
# 4. The local Taylor model T_i
# ---------------------------------------------------------------------------


def taylor_model(q, q_ref, energy, gradient, hessian, kinds):
    """Equation (4), written from the InterpolationDriver docstring."""

    delta = q - q_ref
    s = delta.copy()
    dihedral, improper = (kinds == 'dihedral'), (kinds == 'improper')
    s[dihedral] = 2.0 * np.sin(0.5 * wrap_angle(delta[dihedral]))
    s[improper] = 2.0 * np.tan(0.5 * wrap_angle(delta[improper]))
    return energy + gradient @ s + 0.5 * s @ hessian @ s


class TestTaylorModel:
    """
    Equations (4) and (5) at large torsion displacements, where
    2 sin(Delta/2), 2 tan(Delta/2) and Delta differ: for Delta = 170 deg,
    2 sin(Delta/2) = 1.99 and Delta = 2.97. Near Delta = 0 all three agree to
    second order, so only large displacements can tell them apart (mutant M3).
    The gradient tests the chain factor c (mutant M5).

    With a single datapoint W = 1 and grad W = 0, so the interpolant is T_i.
    The datapoint gets chosen g_q and H_q, so the test does not depend on the
    transformation of section 3.
    """

    @pytest.mark.parametrize('torsion_section, phi_ref, phi_probe', [
        ('dihedrals', 115.0, 215.0),    # Delta = 100 deg
        ('dihedrals', 115.0, 285.0),    # Delta = 170 deg
        ('dihedrals', 115.0, -64.0),    # Delta = -179 deg
        ('dihedrals', 175.0, -175.0),   # across +-pi: -350 deg wraps to +10 deg
        ('impropers', 115.0, 175.0),    # Delta = 60 deg
        ('impropers', 115.0, 70.0),     # Delta = -45 deg
        ('impropers', 175.0, -175.0),   # across +-pi
    ])
    def test_taylor_model_at_large_torsion_displacement(
            self, torsion_section, phi_ref, phi_probe):

        z_matrix = h2o2_z_matrix(torsion_section)
        kinds = coordinate_kinds(z_matrix)
        settings = im_settings('inverse')
        rng = np.random.default_rng(5)

        reference = h2o2_molecule(phi_ref)
        # the probe also changes two bonds and an angle
        probe = h2o2_molecule(phi_probe, r_oo=2.80, r_oh=1.86,
                              angle_degrees=97.4)
        x_probe = probe.get_coordinates_in_bohr()

        energy = -151.3
        gradient = rng.normal(scale=0.05, size=len(kinds))
        a = rng.normal(scale=0.1, size=(len(kinds), len(kinds)))
        hessian = 0.5 * (a + a.T) + np.eye(len(kinds))

        datapoint = InterpolationDatapoint(z_matrix, ostream=silent())
        datapoint.update_settings(settings)
        datapoint.reset_coordinates(reference.get_coordinates_in_bohr())
        datapoint.energy = energy
        datapoint.internal_gradient = gradient
        datapoint.internal_hessian = hessian
        datapoint.confidence_radius = 1.0
        datapoint.imp_int_coordinates = NO_IMPORTANT

        driver = make_driver(z_matrix, settings, [datapoint], 4)
        msi_energy, msi_gradient = interpolate(driver, probe.get_labels(),
                                               x_probe)

        q_ref = coordinate_values(reference, z_matrix, 'inverse')

        def expected_energy(x):
            return taylor_model(
                coordinate_values(moved(probe, x), z_matrix, 'inverse'), q_ref,
                energy, gradient, hessian, kinds)

        expected_gradient = derivative(expected_energy, x_probe).reshape(-1, 3)

        assert abs(msi_energy - expected_energy(x_probe)) < 1.0e-10
        assert np.max(np.abs(msi_gradient - expected_gradient)) < 1.0e-8


# ---------------------------------------------------------------------------
# 5. The Shepard weights
# ---------------------------------------------------------------------------


class TestShepardWeights:
    """Equation (6), with p != q, so that both terms of u^p + u^q matter."""

    def test_weights_follow_kabsch_distance(self):

        molecule = Molecule.read_xyz_string(WATER_XYZ)
        labels = molecule.get_labels()
        x = molecule.get_coordinates_in_bohr()
        rng = np.random.default_rng(6)
        p, q = 6, 2
        settings = im_settings('inverse', p, q)
        potential = ReferencePotential(WATER_Z_MATRIX, x, '1/r', torsions=False)

        # the datapoints are displaced, rotated and translated copies of x, so
        # the distance is wrong without the Kabsch alignment
        radii = [0.8, 1.3]
        datapoints = []
        for radius in radii:
            y = ((x + rng.normal(scale=0.15, size=x.shape)) @
                 random_rotation(rng).T + rng.normal(size=3))
            datapoints.append(make_datapoint(WATER_Z_MATRIX, settings, y,
                                             *potential.compute(y),
                                             radius=radius))
        driver = make_driver(WATER_Z_MATRIX, settings, datapoints, 3)
        interpolate(driver, labels, x)

        u = np.array([
            kabsch_distance(x, dp.cartesian_coordinates)**2 / dp.confidence_radius**2
            for dp in datapoints
        ])
        w = 1.0 / (u**p + u**q)

        assert np.max(np.abs(driver.raw_weights_array / w - 1.0)) < 1.0e-10
        assert np.max(np.abs(driver.normalized_weights_array -
                             w / w.sum())) < 1.0e-12

    def test_target_customised_gate(self):
        """
        The multiplicative gate exp(-lambda A), A = mean_k z_k^2, with
        z = Delta(1/r) r_ref^2 / sigma_r for an inverse bond and
        z = 2 sin(Delta/2) / (2 sin(sigma/2)) for a dihedral, lambda = 0.7,
        sigma_r = 0.04 Angstrom and sigma = 12 deg (InterpolationDriver
        defaults).
        """

        z_matrix = h2o2_z_matrix('dihedrals')
        p, q = 6, 2
        settings = im_settings('inverse', p, q, tc=True)
        reference = h2o2_molecule(115.0)
        probe = h2o2_molecule(135.0, r_oo=2.80)
        x_ref = reference.get_coordinates_in_bohr()
        x_probe = probe.get_coordinates_in_bohr()
        potential = ReferencePotential(z_matrix, x_ref, '1/r')

        datapoint = make_datapoint(z_matrix, settings, x_ref,
                                   *potential.compute(x_ref), radius=1.0,
                                   important=H2O2_IMPORTANT)
        driver = make_driver(z_matrix, settings, [datapoint], 4)
        interpolate(driver, probe.get_labels(), x_probe)

        u = kabsch_distance(x_probe, x_ref)**2
        r = probe.get_distance([1, 2], 'bohr')
        r_ref = reference.get_distance([1, 2], 'bohr')
        sigma_r = 0.04 / bohr_in_angstrom()
        z_bond = (1.0 / r - 1.0 / r_ref) * r_ref**2 / sigma_r
        delta = wrap_angle(probe.get_dihedral([3, 1, 2, 4], 'radian') -
                           reference.get_dihedral([3, 1, 2, 4], 'radian'))
        z_dihedral = 2.0 * np.sin(0.5 * delta) / (2.0 * np.sin(np.deg2rad(6.0)))
        gate = np.exp(-0.7 * 0.5 * (z_bond**2 + z_dihedral**2))
        expected = gate / (u**p + u**q)

        assert abs(driver.raw_weights_array[0] / expected - 1.0) < 1.0e-10

    @pytest.mark.parametrize('tc', [False, True])
    def test_radius_derivatives(self, tc):
        """
        dw/dR and d(grad w)/dR, which IMTrustRadiusOptimizer uses to fit the
        confidence radii, are the derivatives of w and grad w with respect to
        the radius R.
        """

        z_matrix = h2o2_z_matrix('dihedrals')
        settings = im_settings('inverse', 6, 2, tc=tc)
        reference = h2o2_molecule(115.0)
        probe = h2o2_molecule(135.0, r_oo=2.80)
        x_ref = reference.get_coordinates_in_bohr()
        potential = ReferencePotential(z_matrix, x_ref, '1/r')

        datapoint = make_datapoint(z_matrix, settings, x_ref,
                                   *potential.compute(x_ref), radius=0.9,
                                   important=H2O2_IMPORTANT)
        driver = make_driver(z_matrix, settings, [datapoint], 4)
        interpolate(driver, probe.get_labels(),
                    probe.get_coordinates_in_bohr())
        driver.calc_optim_trust_radius = True

        def weight_terms(radius):
            # cartesian_distance returns
            # (d, -, 1/w, grad w, v, -, dw/dR, d(grad w)/dR)
            datapoint.confidence_radius = radius
            return driver.cartesian_distance(datapoint,
                                             store_alpha_gradients=False)

        terms = weight_terms(0.9)
        analytic = np.concatenate([[terms[6]], np.ravel(terms[7])])

        def weight_and_gradient(radius):
            terms = weight_terms(radius[0])
            return np.concatenate([[1.0 / terms[2]], np.ravel(terms[3])])

        expected = derivative(weight_and_gradient, np.array([0.9]), h=1.0e-4)[:, 0]

        assert np.max(np.abs(analytic - expected)) < 1.0e-8 * np.max(np.abs(expected))


# ---------------------------------------------------------------------------
# 6. The interpolated energy and gradient
# ---------------------------------------------------------------------------


class TestInterpolation:
    """Equation (7) on a redundant z-matrix (ethanol) with three datapoints."""

    @pytest.mark.parametrize('case', SETTINGS_CASES, ids=SETTINGS_IDS)
    def test_gradient_is_derivative_of_energy(self, case):

        labels, x0, _, _, driver = reference_database(ETHANOL_XYZ,
                                                      ETHANOL_Z_MATRIX, case)
        x = x0 + np.random.default_rng(8).normal(scale=0.06, size=x0.shape)

        _, gradient = interpolate(driver, labels, x)
        fd_gradient = derivative(lambda y: interpolate(driver, labels, y)[0], x)

        assert np.max(np.abs(gradient.reshape(-1) - fd_gradient)) < 1.0e-8

    def test_invariance_under_rigid_motion(self):

        case = ('inverse', True, 6, 2, True)
        labels, x0, _, _, driver = reference_database(ETHANOL_XYZ,
                                                      ETHANOL_Z_MATRIX, case)
        rng = np.random.default_rng(9)
        x = x0 + rng.normal(scale=0.06, size=x0.shape)
        rotation = random_rotation(rng)

        energy, gradient = interpolate(driver, labels, x)
        moved_energy, moved_gradient = interpolate(
            driver, labels, x @ rotation.T + np.array([1.0, -0.5, 0.2]))

        assert abs(moved_energy - energy) < 1.0e-10
        assert np.max(np.abs(moved_gradient - gradient @ rotation.T)) < 1.0e-10
        # no net force and no net torque
        assert np.max(np.abs(gradient.sum(axis=0))) < 1.0e-10
        assert np.max(np.abs(np.cross(x, gradient).sum(axis=0))) < 1.0e-10


# ---------------------------------------------------------------------------
# 7. Fidelity: the interpolant reproduces the reference potential
# ---------------------------------------------------------------------------


class TestFidelity:
    """
    What the construction relies on: the interpolant reproduces the data it
    was built from, it is exact where the method must be exact, and its error
    near a datapoint has the order the method promises.
    """

    @pytest.mark.parametrize('case', SETTINGS_CASES, ids=SETTINGS_IDS)
    def test_reproduces_datapoints(self, case):
        """
        At every datapoint, E = E_i, grad E = g_i and the Hessian is H_i.

        The Hessian of the interpolant is the derivative of its gradient. With
        q = 2 the weights of the other datapoints vanish as d^4 at the
        datapoint, so they do not contribute. The step is 2e-3 bohr: the own
        weight of the datapoint grows as 1/d^4, and smaller steps amplify the
        round-off of the normalisation (1e-8 at 1e-4 bohr, 4e-11 at 2e-3).
        """

        labels, _, potential, datapoints, driver = reference_database(
            ETHANOL_XYZ, ETHANOL_Z_MATRIX, case)

        for datapoint in datapoints:
            x = datapoint.cartesian_coordinates
            ref_energy, ref_gradient, ref_hessian = potential.compute(x)

            energy, gradient = interpolate(driver, labels, x)
            hessian = derivative(lambda y: interpolate(driver, labels, y)[1],
                                 x, h=2.0e-3).reshape(ref_hessian.shape)

            assert abs(energy - ref_energy) < 1.0e-10
            assert np.max(np.abs(gradient - ref_gradient)) < 1.0e-10
            assert np.max(np.abs(hessian - ref_hessian)) < 1.0e-8

    @pytest.mark.parametrize('case', SETTINGS_CASES, ids=SETTINGS_IDS)
    def test_exact_for_quadratic_potential(self, case):
        """
        Water: 3 internal coordinates for 3N - 6 = 3, and a reference
        potential quadratic in them. Every T_i is then exactly V, so
        E = sum_i W_i V = V and grad E = grad V everywhere, for any weights.
        """

        labels, x0, potential, _, driver = reference_database(
            WATER_XYZ, WATER_Z_MATRIX, case, scale=0.1)
        rng = np.random.default_rng(10)

        for _ in range(5):
            x = ((x0 + rng.normal(scale=0.15, size=x0.shape)) @
                 random_rotation(rng).T + rng.normal(size=3))
            energy, gradient = interpolate(driver, labels, x)
            ref_energy, ref_gradient, _ = potential.compute(x)

            assert abs(energy - ref_energy) < 1.0e-10
            assert np.max(np.abs(gradient - ref_gradient)) < 1.0e-10

    @pytest.mark.parametrize('case', [
        ('inverse', False, 2, 2, False),
        ('plain', True, 2, 2, False),
        ('eq', False, 2, 2, False),
    ], ids=['inverse', 'plain-mass_weighted', 'eq'])
    def test_error_is_third_order_near_datapoint(self, case):
        """
        Ethanol, one datapoint: |E - V| = O(t^3) and |grad E - grad V| = O(t^2)
        along x_i + t u. Every halving of t must divide the errors by 2^3 and
        2^2; the fitted orders must lie within 0.2 of 3 and 2. A Hessian
        transformation without the B2 term (mutant M1) reduces the energy
        order to 2.
        """

        labels, _, potential, datapoints, driver = reference_database(
            ETHANOL_XYZ, ETHANOL_Z_MATRIX, case, n_points=1)
        x = datapoints[0].cartesian_coordinates
        direction = np.random.default_rng(11).normal(size=x.shape)
        direction /= np.linalg.norm(direction)

        energy_errors, gradient_errors = [], []
        for t in [0.04, 0.02, 0.01, 0.005]:
            energy, gradient = interpolate(driver, labels, x + t * direction)
            ref_energy, ref_gradient, _ = potential.compute(x + t * direction)
            energy_errors.append(abs(energy - ref_energy))
            gradient_errors.append(np.max(np.abs(gradient - ref_gradient)))

        energy_orders = np.log2(np.array(energy_errors[:-1]) / energy_errors[1:])
        gradient_orders = np.log2(
            np.array(gradient_errors[:-1]) / gradient_errors[1:])

        assert np.all(np.abs(energy_orders - 3.0) < 0.2)
        assert np.all(np.abs(gradient_orders - 2.0) < 0.2)


# ---------------------------------------------------------------------------
# 8. Force units in the dynamics
# ---------------------------------------------------------------------------

# 1 Eh/bohr in kJ/mol/nm (CODATA 2018: Eh = 2625.4996394798 kJ/mol,
# a0 = 0.0529177210903 nm); independent of the veloxchem constants
HARTREE_PER_BOHR_IN_KJ_PER_MOL_NM = 2625.4996394798 / 0.0529177210903

FIXED_GRADIENT = np.array([[0.010, -0.020, 0.030],
                           [-0.005, 0.004, 0.001],
                           [-0.005, 0.016, -0.031]])


class _FixedGradientDynamics(OpenMMIMDynamics):
    """
    OpenMMIMDynamics whose interpolated gradient is FIXED_GRADIENT
    (Hartree/bohr). It records the positions passed to update_gradient.
    """

    def __init__(self):
        super().__init__(ostream=silent())
        self.received_positions = []

    def update_gradient(self, new_positions):
        self.received_positions.append(np.array(new_positions))
        return FIXED_GRADIENT


class _FixedGradientCollector(IMDatabasePointCollecter):
    """
    IMDatabasePointCollecter whose interpolated gradient is FIXED_GRADIENT
    (Hartree/bohr). It records the positions passed to update_gradient.
    """

    def __init__(self):
        super().__init__(ostream=silent())
        self.received_positions = []

    def update_gradient(self, new_positions):
        self.received_positions.append(np.array(new_positions))
        return FIXED_GRADIENT


class TestForceUnits:
    """
    Equation (8). The QM forces of both MD drivers go through a
    CustomExternalForce whose per-particle parameters are the force in
    kJ/mol/nm. The drivers get a fixed gradient (the classes above), and the
    force OpenMM reports is compared with
    -FIXED_GRADIENT * HARTREE_PER_BOHR_IN_KJ_PER_MOL_NM.
    """

    positions_nm = np.array([[0.0, 0.0, 0.0117],
                             [0.0, 0.0757, -0.0467],
                             [0.0, -0.0757, -0.0467]])

    def openmm_system(self):
        """
        A water topology and a system with the QM force of the MD drivers,
        and its context on the Reference platform.
        """

        mm = pytest.importorskip('openmm')
        app = pytest.importorskip('openmm.app')
        unit = pytest.importorskip('openmm.unit')

        topology = app.Topology()
        residue = topology.addResidue('MOL', topology.addChain())
        for symbol in ('O', 'H', 'H'):
            topology.addAtom(symbol, app.Element.getBySymbol(symbol), residue)

        system = mm.System()
        for mass in (16.0, 1.0, 1.0):
            system.addParticle(mass)
        qm_force = mm.CustomExternalForce('-fx*x-fy*y-fz*z')
        for name in ('fx', 'fy', 'fz'):
            qm_force.addPerParticleParameter(name)
        for i in range(3):
            qm_force.addParticle(i, [0.0, 0.0, 0.0])
        system.addForce(qm_force)

        context = mm.Context(system, mm.VerletIntegrator(0.001),
                             mm.Platform.getPlatformByName('Reference'))
        context.setPositions(self.positions_nm * unit.nanometer)
        return topology, system, context

    def openmm_forces(self, context):

        unit = pytest.importorskip('openmm.unit')
        state = context.getState(getForces=True)
        return state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer)

    def test_dynamics_force_conversion(self):

        _, system, context = self.openmm_system()

        dynamics = _FixedGradientDynamics()
        dynamics.system = system
        dynamics.qm_force_index = 0
        dynamics.qm_atoms = [0, 1, 2]

        dynamics.update_forces(context)

        expected = -FIXED_GRADIENT * HARTREE_PER_BOHR_IN_KJ_PER_MOL_NM
        forces = self.openmm_forces(context)

        # the interpolation receives the positions in nm
        assert np.max(np.abs(dynamics.received_positions[0] -
                             self.positions_nm)) < 1.0e-12
        assert np.max(np.abs(forces - expected)) / np.max(np.abs(expected)) < 1.0e-10

    def test_collector_force_conversion(self):

        topology, system, context = self.openmm_system()

        collector = _FixedGradientCollector()
        collector.topology = topology
        collector.system = system
        collector.qm_force_index = 0
        collector.qm_atoms = [0, 1, 2]
        # the charge and multiplicity of the QM molecule of each step
        collector.molecule = Molecule.read_xyz_string(WATER_XYZ)
        # per-step bookkeeping of the collector, set so that this step only
        # applies the forces: skipping_value = 1 skips the QM-check gate, and
        # point_checker stays away from 0, which would reset the trajectory
        collector.skipping_value = 1
        collector.point_checker = 1
        collector.duration = 10**9
        collector.last_point_added = 0
        collector.cycle_iteration = 0
        collector.unadded_cycles = 0
        collector.current_state = 0
        collector.state_specific_molecules = {0: []}

        collector.update_forces(context)

        expected = -FIXED_GRADIENT * HARTREE_PER_BOHR_IN_KJ_PER_MOL_NM
        forces = self.openmm_forces(context)
        step_molecule = collector.state_specific_molecules[0][0]

        assert np.max(np.abs(collector.received_positions[0] -
                             self.positions_nm)) < 1.0e-12
        # the QM molecule of the step is built in Angstrom
        assert np.max(np.abs(step_molecule.get_coordinates_in_angstrom() -
                             10.0 * self.positions_nm)) < 1.0e-12
        assert np.max(np.abs(forces - expected)) / np.max(np.abs(expected)) < 1.0e-10
