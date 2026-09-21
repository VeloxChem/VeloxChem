"""The range separated RI-JK path, against the four center way.

    python rs_rijk_check.py

Checks the whole path a hybrid range separated functional takes through the simd
resolution of the identity: the two metrics, the two sets of B vectors, the Fock
matrices they build for a closed and an open shell, the limits of omega, the
refusals, and a converged SCF energy. Prints a line per case and exits non-zero if
anything is outside tolerance. It runs in a few minutes.

It belongs here rather than in `tests` for the reason
`rs_integral_check.py` does: the simd drivers have no test coverage, and what this
compares against -- the four center driver -- is slow enough that a suite would not
want it. See that file for the integrals themselves, which this takes as given.

**What the comparison can and cannot say.** The resolution of the identity does not
reproduce four center integrals exactly and is not meant to; a fitting set carries an
error of its own. So a fixed tolerance would be measuring the basis and not the code.
Every range separated case here is therefore barred against the **plain hybrid of the
same molecule and basis**, B3LYP or a bare exchange fraction, computed in the same
run. The question that has an answer is whether fitting the attenuated operator is as
good as fitting the plain one, and it is: every ratio below is between 0.98 and 1.35.

**The composition.** `_get_2e_fock_build_params` gives `exchange_scaling_factor` as
alpha + beta and `erf_k_coef` as -beta, and the Fock matrix is

    F = 2 J[D] - (alpha + beta) K[C] - erf_k_coef K_erf[C]

for a closed shell, and the same per spin with an undoubled Coulomb of the total
density for an open one. For CAM-B3LYP that is -0.19 K - 0.46 K_erf, which is 0.19 of
exact exchange at short range and 0.65 at long. Nothing in the drivers knows that;
they take two factors and subtract two exchanges.

**The two ends of omega.** At a large omega the attenuated operator becomes 1/r, so
the build must collapse onto a plain hybrid whose exchange fraction is alpha alone --
the beta terms cancel. That check is here. The other end is not: `prepare` refuses a
vanishing omega, whose attenuated metric is the zero matrix, and no functional asks
for one.
"""
import subprocess
import sys
import textwrap

import numpy as np
import veloxchem as vlx
from veloxchem import FockDriver, MolecularBasis, T4CScreener
from veloxchem.veloxchemlib import (PackedMatrix, SimdRIFockDriver,
                                    SimdRIJKFockDriver, make_matrix, mat_t,
                                    parse_xc_func, rimode)

WATER = """3

O   0.0000  0.0000  0.1173
H   0.0000  0.7572 -0.4692
H   0.0000 -0.7572 -0.4692
"""

CAFFEINE = "/Users/rinkevic/Development/VeloxChem/benchmarks/geometries/caffeine.xyz"

AUX = "def2-universal-jkfit"
THRESHOLD = 1.0e-14
BUDGET = 8 * 1024**3
OMEGAS = (0.2, 0.33, 3.0)
FUNCTIONALS = ("CAM-B3LYP", "LRC-WPBEH", "WB97X-D4")


def geometry(name, charge=0, multiplicity=1):
    mol = (vlx.Molecule.read_xyz_string(WATER) if name == "water"
           else vlx.Molecule.read_xyz_file(CAFFEINE))
    mol.set_charge(charge)
    mol.set_multiplicity(multiplicity)
    return mol


def parameters(name):
    """What the scf driver passes a build, for one functional."""
    xc = parse_xc_func(name)

    if xc.is_range_separated():
        return (xc.get_rs_alpha() + xc.get_rs_beta(), -xc.get_rs_beta(),
                xc.get_rs_omega(), xc.get_rs_alpha())

    return xc.get_frac_exact_exchange(), 0.0, 0.0, None


class Case:
    """One molecule and basis, with everything the checks below share.

    The density is taken from a converged calculation of the resolution of the
    identity rather than of the four center integrals. It is only a density to build
    from, the same one on both sides of every comparison, so which way converged it
    does not enter -- and the four center way would spend ten minutes on caffeine to
    reach a density which serves no better.
    """

    def __init__(self, name, orbital, charge=0, multiplicity=1):
        self.label = f"{name} {orbital}"
        self.mol = geometry(name, charge, multiplicity)
        self.bas = MolecularBasis.read(self.mol, orbital.upper(), ostream=None)
        self.aux = MolecularBasis.read(self.mol, AUX.upper(), ostream=None)
        self.nao = self.bas.get_dimensions_of_basis()

        scf = vlx.ScfRestrictedDriver()
        scf.ostream.mute()
        scf.ri_jk, scf.ri_jk_simd = True, True
        scf.ri_auxiliary_basis, scf.ri_mode = AUX, "in_memory"
        scf.compute(self.mol, self.bas)

        nocc = self.mol.number_of_alpha_electrons()
        orbitals = scf.mol_orbs.alpha_to_numpy()

        # NOTE: a second set of coefficients one orbital short, which is a cation's
        # beta occupation. The open shell build is handed two sets of different
        # widths, and a check which handed it the same set twice would not tell a
        # spin loop which ran twice on alpha from one which ran on each.
        self.ca = np.ascontiguousarray(orbitals[:, :nocc])
        self.cb = np.ascontiguousarray(orbitals[:, :nocc - 1])

        self.da, self.db = self.ca @ self.ca.T, self.cb @ self.cb.T
        self.dm = self.da

        self.coeff_a = PackedMatrix(self.nao, nocc, mat_t.general)
        self.coeff_a.from_numpy(self.ca)
        self.coeff_b = PackedMatrix(self.nao, nocc - 1, mat_t.general)
        self.coeff_b.from_numpy(self.cb)

        self.packed = PackedMatrix(self.nao, self.nao, mat_t.symmetric)
        self.packed.from_numpy(np.ascontiguousarray(self.dm))
        self.packed_ab = PackedMatrix(self.nao, self.nao, mat_t.symmetric)
        self.packed_ab.from_numpy(np.ascontiguousarray(self.da + self.db))

        self.den = self._matrix(self.dm)
        self.den_a, self.den_b = self._matrix(self.da), self._matrix(self.db)
        self.den_ab = self._matrix(self.da + self.db)

        self.screener = T4CScreener()
        self.screener.partition(self.bas, self.mol, "eri")
        self.fock_drv = FockDriver()

    def _matrix(self, values):
        m = make_matrix(self.bas, mat_t.symmetric)
        m.set_values(values)
        return m

    def four_centre(self, kind, density, factor, omega=0.0):
        return self.fock_drv.compute(self.screener, density, kind, factor, omega,
                                     14).to_numpy()

    def prepared(self, omega):
        drv = SimdRIJKFockDriver()
        drv.prepare(self.mol, self.bas, self.aux, THRESHOLD, BUDGET, 1.0e-12, False,
                    rimode.in_memory, [], PackedMatrix(), 1, omega, PackedMatrix())
        return drv


def relative(got, want):
    return float(np.abs(got - want).max()) / max(float(np.abs(want).max()), 1.0)


def check_metrics(case, fails):
    """Each inverted metric against the inverse formed in numpy.

    Restricted to the directions above a relative cut. The attenuated metric is
    numerically singular by construction -- the transform of erf(omega r)/r carries a
    Gaussian factor -- so a whole-space comparison reports 1e-4 to 1e+6 on a metric
    which serves the exchange perfectly, and what it is reporting is a near null space
    the attenuated integrals are themselves zero in.
    """
    from veloxchem.veloxchemlib import (SimdTwoCenterElectronRepulsionDriver,
                                        SimdTwoCenterElectronRepulsionRsDriver)

    drv = SimdRIJKFockDriver()
    plain_only, _ = drv.make_metric(case.mol, case.aux, 1.0e-12, False,
                                    rimode.in_memory)
    plain_only = plain_only.to_numpy()

    v = SimdTwoCenterElectronRepulsionDriver().compute(case.mol, case.aux).to_numpy()

    for omega in OMEGAS:
        m, m_erf = drv.make_metric_rs(case.mol, case.aux, 1.0e-12, False,
                                      rimode.in_memory, omega)
        m, m_erf = m.to_numpy(), m_erf.to_numpy()

        v_c, v_erf = SimdTwoCenterElectronRepulsionRsDriver().compute(
            case.mol, case.aux, omega)
        v_c, v_erf = v_c.to_numpy(), v_erf.to_numpy()

        def inverted(metric, target, cut=1.0e-8):
            w, u = np.linalg.eigh(target)
            keep = w > cut * w[-1]
            got = u[:, keep].T @ metric.T @ metric @ u[:, keep]
            want = np.diag(1.0 / w[keep])
            return float(np.abs(got - want).max()) / float(np.abs(want).max())

        same = float(np.abs(m - plain_only).max())
        operator = float(np.abs(v_c - v).max())
        plain_ok, erf_ok = inverted(m, v), inverted(m_erf, v_erf)

        ok = (same == 0.0 and operator == 0.0 and plain_ok < 1e-7 and erf_ok < 1e-7)
        if not ok:
            fails.append(f"metrics {case.label} omega {omega}")

        print(f"  {'ok  ' if ok else 'FAIL'} metrics      {case.label:20s} omega {omega:<5g} "
              f" plain half same {same:.1e}  inverted {plain_ok:.1e} / {erf_ok:.1e}")


def check_b_vectors(case, fails):
    """Each exchange assembled from its own B vectors, against four centres."""
    ri = SimdRIFockDriver()
    jk = SimdRIJKFockDriver()
    naux = case.aux.get_dimensions_of_basis()

    def exchange(bq):
        w = ri.compute_w_vectors(bq, case.bas, case.aux, case.coeff_a, 0, naux)
        k = PackedMatrix(case.nao, case.nao, mat_t.symmetric)
        k.zero()
        ri.compute_exchange_matrix(w, k, 1.0)
        return k.to_numpy()

    exact = case.four_centre("kx", case.den, 1.0)

    metric, _ = jk.make_metric(case.mol, case.aux, 1.0e-12, False, rimode.in_memory)
    bq_plain = ri.compute_bq_vectors(case.mol, case.bas, case.aux, metric, THRESHOLD)
    control = relative(exchange(bq_plain), exact)

    for omega in OMEGAS:
        m, m_erf = jk.make_metric_rs(case.mol, case.aux, 1.0e-12, False,
                                     rimode.in_memory, omega)
        bq, bq_erf = ri.compute_bq_vectors_rs(case.mol, case.bas, case.aux, m, m_erf,
                                              THRESHOLD, omega)

        same = max(
            float(np.abs(np.asarray(bq.block_to_numpy(i))
                         - np.asarray(bq_plain.block_to_numpy(i))).max())
            for i in range(bq.number_of_blocks()))

        got = relative(exchange(bq_erf),
                       case.four_centre("kx_rs", case.den, 1.0, omega))

        ok = (same == 0.0 and got < 2.0 * control + 1.0e-9)
        if not ok:
            fails.append(f"B vectors {case.label} omega {omega}")

        print(f"  {'ok  ' if ok else 'FAIL'} B vectors    {case.label:20s} omega {omega:<5g} "
              f" plain half same {same:.1e}  K {control:.2e}  K_erf {got:.2e}")


def check_fock(case, fails):
    """Whole Fock matrices, closed and open shell, against four centres."""
    control = {}

    for name in ("B3LYP",) + FUNCTIONALS:
        a_x, erf_k_coef, omega, _ = parameters(name)

        drv = case.prepared(omega)

        closed = relative(
            drv.compute(case.packed, case.coeff_a, a_x, erf_k_coef).to_numpy(),
            case.four_centre("2jkx", case.den, a_x)
            - (case.four_centre("kx_rs", case.den, erf_k_coef, omega)
               if omega > 0.0 else 0.0))

        got_a, got_b = drv.compute(case.packed_ab, case.coeff_a, case.coeff_b, a_x,
                                   erf_k_coef)

        j = case.four_centre("j", case.den_ab, 0.0)
        want = []
        for density in (case.den_a, case.den_b):
            f = j - case.four_centre("kx", density, a_x)
            if omega > 0.0:
                f -= case.four_centre("kx_rs", density, erf_k_coef, omega)
            want.append(f)

        opened = max(relative(got_a.to_numpy(), want[0]),
                     relative(got_b.to_numpy(), want[1]))

        if name == "B3LYP":
            control["closed"], control["open"] = closed, opened
            ok = True
            note = "control"
        else:
            ok = (closed < 3.0 * control["closed"] + 1e-9 and
                  opened < 3.0 * control["open"] + 1e-9)
            note = f"{closed / control['closed']:.2f} x control"

        if not ok:
            fails.append(f"fock {case.label} {name}")

        print(f"  {'ok  ' if ok else 'FAIL'} fock         {case.label:20s} {name:10s} "
              f" closed {closed:.2e}  open {opened:.2e}  ({note})")


def check_large_omega(case, fails):
    """At a large omega the beta terms cancel and a plain hybrid is left.

    erf(omega r)/r goes to 1/r, so K_erf goes to K and
    -(alpha + beta) K - (-beta) K_erf becomes -alpha K. The right hand side is built
    from the plain B vectors alone, so this says the attenuated half of the path
    lands exactly where the plain half already is.
    """
    for name in FUNCTIONALS:
        a_x, erf_k_coef, _, alpha = parameters(name)

        got = case.prepared(1.0e6).compute(case.packed, case.coeff_a, a_x,
                                           erf_k_coef).to_numpy()

        want = case.prepared(0.0).compute(case.packed, case.coeff_a, alpha).to_numpy()

        apart = relative(got, want)

        ok = apart < 1.0e-9
        if not ok:
            fails.append(f"large omega {case.label} {name}")

        print(f"  {'ok  ' if ok else 'FAIL'} large omega  {case.label:20s} {name:10s} "
              f" against a plain hybrid at a_x {alpha:.3f}: {apart:.2e}")


def check_scf(fails):
    """A converged energy, closed shell and open, against the four center way."""
    control = {}

    for charge, multiplicity in ((0, 1), (1, 2)):
        mol = geometry("water", charge, multiplicity)
        bas = MolecularBasis.read(mol, "DEF2-SVP", ostream=None)

        for name in ("B3LYP",) + FUNCTIONALS:
            energies = {}
            for label, ri in (("four centre", False), ("RI-JK simd", True)):
                scf = (vlx.ScfRestrictedDriver() if multiplicity == 1
                       else vlx.ScfUnrestrictedDriver())
                scf.ostream.mute()
                scf.xcfun, scf.conv_thresh, scf.max_iter = name, 1.0e-8, 150
                if ri:
                    scf.ri_jk, scf.ri_jk_simd = True, True
                    scf.ri_auxiliary_basis, scf.ri_mode = AUX, "in_memory"
                result = scf.compute(mol, bas)
                energies[label] = result["scf_energy"] if result else None

            if None in energies.values():
                fails.append(f"scf water q{charge} {name} did not converge")
                print(f"  FAIL scf          water q{charge} m{multiplicity} {name:10s} "
                      " did not converge")
                continue

            delta = energies["RI-JK simd"] - energies["four centre"]

            if name == "B3LYP":
                control[charge] = abs(delta)
                ok, note = True, "control"
            else:
                ok = abs(delta) < 3.0 * control[charge] + 1.0e-8
                note = f"{abs(delta) / control[charge]:.2f} x control"

            if not ok:
                fails.append(f"scf water q{charge} {name}")

            print(f"  {'ok  ' if ok else 'FAIL'} scf          water q{charge} m{multiplicity} "
                  f"{name:10s}  four centre {energies['four centre']:.8f}  "
                  f"delta {delta:+.2e}  ({note})")


def check_refusals(fails):
    """What the path must not do quietly, each in a process of its own.

    A critical assertion aborts rather than raising, so a refusal cannot be caught
    in the calling process and is checked by requiring that a child dies saying it.
    """
    setup = textwrap.dedent(f"""
        import numpy as np
        import veloxchem as vlx
        from veloxchem import MolecularBasis
        from veloxchem.veloxchemlib import (SimdRIJKFockDriver, PackedMatrix,
                                            rimode, mat_t)
        mol = vlx.Molecule.read_xyz_string({WATER!r})
        bas = MolecularBasis.read(mol, "DEF2-SVP", ostream=None)
        aux = MolecularBasis.read(mol, "DEF2-UNIVERSAL-JKFIT", ostream=None)
        d = SimdRIJKFockDriver()
        drv = SimdRIJKFockDriver()
        BIG = 8 * 1024**3
    """)

    cases = {
        "the direct way asked for with an omega":
            "drv.prepare(mol, bas, aux, 1e-12, BIG, 1e-12, False, rimode.direct,"
            " [], PackedMatrix(), 1, 0.33)",
        "the automatic way falling to direct":
            "drv.prepare(mol, bas, aux, 1e-12, 1000, 1e-12, False, rimode.automatic,"
            " [], PackedMatrix(), 1, 0.33)",
        "an attenuated metric without an omega":
            "drv.prepare(mol, bas, aux, 1e-12, BIG, 1e-12, False, rimode.in_memory,"
            " [], PackedMatrix(), 1, 0.0,"
            " d.make_metric_rs(mol, aux, 1e-12, False, rimode.in_memory, 0.33)[1])",
        "a negative omega":
            "drv.prepare(mol, bas, aux, 1e-12, BIG, 1e-12, False, rimode.in_memory,"
            " [], PackedMatrix(), 1, -0.3)",
        "an attenuated exchange from a plain driver":
            "drv.prepare(mol, bas, aux, 1e-12, BIG, 1e-12, False, rimode.in_memory);"
            "n = bas.get_dimensions_of_basis();"
            "p = PackedMatrix(n, n, mat_t.symmetric); p.from_numpy(np.eye(n));"
            "c = PackedMatrix(n, 1, mat_t.general);"
            "c.from_numpy(np.ascontiguousarray(np.eye(n)[:, :1]));"
            "drv.compute(p, c, 0.2, 0.46)",
        "the metrics of a vanishing omega":
            "d.make_metric_rs(mol, aux, 1e-12, False, rimode.in_memory, 0.0)",
    }

    for what, call in cases.items():
        out = subprocess.run([sys.executable, "-c", setup + call],
                             capture_output=True, text=True)

        said = [l.strip() for l in (out.stdout + out.stderr).splitlines()
                if l.strip().startswith("RIJKFockDriver:")]

        ok = (out.returncode != 0) and bool(said)
        if not ok:
            fails.append(f"refusal: {what}")

        print(f"  {'ok  ' if ok else 'FAIL'} refusal      {what:42s} "
              f"{said[-1][15:85] if said else 'was not refused'}")


def check_memory(fails):
    """The figure the budget is checked against is the figure that will be held."""
    case_mol = geometry("water")
    bas = MolecularBasis.read(case_mol, "DEF2-SVP", ostream=None)
    aux = MolecularBasis.read(case_mol, AUX.upper(), ostream=None)

    drv = SimdRIJKFockDriver()
    one = drv.required_memory(case_mol, bas, aux, 1.0e-12)
    two = drv.required_memory(case_mol, bas, aux, 1.0e-12, [], True)

    ok = (two == 2 * one)
    if not ok:
        fails.append("memory does not double")

    print(f"  {'ok  ' if ok else 'FAIL'} memory       two sets against one: "
          f"{one / 1024**2:.2f} MB -> {two / 1024**2:.2f} MB")


def main():
    fails = []

    cases = [Case("water", "def2-svp"), Case("water", "def2-tzvp"),
             Case("caffeine", "def2-svp")]

    for case in cases:
        check_metrics(case, fails)

    for case in cases:
        check_b_vectors(case, fails)

    for case in cases:
        check_fock(case, fails)

    for case in cases:
        check_large_omega(case, fails)

    check_memory(fails)
    check_refusals(fails)
    check_scf(fails)

    print("\nall checks passed" if not fails else "\nFAILURES:")
    for failure in fails:
        print("  " + failure)

    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
