"""The range separated two-center Coulomb gradient kernels, against pyscf.

    python rs_gradient_check.py

The kernels compute the derivative with respect to the center of the bra of both
`(l|1/r12|l')` and `(l|erf(omega r12)/r12|l')` in one call, and the driver contracts
each with its own weighted density. This checks the Coulomb half against the
unattenuated gradient driver and the attenuated half against pyscf, which is the only
outside reference there is for it.

**The Coulomb half is not bit-identical to the unattenuated driver, and should not be
expected to be.** It is for the *integrals*, where the range separated kernel computes
the same values and they are copied out. It is not for the *gradient*, because the
driver accumulates two operators into two partials in one loop where the unattenuated
one accumulates into a single partial, and a sum re-associated is a sum rounded
differently. The residual measured is two to four times the machine epsilon relative,
and it appears with the same weighted density in both slots, so it is the loop and not
the weights. The bar below is a few epsilon rather than zero.

**Why the weighted density is the metric itself.** The obvious check -- a random
Omega contracted with both codes' derivatives -- does not work, because the two codes
do not order the auxiliary functions the same way and a matrix written in one
ordering means something else in the other. Letting each code build Omega from its
**own** integrals removes the question: the contraction is then

    sum over P and Q of V_PQ dV_PQ/dA  =  half the derivative of the squared
                                          Frobenius norm of V

which no relabelling of the functions can move. A first attempt with a random Omega
gave ratios of -0.75, -1.25 and -7.46 between the two codes, which looks like a
broken kernel and is a broken comparison.

**The convention**, calibrated against the unattenuated driver where both sides are
known good: what the driver returns is

    -2 * sum over P on the atom and all Q of  Omega_PQ * pyscf's int2c2e_ip1

to twelve figures. pyscf's ip1 is the derivative of the bra function, which is the
negative of the derivative with respect to the atom it sits on; the factor of two is
the pair counted from both of its centers.
"""
import sys

import numpy as np
import veloxchem as vlx
from pyscf import df, gto

from veloxchem.veloxchemlib import (PackedMatrix,
                                    SimdTwoCenterElectronRepulsionDriver,
                                    SimdTwoCenterElectronRepulsionGradientDriver,
                                    SimdTwoCenterElectronRepulsionGradientRsDriver,
                                    SimdTwoCenterElectronRepulsionRsDriver, mat_t)

OMEGAS = (0.2, 0.33, 3.0)
TOLERANCE = 1.0e-9

# NOTE: a few machine epsilon, for the reason the docstring gives.
SAME_OPERATOR = 1.0e-13

# NOTE: molecules of distinct elements, and one of three atoms so that a pair of
# atoms is not the whole molecule and the gradient of an atom takes contributions
# from more than one pair.
CASES = [(("H", (0.0, 0.0, 0.0)), ("F", (0.0, 0.0, 0.917))),
         (("Li", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, 1.595))),
         (("C", (0.0, 0.0, 0.0)), ("O", (0.0, 0.0, 1.128))),
         (("O", (0.0, 0.0, 0.1173)), ("H", (0.0, 0.7572, -0.4692)),
          ("H", (0.0, -0.7572, -0.4692)))]

FITTING = "def2-universal-jkfit"

fails = []


def as_xyz(geometry):
    return f"{len(geometry)}\n\n" + "\n".join(
        f"{a}  {x} {y} {z}" for a, (x, y, z) in geometry)


def packed(values):
    n = values.shape[0]
    m = PackedMatrix(n, n, mat_t.symmetric)
    m.from_numpy(np.ascontiguousarray(values))
    return m


def reference(auxmol, weights, natoms, omega=None):
    """What pyscf gives for the same contraction, in pyscf's own ordering."""
    if omega is None:
        ip1 = auxmol.intor("int2c2e_ip1")
    else:
        with auxmol.with_range_coulomb(omega):
            ip1 = auxmol.intor("int2c2e_ip1")

    slices = auxmol.aoslice_by_atom()

    out = np.zeros((natoms, 3))

    for atom in range(natoms):
        first, last = slices[atom][2], slices[atom][3]
        for c in range(3):
            out[atom, c] = -2.0 * np.einsum("pq,pq->", weights[first:last, :],
                                            ip1[c, first:last, :])

    return out


def relative(got, want):
    return float(np.abs(got - want).max()) / max(float(np.abs(want).max()), 1.0)


for geometry in CASES:
    label = "".join(a for a, _ in geometry)

    mol = vlx.Molecule.read_xyz_string(as_xyz(geometry))
    aux = vlx.MolecularBasis.read(mol, FITTING.upper(), ostream=None)
    natoms = len(geometry)

    pm = gto.M(atom=[[a, xyz] for a, xyz in geometry], basis="def2-svp",
               unit="Angstrom")
    am = df.addons.make_auxmol(pm, FITTING)

    v_plain = SimdTwoCenterElectronRepulsionDriver().compute(mol, aux).to_numpy()
    v_pyscf = am.intor("int2c2e")

    plain_grad = np.array(
        SimdTwoCenterElectronRepulsionGradientDriver().compute(
            mol, aux, packed(v_plain)).to_numpy())

    plain_vs_pyscf = relative(plain_grad, reference(am, v_pyscf, natoms))

    for omega in OMEGAS:
        coulomb, attenuated = SimdTwoCenterElectronRepulsionRsDriver().compute(
            mol, aux, omega)
        coulomb, attenuated = coulomb.to_numpy(), attenuated.to_numpy()

        got_c, got_a = SimdTwoCenterElectronRepulsionGradientRsDriver().compute(
            mol, aux, packed(v_plain), packed(attenuated), omega)
        got_c, got_a = np.array(got_c.to_numpy()), np.array(got_a.to_numpy())

        # the Coulomb half is the same operator by the same kernels, so it differs
        # from the unattenuated driver only by how the sum was accumulated
        exact = relative(got_c, plain_grad)

        with am.with_range_coulomb(omega):
            v_erf_pyscf = am.intor("int2c2e")

        att_vs_pyscf = relative(got_a, reference(am, v_erf_pyscf, natoms, omega))

        ok = (exact < SAME_OPERATOR and plain_vs_pyscf < TOLERANCE
              and att_vs_pyscf < TOLERANCE)
        if not ok:
            fails.append(f"{label} omega {omega}")

        print(f"  {'ok  ' if ok else 'FAIL'} {label:4s} omega {omega:<5g}  "
              f"coulomb half vs plain driver {exact:.1e}   plain vs pyscf {plain_vs_pyscf:.2e}   "
              f"attenuated vs pyscf {att_vs_pyscf:.2e}")

# NOTE: at a large omega the attenuated operator is 1/r and its derivative is the
# Coulomb one. Both sides are veloxchem here, so nothing of the fitting or the
# ordering enters and the two must agree to rounding.
mol = vlx.Molecule.read_xyz_string(as_xyz(CASES[2]))
aux = vlx.MolecularBasis.read(mol, FITTING.upper(), ostream=None)
v = SimdTwoCenterElectronRepulsionDriver().compute(mol, aux).to_numpy()
got_c, got_a = SimdTwoCenterElectronRepulsionGradientRsDriver().compute(
    mol, aux, packed(v), packed(v), 1.0e6)
apart = float(np.abs(np.array(got_a.to_numpy()) - np.array(got_c.to_numpy())).max())
apart /= max(float(np.abs(np.array(got_c.to_numpy())).max()), 1.0)

ok = apart < 1.0e-10
if not ok:
    fails.append("the large omega limit")
print(f"\n  {'ok  ' if ok else 'FAIL'} at a large omega the attenuated gradient is the "
      f"Coulomb one: {apart:.2e}")

print("\nall checks passed" if not fails else "\nFAILURES: " + ", ".join(fails))
sys.exit(1 if fails else 0)
