"""The range separated three-center Coulomb gradient kernels, against pyscf.

    python rs_t3c_gradient_check.py

The kernels compute the derivative with respect to the two centers on bra side of
both `(ab|1/r12|c)` and `(ab|erf(omega r12)/r12|c)` in one call, six components an
element, and the driver returns one tensor per operator on one sparsity pattern.

Three things are checked, and it is worth saying what each of them can and cannot
catch, because no one of them is sufficient.

**The Coulomb half against the unattenuated gradient driver, element for element.**
This is exact -- the driver copies values rather than contracting them, so unlike the
two-center gradient there is no re-association to round differently. It settles which
three of the six blocks belong to which operator, and it settles the assignment of
d/dA against d/dB for the Coulomb half. It says nothing about the attenuated half.

**The sum of squares of every component against pyscf.** This is a real cross-code
comparison of the attenuated derivative, and it needs nothing of the tensor's
internal layout: squaring and summing over every element and all six components is
invariant both under the relabelling of the basis functions, which the two codes do
differently, and under any permutation of the components within a combination. That
last invariance is also its limitation: **it cannot catch d/dA swapped with d/dB.**
What rules that out is the Coulomb half above, the two operators being written by one
kernel with one layout.

**The large omega limit.** The attenuated derivative must become the Coulomb one
element for element. Both sides are veloxchem and the layouts are identical, so this
is a straight comparison and it exercises the attenuated kernels rather than the
plain ones.

A fourth check was attempted and abandoned: contracting the derivative with the
integrals to give a per-component scalar. The six components are laid out **per
combination of basis functions**, not as six halves of a block, so a flat reshape
mixes them -- which showed up as a non-zero x and y for a molecule lying along z, a
quantity symmetry forbids. Element-wise comparison per component is possible through
`block.element_offset` and the basis functions of each atom basis, and would be the
way to catch a d/dA and d/dB swap directly.
"""
import sys

import numpy as np
import veloxchem as vlx
from pyscf import df, gto

from veloxchem.veloxchemlib import (SimdThreeCenterElectronRepulsionDriver,
                                    SimdThreeCenterElectronRepulsionGradientDriver,
                                    SimdThreeCenterElectronRepulsionGradientRsDriver)

OMEGAS = (0.2, 0.33, 3.0)
THRESHOLD = 1.0e-14
FITTING = "def2-universal-jkfit"
TOLERANCE = 1.0e-9

CASES = [(("H", (0.0, 0.0, 0.0)), ("F", (0.0, 0.0, 0.917)), "def2-svp"),
         (("Li", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, 1.595)), "def2-svp"),
         (("C", (0.0, 0.0, 0.0)), ("O", (0.0, 0.0, 1.128)), "def2-tzvp"),
         (("O", (0.0, 0.0, 0.1173)), ("H", (0.0, 0.7572, -0.4692)),
          ("H", (0.0, -0.7572, -0.4692)), "def2-svp")]

fails = []


def as_xyz(geometry):
    return f"{len(geometry)}\n\n" + "\n".join(
        f"{a}  {x} {y} {z}" for a, (x, y, z) in geometry)


def squares(tensor):
    """The sum of the squares of every value, counting an atom pair both ways round.

    The pattern keeps the upper triangle of the **atom pairs**, and that is not the
    same as the upper triangle of the atom basis pairs. Water has a block whose two
    sides are the same atom basis, the hydrogens, and which holds three pairs: H1
    with H1, H2 with H2, and H1 with H2. The first two stand for themselves and the
    third stands for two. Doubling whole blocks instead gets a diatomic right, where
    every block is one kind or the other, and a polyatomic wrong -- which is how
    this was found, water failing by 4e-03 on the Coulomb column that pyscf and the
    unattenuated driver already agree on.

    The atom pairs are the fastest varying index of a block, so the values reshape
    to (everything else, pairs) and the weights apply to the columns.
    """
    total = 0.0

    for i in range(tensor.number_of_blocks()):
        block = tensor.block(i)
        npairs = block.number_of_pairs()

        values = np.asarray(tensor.block_to_numpy(i))

        assert values.size % npairs == 0, "the pairs are not the fastest index"

        weights = np.array([2.0 if a != b else 1.0
                            for a, b in zip(block.a_atoms(), block.b_atoms())])

        total += float(np.sum(values.reshape(-1, npairs) ** 2 * weights))

    return total


def elementwise(a, b):
    return max(float(np.abs(np.asarray(a.block_to_numpy(i))
                            - np.asarray(b.block_to_numpy(i))).max())
               for i in range(a.number_of_blocks()))


for case in CASES:
    *atoms, orbital = case
    geometry = list(atoms)
    label = "".join(a for a, _ in geometry)

    mol = vlx.Molecule.read_xyz_string(as_xyz(geometry))
    bas = vlx.MolecularBasis.read(mol, orbital.upper(), ostream=None)
    aux = vlx.MolecularBasis.read(mol, FITTING.upper(), ostream=None)

    pattern = SimdThreeCenterElectronRepulsionDriver().make_pattern(
        mol, bas, aux, THRESHOLD)

    plain = SimdThreeCenterElectronRepulsionGradientDriver().compute(
        pattern, mol, bas, aux)

    pm = gto.M(atom=[[a, xyz] for a, xyz in geometry], basis=orbital,
               unit="Angstrom")
    am = df.addons.make_auxmol(pm, FITTING)

    for omega in OMEGAS:
        coulomb, attenuated = SimdThreeCenterElectronRepulsionGradientRsDriver().compute(
            pattern, mol, bas, aux, omega)

        exact = elementwise(coulomb, plain)

        ip_c = df.incore.aux_e2(pm, am, intor="int3c2e_ip1", comp=3)
        with pm.with_range_coulomb(omega):
            ip_a = df.incore.aux_e2(pm, am, intor="int3c2e_ip1", comp=3)

        ours_c, ours_a = squares(coulomb), squares(attenuated)
        ref_c = 2.0 * float(np.sum(ip_c ** 2))
        ref_a = 2.0 * float(np.sum(ip_a ** 2))

        rel_c = abs(ours_c - ref_c) / abs(ref_c)
        rel_a = abs(ours_a - ref_a) / abs(ref_a)

        ok = (exact == 0.0 and rel_c < TOLERANCE and rel_a < TOLERANCE)
        if not ok:
            fails.append(f"{label} omega {omega}")

        print(f"  {'ok  ' if ok else 'FAIL'} {label:4s} {orbital:9s} omega {omega:<5g}  "
              f"coulomb half exact {exact:.1e}   squares vs pyscf: coulomb {rel_c:.2e}  "
              f"attenuated {rel_a:.2e}")

# NOTE: at a large omega the attenuated operator is 1 / r and its derivative is the
# Coulomb one. Both tensors are veloxchem's and their layouts are identical, so this
# is an element for element comparison and needs nothing of pyscf.
mol = vlx.Molecule.read_xyz_string(as_xyz(list(CASES[2][:-1])))
bas = vlx.MolecularBasis.read(mol, CASES[2][-1].upper(), ostream=None)
aux = vlx.MolecularBasis.read(mol, FITTING.upper(), ostream=None)
pattern = SimdThreeCenterElectronRepulsionDriver().make_pattern(mol, bas, aux, THRESHOLD)
coulomb, attenuated = SimdThreeCenterElectronRepulsionGradientRsDriver().compute(
    mol and pattern, mol, bas, aux, 1.0e6)

biggest = max(float(np.abs(np.asarray(coulomb.block_to_numpy(i))).max())
              for i in range(coulomb.number_of_blocks()))
apart = elementwise(attenuated, coulomb) / max(biggest, 1.0)

ok = apart < 1.0e-10
if not ok:
    fails.append("the large omega limit")
print(f"\n  {'ok  ' if ok else 'FAIL'} at a large omega the attenuated derivative is the "
      f"Coulomb one: {apart:.2e}")

print("\nall checks passed" if not fails else "\nFAILURES: " + ", ".join(fails))
sys.exit(1 if fails else 0)
