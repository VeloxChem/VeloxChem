//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#include "SimdRIJKGradientDriver.hpp"

#include <numeric>
#include <set>
#include <string>

#include "ErrorHandler.hpp"

auto
CSimdRIJKGradientDriver::get_threshold() const -> double
{
    return _threshold;
}

auto
CSimdRIJKGradientDriver::get_block_size() const -> size_t
{
    return _block_size;
}

auto
CSimdRIJKGradientDriver::compute(const CMolecule        &molecule,
                                 const CMolecularBasis  &basis,
                                 const CMolecularBasis  &aux_basis,
                                 const CSparseTensor    &bq_vectors,
                                 const CPackedMatrix    &metric,
                                 const CPackedMatrix    &density,
                                 const CPackedMatrix    &coefficients,
                                 const std::vector<int> &atoms,
                                 const std::vector<int> &aux_atoms) const -> CDenseMatrix
{
    const auto natoms = molecule.number_of_atoms();

    // NOTE: an atom named twice would be accumulated twice, and an atom which is
    // not in the molecule has no row to accumulate into. Both are caught here
    // rather than by the arithmetic, which would give a gradient that is merely
    // wrong.
    for (const auto iatom : atoms)
    {
        errors::assertMsgCritical(
            (iatom >= 0) && (static_cast<size_t>(iatom) < natoms),
            std::string("SimdRIJKGradientDriver.compute: Atom is not an atom of the molecule"));
    }

    errors::assertMsgCritical(
        std::set<int>(atoms.begin(), atoms.end()).size() == atoms.size(),
        std::string("SimdRIJKGradientDriver.compute: An atom is named more than once"));

    errors::assertMsgCritical(
        density.get_type() == mat_t::symmetric,
        std::string("SimdRIJKGradientDriver.compute: The density matrix is expected to be symmetric"));

    auto gradient = CDenseMatrix(natoms, 3);

    gradient.zero();

    if (atoms.empty()) return gradient;

    // TODO: the derivative integrals and their contraction. The kernels which
    // form them are not written yet -- the plain tree has the derivative of the
    // three-center integral with respect to the auxiliary center alone, and none
    // of them are vectorised -- so this driver is a skeleton whose arguments are
    // fixed and whose body is not.
    //
    // What it will do, once they are: form the fitting coefficients from the B
    // vectors and the density, then for each atom of `atoms` accumulate the
    // derivative of the three-center integrals against them for the Coulomb part
    // and against the W matrices of the occupied orbitals for the exchange part,
    // and the derivative of the metric against the coefficients on both sides.
    //
    // The unused arguments are named and not commented out so that the shape of
    // the call does not change when the body arrives.
    (void)basis;
    (void)aux_basis;
    (void)bq_vectors;
    (void)metric;
    (void)coefficients;
    (void)aux_atoms;

    return gradient;
}

auto
CSimdRIJKGradientDriver::compute(const CMolecule       &molecule,
                                 const CMolecularBasis &basis,
                                 const CMolecularBasis &aux_basis,
                                 const CSparseTensor   &bq_vectors,
                                 const CPackedMatrix   &metric,
                                 const CPackedMatrix   &density,
                                 const CPackedMatrix   &coefficients) const -> CDenseMatrix
{
    // NOTE: every atom of the molecule, which is what a gradient usually means.
    // The form which takes a list is for a caller holding a share of them.
    auto atoms = std::vector<int>(molecule.number_of_atoms());

    std::iota(atoms.begin(), atoms.end(), 0);

    return compute(molecule, basis, aux_basis, bq_vectors, metric, density, coefficients, atoms, {});
}
