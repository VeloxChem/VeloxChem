//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#include "SimdRIFockCommon.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>

#include "DenseIndexFunc.hpp"
#include "ErrorHandler.hpp"
#include "Eigen/Dense"
#include "MathLibrary.hpp"
#include "PackedLinearAlgebra.hpp"
#include "SimdT3CDistributor.hpp"
#include "SimdThreeCenterElectronRepulsionDriver.hpp"
#include "SimdTwoCenterElectronRepulsionDriver.hpp"
#include "StringFormat.hpp"

namespace simdri {

namespace {

using prof_clock = std::chrono::steady_clock;

/// @brief The seconds since a mark.
auto
prof_since(const prof_clock::time_point &mark) -> double
{
    return std::chrono::duration<double>(prof_clock::now() - mark).count();
}

}  // namespace

/// @brief Measures what each atom of the auxiliary side of a pattern carries.
/// @param pattern The sparsity pattern to measure.
/// @param natoms The number of atoms of the molecule.
/// @return The memory of the values of each atom, in bytes.
/// @note A block holds as many values for one of its atoms on the auxiliary side as
/// for any other, so its memory divides evenly over them and the memory of an atom is
/// the sum of the shares of the blocks which carry it.
auto
atom_shares(const CTripleSparsityPattern &pattern, const size_t natoms) -> std::vector<double>
{
    std::vector<double> shares(natoms, 0.0);

    for (const auto &block : pattern.blocks())
    {
        const auto &c_atoms = block.c_atoms();

        if (c_atoms.empty()) continue;

        const auto share =
            static_cast<double>(block.number_of_elements() * sizeof(double)) / static_cast<double>(c_atoms.size());

        for (const auto atom : c_atoms) shares[static_cast<size_t>(atom)] += share;
    }

    return shares;
}

/// @brief Gets the dense indices of the auxiliary basis functions of given atoms.
/// @param aux_basis The auxiliary molecular basis.
/// @param atoms The atoms, as their indices in the molecule, or none of them for all
/// of them.
/// @return The indices, in ascending order.
/// @note The dense index runs over the angular momenta of the whole molecule before
/// it runs over the atoms, so the functions of one atom are scattered through it and
/// the functions of a set of atoms are a set rather than a range.
auto
aux_functions_of(const CMolecularBasis &aux_basis, const std::vector<int> &atoms) -> std::vector<size_t>
{
    const auto set_indices = aux_basis.basis_sets_indices();

    const auto natoms = set_indices.size();

    const auto indices = denseidx::index_functions(aux_basis);

    const auto starts = denseidx::make_dense_starts(aux_basis);

    const auto strides = denseidx::make_dense_strides(aux_basis);

    const auto nmoms = static_cast<size_t>(aux_basis.max_angular_momentum() + 1);

    std::vector<int> all_atoms;

    if (atoms.empty())
    {
        all_atoms.reserve(natoms);

        for (size_t atom = 0; atom < natoms; atom++) all_atoms.push_back(static_cast<int>(atom));
    }

    std::vector<size_t> functions;

    for (const auto atom : (atoms.empty() ? all_atoms : atoms))
    {
        const auto index = static_cast<size_t>(atom);

        for (const auto [lc, kc] : indices[static_cast<size_t>(set_indices[index])])
        {
            const auto lval = static_cast<size_t>(lc);

            for (size_t mc = 0; mc < static_cast<size_t>(2 * lc + 1); mc++)
            {
                functions.push_back(starts[index * nmoms + lval] + kc + mc * strides[lval]);
            }
        }
    }

    std::sort(functions.begin(), functions.end());

    return functions;
}

/// @brief Inverts a metric by the route asked for, falling back where there is no
/// Cholesky factor to be had.
/// @param two_center The metric of the fitting basis.
/// @param metric_threshold The threshold below which a direction is dropped.
/// @param use_inverse_square_root Whether to invert the square root of the metric.
/// @param what What the metric is, named in the warning a fallback prints, so that a
/// caller inverting two of them says which of the two it was.
/// @return The inverted metric.
auto
invert_metric(const CPackedMatrix &two_center,
              const double         metric_threshold,
              const bool           use_inverse_square_root,
              const std::string   &what) -> CPackedMatrix
{
    if (use_inverse_square_root) return packlin::inverse_square_root(two_center, metric_threshold);

    try
    {
        return packlin::cholesky_inverse(two_center);
    }
    catch (const std::runtime_error &)
    {
        errors::msg(std::string("RIJKFockDriver: The ") + what +
                        std::string(" has no Cholesky factor, so its square root is inverted instead. This is a "
                                    "nearly linearly dependent fitting basis."),
                    "Warning");

        return packlin::inverse_square_root(two_center, metric_threshold);
    }
}

/// @brief Forms the metric a way of building asks for, and the way it is for.
/// @param molecule The molecule to compute the metric of.
/// @param aux_basis The auxiliary molecular basis.
/// @param metric_threshold The threshold below which a direction is dropped.
/// @param use_inverse_square_root Whether to invert the square root of the metric.
/// @param mode The way of building the metric is for.
/// @param two_center_time Where to add the time of the two-center integrals, if
/// anywhere.
/// @param metric_time Where to add the time of the inversion, if anywhere.
/// @return The metric, and the way it is for.
auto
form_metric(const CMolecule       &molecule,
            const CMolecularBasis &aux_basis,
            const double           metric_threshold,
            const bool             use_inverse_square_root,
            const rimode           mode,
            double                *two_center_time,
            double                *metric_time) -> std::pair<CPackedMatrix, rimode>
{
    const auto mark_two_center = prof_clock::now();

    const auto two_center = CSimdTwoCenterElectronRepulsionDriver().compute(molecule, aux_basis);

    if (two_center_time) *two_center_time += prof_since(mark_two_center);

    const auto mark_metric = prof_clock::now();

    // NOTE: both forms of the metric close the resolution of the identity, and the
    // Cholesky factor costs an order of magnitude less, so it is tried first. A
    // fitting basis which is close to linearly dependent has none, and the square
    // root is inverted in its place, dropping the directions which carry nothing.

    // NOTE: the direct way solves with the factor rather than multiplying by its
    // inverse, so it is the factor which is kept. The inverted square root is
    // taken for a metric which has no factor, in either way, as the B vectors
    // formed with it close the same sum.

    if (mode == rimode::direct)
    {
        if (!use_inverse_square_root)
        {
            try
            {
                auto factor = packlin::cholesky_factor(two_center);

                if (metric_time) *metric_time += prof_since(mark_metric);

                return {std::move(factor), rimode::direct};
            }
            catch (const std::runtime_error &)
            {
                errors::msg(std::string("RIJKFockDriver: The metric of the fitting basis has no Cholesky factor, so "
                                        "its square root is inverted instead and multiplied by. This is a nearly "
                                        "linearly dependent fitting basis."),
                            "Warning");
            }
        }

        // NOTE: the direct way solves the Cholesky factor against the half
        // transformed integrals where it has one, and multiplies by the inverted
        // square root where that is what was asked for, or where there is no factor
        // to be had. The two close the same sum: solving the factor gives B with
        // B^T B equal to A^T V^-1 A, and so does multiplying by the root, the root
        // being its own transpose. The root costs twice the arithmetic and is a
        // product rather than a substitution, which is the trade. Falling back to
        // holding the B vectors instead, which is what this did, asks for the
        // memory the direct way was chosen for want of.

        auto root = packlin::inverse_square_root(two_center, metric_threshold);

        if (metric_time) *metric_time += prof_since(mark_metric);

        return {std::move(root), rimode::direct};
    }

    auto metric = invert_metric(two_center, metric_threshold, use_inverse_square_root,
                                "metric of the fitting basis");

    if (metric_time) *metric_time += prof_since(mark_metric);

    return {std::move(metric), rimode::in_memory};
}

auto
pattern_memory(const CMolecule        &molecule,
               const CMolecularBasis  &basis,
               const CMolecularBasis  &aux_basis,
               const double            threshold,
               const std::vector<int> &aux_atoms) -> size_t
{
    // NOTE: the pattern of what a sweep leaves is the pattern of the three-center
    // integrals. The metric is dense and the transformation of the auxiliary side
    // keeps every atom which survives, so neither of them makes it sparser.

    const CSimdThreeCenterElectronRepulsionDriver eri_drv;

    const auto pattern = aux_atoms.empty() ? eri_drv.make_pattern(molecule, basis, aux_basis, threshold)
                                           : eri_drv.make_pattern(molecule, basis, aux_basis, threshold, aux_atoms);

    size_t nvalues = 0;

    for (size_t i = 0; i < static_cast<size_t>(pattern.number_of_blocks()); i++)
    {
        nvalues += pattern.block(i).number_of_elements();
    }

    return nvalues * sizeof(double);
}

auto
aux_atom_weights(const CMolecule       &molecule,
                 const CMolecularBasis &basis,
                 const CMolecularBasis &aux_basis,
                 const double           threshold) -> std::vector<double>
{
    const auto pattern = CSimdThreeCenterElectronRepulsionDriver().make_pattern(molecule, basis, aux_basis, threshold);

    return atom_shares(pattern, static_cast<size_t>(molecule.number_of_atoms()));
}

auto
integrals_of_part(const CTripleSparsityPattern &pattern,
                  const CMolecule              &molecule,
                  const CMolecularBasis        &basis,
                  const CMolecularBasis        &aux_basis) -> CSparseTensor
{
    auto integrals = CSparseTensor(pattern);

    integrals.allocate();

    auto distributor = CSimdT3CDistributor<CSparseTensor>(&integrals);

    CSimdThreeCenterElectronRepulsionDriver().compute(pattern, molecule, basis, aux_basis, distributor);

    return integrals;
}

auto
make_parts(const CMolecule              &molecule,
           const CMolecularBasis        &basis,
           const CMolecularBasis        &aux_basis,
           const double                  threshold,
           const CTripleSparsityPattern &pattern,
           const size_t                  min_parts,
           const size_t                  budget) -> std::vector<CTripleSparsityPattern>
{
    const auto natoms = static_cast<size_t>(molecule.number_of_atoms());

    const auto shares = simdri::atom_shares(pattern, natoms);

    // NOTE: the parts are cut at whichever is the smaller of what the memory allows
    // and an equal division into the number asked for. A machine with memory to spare
    // gives one part, and the caller which divides the Coulomb pass over the ranks of
    // a communicator then has one part for all of them: one rank sweeps the integrals
    // a second time and the others wait. Cutting finer costs nothing, as the parts are
    // a division of the same atoms and their integrals are the same integrals however
    // they are grouped.

    const auto total = std::accumulate(shares.begin(), shares.end(), 0.0);

    const auto by_memory = static_cast<double>(budget / 2);

    const auto by_parts = total / static_cast<double>(std::max(min_parts, size_t{1}));

    const auto cut = std::min(by_memory, by_parts);

    // NOTE: the atoms are gathered in the order they are given until the integrals
    // of a part reach the budget. An atom whose own integrals are above it is a part
    // of its own, as there is nothing smaller to divide.

    const CSimdThreeCenterElectronRepulsionDriver eri_drv;

    std::vector<CTripleSparsityPattern> parts;

    std::vector<int> atoms;

    double memory = 0.0;

    for (size_t atom = 0; atom < natoms; atom++)
    {
        const auto share = shares[atom];

        if (share <= 0.0) continue;

        if ((!atoms.empty()) && ((memory + share) > cut))
        {
            parts.push_back(eri_drv.make_pattern(molecule, basis, aux_basis, threshold, atoms));

            atoms.clear();

            memory = 0.0;
        }

        atoms.push_back(static_cast<int>(atom));

        memory += share;
    }

    if (!atoms.empty()) parts.push_back(eri_drv.make_pattern(molecule, basis, aux_basis, threshold, atoms));

    return parts;
}


auto
invert_metric_full(const CPackedMatrix &two_center, const double metric_threshold) -> CPackedMatrix
{
    // NOTE: the inverse itself and not a factor of it. A Coulomb only fitting applies
    // the metric to a vector of one value per auxiliary function, twice a build, so
    // there is nothing to be gained from a factor and one multiply is the whole of it.
    // The driver which applies the metric to a tensor wants a factor, and asks for one.

    // NOTE: the directions the fitting basis does not really span are dropped rather
    // than inverted, on the same threshold and by the same construction the inverted
    // square root uses. A plain factorization would invert them: a metric of
    // def2-universal-jkfit on thirty-two waters has thirty-three eigenvalues under
    // 1e-8 of its largest, and dividing by those multiplies numerical noise into
    // every fitting coefficient the driver forms.

    return packlin::pseudo_inverse(two_center, metric_threshold);
}

}  // namespace simdri
