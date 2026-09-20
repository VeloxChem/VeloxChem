//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#ifndef SimdRIFockCommon_hpp
#define SimdRIFockCommon_hpp

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "PackedMatrix.hpp"
#include "SparseTensor.hpp"
#include "TripleSparsityPattern.hpp"

/// @brief How a driver of the resolution of the identity forms its matrices.
/// rimode::in_memory - what a sweep of the integrals leaves is formed once and held
/// rimode::direct - the integrals are formed again on every call
/// @note It is here and not with either driver because both of them have the choice
/// and the calculation makes it once for whichever one it is about to use.
enum class rimode
{
    automatic,
    in_memory,
    direct
};

/// @brief What the drivers of the resolution of the identity share.
///
/// @note The Coulomb only driver and the one which also forms the exchange do the
/// same things to the auxiliary basis: they measure what each of its atoms carries
/// so the ranks can be given equal shares, they name the functions of a set of
/// atoms, and they form and invert the metric. Those belong to neither of them.
///
/// @note What is **not** here is anything which reads a driver's own state. The
/// sweeps of the integrals and the fitting differ between the two -- one inverts the
/// metric outright and applies it to a vector, the other holds a factor of it and
/// applies that to a tensor -- so they stay with their drivers.
namespace simdri {

/// @brief Measures what each atom of the auxiliary side of a pattern carries.
/// @param pattern The sparsity pattern to measure.
/// @param natoms The number of atoms of the molecule.
/// @return The memory of the values of each atom, in bytes.
auto atom_shares(const CTripleSparsityPattern &pattern, const size_t natoms) -> std::vector<double>;

/// @brief Gets the dense indices of the auxiliary basis functions of given atoms.
/// @param aux_basis The auxiliary molecular basis.
/// @param atoms The atoms, as their indices in the molecule, or none of them for all
/// of them.
/// @return The indices, in ascending order.
auto aux_functions_of(const CMolecularBasis &aux_basis, const std::vector<int> &atoms) -> std::vector<size_t>;

/// @brief Inverts the metric of the auxiliary basis.
/// @param two_center The two-center integrals of the auxiliary basis.
/// @param metric_threshold The threshold of the linear dependence.
/// @param use_inverse_square_root Whether to invert the square root of the metric.
/// @param what What the metric is, named in the warning a fallback prints.
/// @return The inverted metric.
/// @brief The memory the values of a pattern of the three-center integrals take.
/// @param molecule The molecule.
/// @param basis The molecular basis.
/// @param aux_basis The auxiliary molecular basis.
/// @param threshold The screening threshold.
/// @param aux_atoms The atoms of the auxiliary side to count, or none of them for
/// the whole molecule.
/// @return The memory of one tensor of that pattern, in bytes.
/// @note One tensor. A driver which holds more than one of them -- a hybrid range
/// separated functional holds two, on the same pattern -- multiplies this itself, so
/// that what it answers is the memory it will really hold.
/// @note The memory answered is that of the atoms asked for, which is this rank's
/// when the auxiliary basis is divided over a communicator. Answering the whole
/// molecule's would put every rank on the direct way for a calculation each of them
/// holds a fitting share of.
auto pattern_memory(const CMolecule        &molecule,
                    const CMolecularBasis  &basis,
                    const CMolecularBasis  &aux_basis,
                    const double            threshold,
                    const std::vector<int> &aux_atoms) -> size_t;

/// @brief Measures what each atom of the auxiliary basis carries, for dividing them.
/// @return The memory of the values of each atom, in bytes.
auto aux_atom_weights(const CMolecule       &molecule,
                      const CMolecularBasis &basis,
                      const CMolecularBasis &aux_basis,
                      const double           threshold) -> std::vector<double>;

/// @brief The three-center integrals of one part of the auxiliary basis.
/// @param pattern The sparsity pattern of the part.
/// @param molecule The molecule.
/// @param basis The molecular basis on a and b sides.
/// @param aux_basis The auxiliary molecular basis.
/// @return The integrals, allocated and computed.
/// @note Both sweeps of the direct way need this and neither needs anything else of
/// the other: the first contracts what comes back against a density and the second
/// against the coefficients of the fitting. Which of the two a driver does, and with
/// what factor, stays with the driver.
auto integrals_of_part(const CTripleSparsityPattern &pattern,
                       const CMolecule              &molecule,
                       const CMolecularBasis        &basis,
                       const CMolecularBasis        &aux_basis) -> CSparseTensor;

/// @brief Divides the auxiliary basis into the parts a sweep takes one at a time.
/// @param pattern The sparsity pattern of the whole auxiliary basis.
/// @param min_parts The fewest parts to make, so that the ranks of a communicator
/// have one each at least.
/// @param budget The memory a part's integrals may take, in bytes.
/// @return The sparsity pattern of each part.
/// @note The parts are cut at whichever is the smaller of what the memory allows and
/// an equal division into the number asked for. Cutting finer costs nothing: the
/// parts are a division of the same atoms and their integrals are the same integrals
/// however they are grouped.
auto make_parts(const CMolecule              &molecule,
                const CMolecularBasis        &basis,
                const CMolecularBasis        &aux_basis,
                const double                  threshold,
                const CTripleSparsityPattern &pattern,
                const size_t                  min_parts,
                const size_t                  budget) -> std::vector<CTripleSparsityPattern>;

/// @brief Inverts the metric of the auxiliary basis outright.
/// @param two_center The two-center integrals of the auxiliary basis.
/// @param metric_threshold The threshold of the linear dependence.
/// @return The inverse, as a symmetric matrix.
/// @note This is what a Coulomb only fitting wants. It applies the metric to a
/// vector of one value per auxiliary function, twice a build, so the inverse itself
/// is the cheap thing to hold and one multiply is the cheap thing to do. A driver
/// which applies the metric to a **tensor** wants a factor of it instead, which is
/// what invert_metric answers, because a factor halves the work of that.
auto invert_metric_full(const CPackedMatrix &two_center, const double metric_threshold) -> CPackedMatrix;

auto invert_metric(const CPackedMatrix &two_center,
                   const double         metric_threshold,
                   const bool           use_inverse_square_root,
                   const std::string   &what) -> CPackedMatrix;

/// @brief Forms the metric of the auxiliary basis, and the way it is for.
/// @return The metric, and the way it is for.
auto form_metric(const CMolecule       &molecule,
                 const CMolecularBasis &aux_basis,
                 const double           metric_threshold,
                 const bool             use_inverse_square_root,
                 const rimode           mode,
                 double                *two_center_time,
                 double                *metric_time) -> std::pair<CPackedMatrix, rimode>;

}  // namespace simdri

#endif /* SimdRIFockCommon_hpp */
