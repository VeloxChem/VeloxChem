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
