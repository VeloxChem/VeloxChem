//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#ifndef SimdThreeCenterElectronRepulsionGeom100Func_hpp
#define SimdThreeCenterElectronRepulsionGeom100Func_hpp

#include <cstddef>
#include <string>

#include "BasisFunction.hpp"
#include "ErrorHandler.hpp"
#include "SimdMatrix.hpp"

namespace simdt3cerigrad {  // simdt3cerigrad namespace

/// @brief The number of components of the derivative a combination writes.
/// @note Three Cartesian directions for each of the two centers on bra side. The
/// center of the auxiliary function is not differentiated: the three derivatives
/// of an integral sum to zero, so the caller takes the third as the negative of
/// the other two.
inline constexpr size_t number_of_components = 6;

/// @brief Computes the derivative of the three-center electron repulsion
/// integrals with respect to the positions of the two atoms on bra side.
/// @param values The values of the derivative, one row per component and pair of
/// angular components, each of one value per atom pair.
/// @param nvalues The number of atom pairs.
/// @param natoms The number of atoms on the auxiliary side of the block.
/// @param bra The basis function of the first center on bra side.
/// @param ket The basis function of the second center on bra side.
/// @param aux The basis function of the auxiliary center.
/// @param coordinates The coordinates of the atom pairs.
/// @param c_coordinates The coordinates of the atoms on the auxiliary side.
/// @param buffer The scratch of the combination of basis functions.
/// @param threshold The screening threshold of the pattern.
/// @note The six components of a combination are consecutive, the whole run of a
/// component before the next, which is how the offsets of a block scale by six.
/// @note **The kernels are not written yet.** Every combination refuses, so that
/// a caller which reaches one is told rather than handed zeros: a gradient of
/// zeros is a gradient which looks converged everywhere.
inline auto
compute_electron_repulsion_geom_100(double               *values,
                                    const size_t          nvalues,
                                    const size_t          natoms,
                                    const CBasisFunction &bra,
                                    const CBasisFunction &ket,
                                    const CBasisFunction &aux,
                                    const CSimdMatrix    &coordinates,
                                    const CSimdMatrix    &c_coordinates,
                                    CSimdMatrix          &buffer,
                                    const double          threshold) -> void
{
    errors::assertMsgCritical(
        false,
        std::string("SimdThreeCenterElectronRepulsionGeom100Func: No kernel for the combination of angular momenta ") +
            std::to_string(bra.get_angular_momentum()) + std::string(", ") + std::to_string(ket.get_angular_momentum()) +
            std::string(" and ") + std::to_string(aux.get_angular_momentum()) +
            std::string(". The derivative kernels are not written yet"));

    (void)values;
    (void)nvalues;
    (void)natoms;
    (void)coordinates;
    (void)c_coordinates;
    (void)buffer;
    (void)threshold;
}

/// @brief Gets the number of rows of the buffer a combination needs.
/// @note A placeholder until the kernels are generated, which will bring their
/// own table. The derivative raises the momentum of the differentiated center by
/// one and writes six components, so it is bounded by six times the rows of the
/// integral one shell higher; it is never read until a kernel exists to read it.
inline auto
number_of_buffer_rows(const int a_angular_momentum, const int b_angular_momentum, const int c_angular_momentum) -> size_t
{
    errors::assertMsgCritical((a_angular_momentum >= 0) && (a_angular_momentum < 7) && (b_angular_momentum >= 0) &&
                                  (b_angular_momentum < 7) && (c_angular_momentum >= 0) && (c_angular_momentum < 9),
                              std::string("SimdThreeCenterElectronRepulsionGeom100Func: Angular momentum is out of range"));

    const auto la = static_cast<size_t>(a_angular_momentum) + 1;

    const auto lb = static_cast<size_t>(b_angular_momentum) + 1;

    const auto lc = static_cast<size_t>(c_angular_momentum);

    return number_of_components * (la + 1) * (la + 2) * (lb + 1) * (lb + 2) * (lc + 1) * (lc + 2) / 8 + 64;
}

}  // namespace simdt3cerigrad

#endif /* SimdThreeCenterElectronRepulsionGeom100Func_hpp */
