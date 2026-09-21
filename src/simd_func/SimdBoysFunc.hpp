//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#ifndef SimdBoysFunc_hpp
#define SimdBoysFunc_hpp

#include <cstddef>
#include <initializer_list>

#include "BasisFunction.hpp"
#include "SimdMatrix.hpp"

namespace simdfunc {  // simdfunc namespace

// NOTE: this header holds two layers. compute_boys_values below is the Boys
// function itself and knows nothing of any operator: it reads the arguments from a
// row of a buffer and writes the values to the rows which follow. The two forms
// under it are the wrapper the two-center electron repulsion kernels call, which
// forms the argument from the geometry and the exponents of a pair of primitives
// and scales the values by the prefactor of the integral. An operator whose
// argument or prefactor differs wants its own wrapper and the same core.

/// @brief Computes the values of Boys function of every order up to the requested
/// one, from the arguments already in the buffer.
/// @param buffer The buffer holding the arguments in the row at target.
/// @param target The row of the buffer holding the arguments. The values follow it,
/// so the routine reads one row and writes order + 1 of them.
/// @param order The highest order of Boys function to compute.
/// @param ncols The number of arguments to compute the values of.
/// @note The rows from target + 1 to target + order + 1 hold the values of order
/// zero to order. The arguments are not modified.
/// @note The values are computed from the Taylor expansion on a grid of arguments
/// below the end of the grid and from the asymptotic expansion above it, with the
/// remaining orders following from the downward and the upward recursion
/// respectively. The relative accuracy is better than 1.0e-14 for every order up to
/// the highest supported one and for every argument.
auto compute_boys_values(CSimdMatrix &buffer, const size_t target, const size_t order, const size_t ncols) -> void;

/// @brief Computes the values of Boys function of the requested orders alone, from
/// the arguments already in the buffer.
/// @param buffer The buffer holding the arguments in the row at target.
/// @param target The row of the buffer holding the arguments. The values follow it,
/// one row per requested order and in the order they are requested.
/// @param orders The orders of Boys function to compute.
/// @param ncols The number of arguments to compute the values of.
/// @note A recursion carries the values from the order the expansion gives to the
/// order which is wanted, so the orders below the highest requested one are computed
/// whether they are asked for or not. What this form saves is the rows of the
/// buffer, not the arithmetic: a caller which needs two orders out of twelve holds
/// three rows rather than fourteen.
auto compute_boys_values(CSimdMatrix &buffer, const size_t target, const std::initializer_list<size_t> orders, const size_t ncols)
    -> void;

/// @brief Computes the values of Boys function of every order up to the requested
/// one for one pair of primitives of a two-center electron repulsion integral.
/// @param buffer The buffer of the combination of basis functions.
/// @param coordinates The coordinates of the atom pairs, whose row nine holds the
/// squared distance of the atom pair.
/// @param target The row of the buffer to write the argument to, with the values
/// following it.
/// @param order The highest order of Boys function to compute.
/// @param ncols The number of atom pairs the pair of primitives reaches.
/// @param fj The prefactor of the integral, which every value is scaled by.
/// @param mu The factor the squared distance is scaled by to give the argument.
auto compute_full_boys_function(CSimdMatrix       &buffer,
                                const CSimdMatrix &coordinates,
                                const size_t       target,
                                const size_t       order,
                                const size_t       ncols,
                                const double       fj,
                                const double       mu) -> void;

/// @brief Computes the values of Boys function of the requested orders alone for one
/// pair of primitives of a two-center electron repulsion integral.
/// @param buffer The buffer of the combination of basis functions.
/// @param coordinates The coordinates of the atom pairs, whose row nine holds the
/// squared distance of the atom pair.
/// @param target The row of the buffer to write the argument to, with the values
/// following it.
/// @param orders The orders of Boys function to compute.
/// @param ncols The number of atom pairs the pair of primitives reaches.
/// @param fj The prefactor of the integral, which every value is scaled by.
/// @param mu The factor the squared distance is scaled by to give the argument.
auto compute_boys_function(CSimdMatrix                        &buffer,
                           const CSimdMatrix                  &coordinates,
                           const size_t                        target,
                           const std::initializer_list<size_t> orders,
                           const size_t                        ncols,
                           const double                        fj,
                           const double                        mu) -> void;

/// @brief Computes the values of the attenuated Boys function of every order up to
/// the requested one for one pair of primitives of a two-center electron repulsion
/// integral of the operator erf(omega r) / r.
/// @param buffer The buffer of the combination of basis functions.
/// @param coordinates The coordinates of the atom pairs, whose row nine holds the
/// squared distance of the atom pair.
/// @param target The row of the buffer to write the argument to, with the values
/// following it.
/// @param order The highest order to compute.
/// @param ncols The number of atom pairs the pair of primitives reaches.
/// @param fj The prefactor of the integral.
/// @param mu The factor the squared distance is scaled by for the plain operator.
/// @param omega The range separation parameter.
/// @note The attenuated function is the plain one with its argument scaled by
/// theta squared and its value of order m scaled by theta to the 2m + 1, where
/// theta squared is omega squared over omega squared plus mu. It is therefore the
/// same expansion and the same accuracy; what the attenuation changes is where the
/// expansion is evaluated and what each order is weighted by.
/// @note Sending omega to infinity returns the plain function exactly, and sending
/// it to zero returns zeros, which is what erf(0 r) / r is. Both are checked.
auto compute_full_erf_boys_function(CSimdMatrix       &buffer,
                                    const CSimdMatrix &coordinates,
                                    const size_t       target,
                                    const size_t       order,
                                    const size_t       ncols,
                                    const double       fj,
                                    const double       mu,
                                    const double       omega) -> void;

/// @brief Computes the values of the attenuated Boys function of the requested
/// orders alone for one pair of primitives of a two-center electron repulsion
/// integral of the operator erf(omega r) / r.
/// @param buffer The buffer of the combination of basis functions.
/// @param coordinates The coordinates of the atom pairs, whose row nine holds the
/// squared distance of the atom pair.
/// @param target The row of the buffer to write the argument to, with the values
/// following it.
/// @param orders The orders to compute.
/// @param ncols The number of atom pairs the pair of primitives reaches.
/// @param fj The prefactor of the integral.
/// @param mu The factor the squared distance is scaled by for the plain operator.
/// @param omega The range separation parameter.
/// @note The rows are packed consecutively while the orders they stand for need not
/// be, so each row is weighted by the order it stands for and not by where it sits.
auto compute_erf_boys_function(CSimdMatrix                        &buffer,
                               const CSimdMatrix                  &coordinates,
                               const size_t                        target,
                               const std::initializer_list<size_t> orders,
                               const size_t                        ncols,
                               const double                        fj,
                               const double                        mu,
                               const double                        omega) -> void;

/// @brief Computes the values of Boys function of every order up to the requested
/// one for one triple of primitives of a three-center electron repulsion integral.
/// @param buffer The buffer of the combination of basis functions, holding the
/// displacement of the Gaussian product center from the atom on the ket side in the
/// three rows at pc.
/// @param coordinates The coordinates of the atom pairs, whose row nine holds the
/// squared distance of the atom pair.
/// @param target The row of the buffer to write the argument to, with the values
/// following it.
/// @param pc The first of the three rows holding that displacement, as
/// simdfunc::compute_pc wrote them.
/// @param order The highest order of Boys function to compute.
/// @param ncols The number of atom pairs the triple of primitives reaches.
/// @param fj The prefactor of the integral, which every value is scaled by.
/// @param mu The factor the squared distance of the atom pair is scaled by, whose
/// exponential every value is scaled by as well.
/// @param fq The factor the squared displacement is scaled by to give the argument.
/// @note The scaling is not one number here, as it is for the two-center form. The
/// pair of primitives contributes exp(-mu R_AB^2), which varies with the atom pair,
/// so the values are scaled column by column and not by fj alone.
/// @brief Computes the values of Boys function of the requested orders alone for one
/// triple of primitives of a three-center electron repulsion integral.
/// @param buffer The buffer of the combination of basis functions, holding the
/// displacement of the Gaussian product center from the atom on the ket side in the
/// three rows at pc.
/// @param coordinates The coordinates of the atom pairs, whose row nine holds the
/// squared distance of the atom pair.
/// @param target The row of the buffer to write the argument to, with the values
/// following it, one row per requested order and in the order they are requested.
/// @param pc The first of the three rows holding that displacement, as
/// simdfunc::compute_pc wrote them.
/// @param orders The orders of Boys function to compute.
/// @param ncols The number of atom pairs the triple of primitives reaches.
/// @param fj The prefactor of the integral, which every value is scaled by.
/// @param mu The factor the squared distance of the atom pair is scaled by, whose
/// exponential every value is scaled by as well.
/// @param fq The factor the squared displacement is scaled by to give the argument.
/// @note What this form saves against the full one is the rows of the buffer and not
/// the arithmetic, as the orders below the highest requested one are computed whether
/// they are asked for or not.
auto compute_t3c_boys_function(CSimdMatrix                        &buffer,
                               const CSimdMatrix                  &coordinates,
                               const size_t                        target,
                               const size_t                        pc,
                               const std::initializer_list<size_t> orders,
                               const size_t                        ncols,
                               const double                        fj,
                               const size_t                        pair,
                               const double                        fq) -> void;

auto compute_full_t3c_boys_function(CSimdMatrix       &buffer,
                                    const CSimdMatrix &coordinates,
                                    const size_t       target,
                                    const size_t       pc,
                                    const size_t       order,
                                    const size_t       ncols,
                                    const double       fj,
                                    const size_t       pair,
                                    const double       fq) -> void;

/// @brief Computes the values of the attenuated Boys function of every order up to
/// the requested one for one triple of primitives of a three-center electron
/// repulsion integral of the operator erf(omega r) / r.
/// @param fq The reduced exponent of the pair and the auxiliary primitive, which the
/// argument is formed from and which the attenuation is taken against.
/// @param omega The range separation parameter.
/// @note The attenuation is the same as in the two-center case: the argument is
/// scaled by theta squared and the value of order m by theta to the 2m + 1, with
/// theta squared being omega squared over omega squared plus fq. What differs is
/// only that the argument is the scaled displacement of the product center and that
/// the values carry the exponential of the pair as well.
auto compute_full_t3c_erf_boys_function(CSimdMatrix       &buffer,
                                        const CSimdMatrix &coordinates,
                                        const size_t       target,
                                        const size_t       pc,
                                        const size_t       order,
                                        const size_t       ncols,
                                        const double       fj,
                                        const size_t       pair,
                                        const double       fq,
                                        const double       omega) -> void;

/// @brief Computes the values of the attenuated Boys function of the requested
/// orders alone for one triple of primitives of a three-center electron repulsion
/// integral of the operator erf(omega r) / r.
/// @note The rows are packed consecutively while the orders they stand for need not
/// be, so each row is weighted by the order it stands for and not by where it sits.
auto compute_t3c_erf_boys_function(CSimdMatrix                        &buffer,
                                   const CSimdMatrix                  &coordinates,
                                   const size_t                        target,
                                   const size_t                        pc,
                                   const std::initializer_list<size_t> orders,
                                   const size_t                        ncols,
                                   const double                        fj,
                                   const size_t                        pair,
                                   const double                        fq,
                                   const double                        omega) -> void;

/// @brief Computes the values of Boys function of the requested orders alone for one
/// pair of primitives and one point charge of a nuclear attraction integral.
/// @param buffer The buffer of the combination of basis functions, holding the
/// displacement of the Gaussian product center from the charge in the three rows at
/// pc.
/// @param coordinates The coordinates of the atom pairs, whose row nine holds the
/// squared distance of the atom pair.
/// @param target The row of the buffer to write the argument to, with the values
/// following it, one row per requested order and in the order they are requested.
/// @param pc The first of the three rows holding that displacement, as
/// simdfunc::compute_pc wrote them.
/// @param orders The orders of Boys function to compute.
/// @param ncols The number of atom pairs the pair of primitives reaches.
/// @param fz The prefactor of the integral, 2 pi / p times the normalisations and the
/// charge, which every value is scaled by.
/// @param mu The factor the squared distance of the atom pair is scaled by, whose
/// exponential every value is scaled by as well.
/// @param p The exponent sum of the pair, which the squared displacement is scaled by
/// to give the argument.
/// @note The argument is p PC^2 where the electron repulsion's is mu R_AB^2, which is
/// the whole of the difference: the operator is anchored at the charge rather than
/// carried between two distributions. The scaling is the same shape as the
/// three-center form's, the pair having been collapsed onto one center either way.
/// @note What this form saves against the full one is the rows of the buffer and not
/// the arithmetic, as the orders below the highest requested one are computed whether
/// they are asked for or not.
auto compute_npot_boys_function(CSimdMatrix                        &buffer,
                                const CSimdMatrix                  &coordinates,
                                const size_t                        target,
                                const size_t                        pc,
                                const std::initializer_list<size_t> orders,
                                const size_t                        ncols,
                                const double                        fz,
                                const size_t                        pair_exp,
                                const double                        p) -> void;

auto compute_full_npot_boys_function(CSimdMatrix       &buffer,
                                     const CSimdMatrix &coordinates,
                                     const size_t       target,
                                     const size_t       pc,
                                     const size_t       order,
                                     const size_t       ncols,
                                     const double       fz,
                                     const size_t       pair_exp,
                                     const double       p) -> void;

/// @brief Writes the exponential the pair of primitives contributes,
/// exp(-mu AB^2), into one row.
/// @param buffer The buffer of the combination of basis functions.
/// @param coordinates The coordinates of the atom pairs, whose row nine holds the
/// squared distance of the atom pair.
/// @param target The row to write it to.
/// @param ncols The number of atom pairs the pair of primitives reaches.
/// @param mu The factor the squared distance of the atom pair is scaled by.
/// @note It depends on the pair of primitives and on the atom pair and on neither
/// the order of the Boys function nor the charge, so an anchored operator forms it
/// once above its loop over the charges and every call inside that loop reads it.
/// Evaluating it where it was read instead cost 17 to 26 per cent of a nuclear
/// attraction call, rising with the molecule, the loop over the charges being the
/// molecule.
/// @brief Forms the exponential every pair of primitives of the bra contributes,
/// exp(-mu AB^2), for all of them at once.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param coordinates The coordinates of the atom pairs, whose row nine holds the
/// squared distance of the atom pair.
/// @param ncols The number of atom pairs the widest pair of primitives reaches,
/// which is the stride of the scratch.
/// @note For a three-center form, whose loop over the atoms on the ket side stands
/// outside the pair of primitives. A row of the buffer would be refilled for every
/// one of those atoms; this is filled once for the call, and the pair is named
/// rather than the row. It costs nprim_a * nprim_b * ncols of scratch a thread,
/// which is what buys the whole repetition over the atoms.
/// @note The scratch is held per thread and grown rather than allocated per call.
/// **compute_t3c_boys_function and its full form read it, so it has to be filled
/// before them**; a kernel that does not is reading the last call's pairs.
auto compute_pair_exponents(const CBasisFunction &bra,
                            const CBasisFunction &ket,
                            const CSimdMatrix    &coordinates,
                            const size_t          ncols) -> void;

auto compute_pair_exponent(CSimdMatrix       &buffer,
                           const CSimdMatrix &coordinates,
                           const size_t       target,
                           const size_t       ncols,
                           const double       mu) -> void;

}  // namespace simdfunc

#endif /* SimdBoysFunc_hpp */
