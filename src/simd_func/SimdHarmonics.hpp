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



#ifndef SimdHarmonics_hpp
#define SimdHarmonics_hpp

#include "SimdMatrix.hpp"

namespace simdfunc {  // simdfunc namespace

/// @brief Creates the real solid harmonics of angular momentum one of the vectors
/// between the atoms of the atom pairs.
/// @param coordinates The coordinates of the atom pairs, as seven rows holding
/// the coordinates of the atoms on bra side in rows zero to two and of the atoms
/// on ket side in rows three to five.
/// @return The matrix of three rows and as many columns as the coordinates, with
/// the row of index m + 1 holding the solid harmonic of order m.
/// @note The vector between the atoms is taken from the atom on ket side to the
/// atom on bra side, as the squared distance of the coordinates is, so the sign
/// of the solid harmonics of odd angular momenta follows that convention.
/// @note The solid harmonics are ordered and normalized as the spherical
/// components of a basis function of the same angular momentum, i.e. the row of
/// index zero holds y, the row of index one holds z and the row of index two
/// holds x. The padding of the rows is left undefined, as the atom pairs beyond
/// the number of columns are not part of the sparsity pattern.
auto make_p_solid_harmonics(const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum 2 of the vectors
/// between the atoms of the atom pairs.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 5 rows and as many columns as the coordinates, with
/// the row of index m + 2 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570. The rows of all orders are formed in a
/// single loop, so that every value entering the recursion is loaded once.
auto make_d_solid_harmonics(const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum 3 of the vectors
/// between the atoms of the atom pairs.
/// @param d_harmonics The solid harmonics of angular momentum 2.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 7 rows and as many columns as the coordinates, with
/// the row of index m + 3 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570. The rows of all orders are formed in a
/// single loop, so that every value entering the recursion is loaded once.
auto make_f_solid_harmonics(const CSimdMatrix &d_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum 4 of the vectors
/// between the atoms of the atom pairs.
/// @param f_harmonics The solid harmonics of angular momentum 3.
/// @param d_harmonics The solid harmonics of angular momentum 2.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 9 rows and as many columns as the coordinates, with
/// the row of index m + 4 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570. The rows of all orders are formed in a
/// single loop, so that every value entering the recursion is loaded once.
auto make_g_solid_harmonics(const CSimdMatrix &f_harmonics, const CSimdMatrix &d_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum 5 of the vectors
/// between the atoms of the atom pairs.
/// @param g_harmonics The solid harmonics of angular momentum 4.
/// @param f_harmonics The solid harmonics of angular momentum 3.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 11 rows and as many columns as the coordinates, with
/// the row of index m + 5 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570. The rows of all orders are formed in a
/// single loop, so that every value entering the recursion is loaded once.
auto make_h_solid_harmonics(const CSimdMatrix &g_harmonics, const CSimdMatrix &f_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum 6 of the vectors
/// between the atoms of the atom pairs.
/// @param h_harmonics The solid harmonics of angular momentum 5.
/// @param g_harmonics The solid harmonics of angular momentum 4.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 13 rows and as many columns as the coordinates, with
/// the row of index m + 6 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570. The rows of all orders are formed in a
/// single loop, so that every value entering the recursion is loaded once.
auto make_i_solid_harmonics(const CSimdMatrix &h_harmonics, const CSimdMatrix &g_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum 7 of the vectors
/// between the atoms of the atom pairs.
/// @param i_harmonics The solid harmonics of angular momentum 6.
/// @param h_harmonics The solid harmonics of angular momentum 5.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 15 rows and as many columns as the coordinates, with
/// the row of index m + 7 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570. The rows of all orders are formed in a
/// single loop, so that every value entering the recursion is loaded once.
auto make_k_solid_harmonics(const CSimdMatrix &i_harmonics, const CSimdMatrix &h_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum 8 of the vectors
/// between the atoms of the atom pairs.
/// @param k_harmonics The solid harmonics of angular momentum 7.
/// @param i_harmonics The solid harmonics of angular momentum 6.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 17 rows and as many columns as the coordinates, with
/// the row of index m + 8 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570. The rows of all orders are formed in a
/// single loop, so that every value entering the recursion is loaded once.
auto make_l_solid_harmonics(const CSimdMatrix &k_harmonics, const CSimdMatrix &i_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

}  // namespace simdfunc

#endif /* SimdHarmonics_hpp */
