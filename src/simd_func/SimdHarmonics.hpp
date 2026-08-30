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

#include <vector>

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
/// @brief Creates the real solid harmonics of all angular momenta up to the given
/// one of the vectors between the atoms of the atom pairs.
/// @param coordinates The coordinates of the atom pairs, as seven rows.
/// @param lmax The highest angular momentum to create the harmonics of.
/// @return The vector of matrices, with the element of index l - 1 holding the
/// harmonics of angular momentum l as 2 l + 1 rows. The vector is empty for a
/// highest angular momentum below one, as the harmonic of angular momentum zero
/// is one for every atom pair and is not stored.
/// @note The harmonics of an angular momentum are formed from those of the two
/// angular momenta below it, so all of them are created on the way up and kept.
auto make_solid_harmonics(const CSimdMatrix &coordinates, const int lmax) -> std::vector<CSimdMatrix>;

auto make_p_solid_harmonics(const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum two of the
/// vectors between the atoms of the atom pairs.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 5 rows and as many columns as the coordinates, with
/// the row of index m + 2 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570.
auto make_d_solid_harmonics(const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum three of the
/// vectors between the atoms of the atom pairs.
/// @param d_harmonics The solid harmonics of angular momentum two.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 7 rows and as many columns as the coordinates, with
/// the row of index m + 3 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570.
auto make_f_solid_harmonics(const CSimdMatrix &d_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum four of the
/// vectors between the atoms of the atom pairs.
/// @param f_harmonics The solid harmonics of angular momentum three.
/// @param d_harmonics The solid harmonics of angular momentum two.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 9 rows and as many columns as the coordinates, with
/// the row of index m + 4 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570.
auto make_g_solid_harmonics(const CSimdMatrix &f_harmonics, const CSimdMatrix &d_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum five of the
/// vectors between the atoms of the atom pairs.
/// @param g_harmonics The solid harmonics of angular momentum four.
/// @param f_harmonics The solid harmonics of angular momentum three.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 11 rows and as many columns as the coordinates, with
/// the row of index m + 5 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570.
auto make_h_solid_harmonics(const CSimdMatrix &g_harmonics, const CSimdMatrix &f_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum six of the
/// vectors between the atoms of the atom pairs.
/// @param h_harmonics The solid harmonics of angular momentum five.
/// @param g_harmonics The solid harmonics of angular momentum four.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 13 rows and as many columns as the coordinates, with
/// the row of index m + 6 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570.
auto make_i_solid_harmonics(const CSimdMatrix &h_harmonics, const CSimdMatrix &g_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum seven of the
/// vectors between the atoms of the atom pairs.
/// @param i_harmonics The solid harmonics of angular momentum six.
/// @param h_harmonics The solid harmonics of angular momentum five.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 15 rows and as many columns as the coordinates, with
/// the row of index m + 7 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570.
auto make_k_solid_harmonics(const CSimdMatrix &i_harmonics, const CSimdMatrix &h_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum eight of the
/// vectors between the atoms of the atom pairs.
/// @param k_harmonics The solid harmonics of angular momentum seven.
/// @param i_harmonics The solid harmonics of angular momentum six.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 17 rows and as many columns as the coordinates, with
/// the row of index m + 8 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570.
auto make_l_solid_harmonics(const CSimdMatrix &k_harmonics, const CSimdMatrix &i_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum nine of the
/// vectors between the atoms of the atom pairs.
/// @param l_harmonics The solid harmonics of angular momentum eight.
/// @param k_harmonics The solid harmonics of angular momentum seven.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 19 rows and as many columns as the coordinates, with
/// the row of index m + 9 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570.
auto make_m_solid_harmonics(const CSimdMatrix &l_harmonics, const CSimdMatrix &k_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum ten of the
/// vectors between the atoms of the atom pairs.
/// @param m_harmonics The solid harmonics of angular momentum nine.
/// @param l_harmonics The solid harmonics of angular momentum eight.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 21 rows and as many columns as the coordinates, with
/// the row of index m + 10 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570.
auto make_n_solid_harmonics(const CSimdMatrix &m_harmonics, const CSimdMatrix &l_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum eleven of the
/// vectors between the atoms of the atom pairs.
/// @param n_harmonics The solid harmonics of angular momentum ten.
/// @param m_harmonics The solid harmonics of angular momentum nine.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 23 rows and as many columns as the coordinates, with
/// the row of index m + 11 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570.
auto make_o_solid_harmonics(const CSimdMatrix &n_harmonics, const CSimdMatrix &m_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

/// @brief Creates the real solid harmonics of angular momentum twelve of the
/// vectors between the atoms of the atom pairs.
/// @param o_harmonics The solid harmonics of angular momentum eleven.
/// @param n_harmonics The solid harmonics of angular momentum ten.
/// @param p_harmonics The solid harmonics of angular momentum one, which supply
/// the components of the vector between the atoms.
/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last
/// row supplies the squared distance of the atom pair.
/// @return The matrix of 25 rows and as many columns as the coordinates, with
/// the row of index m + 12 holding the solid harmonic of order m.
/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of
/// J. Chem. Theory Comput. 2020, 16, 2570.
auto make_q_solid_harmonics(const CSimdMatrix &o_harmonics, const CSimdMatrix &n_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix;

}  // namespace simdfunc

#endif /* SimdHarmonics_hpp */
