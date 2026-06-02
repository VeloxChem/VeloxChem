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

#ifndef newints_CauchySchwarzPairData_hpp
#define newints_CauchySchwarzPairData_hpp

#include <cstddef>
#include <vector>

namespace newints {

/// @brief A significant orbital shell pair and its Cauchy-Schwarz factor.
struct ShellPair
{
    /// @brief First contracted-shell index (atom-major), canonical i >= j.
    int i;
    /// @brief Second contracted-shell index (atom-major).
    int j;
    /// @brief Q_ij = max over (mu in shell i, nu in shell j) of sqrt((mu nu|mu nu)).
    double q;
};

/// @brief Cauchy-Schwarz factors for unique orbital shell pairs (the bra/ket side of
/// three- and four-center electron-repulsion integrals).
///
/// For each unique contracted-shell pair (i >= j) the factor is
///   Q_ij = max over (mu in shell i, nu in shell j) of sqrt((mu nu|mu nu)),
/// which bounds later integrals by Cauchy-Schwarz: |(mu nu|lambda sigma)| <= Q_ij Q_kl
/// (four-center) and |(mu nu|P)| <= Q_ij Q_P (three-center, with Q_P from
/// CauchySchwarzData).
///
/// Because the charge distribution mu nu decays with the mu-nu separation, the vast
/// majority of shell pairs are negligible, so only SIGNIFICANT pairs (Q_ij >= the build
/// threshold) are stored, as a sparse list -- O(N_shells) entries rather than the dense
/// O(N_shells^2) triangle. Shell indices are the same atom-major contracted-shell indices
/// as MolecularBasisOutline / CauchySchwarzData.
///
/// This is a pure storage container; the factors are produced by a separate builder
/// (see ShellPairSchwarzDriver), since (mu nu|mu nu) requires four-center ERI-diagonal
/// kernels.
class CauchySchwarzPairData
{
   public:
    /// @brief Creates empty screening data.
    CauchySchwarzPairData() = default;

    /// @brief Creates screening data from a precomputed list of significant pairs.
    /// @param nshells The total number of contracted shells (the index range of i, j).
    /// @param threshold The build threshold below which pairs were discarded.
    /// @param pairs The significant pairs (i >= j); sorted on construction.
    CauchySchwarzPairData(std::size_t nshells, double threshold, std::vector<ShellPair> pairs);

    /// @brief Gets the significant shell pairs (i >= j), sorted by (i, j).
    /// @return The list of significant pairs.
    auto pairs() const -> const std::vector<ShellPair> & { return pairs_; }

    /// @brief Gets the number of significant shell pairs.
    /// @return The number of stored pairs.
    auto number_of_pairs() const -> std::size_t { return pairs_.size(); }

    /// @brief Gets the total number of contracted shells.
    /// @return The shell-index range.
    auto number_of_shells() const -> std::size_t { return nshells_; }

    /// @brief Gets the build threshold used to filter significant pairs.
    /// @return The threshold.
    auto threshold() const -> double { return threshold_; }

    /// @brief Gets the Cauchy-Schwarz factor Q_ij of a shell pair.
    /// @param i The first shell index.
    /// @param j The second shell index (order does not matter; canonicalized to i >= j).
    /// @return Q_ij if the pair is significant, 0 otherwise.
    auto q(std::size_t i, std::size_t j) const -> double;

   private:
    /// @brief Total number of contracted shells.
    std::size_t nshells_ = 0;
    /// @brief Build threshold below which pairs were discarded.
    double threshold_ = 0.0;
    /// @brief Significant pairs (i >= j), sorted lexicographically by (i, j).
    std::vector<ShellPair> pairs_;
};

}  // namespace newints

#endif /* newints_CauchySchwarzPairData_hpp */
