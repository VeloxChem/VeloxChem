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

#ifndef newints_CauchySchwarzData_hpp
#define newints_CauchySchwarzData_hpp

#include <cstddef>
#include <vector>

class CMolecularBasis;

namespace newints {

/// @brief Per-shell Cauchy-Schwarz factors for the two-center Coulomb metric.
///
/// For each contracted shell P it stores the single number
///   Q_P = sqrt((P|P)),
/// where (P|P) is the two-center Coulomb self-integral. These bound later integrals
/// by the Cauchy-Schwarz inequality, e.g. |(P|Q)| <= Q_P Q_Q for the two-center
/// metric and |(mu nu|P)| <= sqrt((mu nu|mu nu)) Q_P for three-center RI integrals.
///
/// One double is stored per contracted shell. Because the concentric block (P|P) is
/// diagonal in m and independent of m for solid-harmonic Gaussians, the per-component
/// self-integrals are all equal, so this single value is the exact max over the
/// shell's components -- not a lossy approximation.
///
/// The values are laid out in the same atom-major shell order as MolecularBasisOutline
/// (its indices()/angular_momenta()/atom_indices()), so shell s maps to the SparseMatrix
/// key outline.indices()[s].
class CauchySchwarzData
{
   public:
    /// @brief Creates empty screening data.
    CauchySchwarzData() = default;

    /// @brief Builds Q_P = sqrt((P|P)) for every contracted shell of the basis from
    /// the concentric two-center Coulomb self-integral.
    /// @param basis The molecular basis.
    explicit CauchySchwarzData(const CMolecularBasis &basis);

    /// @brief Gets the Cauchy-Schwarz factor Q_P of a shell.
    /// @param shell The atom-major contracted-shell index.
    /// @return Q_P = sqrt((P|P)).
    auto q(const std::size_t shell) const -> double { return q_[shell]; }

    /// @brief Gets all Cauchy-Schwarz factors in atom-major shell order.
    /// @return The vector of Q_P values.
    auto values() const -> const std::vector<double> & { return q_; }

    /// @brief Gets the number of contracted shells.
    /// @return The number of stored factors.
    auto size() const -> std::size_t { return q_.size(); }

   private:
    /// @brief Q_P = sqrt((P|P)) per contracted shell, atom-major order.
    std::vector<double> q_;
};

}  // namespace newints

#endif /* newints_CauchySchwarzData_hpp */
