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



#ifndef SimdThreeCenterElectronRepulsionDriver_hpp
#define SimdThreeCenterElectronRepulsionDriver_hpp

#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "AtomBasisTripleSparsity.hpp"
#include "SimdMatrix.hpp"
#include "SparseTensor.hpp"

/// @brief Class CSimdThreeCenterElectronRepulsionDriver computes the three-center
/// electron repulsion integrals of a molecular basis on the a and b sides and an
/// auxiliary molecular basis on the c side, and stores them in a sparse tensor.
/// @note The integrals are screened, unlike the two-center ones. The bound of a
/// combination of basis functions falls with the separation of the atom pair on
/// the a and b sides, so a pair of distant atoms contributes nothing to the
/// auxiliary functions it does not reach. What survives is described by
/// CSparseTensor, which is sparse along the a and b sides and dense along the c
/// side, and the driver fills the values blocks it lays out.
/// @note The threshold is an argument of compute rather than a property of the
/// driver, as the tensor of one molecule is built at several thresholds and the
/// bound is what decides its storage.
/// @note The atoms on the c side may be given, in which case the tensor holds
/// the part of the whole which those atoms carry. This is how a tensor too large
/// for the memory is formed in batches: the caller divides the atoms and the
/// parts are computed one after another.
class CSimdThreeCenterElectronRepulsionDriver
{
   public:
    /// @brief The default constructor.
    CSimdThreeCenterElectronRepulsionDriver() = default;

    /// @brief Computes the three-center electron repulsion integrals of a
    /// molecular basis and an auxiliary molecular basis, for all atoms on c side.
    /// @param molecule The molecule to compute the integrals of.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis on c side.
    /// @param threshold The screening threshold.
    /// @return The symmetric sparse tensor of the integrals.
    auto compute(const CMolecule       &molecule,
                 const CMolecularBasis &basis,
                 const CMolecularBasis &aux_basis,
                 const double           threshold) const -> CSparseTensor;

    /// @brief Computes the three-center electron repulsion integrals of a
    /// molecular basis and an auxiliary molecular basis, for the given atoms on
    /// c side.
    /// @param molecule The molecule to compute the integrals of.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis on c side.
    /// @param threshold The screening threshold.
    /// @param atoms The atoms on c side to compute the integrals of, as their
    /// indices in the molecule and without repetition.
    /// @return The symmetric sparse tensor of the integrals.
    auto compute(const CMolecule        &molecule,
                 const CMolecularBasis  &basis,
                 const CMolecularBasis  &aux_basis,
                 const double            threshold,
                 const std::vector<int> &atoms) const -> CSparseTensor;

   private:
    /// @brief Computes the integrals of the blocks of a sparse tensor.
    /// @param tensor The sparse tensor to compute the integrals of.
    /// @param molecule The molecule to take the atomic coordinates from.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis on c side.
    auto _compute_blocks(CSparseTensor         &tensor,
                         const CMolecule       &molecule,
                         const CMolecularBasis &basis,
                         const CMolecularBasis &aux_basis) const -> void;

};

#endif /* SimdThreeCenterElectronRepulsionDriver_hpp */
