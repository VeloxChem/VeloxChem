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


#ifndef SimdKineticEnergyDriver_hpp
#define SimdKineticEnergyDriver_hpp

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "SparseMatrix.hpp"

/// @brief Class CSimdKineticEnergyDriver computes the two-center overlap integrals of
/// a molecular basis and stores them in a sparse matrix, using the sparsity
/// patterns of the atom basis pair groups to skip the atom pairs and the
/// combinations of basis functions whose integrals are below the threshold.
class CSimdKineticEnergyDriver
{
   public:
    /// @brief The constructor with screening threshold.
    /// @param threshold The screening threshold of the integrals.
    explicit CSimdKineticEnergyDriver(const double threshold = 1.0e-14)

        : _threshold(threshold)
    {
    }

    /// @brief Computes the kinetic energy matrix of a molecular basis.
    /// @param molecule The molecule to compute the kinetic energy matrix of.
    /// @param basis The molecular basis on bra and ket sides.
    /// @return The symmetric sparse kinetic energy matrix.
    auto compute(const CMolecule &molecule, const CMolecularBasis &basis) const -> CSparseMatrix;

    /// @brief Computes the kinetic energy matrix of a pair of molecular bases.
    /// @param molecule The molecule to compute the kinetic energy matrix of.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    /// @return The general sparse kinetic energy matrix.
    auto compute(const CMolecule &molecule, const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis) const
        -> CSparseMatrix;

    /// @brief Gets screening threshold of the integrals.
    /// @return The screening threshold.
    auto
    get_threshold() const -> double
    {
        return _threshold;
    }

   private:
    /// @brief Computes the integrals of the off-diagonal blocks of a sparse matrix.
    /// @param matrix The sparse matrix to compute the integrals of.
    /// @param molecule The molecule to compute the integrals of.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    auto _compute_pair_blocks(CSparseMatrix         &matrix,
                              const CMolecule       &molecule,
                              const CMolecularBasis &bra_basis,
                              const CMolecularBasis &ket_basis) const -> void;

    /// @brief Computes the integrals of the diagonal blocks of a sparse matrix.
    /// @param matrix The sparse matrix to compute the integrals of.
    /// @param molecule The molecule to compute the integrals of.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    auto _compute_diagonal_blocks(CSparseMatrix         &matrix,
                                  const CMolecule       &molecule,
                                  const CMolecularBasis &bra_basis,
                                  const CMolecularBasis &ket_basis) const -> void;

    /// @brief The screening threshold of the integrals.
    double _threshold;
};

#endif /* SimdKineticEnergyDriver_hpp */
