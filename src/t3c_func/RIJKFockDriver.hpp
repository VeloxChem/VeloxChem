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

#ifndef RIJKFockDriver_hpp
#define RIJKFockDriver_hpp

#include "SubMatrix.hpp"
#include "Matrix.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "T3FlatBuffer.hpp"

#include <string>

/// Class CRIJKFockDriver provides methods for computing Coulomb/Exchange Fock matrices
/// using three center electron repulsion integrals.
class CRIJKFockDriver
{
   public:
    /// Creates a Fock matrices  driver.
    CRIJKFockDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The Fock matrices driver to be copied.
    CRIJKFockDriver(const CRIJKFockDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The Fock matrices driver  to be moved.
    CRIJKFockDriver(CRIJKFockDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CRIJKFockDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The Fock matrices driver to be copy assigned.
    /// @return The assigned Fock matrices driver.
    auto operator=(const CRIJKFockDriver &other) -> CRIJKFockDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The Fock matrices driver to be move assigned.
    /// @return The assigned Fock matrices driver .
    auto operator=(CRIJKFockDriver &&other) noexcept -> CRIJKFockDriver & = delete;

    /// @brief The equality operator.
    /// @param other The Fock matrices driver  to be compared.
    /// @return True if Fock matrices drivers  are equal, False otherwise.
    auto operator==(const CRIJKFockDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The Fock matrices driver to be compared.
    /// @return True if Fock matrices drivers  are not equal, False otherwise.
    auto operator!=(const CRIJKFockDriver &other) const -> bool = delete;
    
    /// @brief Computes B^Q vectors in distributed form.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxilary molecular  basis.
    /// @param metric The B^Q vector transformation metric.
    /// @param rank The rank of MPI process to store batch of B^Q vectors.
    /// @param nodes The number of MPI nodes in communicator.
    auto compute_bq_vectors(const CMolecule&        molecule,
                            const CMolecularBasis&  basis,
                            const CMolecularBasis&  aux_basis,
                            const CSubMatrix&       metric,
                            const int               rank,
                            const int               nodes) -> void;
    
    /// @brief Computes Coulomb Fock matrix for given density.
    /// @param density The density matrix to construct Fock matrix.
    /// @param label The label of Fock matrix type.
    /// @return The Fock matrix.
    auto compute_j_fock(const CMatrix     &density,
                        const std::string &label) const -> CMatrix;
    
    /// @brief Computes Coulomb Fock matrix for given density.
    /// @param molorbs The occupied molecular orbitals to construct Fock matrix.
    /// @return The Fock matrix in submatrix storage.
    auto compute_k_fock(const CSubMatrix &molorbs) const -> CSubMatrix;
    
    private:
    
    /// @brief The distributed B^Q vectors.
    CT3FlatBuffer<double> _bq_vectors;
    
    /// @brief Computes M^P vector  for given density.
    /// @param density The density matrix to construct Fock matrix.
    /// @return The M^P vectors.
    auto _comp_m_vector(const CMatrix &density) const -> std::vector<double>;
    
    /// @brief Computes J vector  from given M^P vector.
    /// @param mvector The M^P vector.
    /// @return The computed J vector.
    auto _comp_j_vector(const std::vector<double>& mvector) const -> std::vector<double>;
};

#endif /* RIJKFockDriver_hpp */
