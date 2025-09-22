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

#ifndef RIFockDriver_hpp
#define RIFockDriver_hpp

#include <vector>
#include <map>

#include "SubMatrix.hpp"
#include "Matrix.hpp"
#include "T3FlatBuffer.hpp"
#include "T3RectFlatBuffer.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"

/// Class CRIFockDriver provides methods for computing Coulomb Fock matrices
/// using three center electron repulsion integrals.
class CRIFockDriver
{
   public:
    /// Creates a Fock matrices  driver.
    CRIFockDriver();
    
    /// Creates a Fock matrices  driver.
    /// @param j_metric The metric matrix for J fitting.
    CRIFockDriver(const CSubMatrix& j_metric);
    
    /// Creates a Fock matrices  driver.
    /// @param j_metric The metric matrix for J fitting.
    /// @param k_metric The metric matrix for K fitting.
    CRIFockDriver(const CSubMatrix& j_metric,
                  const CSubMatrix& k_metric);

    /// @brief The default copy constructor.
    /// @param other The Fock matrices driver to be copied.
    CRIFockDriver(const CRIFockDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The Fock matrices driver  to be moved.
    CRIFockDriver(CRIFockDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CRIFockDriver();

    /// @brief The default copy assignment operator.
    /// @param other The Fock matrices driver to be copy assigned.
    /// @return The assigned Fock matrices driver.
    auto operator=(const CRIFockDriver &other) -> CRIFockDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The Fock matrices driver to be move assigned.
    /// @return The assigned Fock matrices driver .
    auto operator=(CRIFockDriver &&other) noexcept -> CRIFockDriver & = delete;

    /// @brief The equality operator.
    /// @param other The Fock matrices driver  to be compared.
    /// @return True if Fock matrices drivers  are equal, False otherwise.
    auto operator==(const CRIFockDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The Fock matrices driver to be compared.
    /// @return True if Fock matrices drivers  are not equal, False otherwise.
    auto operator!=(const CRIFockDriver &other) const -> bool = delete;
    
    /// @brief Computes three center electron repulsion integral buffers.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxilary molecular  basis.
    auto prepare_buffers(const CMolecule       &molecule,
                         const CMolecularBasis &basis,
                         const CMolecularBasis &aux_basis) -> void;
    
    /// @brief Computes three center electron repulsion integral buffers.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxilary molecular  basis.
    /// @param atoms The vector of atoms to compute three-center electron repulsion integrals.
    auto prepare_buffers(const CMolecule&        molecule,
                         const CMolecularBasis&  basis,
                         const CMolecularBasis&  aux_basis,
                         const std::vector<int>& atoms) -> void;
    
    /// @brief Computes Fock matrix for given density.
    /// @param density The density matrix to construct Fock matrix.
    /// @param label The label of Fock matrix type.
    /// @return The Fock matrix.
    auto compute(const CMatrix     &density,
                 const std::string &label) const -> CMatrix;
    
    /// @brief Computes Fock matrix for given density.
    /// @param density The density matrix to construct Fock matrix.
    /// @param gvector The Gamma vector.
    /// @param label The label of Fock matrix type.
    /// @return The Fock matrix.
    auto compute(const CMatrix     &density,
                 const std::vector<double>& gvector,
                 const std::string &label) const -> CMatrix;
    
    /// @brief Computes local Fock matrix for given density.
    /// @param density The density matrix to construct Fock matrix.
    /// @param gvector The Gamma vector.
    /// @param label The label of Fock matrix type.
    /// @return The Fock matrix.
    auto local_compute(const CMatrix     &density,
                       const std::vector<double>& gvector,
                       const std::string &label) const -> CMatrix;
    
    /// @brief Computes transformed Gamma vector with J metric for given density.
    /// @param density The density matrix to construct Fock matrix.
    /// @return The transformed Gamma vector.
    auto compute_bq_vector(const CMatrix &density) const -> std::vector<double>;
    
    /// @brief Computes local transformed Gamma vector with J metric for given density.
    /// @param density The density matrix to construct Fock matrix.
    /// @return The transformed Gamma vector.
    auto compute_local_bq_vector(const CMatrix &density) const -> std::vector<double>;
    
    /// @brief Computes transformed Bq vector with K metric for given similarity transformed MOs.
    /// @param lambda_p The MOs augmented by particle single excitations.
    /// @param lambda_h The MOs augmented by hole single excitations.
    /// @return The transformed Gamma vector.
    auto compute_bq_vector(const CSubMatrix& lambda_p, const CSubMatrix& lambda_h) const -> CT3RectFlatBuffer<double>;
    
    /// @brief Gets mask indices of distributed auxilary AOs.
    /// @return The mask indices of distributed auxilary AOs.
    auto mask_indices() const -> std::map<size_t, size_t> { return _eri_buffer.mask_indices(); };
    
    private:
    /// @brief Pointer to metric matrix for J fitting.
    CSubMatrix* _j_metric;
    
    /// @brief Pointer to metric matrix for K fitting.
    CSubMatrix* _k_metric;
    
    /// @brief Three center electron repulsion integrals buffer.
    CT3FlatBuffer<double> _eri_buffer;
    
    /// @brief Computes Gamma vector  for given density.
    /// @param density The density matrix to construct Fock matrix.
    /// @return The Gamma vector.
    auto _comp_gamma_vector(const CMatrix &density) const -> std::vector<double>;
    
    /// @brief Transforms Gamma vector  with J metric.
    /// @param gvector The Gamma vector.
    /// @return The transformed Gamma vector.
    auto _trafo_gamma_vector(const std::vector<double>& gvector) const -> std::vector<double>;
    
    /// @brief Transforms local Gamma vector  with J metric.
    /// @param gvector The Gamma vector.
    /// @return The transformed Gamma vector.
    auto _trafo_local_gamma_vector(const std::vector<double>& gvector) const -> std::vector<double>;
    
    /// @brief Computes J vector  for given transformed Gamma vector.
    /// @param gvector The transformed Gamma vector.
    /// @return The computed J vector.
    auto _comp_j_vector(const std::vector<double>& gvector) const -> std::vector<double>;
    
    /// @brief Computes local J vector  for given transformed Gamma vector.
    /// @param gvector The transformed Gamma vector.
    /// @return The computed J vector.
    auto _comp_local_j_vector(const std::vector<double>& gvector) const -> std::vector<double>;
};

#endif /* RIFockDriver_hpp */
