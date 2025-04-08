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

#ifndef RIFockGradDriver_hpp
#define RIFockGradDriver_hpp

#include <vector>

#include "SubMatrix.hpp"
#include "Matrix.hpp"
#include "T3FlatBuffer.hpp"
#include "Point.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"

/// Class CRIFockgradDriver provides methods for computing Coulomb Fock gradient
/// using three center electron repulsion integrals.
class CRIFockGradDriver
{
   public:
    /// Creates a Fock matrices  driver.
    CRIFockGradDriver() = default;
    
    /// @brief The default copy constructor.
    /// @param other The Fock matrices driver to be copied.
    CRIFockGradDriver(const CRIFockGradDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The Fock matrices driver  to be moved.
    CRIFockGradDriver(CRIFockGradDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CRIFockGradDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The Fock matrices driver to be copy assigned.
    /// @return The assigned Fock matrices driver.
    auto operator=(const CRIFockGradDriver &other) -> CRIFockGradDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The Fock matrices driver to be move assigned.
    /// @return The assigned Fock matrices driver .
    auto operator=(CRIFockGradDriver &&other) noexcept -> CRIFockGradDriver & = delete;

    /// @brief The equality operator.
    /// @param other The Fock matrices driver  to be compared.
    /// @return True if Fock matrices drivers  are equal, False otherwise.
    auto operator==(const CRIFockGradDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The Fock matrices driver to be compared.
    /// @return True if Fock matrices drivers  are not equal, False otherwise.
    auto operator!=(const CRIFockGradDriver &other) const -> bool = delete;
    
    /// @brief Computes Fock matrix gradient.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxilary molecular basis for fiting of four-center repulsion integrals.
    /// @param molecule The molecule.
    /// @param gamma The transformed Gamma vector.
    /// @param density The density matrix to construct Fock matrix.
    /// @param iatom The index of requested atom.
    /// @return The Fock contribution to atom's gradient.
    auto compute(const CMolecularBasis&     basis,
                 const CMolecularBasis&     aux_basis,
                 const CMolecule&           molecule,
                 const std::vector<double>& gamma,
                 const CMatrix&             density,
                 const int                  iatom) const -> TPoint<double>;
    
    
    
    /// @brief Computes Fock matrix gradient.
    /// @param screener The screener with basis function pairs data.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxilary molecular basis for fiting of four-center repulsion integrals.
    /// @param molecule The molecule.
    /// @param bra_gamma The transformed Gamma vector on bra side.
    /// @param ket_gamma The transformed Gamma vector on ket side.
    /// @param bra_density The density matrix to construct Fock matrix on bra side.
    /// @param ket_density The density matrix to construct Fock matrix on ket side.
    /// @param iatom The index of requested atom.
    /// @return The Fock contribution to atom's gradient.
    auto direct_compute(const CT4CScreener&        screener,
                        const CMolecularBasis&     basis,
                        const CMolecularBasis&     aux_basis,
                        const CMolecule&           molecule,
                        const std::vector<double>& bra_gamma,
                        const std::vector<double>& ket_gamma,
                        const CMatrix&             bra_density,
                        const CMatrix&             ket_density,
                        const int                  iatom,
                        const int                  ithreshold) const -> TPoint<double>;
    
    /// @brief Computes Fock matrix gradient.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxilary molecular basis for fiting of four-center repulsion integrals.
    /// @param molecule The molecule.
    /// @param gamma The transformed Gamma vector.
    /// @param density The density matrix to construct Fock matrix.
    /// @param atoms The indices of requested atoms.
    /// @return The Fock contribution to atom's gradient.
    auto compute(const CMolecularBasis&     basis,
                 const CMolecularBasis&     aux_basis,
                 const CMolecule&           molecule,
                 const std::vector<double>& gamma,
                 const CMatrix&             density,
                 const std::vector<int>     atoms) const -> std::vector<TPoint<double>>;
    
    private:
    
    /// @brief Computes Coulomb contribution to atom's gradient from resolution of identity integrals.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxilary molecular basis for fiting of four-center repulsion integrals.
    /// @param molecule The molecule.
    /// @param gamma The transformed Gamma vector.
    /// @param density The density matrix to construct Fock matrix.
    /// @param iatom The index of requested atom.
    /// @return The Fock contribution to atom's gradient.
    auto _comp_eri_grad(const CMolecularBasis&     basis,
                        const CMolecularBasis&     aux_basis,
                        const CMolecule&           molecule,
                        const std::vector<double>& gamma,
                        const CMatrix&             density,
                        const int                  iatom) const -> std::array<double, 3>;
};

#endif /* RIFockGradDriver_hpp */
