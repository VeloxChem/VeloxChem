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

#ifndef ThreeCenterElectronRepulsionGeom0X0Driver_hpp
#define ThreeCenterElectronRepulsionGeom0X0Driver_hpp

#include <vector>
#include <ranges>

#include "GtoFunc.hpp"
#include "GtoPairBlockFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"
#include "TensorComponents.hpp"
#include "T3CUtils.hpp"
#include "T3RectFlatBuffer.hpp"
#include "T3CGeom0X0Distributor.hpp"
#include "ThreeCenterElectronRepulsionGeom010Func.hpp"

#include <iostream>

/// @brief Class  CThreeCenterElectronRepulsionGeom0X0Driver provides methods for computing arbitrary order three-center
/// electron repulsion integral derivatives with respect bra side.
template <int N>
class CThreeCenterElectronRepulsionGeom0X0Driver
{
   public:
    /// @brief Creates an electron repulsion derivative integrals driver.
    CThreeCenterElectronRepulsionGeom0X0Driver() = default;

    /// @brief The default copy constructor.
    /// @param other The overlap derivative integrals driver to be copied.
    CThreeCenterElectronRepulsionGeom0X0Driver(const CThreeCenterElectronRepulsionGeom0X0Driver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The overlap derivative integrals driver  to be moved.
    CThreeCenterElectronRepulsionGeom0X0Driver(CThreeCenterElectronRepulsionGeom0X0Driver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CThreeCenterElectronRepulsionGeom0X0Driver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The overlap derivative integrals driver to be copy assigned.
    /// @return The assigned overlap derivative integrals driver.
    auto operator=(const CThreeCenterElectronRepulsionGeom0X0Driver &other) -> CThreeCenterElectronRepulsionGeom0X0Driver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The overlap derivative integrals driver to be move assigned.
    /// @return The assigned overlap derivative integrals driver .
    auto operator=(CThreeCenterElectronRepulsionGeom0X0Driver &&other) noexcept -> CThreeCenterElectronRepulsionGeom0X0Driver & = delete;

    /// @brief The equality operator.
    /// @param other The overlap derivative integrals driver  to be compared.
    /// @return True if overlap derivative integrals drivers  are equal, False otherwise.
    auto operator==(const CThreeCenterElectronRepulsionGeom0X0Driver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The overlap derivative integrals driver to be compared.
    /// @return True if overlap derivative integrals drivers  are not equal, False otherwise.
    auto operator!=(const CThreeCenterElectronRepulsionGeom0X0Driver &other) const -> bool = delete;

    /// @brief Computes electron repulsion matrix derivative fro molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxilary molecular basis for fiting of four-center repulsion integrals.
    /// @param molecule The molecule.
    /// @param iatom The index of atom.
    /// @return The nuclear potential matrix.
    auto compute(const CMolecularBasis &basis,
                 const CMolecularBasis &aux_basis,
                 const CMolecule       &molecule,
                 const int             iatom) const -> CT3RectFlatBuffer<double>;
};

template <int N>
auto
CThreeCenterElectronRepulsionGeom0X0Driver<N>::compute(const CMolecularBasis &basis,
                                                       const CMolecularBasis &aux_basis,
                                                       const CMolecule       &molecule,
                                                       const int             iatom) const -> CT3RectFlatBuffer<double>
{
    // set up GTOs data
    
    const auto bra_gto_blocks = gtofunc::make_gto_blocks(aux_basis, molecule);
    
    const auto ketr_gto_blocks = gtofunc::make_gto_blocks(basis, molecule, {iatom, });
    
    const auto ketf_gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto gto_pair_blocks = gtofunc::make_gto_pair_blocks(ketr_gto_blocks, ketf_gto_blocks);
    
    // set up composite flat tensor for integrals
    
    std::vector<size_t> aux_indices(aux_basis.dimensions_of_basis());
    
    std::iota(aux_indices.begin(), aux_indices.end(), size_t{0});
    
    const auto red_indices = t3cfunc::unique_indices(ketr_gto_blocks);
    
    std::vector<size_t> full_indices(basis.dimensions_of_basis());
    
    std::iota(full_indices.begin(), full_indices.end(), size_t{0});
    
    const auto mask_indices = t3cfunc::mask_indices(ketr_gto_blocks);
    
    CT3RectFlatBuffer<double> buffer(aux_indices, mask_indices, basis.dimensions_of_basis(), 3);
    
    // set distributor
    
    CT3CGeom0X0Distributor distributor(&buffer);
    
    // main compute loop
    
    std::ranges::for_each(bra_gto_blocks, [&](const auto& gblock) {
        auto bra_range = std::pair<size_t, size_t>{size_t{0}, gblock.number_of_basis_functions()};
        std::ranges::for_each(gto_pair_blocks, [&](const auto& gp_pairs) {
            if constexpr (N == 1)
            {
                t3cerifunc::compute_geom_010(distributor, gblock, gp_pairs, bra_range);
            }
        });
    });
    
    return buffer;
}

#endif /* ThreeCenterElectronRepulsionGeom0X0Driver_hpp */
