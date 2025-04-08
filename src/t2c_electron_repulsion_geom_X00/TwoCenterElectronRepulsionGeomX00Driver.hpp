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

#ifndef TwoCenterElectronRepulsionGeomX00Driver_hpp
#define TwoCenterElectronRepulsionGeomX00Driver_hpp

#include <vector>

#include "GtoFunc.hpp"
#include "Matrices.hpp"
#include "MatricesFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "OpenMPFunc.hpp"
#include "Point.hpp"
#include "T2CDistributor.hpp"
#include "T2CGeom10SumDistributor.hpp"
#include "T2CGeom10SumTwoDistributor.hpp"
#include "TwoCenterElectronRepulsionGeom100Func.hpp"

/// @brief Class  CTwoCenterElectronRepulsionGeomX00Driver provides methods for computing arbitrary order two-center
/// electron repulsion integral derivatives with respect bra side.
template <int N>
class CTwoCenterElectronRepulsionGeomX00Driver
{
   public:
    /// @brief Creates an electron repulsion derivative integrals driver.
    CTwoCenterElectronRepulsionGeomX00Driver() = default;

    /// @brief The default copy constructor.
    /// @param other The overlap derivative integrals driver to be copied.
    CTwoCenterElectronRepulsionGeomX00Driver(const CTwoCenterElectronRepulsionGeomX00Driver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The overlap derivative integrals driver  to be moved.
    CTwoCenterElectronRepulsionGeomX00Driver(CTwoCenterElectronRepulsionGeomX00Driver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CTwoCenterElectronRepulsionGeomX00Driver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The overlap derivative integrals driver to be copy assigned.
    /// @return The assigned overlap derivative integrals driver.
    auto operator=(const CTwoCenterElectronRepulsionGeomX00Driver &other) -> CTwoCenterElectronRepulsionGeomX00Driver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The overlap derivative integrals driver to be move assigned.
    /// @return The assigned overlap derivative integrals driver .
    auto operator=(CTwoCenterElectronRepulsionGeomX00Driver &&other) noexcept -> CTwoCenterElectronRepulsionGeomX00Driver & = delete;

    /// @brief The equality operator.
    /// @param other The overlap derivative integrals driver  to be compared.
    /// @return True if overlap derivative integrals drivers  are equal, False otherwise.
    auto operator==(const CTwoCenterElectronRepulsionGeomX00Driver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The overlap derivative integrals driver to be compared.
    /// @return True if overlap derivative integrals drivers  are not equal, False otherwise.
    auto operator!=(const CTwoCenterElectronRepulsionGeomX00Driver &other) const -> bool = delete;

    /// @brief Computes electron repulsion matrix derivative fro molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param iatom The index of atom.
    /// @return The nuclear potential matrix.
    auto compute(const CMolecularBasis             &basis,
                 const CMolecule                   &molecule,
                 const int                         iatom) const -> CMatrices;
    
    /// @brief Computes electron repulsion matrix derivative fro molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param gamma The transformed Gamma vector.
    /// @param iatom The index of atom.
    /// @return The electron repulsion gradient.
    auto compute(const CMolecularBasis             &basis,
                 const CMolecule                   &molecule,
                 const std::vector<double>&        gamma,
                 const int                         iatom) const -> std::array<double, 3>;
    
    /// @brief Computes electron repulsion matrix derivative fro molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param bra_gamma The transformed Gamma vector on bra side.
    /// @param ket_gamma The transformed Gamma vector on ket side.
    /// @param iatom The index of atom.
    /// @return The electron repulsion gradient.
    auto compute(const CMolecularBasis             &basis,
                 const CMolecule                   &molecule,
                 const std::vector<double>&        bra_gamma,
                 const std::vector<double>&        ket_gamma,
                 const int                         iatom) const -> std::array<double, 3>;
};

template <int N>
auto
CTwoCenterElectronRepulsionGeomX00Driver<N>::compute(const CMolecularBasis             &basis,
                                                     const CMolecule                   &molecule,
                                                     const int                         iatom) const -> CMatrices
{
    // set up operator derivatives matrices

    auto eri_mats = matfunc::make_matrices(std::array<int, 1>{N}, basis, mat_t::general);

    eri_mats.zero();

    // prepare pointers for OMP parallel region

    auto ptr_basis = &basis;

    auto ptr_molecule = &molecule;

    auto ptr_eri_mats = &eri_mats;

    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_basis, ptr_molecule, ptr_eri_mats, iatom)
    {
#pragma omp single nowait
        {
            const auto bra_gto_blocks = gtofunc::make_gto_blocks(*ptr_basis, *ptr_molecule, {iatom, });
            
            const auto ket_gto_blocks = gtofunc::make_gto_blocks(*ptr_basis, *ptr_molecule);

            const auto tasks = omp::make_work_tasks(bra_gto_blocks, ket_gto_blocks);

            std::ranges::for_each(std::ranges::reverse_view(tasks), [&](const auto& task) {
                auto bra_gtos    = bra_gto_blocks[task[0]];
                auto ket_gtos    = ket_gto_blocks[task[1]];
                auto bra_indices = std::pair<size_t, size_t>{task[2], task[3]};
                auto ket_indices = std::pair<size_t, size_t>{task[4], task[5]};
#pragma omp task firstprivate(bra_gtos, ket_gtos, bra_indices, ket_indices)
                {
                    CT2CDistributor<CMatrices> distributor(ptr_eri_mats);
                    if constexpr (N == 1)
                    {
                        t2cerifunc::compute_geom_100(distributor, bra_gtos, ket_gtos, bra_indices, ket_indices, false);
                    }
                }
            });
        }
    }

    return eri_mats;
}

template <int N>
auto
CTwoCenterElectronRepulsionGeomX00Driver<N>::compute(const CMolecularBasis             &basis,
                                                     const CMolecule                   &molecule,
                                                     const std::vector<double>&        gamma,
                                                     const int                         iatom) const -> std::array<double, 3>
{
    // set up bra and ket GTOs blocks
    
    const auto bra_gto_blocks = gtofunc::make_gto_blocks(basis, molecule, {iatom, });
    
    const auto ket_gto_blocks = gtofunc::make_gto_blocks(basis, molecule);
    
    // set up tasks

    const auto tasks = omp::make_work_tasks(bra_gto_blocks, ket_gto_blocks);
    
    // set up accumulation vector
    
    auto gvec_xyz = std::vector<std::array<double, 3>>(tasks.size(), {0.0, 0.0, 0.0});
    
    // prepare pointers for OMP parallel region
    
    auto ptr_bra_gto_blocks = bra_gto_blocks.data();
    
    auto ptr_ket_gto_blocks = ket_gto_blocks.data();
    
    auto ptr_gamma = gamma.data();
    
    auto ptr_tasks = &tasks;
    
    auto ptr_gvec_xyz = gvec_xyz.data();

    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_bra_gto_blocks, ptr_ket_gto_blocks, ptr_tasks, ptr_gvec_xyz)
    {
#pragma omp single nowait
        {
            size_t index = 0;
            std::ranges::for_each(std::ranges::reverse_view(*ptr_tasks), [&](const auto& task) {
                auto bra_idx     = task[0];
                auto ket_idx     = task[1];
                auto bra_indices = std::pair<size_t, size_t>{task[2], task[3]};
                auto ket_indices = std::pair<size_t, size_t>{task[4], task[5]};
#pragma omp task firstprivate(bra_idx, ket_idx, bra_indices, ket_indices)
                {
                    auto bra_gtos = ptr_bra_gto_blocks[bra_idx];
                    auto ket_gtos = ptr_ket_gto_blocks[ket_idx];
                    CT2CGeom10SumDistributor distributor(ptr_gvec_xyz[index].data(), ptr_gamma);
                    if constexpr (N == 1)
                    {
                        t2cerifunc::compute_geom_100(distributor, bra_gtos, ket_gtos, bra_indices, ket_indices, false);
                    }
                }
                index++;
            });
        }
    }

    // sum up contributions from accumulation vector
    
    auto g_xyz = std::array<double, 3>({0.0, 0.0, 0.0});
    
    for (const auto& t_xyz : gvec_xyz)
    {
        g_xyz[0] += t_xyz[0];
        
        g_xyz[1] += t_xyz[1];
        
        g_xyz[2] += t_xyz[2];
    }
    
    return g_xyz;
}

template <int N>
auto
CTwoCenterElectronRepulsionGeomX00Driver<N>::compute(const CMolecularBasis             &basis,
                                                     const CMolecule                   &molecule,
                                                     const std::vector<double>&        bra_gamma,
                                                     const std::vector<double>&        ket_gamma,
                                                     const int                         iatom) const -> std::array<double, 3>
{
    // set up bra and ket GTOs blocks
    
    const auto bra_gto_blocks = gtofunc::make_gto_blocks(basis, molecule, {iatom, });
    
    const auto ket_gto_blocks = gtofunc::make_gto_blocks(basis, molecule);
    
    // set up tasks

    const auto tasks = omp::make_work_tasks(bra_gto_blocks, ket_gto_blocks);
    
    // set up accumulation vector
    
    auto gvec_xyz = std::vector<std::array<double, 3>>(tasks.size(), {0.0, 0.0, 0.0});
    
    // prepare pointers for OMP parallel region
    
    auto ptr_bra_gto_blocks = bra_gto_blocks.data();
    
    auto ptr_ket_gto_blocks = ket_gto_blocks.data();
    
    auto ptr_bra_gamma = bra_gamma.data();
    
    auto ptr_ket_gamma = ket_gamma.data();
    
    auto ptr_tasks = &tasks;
    
    auto ptr_gvec_xyz = gvec_xyz.data();

    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_bra_gto_blocks, ptr_ket_gto_blocks, ptr_tasks, ptr_gvec_xyz)
    {
#pragma omp single nowait
        {
            size_t index = 0;
            std::ranges::for_each(std::ranges::reverse_view(*ptr_tasks), [&](const auto& task) {
                auto bra_idx     = task[0];
                auto ket_idx     = task[1];
                auto bra_indices = std::pair<size_t, size_t>{task[2], task[3]};
                auto ket_indices = std::pair<size_t, size_t>{task[4], task[5]};
#pragma omp task firstprivate(bra_idx, ket_idx, bra_indices, ket_indices)
                {
                    auto bra_gtos = ptr_bra_gto_blocks[bra_idx];
                    auto ket_gtos = ptr_ket_gto_blocks[ket_idx];
                    CT2CGeom10SumTwoDistributor distributor(ptr_gvec_xyz[index].data(), ptr_bra_gamma, ptr_ket_gamma);
                    if constexpr (N == 1)
                    {
                        t2cerifunc::compute_geom_100(distributor, bra_gtos, ket_gtos, bra_indices, ket_indices, false);
                    }
                }
                index++;
            });
        }
    }

    // sum up contributions from accumulation vector
    
    auto g_xyz = std::array<double, 3>({0.0, 0.0, 0.0});
    
    for (const auto& t_xyz : gvec_xyz)
    {
        g_xyz[0] += t_xyz[0];
        
        g_xyz[1] += t_xyz[1];
        
        g_xyz[2] += t_xyz[2];
    }
    
    return g_xyz;
}



#endif /* TwoCenterElectronRepulsionX00Driver_hpp */
