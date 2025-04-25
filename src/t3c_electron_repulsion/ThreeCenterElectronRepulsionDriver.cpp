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

#include "ThreeCenterElectronRepulsionDriver.hpp"

#include <algorithm>

#include "GtoFunc.hpp"
#include "GtoPairBlockFunc.hpp"
#include "OpenMPFunc.hpp"
#include "ThreeCenterElectronRepulsionFunc.hpp"
#include "T3CDistributor.hpp"
#include "T3CLocalDistributor.hpp"

#include <iostream>

auto
CThreeCenterElectronRepulsionDriver::compute(const CMolecularBasis &basis,
                                             const CMolecularBasis &aux_basis,
                                             const CMolecule &molecule) const -> CT3FlatBuffer<double>
{
    std::vector<size_t> aux_indices(aux_basis.dimensions_of_basis());
    
    std::iota(aux_indices.begin(), aux_indices.end(), size_t{0});
    
    CT3FlatBuffer<double> buffer(aux_indices, basis.dimensions_of_basis());
    
    // set up basis function pairs blocks

    const auto gto_pair_blocks = gtofunc::make_gto_pair_blocks(basis, molecule);
    
    const auto aux_gto_blocks = gtofunc::make_gto_blocks(aux_basis, molecule);
    
    // prepare pointers for OMP parallel region
    
    auto ptr_gto_pair_blocks = &gto_pair_blocks;
    
    auto ptr_aux_gto_blocks = &aux_gto_blocks;
    
    auto ptr_buffer = &buffer;
    
    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_gto_pair_blocks, ptr_aux_gto_blocks, ptr_buffer)
    {
#pragma omp single nowait
        {
            const auto tasks = omp::make_aux_work_tasks(*ptr_aux_gto_blocks);

            std::ranges::for_each(std::ranges::reverse_view(tasks), [&](const auto& task) {
                auto aux_idx = task[0];
                auto bra_range = std::pair<size_t, size_t>{task[1], task[2]};
                if (const auto nblocks = ptr_gto_pair_blocks->size(); nblocks > 0)
                {
                    for (size_t i = 0; i < nblocks; i++)
                    {
#pragma omp task firstprivate(aux_idx, bra_range, i)
                        {
                            auto aux_gtos = ptr_aux_gto_blocks->at(aux_idx);
                            
                            const auto gp_pairs = ptr_gto_pair_blocks->at(i);
                            
                            CT3CDistributor distributor(ptr_buffer);
                            
                            t3cerifunc::compute(distributor, aux_gtos, gp_pairs, bra_range);
                        }
                    }
                }
            });
        }
    }

    return buffer;
}

auto
CThreeCenterElectronRepulsionDriver::compute(const CMolecularBasis &basis,
                                             const CMolecularBasis &aux_basis,
                                             const CMolecule &molecule,
                                             const std::vector<int>& atoms) const -> CT3FlatBuffer<double>
{
    // set up GTOs data
    
    const auto aux_gto_blocks = gtofunc::make_gto_blocks(aux_basis, molecule, atoms);
    
    const auto gto_pair_blocks = gtofunc::make_gto_pair_blocks(basis, molecule);
    
    // set up composite flat tensor for integrals
    
    const auto mask_indices = t3cfunc::mask_indices(aux_gto_blocks);
    
    CT3FlatBuffer<double> buffer(mask_indices, basis.dimensions_of_basis());
    
    // prepare pointers for OMP parallel region
    
    auto ptr_gto_pair_blocks = &gto_pair_blocks;
    
    auto ptr_aux_gto_blocks = &aux_gto_blocks;
    
    auto ptr_buffer = &buffer;
    
    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();
    
#pragma omp parallel shared(ptr_gto_pair_blocks, ptr_aux_gto_blocks, ptr_buffer)
    {
#pragma omp single nowait
        {
            const auto tasks = omp::make_aux_work_tasks(*ptr_aux_gto_blocks);

            std::ranges::for_each(std::ranges::reverse_view(tasks), [&](const auto& task) {
                auto aux_idx = task[0];
                auto bra_range = std::pair<size_t, size_t>{task[1], task[2]};
                if (const auto nblocks = ptr_gto_pair_blocks->size(); nblocks > 0)
                {
                    for (size_t i = 0; i < nblocks; i++)
                    {
#pragma omp task firstprivate(aux_idx, bra_range, i)
                        {
                            auto aux_gtos = ptr_aux_gto_blocks->at(aux_idx);
                            
                            const auto gp_pairs = ptr_gto_pair_blocks->at(i);
                            
                            CT3CLocalDistributor distributor(ptr_buffer);
                            
                            t3cerifunc::compute(distributor, aux_gtos, gp_pairs, bra_range);
                        }
                    }
                }
            });
        }
    }
    
    
    return buffer;
}
