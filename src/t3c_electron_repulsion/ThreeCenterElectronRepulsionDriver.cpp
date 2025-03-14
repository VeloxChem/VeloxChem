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
