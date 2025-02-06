#include "ThreeCenterElectronRepulsionDriver.hpp"

#include <algorithm>

#include "GtoFunc.hpp"
#include "GtoPairBlockFunc.hpp"
#include "OpenMPFunc.hpp"
#include "ThreeCenterElectronRepulsionFunc.hpp"
#include "T3CDistributor.hpp"

#include <iostream>

auto
CThreeCenterElectronRepulsionDriver::compute(const CMolecularBasis &basis, const CMolecularBasis &aux_basis, const CMolecule &molecule) const -> CT3FlatBuffer<double>
{
    std::vector<size_t> aux_indices(aux_basis.dimensions_of_basis());
    
    std::iota(aux_indices.begin(), aux_indices.end(), size_t{0});
    
    CT3FlatBuffer<double> buffer(aux_indices, basis.dimensions_of_basis());
    
    // set up basis function pairs blocks

    const auto gto_pair_blocks = gtofunc::make_gto_pair_blocks(basis, molecule);
    
    const auto aux_gto_blocks = gtofunc::make_gto_blocks(aux_basis, molecule);
    
    // prepare pointers for OMP parallel region

    auto ptr_aux_basis = &aux_basis;

    auto ptr_molecule = &molecule;
    
    auto ptr_gto_pair_blocks = &gto_pair_blocks;
    
    auto ptr_aux_gto_blocks = aux_gto_blocks.data();
    
    auto ptr_buffer = &buffer;
    
    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_aux_basis, ptr_molecule, ptr_gto_pair_blocks, ptr_aux_gto_blocks, ptr_buffer)
    {
#pragma omp single nowait
        {
            const auto tasks = omp::make_aux_work_tasks(aux_gto_blocks);

            std::ranges::for_each(std::ranges::reverse_view(tasks), [&](const auto& task) {
                ///std::cout << "Task : " << task[0] << " , " << task[1] << " , " << task[2] << std::endl;
                auto aux_gtos    = ptr_aux_gto_blocks[task[0]];
                auto aux_indices = std::pair<size_t, size_t>{task[1], task[2]};
#pragma omp task firstprivate(aux_gtos, aux_indices) shared(ptr_gto_pair_blocks)
                {
                    if (const auto nblocks = ptr_gto_pair_blocks->size(); nblocks > 0)
                    {
                        for (size_t i = 0; i < nblocks; i++)
                        {
                            const auto gp_pairs = ptr_gto_pair_blocks->at(i);
                            
                            CT3CDistributor distributor(ptr_buffer);
                            
                            t3cerifunc::compute(distributor, aux_gtos, gp_pairs, aux_indices);
                        }
                    }
                }
            });
        }
    }

    return buffer;
}
