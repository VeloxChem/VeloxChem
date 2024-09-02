#include "T4CScreener.hpp"

#include "DiagonalElectronRepulsionFunc.hpp"
#include "GtoPairBlockFunc.hpp"
#include "OpenMPFunc.hpp"
#include "T4CDiagonalDistributor.hpp"

auto
CT4CScreener::partition(const CMolecularBasis& basis, const CMolecule& molecule, const std::string& label) -> void
{
    const auto gto_pair_blocks = gtofunc::make_gto_pair_blocks(basis, molecule);

    if (const auto nblocks = gto_pair_blocks.size(); nblocks > 0)
    {
        // set up max. integral values

        std::vector<std::vector<double>> max_values;

        for (size_t i = 0; i < nblocks; i++)
        {
            max_values.push_back(std::vector<double>(gto_pair_blocks[i].number_of_contracted_pairs(), 0.0));
        }

        // set up work groups

        const auto work_groups = omp::make_diag_work_group(gto_pair_blocks);

        // prepare pointers for OMP parallel region

        auto ptr_gto_pair_blocks = gto_pair_blocks.data();

        auto ptr_work_groups = work_groups.data();

        auto ptr_max_values = max_values.data();

        // execute OMP tasks with static scheduling

        omp::set_static_scheduler();

        const auto ntasks = work_groups.size();

#pragma omp parallel num_threads(ntasks) shared(ntasks, ptr_gto_pair_blocks, ptr_work_groups, ptr_max_values)
        {
#pragma omp single nowait
            {
                for (size_t i = 0; i < ntasks; i++)
                {
#pragma omp task firstprivate(i)
                    {
                        for (const auto& wtask : ptr_work_groups[i])
                        {
                            const auto gto_pair_block = ptr_gto_pair_blocks[wtask[0]];

                            CT4CDiagonalDistributor distributor(ptr_max_values[wtask[0]].data());

                            if (label == "eri")
                            {
                                erifunc::diag_compute<CT4CDiagonalDistributor>(&distributor, gto_pair_block, {wtask[1], wtask[2]});
                            }
                        }
                    }
                }
            }
        }

        //        for (size_t i = 0; i < max_values.size(); i++)
        //        {
        //            std::cout << "*** GTO pair block : " << gto_pair_blocks[i].angular_momentums()[0] << gto_pair_blocks[i].angular_momentums()[1]
        //            << std::endl;
        //
        //            for (size_t j = 0; j < max_values[i].size(); j++)
        //            {
        //                std::cout << max_values[i][j] << " ";
        //            }
        //
        //            std::cout << std::endl;
        //        }

        // apply Cauchy–Schwarz partitioning

        for (size_t i = 0; i < nblocks; i++)
        {
            _gto_pair_blocks.push_back(CBlockedGtoPairBlock(gto_pair_blocks[i], max_values[i]));
        }
    }
}

auto
CT4CScreener::gto_pair_blocks() const -> std::vector<CBlockedGtoPairBlock>
{
    return _gto_pair_blocks;
}
