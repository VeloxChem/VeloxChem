#include "FockDriver.hpp"

#include "BatchFunc.hpp"
#include "ElectronRepulsionFunc.hpp"
#include "GtoBlock.hpp"
#include "GtoPairBlock.hpp"
#include "GtoPairBlockFunc.hpp"
#include "OpenMPFunc.hpp"
#include "T4CMatricesDistributor.hpp"
#include "T4CMatrixDistributor.hpp"
#include "T4COrderedMatrixDistributor.hpp"

auto
CFockDriver::compute(const CMolecularBasis& basis,
                     const CMolecule&       molecule,
                     const CMatrix&         density,
                     const std::string&     label,
                     const double           exchange_factor) const -> CMatrix
{
    // set up Fock matrix

    auto fock_mat = CMatrix(density);

    fock_mat.zero();

    // set up distributor

    auto distributor = CT4CMatrixDistributor(fock_mat, density, label, exchange_factor);

    // set up work groups

    const auto gto_pair_blocks = gtofunc::make_gto_pair_blocks(basis, molecule);

    const auto work_groups = omp::make_work_group(gto_pair_blocks);

    // prepare pointers for OMP parallel region

    auto ptr_gto_pair_blocks = gto_pair_blocks.data();

    auto ptr_work_groups = work_groups.data();

    auto ptr_distributor = &distributor;

    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

    const auto ntasks = work_groups.size();

#pragma omp parallel num_threads(ntasks) shared(ntasks, ptr_gto_pair_blocks, ptr_work_groups, ptr_distributor)
    {
#pragma omp single nowait
        {
            for (size_t i = 0; i < ntasks; i++)
            {
#pragma omp task firstprivate(i)
                {
                    for (const auto& wtask : ptr_work_groups[i])
                    {
                        const std::array<int, 2> bra_indices{wtask[2], wtask[3]};

                        const std::array<int, 2> ket_indices{wtask[4], wtask[5]};

                        if ((wtask[0] == wtask[1]) && (bra_indices == ket_indices))
                        {
                            const auto gto_pair_block = ptr_gto_pair_blocks[wtask[0]];

                            erifunc::compute<CT4CMatrixDistributor>(ptr_distributor, gto_pair_block, bra_indices);
                        }
                        else
                        {
                            const auto bra_gto_pair_block = ptr_gto_pair_blocks[wtask[0]];

                            const auto ket_gto_pair_block = ptr_gto_pair_blocks[wtask[1]];

                            erifunc::compute<CT4CMatrixDistributor>(
                                ptr_distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
                        }
                    }
                }
            }
        }
    }

    fock_mat.symmetrize();

    return fock_mat;
}

auto
CFockDriver::compute(const CT4CScreener& screener,
                     const CMatrix&      density,
                     const std::string&  label,
                     const double        exchange_factor,
                     const int           ithreshold) const -> CMatrix
{
    // set up Fock matrix

    auto fock_mat = CMatrix(density);

    fock_mat.zero();

    // set up distributor

    auto distributor = CT4CMatrixDistributor(fock_mat, density, label, exchange_factor);

    // set up work groups

    const auto gto_pair_blocks = screener.gto_pair_blocks();

    const auto work_tasks = omp::make_work_group(gto_pair_blocks, ithreshold);

    // prepare pointers for OMP parallel region

    auto ptr_gto_pair_blocks = gto_pair_blocks.data();

    auto ptr_work_tasks = work_tasks.data();

    auto ptr_distributor = &distributor;

    // execute OMP tasks with static scheduling

    const auto ntasks = work_tasks.size();

#pragma omp parallel shared(ntasks, ptr_gto_pair_blocks, ptr_work_tasks, ptr_distributor)
    {
#pragma omp single nowait
        {
            for (size_t i = 0; i < ntasks; i++)
            {
#pragma omp task firstprivate(i)
                {
                    const auto wtask = ptr_work_tasks[i];

                    const std::array<int, 2> bra_indices{wtask[4], wtask[5]};

                    const std::array<int, 2> ket_indices{wtask[6], wtask[7]};

                    if ((wtask[0] == wtask[1]) && (wtask[2] == wtask[3]) && (bra_indices == ket_indices))
                    {
                        const auto gto_pair_block = ptr_gto_pair_blocks[wtask[0]].gto_pair_block(wtask[2]);

                        erifunc::compute<CT4CMatrixDistributor>(ptr_distributor, gto_pair_block, bra_indices);
                    }
                    else
                    {
                        const auto bra_gto_pair_block = ptr_gto_pair_blocks[wtask[0]].gto_pair_block(wtask[2]);

                        const auto ket_gto_pair_block = ptr_gto_pair_blocks[wtask[1]].gto_pair_block(wtask[3]);

                        erifunc::compute<CT4CMatrixDistributor>(ptr_distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
                    }
                }
            }
        }
    }

    fock_mat.symmetrize();

    return fock_mat;
}

auto
CFockDriver::compute(const CT4CScreener&             screener,
                     const CMatrices&                densities,
                     const std::vector<std::string>& labels,
                     const std::vector<double>&      exchange_factors,
                     const int                       ithreshold) const -> CMatrices
{
    // set up Fock matrix

    auto fock_mats = CMatrices(densities);

    fock_mats.zero();

    // set up distributor

    auto distributor = CT4CMatricesDistributor(fock_mats, densities, labels, exchange_factors);

    // set up work groups

    const auto gto_pair_blocks = screener.gto_pair_blocks();

    const auto work_tasks = omp::make_work_group(gto_pair_blocks, ithreshold);

    // prepare pointers for OMP parallel region

    auto ptr_gto_pair_blocks = gto_pair_blocks.data();

    auto ptr_work_tasks = work_tasks.data();

    auto ptr_distributor = &distributor;

    // execute OMP tasks with static scheduling

    const auto ntasks = work_tasks.size();

#pragma omp parallel shared(ntasks, ptr_gto_pair_blocks, ptr_work_tasks, ptr_distributor)
    {
#pragma omp single nowait
        {
            for (size_t i = 0; i < ntasks; i++)
            {
#pragma omp task firstprivate(i)
                {
                    const auto wtask = ptr_work_tasks[i];

                    const std::array<int, 2> bra_indices{wtask[4], wtask[5]};

                    const std::array<int, 2> ket_indices{wtask[6], wtask[7]};

                    if ((wtask[0] == wtask[1]) && (wtask[2] == wtask[3]) && (bra_indices == ket_indices))
                    {
                        const auto gto_pair_block = ptr_gto_pair_blocks[wtask[0]].gto_pair_block(wtask[2]);

                        erifunc::compute<CT4CMatricesDistributor>(ptr_distributor, gto_pair_block, bra_indices);
                    }
                    else
                    {
                        const auto bra_gto_pair_block = ptr_gto_pair_blocks[wtask[0]].gto_pair_block(wtask[2]);

                        const auto ket_gto_pair_block = ptr_gto_pair_blocks[wtask[1]].gto_pair_block(wtask[3]);

                        erifunc::compute<CT4CMatricesDistributor>(ptr_distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
                    }
                }
            }
        }
    }

    fock_mats.symmetrize();

    return fock_mats;
}

auto
CFockDriver::ordered_compute(const CT4CScreener& screener,
                             const CMatrix&      density,
                             const std::string&  label,
                             const double        exchange_factor,
                             const int           ithreshold) const -> CMatrix
{
    // set up Fock matrix

    auto fock_mat = CMatrix(density);

    fock_mat.zero();

    // set up work groups

    const auto gto_pair_blocks = screener.gto_pair_blocks();

    const auto work_tasks = omp::make_ordered_work_group(gto_pair_blocks, ithreshold);

    // prepare pointers for OMP parallel region

    auto ptr_gto_pair_blocks = gto_pair_blocks.data();

    auto ptr_work_tasks = work_tasks.data();

    auto ptr_density = density.pointer();

    auto ptr_fock = fock_mat.pointer();

    // execute OMP tasks with static scheduling

    const auto ntasks = work_tasks.size();

#pragma omp parallel shared(ntasks, ptr_gto_pair_blocks, ptr_work_tasks, ptr_fock, ptr_density)
    {
#pragma omp single nowait
        {
            for (size_t i = 0; i < ntasks; i++)
            {
#pragma omp task firstprivate(i)
                {
                    // set up distributor

                    const int max_width = 256;

                    const auto wtask = ptr_work_tasks[i];

                    if ((wtask[0] == wtask[1]) && (wtask[2] == wtask[3]))
                    {
                        const auto gto_pair_block = ptr_gto_pair_blocks[wtask[0]].gto_pair_block(wtask[2]);

                        auto distributor = CT4COrderedMatrixDistributor(ptr_fock, ptr_density, gto_pair_block, label, exchange_factor);

                        const auto bra_size = gto_pair_block.number_of_contracted_pairs();

                        const auto bra_tasks = batch::reduced_number_of_batches(bra_size, max_width);

                        for (int k = 0; k < bra_tasks; k++)
                        {
                            const auto bra_first = batch::batch_index(k, bra_size, bra_tasks);

                            const auto bra_last = batch::batch_index(k + 1, bra_size, bra_tasks);

                            for (int l = k; l < bra_tasks; l++)
                            {
                                if (l == k)
                                {
                                    erifunc::compute<CT4COrderedMatrixDistributor>(&distributor, gto_pair_block, {bra_first, bra_last});

                                    // std::cout << "Diagonal term: " << bra_first << " , " << bra_last << " for " <<
                                    // gto_pair_block.angular_momentums()[0] << " , " <<  gto_pair_block.angular_momentums()[1] << " Index : " <<
                                    // wtask[0] << " , " << wtask[1] << " , " << wtask[2] << " , " << wtask[3]  << std::endl;
                                }
                                else
                                {
                                    const auto ket_first = batch::batch_index(l, bra_size, bra_tasks);

                                    const auto ket_last = batch::batch_index(l + 1, bra_size, bra_tasks);

                                    // std::cout << "Diagonal off-term : " << bra_first << " , " << bra_last << " , " << ket_first << " , " <<
                                    // ket_last <<  " for " << gto_pair_block.angular_momentums()[0] << " , " << gto_pair_block.angular_momentums()[1]
                                    // << " , "  << gto_pair_block.angular_momentums()[0] << " , " <<  gto_pair_block.angular_momentums()[1] << "
                                    // Index : " << wtask[0] << " , " << wtask[1] << " , " << wtask[2] << " , " << wtask[3]  << std::endl;

                                    erifunc::compute<CT4COrderedMatrixDistributor>(
                                        &distributor, gto_pair_block, gto_pair_block, {bra_first, bra_last}, {ket_first, ket_last});
                                }
                            }
                        }

                        distributor.accumulate(gto_pair_block);
                    }
                    else
                    {
                        const auto bra_gto_pair_block = ptr_gto_pair_blocks[wtask[0]].gto_pair_block(wtask[2]);

                        const auto ket_gto_pair_block = ptr_gto_pair_blocks[wtask[1]].gto_pair_block(wtask[3]);

                        auto distributor =
                            CT4COrderedMatrixDistributor(ptr_fock, ptr_density, bra_gto_pair_block, ket_gto_pair_block, label, exchange_factor);

                        const auto bra_size = bra_gto_pair_block.number_of_contracted_pairs();

                        const auto ket_size = ket_gto_pair_block.number_of_contracted_pairs();

                        const auto bra_tasks = batch::reduced_number_of_batches(bra_size, max_width);

                        const auto ket_tasks = batch::reduced_number_of_batches(ket_size, max_width);

                        for (int k = 0; k < bra_tasks; k++)
                        {
                            const auto bra_first = batch::batch_index(k, bra_size, bra_tasks);

                            const auto bra_last = batch::batch_index(k + 1, bra_size, bra_tasks);

                            for (int l = 0; l < ket_tasks; l++)
                            {
                                const auto ket_first = batch::batch_index(l, ket_size, ket_tasks);

                                const auto ket_last = batch::batch_index(l + 1, ket_size, ket_tasks);

                                // std::cout << "Off-Diagonal term : " << bra_first << " , " << bra_last << " , " << ket_first << " , " << ket_last <<
                                // " for " << bra_gto_pair_block.angular_momentums()[0] << " , " <<  bra_gto_pair_block.angular_momentums()[1] << " ,
                                // "  << ket_gto_pair_block.angular_momentums()[0] << " , " <<  ket_gto_pair_block.angular_momentums()[1] << " Index :
                                // " << wtask[0] << " , " << wtask[1]  << " , " << wtask[2] << " , " << wtask[3] <<  std::endl;

                                erifunc::compute<CT4COrderedMatrixDistributor>(
                                    &distributor, bra_gto_pair_block, ket_gto_pair_block, {bra_first, bra_last}, {ket_first, ket_last});
                            }
                        }

                        distributor.accumulate(bra_gto_pair_block, ket_gto_pair_block);
                    }
                }
            }
        }
    }

    fock_mat.symmetrize();

    return fock_mat;
}
