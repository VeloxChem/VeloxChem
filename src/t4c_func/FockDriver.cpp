#include "FockDriver.hpp"

#include <algorithm>
#include <ranges>

#include "CustomViews.hpp"
#include "ElectronRepulsionFunc.hpp"
#include "GtoBlock.hpp"
#include "GtoPairBlock.hpp"
#include "GtoPairBlockFunc.hpp"
#include "OpenMPFunc.hpp"
#include "T4CMatricesDistributor.hpp"
#include "T4CMatrixDistributor.hpp"

auto
CFockDriver::compute(const CMolecularBasis& basis,
                     const CMolecule&       molecule,
                     const CMatrix&         density,
                     const std::string&     label,
                     const double           exchange_factor,
                     const double           omega) const -> CMatrix
{
    // set up Fock matrix

    auto fock_mat = CMatrix(density);

    fock_mat.zero();

    // set up basis function pairs blocks

    const auto gto_pair_blocks = gtofunc::make_gto_pair_blocks(basis, molecule);

    // prepare pointers for OMP parallel region

    auto ptr_gto_pair_blocks = &gto_pair_blocks;

    auto ptr_density = &density;

    auto ptr_fock = &fock_mat;

    // execute OMP tasks with static scheduling

#pragma omp parallel shared(ptr_gto_pair_blocks, ptr_density, ptr_fock, label, exchange_factor, omega)
    {
#pragma omp single nowait
        {
            const auto nblocks = ptr_gto_pair_blocks->size();

            auto ptr_gto_pairs_data = ptr_gto_pair_blocks->data();

            std::ranges::for_each(views::triangular(nblocks) | std::views::reverse, [&](const auto& index) {
                const size_t i = index.first;
                const size_t j = index.second;
#pragma omp task firstprivate(i, j)
                {
                    auto                  bra_gpairs = ptr_gto_pairs_data[i];
                    auto                  ket_gpairs = ptr_gto_pairs_data[j];
                    CT4CMatrixDistributor distributor(ptr_fock, ptr_density, label, exchange_factor, omega);
                    distributor.set_indices(bra_gpairs, ket_gpairs);
                    auto bra_range = std::pair<size_t, size_t>(0, bra_gpairs.number_of_contracted_pairs());
                    auto ket_range = std::pair<size_t, size_t>(0, ket_gpairs.number_of_contracted_pairs());
                    erifunc::compute<CT4CMatrixDistributor>(distributor, bra_gpairs, ket_gpairs, bra_range, ket_range, i == j);
                    distributor.accumulate(bra_gpairs, ket_gpairs);
                }
            });
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
                     const double        omega,
                     const int           ithreshold) const -> CMatrix
{
    auto bsfac = _determine_block_size_factor(_get_nao(density));

    // set up Fock matrix

    auto fock_mat = CMatrix(density);

    fock_mat.zero();

    // prepare pointers for OMP parallel region

    auto ptr_screener = &screener;

    auto ptr_density = &density;

    auto ptr_fock = &fock_mat;

    // execute OMP tasks with static scheduling

#pragma omp parallel shared(ptr_screener, ptr_density, ptr_fock, label, exchange_factor, omega, ithreshold)
    {
#pragma omp single nowait
        {
            auto gto_pair_blocks = ptr_screener->gto_pair_blocks();

            const auto work_tasks = omp::make_work_group(gto_pair_blocks, ithreshold, bsfac);

            std::ranges::for_each(std::views::reverse(work_tasks), [&](const auto& task) {
                const auto bra_gpairs = gto_pair_blocks[task[0]].gto_pair_block(static_cast<int>(task[2]));
                const auto ket_gpairs = gto_pair_blocks[task[1]].gto_pair_block(static_cast<int>(task[3]));
                const auto bra_range  = std::pair<size_t, size_t>{task[4], task[5]};
                const auto ket_range  = std::pair<size_t, size_t>{task[6], task[7]};
                const bool diagonal   = (task[0] == task[1]) && (task[2] == task[3]) && (bra_range == ket_range);
#pragma omp task firstprivate(bra_gpairs, ket_gpairs, bra_range, ket_range, diagonal)
                {
                    CT4CMatrixDistributor distributor(ptr_fock, ptr_density, label, exchange_factor, omega);
                    distributor.set_indices(bra_gpairs, ket_gpairs);
                    erifunc::compute<CT4CMatrixDistributor>(distributor, bra_gpairs, ket_gpairs, bra_range, ket_range, diagonal);
                    distributor.accumulate(bra_gpairs, ket_gpairs);
                }
            });
        }
    }

    fock_mat.symmetrize();

    return fock_mat;
}

auto
CFockDriver::compute(const CT4CScreener&             screener,
                     const CMatrices&                densities,
                     const std::vector<std::string>& labels,
                     const double                    exchange_factor,
                     const double                    omega,
                     const int                       ithreshold) const -> CMatrices
{
    auto bsfac = _determine_block_size_factor(_get_nao(densities));

    // set up Fock matrices

    auto fock_mats = CMatrices(densities);

    fock_mats.zero();

    // prepare pointers for OMP parallel region

    auto ptr_screener = &screener;

    auto ptr_densities = &densities;

    auto ptr_focks = &fock_mats;

    auto ptr_labels = &labels;

    // execute OMP tasks with static scheduling

#pragma omp parallel shared(ptr_screener, ptr_densities, ptr_focks, ptr_labels, exchange_factor, omega, ithreshold)
    {
#pragma omp single nowait
        {
            auto gto_pair_blocks = ptr_screener->gto_pair_blocks();

            const auto work_tasks = omp::make_work_group(gto_pair_blocks, ithreshold, bsfac);

            std::ranges::for_each(std::views::reverse(work_tasks), [&](const auto& task) {
                const auto bra_gpairs = gto_pair_blocks[task[0]].gto_pair_block(static_cast<int>(task[2]));
                const auto ket_gpairs = gto_pair_blocks[task[1]].gto_pair_block(static_cast<int>(task[3]));
                const auto bra_range  = std::pair<size_t, size_t>{task[4], task[5]};
                const auto ket_range  = std::pair<size_t, size_t>{task[6], task[7]};
                const bool diagonal   = (task[0] == task[1]) && (task[2] == task[3]) && (bra_range == ket_range);
#pragma omp task firstprivate(bra_gpairs, ket_gpairs, bra_range, ket_range, diagonal)
                {
                    CT4CMatricesDistributor distributor(ptr_focks, ptr_densities, *ptr_labels, exchange_factor, omega);
                    distributor.set_indices(bra_gpairs, ket_gpairs);
                    erifunc::compute<CT4CMatricesDistributor>(distributor, bra_gpairs, ket_gpairs, bra_range, ket_range, diagonal);
                    distributor.accumulate(bra_gpairs, ket_gpairs);
                }
            });
        }
    }

    fock_mats.symmetrize();

    return fock_mats;
}

auto
CFockDriver::compute(const CT4CScreener& screener,
                     const int           rank,
                     const int           nodes,
                     const CMatrix&      density,
                     const std::string&  label,
                     const double        exchange_factor,
                     const double        omega,
                     const int           ithreshold) const -> CMatrix
{
    auto bsfac = _determine_block_size_factor(_get_nao(density));

    // set up Fock matrix

    auto fock_mat = CMatrix(density);

    fock_mat.zero();

    // prepare pointers for OMP parallel region

    auto ptr_screener = &screener;

    auto ptr_density = &density;

    auto ptr_fock = &fock_mat;

    // execute OMP tasks with static scheduling

#pragma omp parallel shared(ptr_screener, rank, nodes, ptr_density, ptr_fock, label, exchange_factor, omega, ithreshold)
    {
#pragma omp single nowait
        {
            auto gto_pair_blocks = ptr_screener->gto_pair_blocks();

            const auto work_tasks = omp::make_work_group(gto_pair_blocks, ithreshold, bsfac);

            const auto red_tasks = omp::partition_tasks(work_tasks, rank, nodes);

            std::ranges::for_each(std::views::reverse(red_tasks), [&](const auto& task) {
                const auto bra_gpairs = gto_pair_blocks[task[0]].gto_pair_block(static_cast<int>(task[2]));
                const auto ket_gpairs = gto_pair_blocks[task[1]].gto_pair_block(static_cast<int>(task[3]));
                const auto bra_range  = std::pair<size_t, size_t>{task[4], task[5]};
                const auto ket_range  = std::pair<size_t, size_t>{task[6], task[7]};
                const bool diagonal   = (task[0] == task[1]) && (task[2] == task[3]) && (bra_range == ket_range);
#pragma omp task firstprivate(bra_gpairs, ket_gpairs, bra_range, ket_range, diagonal)
                {
                    CT4CMatrixDistributor distributor(ptr_fock, ptr_density, label, exchange_factor, omega);
                    distributor.set_indices(bra_gpairs, ket_gpairs);
                    erifunc::compute<CT4CMatrixDistributor>(distributor, bra_gpairs, ket_gpairs, bra_range, ket_range, diagonal);
                    distributor.accumulate(bra_gpairs, ket_gpairs);
                }
            });
        }
    }

    fock_mat.symmetrize();

    return fock_mat;
}

auto
CFockDriver::compute(const CT4CScreener&             screener,
                     const int                       rank,
                     const int                       nodes,
                     const CMatrices&                densities,
                     const std::vector<std::string>& labels,
                     const double                    exchange_factor,
                     const double                    omega,
                     const int                       ithreshold) const -> CMatrices
{
    auto bsfac = _determine_block_size_factor(_get_nao(densities));

    // set up Fock matrices

    auto fock_mats = CMatrices(densities);

    fock_mats.zero();

    // prepare pointers for OMP parallel region

    auto ptr_screener = &screener;

    auto ptr_densities = &densities;

    auto ptr_focks = &fock_mats;

    auto ptr_labels = &labels;

    // execute OMP tasks with static scheduling

#pragma omp parallel shared(ptr_screener, ptr_densities, ptr_focks, ptr_labels, exchange_factor, omega, ithreshold)
    {
#pragma omp single nowait
        {
            auto gto_pair_blocks = ptr_screener->gto_pair_blocks();

            const auto work_tasks = omp::make_work_group(gto_pair_blocks, ithreshold, bsfac);

            const auto red_tasks = omp::partition_tasks(work_tasks, rank, nodes);

            std::ranges::for_each(std::views::reverse(red_tasks), [&](const auto& task) {
                const auto bra_gpairs = gto_pair_blocks[task[0]].gto_pair_block(static_cast<int>(task[2]));
                const auto ket_gpairs = gto_pair_blocks[task[1]].gto_pair_block(static_cast<int>(task[3]));
                const auto bra_range  = std::pair<size_t, size_t>{task[4], task[5]};
                const auto ket_range  = std::pair<size_t, size_t>{task[6], task[7]};
                const bool diagonal   = (task[0] == task[1]) && (task[2] == task[3]) && (bra_range == ket_range);
#pragma omp task firstprivate(bra_gpairs, ket_gpairs, bra_range, ket_range, diagonal)
                {
                    CT4CMatricesDistributor distributor(ptr_focks, ptr_densities, *ptr_labels, exchange_factor, omega);
                    distributor.set_indices(bra_gpairs, ket_gpairs);
                    erifunc::compute<CT4CMatricesDistributor>(distributor, bra_gpairs, ket_gpairs, bra_range, ket_range, diagonal);
                    distributor.accumulate(bra_gpairs, ket_gpairs);
                }
            });
        }
    }

    fock_mats.symmetrize();

    return fock_mats;
}

auto
CFockDriver::set_block_size_factor(const int factor) -> void
{
    _block_size_factor = factor;
}

auto
CFockDriver::_get_nao(const CMatrix& mat) const -> int
{
    return mat.number_of_rows();
}

auto
CFockDriver::_get_nao(const CMatrices& mats) const -> int
{
    auto keys = mats.keys();

    auto mat_ptr = mats.matrix(keys[0]);

    return mat_ptr->number_of_rows();
}

auto
CFockDriver::_determine_block_size_factor(const int nao) const -> int
{
    if (nao < 300)
    {
        return 32 * _block_size_factor;
    }
    else if (nao < 1500)
    {
        return 16 * _block_size_factor;
    }
    else if (nao < 4500)
    {
        return 8 * _block_size_factor;
    }
    else if (nao < 6000)
    {
        return 4 * _block_size_factor;
    }
    else if (nao < 10000)
    {
        return 2 * _block_size_factor;
    }
    else
    {
        return _block_size_factor;
    }
}
