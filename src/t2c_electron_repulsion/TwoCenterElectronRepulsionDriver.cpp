#include "TwoCenterElectronRepulsionDriver.hpp"

#include <algorithm>
#include <iostream>
#include <ranges>
#include <utility>

#include "GtoFunc.hpp"
#include "TwoCenterElectronRepulsionFunc.hpp"
#include "MatrixFunc.hpp"
#include "OpenMPFunc.hpp"
#include "T2CDistributor.hpp"

auto
CTwoCenterElectronRepulsionDriver::compute(const CMolecularBasis& basis, const CMolecule& molecule) const -> CMatrix
{
    // set up electron repulsion matrix

    auto eri_mat = matfunc::make_matrix(basis, mat_t::symmetric);

    eri_mat.zero();

    // prepare pointers for OMP parallel region

    auto ptr_basis = &basis;

    auto ptr_molecule = &molecule;

    auto ptr_eri_mat = &eri_mat;

    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_basis, ptr_molecule, ptr_eri_mat)
    {
#pragma omp single nowait
        {
            const auto gto_blocks = gtofunc::make_gto_blocks(*ptr_basis, *ptr_molecule);

            const auto tasks = omp::make_work_tasks(gto_blocks);

            std::ranges::for_each(std::ranges::reverse_view(tasks), [&](const auto& task) {
                auto bra_gtos    = gto_blocks[task[0]];
                auto ket_gtos    = gto_blocks[task[1]];
                auto bra_indices = std::pair<size_t, size_t>{task[2], task[3]};
                auto ket_indices = std::pair<size_t, size_t>{task[4], task[5]};
                bool bkequal     = (task[0] == task[1]) && (task[2] == task[4]) && (task[3] == task[5]);
#pragma omp task firstprivate(bra_gtos, ket_gtos, bra_indices, ket_indices, bkequal)
                {
                    CT2CDistributor<CMatrix> distributor(ptr_eri_mat);
                    t2cerifunc::compute(distributor, bra_gtos, ket_gtos, bra_indices, ket_indices, bkequal);
                }
            });
        }
    }

    return eri_mat;
}
