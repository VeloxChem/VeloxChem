#include "ThreeCenterOverlapDriver.hpp"

#include <algorithm>
#include <iostream>
#include <ranges>
#include <utility>

#include "BasisFunction.hpp"
#include "GtoFunc.hpp"
#include "MatrixFunc.hpp"
#include "OpenMPFunc.hpp"
#include "ThreeCenterOverlapFunc.hpp"
#include "T2CDistributor.hpp"

auto
CThreeCenterOverlapDriver::compute(const std::vector<double>         &exponents,
                                   const std::vector<double>         &factors,
                                   const std::vector<TPoint<double>> &coordinates,
                                   const CMolecularBasis             &basis,
                                   const CMolecule                   &molecule) const -> CMatrix
{
    // set up overlap matrix

    auto ovl_mat = matfunc::make_matrix(basis, mat_t::symmetric);

    ovl_mat.zero();
    
    // prepare external Gaussians data
    
    auto gnorms = exponents;
    
    gnorms.insert(gnorms.end(), factors.cbegin(), factors.cend());
    
    // prepare pointers for OMP parallel region

    auto ptr_basis = &basis;

    auto ptr_molecule = &molecule;

    auto ptr_ovl_mat = &ovl_mat;

    auto ptr_coordinates = &coordinates;

    auto ptr_gnorms = &gnorms;

    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_basis, ptr_molecule, ptr_ovl_mat, ptr_coordinates, ptr_gnorms)
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
                    CT2CDistributor<CMatrix> distributor(ptr_ovl_mat, *ptr_coordinates, *ptr_gnorms);
                    t3ovlfunc::compute(distributor, bra_gtos, ket_gtos, bra_indices, ket_indices, bkequal);
                }
            });
        }
    }

    return ovl_mat;
}
