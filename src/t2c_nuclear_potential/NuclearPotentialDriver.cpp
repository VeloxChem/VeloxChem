#include "NuclearPotentialDriver.hpp"

#include <algorithm>
#include <iostream>
#include <ranges>
#include <utility>

#include "GtoFunc.hpp"
#include "MatrixFunc.hpp"
#include "NuclearPotentialFunc.hpp"
#include "OpenMPFunc.hpp"
#include "DenseMatrixDistributor.hpp"
#include "T2CDistributor.hpp"

auto
CNuclearPotentialDriver::compute(const CMolecularBasis& basis, const CMolecule& molecule) const -> CMatrix
{
    auto charges = molecule.charges();

    auto coordinates = molecule.coordinates("au");

    return compute(charges, coordinates, basis, molecule);
}

auto
CNuclearPotentialDriver::compute(const std::vector<double>&         charges,
                                 const std::vector<TPoint<double>>& coordinates,
                                 const CMolecularBasis&             basis,
                                 const CMolecule&                   molecule) const -> CMatrix
{
    // set up nuclear potential matrix

    auto npot_mat = matfunc::make_matrix(basis, mat_t::symmetric);

    npot_mat.zero();

    // prepare pointers for OMP parallel region

    auto ptr_charges = &charges;

    auto ptr_coordinates = &coordinates;

    auto ptr_basis = &basis;

    auto ptr_molecule = &molecule;

    auto ptr_npot_mat = &npot_mat;

    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_charges, ptr_coordinates, ptr_basis, ptr_molecule, ptr_npot_mat)
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
                    CT2CDistributor<CMatrix> distributor(ptr_npot_mat, *ptr_coordinates, *ptr_charges);
                    npotfunc::compute(distributor, bra_gtos, ket_gtos, bra_indices, ket_indices, bkequal);
                }
            });
        }
    }

    return npot_mat;
}

auto
CNuclearPotentialDriver::compute(CDenseMatrix&                 gmatrix,
                                 const std::vector<CGtoBlock>& gto_blocks,
                                 const CDenseMatrix&           fmatrix,
                                 const size_t                  gindex,
                                 const size_t                  naos,
                                 const double                  gpoint_x,
                                 const double                  gpoint_y,
                                 const double                  gpoint_z,
                                 const double                  gpoint_w) const -> void
{
    std::vector<double> charges({1.0, });
    
    std::vector<TPoint<double>> coordinates({TPoint<double>({gpoint_x, gpoint_y, gpoint_z}),});
    
    if (const auto nblocks = gto_blocks.size(); nblocks > 0)
    {
        for (size_t i = 0; i < nblocks; i++)
        {
            auto bra_indices = std::pair<size_t, size_t>{size_t{0}, gto_blocks[i].number_of_basis_functions()};
            
            CDenseMatrixDistributor distributor_ii(&gmatrix, coordinates, charges, &fmatrix, gpoint_w);
            
            npotfunc::compute(distributor_ii, gto_blocks[i], gto_blocks[i], bra_indices, bra_indices, true);
            
            for (size_t j = i + 1; j < nblocks; j++)
            {
                auto ket_indices = std::pair<size_t, size_t>{size_t{0}, gto_blocks[j].number_of_basis_functions()};
                
                CDenseMatrixDistributor distributor_ij(&gmatrix, coordinates, charges, &fmatrix, gpoint_w);
                
                npotfunc::compute(distributor_ij, gto_blocks[i], gto_blocks[j], bra_indices, ket_indices, false);
            }
        }
    }
}
