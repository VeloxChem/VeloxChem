#include "OverlapDriver.hpp"

#include "GtoFunc.hpp"
#include "MatrixFunc.hpp"
#include "OpenMPFunc.hpp"
#include "OverlapFunc.hpp"

auto
COverlapDriver::compute(const CMolecularBasis& basis, const CMolecule& molecule) const -> CMatrix
{
    auto ovl_matrix = matfunc::makeMatrix(basis, mat_t::symm);

    ovl_matrix.zero();

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto work_groups = omp::makeWorkGroup(gto_blocks);

    // prepare pointers for OMP parallel region

    auto ptr_gto_blocks = gto_blocks.data();

    auto ptr_work_groups = work_groups.data();

    auto ptr_ovl_matrix = &ovl_matrix;

    // execute OMP tasks with static scheduling

    omp::setStaticScheduler();

    const auto ntasks = work_groups.size();

#pragma omp parallel num_threads(ntasks) shared(ntasks, ptr_gto_blocks, ptr_work_groups, ptr_ovl_matrix)
    {
#pragma omp single nowait
        {
            for (size_t i = 0; i < ntasks; i++)
            {
#pragma omp task firstprivate(i)
                {
                    for (const auto& task : ptr_work_groups[i])
                    {
                        const auto mat_type = ptr_ovl_matrix->getType();

                        if (task[0] == task[1])
                        {
                            const auto gto_block = ptr_gto_blocks[task[0]];

                            const auto angmom = gto_block.getAngularMomentum();

                            auto ptr_submatrix = ptr_ovl_matrix->getSubMatrix({angmom, angmom});

                            ovlfunc::compute(ptr_submatrix, gto_block, angmom, task[2], task[3]);
                        }
                        else
                        {
                            const auto bra_gto_block = ptr_gto_blocks[task[0]];

                            const auto ket_gto_block = ptr_gto_blocks[task[1]];

                            const auto bra_angmom = bra_gto_block.getAngularMomentum();

                            const auto ket_angmom = ket_gto_block.getAngularMomentum();

                            auto ptr_submatrix = ptr_ovl_matrix->getSubMatrix({bra_angmom, ket_angmom});

                            const auto ang_order = ptr_ovl_matrix->isAngularOrder({bra_angmom, ket_angmom});

                            ovlfunc::compute(
                                ptr_submatrix, bra_gto_block, ket_gto_block, bra_angmom, ket_angmom, ang_order, task[2], task[3], mat_type);
                        }
                    }
                }
            }
        }
    }

    return ovl_matrix;
}
