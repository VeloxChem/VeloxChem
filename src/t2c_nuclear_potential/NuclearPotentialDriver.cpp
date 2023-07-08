#include "NuclearPotentialDriver.hpp"

#include "GtoFunc.hpp"
#include "Matrix.hpp"
#include "MatrixFunc.hpp"
#include "MatrixType.hpp"
#include "NuclearPotentialFunc.hpp"
#include "OpenMPFunc.hpp"

auto
CNuclearPotentialDriver::compute(const CMolecularBasis&       basis,
                                 const CMolecule&             molecule,
                                 const std::vector<double>&   charges,
                                 const std::vector<TPoint3D>& points) const -> CMatrices
{
    CMatrices matrices;

    if (const auto natoms = static_cast<int64_t>(charges.size()); natoms > 0)
    {
        // set up matrices

        for (int64_t i = 0; i < natoms; i++)
        {
            matrices.add(matfunc::makeMatrix(basis, mat_t::symm), i);
        }

        matrices.zero();

        // set up parallelization data

        const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

        const auto work_groups = omp::makeWorkGroup(gto_blocks);

        // prepare pointers for OMP parallel region

        auto ptr_gto_blocks = gto_blocks.data();

        auto ptr_work_groups = work_groups.data();

        auto ptr_matrices = &matrices;

        auto ptr_charges = charges.data();

        auto ptr_points = points.data();

        // execute OMP tasks with static scheduling

        omp::setStaticScheduler();

        const auto ntasks = work_groups.size();

        // compute nuclear potential for each external point

        for (int64_t i = 0; i < natoms; i++)
        {
#pragma omp parallel num_threads(ntasks) shared(ntasks, ptr_gto_blocks, ptr_work_groups, ptr_matrices, ptr_charges, ptr_points) firstprivate(i)
            {
#pragma omp single nowait
                {
                    for (size_t j = 0; j < ntasks; j++)
                    {
#pragma omp task firstprivate(i)
                        {
                            for (const auto& task : ptr_work_groups[j])
                            {
                                auto ptr_matrix = ptr_matrices->getMatrix("X");

                                if (task[0] == task[1])
                                {
                                    const auto gto_block = ptr_gto_blocks[task[0]];

                                    const auto angmom = gto_block.getAngularMomentum();

                                    auto ptr_submatrix = ptr_matrix->getSubMatrix({angmom, angmom});

                                    npotfunc::compute(ptr_submatrix, ptr_charges[i], ptr_points[i], gto_block, angmom, task[2], task[3]);
                                }
                                else
                                {
                                    const auto bra_gto_block = ptr_gto_blocks[task[0]];

                                    const auto ket_gto_block = ptr_gto_blocks[task[1]];

                                    const auto bra_angmom = bra_gto_block.getAngularMomentum();

                                    const auto ket_angmom = ket_gto_block.getAngularMomentum();

                                    auto ptr_submatrix = ptr_matrix->getSubMatrix({bra_angmom, ket_angmom});

                                    const auto ang_order = ptr_matrix->isAngularOrder({bra_angmom, ket_angmom});

                                    npotfunc::compute(ptr_submatrix,
                                                      ptr_charges[i],
                                                      ptr_points[i],
                                                      bra_gto_block,
                                                      ket_gto_block,
                                                      bra_angmom,
                                                      ket_angmom,
                                                      ang_order,
                                                      task[2],
                                                      task[3],
                                                      mat_t::symm);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return matrices;
}
