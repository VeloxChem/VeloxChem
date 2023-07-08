#include "QuadrupoleDriver.hpp"

#include "GtoFunc.hpp"
#include "Matrix.hpp"
#include "MatrixFunc.hpp"
#include "MatrixType.hpp"
#include "OpenMPFunc.hpp"
#include "QuadrupoleFunc.hpp"

auto
CQuadrupoleDriver::compute(const CMolecularBasis& basis, const CMolecule& molecule, const TPoint3D& point) const -> CMatrices
{
    CMatrices quad_matrix;

    quad_matrix.add(matfunc::makeMatrix(basis, mat_t::symm), "xx");

    quad_matrix.add(matfunc::makeMatrix(basis, mat_t::symm), "xy");

    quad_matrix.add(matfunc::makeMatrix(basis, mat_t::symm), "xz");

    quad_matrix.add(matfunc::makeMatrix(basis, mat_t::symm), "yy");

    quad_matrix.add(matfunc::makeMatrix(basis, mat_t::symm), "yz");

    quad_matrix.add(matfunc::makeMatrix(basis, mat_t::symm), "zz");

    quad_matrix.zero();

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto work_groups = omp::makeWorkGroup(gto_blocks);

    // prepare pointers for OMP parallel region

    auto ptr_gto_blocks = gto_blocks.data();

    auto ptr_work_groups = work_groups.data();

    auto ptr_quad_matrix = &quad_matrix;

    // execute OMP tasks with static scheduling

    omp::setStaticScheduler();

    const auto ntasks = work_groups.size();

#pragma omp parallel num_threads(ntasks) shared(ntasks, ptr_gto_blocks, ptr_work_groups, ptr_quad_matrix, point)
    {
#pragma omp single nowait
        {
            for (size_t i = 0; i < ntasks; i++)
            {
#pragma omp task firstprivate(i)
                {
                    for (const auto& task : ptr_work_groups[i])
                    {
                        auto ptr_matrix_xx = ptr_quad_matrix->getMatrix("XX");

                        auto ptr_matrix_xy = ptr_quad_matrix->getMatrix("XY");

                        auto ptr_matrix_xz = ptr_quad_matrix->getMatrix("XZ");

                        auto ptr_matrix_yy = ptr_quad_matrix->getMatrix("YY");

                        auto ptr_matrix_yz = ptr_quad_matrix->getMatrix("YZ");

                        auto ptr_matrix_zz = ptr_quad_matrix->getMatrix("ZZ");

                        if (task[0] == task[1])
                        {
                            const auto gto_block = ptr_gto_blocks[task[0]];

                            const auto angmom = gto_block.getAngularMomentum();

                            auto ptr_submatrix_xx = ptr_matrix_xx->getSubMatrix({angmom, angmom});

                            auto ptr_submatrix_xy = ptr_matrix_xy->getSubMatrix({angmom, angmom});

                            auto ptr_submatrix_xz = ptr_matrix_xz->getSubMatrix({angmom, angmom});

                            auto ptr_submatrix_yy = ptr_matrix_yy->getSubMatrix({angmom, angmom});

                            auto ptr_submatrix_yz = ptr_matrix_yz->getSubMatrix({angmom, angmom});

                            auto ptr_submatrix_zz = ptr_matrix_zz->getSubMatrix({angmom, angmom});

                            quadfunc::compute(ptr_submatrix_xx,
                                              ptr_submatrix_xy,
                                              ptr_submatrix_xz,
                                              ptr_submatrix_yy,
                                              ptr_submatrix_yz,
                                              ptr_submatrix_zz,
                                              point,
                                              gto_block,
                                              angmom,
                                              task[2],
                                              task[3]);
                        }
                        else
                        {
                            const auto bra_gto_block = ptr_gto_blocks[task[0]];

                            const auto ket_gto_block = ptr_gto_blocks[task[1]];

                            const auto bra_angmom = bra_gto_block.getAngularMomentum();

                            const auto ket_angmom = ket_gto_block.getAngularMomentum();

                            auto ptr_submatrix_xx = ptr_matrix_xx->getSubMatrix({bra_angmom, ket_angmom});

                            auto ptr_submatrix_xy = ptr_matrix_xy->getSubMatrix({bra_angmom, ket_angmom});

                            auto ptr_submatrix_xz = ptr_matrix_xz->getSubMatrix({bra_angmom, ket_angmom});

                            auto ptr_submatrix_yy = ptr_matrix_yy->getSubMatrix({bra_angmom, ket_angmom});

                            auto ptr_submatrix_yz = ptr_matrix_yz->getSubMatrix({bra_angmom, ket_angmom});

                            auto ptr_submatrix_zz = ptr_matrix_zz->getSubMatrix({bra_angmom, ket_angmom});

                            const auto ang_order = ptr_matrix_xx->isAngularOrder({bra_angmom, ket_angmom});

                            quadfunc::compute(ptr_submatrix_xx,
                                              ptr_submatrix_xy,
                                              ptr_submatrix_xz,
                                              ptr_submatrix_yy,
                                              ptr_submatrix_yz,
                                              ptr_submatrix_zz,
                                              point,
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

    return quad_matrix;
}
