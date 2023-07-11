#include "OctupoleDriver.hpp"

#include "GtoFunc.hpp"
#include "Matrix.hpp"
#include "MatrixFunc.hpp"
#include "MatrixType.hpp"
#include "OctupoleFunc.hpp"
#include "OpenMPFunc.hpp"
#include "MatricesFunc.hpp"

auto
COctupoleDriver::compute(const CMolecularBasis& basis, const CMolecule& molecule, const TPoint3D& point) const -> CMatrices
{
    // set up octupoles matrix
    
    auto octu_matrix = matfunc::makeMatrices(3, basis, mat_t::symm);

    octu_matrix.zero();

    // set up work groups 
    
    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto work_groups = omp::makeWorkGroup(gto_blocks);

    // prepare pointers for OMP parallel region

    auto ptr_gto_blocks = gto_blocks.data();

    auto ptr_work_groups = work_groups.data();

    auto ptr_octu_matrix = &octu_matrix;

    // execute OMP tasks with static scheduling

    omp::setStaticScheduler();

    const auto ntasks = work_groups.size();

#pragma omp parallel num_threads(ntasks) shared(ntasks, ptr_gto_blocks, ptr_work_groups, ptr_octu_matrix, point)
    {
#pragma omp single nowait
        {
            for (size_t i = 0; i < ntasks; i++)
            {
#pragma omp task firstprivate(i)
                {
                    for (const auto& task : ptr_work_groups[i])
                    {
                        auto ptr_matrix_xxx = ptr_octu_matrix->getMatrix("XXX");

                        auto ptr_matrix_xxy = ptr_octu_matrix->getMatrix("XXY");

                        auto ptr_matrix_xxz = ptr_octu_matrix->getMatrix("XXZ");

                        auto ptr_matrix_xyy = ptr_octu_matrix->getMatrix("XYY");

                        auto ptr_matrix_xyz = ptr_octu_matrix->getMatrix("XYZ");

                        auto ptr_matrix_xzz = ptr_octu_matrix->getMatrix("XZZ");

                        auto ptr_matrix_yyy = ptr_octu_matrix->getMatrix("YYY");

                        auto ptr_matrix_yyz = ptr_octu_matrix->getMatrix("YYZ");

                        auto ptr_matrix_yzz = ptr_octu_matrix->getMatrix("YZZ");

                        auto ptr_matrix_zzz = ptr_octu_matrix->getMatrix("ZZZ");

                        if (task[0] == task[1])
                        {
                            const auto gto_block = ptr_gto_blocks[task[0]];

                            const auto angmom = gto_block.getAngularMomentum();

                            auto ptr_submatrix_xxx = ptr_matrix_xxx->getSubMatrix({angmom, angmom});

                            auto ptr_submatrix_xxy = ptr_matrix_xxy->getSubMatrix({angmom, angmom});

                            auto ptr_submatrix_xxz = ptr_matrix_xxz->getSubMatrix({angmom, angmom});

                            auto ptr_submatrix_xyy = ptr_matrix_xyy->getSubMatrix({angmom, angmom});

                            auto ptr_submatrix_xyz = ptr_matrix_xyz->getSubMatrix({angmom, angmom});

                            auto ptr_submatrix_xzz = ptr_matrix_xzz->getSubMatrix({angmom, angmom});

                            auto ptr_submatrix_yyy = ptr_matrix_yyy->getSubMatrix({angmom, angmom});

                            auto ptr_submatrix_yyz = ptr_matrix_yyz->getSubMatrix({angmom, angmom});

                            auto ptr_submatrix_yzz = ptr_matrix_yzz->getSubMatrix({angmom, angmom});

                            auto ptr_submatrix_zzz = ptr_matrix_zzz->getSubMatrix({angmom, angmom});

                            octufunc::compute(ptr_submatrix_xxx,
                                              ptr_submatrix_xxy,
                                              ptr_submatrix_xxz,
                                              ptr_submatrix_xyy,
                                              ptr_submatrix_xyz,
                                              ptr_submatrix_xzz,
                                              ptr_submatrix_yyy,
                                              ptr_submatrix_yyz,
                                              ptr_submatrix_yzz,
                                              ptr_submatrix_zzz,
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

                            auto ptr_submatrix_xxx = ptr_matrix_xxx->getSubMatrix({bra_angmom, ket_angmom});

                            auto ptr_submatrix_xxy = ptr_matrix_xxy->getSubMatrix({bra_angmom, ket_angmom});

                            auto ptr_submatrix_xxz = ptr_matrix_xxz->getSubMatrix({bra_angmom, ket_angmom});

                            auto ptr_submatrix_xyy = ptr_matrix_xyy->getSubMatrix({bra_angmom, ket_angmom});

                            auto ptr_submatrix_xyz = ptr_matrix_xyz->getSubMatrix({bra_angmom, ket_angmom});

                            auto ptr_submatrix_xzz = ptr_matrix_xzz->getSubMatrix({bra_angmom, ket_angmom});

                            auto ptr_submatrix_yyy = ptr_matrix_yyy->getSubMatrix({bra_angmom, ket_angmom});

                            auto ptr_submatrix_yyz = ptr_matrix_yyz->getSubMatrix({bra_angmom, ket_angmom});

                            auto ptr_submatrix_yzz = ptr_matrix_yzz->getSubMatrix({bra_angmom, ket_angmom});

                            auto ptr_submatrix_zzz = ptr_matrix_zzz->getSubMatrix({bra_angmom, ket_angmom});

                            const auto ang_order = ptr_matrix_xxx->isAngularOrder({bra_angmom, ket_angmom});

                            octufunc::compute(ptr_submatrix_xxx,
                                              ptr_submatrix_xxy,
                                              ptr_submatrix_xxz,
                                              ptr_submatrix_xyy,
                                              ptr_submatrix_xyz,
                                              ptr_submatrix_xzz,
                                              ptr_submatrix_yyy,
                                              ptr_submatrix_yyz,
                                              ptr_submatrix_yzz,
                                              ptr_submatrix_zzz,
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

    return octu_matrix;
}
