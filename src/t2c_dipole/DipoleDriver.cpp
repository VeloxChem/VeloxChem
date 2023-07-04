#include "DipoleDriver.hpp"

#include "Matrix.hpp"
#include "MatrixType.hpp"
#include "GtoFunc.hpp"
#include "MatrixFunc.hpp"
#include "OpenMPFunc.hpp"
#include "DipoleFunc.hpp"

auto
CDipoleDriver::compute(const CMolecularBasis& basis,
                       const CMolecule& molecule,
                       const TPoint3D& point) const -> CMatrices
{
    CMatrices dip_matrix;
    
    dip_matrix.add(matfunc::makeMatrix(basis, mat_t::symm), "x");
    
    dip_matrix.add(matfunc::makeMatrix(basis, mat_t::symm), "y");
    
    dip_matrix.add(matfunc::makeMatrix(basis, mat_t::symm), "z");
    
    dip_matrix.zero();

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto work_groups = omp::makeWorkGroup(gto_blocks);

    // prepare pointers for OMP parallel region

    auto ptr_gto_blocks = gto_blocks.data();

    auto ptr_work_groups = work_groups.data();

    auto ptr_dip_matrix = &dip_matrix;

    // execute OMP tasks with static scheduling

    omp::setStaticScheduler();

    const auto ntasks = work_groups.size();

#pragma omp parallel num_threads(ntasks) shared(ntasks, ptr_gto_blocks, ptr_work_groups, ptr_dip_matrix, point)
    {
#pragma omp single nowait
        {
            for (size_t i = 0; i < ntasks; i++)
            {
#pragma omp task firstprivate(i)
                {
                    for (const auto& task : ptr_work_groups[i])
                    {
                        const auto mat_type = mat_t::symm;

                        if (task[0] == task[1])
                        {
                            const auto gto_block = ptr_gto_blocks[task[0]];

                            const auto angmom = gto_block.getAngularMomentum();
                            
                            auto ptr_matrix_x = ptr_dip_matrix->getMatrix("X");
                            
                            auto ptr_matrix_y = ptr_dip_matrix->getMatrix("Y");
                            
                            auto ptr_matrix_z = ptr_dip_matrix->getMatrix("Z");
                        
                            auto ptr_submatrix_x = ptr_matrix_x->getSubMatrix({angmom, angmom});
                            
                            auto ptr_submatrix_y = ptr_matrix_y->getSubMatrix({angmom, angmom});
                            
                            auto ptr_submatrix_z = ptr_matrix_z->getSubMatrix({angmom, angmom});

                            dipfunc::compute(ptr_submatrix_x, ptr_submatrix_y, ptr_submatrix_z, point, gto_block, angmom, task[2], task[3]);
                        }
                        else
                        {
                            const auto bra_gto_block = ptr_gto_blocks[task[0]];

                            const auto ket_gto_block = ptr_gto_blocks[task[1]];

                            const auto bra_angmom = bra_gto_block.getAngularMomentum();

                            const auto ket_angmom = ket_gto_block.getAngularMomentum();

                            auto ptr_matrix_x = ptr_dip_matrix->getMatrix("X");
                            
                            auto ptr_matrix_y = ptr_dip_matrix->getMatrix("Y");
                            
                            auto ptr_matrix_z = ptr_dip_matrix->getMatrix("Z");
                        
                            auto ptr_submatrix_x = ptr_matrix_x->getSubMatrix({bra_angmom, ket_angmom});
                            
                            auto ptr_submatrix_y = ptr_matrix_y->getSubMatrix({bra_angmom, ket_angmom});
                            
                            auto ptr_submatrix_z = ptr_matrix_z->getSubMatrix({bra_angmom, ket_angmom});

                            const auto ang_order = ptr_matrix_x->isAngularOrder({bra_angmom, ket_angmom});

                            dipfunc::compute(
                                             ptr_submatrix_x, ptr_submatrix_y, ptr_submatrix_z, point, bra_gto_block, ket_gto_block, bra_angmom, ket_angmom, ang_order, task[2], task[3], mat_type);
                        }
                    }
                }
            }
        }
    }

    return dip_matrix;
}
