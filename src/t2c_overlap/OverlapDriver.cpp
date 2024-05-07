//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

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


auto
COverlapDriver::compute(const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis, const CMolecule& molecule) const -> CMatrix
{
    auto ovl_matrix = matfunc::makeMatrix(bra_basis, ket_basis);

    ovl_matrix.zero();
    
    const auto bra_gto_blocks = gtofunc::makeGtoBlocks(bra_basis, molecule);
    
    const auto ket_gto_blocks = gtofunc::makeGtoBlocks(ket_basis, molecule);
    
    const auto work_groups = omp::makeWorkGroup(bra_gto_blocks, ket_gto_blocks);
    
    // prepare pointers for OMP parallel region

    auto ptr_bra_gto_blocks = bra_gto_blocks.data();
    
    auto ptr_ket_gto_blocks = ket_gto_blocks.data();

    auto ptr_work_groups = work_groups.data();

    auto ptr_ovl_matrix = &ovl_matrix;

    // execute OMP tasks with static scheduling

    omp::setStaticScheduler();

    const auto ntasks = work_groups.size();
    
#pragma omp parallel num_threads(ntasks) shared(ntasks, ptr_bra_gto_blocks, ptr_ket_gto_blocks, ptr_work_groups, ptr_ovl_matrix)
    {
#pragma omp single nowait
        {
            for (size_t i = 0; i < ntasks; i++)
            {
#pragma omp task firstprivate(i)
                {
                    for (const auto& task : ptr_work_groups[i])
                    {
                        const auto bra_gto_block = ptr_bra_gto_blocks[task[0]];

                        const auto ket_gto_block = ptr_ket_gto_blocks[task[1]];

                        const auto bra_angmom = bra_gto_block.getAngularMomentum();

                        const auto ket_angmom = ket_gto_block.getAngularMomentum();

                        auto ptr_submatrix = ptr_ovl_matrix->getSubMatrix({bra_angmom, ket_angmom});

                        ovlfunc::compute(ptr_submatrix, bra_gto_block, ket_gto_block, bra_angmom, ket_angmom, true, task[2], task[3], mat_t::gen);
                    }
                }
            }
        }
    }
    
    return ovl_matrix;
}
