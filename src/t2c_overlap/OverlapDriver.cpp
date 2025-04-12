//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
