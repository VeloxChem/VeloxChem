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

#include "XCIntegratorForPDFT.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>

#include "DenseMatrix.hpp"
#include "DftSubMatrix.hpp"
#include "ErrorHandler.hpp"
#include "FunctionalParser.hpp"
#include "GtoFunc.hpp"
#include "GtoValues.hpp"
#include "MathFunc.hpp"
#include "MultiTimer.hpp"
#include "PairDensityGridGenerator.hpp"
#include "Prescreener.hpp"
#include "SerialDenseLinearAlgebra.hpp"
#include "StringFormat.hpp"
#include "Timer.hpp"
#include "XCIntegratorForGGA.hpp"
#include "XCIntegratorForLDA.hpp"

namespace xcintpdft {  // xcintpdft namespace

void
integrateVxcPDFTForLDA(CAOKohnShamMatrix&              aoFockMatrix,
                       CDense4DTensor&                 tensorWxc,
                       const CMolecule&                molecule,
                       const CMolecularBasis&          basis,
                       const double*                   densityMatrixPointer,
                       const CDenseMatrix&             twoBodyDensityMatrix,
                       const CDenseMatrix&             activeMOs,
                       const CMolecularGrid&           molecularGrid,
                       const double                    screeningThresholdForGTOValues,
                       const CXCPairDensityFunctional& xcFunctional,
                       const double                    rs_omega)
{
    CMultiTimer timer;

    timer.start("Total timing");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // Set up Fock matrix

    aoFockMatrix.zero();

    // initial values for XC energy and number of electrons

    double nele = 0.0, xcene = 0.0;

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    // set up number of grid blocks

    const auto n_boxes = counts.size();

    const auto n_gto_blocks = gto_blocks.size();

    // set up pointers to OMP data

    auto ptr_counts = counts.data();

    auto ptr_displacements = displacements.data();

    auto ptr_gto_blocks = gto_blocks.data();

    auto ptr_aoFockMatrix = &aoFockMatrix;
    auto ptr_tensorWxc = &tensorWxc;

    auto ptr_twoBodyDensityMatrix = &twoBodyDensityMatrix;
    auto ptr_activeMOs = &activeMOs;

    auto ptr_xcFunctional = &xcFunctional;

#pragma omp parallel shared(ptr_counts, ptr_displacements, xcoords, ycoords, zcoords, \
                            ptr_gto_blocks, densityMatrixPointer, ptr_aoFockMatrix, ptr_tensorWxc, \
                            ptr_twoBodyDensityMatrix, ptr_activeMOs, ptr_xcFunctional, \
                            n_boxes, n_gto_blocks, naos, nele, xcene)
    {

#pragma omp single nowait
    {

    for (size_t box_id = 0; box_id < n_boxes; box_id++)
    {

    #pragma omp task firstprivate(box_id)
    {
        auto thread_id = omp_get_thread_num();

        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // prescreening

        omptimers[thread_id].start("GTO pre-screening");

        std::vector<std::vector<int>> cgto_mask_blocks, pre_ao_inds_blocks;

        std::vector<int> aoinds;

        cgto_mask_blocks.reserve(n_gto_blocks);

        pre_ao_inds_blocks.reserve(n_gto_blocks);

        aoinds.reserve(naos); 

        for (size_t i = 0; i < n_gto_blocks; i++)
        {
            // 0th order GTO derivative
            auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(ptr_gto_blocks[i], 0, screeningThresholdForGTOValues, boxdim);

            cgto_mask_blocks.push_back(cgto_mask);

            pre_ao_inds_blocks.push_back(pre_ao_inds);

            for (const auto nu : pre_ao_inds)
            {
                aoinds.push_back(nu);
            }
        }

        const auto aocount = static_cast<int>(aoinds.size());

        omptimers[thread_id].stop("GTO pre-screening");

        if (aocount > 0)
        {
            omptimers[thread_id].start("Density matrix slicing");

            auto sub_dens_mat_a = dftsubmat::getSubDensityMatrix(densityMatrixPointer, aoinds, naos);

            auto sub_active_mos = dftsubmat::getSubMatrixByColumnSlicing(*ptr_activeMOs, aoinds, naos);

            omptimers[thread_id].stop("Density matrix slicing");

            // GTO values on grid points

            omptimers[thread_id].start("gtoeval");

            CDenseMatrix mat_chi(aocount, npoints);

            const auto grid_x_ptr = xcoords + gridblockpos;
            const auto grid_y_ptr = ycoords + gridblockpos;
            const auto grid_z_ptr = zcoords + gridblockpos;

            std::vector<double> grid_x(grid_x_ptr, grid_x_ptr + npoints);
            std::vector<double> grid_y(grid_y_ptr, grid_y_ptr + npoints);
            std::vector<double> grid_z(grid_z_ptr, grid_z_ptr + npoints);

            // go through GTO blocks

            for (size_t i_block = 0, idx = 0; i_block < n_gto_blocks; i_block++)
            {
                const auto& gto_block = ptr_gto_blocks[i_block];

                const auto& cgto_mask = cgto_mask_blocks[i_block];

                const auto& pre_ao_inds = pre_ao_inds_blocks[i_block];

                auto cmat = gtoval::get_gto_values_for_lda(gto_block, grid_x, grid_y, grid_z, cgto_mask);

                if (cmat.is_empty()) continue;

                auto submat_ptr = cmat.sub_matrix({0, 0});

                auto submat_data = submat_ptr->data();

                for (int nu = 0; nu < static_cast<int>(pre_ao_inds.size()); nu++, idx++)
                {
                    std::memcpy(mat_chi.row(idx), submat_data + nu * npoints, npoints * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");

            omptimers[thread_id].start("Generate density grid");

            auto local_xcfunc = CXCPairDensityFunctional(*ptr_xcFunctional);

            std::vector<double> local_weights_data(weights + gridblockpos, weights + gridblockpos + npoints);

            std::vector<double> rho_data(2 * npoints);

            std::vector<double> exc_data(1 * npoints);
            std::vector<double> vrho_data(2 * npoints);

            auto local_weights = local_weights_data.data();

            auto rho  = rho_data.data();

            auto exc  = exc_data.data();
            auto vrho = vrho_data.data();

            // generate density and on-top pair density on the grid

            pairdengridgen::serialGeneratePairDensityForLDA(rho, mat_chi, sub_dens_mat_a, sub_active_mos, *ptr_twoBodyDensityMatrix);

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            local_xcfunc.compute_exc_vxc_for_plda(npoints, rho, exc, vrho, rs_omega);

            omptimers[thread_id].stop("XC functional eval.");

            // compute partial contribution to Vxc and Wxc

            auto partial_mat_Vxc = xcintpdft::integratePartialVxcFockForLDA(
                    local_weights, mat_chi, vrho, omptimers[thread_id]);

            auto partial_tensorWxc = xcintpdft::integratePartialWxcFockForPLDA(
                    local_weights, mat_chi, sub_active_mos, vrho, omptimers[thread_id]);

            omptimers[thread_id].start("Vxc and Wxc dist.");

            double local_nele = 0.0, local_xcene = 0.0;

            for (int g = 0; g < npoints; g++)
            {
                auto rho_total = rho[2 * g + 0];

                local_nele += local_weights[g] * rho_total;

                local_xcene += local_weights[g] * exc[g] * rho_total;
            }

            #pragma omp critical
            {
                nele += local_nele;

                xcene += local_xcene;

                dftsubmat::distributeSubMatrixToKohnSham(*ptr_aoFockMatrix, partial_mat_Vxc, aoinds);

                dftsubmat::distributeSubmatrixTo4DTensor(*ptr_tensorWxc, partial_tensorWxc, aoinds);
            }

            omptimers[thread_id].stop("Vxc and Wxc dist.");
        }
    }
    }
    }
    }

    aoFockMatrix.setNumberOfElectrons(nele);

    aoFockMatrix.setExchangeCorrelationEnergy(xcene);

    timer.stop("Total timing");

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;
    // std::cout << "OpenMP timing" << std::endl;
    // for (int thread_id = 0; thread_id < nthreads; thread_id++)
    // {
    //     std::cout << "Thread " << thread_id << std::endl;
    //     std::cout << omptimers[thread_id].getSummary() << std::endl;
    // }
}

void
integrateVxcPDFTForGGA(CAOKohnShamMatrix&              aoFockMatrix,
                       CDense4DTensor&                 tensorWxc,
                       const CMolecule&                molecule,
                       const CMolecularBasis&          basis,
                       const double*                   densityMatrixPointer,
                       const CDenseMatrix&             twoBodyDensityMatrix,
                       const CDenseMatrix&             activeMOs,
                       const CMolecularGrid&           molecularGrid,
                       const double                    screeningThresholdForGTOValues,
                       const CXCPairDensityFunctional& xcFunctional,
                       const double                    rs_omega)
{
    CMultiTimer timer;

    timer.start("Total timing");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // Set up Fock matrix

    aoFockMatrix.zero();

    // initial values for XC energy and number of electrons

    double nele = 0.0, xcene = 0.0;

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    // set up number of grid blocks

    const auto n_boxes = counts.size();

    const auto n_gto_blocks = gto_blocks.size();

    // set up pointers to OMP data

    auto ptr_counts = counts.data();

    auto ptr_displacements = displacements.data();

    auto ptr_gto_blocks = gto_blocks.data();

    auto ptr_aoFockMatrix = &aoFockMatrix;
    auto ptr_tensorWxc = &tensorWxc;

    auto ptr_twoBodyDensityMatrix = &twoBodyDensityMatrix;
    auto ptr_activeMOs = &activeMOs;

    auto ptr_xcFunctional = &xcFunctional;

#pragma omp parallel shared(ptr_counts, ptr_displacements, xcoords, ycoords, zcoords, \
                            ptr_gto_blocks, densityMatrixPointer, ptr_aoFockMatrix, ptr_tensorWxc, \
                            ptr_twoBodyDensityMatrix, ptr_activeMOs, ptr_xcFunctional, \
                            n_boxes, n_gto_blocks, naos, nele, xcene)
    {

#pragma omp single nowait
    {

    for (size_t box_id = 0; box_id < n_boxes; box_id++)
    {

    #pragma omp task firstprivate(box_id)
    {
        auto thread_id = omp_get_thread_num();

        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // prescreening

        omptimers[thread_id].start("GTO pre-screening");

        std::vector<std::vector<int>> cgto_mask_blocks, pre_ao_inds_blocks;

        std::vector<int> aoinds;

        cgto_mask_blocks.reserve(n_gto_blocks);

        pre_ao_inds_blocks.reserve(n_gto_blocks);

        aoinds.reserve(naos); 

        for (size_t i = 0; i < n_gto_blocks; i++)
        {
            // 1st order GTO derivative
            auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(ptr_gto_blocks[i], 1, screeningThresholdForGTOValues, boxdim);

            cgto_mask_blocks.push_back(cgto_mask);

            pre_ao_inds_blocks.push_back(pre_ao_inds);

            for (const auto nu : pre_ao_inds)
            {
                aoinds.push_back(nu);
            }
        }

        const auto aocount = static_cast<int>(aoinds.size());

        omptimers[thread_id].stop("GTO pre-screening");

        if (aocount > 0)
        {
            omptimers[thread_id].start("Density matrix slicing");

            auto sub_dens_mat_a = dftsubmat::getSubDensityMatrix(densityMatrixPointer, aoinds, naos);

            auto sub_active_mos = dftsubmat::getSubMatrixByColumnSlicing(*ptr_activeMOs, aoinds, naos);

            omptimers[thread_id].stop("Density matrix slicing");

            // GTO values on grid points

            omptimers[thread_id].start("gtoeval");

            CDenseMatrix mat_chi(aocount, npoints);
            CDenseMatrix mat_chi_x(aocount, npoints);
            CDenseMatrix mat_chi_y(aocount, npoints);
            CDenseMatrix mat_chi_z(aocount, npoints);

            const auto grid_x_ptr = xcoords + gridblockpos;
            const auto grid_y_ptr = ycoords + gridblockpos;
            const auto grid_z_ptr = zcoords + gridblockpos;

            std::vector<double> grid_x(grid_x_ptr, grid_x_ptr + npoints);
            std::vector<double> grid_y(grid_y_ptr, grid_y_ptr + npoints);
            std::vector<double> grid_z(grid_z_ptr, grid_z_ptr + npoints);

            // go through GTO blocks

            for (size_t i_block = 0, idx = 0; i_block < n_gto_blocks; i_block++)
            {
                const auto& gto_block = ptr_gto_blocks[i_block];

                const auto& cgto_mask = cgto_mask_blocks[i_block];

                const auto& pre_ao_inds = pre_ao_inds_blocks[i_block];

                auto cmat = gtoval::get_gto_values_for_gga(gto_block, grid_x, grid_y, grid_z, cgto_mask);

                if (cmat.is_empty()) continue;

                auto submat_0_ptr = cmat.sub_matrix({0, 0});
                auto submat_x_ptr = cmat.sub_matrix({1, 0});
                auto submat_y_ptr = cmat.sub_matrix({1, 1});
                auto submat_z_ptr = cmat.sub_matrix({1, 2});

                auto submat_0_data = submat_0_ptr->data();
                auto submat_x_data = submat_x_ptr->data();
                auto submat_y_data = submat_y_ptr->data();
                auto submat_z_data = submat_z_ptr->data();

                for (int nu = 0; nu < static_cast<int>(pre_ao_inds.size()); nu++, idx++)
                {
                    std::memcpy(mat_chi.row(idx), submat_0_data + nu * npoints, npoints * sizeof(double));
                    std::memcpy(mat_chi_x.row(idx), submat_x_data + nu * npoints, npoints * sizeof(double));
                    std::memcpy(mat_chi_y.row(idx), submat_y_data + nu * npoints, npoints * sizeof(double));
                    std::memcpy(mat_chi_z.row(idx), submat_z_data + nu * npoints, npoints * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");

            omptimers[thread_id].start("Generate density grid");

            auto local_xcfunc = CXCPairDensityFunctional(*ptr_xcFunctional);

            std::vector<double> local_weights_data(weights + gridblockpos, weights + gridblockpos + npoints);

            std::vector<double> rho_data(2 * npoints);
            std::vector<double> rhograd_data(2 * 3 * npoints);
            std::vector<double> sigma_data(3 * npoints);

            std::vector<double> exc_data(1 * npoints);
            std::vector<double> vrho_data(2 * npoints);
            std::vector<double> vsigma_data(3 * npoints);

            auto local_weights = local_weights_data.data();

            auto rho     = rho_data.data();
            auto rhograd = rhograd_data.data();
            auto sigma   = sigma_data.data();

            auto exc    = exc_data.data();
            auto vrho   = vrho_data.data();
            auto vsigma = vsigma_data.data();

            // generate density and on-top pair density on the grid

            pairdengridgen::serialGeneratePairDensityForGGA(
                rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat_a, sub_active_mos, *ptr_twoBodyDensityMatrix);

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            local_xcfunc.compute_exc_vxc_for_pgga(npoints, rho, sigma, exc, vrho, vsigma, rs_omega);

            omptimers[thread_id].stop("XC functional eval.");

            // compute partial contribution to Vxc and Wxc

            // TODO (MGD) gradient not correct for vsigma[1] and vsigma[2]

            auto partial_mat_Vxc = xcintpdft::integratePartialVxcFockForGGA(
                    local_weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rhograd, vrho, vsigma, omptimers[thread_id]);

            auto partial_tensorWxc = xcintpdft::integratePartialWxcFockForPLDA(
                    local_weights, mat_chi, sub_active_mos, vrho, omptimers[thread_id]);

            omptimers[thread_id].start("Vxc and Wxc dist.");

            double local_nele = 0.0, local_xcene = 0.0;

            for (int g = 0; g < npoints; g++)
            {
                auto rho_total = rho[2 * g + 0];

                local_nele += local_weights[g] * rho_total;

                local_xcene += local_weights[g] * exc[g] * rho_total;
            }

            #pragma omp critical
            {
                nele += local_nele;

                xcene += local_xcene;

                dftsubmat::distributeSubMatrixToKohnSham(*ptr_aoFockMatrix, partial_mat_Vxc, aoinds);

                dftsubmat::distributeSubmatrixTo4DTensor(*ptr_tensorWxc, partial_tensorWxc, aoinds);
            }

            omptimers[thread_id].stop("Vxc and Wxc dist.");
        }
    }
    }
    }
    }

    aoFockMatrix.setNumberOfElectrons(nele);

    aoFockMatrix.setExchangeCorrelationEnergy(xcene);

    timer.stop("Total timing");

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;
    // std::cout << "OpenMP timing" << std::endl;
    // for (int thread_id = 0; thread_id < nthreads; thread_id++)
    // {
    //     std::cout << "Thread " << thread_id << std::endl;
    //     std::cout << omptimers[thread_id].getSummary() << std::endl;
    // }
}

auto
integratePartialVxcFockForLDA(const double* weights, const CDenseMatrix& gtoValues, const double* vrho, CMultiTimer& timer) -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    // GTO values on grid points

    timer.start("Vxc matrix G");

    auto chi_val = gtoValues.values();

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);

    auto G_val = mat_G.values();

    {
        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd
            for (int g = 0; g < npoints; g++)
            {
                G_val[nu_offset + g] = weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];
            }
        }
    }

    timer.stop("Vxc matrix G");

    timer.start("Vxc matrix matmul");

    auto mat_Vxc = sdenblas::serialMultABt(gtoValues, mat_G);

    timer.stop("Vxc matrix matmul");

    return mat_Vxc;
}

auto
integratePartialVxcFockForGGA(const double*       weights,
                              const CDenseMatrix& gtoValues,
                              const CDenseMatrix& gtoValuesX,
                              const CDenseMatrix& gtoValuesY,
                              const CDenseMatrix& gtoValuesZ,
                              const double*       rhograd,
                              const double*       vrho,
                              const double*       vsigma,
                              CMultiTimer&        timer) -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    timer.start("Vxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);
    CDenseMatrix mat_G_gga(naos, npoints);

    auto G_val     = mat_G.values();
    auto G_gga_val = mat_G_gga.values();

    auto chi_val   = gtoValues.values();
    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    {
        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd
            for (int g = 0; g < npoints; g++)
            {
                auto vx = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                auto vy = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                auto vz = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                G_val[nu_offset + g] = weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];

                G_gga_val[nu_offset + g] =
                    weights[g] * (vx * chi_x_val[nu_offset + g] + vy * chi_y_val[nu_offset + g] + vz * chi_z_val[nu_offset + g]);
            }
        }
    }

    timer.stop("Vxc matrix G");

    // Note that we use matrix-matrix multiplication only once, and symmetrize
    // the result. This is because the density matrix is symmetric, and the
    // Kohn-Sham matrix from mat_G is also symmetric. Formally only the
    // mat_G_gga contribution should be symmetrized.

    timer.start("Vxc matrix matmul");

    auto mat_Vxc = sdenblas::serialMultABt(gtoValues, sdenblas::serialAddAB(mat_G, mat_G_gga, 2.0));

    mat_Vxc.symmetrizeAndScale(0.5);

    timer.stop("Vxc matrix matmul");

    return mat_Vxc;
}

CDenseMatrix
integratePartialWxcFockForPLDA(const double*       weights,
                               const CDenseMatrix& gtoValues,
                               const CDenseMatrix& activeMOs,
                               const double*       vrho,
                               CMultiTimer&        timer)
{
    const auto npoints = gtoValues.getNumberOfColumns();

    timer.start("Wxc matrix");

    auto nActive = activeMOs.getNumberOfRows();

    CDenseMatrix mos_on_grid;

    if (nActive > 0)
    {
        mos_on_grid = sdenblas::serialMultAB(activeMOs, gtoValues);
    }

    // created empty partial mat_W
    CDenseMatrix matrixWxc(nActive * nActive * nActive, npoints);
    auto         W_val = matrixWxc.values();

    {
        for (int j = 0; j < nActive; j++)
        {
            auto mo_j = mos_on_grid.row(j);

            for (int k = 0; k < nActive; k++)
            {
                auto jk   = j * nActive + k;
                auto mo_k = mos_on_grid.row(k);

                for (int l = 0; l < nActive; l++)
                {
                    auto jkl  = jk * nActive + l;
                    auto mo_l = mos_on_grid.row(l);

                    #pragma omp simd
                    for (int g = 0; g < npoints; g++)
                    {
                        W_val[jkl * npoints + g] = weights[g] * vrho[2 * g + 1] * mo_j[g] * mo_k[g] * mo_l[g];
                    }
                }
            }
        }
    }

    auto tensorWxc = sdenblas::serialMultABt(gtoValues, matrixWxc);

    timer.stop("Wxc matrix");

    return tensorWxc;
}

}  // namespace xcintpdft
