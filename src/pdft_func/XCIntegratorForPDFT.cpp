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

#include "XCIntegratorForPDFT.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>

#include "DenseLinearAlgebra.hpp"
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
                       const CXCPairDensityFunctional& xcFunctional)
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // Set up Fock matrix

    aoFockMatrix.zero();

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    // density and functional derivatives

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(2 * max_npoints_per_box);

    std::vector<double> exc_data(1 * max_npoints_per_box);
    std::vector<double> vrho_data(2 * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();

    auto exc  = exc_data.data();
    auto vrho = vrho_data.data();

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

    timer.stop("Preparation");

    for (int box_id = 0; box_id < counts.size(); box_id++)
    {
        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // prescreening

        timer.start("GTO pre-screening");

        std::vector<std::vector<int>> cgto_mask_blocks, pre_ao_inds_blocks;

        std::vector<int> aoinds;

        for (const auto& gto_block : gto_blocks)
        {
            // 0th order GTO derivative
            auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(gto_block, 0, screeningThresholdForGTOValues, boxdim);

            cgto_mask_blocks.push_back(cgto_mask);

            pre_ao_inds_blocks.push_back(pre_ao_inds);

            for (const auto nu : pre_ao_inds)
            {
                aoinds.push_back(nu);
            }
        }

        const auto aocount = static_cast<int>(aoinds.size());

        timer.stop("GTO pre-screening");

        if (aocount == 0) continue;

        // GTO values on grid points

        timer.start("OMP GTO evaluation");

        CDenseMatrix mat_chi(aocount, npoints);

#pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            omptimers[thread_id].start("gtoeval");

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            const auto grid_x_ptr = xcoords + gridblockpos + grid_batch_offset;
            const auto grid_y_ptr = ycoords + gridblockpos + grid_batch_offset;
            const auto grid_z_ptr = zcoords + gridblockpos + grid_batch_offset;

            std::vector<double> grid_x(grid_x_ptr, grid_x_ptr + grid_batch_size);
            std::vector<double> grid_y(grid_y_ptr, grid_y_ptr + grid_batch_size);
            std::vector<double> grid_z(grid_z_ptr, grid_z_ptr + grid_batch_size);

            // go through GTO blocks

            for (size_t i_block = 0, idx = 0; i_block < gto_blocks.size(); i_block++)
            {
                const auto& gto_block = gto_blocks[i_block];

                const auto& cgto_mask = cgto_mask_blocks[i_block];

                const auto& pre_ao_inds = pre_ao_inds_blocks[i_block];

                auto cmat = gtoval::get_gto_values_for_lda(gto_block, grid_x, grid_y, grid_z, cgto_mask);

                if (cmat.is_empty()) continue;

                auto submat_ptr = cmat.sub_matrix({0, 0});

                auto submat_data = submat_ptr->data();

                for (int nu = 0; nu < static_cast<int>(pre_ao_inds.size()); nu++, idx++)
                {
                    std::memcpy(mat_chi.row(idx) + grid_batch_offset, submat_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");
        }

        timer.stop("OMP GTO evaluation");

        // generate sub density matrix and MO coefficients

        timer.start("Density matrix slicing");

        auto sub_dens_mat_a = dftsubmat::getSubDensityMatrix(densityMatrixPointer, aoinds, naos);

        auto sub_active_mos = dftsubmat::getSubMatrixByColumnSlicing(activeMOs, aoinds, naos);

        timer.stop("Density matrix slicing");

        // generate density and on-top pair density on the grid

        pairdengridgen::generatePairDensityForLDA(rho, mat_chi, sub_dens_mat_a, sub_active_mos, twoBodyDensityMatrix, timer);

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_plda(npoints, rho, exc, vrho);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        auto partial_mat_Vxc = xcintlda::integratePartialVxcFockForLDA(local_weights, mat_chi, vrho, timer);

        timer.start("Vxc matrix dist.");

        dftsubmat::distributeSubMatrixToKohnSham(aoFockMatrix, partial_mat_Vxc, aoinds);

        timer.stop("Vxc matrix dist.");

        auto partial_tensorWxc = xcintpdft::integratePartialWxcFockForPLDA(local_weights, mat_chi, sub_active_mos, vrho, timer);

        timer.start("Wxc matrix dist.");

        dftsubmat::distributeSubmatrixTo4DTensor(tensorWxc, partial_tensorWxc, aoinds);

        timer.stop("Wxc matrix dist.");

        // compute partial contribution to XC energy

        timer.start("XC energy");

        for (int g = 0; g < npoints; g++)
        {
            auto rho_total = rho[2 * g + 0];

            nele += local_weights[g] * rho_total;

            xcene += local_weights[g] * exc[g] * rho_total;
        }

        timer.stop("XC energy");
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
                       const CXCPairDensityFunctional& xcFunctional)
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // Set up Fock matrix

    aoFockMatrix.zero();

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    // density and functional derivatives

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(2 * max_npoints_per_box);
    std::vector<double> rhograd_data(2 * 3 * max_npoints_per_box);
    std::vector<double> sigma_data(3 * max_npoints_per_box);

    std::vector<double> exc_data(1 * max_npoints_per_box);
    std::vector<double> vrho_data(2 * max_npoints_per_box);
    std::vector<double> vsigma_data(3 * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();

    auto exc    = exc_data.data();
    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();

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

    timer.stop("Preparation");

    for (int box_id = 0; box_id < counts.size(); box_id++)
    {
        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // prescreening

        timer.start("GTO pre-screening");

        std::vector<std::vector<int>> cgto_mask_blocks, pre_ao_inds_blocks;

        std::vector<int> aoinds;

        for (const auto& gto_block : gto_blocks)
        {
            // 1st order GTO derivative
            auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(gto_block, 1, screeningThresholdForGTOValues, boxdim);

            cgto_mask_blocks.push_back(cgto_mask);

            pre_ao_inds_blocks.push_back(pre_ao_inds);

            for (const auto nu : pre_ao_inds)
            {
                aoinds.push_back(nu);
            }
        }

        const auto aocount = static_cast<int>(aoinds.size());

        timer.stop("GTO pre-screening");

        if (aocount == 0) continue;

        // GTO values on grid points

        timer.start("OMP GTO evaluation");

        CDenseMatrix mat_chi(aocount, npoints);
        CDenseMatrix mat_chi_x(aocount, npoints);
        CDenseMatrix mat_chi_y(aocount, npoints);
        CDenseMatrix mat_chi_z(aocount, npoints);

#pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            omptimers[thread_id].start("gtoeval");

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            const auto grid_x_ptr = xcoords + gridblockpos + grid_batch_offset;
            const auto grid_y_ptr = ycoords + gridblockpos + grid_batch_offset;
            const auto grid_z_ptr = zcoords + gridblockpos + grid_batch_offset;

            std::vector<double> grid_x(grid_x_ptr, grid_x_ptr + grid_batch_size);
            std::vector<double> grid_y(grid_y_ptr, grid_y_ptr + grid_batch_size);
            std::vector<double> grid_z(grid_z_ptr, grid_z_ptr + grid_batch_size);

            // go through GTO blocks

            for (size_t i_block = 0, idx = 0; i_block < gto_blocks.size(); i_block++)
            {
                const auto& gto_block = gto_blocks[i_block];

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
                    std::memcpy(mat_chi.row(idx) + grid_batch_offset, submat_0_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_x.row(idx) + grid_batch_offset, submat_x_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_y.row(idx) + grid_batch_offset, submat_y_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_z.row(idx) + grid_batch_offset, submat_z_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");
        }

        timer.stop("OMP GTO evaluation");

        // generate sub density matrix and MO coefficients

        timer.start("Density matrix slicing");

        auto sub_dens_mat_a = dftsubmat::getSubDensityMatrix(densityMatrixPointer, aoinds, naos);

        auto sub_active_mos = dftsubmat::getSubMatrixByColumnSlicing(activeMOs, aoinds, naos);

        timer.stop("Density matrix slicing");

        // generate density and on-top pair density on the grid

        pairdengridgen::generatePairDensityForGGA(
            rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat_a, sub_active_mos, twoBodyDensityMatrix, timer);

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_pgga(npoints, rho, sigma, exc, vrho, vsigma);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // Compute Vxc matrix

        // TODO (MGD) gradient not correct for vsigma[1] and vsigma[2]

        auto partial_mat_Vxc =
            xcintgga::integratePartialVxcFockForGGA(local_weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rhograd, vrho, vsigma, timer);

        timer.start("Vxc matrix dist.");

        dftsubmat::distributeSubMatrixToKohnSham(aoFockMatrix, partial_mat_Vxc, aoinds);

        timer.stop("Vxc matrix dist.");

        // Compute Wxc tensor

        auto partial_tensorWxc = xcintpdft::integratePartialWxcFockForPLDA(local_weights, mat_chi, sub_active_mos, vrho, timer);

        timer.start("Wxc matrix dist.");

        dftsubmat::distributeSubmatrixTo4DTensor(tensorWxc, partial_tensorWxc, aoinds);

        timer.stop("Wxc matrix dist.");

        // compute partial contribution to XC energy

        timer.start("XC energy");

        for (int g = 0; g < npoints; g++)
        {
            auto rho_total = rho[2 * g + 0];

            nele += local_weights[g] * rho_total;

            xcene += local_weights[g] * exc[g] * rho_total;
        }

        timer.stop("XC energy");
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
        mos_on_grid = denblas::multAB(activeMOs, gtoValues);
    }

    // created empty partial mat_W
    CDenseMatrix matrixWxc(nActive * nActive * nActive, npoints);
    auto         W_val = matrixWxc.values();

#pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();
        auto nthreads  = omp_get_max_threads();

        auto grid_batch_size   = mathfunc::batch_size(npoints, thread_id, nthreads);
        auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

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
                    for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                    {
                        W_val[jkl * npoints + g] = weights[g] * vrho[2 * g + 1] * mo_j[g] * mo_k[g] * mo_l[g];
                    }
                }
            }
        }
    }

    auto tensorWxc = denblas::multABt(gtoValues, matrixWxc);

    timer.stop("Wxc matrix");

    return tensorWxc;
}

}  // namespace xcintpdft
