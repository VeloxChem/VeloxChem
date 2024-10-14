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

#include "XCIntegratorForGGA.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <sstream>

#include "DenseLinearAlgebra.hpp"
#include "DenseMatrix.hpp"
#include "DensityGridGenerator.hpp"
#include "DftSubMatrix.hpp"
#include "ErrorHandler.hpp"
#include "GtoFunc.hpp"
#include "GtoValues.hpp"
#include "MathFunc.hpp"
#include "MultiTimer.hpp"
#include "Prescreener.hpp"
#include "StringFormat.hpp"

namespace xcintgga {  // xcintgga namespace

auto
integrateVxcFockForGGA(const CMolecule&                  molecule,
                       const CMolecularBasis&            basis,
                       const std::vector<const double*>& gsDensityPointers,
                       const CMolecularGrid&             molecularGrid,
                       const double                      screeningThresholdForGTOValues,
                       const CXCFunctional&              xcFunctional,
                       const std::string&                flag) -> CAOKohnShamMatrix
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // Kohn-Sham matrix

    bool closedshell = (format::upper_case(flag) == std::string("CLOSEDSHELL"));

    CAOKohnShamMatrix mat_Vxc(naos, naos, closedshell);

    mat_Vxc.zero();

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    // density and functional derivatives

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhograd_data(dim->rho * 3 * max_npoints_per_box);
    std::vector<double> sigma_data(dim->sigma * max_npoints_per_box);

    std::vector<double> exc_data(dim->zk * max_npoints_per_box);
    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> vsigma_data(dim->vsigma * max_npoints_per_box);

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

    for (size_t box_id = 0; box_id < counts.size(); box_id++)
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

        // generate sub density matrix and density grid

        if (closedshell)
        {
            timer.start("Density matrix slicing");

            auto sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

            timer.stop("Density matrix slicing");

            dengridgen::generateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat, timer);
        }
        else
        {
            timer.start("Density matrix slicing");

            auto sub_dens_mat_a = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);
            auto sub_dens_mat_b = dftsubmat::getSubDensityMatrix(gsDensityPointers[1], aoinds, naos);

            dengridgen::generateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat_a, sub_dens_mat_b, timer);

            timer.stop("Density matrix slicing");
        }

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_gga(npoints, rho, sigma, exc, vrho, vsigma);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        if (closedshell)
        {
            auto partial_mat_Vxc =
                integratePartialVxcFockForGGA(local_weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rhograd, vrho, vsigma, timer);

            timer.start("Vxc matrix dist.");

            dftsubmat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc, aoinds);

            timer.stop("Vxc matrix dist.");
        }
        else
        {
            auto partial_mat_Vxc_ab =
                integratePartialVxcFockForGGAOpenShell(local_weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rhograd, vrho, vsigma, timer);

            timer.start("Vxc matrix dist.");

            dftsubmat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc_ab, aoinds);

            timer.stop("Vxc matrix dist.");
        }

        // compute partial contribution to XC energy

        timer.start("XC energy");

        for (int g = 0; g < npoints; g++)
        {
            auto rho_total = rho[2 * g + 0] + rho[2 * g + 1];

            nele += local_weights[g] * rho_total;

            xcene += local_weights[g] * exc[g] * rho_total;
        }

        timer.stop("XC energy");
    }

    mat_Vxc.setNumberOfElectrons(nele);

    mat_Vxc.setExchangeCorrelationEnergy(xcene);

    timer.stop("Total timing");

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

#pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

    auto mat_Vxc = denblas::multABt(gtoValues, denblas::addAB(mat_G, mat_G_gga, 2.0));

    mat_Vxc.symmetrizeAndScale(0.5);

    timer.stop("Vxc matrix matmul");

    return mat_Vxc;
}

auto
integratePartialVxcFockForGGAOpenShell(const double*       weights,
                                       const CDenseMatrix& gtoValues,
                                       const CDenseMatrix& gtoValuesX,
                                       const CDenseMatrix& gtoValuesY,
                                       const CDenseMatrix& gtoValuesZ,
                                       const double*       rhograd,
                                       const double*       vrho,
                                       const double*       vsigma,
                                       CMultiTimer&        timer) -> std::vector<CDenseMatrix>
{
    const auto npoints = gtoValues.getNumberOfColumns();

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G_a(naos, npoints);
    CDenseMatrix mat_G_b(naos, npoints);

    CDenseMatrix mat_G_a_gga(naos, npoints);
    CDenseMatrix mat_G_b_gga(naos, npoints);

    auto G_a_val = mat_G_a.values();
    auto G_b_val = mat_G_b.values();

    auto G_a_gga_val = mat_G_a_gga.values();
    auto G_b_gga_val = mat_G_b_gga.values();

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

#pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                auto vxa = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                auto vya = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                auto vza = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                auto vxb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 3] + vsigma[3 * g + 1] * rhograd[6 * g + 0];
                auto vyb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 4] + vsigma[3 * g + 1] * rhograd[6 * g + 1];
                auto vzb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 5] + vsigma[3 * g + 1] * rhograd[6 * g + 2];

                G_a_val[nu_offset + g] = weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];
                G_b_val[nu_offset + g] = weights[g] * vrho[2 * g + 1] * chi_val[nu_offset + g];

                G_a_gga_val[nu_offset + g] =
                    weights[g] * (vxa * chi_x_val[nu_offset + g] + vya * chi_y_val[nu_offset + g] + vza * chi_z_val[nu_offset + g]);
                G_b_gga_val[nu_offset + g] =
                    weights[g] * (vxb * chi_x_val[nu_offset + g] + vyb * chi_y_val[nu_offset + g] + vzb * chi_z_val[nu_offset + g]);
            }
        }
    }

    timer.stop("Vxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    // Note that we use matrix-matrix multiplication only once, and symmetrize
    // the result. This is because the density matrix is symmetric, and the
    // Kohn-Sham matrix from mat_G is also symmetric. Formally only the
    // mat_G_gga contribution should be symmetrized.

    timer.start("Vxc matrix matmul");

    auto mat_Vxc_a = denblas::multABt(gtoValues, denblas::addAB(mat_G_a, mat_G_a_gga, 2.0));
    auto mat_Vxc_b = denblas::multABt(gtoValues, denblas::addAB(mat_G_b, mat_G_b_gga, 2.0));

    mat_Vxc_a.symmetrizeAndScale(0.5);
    mat_Vxc_b.symmetrizeAndScale(0.5);

    timer.stop("Vxc matrix matmul");

    return std::vector<CDenseMatrix>{mat_Vxc_a, mat_Vxc_b};
}

auto
integrateFxcFockForGGA(const std::vector<double*>&       aoFockPointers,
                       const CMolecule&                  molecule,
                       const CMolecularBasis&            basis,
                       const std::vector<const double*>& rwDensityPointers,
                       const std::vector<const double*>& gsDensityPointers,
                       const CMolecularGrid&             molecularGrid,
                       const double                      screeningThresholdForGTOValues,
                       const CXCFunctional&              xcFunctional) -> void
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    // density and functional derivatives

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhograd_data(dim->rho * 3 * max_npoints_per_box);
    std::vector<double> sigma_data(dim->sigma * max_npoints_per_box);

    std::vector<double> rhow_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhowgrad_data(dim->rho * 3 * max_npoints_per_box);

    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> vsigma_data(dim->vsigma * max_npoints_per_box);

    std::vector<double> v2rho2_data(dim->v2rho2 * max_npoints_per_box);
    std::vector<double> v2rhosigma_data(dim->v2rhosigma * max_npoints_per_box);
    std::vector<double> v2sigma2_data(dim->v2sigma2 * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();

    auto rhow     = rhow_data.data();
    auto rhowgrad = rhowgrad_data.data();

    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();

    auto v2rho2     = v2rho2_data.data();
    auto v2rhosigma = v2rhosigma_data.data();
    auto v2sigma2   = v2sigma2_data.data();

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

        // pre-screening

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

        // generate sub density matrix and density grid

        timer.start("Density matrix slicing");

        auto sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat, timer);

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

        xcFunctional.compute_fxc_for_gga(npoints, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // go through rhow density matrices

        for (size_t idensity = 0; idensity < rwDensityPointers.size(); idensity++)
        {
            // generate sub density matrix

            timer.start("Density matrix slicing");

            auto sub_dens_mat = dftsubmat::getSubDensityMatrix(rwDensityPointers[idensity], aoinds, naos);

            timer.stop("Density matrix slicing");

            // generate density grid

            dengridgen::generateDensityForGGA(rhow, rhowgrad, nullptr, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat, timer);

            // compute partial contribution to Fxc matrix

            auto partial_mat_Fxc = integratePartialFxcFockForGGA(xcFunctional,
                                                                 local_weights,
                                                                 mat_chi,
                                                                 mat_chi_x,
                                                                 mat_chi_y,
                                                                 mat_chi_z,
                                                                 rhow,
                                                                 rhograd,
                                                                 rhowgrad,
                                                                 vsigma,
                                                                 v2rho2,
                                                                 v2rhosigma,
                                                                 v2sigma2,
                                                                 timer);

            // distribute partial Fxc to full Fock matrix

            timer.start("Fxc matrix dist.");

            dftsubmat::distributeSubMatrixToFock(aoFockPointers, idensity, partial_mat_Fxc, aoinds, naos);

            timer.stop("Fxc matrix dist.");
        }
    }

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
integratePartialFxcFockForGGA(const CXCFunctional& xcFunctional,
                              const double*        weights,
                              const CDenseMatrix&  gtoValues,
                              const CDenseMatrix&  gtoValuesX,
                              const CDenseMatrix&  gtoValuesY,
                              const CDenseMatrix&  gtoValuesZ,
                              const double*        rhow,
                              const double*        rhograd,
                              const double*        rhowgrad,
                              const double*        vsigma,
                              const double*        v2rho2,
                              const double*        v2rhosigma,
                              const double*        v2sigma2,
                              CMultiTimer&         timer) -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Fxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);
    CDenseMatrix mat_G_gga(naos, npoints);

    auto G_val     = mat_G.values();
    auto G_gga_val = mat_G_gga.values();

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

#pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                double w = weights[g];

                auto rhow_a = rhow[2 * g + 0];
                auto rhow_b = rhow[2 * g + 1];

                auto rhowgrad_ax = rhowgrad[6 * g + 0];
                auto rhowgrad_ay = rhowgrad[6 * g + 1];
                auto rhowgrad_az = rhowgrad[6 * g + 2];
                auto rhowgrad_bx = rhowgrad[6 * g + 3];
                auto rhowgrad_by = rhowgrad[6 * g + 4];
                auto rhowgrad_bz = rhowgrad[6 * g + 5];

                auto rhograd_ax = rhograd[6 * g + 0];
                auto rhograd_ay = rhograd[6 * g + 1];
                auto rhograd_az = rhograd[6 * g + 2];
                auto rhograd_bx = rhograd[6 * g + 3];
                auto rhograd_by = rhograd[6 * g + 4];
                auto rhograd_bz = rhograd[6 * g + 5];

                auto grhow_grho_aa = 2.0 * (rhowgrad_ax * rhograd_ax + rhowgrad_ay * rhograd_ay + rhowgrad_az * rhograd_az);

                auto grhow_grho_bb = 2.0 * (rhowgrad_bx * rhograd_bx + rhowgrad_by * rhograd_by + rhowgrad_bz * rhograd_bz);

                auto grhow_grho_ab = (rhowgrad_ax * rhograd_bx + rhowgrad_ay * rhograd_by + rhowgrad_az * rhograd_bz +

                                      rhowgrad_bx * rhograd_ax + rhowgrad_by * rhograd_ay + rhowgrad_bz * rhograd_az);

                // functional derivatives

                // first-order

                auto vsigma_a = vsigma[dim->vsigma * g + 0];
                auto vsigma_c = vsigma[dim->vsigma * g + 1];

                // second-order

                auto v2rho2_aa = v2rho2[dim->v2rho2 * g + 0];
                auto v2rho2_ab = v2rho2[dim->v2rho2 * g + 1];

                auto v2rhosigma_aa = v2rhosigma[dim->v2rhosigma * g + 0];
                auto v2rhosigma_ac = v2rhosigma[dim->v2rhosigma * g + 1];
                auto v2rhosigma_ab = v2rhosigma[dim->v2rhosigma * g + 2];
                auto v2rhosigma_ba = v2rhosigma[dim->v2rhosigma * g + 3];
                auto v2rhosigma_bc = v2rhosigma[dim->v2rhosigma * g + 4];

                auto v2sigma2_aa = v2sigma2[dim->v2sigma2 * g + 0];
                auto v2sigma2_ac = v2sigma2[dim->v2sigma2 * g + 1];
                auto v2sigma2_ab = v2sigma2[dim->v2sigma2 * g + 2];
                auto v2sigma2_cc = v2sigma2[dim->v2sigma2 * g + 3];
                auto v2sigma2_cb = v2sigma2[dim->v2sigma2 * g + 4];

                // scalar contribution

                double f_0 = v2rho2_aa * rhow_a + v2rho2_ab * rhow_b + v2rhosigma_aa * grhow_grho_aa + v2rhosigma_ac * grhow_grho_ab +
                             v2rhosigma_ab * grhow_grho_bb;

                G_val[nu_offset + g] = w * f_0 * chi_val[nu_offset + g];

                // vector contribution

                double f_aa = v2rhosigma_aa * rhow_a + v2rhosigma_ba * rhow_b + v2sigma2_aa * grhow_grho_aa + v2sigma2_ac * grhow_grho_ab +
                              v2sigma2_ab * grhow_grho_bb;

                double f_ab = v2rhosigma_ac * rhow_a + v2rhosigma_bc * rhow_b + v2sigma2_ac * grhow_grho_aa + v2sigma2_cc * grhow_grho_ab +
                              v2sigma2_cb * grhow_grho_bb;

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                xcomp += 2.0 * f_aa * rhograd_ax + f_ab * rhograd_bx;
                ycomp += 2.0 * f_aa * rhograd_ay + f_ab * rhograd_by;
                zcomp += 2.0 * f_aa * rhograd_az + f_ab * rhograd_bz;

                xcomp += 2.0 * vsigma_a * rhowgrad_ax + vsigma_c * rhowgrad_bx;
                ycomp += 2.0 * vsigma_a * rhowgrad_ay + vsigma_c * rhowgrad_by;
                zcomp += 2.0 * vsigma_a * rhowgrad_az + vsigma_c * rhowgrad_bz;

                G_gga_val[nu_offset + g] =
                    w * (xcomp * chi_x_val[nu_offset + g] + ycomp * chi_y_val[nu_offset + g] + zcomp * chi_z_val[nu_offset + g]);
            }
        }
    }

    timer.stop("Fxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Fxc matrix matmul");

    auto mat_Fxc = denblas::multABt(gtoValues, mat_G);

    auto mat_Fxc_gga = denblas::multABt(gtoValues, mat_G_gga);

    mat_Fxc_gga.symmetrize();  // matrix + matrix.T

    mat_Fxc = denblas::addAB(mat_Fxc, mat_Fxc_gga, 1.0);

    timer.stop("Fxc matrix matmul");

    return mat_Fxc;
}

auto
integrateKxcFockForGGA(const std::vector<double*>& aoFockPointers,
                       const CMolecule&            molecule,
                       const CMolecularBasis&      basis,
                       const CAODensityMatrix&     rwDensityMatrix,
                       const CAODensityMatrix&     rw2DensityMatrix,
                       const CAODensityMatrix&     gsDensityMatrix,
                       const CMolecularGrid&       molecularGrid,
                       const double                screeningThresholdForGTOValues,
                       const CXCFunctional&        xcFunctional,
                       const std::string&          quadMode) -> void
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    // density and functional derivatives

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhograd_data(dim->rho * 3 * max_npoints_per_box);
    std::vector<double> sigma_data(dim->sigma * max_npoints_per_box);

    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> vsigma_data(dim->vsigma * max_npoints_per_box);

    std::vector<double> v2rho2_data(dim->v2rho2 * max_npoints_per_box);
    std::vector<double> v2rhosigma_data(dim->v2rhosigma * max_npoints_per_box);
    std::vector<double> v2sigma2_data(dim->v2sigma2 * max_npoints_per_box);

    std::vector<double> v3rho3_data(dim->v3rho3 * max_npoints_per_box);
    std::vector<double> v3rho2sigma_data(dim->v3rho2sigma * max_npoints_per_box);
    std::vector<double> v3rhosigma2_data(dim->v3rhosigma2 * max_npoints_per_box);
    std::vector<double> v3sigma3_data(dim->v3sigma3 * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();

    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();

    auto v2rho2     = v2rho2_data.data();
    auto v2rhosigma = v2rhosigma_data.data();
    auto v2sigma2   = v2sigma2_data.data();

    auto v3rho3      = v3rho3_data.data();
    auto v3rho2sigma = v3rho2sigma_data.data();
    auto v3rhosigma2 = v3rhosigma2_data.data();
    auto v3sigma3    = v3sigma3_data.data();

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

        // pre-screening

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

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityMatrix.alphaDensity(0), aoinds, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat, timer);

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto rw_sub_dens_mat = dftsubmat::getSubAODensityMatrix(rwDensityMatrix, aoinds);

        auto rw2_sub_dens_mat = dftsubmat::getSubAODensityMatrix(rw2DensityMatrix, aoinds);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto rwdengrid = dengridgen::generateDensityGridForGGA(mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw_sub_dens_mat, xcfuntype, timer);

        auto rw2dengrid =
            dengridgen::generateDensityGridForGGA(mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw2_sub_dens_mat, xcfuntype, timer);

        // compute perturbed density

        timer.start("Density grid quad");

        auto numdens_rw2 = rw2DensityMatrix.getNumberOfDensityMatrices();

        CDensityGridQuad rwdengridquad(npoints, numdens_rw2, xcfuntype, dengrid::ab);

        rwdengridquad.DensityProd(rwdengrid, xcfuntype, numdens_rw2, quadMode);

        timer.stop("Density grid quad");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

        xcFunctional.compute_fxc_for_gga(npoints, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

        xcFunctional.compute_kxc_for_gga(npoints, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // go through density matrices

        for (int idensity = 0; idensity < numdens_rw2; idensity++)
        {
            // compute partial contribution to Kxc matrix

            auto partial_mat_Kxc = integratePartialKxcFockForGGA(xcFunctional,
                                                                 local_weights,
                                                                 mat_chi,
                                                                 mat_chi_x,
                                                                 mat_chi_y,
                                                                 mat_chi_z,
                                                                 rhograd,
                                                                 vsigma,
                                                                 v2rho2,
                                                                 v2rhosigma,
                                                                 v2sigma2,
                                                                 v3rho3,
                                                                 v3rho2sigma,
                                                                 v3rhosigma2,
                                                                 v3sigma3,
                                                                 rwdengridquad,
                                                                 rw2dengrid,
                                                                 idensity,
                                                                 timer);

            // distribute partial Kxc to full Fock matrix

            timer.start("Kxc matrix dist.");

            dftsubmat::distributeSubMatrixToFock(aoFockPointers, idensity, partial_mat_Kxc, aoinds, naos);

            timer.stop("Kxc matrix dist.");
        }
    }

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
integratePartialKxcFockForGGA(const CXCFunctional&    xcFunctional,
                              const double*           weights,
                              const CDenseMatrix&     gtoValues,
                              const CDenseMatrix&     gtoValuesX,
                              const CDenseMatrix&     gtoValuesY,
                              const CDenseMatrix&     gtoValuesZ,
                              const double*           rhograd,
                              const double*           vsigma,
                              const double*           v2rho2,
                              const double*           v2rhosigma,
                              const double*           v2sigma2,
                              const double*           v3rho3,
                              const double*           v3rho2sigma,
                              const double*           v3rhosigma2,
                              const double*           v3sigma3,
                              const CDensityGridQuad& rwDensityGridQuad,
                              const CDensityGrid&     rw2DensityGrid,
                              const int           iFock,
                              CMultiTimer&            timer) -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    timer.start("Kxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    // pointers to perturbed density gradient norms

    auto gam = rwDensityGridQuad.gam(iFock);

    auto gamx = rwDensityGridQuad.gamX(iFock);
    auto gamy = rwDensityGridQuad.gamY(iFock);
    auto gamz = rwDensityGridQuad.gamZ(iFock);

    auto gamxx = rwDensityGridQuad.gamXX(iFock);
    auto gamxy = rwDensityGridQuad.gamXY(iFock);
    auto gamxz = rwDensityGridQuad.gamXZ(iFock);

    auto gamyx = rwDensityGridQuad.gamYX(iFock);
    auto gamyy = rwDensityGridQuad.gamYY(iFock);
    auto gamyz = rwDensityGridQuad.gamYZ(iFock);

    auto gamzx = rwDensityGridQuad.gamZX(iFock);
    auto gamzy = rwDensityGridQuad.gamZY(iFock);
    auto gamzz = rwDensityGridQuad.gamZZ(iFock);

    auto rhow12a = rw2DensityGrid.alphaDensity(iFock);

    auto gradw12a_x = rw2DensityGrid.alphaDensityGradientX(iFock);
    auto gradw12a_y = rw2DensityGrid.alphaDensityGradientY(iFock);
    auto gradw12a_z = rw2DensityGrid.alphaDensityGradientZ(iFock);

    timer.stop("Kxc matrix prep.");

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Kxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);
    CDenseMatrix mat_G_gga(naos, npoints);

    auto G_val     = mat_G.values();
    auto G_gga_val = mat_G_gga.values();

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd 
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                double w = weights[g];

                double rxw12a = gradw12a_x[g];
                double ryw12a = gradw12a_y[g];
                double rzw12a = gradw12a_z[g];

                double grada_x_g = rhograd[6 * g + 0];
                double grada_y_g = rhograd[6 * g + 1];
                double grada_z_g = rhograd[6 * g + 2];

                double l2contract   = grada_x_g * rxw12a + grada_y_g * ryw12a + grada_z_g * rzw12a;
                double l5contract_x = grada_x_g * l2contract;
                double l5contract_y = grada_y_g * l2contract;
                double l5contract_z = grada_z_g * l2contract;
                double q2contract   = grada_x_g * gamx[g] + grada_y_g * gamy[g] + grada_z_g * gamz[g];
                double q3contract   = grada_x_g * grada_x_g * gamxx[g] + grada_x_g * grada_y_g * gamxy[g] + grada_x_g * grada_z_g * gamxz[g] +
                                    grada_y_g * grada_x_g * gamyx[g] + grada_y_g * grada_y_g * gamyy[g] + grada_y_g * grada_z_g * gamyz[g] +
                                    grada_z_g * grada_x_g * gamzx[g] + grada_z_g * grada_y_g * gamzy[g] + grada_z_g * grada_z_g * gamzz[g];

                double q4contract   = gamxx[g] + gamyy[g] + gamzz[g];
                double q7contract_x = grada_x_g * (grada_x_g * gamx[g] + grada_y_g * gamy[g] + grada_z_g * gamz[g]);
                double q7contract_y = grada_y_g * (grada_x_g * gamx[g] + grada_y_g * gamy[g] + grada_z_g * gamz[g]);
                double q7contract_z = grada_z_g * (grada_x_g * gamx[g] + grada_y_g * gamy[g] + grada_z_g * gamz[g]);

                double q8contract_x = grada_x_g * gamxx[g] + grada_y_g * gamxy[g] + grada_z_g * gamxz[g];
                double q8contract_y = grada_x_g * gamyx[g] + grada_y_g * gamyy[g] + grada_z_g * gamyz[g];
                double q8contract_z = grada_x_g * gamzx[g] + grada_y_g * gamzy[g] + grada_z_g * gamzz[g];

                double q9contract_x = grada_x_g * q3contract;
                double q9contract_y = grada_y_g * q3contract;
                double q9contract_z = grada_z_g * q3contract;

                double q10contract_x = grada_x_g * gamxx[g] + grada_y_g * gamyx[g] + grada_z_g * gamzx[g];
                double q10contract_y = grada_x_g * gamxy[g] + grada_y_g * gamyy[g] + grada_z_g * gamzy[g];
                double q10contract_z = grada_x_g * gamxz[g] + grada_y_g * gamyz[g] + grada_z_g * gamzz[g];

                double q11contract_x = grada_x_g * gamxx[g] + grada_x_g * gamyy[g] + grada_x_g * gamzz[g];
                double q11contract_y = grada_y_g * gamxx[g] + grada_y_g * gamyy[g] + grada_y_g * gamzz[g];
                double q11contract_z = grada_z_g * gamxx[g] + grada_z_g * gamyy[g] + grada_z_g * gamzz[g];

                // functional derivatives

                // first-order

                auto vsigma_a = vsigma[dim->vsigma * g + 0];
                auto vsigma_c = vsigma[dim->vsigma * g + 1];

                // second-order

                auto v2rho2_aa = v2rho2[dim->v2rho2 * g + 0];
                auto v2rho2_ab = v2rho2[dim->v2rho2 * g + 1];

                auto v2rhosigma_aa = v2rhosigma[dim->v2rhosigma * g + 0];
                auto v2rhosigma_ac = v2rhosigma[dim->v2rhosigma * g + 1];
                auto v2rhosigma_ab = v2rhosigma[dim->v2rhosigma * g + 2];
                auto v2rhosigma_ba = v2rhosigma[dim->v2rhosigma * g + 3];
                auto v2rhosigma_bc = v2rhosigma[dim->v2rhosigma * g + 4];

                auto v2sigma2_aa = v2sigma2[dim->v2sigma2 * g + 0];
                auto v2sigma2_ac = v2sigma2[dim->v2sigma2 * g + 1];
                auto v2sigma2_ab = v2sigma2[dim->v2sigma2 * g + 2];
                auto v2sigma2_cc = v2sigma2[dim->v2sigma2 * g + 3];
                auto v2sigma2_cb = v2sigma2[dim->v2sigma2 * g + 4];

                // third-order

                auto v3rho3_aaa = v3rho3[dim->v3rho3 * g + 0];
                auto v3rho3_aab = v3rho3[dim->v3rho3 * g + 1];
                auto v3rho3_abb = v3rho3[dim->v3rho3 * g + 2];

                auto v3rho2sigma_aaa = v3rho2sigma[dim->v3rho2sigma * g + 0];
                auto v3rho2sigma_aac = v3rho2sigma[dim->v3rho2sigma * g + 1];
                auto v3rho2sigma_aab = v3rho2sigma[dim->v3rho2sigma * g + 2];
                auto v3rho2sigma_aba = v3rho2sigma[dim->v3rho2sigma * g + 3];
                auto v3rho2sigma_abc = v3rho2sigma[dim->v3rho2sigma * g + 4];
                auto v3rho2sigma_abb = v3rho2sigma[dim->v3rho2sigma * g + 5];
                auto v3rho2sigma_bba = v3rho2sigma[dim->v3rho2sigma * g + 6];
                auto v3rho2sigma_bbc = v3rho2sigma[dim->v3rho2sigma * g + 7];

                auto v3rhosigma2_aaa = v3rhosigma2[dim->v3rhosigma2 * g + 0];
                auto v3rhosigma2_aac = v3rhosigma2[dim->v3rhosigma2 * g + 1];
                auto v3rhosigma2_aab = v3rhosigma2[dim->v3rhosigma2 * g + 2];
                auto v3rhosigma2_acc = v3rhosigma2[dim->v3rhosigma2 * g + 3];
                auto v3rhosigma2_acb = v3rhosigma2[dim->v3rhosigma2 * g + 4];
                auto v3rhosigma2_abb = v3rhosigma2[dim->v3rhosigma2 * g + 5];
                auto v3rhosigma2_baa = v3rhosigma2[dim->v3rhosigma2 * g + 6];
                auto v3rhosigma2_bac = v3rhosigma2[dim->v3rhosigma2 * g + 7];
                auto v3rhosigma2_bab = v3rhosigma2[dim->v3rhosigma2 * g + 8];
                auto v3rhosigma2_bcc = v3rhosigma2[dim->v3rhosigma2 * g + 9];
                auto v3rhosigma2_bcb = v3rhosigma2[dim->v3rhosigma2 * g + 10];

                auto v3sigma3_aaa = v3sigma3[dim->v3sigma3 * g + 0];
                auto v3sigma3_aac = v3sigma3[dim->v3sigma3 * g + 1];
                auto v3sigma3_aab = v3sigma3[dim->v3sigma3 * g + 2];
                auto v3sigma3_acc = v3sigma3[dim->v3sigma3 * g + 3];
                auto v3sigma3_acb = v3sigma3[dim->v3sigma3 * g + 4];
                auto v3sigma3_abb = v3sigma3[dim->v3sigma3 * g + 5];
                auto v3sigma3_ccc = v3sigma3[dim->v3sigma3 * g + 6];
                auto v3sigma3_ccb = v3sigma3[dim->v3sigma3 * g + 7];
                auto v3sigma3_cbb = v3sigma3[dim->v3sigma3 * g + 8];

                // functional derivatives
                double rr  = (v2rho2_aa + v2rho2_ab);
                double rrr = (v3rho3_aaa + 2.0 * v3rho3_aab + v3rho3_abb);
                double rx  = (2.0 * v2rhosigma_ac + 2.0 * v2rhosigma_ab + 2.0 * v2rhosigma_aa);
                double rxr = (2.0 * v3rho2sigma_abc + 2.0 * v3rho2sigma_abb + 2.0 * v3rho2sigma_aba + 2.0 * v3rho2sigma_aac + 2.0 * v3rho2sigma_aab +
                              2.0 * v3rho2sigma_aaa);
                double rxx = (4.0 * v3rhosigma2_acc + 8.0 * v3rhosigma2_acb + 4.0 * v3rhosigma2_abb + 8.0 * v3rhosigma2_aac + 8.0 * v3rhosigma2_aab +
                              4.0 * v3rhosigma2_aaa);
                double x   = vsigma_c + 2.0 * vsigma_a;
                double xr  = v2rhosigma_bc + 2.0 * v2rhosigma_ba + v2rhosigma_ac + 2.0 * v2rhosigma_aa;
                double xx  = 2.0 * v2sigma2_cc + 2.0 * v2sigma2_cb + 6.0 * v2sigma2_ac + 4.0 * v2sigma2_ab + 4.0 * v2sigma2_aa;
                double xrr =
                    v3rho2sigma_bbc + 2.0 * v3rho2sigma_bba + 2.0 * v3rho2sigma_abc + 4.0 * v3rho2sigma_aba + v3rho2sigma_aac + 2.0 * v3rho2sigma_aaa;
                double xxr = 2.0 * v3rhosigma2_bcc + 2.0 * v3rhosigma2_bcb + 6.0 * v3rhosigma2_bac + 4.0 * v3rhosigma2_bab + 4.0 * v3rhosigma2_baa +
                             2.0 * v3rhosigma2_acc + 2.0 * v3rhosigma2_acb + 6.0 * v3rhosigma2_aac + 4.0 * v3rhosigma2_aab + 4.0 * v3rhosigma2_aaa;
                double xxx = 4.0 * v3sigma3_ccc + 8.0 * v3sigma3_ccb + 4.0 * v3sigma3_cbb + 16.0 * v3sigma3_acc + 24.0 * v3sigma3_acb +
                             8.0 * v3sigma3_abb + 20.0 * v3sigma3_aac + 16.0 * v3sigma3_aab + 8.0 * v3sigma3_aaa;

                // Scalar contribution

                double prefac = 0.0;

                // vxc 1 contributions

                prefac += rr * rhow12a[g]  // l1
                          + rx * l2contract;

                // vxc 2 contributions

                prefac += rrr * gam[g]  // q1
                          + rxr * q2contract + rxx * q3contract + rx * q4contract;

                G_val[nu_offset + g] = w * prefac * chi_val[nu_offset + g];

                // vector contribution

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                // vxc 1 contributions

                xcomp += xr * grada_x_g * rhow12a[g]  // l3
                         + x * rxw12a                 // l4
                         + xx * l5contract_x;

                ycomp += xr * grada_y_g * rhow12a[g]  // l3
                         + x * ryw12a                 // l4
                         + xx * l5contract_y;

                zcomp += xr * grada_z_g * rhow12a[g]  // l3
                         + x * rzw12a                 // l4
                         + xx * l5contract_z;

                // vxc 2 contributions

                xcomp += xrr * grada_x_g * gam[g]  // q5
                         + xr * gamx[g]            // q6
                         + xxr * q7contract_x + xx * (q8contract_x + q10contract_x + q11contract_x) + xxx * q9contract_x;

                ycomp += xrr * grada_y_g * gam[g]  // q5
                         + xr * gamy[g]            // q6
                         + xxr * q7contract_y + xx * (q8contract_y + q10contract_y + q11contract_y) + xxx * q9contract_y;

                zcomp += xrr * grada_z_g * gam[g]  // q5
                         + xr * gamz[g]            // q6
                         + xxr * q7contract_z + xx * (q8contract_z + q10contract_z + q11contract_z) + xxx * q9contract_z;

                G_gga_val[nu_offset + g] =
                    w * (xcomp * chi_x_val[nu_offset + g] + ycomp * chi_y_val[nu_offset + g] + zcomp * chi_z_val[nu_offset + g]);
            }
        }
    }

    timer.stop("Kxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Kxc matrix matmul");

    auto mat_Kxc = denblas::multABt(gtoValues, mat_G);

    auto mat_Kxc_gga = denblas::multABt(gtoValues, mat_G_gga);

    mat_Kxc_gga.symmetrize();  // matrix + matrix.T

    mat_Kxc = denblas::addAB(mat_Kxc, mat_Kxc_gga, 1.0);

    timer.stop("Kxc matrix matmul");

    return mat_Kxc;
}

}  // namespace xcintgga
