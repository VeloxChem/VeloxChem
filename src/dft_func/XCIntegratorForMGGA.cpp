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

#include "XCIntegratorForMGGA.hpp"

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

namespace xcintmgga {  // xcintmgga namespace

auto
integrateVxcFockForMGGA(const CMolecule&                  molecule,
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

    auto       mggafunc = xcFunctional.getFunctionalPointerToMetaGgaComponent();
    const auto dim      = &(mggafunc->dim);

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhograd_data(dim->rho * 3 * max_npoints_per_box);
    std::vector<double> sigma_data(dim->sigma * max_npoints_per_box);
    std::vector<double> lapl_data(dim->lapl * max_npoints_per_box);
    std::vector<double> tau_data(dim->tau * max_npoints_per_box);

    std::vector<double> exc_data(dim->zk * max_npoints_per_box);
    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> vsigma_data(dim->vsigma * max_npoints_per_box);
    std::vector<double> vlapl_data(dim->vlapl * max_npoints_per_box);
    std::vector<double> vtau_data(dim->vtau * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();
    auto lapl    = lapl_data.data();
    auto tau     = tau_data.data();

    auto exc    = exc_data.data();
    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();
    auto vlapl  = vlapl_data.data();
    auto vtau   = vtau_data.data();

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

            dengridgen::generateDensityForMGGA(rho, rhograd, sigma, lapl, tau, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat, timer);
        }
        else
        {
            timer.start("Density matrix slicing");

            auto sub_dens_mat_a = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);
            auto sub_dens_mat_b = dftsubmat::getSubDensityMatrix(gsDensityPointers[1], aoinds, naos);

            dengridgen::generateDensityForMGGA(
                rho, rhograd, sigma, lapl, tau, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat_a, sub_dens_mat_b, timer);

            timer.stop("Density matrix slicing");
        }

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_mgga(npoints, rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        if (closedshell)
        {
            auto partial_mat_Vxc =
                integratePartialVxcFockForMGGA(local_weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rhograd, vrho, vsigma, vlapl, vtau, timer);

            timer.start("Vxc matrix dist.");

            dftsubmat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc, aoinds);

            timer.stop("Vxc matrix dist.");
        }
        else
        {
            auto partial_mat_Vxc_ab = integratePartialVxcFockForMGGAOpenShell(
                local_weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rhograd, vrho, vsigma, vlapl, vtau, timer);

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
integratePartialVxcFockForMGGA(const double*       weights,
                               const CDenseMatrix& gtoValues,
                               const CDenseMatrix& gtoValuesX,
                               const CDenseMatrix& gtoValuesY,
                               const CDenseMatrix& gtoValuesZ,
                               const double*       rhograd,
                               const double*       vrho,
                               const double*       vsigma,
                               const double*       vlapl,
                               const double*       vtau,
                               CMultiTimer&        timer) -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    timer.start("Vxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    // LDA contribution
    CDenseMatrix mat_G(naos, npoints);

    // GGA contribution
    CDenseMatrix mat_G_gga(naos, npoints);

    // tau contribution
    CDenseMatrix mat_G_gga_x(naos, npoints);
    CDenseMatrix mat_G_gga_y(naos, npoints);
    CDenseMatrix mat_G_gga_z(naos, npoints);

    auto G_val = mat_G.values();

    auto G_gga_val = mat_G_gga.values();

    auto G_gga_x_val = mat_G_gga_x.values();
    auto G_gga_y_val = mat_G_gga_y.values();
    auto G_gga_z_val = mat_G_gga_z.values();

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
                auto vx = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                auto vy = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                auto vz = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                // LDA contribution
                G_val[nu_offset + g] = weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];

                // GGA contribution (will be scaled by 2 later)
                G_gga_val[nu_offset + g] =
                    weights[g] * (vx * chi_x_val[nu_offset + g] + vy * chi_y_val[nu_offset + g] + vz * chi_z_val[nu_offset + g]);

                // tau contribution (will be scaled by 0.5 later)
                G_gga_x_val[nu_offset + g] = weights[g] * vtau[2 * g + 0] * chi_x_val[nu_offset + g];
                G_gga_y_val[nu_offset + g] = weights[g] * vtau[2 * g + 0] * chi_y_val[nu_offset + g];
                G_gga_z_val[nu_offset + g] = weights[g] * vtau[2 * g + 0] * chi_z_val[nu_offset + g];
            }
        }
    }

    timer.stop("Vxc matrix G");

    // Note that we use matrix-matrix multiplication only once, and symmetrize
    // the result. This is because the density matrix is symmetric, and the
    // Kohn-Sham matrix from mat_G is also symmetric. Formally only the
    // mat_G_gga contribution should be symmetrized.

    timer.start("Vxc matrix matmul");

    // LDA and GGA contribution
    auto mat_Vxc = denblas::multABt(gtoValues, denblas::addAB(mat_G, mat_G_gga, 2.0));

    // tau contribution
    auto mat_Vxc_x = denblas::multABt(gtoValuesX, mat_G_gga_x);
    auto mat_Vxc_y = denblas::multABt(gtoValuesY, mat_G_gga_y);
    auto mat_Vxc_z = denblas::multABt(gtoValuesZ, mat_G_gga_z);

    mat_Vxc = denblas::addAB(mat_Vxc, mat_Vxc_x, 0.5);
    mat_Vxc = denblas::addAB(mat_Vxc, mat_Vxc_y, 0.5);
    mat_Vxc = denblas::addAB(mat_Vxc, mat_Vxc_z, 0.5);

    mat_Vxc.symmetrizeAndScale(0.5);

    timer.stop("Vxc matrix matmul");

    return mat_Vxc;
}

auto
integratePartialVxcFockForMGGAOpenShell(const double*       weights,
                                        const CDenseMatrix& gtoValues,
                                        const CDenseMatrix& gtoValuesX,
                                        const CDenseMatrix& gtoValuesY,
                                        const CDenseMatrix& gtoValuesZ,
                                        const double*       rhograd,
                                        const double*       vrho,
                                        const double*       vsigma,
                                        const double*       vlapl,
                                        const double*       vtau,
                                        CMultiTimer&        timer) -> std::vector<CDenseMatrix>
{
    const auto npoints = gtoValues.getNumberOfColumns();

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    // LDA contribution
    CDenseMatrix mat_G_a(naos, npoints);
    CDenseMatrix mat_G_b(naos, npoints);

    // GGA contribution
    CDenseMatrix mat_G_a_gga(naos, npoints);
    CDenseMatrix mat_G_b_gga(naos, npoints);

    // tau contribution
    CDenseMatrix mat_G_a_gga_x(naos, npoints);
    CDenseMatrix mat_G_a_gga_y(naos, npoints);
    CDenseMatrix mat_G_a_gga_z(naos, npoints);

    CDenseMatrix mat_G_b_gga_x(naos, npoints);
    CDenseMatrix mat_G_b_gga_y(naos, npoints);
    CDenseMatrix mat_G_b_gga_z(naos, npoints);

    auto G_a_val = mat_G_a.values();
    auto G_b_val = mat_G_b.values();

    auto G_a_gga_val = mat_G_a_gga.values();
    auto G_b_gga_val = mat_G_b_gga.values();

    auto G_a_gga_x_val = mat_G_a_gga_x.values();
    auto G_a_gga_y_val = mat_G_a_gga_y.values();
    auto G_a_gga_z_val = mat_G_a_gga_z.values();

    auto G_b_gga_x_val = mat_G_b_gga_x.values();
    auto G_b_gga_y_val = mat_G_b_gga_y.values();
    auto G_b_gga_z_val = mat_G_b_gga_z.values();

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

                // LDA contribution
                G_a_val[nu_offset + g] = weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];
                G_b_val[nu_offset + g] = weights[g] * vrho[2 * g + 1] * chi_val[nu_offset + g];

                // GGA contribution (will be scaled by 2 later)
                G_a_gga_val[nu_offset + g] =
                    weights[g] * (vxa * chi_x_val[nu_offset + g] + vya * chi_y_val[nu_offset + g] + vza * chi_z_val[nu_offset + g]);
                G_b_gga_val[nu_offset + g] =
                    weights[g] * (vxb * chi_x_val[nu_offset + g] + vyb * chi_y_val[nu_offset + g] + vzb * chi_z_val[nu_offset + g]);

                // TODO implement Laplacian dependence

                // tau contribution (will be scaled by 0.5 later)
                G_a_gga_x_val[nu_offset + g] = weights[g] * vtau[2 * g + 0] * chi_x_val[nu_offset + g];
                G_a_gga_y_val[nu_offset + g] = weights[g] * vtau[2 * g + 0] * chi_y_val[nu_offset + g];
                G_a_gga_z_val[nu_offset + g] = weights[g] * vtau[2 * g + 0] * chi_z_val[nu_offset + g];

                G_b_gga_x_val[nu_offset + g] = weights[g] * vtau[2 * g + 1] * chi_x_val[nu_offset + g];
                G_b_gga_y_val[nu_offset + g] = weights[g] * vtau[2 * g + 1] * chi_y_val[nu_offset + g];
                G_b_gga_z_val[nu_offset + g] = weights[g] * vtau[2 * g + 1] * chi_z_val[nu_offset + g];
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

    // LDA and GGA contribution
    auto mat_Vxc_a = denblas::multABt(gtoValues, denblas::addAB(mat_G_a, mat_G_a_gga, 2.0));
    auto mat_Vxc_b = denblas::multABt(gtoValues, denblas::addAB(mat_G_b, mat_G_b_gga, 2.0));

    // tau contribution
    auto mat_Vxc_a_x = denblas::multABt(gtoValuesX, mat_G_a_gga_x);
    auto mat_Vxc_a_y = denblas::multABt(gtoValuesY, mat_G_a_gga_y);
    auto mat_Vxc_a_z = denblas::multABt(gtoValuesZ, mat_G_a_gga_z);

    auto mat_Vxc_b_x = denblas::multABt(gtoValuesX, mat_G_b_gga_x);
    auto mat_Vxc_b_y = denblas::multABt(gtoValuesY, mat_G_b_gga_y);
    auto mat_Vxc_b_z = denblas::multABt(gtoValuesZ, mat_G_b_gga_z);

    mat_Vxc_a = denblas::addAB(mat_Vxc_a, mat_Vxc_a_x, 0.5);
    mat_Vxc_a = denblas::addAB(mat_Vxc_a, mat_Vxc_a_y, 0.5);
    mat_Vxc_a = denblas::addAB(mat_Vxc_a, mat_Vxc_a_z, 0.5);

    mat_Vxc_b = denblas::addAB(mat_Vxc_b, mat_Vxc_b_x, 0.5);
    mat_Vxc_b = denblas::addAB(mat_Vxc_b, mat_Vxc_b_y, 0.5);
    mat_Vxc_b = denblas::addAB(mat_Vxc_b, mat_Vxc_b_z, 0.5);

    mat_Vxc_a.symmetrizeAndScale(0.5);
    mat_Vxc_b.symmetrizeAndScale(0.5);

    timer.stop("Vxc matrix matmul");

    return std::vector<CDenseMatrix>{mat_Vxc_a, mat_Vxc_b};
}

auto
integrateFxcFockForMGGA(const std::vector<double*>&       aoFockPointers,
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

    auto       mggafunc = xcFunctional.getFunctionalPointerToMetaGgaComponent();
    const auto dim      = &(mggafunc->dim);

    std::vector<double> local_weights_data(max_npoints_per_box);

    // ground-state
    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhograd_data(dim->rho * 3 * max_npoints_per_box);
    std::vector<double> sigma_data(dim->sigma * max_npoints_per_box);
    std::vector<double> lapl_data(dim->lapl * max_npoints_per_box);
    std::vector<double> tau_data(dim->tau * max_npoints_per_box);

    // perturbed
    std::vector<double> rhow_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhowgrad_data(dim->rho * 3 * max_npoints_per_box);
    std::vector<double> laplw_data(dim->lapl * max_npoints_per_box);
    std::vector<double> tauw_data(dim->tau * max_npoints_per_box);

    // First-order
    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> vsigma_data(dim->vsigma * max_npoints_per_box);
    std::vector<double> vlapl_data(dim->vlapl * max_npoints_per_box);
    std::vector<double> vtau_data(dim->vtau * max_npoints_per_box);

    // Second-order
    std::vector<double> v2rho2_data(dim->v2rho2 * max_npoints_per_box);
    std::vector<double> v2rhosigma_data(dim->v2rhosigma * max_npoints_per_box);
    std::vector<double> v2rholapl_data(dim->v2rholapl * max_npoints_per_box);
    std::vector<double> v2rhotau_data(dim->v2rhotau * max_npoints_per_box);
    std::vector<double> v2sigma2_data(dim->v2sigma2 * max_npoints_per_box);
    std::vector<double> v2sigmalapl_data(dim->v2sigmalapl * max_npoints_per_box);
    std::vector<double> v2sigmatau_data(dim->v2sigmatau * max_npoints_per_box);
    std::vector<double> v2lapl2_data(dim->v2lapl2 * max_npoints_per_box);
    std::vector<double> v2lapltau_data(dim->v2lapltau * max_npoints_per_box);
    std::vector<double> v2tau2_data(dim->v2tau2 * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    // Ground-state
    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();
    auto lapl    = lapl_data.data();
    auto tau     = tau_data.data();

    // Perturbed
    auto rhow     = rhow_data.data();
    auto rhowgrad = rhowgrad_data.data();
    auto laplw    = laplw_data.data();
    auto tauw     = tauw_data.data();

    // First-order
    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();
    auto vlapl  = vlapl_data.data();
    auto vtau   = vtau_data.data();

    // Second-order
    auto v2rho2      = v2rho2_data.data();
    auto v2rhosigma  = v2rhosigma_data.data();
    auto v2rholapl   = v2rholapl_data.data();
    auto v2rhotau    = v2rhotau_data.data();
    auto v2sigma2    = v2sigma2_data.data();
    auto v2sigmalapl = v2sigmalapl_data.data();
    auto v2sigmatau  = v2sigmatau_data.data();
    auto v2lapl2     = v2lapl2_data.data();
    auto v2lapltau   = v2lapltau_data.data();
    auto v2tau2      = v2tau2_data.data();

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

        // pre-screening of GTOs

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

        dengridgen::generateDensityForMGGA(rho, rhograd, sigma, lapl, tau, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat, timer);

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_mgga(npoints, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau);

        xcFunctional.compute_fxc_for_mgga(
            npoints, rho, sigma, lapl, tau, v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2);

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

            dengridgen::generateDensityForMGGA(rhow, rhowgrad, sigma, laplw, tauw, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat, timer);

            // compute partial contribution to Fxc matrix

            auto partial_mat_Fxc = integratePartialFxcFockForMGGA(xcFunctional,
                                                                  local_weights,
                                                                  mat_chi,
                                                                  mat_chi_x,
                                                                  mat_chi_y,
                                                                  mat_chi_z,
                                                                  rhow,
                                                                  rhograd,
                                                                  rhowgrad,
                                                                  tauw,
                                                                  laplw,
                                                                  vrho,
                                                                  vsigma,
                                                                  vlapl,
                                                                  vtau,
                                                                  v2rho2,
                                                                  v2lapl2,
                                                                  v2tau2,
                                                                  v2rholapl,
                                                                  v2rhotau,
                                                                  v2lapltau,
                                                                  v2rhosigma,
                                                                  v2sigmalapl,
                                                                  v2sigmatau,
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
integratePartialFxcFockForMGGA(const CXCFunctional& xcFunctional,
                               const double*        weights,
                               const CDenseMatrix&  gtoValues,
                               const CDenseMatrix&  gtoValuesX,
                               const CDenseMatrix&  gtoValuesY,
                               const CDenseMatrix&  gtoValuesZ,
                               const double*        rhow,
                               const double*        rhograd,
                               const double*        rhowgrad,
                               const double*        tauw,
                               const double*        laplw,
                               const double*        vrho,
                               const double*        vsigma,
                               const double*        vlapl,
                               const double*        vtau,
                               const double*        v2rho2,
                               const double*        v2lapl2,
                               const double*        v2tau2,
                               const double*        v2rholapl,
                               const double*        v2rhotau,
                               const double*        v2lapltau,
                               const double*        v2rhosigma,
                               const double*        v2sigmalapl,
                               const double*        v2sigmatau,
                               const double*        v2sigma2,
                               CMultiTimer&         timer) -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Fxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    // LDA contribution
    CDenseMatrix mat_G(naos, npoints);

    // GGA contribution
    CDenseMatrix mat_G_gga(naos, npoints);

    // tau contribution
    CDenseMatrix mat_G_gga_x(naos, npoints);
    CDenseMatrix mat_G_gga_y(naos, npoints);
    CDenseMatrix mat_G_gga_z(naos, npoints);

    auto G_val = mat_G.values();

    auto G_gga_val = mat_G_gga.values();

    auto G_gga_x_val = mat_G_gga_x.values();
    auto G_gga_y_val = mat_G_gga_y.values();
    auto G_gga_z_val = mat_G_gga_z.values();

    auto chi_val   = gtoValues.values();
    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    auto       mggafunc = xcFunctional.getFunctionalPointerToMetaGgaComponent();
    const auto dim      = &(mggafunc->dim);

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

                // ground-state gardients
                double grada_x_g = rhograd[6 * g + 0];
                double grada_y_g = rhograd[6 * g + 1];
                double grada_z_g = rhograd[6 * g + 2];

                // perturbed density and its gradients
                double rwa   = rhow[2 * g + 0];
                double tauwa = tauw[2 * g + 0];
                // double laplwa = laplw[2 * g + 0];

                double rwa_x = rhowgrad[6 * g + 0];
                double rwa_y = rhowgrad[6 * g + 1];
                double rwa_z = rhowgrad[6 * g + 2];

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

                // auto v2rholapl_aa = v2rholapl[dim->v2rholapl * g + 0];
                // auto v2rholapl_ab = v2rholapl[dim->v2rholapl * g + 1];

                auto v2rhotau_aa = v2rhotau[dim->v2rhotau * g + 0];
                auto v2rhotau_ab = v2rhotau[dim->v2rhotau * g + 1];
                auto v2rhotau_ba = v2rhotau[dim->v2rhotau * g + 2];

                auto v2sigma2_aa = v2sigma2[dim->v2sigma2 * g + 0];
                auto v2sigma2_ac = v2sigma2[dim->v2sigma2 * g + 1];
                auto v2sigma2_ab = v2sigma2[dim->v2sigma2 * g + 2];
                auto v2sigma2_cc = v2sigma2[dim->v2sigma2 * g + 3];
                auto v2sigma2_cb = v2sigma2[dim->v2sigma2 * g + 4];

                // auto v2sigmalapl_aa = v2sigmalapl[dim->v2sigmalapl * g + 0];
                // auto v2sigmalapl_ab = v2sigmalapl[dim->v2sigmalapl * g + 1];
                // auto v2sigmalapl_ca = v2sigmalapl[dim->v2sigmalapl * g + 2];
                // auto v2sigmalapl_cb = v2sigmalapl[dim->v2sigmalapl * g + 3];

                auto v2sigmatau_aa = v2sigmatau[dim->v2sigmatau * g + 0];
                auto v2sigmatau_ab = v2sigmatau[dim->v2sigmatau * g + 1];
                auto v2sigmatau_ca = v2sigmatau[dim->v2sigmatau * g + 2];
                auto v2sigmatau_cb = v2sigmatau[dim->v2sigmatau * g + 3];
                auto v2sigmatau_ba = v2sigmatau[dim->v2sigmatau * g + 4];

                // auto v2lapl2_aa = v2lapl2[dim->v2lapl2 * g + 0];
                // auto v2lapl2_ab = v2lapl2[dim->v2lapl2 * g + 1];

                // auto v2lapltau_aa = v2lapltau[dim->v2lapltau * g + 0];
                // auto v2lapltau_ba = v2lapltau[dim->v2lapltau * g + 2];

                auto v2tau2_aa = v2tau2[dim->v2tau2 * g + 0];
                auto v2tau2_ab = v2tau2[dim->v2tau2 * g + 1];

                // sums of functional derivatives that can be used in the restricted case

                // first-order
                double x = vsigma_c + 2.0 * vsigma_a;
                // second-order
                double rr = v2rho2_aa + v2rho2_ab;
                double rx = 2.0 * v2rhosigma_ac + 2.0 * v2rhosigma_ab + 2.0 * v2rhosigma_aa;
                double rt = v2rhotau_aa + v2rhotau_ab;
                // double rl = v2rholapl_aa + v2rholapl_ab;

                // sigma and gamma
                double xr = v2rhosigma_bc + 2.0 * v2rhosigma_ba + v2rhosigma_ac + 2.0 * v2rhosigma_aa;
                double xt = v2sigmatau_cb + 2.0 * v2sigmatau_ab + v2sigmatau_ca + 2.0 * v2sigmatau_aa;
                // double xl = v2sigmalapl_cb + 2.0*v2sigmalapl_ab + v2sigmalapl_ca + 2.0 * v2sigmalapl_aa;
                double xx = 2.0 * v2sigma2_cc + 2.0 * v2sigma2_cb + 6.0 * v2sigma2_ac + 4.0 * v2sigma2_ab + 4.0 * v2sigma2_aa;

                // tau
                double tt = v2tau2_aa + v2tau2_ab;
                double tx = 2.0 * v2sigmatau_ca + 2.0 * v2sigmatau_ba + 2.0 * v2sigmatau_aa;
                double tr = v2rhotau_aa + v2rhotau_ba;
                // double tl = v2lapltau_aa + v2lapltau_ba;

                // lapl
                // double ll = v2lapl2_aa + v2lapl2_ab;
                // double lx = 2.0 * v2sigmalapl_ca + 2.0 * v2sigmalapl_ba + 2.0 * v2sigmalapl_aa;
                // double lr = v2rholapl_aa + v2rholapl_ba;
                // double lt = v2lapltau_aa + v2lapltau_ab;

                // contraction of perturbed density for restricted case
                double contract = grada_x_g * rwa_x + grada_y_g * rwa_y + grada_z_g * rwa_z;

                // rho-operator
                double r_0 = rr * rwa + rx * contract + rt * tauwa;
                //+ rl * laplwa;

                G_val[nu_offset + g] = w * r_0 * chi_val[nu_offset + g];

                // GGA contribution (will be scaled by 2 later)

                // vector contribution

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                // grad-operator

                // xcomp +=  grada_x_g * ( xr * rwa + xt * tauwa + xl * laplwa)
                xcomp += grada_x_g * (xr * rwa + xt * tauwa) + x * rwa_x + xx * grada_x_g * contract;

                // ycomp += grada_y_g * ( xr * rwa + xt * tauwa + xl * laplwa)
                ycomp += grada_y_g * (xr * rwa + xt * tauwa) + x * rwa_y + xx * grada_y_g * contract;

                // zcomp += grada_z_g * ( xr * rwa + xt * tauwa + xl * laplwa)
                zcomp += grada_z_g * (xr * rwa + xt * tauwa) + x * rwa_z + xx * grada_z_g * contract;

                G_gga_val[nu_offset + g] =
                    w * (xcomp * chi_x_val[nu_offset + g] + ycomp * chi_y_val[nu_offset + g] + zcomp * chi_z_val[nu_offset + g]);

                // TODO implement Laplacian dependence

                // double lap_0 =    lr * rwa
                //                 + lx * contract
                //                 + lt * tauwa
                //                 + ll * laplwa;
                // G_gga_val[nu_offset + g] += w * lap_0 * (chi_xx_val[nu_offset + g] +
                //                                          chi_yy_val[nu_offset + g] +
                //                                          chi_zz_val[nu_offset + g]);

                // tau contribution (will be scaled by 0.5 later)
                double tau_0 = tr * rwa + tx * contract + tt * tauwa;
                //+ tl * laplwa;

                G_gga_x_val[nu_offset + g] = w * tau_0 * chi_x_val[nu_offset + g];
                G_gga_y_val[nu_offset + g] = w * tau_0 * chi_y_val[nu_offset + g];
                G_gga_z_val[nu_offset + g] = w * tau_0 * chi_z_val[nu_offset + g];
            }
        }
    }

    timer.stop("Fxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Fxc matrix matmul");

    // LDA and GGA contribution
    auto mat_Fxc = denblas::multABt(gtoValues, mat_G);

    auto mat_Fxc_gga = denblas::multABt(gtoValues, mat_G_gga);

    mat_Fxc_gga.symmetrize();  // (matrix + matrix.T)

    mat_Fxc = denblas::addAB(mat_Fxc, mat_Fxc_gga, 1.0);

    // tau contribution
    auto mat_Fxc_x = denblas::multABt(gtoValuesX, mat_G_gga_x);
    auto mat_Fxc_y = denblas::multABt(gtoValuesY, mat_G_gga_y);
    auto mat_Fxc_z = denblas::multABt(gtoValuesZ, mat_G_gga_z);

    mat_Fxc = denblas::addAB(mat_Fxc, mat_Fxc_x, 0.5);
    mat_Fxc = denblas::addAB(mat_Fxc, mat_Fxc_y, 0.5);
    mat_Fxc = denblas::addAB(mat_Fxc, mat_Fxc_z, 0.5);

    timer.stop("Fxc matrix matmul");

    return mat_Fxc;
}

auto
integrateKxcFockForMGGA(const std::vector<double*>& aoFockPointers,
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

    std::vector<double> local_weights_data(max_npoints_per_box);

    auto       mggafunc = xcFunctional.getFunctionalPointerToMetaGgaComponent();
    const auto dim      = &(mggafunc->dim);

    // Input
    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhograd_data(dim->rho * 3 * max_npoints_per_box);
    std::vector<double> sigma_data(dim->sigma * max_npoints_per_box);
    std::vector<double> lapl_data(dim->lapl * max_npoints_per_box);
    std::vector<double> tau_data(dim->tau * max_npoints_per_box);

    // First-order
    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> vsigma_data(dim->vsigma * max_npoints_per_box);
    std::vector<double> vlapl_data(dim->vlapl * max_npoints_per_box);
    std::vector<double> vtau_data(dim->vtau * max_npoints_per_box);

    // Second-order
    std::vector<double> v2rho2_data(dim->v2rho2 * max_npoints_per_box);
    std::vector<double> v2rhosigma_data(dim->v2rhosigma * max_npoints_per_box);
    std::vector<double> v2rholapl_data(dim->v2rholapl * max_npoints_per_box);
    std::vector<double> v2rhotau_data(dim->v2rhotau * max_npoints_per_box);
    std::vector<double> v2sigma2_data(dim->v2sigma2 * max_npoints_per_box);
    std::vector<double> v2sigmalapl_data(dim->v2sigmalapl * max_npoints_per_box);
    std::vector<double> v2sigmatau_data(dim->v2sigmatau * max_npoints_per_box);
    std::vector<double> v2lapl2_data(dim->v2lapl2 * max_npoints_per_box);
    std::vector<double> v2lapltau_data(dim->v2lapltau * max_npoints_per_box);
    std::vector<double> v2tau2_data(dim->v2tau2 * max_npoints_per_box);

    // Third-order
    std::vector<double> v3rho3_data(dim->v3rho3 * max_npoints_per_box);
    std::vector<double> v3rho2sigma_data(dim->v3rho2sigma * max_npoints_per_box);
    std::vector<double> v3rho2lapl_data(dim->v3rho2lapl * max_npoints_per_box);
    std::vector<double> v3rho2tau_data(dim->v3rho2tau * max_npoints_per_box);
    std::vector<double> v3rhosigma2_data(dim->v3rhosigma2 * max_npoints_per_box);
    std::vector<double> v3rhosigmalapl_data(dim->v3rhosigmalapl * max_npoints_per_box);
    std::vector<double> v3rhosigmatau_data(dim->v3rhosigmatau * max_npoints_per_box);
    std::vector<double> v3rholapl2_data(dim->v3rholapl2 * max_npoints_per_box);
    std::vector<double> v3rholapltau_data(dim->v3rholapltau * max_npoints_per_box);
    std::vector<double> v3rhotau2_data(dim->v3rhotau2 * max_npoints_per_box);
    std::vector<double> v3sigma3_data(dim->v3sigma3 * max_npoints_per_box);
    std::vector<double> v3sigma2lapl_data(dim->v3sigma2lapl * max_npoints_per_box);
    std::vector<double> v3sigma2tau_data(dim->v3sigma2tau * max_npoints_per_box);
    std::vector<double> v3sigmalapl2_data(dim->v3sigmalapl2 * max_npoints_per_box);
    std::vector<double> v3sigmalapltau_data(dim->v3sigmalapltau * max_npoints_per_box);
    std::vector<double> v3sigmatau2_data(dim->v3sigmatau2 * max_npoints_per_box);
    std::vector<double> v3lapl3_data(dim->v3lapl3 * max_npoints_per_box);
    std::vector<double> v3lapl2tau_data(dim->v3lapl2tau * max_npoints_per_box);
    std::vector<double> v3lapltau2_data(dim->v3lapltau2 * max_npoints_per_box);
    std::vector<double> v3tau3_data(dim->v3tau3 * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    // Input
    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();
    auto lapl    = lapl_data.data();
    auto tau     = tau_data.data();

    // First-order
    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();
    auto vlapl  = vlapl_data.data();
    auto vtau   = vtau_data.data();

    // Second-order
    auto v2rho2      = v2rho2_data.data();
    auto v2rhosigma  = v2rhosigma_data.data();
    auto v2rholapl   = v2rholapl_data.data();
    auto v2rhotau    = v2rhotau_data.data();
    auto v2sigma2    = v2sigma2_data.data();
    auto v2sigmalapl = v2sigmalapl_data.data();
    auto v2sigmatau  = v2sigmatau_data.data();
    auto v2lapl2     = v2lapl2_data.data();
    auto v2lapltau   = v2lapltau_data.data();
    auto v2tau2      = v2tau2_data.data();

    // Third-order
    auto v3rho3         = v3rho3_data.data();
    auto v3rho2sigma    = v3rho2sigma_data.data();
    auto v3rho2lapl     = v3rho2lapl_data.data();
    auto v3rho2tau      = v3rho2tau_data.data();
    auto v3rhosigma2    = v3rhosigma2_data.data();
    auto v3rhosigmalapl = v3rhosigmalapl_data.data();
    auto v3rhosigmatau  = v3rhosigmatau_data.data();
    auto v3rholapl2     = v3rholapl2_data.data();
    auto v3rholapltau   = v3rholapltau_data.data();
    auto v3rhotau2      = v3rhotau2_data.data();
    auto v3sigma3       = v3sigma3_data.data();
    auto v3sigma2lapl   = v3sigma2lapl_data.data();
    auto v3sigma2tau    = v3sigma2tau_data.data();
    auto v3sigmalapl2   = v3sigmalapl2_data.data();
    auto v3sigmalapltau = v3sigmalapltau_data.data();
    auto v3sigmatau2    = v3sigmatau2_data.data();
    auto v3lapl3        = v3lapl3_data.data();
    auto v3lapl2tau     = v3lapl2tau_data.data();
    auto v3lapltau2     = v3lapltau2_data.data();
    auto v3tau3         = v3tau3_data.data();

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

        // pre-screening of GTOs

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

        auto sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityMatrix.alphaDensity(0), aoinds, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForMGGA(rho, rhograd, sigma, lapl, tau, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat, timer);

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto rw_sub_dens_mat = dftsubmat::getSubAODensityMatrix(rwDensityMatrix, aoinds);

        auto rw2_sub_dens_mat = dftsubmat::getSubAODensityMatrix(rw2DensityMatrix, aoinds);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto rwdengrid = dengridgen::generateDensityGridForMGGA(mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw_sub_dens_mat, xcfuntype, timer);

        auto rw2dengrid =
            dengridgen::generateDensityGridForMGGA(mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw2_sub_dens_mat, xcfuntype, timer);

        // compute perturbed density

        timer.start("Density grid quad");

        auto numdens_rw2 = rw2DensityMatrix.getNumberOfDensityMatrices();

        CDensityGridQuad rwdengridquad(npoints, numdens_rw2, xcfuntype, dengrid::ab);

        rwdengridquad.DensityProd(rwdengrid, xcfuntype, numdens_rw2, quadMode);

        timer.stop("Density grid quad");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_mgga(npoints, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau);

        xcFunctional.compute_fxc_for_mgga(
            npoints, rho, sigma, lapl, tau, v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2);

        xcFunctional.compute_kxc_for_mgga(npoints,
                                          rho,
                                          sigma,
                                          lapl,
                                          tau,
                                          v3rho3,
                                          v3rho2sigma,
                                          v3rho2lapl,
                                          v3rho2tau,
                                          v3rhosigma2,
                                          v3rhosigmalapl,
                                          v3rhosigmatau,
                                          v3rholapl2,
                                          v3rholapltau,
                                          v3rhotau2,
                                          v3sigma3,
                                          v3sigma2lapl,
                                          v3sigma2tau,
                                          v3sigmalapl2,
                                          v3sigmalapltau,
                                          v3sigmatau2,
                                          v3lapl3,
                                          v3lapl2tau,
                                          v3lapltau2,
                                          v3tau3);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // go through density matrices

        for (int idensity = 0; idensity < numdens_rw2; idensity++)
        {
            // compute partial contribution to Kxc matrix

            auto partial_mat_Kxc = integratePartialKxcFockForMGGA(xcFunctional,
                                                                  local_weights,
                                                                  mat_chi,
                                                                  mat_chi_x,
                                                                  mat_chi_y,
                                                                  mat_chi_z,
                                                                  rhograd,
                                                                  vsigma,
                                                                  v2rho2,
                                                                  v2rhosigma,
                                                                  v2rholapl,
                                                                  v2rhotau,
                                                                  v2sigma2,
                                                                  v2sigmalapl,
                                                                  v2sigmatau,
                                                                  v2lapl2,
                                                                  v2lapltau,
                                                                  v2tau2,
                                                                  v3rho3,
                                                                  v3rho2sigma,
                                                                  v3rho2lapl,
                                                                  v3rho2tau,
                                                                  v3rhosigma2,
                                                                  v3rhosigmalapl,
                                                                  v3rhosigmatau,
                                                                  v3rholapl2,
                                                                  v3rholapltau,
                                                                  v3rhotau2,
                                                                  v3sigma3,
                                                                  v3sigma2lapl,
                                                                  v3sigma2tau,
                                                                  v3sigmalapl2,
                                                                  v3sigmalapltau,
                                                                  v3sigmatau2,
                                                                  v3lapl3,
                                                                  v3lapl2tau,
                                                                  v3lapltau2,
                                                                  v3tau3,
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
integratePartialKxcFockForMGGA(const CXCFunctional&    xcFunctional,
                               const double*           weights,
                               const CDenseMatrix&     gtoValues,
                               const CDenseMatrix&     gtoValuesX,
                               const CDenseMatrix&     gtoValuesY,
                               const CDenseMatrix&     gtoValuesZ,
                               const double*           rhograd,
                               const double*           vsigma,
                               const double*           v2rho2,
                               const double*           v2rhosigma,
                               const double*           v2rholapl,
                               const double*           v2rhotau,
                               const double*           v2sigma2,
                               const double*           v2sigmalapl,
                               const double*           v2sigmatau,
                               const double*           v2lapl2,
                               const double*           v2lapltau,
                               const double*           v2tau2,
                               const double*           v3rho3,
                               const double*           v3rho2sigma,
                               const double*           v3rho2lapl,
                               const double*           v3rho2tau,
                               const double*           v3rhosigma2,
                               const double*           v3rhosigmalapl,
                               const double*           v3rhosigmatau,
                               const double*           v3rholapl2,
                               const double*           v3rholapltau,
                               const double*           v3rhotau2,
                               const double*           v3sigma3,
                               const double*           v3sigma2lapl,
                               const double*           v3sigma2tau,
                               const double*           v3sigmalapl2,
                               const double*           v3sigmalapltau,
                               const double*           v3sigmatau2,
                               const double*           v3lapl3,
                               const double*           v3lapl2tau,
                               const double*           v3lapltau2,
                               const double*           v3tau3,
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

    // pointers to perturbed densities

    auto gam    = rwDensityGridQuad.gam(iFock);
    auto rt_gam = rwDensityGridQuad.rt_gam(iFock);
    auto rl_gam = rwDensityGridQuad.rl_gam(iFock);
    auto tt_gam = rwDensityGridQuad.tt_gam(iFock);
    auto tl_gam = rwDensityGridQuad.tl_gam(iFock);
    auto ll_gam = rwDensityGridQuad.ll_gam(iFock);

    auto gamx    = rwDensityGridQuad.gamX(iFock);
    auto gamy    = rwDensityGridQuad.gamY(iFock);
    auto gamz    = rwDensityGridQuad.gamZ(iFock);
    auto st_gamx = rwDensityGridQuad.st_gamX(iFock);
    auto st_gamy = rwDensityGridQuad.st_gamY(iFock);
    auto st_gamz = rwDensityGridQuad.st_gamZ(iFock);
    auto sl_gamx = rwDensityGridQuad.sl_gamX(iFock);
    auto sl_gamy = rwDensityGridQuad.sl_gamY(iFock);
    auto sl_gamz = rwDensityGridQuad.sl_gamZ(iFock);

    auto gamxx = rwDensityGridQuad.gamXX(iFock);
    auto gamxy = rwDensityGridQuad.gamXY(iFock);
    auto gamxz = rwDensityGridQuad.gamXZ(iFock);
    auto gamyx = rwDensityGridQuad.gamYX(iFock);
    auto gamyy = rwDensityGridQuad.gamYY(iFock);
    auto gamyz = rwDensityGridQuad.gamYZ(iFock);
    auto gamzx = rwDensityGridQuad.gamZX(iFock);
    auto gamzy = rwDensityGridQuad.gamZY(iFock);
    auto gamzz = rwDensityGridQuad.gamZZ(iFock);

    auto rhow12a    = rw2DensityGrid.alphaDensity(iFock);
    auto tauw12a    = rw2DensityGrid.alphaDensitytau(iFock);
    auto laplw12a   = rw2DensityGrid.alphaDensitylapl(iFock);
    auto gradw12a_x = rw2DensityGrid.alphaDensityGradientX(iFock);
    auto gradw12a_y = rw2DensityGrid.alphaDensityGradientY(iFock);
    auto gradw12a_z = rw2DensityGrid.alphaDensityGradientZ(iFock);

    timer.stop("Kxc matrix prep.");

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Kxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    // LDA contribution
    CDenseMatrix mat_G(naos, npoints);

    // GGA contribution
    CDenseMatrix mat_G_gga(naos, npoints);

    // tau contribution
    CDenseMatrix mat_G_gga_x(naos, npoints);
    CDenseMatrix mat_G_gga_y(naos, npoints);
    CDenseMatrix mat_G_gga_z(naos, npoints);

    auto G_val = mat_G.values();

    auto G_gga_val = mat_G_gga.values();

    auto G_gga_x_val = mat_G_gga_x.values();
    auto G_gga_y_val = mat_G_gga_y.values();
    auto G_gga_z_val = mat_G_gga_z.values();

    auto       mggafunc = xcFunctional.getFunctionalPointerToMetaGgaComponent();
    const auto dim      = &(mggafunc->dim);

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

                // auto v2rholapl_aa = v2rholapl[dim->v2rholapl * g + 0];
                // auto v2rholapl_ab = v2rholapl[dim->v2rholapl * g + 1];

                auto v2rhotau_aa = v2rhotau[dim->v2rhotau * g + 0];
                auto v2rhotau_ab = v2rhotau[dim->v2rhotau * g + 1];
                auto v2rhotau_ba = v2rhotau[dim->v2rhotau * g + 2];

                auto v2sigma2_aa = v2sigma2[dim->v2sigma2 * g + 0];
                auto v2sigma2_ac = v2sigma2[dim->v2sigma2 * g + 1];
                auto v2sigma2_ab = v2sigma2[dim->v2sigma2 * g + 2];
                auto v2sigma2_cc = v2sigma2[dim->v2sigma2 * g + 3];
                auto v2sigma2_cb = v2sigma2[dim->v2sigma2 * g + 4];

                // auto v2sigmalapl_aa = v2sigmalapl[dim->v2sigmalapl * g + 0];
                // auto v2sigmalapl_ab = v2sigmalapl[dim->v2sigmalapl * g + 1];
                // auto v2sigmalapl_ca = v2sigmalapl[dim->v2sigmalapl * g + 2];
                // auto v2sigmalapl_cb = v2sigmalapl[dim->v2sigmalapl * g + 3];

                auto v2sigmatau_aa = v2sigmatau[dim->v2sigmatau * g + 0];
                auto v2sigmatau_ab = v2sigmatau[dim->v2sigmatau * g + 1];
                auto v2sigmatau_ca = v2sigmatau[dim->v2sigmatau * g + 2];
                auto v2sigmatau_cb = v2sigmatau[dim->v2sigmatau * g + 3];
                auto v2sigmatau_ba = v2sigmatau[dim->v2sigmatau * g + 4];

                // auto v2lapl2_aa = v2lapl2[dim->v2lapl2 * g + 0];
                // auto v2lapl2_ab = v2lapl2[dim->v2lapl2 * g + 1];

                // auto v2lapltau_aa = v2lapltau[dim->v2lapltau * g + 0];
                // auto v2lapltau_ba = v2lapltau[dim->v2lapltau * g + 2];

                auto v2tau2_aa = v2tau2[dim->v2tau2 * g + 0];
                auto v2tau2_ab = v2tau2[dim->v2tau2 * g + 1];

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

                // auto v3rho2lapl_aaa = v3rho2lapl[dim->v3rho2lapl * g + 0];
                // auto v3rho2lapl_aab = v3rho2lapl[dim->v3rho2lapl * g + 1];
                // auto v3rho2lapl_aba = v3rho2lapl[dim->v3rho2lapl * g + 2];
                // auto v3rho2lapl_abb = v3rho2lapl[dim->v3rho2lapl * g + 3];

                auto v3rho2tau_aaa = v3rho2tau[dim->v3rho2tau * g + 0];
                auto v3rho2tau_aab = v3rho2tau[dim->v3rho2tau * g + 1];
                auto v3rho2tau_aba = v3rho2tau[dim->v3rho2tau * g + 2];
                auto v3rho2tau_abb = v3rho2tau[dim->v3rho2tau * g + 3];
                auto v3rho2tau_bba = v3rho2tau[dim->v3rho2tau * g + 4];

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

                // auto v3rhosigmalapl_aaa = v3rhosigmalapl[dim->v3rhosigmalapl * g + 0];
                // auto v3rhosigmalapl_aab = v3rhosigmalapl[dim->v3rhosigmalapl * g + 1];
                // auto v3rhosigmalapl_aca = v3rhosigmalapl[dim->v3rhosigmalapl * g + 2];
                // auto v3rhosigmalapl_acb = v3rhosigmalapl[dim->v3rhosigmalapl * g + 3];
                // auto v3rhosigmalapl_aba = v3rhosigmalapl[dim->v3rhosigmalapl * g + 4];
                // auto v3rhosigmalapl_abb = v3rhosigmalapl[dim->v3rhosigmalapl * g + 5];
                // auto v3rhosigmalapl_baa = v3rhosigmalapl[dim->v3rhosigmalapl * g + 6];
                // auto v3rhosigmalapl_bab = v3rhosigmalapl[dim->v3rhosigmalapl * g + 7];
                // auto v3rhosigmalapl_bca = v3rhosigmalapl[dim->v3rhosigmalapl * g + 8];
                // auto v3rhosigmalapl_bcb = v3rhosigmalapl[dim->v3rhosigmalapl * g + 9];

                auto v3rhosigmatau_aaa = v3rhosigmatau[dim->v3rhosigmatau * g + 0];
                auto v3rhosigmatau_aab = v3rhosigmatau[dim->v3rhosigmatau * g + 1];
                auto v3rhosigmatau_aca = v3rhosigmatau[dim->v3rhosigmatau * g + 2];
                auto v3rhosigmatau_acb = v3rhosigmatau[dim->v3rhosigmatau * g + 3];
                auto v3rhosigmatau_aba = v3rhosigmatau[dim->v3rhosigmatau * g + 4];
                auto v3rhosigmatau_abb = v3rhosigmatau[dim->v3rhosigmatau * g + 5];
                auto v3rhosigmatau_baa = v3rhosigmatau[dim->v3rhosigmatau * g + 6];
                auto v3rhosigmatau_bab = v3rhosigmatau[dim->v3rhosigmatau * g + 7];
                auto v3rhosigmatau_bca = v3rhosigmatau[dim->v3rhosigmatau * g + 8];
                auto v3rhosigmatau_bcb = v3rhosigmatau[dim->v3rhosigmatau * g + 9];
                auto v3rhosigmatau_bba = v3rhosigmatau[dim->v3rhosigmatau * g + 10];

                // auto v3rholapl2_aaa = v3rholapl2[dim->v3rholapl2 * g + 0];
                // auto v3rholapl2_aab = v3rholapl2[dim->v3rholapl2 * g + 1];
                // auto v3rholapl2_abb = v3rholapl2[dim->v3rholapl2 * g + 2];

                // auto v3rholapltau_aaa = v3rholapltau[dim->v3rholapltau * g + 0];
                // auto v3rholapltau_aab = v3rholapltau[dim->v3rholapltau * g + 1];
                // auto v3rholapltau_aba = v3rholapltau[dim->v3rholapltau * g + 2];
                // auto v3rholapltau_abb = v3rholapltau[dim->v3rholapltau * g + 3];
                // auto v3rholapltau_baa = v3rholapltau[dim->v3rholapltau * g + 4];
                // auto v3rholapltau_bba = v3rholapltau[dim->v3rholapltau * g + 6];

                auto v3rhotau2_aaa = v3rhotau2[dim->v3rhotau2 * g + 0];
                auto v3rhotau2_aab = v3rhotau2[dim->v3rhotau2 * g + 1];
                auto v3rhotau2_abb = v3rhotau2[dim->v3rhotau2 * g + 2];
                auto v3rhotau2_baa = v3rhotau2[dim->v3rhotau2 * g + 3];
                auto v3rhotau2_bab = v3rhotau2[dim->v3rhotau2 * g + 4];

                auto v3sigma3_aaa = v3sigma3[dim->v3sigma3 * g + 0];
                auto v3sigma3_aac = v3sigma3[dim->v3sigma3 * g + 1];
                auto v3sigma3_aab = v3sigma3[dim->v3sigma3 * g + 2];
                auto v3sigma3_acc = v3sigma3[dim->v3sigma3 * g + 3];
                auto v3sigma3_acb = v3sigma3[dim->v3sigma3 * g + 4];
                auto v3sigma3_abb = v3sigma3[dim->v3sigma3 * g + 5];
                auto v3sigma3_ccc = v3sigma3[dim->v3sigma3 * g + 6];
                auto v3sigma3_ccb = v3sigma3[dim->v3sigma3 * g + 7];
                auto v3sigma3_cbb = v3sigma3[dim->v3sigma3 * g + 8];

                // auto v3sigma2lapl_aaa = v3sigma2lapl[dim->v3sigma2lapl * g + 0];
                // auto v3sigma2lapl_aab = v3sigma2lapl[dim->v3sigma2lapl * g + 1];
                // auto v3sigma2lapl_aca = v3sigma2lapl[dim->v3sigma2lapl * g + 2];
                // auto v3sigma2lapl_acb = v3sigma2lapl[dim->v3sigma2lapl * g + 3];
                // auto v3sigma2lapl_aba = v3sigma2lapl[dim->v3sigma2lapl * g + 4];
                // auto v3sigma2lapl_abb = v3sigma2lapl[dim->v3sigma2lapl * g + 5];
                // auto v3sigma2lapl_cca = v3sigma2lapl[dim->v3sigma2lapl * g + 6];
                // auto v3sigma2lapl_ccb = v3sigma2lapl[dim->v3sigma2lapl * g + 7];
                // auto v3sigma2lapl_cba = v3sigma2lapl[dim->v3sigma2lapl * g + 8];
                // auto v3sigma2lapl_cbb = v3sigma2lapl[dim->v3sigma2lapl * g + 9];

                auto v3sigma2tau_aaa = v3sigma2tau[dim->v3sigma2tau * g + 0];
                auto v3sigma2tau_aab = v3sigma2tau[dim->v3sigma2tau * g + 1];
                auto v3sigma2tau_aca = v3sigma2tau[dim->v3sigma2tau * g + 2];
                auto v3sigma2tau_acb = v3sigma2tau[dim->v3sigma2tau * g + 3];
                auto v3sigma2tau_aba = v3sigma2tau[dim->v3sigma2tau * g + 4];
                auto v3sigma2tau_abb = v3sigma2tau[dim->v3sigma2tau * g + 5];
                auto v3sigma2tau_cca = v3sigma2tau[dim->v3sigma2tau * g + 6];
                auto v3sigma2tau_ccb = v3sigma2tau[dim->v3sigma2tau * g + 7];
                auto v3sigma2tau_cba = v3sigma2tau[dim->v3sigma2tau * g + 8];
                auto v3sigma2tau_cbb = v3sigma2tau[dim->v3sigma2tau * g + 9];
                auto v3sigma2tau_bba = v3sigma2tau[dim->v3sigma2tau * g + 10];

                // auto v3sigmalapl2_aaa = v3sigmalapl2[dim->v3sigmalapl2 * g + 0];
                // auto v3sigmalapl2_aab = v3sigmalapl2[dim->v3sigmalapl2 * g + 1];
                // auto v3sigmalapl2_abb = v3sigmalapl2[dim->v3sigmalapl2 * g + 2];
                // auto v3sigmalapl2_caa = v3sigmalapl2[dim->v3sigmalapl2 * g + 3];
                // auto v3sigmalapl2_cab = v3sigmalapl2[dim->v3sigmalapl2 * g + 4];
                // auto v3sigmalapl2_cbb = v3sigmalapl2[dim->v3sigmalapl2 * g + 5];

                // auto v3sigmalapltau_aaa = v3sigmalapltau[dim->v3sigmalapltau * g + 0];
                // auto v3sigmalapltau_aab = v3sigmalapltau[dim->v3sigmalapltau * g + 1];
                // auto v3sigmalapltau_aba = v3sigmalapltau[dim->v3sigmalapltau * g + 2];
                // auto v3sigmalapltau_abb = v3sigmalapltau[dim->v3sigmalapltau * g + 3];
                // auto v3sigmalapltau_caa = v3sigmalapltau[dim->v3sigmalapltau * g + 4];
                // auto v3sigmalapltau_cab = v3sigmalapltau[dim->v3sigmalapltau * g + 5];
                // auto v3sigmalapltau_cba = v3sigmalapltau[dim->v3sigmalapltau * g + 6];
                // auto v3sigmalapltau_cbb = v3sigmalapltau[dim->v3sigmalapltau * g + 7];
                // auto v3sigmalapltau_baa = v3sigmalapltau[dim->v3sigmalapltau * g + 8];
                // auto v3sigmalapltau_bba = v3sigmalapltau[dim->v3sigmalapltau * g + 10];

                auto v3sigmatau2_aaa = v3sigmatau2[dim->v3sigmatau2 * g + 0];
                auto v3sigmatau2_aab = v3sigmatau2[dim->v3sigmatau2 * g + 1];
                auto v3sigmatau2_abb = v3sigmatau2[dim->v3sigmatau2 * g + 2];
                auto v3sigmatau2_caa = v3sigmatau2[dim->v3sigmatau2 * g + 3];
                auto v3sigmatau2_cab = v3sigmatau2[dim->v3sigmatau2 * g + 4];
                auto v3sigmatau2_cbb = v3sigmatau2[dim->v3sigmatau2 * g + 5];
                auto v3sigmatau2_baa = v3sigmatau2[dim->v3sigmatau2 * g + 6];
                auto v3sigmatau2_bab = v3sigmatau2[dim->v3sigmatau2 * g + 7];

                // auto v3lapl3_aaa = v3lapl3[dim->v3lapl3 * g + 0];
                // auto v3lapl3_aab = v3lapl3[dim->v3lapl3 * g + 1];
                // auto v3lapl3_abb = v3lapl3[dim->v3lapl3 * g + 2];

                // auto v3lapl2tau_aaa = v3lapl2tau[dim->v3lapl2tau * g + 0];
                // auto v3lapl2tau_aba = v3lapl2tau[dim->v3lapl2tau * g + 2];
                // auto v3lapl2tau_bba = v3lapl2tau[dim->v3lapl2tau * g + 4];

                // auto v3lapltau2_aaa = v3lapltau2[dim->v3lapltau2 * g + 0];
                // auto v3lapltau2_aab = v3lapltau2[dim->v3lapltau2 * g + 1];
                // auto v3lapltau2_baa = v3lapltau2[dim->v3lapltau2 * g + 3];
                // auto v3lapltau2_bab = v3lapltau2[dim->v3lapltau2 * g + 4];

                auto v3tau3_aaa = v3tau3[dim->v3tau3 * g + 0];
                auto v3tau3_aab = v3tau3[dim->v3tau3 * g + 1];
                auto v3tau3_abb = v3tau3[dim->v3tau3 * g + 2];

                // functional derivatives

                // first-order
                double x = vsigma_c + 2.0 * vsigma_a;

                // second-order
                // rho
                double rr = v2rho2_aa + v2rho2_ab;
                double rx = 2.0 * v2rhosigma_ac + 2.0 * v2rhosigma_ab + 2.0 * v2rhosigma_aa;
                double rt = v2rhotau_aa + v2rhotau_ab;
                // double rl = v2rholapl_aa + v2rholapl_ab;

                // sigma and gamma
                double xr = v2rhosigma_bc + 2.0 * v2rhosigma_ba + v2rhosigma_ac + 2.0 * v2rhosigma_aa;
                double xt = v2sigmatau_cb + 2.0 * v2sigmatau_ab + v2sigmatau_ca + 2.0 * v2sigmatau_aa;
                // double xl = v2sigmalapl_cb + 2.0*v2sigmalapl_ab + v2sigmalapl_ca + 2.0 * v2sigmalapl_aa;
                double xx = 2.0 * v2sigma2_cc + 2.0 * v2sigma2_cb + 6.0 * v2sigma2_ac + 4.0 * v2sigma2_ab + 4.0 * v2sigma2_aa;

                // tau
                double tt = v2tau2_aa + v2tau2_ab;
                double tx = 2.0 * v2sigmatau_ca + 2.0 * v2sigmatau_ba + 2.0 * v2sigmatau_aa;
                double tr = v2rhotau_aa + v2rhotau_ba;
                // double tl = v2lapltau_aa + v2lapltau_ba;

                // lapl
                // double ll = v2lapl2_aa + v2lapl2_ab;
                // double lx = 2.0 * v2sigmalapl_ca + 2.0 * v2sigmalapl_ba + 2.0 * v2sigmalapl_aa;
                // double lr = v2rholapl_aa + v2rholapl_ba;
                // double lt = v2lapltau_aa + v2lapltau_ab;

                // Third-oder

                // // sigma and gamma
                double xxx = 4.0 * v3sigma3_ccc + 8.0 * v3sigma3_ccb + 4.0 * v3sigma3_cbb + 16.0 * v3sigma3_acc + 24.0 * v3sigma3_acb +
                             8.0 * v3sigma3_abb + 20.0 * v3sigma3_aac + 16.0 * v3sigma3_aab + 8.0 * v3sigma3_aaa;

                double xxr = 2.0 * v3rhosigma2_bcc + 2.0 * v3rhosigma2_bcb + 6.0 * v3rhosigma2_bac + 4.0 * v3rhosigma2_bab + 4.0 * v3rhosigma2_baa +
                             2.0 * v3rhosigma2_acc + 2.0 * v3rhosigma2_acb + 6.0 * v3rhosigma2_aac + 4.0 * v3rhosigma2_aab + 4.0 * v3rhosigma2_aaa;

                double xrt = v3rhosigmatau_bcb + v3rhosigmatau_bca + 2.0 * v3rhosigmatau_bab + 2.0 * v3rhosigmatau_baa + 2.0 * v3rhosigmatau_aab +
                             2.0 * v3rhosigmatau_aaa + v3rhosigmatau_acb + v3rhosigmatau_aca;

                // double xxl  = 2.0 * v3sigma2lapl_ccb + 2.0 * v3sigma2lapl_cca + 2.0 * v3sigma2lapl_cbb
                //             + 2.0 * v3sigma2lapl_cba + 6.0 * v3sigma2lapl_acb + 6.0 * v3sigma2lapl_aca
                //             + 4.0 * v3sigma2lapl_abb + 4.0 * v3sigma2lapl_aba + 4.0 * v3sigma2lapl_aab + 4.0 * v3sigma2lapl_aaa;
                double xrr =
                    v3rho2sigma_bbc + 2.0 * v3rho2sigma_bba + 2.0 * v3rho2sigma_abc + 4.0 * v3rho2sigma_aba + v3rho2sigma_aac + 2.0 * v3rho2sigma_aaa;

                // double xrl  = v3rhosigmalapl_bcb + v3rhosigmalapl_bca + 2.0 * v3rhosigmalapl_bab
                //             + 2.0 * v3rhosigmalapl_baa + v3rhosigmalapl_acb + v3rhosigmalapl_aca
                //             + 2.0 * v3rhosigmalapl_aab + 2.0 * v3rhosigmalapl_aaa;

                // double xtl  = v3sigmalapltau_cbb + v3sigmalapltau_cba + v3sigmalapltau_cab
                //             + v3sigmalapltau_caa + 2.0 * v3sigmalapltau_abb + 2.0 * v3sigmalapltau_aba
                //             + 2.0 * v3sigmalapltau_aab + 2.0 * v3sigmalapltau_aaa;
                // double xll  = v3sigmalapl2_cbb + 2.0 * v3sigmalapl2_cab + v3sigmalapl2_caa
                //             + 2.0 * v3sigmalapl2_abb + 4.0 * v3sigmalapl2_aab + 2.0 * v3sigmalapl2_aaa;

                double xxt = 2.0 * v3sigma2tau_ccb + 2.0 * v3sigma2tau_cca + 2.0 * v3sigma2tau_cbb + 2.0 * v3sigma2tau_cba + 6.0 * v3sigma2tau_acb +
                             6.0 * v3sigma2tau_aca + 4.0 * v3sigma2tau_abb + 4.0 * v3sigma2tau_aba + 4.0 * v3sigma2tau_aab + 4.0 * v3sigma2tau_aaa;

                double xtt =
                    v3sigmatau2_cbb + 2.0 * v3sigmatau2_cab + v3sigmatau2_caa + 2.0 * v3sigmatau2_abb + 4.0 * v3sigmatau2_aab + 2.0 * v3sigmatau2_aaa;

                // rho
                double rrr = v3rho3_abb + 2.0 * v3rho3_aab + v3rho3_aaa;
                double rrt = v3rho2tau_abb + v3rho2tau_aba + v3rho2tau_aab + v3rho2tau_aaa;
                double rtx = 2.0 * v3rhosigmatau_acb + 2.0 * v3rhosigmatau_aca + 2.0 * v3rhosigmatau_abb + 2.0 * v3rhosigmatau_aba +
                             2.0 * v3rhosigmatau_aab + 2.0 * v3rhosigmatau_aaa;
                // double rrl  = v3rho2lapl_abb + v3rho2lapl_aba + v3rho2lapl_aab + v3rho2lapl_aaa;
                double rrx = 2.0 * v3rho2sigma_abc + 2.0 * v3rho2sigma_abb + 2.0 * v3rho2sigma_aba + 2.0 * v3rho2sigma_aac + 2.0 * v3rho2sigma_aab +
                             2.0 * v3rho2sigma_aaa;
                // double rtl  = v3rholapltau_abb + v3rholapltau_aba + v3rholapltau_aab + v3rholapltau_aaa;

                double rtt = v3rhotau2_abb + 2.0 * v3rhotau2_aab + v3rhotau2_aaa;

                // double rll  = v3rholapl2_abb + 2.0 * v3rholapl2_aab + v3rholapl2_aaa;
                // double rlx  = 2.0 * v3rhosigmalapl_acb + 2.0 * v3rhosigmalapl_aca + 2.0 * v3rhosigmalapl_abb
                //             + 2.0 * v3rhosigmalapl_aba + 2.0 * v3rhosigmalapl_aab + 2.0 * v3rhosigmalapl_aaa;
                double rxx = 4.0 * v3rhosigma2_acc + 8.0 * v3rhosigma2_acb + 4.0 * v3rhosigma2_abb + 8.0 * v3rhosigma2_aac + 8.0 * v3rhosigma2_aab +
                             4.0 * v3rhosigma2_aaa;

                // laplacian
                // double lll  = v3lapl3_abb + 2.0 * v3lapl3_aab + v3lapl3_aaa;
                // double llr  = v3rholapl2_bab + v3rholapl2_baa + v3rholapl2_aab + v3rholapl2_aaa;
                // double llt  = v3lapl2tau_abb + v3lapl2tau_aba + v3lapl2tau_aab + v3lapl2tau_aaa;
                // double llx  = 2.0 * v3sigmalapl2_cab + 2.0 * v3sigmalapl2_caa + 2.0 * v3sigmalapl2_bab
                //            + 2.0 * v3sigmalapl2_baa + 2.0 * v3sigmalapl2_aab + 2.0 * v3sigmalapl2_aaa;
                // double lrr  = v3rho2lapl_bba + 2.0 * v3rho2lapl_aba + v3rho2lapl_aaa;
                // double lrt  = v3rholapltau_bab + v3rholapltau_baa + v3rholapltau_aab + v3rholapltau_aaa;
                // double lrx  = 2.0 * v3rhosigmalapl_bca + 2.0 * v3rhosigmalapl_bba + 2.0 * v3rhosigmalapl_baa
                //            + 2.0 * v3rhosigmalapl_aca + 2.0 * v3rhosigmalapl_aba + 2.0 * v3rhosigmalapl_aaa;
                // double ltt  = v3lapltau2_abb + 2.0 * v3lapltau2_aab + v3lapltau2_aaa;
                // double ltx  = 2.0 * v3sigmalapltau_cab + 2.0 * v3sigmalapltau_caa + 2.0 * v3sigmalapltau_bab
                //            + 2.0 * v3sigmalapltau_baa + 2.0 * v3sigmalapltau_aab + 2.0 * v3sigmalapltau_aaa;
                // double lxx  = 4.0 * v3sigma2lapl_cca + 8.0 * v3sigma2lapl_cba + 4.0 * v3sigma2lapl_bba
                //            + 8.0 * v3sigma2lapl_aca + 8.0 * v3sigma2lapl_aba + 4.0 * v3sigma2lapl_aaa;

                // tau
                double trr = v3rho2tau_bba + 2.0 * v3rho2tau_aba + v3rho2tau_aaa;
                // double ttl  = v3lapltau2_bab + v3lapltau2_baa + v3lapltau2_aab + v3lapltau2_aaa;

                // double trl  = v3rholapltau_bba + v3rholapltau_baa + v3rholapltau_aba + v3rholapltau_aaa;

                // double tll  = v3lapl2tau_bba + 2.0 * v3lapl2tau_aba + v3lapl2tau_aaa;
                // double tlx  = 2.0 * v3sigmalapltau_cba + 2.0 * v3sigmalapltau_caa + 2.0 * v3sigmalapltau_bba
                //             + 2.0 * v3sigmalapltau_baa + 2.0 * v3sigmalapltau_aba + 2.0 * v3sigmalapltau_aaa;

                double ttt = v3tau3_abb + 2.0 * v3tau3_aab + v3tau3_aaa;

                double ttr = v3rhotau2_bab + v3rhotau2_baa + v3rhotau2_aab + v3rhotau2_aaa;

                double trx = 2.0 * v3rhosigmatau_bca + 2.0 * v3rhosigmatau_bba + 2.0 * v3rhosigmatau_baa + 2.0 * v3rhosigmatau_aca +
                             2.0 * v3rhosigmatau_aba + 2.0 * v3rhosigmatau_aaa;

                double txx = 4.0 * v3sigma2tau_cca + 8.0 * v3sigma2tau_cba + 4.0 * v3sigma2tau_bba + 8.0 * v3sigma2tau_aca + 8.0 * v3sigma2tau_aba +
                             4.0 * v3sigma2tau_aaa;

                double ttx = 2.0 * v3sigmatau2_cab + 2.0 * v3sigmatau2_caa + 2.0 * v3sigmatau2_bab + 2.0 * v3sigmatau2_baa + 2.0 * v3sigmatau2_aab +
                             2.0 * v3sigmatau2_aaa;

                double w = weights[g];

                double rxw12a = gradw12a_x[g];
                double ryw12a = gradw12a_y[g];
                double rzw12a = gradw12a_z[g];

                double grada_x_g = rhograd[6 * g + 0];
                double grada_y_g = rhograd[6 * g + 1];
                double grada_z_g = rhograd[6 * g + 2];

                double l2contract = grada_x_g * rxw12a + grada_y_g * ryw12a + grada_z_g * rzw12a;

                double l5contract_x = grada_x_g * l2contract;
                double l5contract_y = grada_y_g * l2contract;
                double l5contract_z = grada_z_g * l2contract;

                double q2contract = grada_x_g * gamx[g] + grada_y_g * gamy[g] + grada_z_g * gamz[g];
                // double sl_q2contract = grada_x_g * sl_gamx[g] + grada_y_g * sl_gamy[g] + grada_z_g * sl_gamz[g];
                double st_q2contract = grada_x_g * st_gamx[g] + grada_y_g * st_gamy[g] + grada_z_g * st_gamz[g];

                double q3contract = grada_x_g * grada_x_g * gamxx[g] + grada_x_g * grada_y_g * gamxy[g] + grada_x_g * grada_z_g * gamxz[g] +
                                    grada_y_g * grada_x_g * gamyx[g] + grada_y_g * grada_y_g * gamyy[g] + grada_y_g * grada_z_g * gamyz[g] +
                                    grada_z_g * grada_x_g * gamzx[g] + grada_z_g * grada_y_g * gamzy[g] + grada_z_g * grada_z_g * gamzz[g];

                double q4contract = gamxx[g] + gamyy[g] + gamzz[g];

                double q7contract_x = grada_x_g * grada_x_g * gamx[g] + grada_x_g * grada_y_g * gamy[g] + grada_x_g * grada_z_g * gamz[g];
                double q7contract_y = grada_y_g * grada_x_g * gamx[g] + grada_y_g * grada_y_g * gamy[g] + grada_y_g * grada_z_g * gamz[g];
                double q7contract_z = grada_z_g * grada_x_g * gamx[g] + grada_z_g * grada_y_g * gamy[g] + grada_z_g * grada_z_g * gamz[g];

                // double sl_q7contract_x =  grada_x_g * grada_x_g * sl_gamx[g] + grada_x_g * grada_y_g * sl_gamy[g] + grada_x_g * grada_z_g *
                // sl_gamz[g]; double sl_q7contract_y =  grada_y_g * grada_x_g * sl_gamx[g] + grada_y_g * grada_y_g * sl_gamy[g] + grada_y_g *
                // grada_z_g * sl_gamz[g]; double sl_q7contract_z =  grada_z_g * grada_x_g * sl_gamx[g] + grada_z_g * grada_y_g * sl_gamy[g] +
                // grada_z_g * grada_z_g * sl_gamz[g];

                double st_q7contract_x = grada_x_g * grada_x_g * st_gamx[g] + grada_x_g * grada_y_g * st_gamy[g] + grada_x_g * grada_z_g * st_gamz[g];
                double st_q7contract_y = grada_y_g * grada_x_g * st_gamx[g] + grada_y_g * grada_y_g * st_gamy[g] + grada_y_g * grada_z_g * st_gamz[g];
                double st_q7contract_z = grada_z_g * grada_x_g * st_gamx[g] + grada_z_g * grada_y_g * st_gamy[g] + grada_z_g * grada_z_g * st_gamz[g];

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

                // Rho operator contributions

                // vxc 1 contributions

                double rho_0 = rr * rhow12a[g] + rx * l2contract + rt * tauw12a[g];
                //+ rl * laplw12a[g];

                // double lap_0 =    lr * rhow12a[g]
                //                 + lx * l2contract
                //                 + lt * tauw12a[g]
                //                 + ll * laplw12a[g];

                double tau_0 = tr * rhow12a[g] + tx * l2contract + tt * tauw12a[g];
                //+ tl * laplw12a[g];

                // vxc 2 contributions

                rho_0 += rrr * gam[g] +
                         rrt * rt_gam[g]
                         //+ rrl * rl_gam[g]
                         //+ rll * ll_gam[g]
                         + rtt * tt_gam[g]
                         //+ rtl * tl_gam[g]
                         + rrx * q2contract
                         //+ rlx * sl_q2contract
                         + rtx * st_q2contract + rxx * q3contract + rx * q4contract;

                // lap_0 += lrr * gam[g]
                //        + lrt * rt_gam[g]
                //        + llr * rl_gam[g]
                //        + lll * ll_gam[g]
                //        + ltt * tt_gam[g]
                //        + llt * tl_gam[g]
                //        + lrx * q2contract
                //        + llx * sl_q2contract
                //        + ltx * st_q2contract
                //        + lxx * q3contract
                //        + lx  * q4contract;

                tau_0 += trr * gam[g] +
                         ttr * rt_gam[g]
                         //+ trl * rl_gam[g]
                         //+ tll * ll_gam[g]
                         + ttt * tt_gam[g]
                         //+ ttl * tl_gam[g]
                         + trx * q2contract
                         //+ tlx * sl_q2contract
                         + ttx * st_q2contract + txx * q3contract + tx * q4contract;

                // Grad operator contributions

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                // vxc 1 contributions

                // xcomp +=  grada_x_g * ( xr * rhow12a[g] + xt * tauw12a[g] + xl * laplw12a[g])
                xcomp += grada_x_g * (xr * rhow12a[g] + xt * tauw12a[g]) + x * rxw12a + xx * l5contract_x;

                // ycomp += grada_y_g * ( xr * rhow12a[g] + xt * tauw12a[g] + xl * laplw12a[g])
                ycomp += grada_y_g * (xr * rhow12a[g] + xt * tauw12a[g]) + x * ryw12a + xx * l5contract_y;

                // zcomp += grada_z_g * ( xr * rhow12a[g] + xt * tauw12a[g] + xl * laplw12a[g])
                zcomp += grada_z_g * (xr * rhow12a[g] + xt * tauw12a[g]) + x * rzw12a + xx * l5contract_z;

                // vxc 2 contributions

                xcomp += xrr * grada_x_g * gam[g] +
                         xrt * grada_x_g * rt_gam[g]
                         //+ xrl * grada_x_g * rl_gam[g]
                         //+ xll * grada_x_g * ll_gam[g]
                         + xtt * grada_x_g * tt_gam[g]
                         //+ xtl * grada_x_g * tl_gam[g]
                         + xr * gamx[g]  // q6
                         //+ xl * sl_gamx[g]
                         + xt * st_gamx[g] +
                         xxr * q7contract_x
                         //+ xxl * sl_q7contract_x
                         + xxt * st_q7contract_x + xx * (q8contract_x + q10contract_x + q11contract_x) + xxx * q9contract_x;

                ycomp += xrr * grada_y_g * gam[g]  // q5
                         + xrt * grada_y_g * rt_gam[g]
                         //+ xrl * grada_y_g * rl_gam[g]
                         //+ xll * grada_y_g * ll_gam[g]
                         + xtt * grada_y_g * tt_gam[g]
                         //+ xtl * grada_y_g * tl_gam[g]
                         + xr * gamy[g]  // q6
                         //+ xl * sl_gamy[g]
                         + xt * st_gamy[g] +
                         xxr * q7contract_y
                         //+ xxl * sl_q7contract_y
                         + xxt * st_q7contract_y + xx * (q8contract_y + q10contract_y + q11contract_y) + xxx * q9contract_y;

                zcomp += xrr * grada_z_g * gam[g]  // q5
                         + xrt * grada_z_g * rt_gam[g]
                         //+ xrl * grada_z_g * rl_gam[g]
                         //+ xll * grada_z_g * ll_gam[g]
                         + xtt * grada_z_g * tt_gam[g]
                         //+ xtl * grada_z_g * tl_gam[g]
                         + xr * gamz[g]  // q6
                         //+ xl * sl_gamz[g]
                         + xt * st_gamz[g] +
                         xxr * q7contract_z
                         //+ xxl * sl_q7contract_z
                         + xxt * st_q7contract_z + xx * (q8contract_z + q10contract_z + q11contract_z) + xxx * q9contract_z;

                G_val[nu_offset + g] = w * rho_0 * chi_val[nu_offset + g];

                G_gga_val[nu_offset + g] =
                    w * (xcomp * chi_x_val[nu_offset + g] + ycomp * chi_y_val[nu_offset + g] + zcomp * chi_z_val[nu_offset + g]);

                // TODO implement Laplacian dependence

                // G_gga_val[nu_offset + g] += w * lap_0 * (chi_xx_val[nu_offset + g] +
                //                                          chi_yy_val[nu_offset + g] +
                //                                          chi_zz_val[nu_offset + g]);

                G_gga_x_val[nu_offset + g] = w * tau_0 * chi_x_val[nu_offset + g];
                G_gga_y_val[nu_offset + g] = w * tau_0 * chi_y_val[nu_offset + g];
                G_gga_z_val[nu_offset + g] = w * tau_0 * chi_z_val[nu_offset + g];
            }
        }
    }

    timer.stop("Kxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Kxc matrix matmul");

    // LDA and GGA contribution
    auto mat_Kxc = denblas::multABt(gtoValues, mat_G);

    auto mat_Kxc_gga = denblas::multABt(gtoValues, mat_G_gga);

    mat_Kxc_gga.symmetrize();  // (matrix + matrix.T)

    mat_Kxc = denblas::addAB(mat_Kxc, mat_Kxc_gga, 1.0);

    // tau contribution
    auto mat_Kxc_x = denblas::multABt(gtoValuesX, mat_G_gga_x);
    auto mat_Kxc_y = denblas::multABt(gtoValuesY, mat_G_gga_y);
    auto mat_Kxc_z = denblas::multABt(gtoValuesZ, mat_G_gga_z);

    mat_Kxc = denblas::addAB(mat_Kxc, mat_Kxc_x, 0.5);
    mat_Kxc = denblas::addAB(mat_Kxc, mat_Kxc_y, 0.5);
    mat_Kxc = denblas::addAB(mat_Kxc, mat_Kxc_z, 0.5);

    timer.stop("Kxc matrix matmul");

    return mat_Kxc;
}

}  // namespace xcintmgga
