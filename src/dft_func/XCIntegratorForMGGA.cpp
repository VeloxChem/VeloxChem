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

}  // namespace xcintmgga
