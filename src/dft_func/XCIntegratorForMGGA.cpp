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
integrateVxcFockForMetaGgaClosedShell(const CMolecule&                  molecule,
                                      const CMolecularBasis&            basis,
                                      const std::vector<const double*>& gsDensityPointers,
                                      const CMolecularGrid&             molecularGrid,
                                      const double                      screeningThresholdForGTOValues,
                                      const CXCFunctional&              xcFunctional) -> CAOKohnShamMatrix
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

    CAOKohnShamMatrix mat_Vxc(naos, naos, std::string("closedshell"));

    mat_Vxc.zero();

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    auto omp_max_npoints = max_npoints_per_box / nthreads;
    if (max_npoints_per_box % nthreads != 0) omp_max_npoints++;

    // density and functional derivatives

    auto       mggafunc = xcFunctional.getFunctionalPointerToMetaGgaComponent();
    const auto dim      = &(mggafunc->dim);

    std::vector<CXCFunctional> omp_xcfuncs(nthreads, CXCFunctional(xcFunctional));

    std::vector<std::vector<double>> omp_local_weights_data(nthreads, std::vector<double>(omp_max_npoints));

    std::vector<std::vector<double>> omp_rho_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));
    std::vector<std::vector<double>> omp_rhograd_data(nthreads, std::vector<double>(dim->rho * 3 * omp_max_npoints));
    std::vector<std::vector<double>> omp_sigma_data(nthreads, std::vector<double>(dim->sigma * omp_max_npoints));
    std::vector<std::vector<double>> omp_lapl_data(nthreads, std::vector<double>(dim->lapl * omp_max_npoints));
    std::vector<std::vector<double>> omp_tau_data(nthreads, std::vector<double>(dim->tau * omp_max_npoints));

    std::vector<std::vector<double>> omp_exc_data(nthreads, std::vector<double>(dim->zk * omp_max_npoints));
    std::vector<std::vector<double>> omp_vrho_data(nthreads, std::vector<double>(dim->vrho * omp_max_npoints));
    std::vector<std::vector<double>> omp_vsigma_data(nthreads, std::vector<double>(dim->vsigma * omp_max_npoints));
    std::vector<std::vector<double>> omp_vlapl_data(nthreads, std::vector<double>(dim->vlapl * omp_max_npoints));
    std::vector<std::vector<double>> omp_vtau_data(nthreads, std::vector<double>(dim->vtau * omp_max_npoints));

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

        timer.start("Density matrix slicing");

        auto sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

        timer.stop("Density matrix slicing");

        // GTO values on grid points

        timer.start("OMP Vxc calc.");

        CDenseMatrix sum_partial_mat_Vxc(aocount, aocount);

        sum_partial_mat_Vxc.zero();

#pragma omp parallel reduction(+ : nele, xcene)
        {
            auto thread_id = omp_get_thread_num();

            omptimers[thread_id].start("gtoeval");

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            CDenseMatrix mat_chi(aocount, grid_batch_size);
            CDenseMatrix mat_chi_x(aocount, grid_batch_size);
            CDenseMatrix mat_chi_y(aocount, grid_batch_size);
            CDenseMatrix mat_chi_z(aocount, grid_batch_size);

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
                    std::memcpy(mat_chi.row(idx), submat_0_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_x.row(idx), submat_x_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_y.row(idx), submat_y_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_z.row(idx), submat_z_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");

            omptimers[thread_id].start("Generate density grid");

            auto local_weights = omp_local_weights_data[thread_id].data();

            auto rho     = omp_rho_data[thread_id].data();
            auto rhograd = omp_rhograd_data[thread_id].data();
            auto sigma   = omp_sigma_data[thread_id].data();
            auto lapl    = omp_lapl_data[thread_id].data();
            auto tau     = omp_tau_data[thread_id].data();

            auto exc    = omp_exc_data[thread_id].data();
            auto vrho   = omp_vrho_data[thread_id].data();
            auto vsigma = omp_vsigma_data[thread_id].data();
            auto vlapl  = omp_vlapl_data[thread_id].data();
            auto vtau   = omp_vtau_data[thread_id].data();

            dengridgen::serialGenerateDensityForMGGA(rho, rhograd, sigma, lapl, tau, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat);

            omptimers[thread_id].stop("Generate density grid");

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_exc_vxc_for_mgga(grid_batch_size, rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            omptimers[thread_id].start("Vxc matrix G");

            // LDA contribution
            CDenseMatrix mat_G(aocount, grid_batch_size);

            // GGA contribution
            CDenseMatrix mat_G_gga(aocount, grid_batch_size);

            // tau contribution
            CDenseMatrix mat_G_gga_x(aocount, grid_batch_size);
            CDenseMatrix mat_G_gga_y(aocount, grid_batch_size);
            CDenseMatrix mat_G_gga_z(aocount, grid_batch_size);

            auto G_val = mat_G.values();

            auto G_gga_val = mat_G_gga.values();

            auto G_gga_x_val = mat_G_gga_x.values();
            auto G_gga_y_val = mat_G_gga_y.values();
            auto G_gga_z_val = mat_G_gga_z.values();

            auto chi_val   = mat_chi.values();

            auto chi_x_val = mat_chi_x.values();
            auto chi_y_val = mat_chi_y.values();
            auto chi_z_val = mat_chi_z.values();

            for (int nu = 0; nu < aocount; nu++)
            {
                auto nu_offset = nu * grid_batch_size;

#pragma omp simd
                for (int g = 0; g < grid_batch_size; g++)
                {
                    auto vx = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                    auto vy = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                    auto vz = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                    // LDA contribution
                    G_val[nu_offset + g] = local_weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];

                    // GGA contribution (will be scaled by 2 later)
                    G_gga_val[nu_offset + g] =
                        local_weights[g] * (vx * chi_x_val[nu_offset + g] + vy * chi_y_val[nu_offset + g] + vz * chi_z_val[nu_offset + g]);

                    // tau contribution (will be scaled by 0.5 later)
                    G_gga_x_val[nu_offset + g] = local_weights[g] * vtau[2 * g + 0] * chi_x_val[nu_offset + g];
                    G_gga_y_val[nu_offset + g] = local_weights[g] * vtau[2 * g + 0] * chi_y_val[nu_offset + g];
                    G_gga_z_val[nu_offset + g] = local_weights[g] * vtau[2 * g + 0] * chi_z_val[nu_offset + g];
                }
            }

            omptimers[thread_id].stop("Vxc matrix G");

            // Note that we use matrix-matrix multiplication only once, and symmetrize
            // the result. This is because the density matrix is symmetric, and the
            // Kohn-Sham matrix from mat_G is also symmetric. Formally only the
            // mat_G_gga contribution should be symmetrized.

            omptimers[thread_id].start("Vxc matmul and symm.");

            // LDA and GGA contribution
            auto partial_mat_Vxc = denblas::serialMultABt(mat_chi, denblas::serialAddAB(mat_G, mat_G_gga, 2.0));

            // tau contribution
            auto partial_mat_Vxc_x = denblas::serialMultABt(mat_chi_x, mat_G_gga_x);
            auto partial_mat_Vxc_y = denblas::serialMultABt(mat_chi_y, mat_G_gga_y);
            auto partial_mat_Vxc_z = denblas::serialMultABt(mat_chi_z, mat_G_gga_z);

            denblas::serialInPlaceAddAB(partial_mat_Vxc, partial_mat_Vxc_x, 0.5);
            denblas::serialInPlaceAddAB(partial_mat_Vxc, partial_mat_Vxc_y, 0.5);
            denblas::serialInPlaceAddAB(partial_mat_Vxc, partial_mat_Vxc_z, 0.5);

            partial_mat_Vxc.symmetrizeAndScale(0.5);

            omptimers[thread_id].stop("Vxc matmul and symm.");

            omptimers[thread_id].start("Vxc local matrix dist.");

#pragma omp critical
            denblas::serialInPlaceAddAB(sum_partial_mat_Vxc, partial_mat_Vxc);

            omptimers[thread_id].stop("Vxc local matrix dist.");

            omptimers[thread_id].start("XC energy and num. elec.");

            double local_nele = 0.0, local_xcene = 0.0;

            for (int g = 0; g < grid_batch_size; g++)
            {
                auto rho_total = rho[2 * g + 0] + rho[2 * g + 1];

                local_nele += local_weights[g] * rho_total;

                local_xcene += local_weights[g] * exc[g] * rho_total;
            }

            nele += local_nele;

            xcene += local_xcene;

            omptimers[thread_id].stop("XC energy and num. elec.");
        }

        timer.stop("OMP Vxc calc.");

        timer.start("Vxc matrix dist.");

        dftsubmat::distributeSubMatrixToKohnSham(mat_Vxc, sum_partial_mat_Vxc, aoinds);

        timer.stop("Vxc matrix dist.");
    }

    mat_Vxc.setNumberOfElectrons(nele);

    mat_Vxc.setExchangeCorrelationEnergy(xcene);

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

    return mat_Vxc;
}

auto
integrateVxcFockForMetaGgaOpenShell(const CMolecule&                  molecule,
                                    const CMolecularBasis&            basis,
                                    const std::vector<const double*>& gsDensityPointers,
                                    const CMolecularGrid&             molecularGrid,
                                    const double                      screeningThresholdForGTOValues,
                                    const CXCFunctional&              xcFunctional) -> CAOKohnShamMatrix
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

    CAOKohnShamMatrix mat_Vxc(naos, naos, std::string("openshell"));

    mat_Vxc.zero();

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    auto omp_max_npoints = max_npoints_per_box / nthreads;
    if (max_npoints_per_box % nthreads != 0) omp_max_npoints++;

    // density and functional derivatives

    auto       mggafunc = xcFunctional.getFunctionalPointerToMetaGgaComponent();
    const auto dim      = &(mggafunc->dim);

    std::vector<CXCFunctional> omp_xcfuncs(nthreads, CXCFunctional(xcFunctional));

    std::vector<std::vector<double>> omp_local_weights_data(nthreads, std::vector<double>(omp_max_npoints));

    std::vector<std::vector<double>> omp_rho_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));
    std::vector<std::vector<double>> omp_rhograd_data(nthreads, std::vector<double>(dim->rho * 3 * omp_max_npoints));
    std::vector<std::vector<double>> omp_sigma_data(nthreads, std::vector<double>(dim->sigma * omp_max_npoints));
    std::vector<std::vector<double>> omp_lapl_data(nthreads, std::vector<double>(dim->lapl * omp_max_npoints));
    std::vector<std::vector<double>> omp_tau_data(nthreads, std::vector<double>(dim->tau * omp_max_npoints));

    std::vector<std::vector<double>> omp_exc_data(nthreads, std::vector<double>(dim->zk * omp_max_npoints));
    std::vector<std::vector<double>> omp_vrho_data(nthreads, std::vector<double>(dim->vrho * omp_max_npoints));
    std::vector<std::vector<double>> omp_vsigma_data(nthreads, std::vector<double>(dim->vsigma * omp_max_npoints));
    std::vector<std::vector<double>> omp_vlapl_data(nthreads, std::vector<double>(dim->vlapl * omp_max_npoints));
    std::vector<std::vector<double>> omp_vtau_data(nthreads, std::vector<double>(dim->vtau * omp_max_npoints));

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

        timer.start("Density matrix slicing");

        auto sub_dens_mat_a = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);
        auto sub_dens_mat_b = dftsubmat::getSubDensityMatrix(gsDensityPointers[1], aoinds, naos);

        timer.stop("Density matrix slicing");

        // GTO values on grid points

        timer.start("OMP Vxc calc.");

        CDenseMatrix sum_partial_mat_Vxc_a(aocount, aocount);
        CDenseMatrix sum_partial_mat_Vxc_b(aocount, aocount);

        sum_partial_mat_Vxc_a.zero();
        sum_partial_mat_Vxc_b.zero();

#pragma omp parallel reduction(+ : nele, xcene)
        {
            auto thread_id = omp_get_thread_num();

            omptimers[thread_id].start("gtoeval");

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            CDenseMatrix mat_chi(aocount, grid_batch_size);
            CDenseMatrix mat_chi_x(aocount, grid_batch_size);
            CDenseMatrix mat_chi_y(aocount, grid_batch_size);
            CDenseMatrix mat_chi_z(aocount, grid_batch_size);

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
                    std::memcpy(mat_chi.row(idx), submat_0_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_x.row(idx), submat_x_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_y.row(idx), submat_y_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_z.row(idx), submat_z_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");

            omptimers[thread_id].start("Generate density grid");

            auto local_weights = omp_local_weights_data[thread_id].data();

            auto rho     = omp_rho_data[thread_id].data();
            auto rhograd = omp_rhograd_data[thread_id].data();
            auto sigma   = omp_sigma_data[thread_id].data();
            auto lapl    = omp_lapl_data[thread_id].data();
            auto tau     = omp_tau_data[thread_id].data();

            auto exc    = omp_exc_data[thread_id].data();
            auto vrho   = omp_vrho_data[thread_id].data();
            auto vsigma = omp_vsigma_data[thread_id].data();
            auto vlapl  = omp_vlapl_data[thread_id].data();
            auto vtau   = omp_vtau_data[thread_id].data();

            dengridgen::serialGenerateDensityForMGGA(
                rho, rhograd, sigma, lapl, tau, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat_a, sub_dens_mat_b);

            omptimers[thread_id].stop("Generate density grid");

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_exc_vxc_for_mgga(grid_batch_size, rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            omptimers[thread_id].start("Vxc matrix G");

            // LDA contribution
            CDenseMatrix mat_G_a(aocount, grid_batch_size);
            CDenseMatrix mat_G_b(aocount, grid_batch_size);

            // GGA contribution
            CDenseMatrix mat_G_a_gga(aocount, grid_batch_size);
            CDenseMatrix mat_G_b_gga(aocount, grid_batch_size);

            // tau contribution
            CDenseMatrix mat_G_a_gga_x(aocount, grid_batch_size);
            CDenseMatrix mat_G_a_gga_y(aocount, grid_batch_size);
            CDenseMatrix mat_G_a_gga_z(aocount, grid_batch_size);

            CDenseMatrix mat_G_b_gga_x(aocount, grid_batch_size);
            CDenseMatrix mat_G_b_gga_y(aocount, grid_batch_size);
            CDenseMatrix mat_G_b_gga_z(aocount, grid_batch_size);

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

            auto chi_val   = mat_chi.values();

            auto chi_x_val = mat_chi_x.values();
            auto chi_y_val = mat_chi_y.values();
            auto chi_z_val = mat_chi_z.values();

            for (int nu = 0; nu < aocount; nu++)
            {
                auto nu_offset = nu * grid_batch_size;

#pragma omp simd
                for (int g = 0; g < grid_batch_size; g++)
                {
                    auto vxa = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                    auto vya = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                    auto vza = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                    auto vxb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 3] + vsigma[3 * g + 1] * rhograd[6 * g + 0];
                    auto vyb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 4] + vsigma[3 * g + 1] * rhograd[6 * g + 1];
                    auto vzb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 5] + vsigma[3 * g + 1] * rhograd[6 * g + 2];

                    // LDA contribution
                    G_a_val[nu_offset + g] = local_weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];
                    G_b_val[nu_offset + g] = local_weights[g] * vrho[2 * g + 1] * chi_val[nu_offset + g];

                    // GGA contribution (will be scaled by 2 later)
                    G_a_gga_val[nu_offset + g] =
                        local_weights[g] * (vxa * chi_x_val[nu_offset + g] + vya * chi_y_val[nu_offset + g] + vza * chi_z_val[nu_offset + g]);
                    G_b_gga_val[nu_offset + g] =
                        local_weights[g] * (vxb * chi_x_val[nu_offset + g] + vyb * chi_y_val[nu_offset + g] + vzb * chi_z_val[nu_offset + g]);

                    // TODO implement Laplacian dependence

                    // tau contribution (will be scaled by 0.5 later)
                    G_a_gga_x_val[nu_offset + g] = local_weights[g] * vtau[2 * g + 0] * chi_x_val[nu_offset + g];
                    G_a_gga_y_val[nu_offset + g] = local_weights[g] * vtau[2 * g + 0] * chi_y_val[nu_offset + g];
                    G_a_gga_z_val[nu_offset + g] = local_weights[g] * vtau[2 * g + 0] * chi_z_val[nu_offset + g];

                    G_b_gga_x_val[nu_offset + g] = local_weights[g] * vtau[2 * g + 1] * chi_x_val[nu_offset + g];
                    G_b_gga_y_val[nu_offset + g] = local_weights[g] * vtau[2 * g + 1] * chi_y_val[nu_offset + g];
                    G_b_gga_z_val[nu_offset + g] = local_weights[g] * vtau[2 * g + 1] * chi_z_val[nu_offset + g];
                }
            }

            omptimers[thread_id].stop("Vxc matrix G");

            // Note that we use matrix-matrix multiplication only once, and symmetrize
            // the result. This is because the density matrix is symmetric, and the
            // Kohn-Sham matrix from mat_G is also symmetric. Formally only the
            // mat_G_gga contribution should be symmetrized.

            omptimers[thread_id].start("Vxc matmul and symm.");

            // LDA and GGA contribution
            auto partial_mat_Vxc_a = denblas::serialMultABt(mat_chi, denblas::serialAddAB(mat_G_a, mat_G_a_gga, 2.0));
            auto partial_mat_Vxc_b = denblas::serialMultABt(mat_chi, denblas::serialAddAB(mat_G_b, mat_G_b_gga, 2.0));

            // tau contribution
            auto partial_mat_Vxc_a_x = denblas::serialMultABt(mat_chi_x, mat_G_a_gga_x);
            auto partial_mat_Vxc_a_y = denblas::serialMultABt(mat_chi_y, mat_G_a_gga_y);
            auto partial_mat_Vxc_a_z = denblas::serialMultABt(mat_chi_z, mat_G_a_gga_z);

            auto partial_mat_Vxc_b_x = denblas::serialMultABt(mat_chi_x, mat_G_b_gga_x);
            auto partial_mat_Vxc_b_y = denblas::serialMultABt(mat_chi_y, mat_G_b_gga_y);
            auto partial_mat_Vxc_b_z = denblas::serialMultABt(mat_chi_z, mat_G_b_gga_z);

            denblas::serialInPlaceAddAB(partial_mat_Vxc_a, partial_mat_Vxc_a_x, 0.5);
            denblas::serialInPlaceAddAB(partial_mat_Vxc_a, partial_mat_Vxc_a_y, 0.5);
            denblas::serialInPlaceAddAB(partial_mat_Vxc_a, partial_mat_Vxc_a_z, 0.5);

            denblas::serialInPlaceAddAB(partial_mat_Vxc_b, partial_mat_Vxc_b_x, 0.5);
            denblas::serialInPlaceAddAB(partial_mat_Vxc_b, partial_mat_Vxc_b_y, 0.5);
            denblas::serialInPlaceAddAB(partial_mat_Vxc_b, partial_mat_Vxc_b_z, 0.5);

            partial_mat_Vxc_a.symmetrizeAndScale(0.5);
            partial_mat_Vxc_b.symmetrizeAndScale(0.5);

            omptimers[thread_id].stop("Vxc matmul and symm.");

            omptimers[thread_id].start("Vxc local matrix dist.");

#pragma omp critical
            {
                denblas::serialInPlaceAddAB(sum_partial_mat_Vxc_a, partial_mat_Vxc_a);
                denblas::serialInPlaceAddAB(sum_partial_mat_Vxc_b, partial_mat_Vxc_b);
            }

            omptimers[thread_id].stop("Vxc local matrix dist.");

            omptimers[thread_id].start("XC energy and num. elec.");

            double local_nele = 0.0, local_xcene = 0.0;

            for (int g = 0; g < grid_batch_size; g++)
            {
                auto rho_total = rho[2 * g + 0] + rho[2 * g + 1];

                local_nele += local_weights[g] * rho_total;

                local_xcene += local_weights[g] * exc[g] * rho_total;
            }

            nele += local_nele;

            xcene += local_xcene;

            omptimers[thread_id].stop("XC energy and num. elec.");
        }

        timer.stop("OMP Vxc calc.");

        timer.start("Vxc matrix dist.");

        dftsubmat::distributeSubMatrixToKohnSham(mat_Vxc, sum_partial_mat_Vxc_a, sum_partial_mat_Vxc_b, aoinds);

        timer.stop("Vxc matrix dist.");
    }

    mat_Vxc.setNumberOfElectrons(nele);

    mat_Vxc.setExchangeCorrelationEnergy(xcene);

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

    for (size_t box_id = 0; box_id < counts.size(); box_id++)
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
                        const std::vector<const double*>& rwDensityPointers,
                        const std::vector<const double*>& rw2DensityPointers,
                        const std::vector<const double*>& gsDensityPointers,
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

    for (size_t box_id = 0; box_id < counts.size(); box_id++)
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

        auto sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForMGGA(rho, rhograd, sigma, lapl, tau, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat, timer);

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto rw_sub_dens_mat = dftsubmat::getSubAODensityMatrix(rwDensityPointers, aoinds, naos);

        auto rw2_sub_dens_mat = dftsubmat::getSubAODensityMatrix(rw2DensityPointers, aoinds, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto rwdengrid = dengridgen::generateDensityGridForMGGA(mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw_sub_dens_mat, xcfuntype, timer);

        auto rw2dengrid =
            dengridgen::generateDensityGridForMGGA(mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw2_sub_dens_mat, xcfuntype, timer);

        // compute perturbed density

        timer.start("Density grid quad");

        auto numdens_rw2 = static_cast<int>(rw2DensityPointers.size());

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

auto
integrateKxcLxcFockForMGGA(const std::vector<double*>& aoFockPointers,
                           const CMolecule&            molecule,
                           const CMolecularBasis&      basis,
                           const std::vector<const double*>& rwDensityPointers,
                           const std::vector<const double*>& rw2DensityPointers,
                           const std::vector<const double*>& rw3DensityPointers,
                           const std::vector<const double*>& gsDensityPointers,
                           const CMolecularGrid&       molecularGrid,
                           const double                screeningThresholdForGTOValues,
                           const CXCFunctional&        xcFunctional,
                           const std::string&          cubeMode) -> void
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

    // Fourth-order
    std::vector<double> v4rho4_data(dim->v4rho4 * max_npoints_per_box);
    std::vector<double> v4rho3sigma_data(dim->v4rho3sigma * max_npoints_per_box);
    std::vector<double> v4rho3lapl_data(dim->v4rho3lapl * max_npoints_per_box);
    std::vector<double> v4rho3tau_data(dim->v4rho3tau * max_npoints_per_box);
    std::vector<double> v4rho2sigma2_data(dim->v4rho2sigma2 * max_npoints_per_box);
    std::vector<double> v4rho2sigmalapl_data(dim->v4rho2sigmalapl * max_npoints_per_box);
    std::vector<double> v4rho2sigmatau_data(dim->v4rho2sigmatau * max_npoints_per_box);
    std::vector<double> v4rho2lapl2_data(dim->v4rho2lapl2 * max_npoints_per_box);
    std::vector<double> v4rho2lapltau_data(dim->v4rho2lapltau * max_npoints_per_box);
    std::vector<double> v4rho2tau2_data(dim->v4rho2tau2 * max_npoints_per_box);
    std::vector<double> v4rhosigma3_data(dim->v4rhosigma3 * max_npoints_per_box);
    std::vector<double> v4rhosigma2lapl_data(dim->v4rhosigma2lapl * max_npoints_per_box);
    std::vector<double> v4rhosigma2tau_data(dim->v4rhosigma2tau * max_npoints_per_box);
    std::vector<double> v4rhosigmalapl2_data(dim->v4rhosigmalapl2 * max_npoints_per_box);
    std::vector<double> v4rhosigmalapltau_data(dim->v4rhosigmalapltau * max_npoints_per_box);
    std::vector<double> v4rhosigmatau2_data(dim->v4rhosigmatau2 * max_npoints_per_box);
    std::vector<double> v4rholapl3_data(dim->v4rholapl3 * max_npoints_per_box);
    std::vector<double> v4rholapl2tau_data(dim->v4rholapl2tau * max_npoints_per_box);
    std::vector<double> v4rholapltau2_data(dim->v4rholapltau2 * max_npoints_per_box);
    std::vector<double> v4rhotau3_data(dim->v4rhotau3 * max_npoints_per_box);
    std::vector<double> v4sigma4_data(dim->v4sigma4 * max_npoints_per_box);
    std::vector<double> v4sigma3lapl_data(dim->v4sigma3lapl * max_npoints_per_box);
    std::vector<double> v4sigma3tau_data(dim->v4sigma3tau * max_npoints_per_box);
    std::vector<double> v4sigma2lapl2_data(dim->v4sigma2lapl2 * max_npoints_per_box);
    std::vector<double> v4sigma2lapltau_data(dim->v4sigma2lapltau * max_npoints_per_box);
    std::vector<double> v4sigma2tau2_data(dim->v4sigma2tau2 * max_npoints_per_box);
    std::vector<double> v4sigmalapl3_data(dim->v4sigmalapl3 * max_npoints_per_box);
    std::vector<double> v4sigmalapl2tau_data(dim->v4sigmalapl2tau * max_npoints_per_box);
    std::vector<double> v4sigmalapltau2_data(dim->v4sigmalapltau2 * max_npoints_per_box);
    std::vector<double> v4sigmatau3_data(dim->v4sigmatau3 * max_npoints_per_box);
    std::vector<double> v4lapl4_data(dim->v4lapl4 * max_npoints_per_box);
    std::vector<double> v4lapl3tau_data(dim->v4lapl3tau * max_npoints_per_box);
    std::vector<double> v4lapl2tau2_data(dim->v4lapl2tau2 * max_npoints_per_box);
    std::vector<double> v4lapltau3_data(dim->v4lapltau3 * max_npoints_per_box);
    std::vector<double> v4tau4_data(dim->v4tau4 * max_npoints_per_box);

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

    // Fourth-order
    auto v4rho4            = v4rho4_data.data();
    auto v4rho3sigma       = v4rho3sigma_data.data();
    auto v4rho3lapl        = v4rho3lapl_data.data();
    auto v4rho3tau         = v4rho3tau_data.data();
    auto v4rho2sigma2      = v4rho2sigma2_data.data();
    auto v4rho2sigmalapl   = v4rho2sigmalapl_data.data();
    auto v4rho2sigmatau    = v4rho2sigmatau_data.data();
    auto v4rho2lapl2       = v4rho2lapl2_data.data();
    auto v4rho2lapltau     = v4rho2lapltau_data.data();
    auto v4rho2tau2        = v4rho2tau2_data.data();
    auto v4rhosigma3       = v4rhosigma3_data.data();
    auto v4rhosigma2lapl   = v4rhosigma2lapl_data.data();
    auto v4rhosigma2tau    = v4rhosigma2tau_data.data();
    auto v4rhosigmalapl2   = v4rhosigmalapl2_data.data();
    auto v4rhosigmalapltau = v4rhosigmalapltau_data.data();
    auto v4rhosigmatau2    = v4rhosigmatau2_data.data();
    auto v4rholapl3        = v4rholapl3_data.data();
    auto v4rholapl2tau     = v4rholapl2tau_data.data();
    auto v4rholapltau2     = v4rholapltau2_data.data();
    auto v4rhotau3         = v4rhotau3_data.data();
    auto v4sigma4          = v4sigma4_data.data();
    auto v4sigma3lapl      = v4sigma3lapl_data.data();
    auto v4sigma3tau       = v4sigma3tau_data.data();
    auto v4sigma2lapl2     = v4sigma2lapl2_data.data();
    auto v4sigma2lapltau   = v4sigma2lapltau_data.data();
    auto v4sigma2tau2      = v4sigma2tau2_data.data();
    auto v4sigmalapl3      = v4sigmalapl3_data.data();
    auto v4sigmalapl2tau   = v4sigmalapl2tau_data.data();
    auto v4sigmalapltau2   = v4sigmalapltau2_data.data();
    auto v4sigmatau3       = v4sigmatau3_data.data();
    auto v4lapl4           = v4lapl4_data.data();
    auto v4lapl3tau        = v4lapl3tau_data.data();
    auto v4lapl2tau2       = v4lapl2tau2_data.data();
    auto v4lapltau3        = v4lapltau3_data.data();
    auto v4tau4            = v4tau4_data.data();

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

        auto sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForMGGA(rho, rhograd, sigma, lapl, tau, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat, timer);

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto rw_sub_dens_mat = dftsubmat::getSubAODensityMatrix(rwDensityPointers, aoinds, naos);

        auto rw2_sub_dens_mat = dftsubmat::getSubAODensityMatrix(rw2DensityPointers, aoinds, naos);

        auto rw3_sub_dens_mat = dftsubmat::getSubAODensityMatrix(rw3DensityPointers, aoinds, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto rwdengrid = dengridgen::generateDensityGridForMGGA(mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw_sub_dens_mat, xcfuntype, timer);

        auto rw2dengrid =
            dengridgen::generateDensityGridForMGGA(mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw2_sub_dens_mat, xcfuntype, timer);

        auto rw3dengrid =
            dengridgen::generateDensityGridForMGGA(mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw3_sub_dens_mat, xcfuntype, timer);

        // compute perturbed density

        timer.start("Density grid quad");

        auto numdens_rw3 = static_cast<int>(rw3DensityPointers.size());

        auto numdens_rw2 = static_cast<int>(rw2DensityPointers.size());

        CDensityGridCubic rwdengridcube(npoints, (numdens_rw2 + numdens_rw3), xcfuntype, dengrid::ab);

        rwdengridcube.DensityProd(rwdengrid, rw2dengrid, xcfuntype, (numdens_rw2 + numdens_rw3), cubeMode);

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

        xcFunctional.compute_lxc_for_mgga(npoints,
                                          rho,
                                          sigma,
                                          lapl,
                                          tau,
                                          v4rho4,
                                          v4rho3sigma,
                                          v4rho3lapl,
                                          v4rho3tau,
                                          v4rho2sigma2,
                                          v4rho2sigmalapl,
                                          v4rho2sigmatau,
                                          v4rho2lapl2,
                                          v4rho2lapltau,
                                          v4rho2tau2,
                                          v4rhosigma3,
                                          v4rhosigma2lapl,
                                          v4rhosigma2tau,
                                          v4rhosigmalapl2,
                                          v4rhosigmalapltau,
                                          v4rhosigmatau2,
                                          v4rholapl3,
                                          v4rholapl2tau,
                                          v4rholapltau2,
                                          v4rhotau3,
                                          v4sigma4,
                                          v4sigma3lapl,
                                          v4sigma3tau,
                                          v4sigma2lapl2,
                                          v4sigma2lapltau,
                                          v4sigma2tau2,
                                          v4sigmalapl3,
                                          v4sigmalapl2tau,
                                          v4sigmalapltau2,
                                          v4sigmatau3,
                                          v4lapl4,
                                          v4lapl3tau,
                                          v4lapl2tau2,
                                          v4lapltau3,
                                          v4tau4);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // go through density matrices

        for (int idensity = 0; idensity < numdens_rw2; idensity++)
        {
            // compute partial contribution to Kxc matrix

            auto partial_mat_Kxc = integratePartialKxcFockForMGGA2(xcFunctional,
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
                                                                   rwdengridcube,
                                                                   rw2dengrid,
                                                                   idensity,
                                                                   timer);

            // distribute partial Kxc to full Fock matrix

            timer.start("Kxc matrix dist.");

            dftsubmat::distributeSubMatrixToFock(aoFockPointers, idensity, partial_mat_Kxc, aoinds, naos);

            timer.stop("Kxc matrix dist.");
        }

        for (int idensity = 0; idensity < numdens_rw3; idensity++)
        {
            // compute partial contribution to Lxc matrix

            auto partial_mat_Lxc = integratePartialLxcFockForMGGA(xcFunctional,
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
                                                                  v4rho4,
                                                                  v4rho3sigma,
                                                                  v4rho3lapl,
                                                                  v4rho3tau,
                                                                  v4rho2sigma2,
                                                                  v4rho2sigmalapl,
                                                                  v4rho2sigmatau,
                                                                  v4rho2lapl2,
                                                                  v4rho2lapltau,
                                                                  v4rho2tau2,
                                                                  v4rhosigma3,
                                                                  v4rhosigma2lapl,
                                                                  v4rhosigma2tau,
                                                                  v4rhosigmalapl2,
                                                                  v4rhosigmalapltau,
                                                                  v4rhosigmatau2,
                                                                  v4rholapl3,
                                                                  v4rholapl2tau,
                                                                  v4rholapltau2,
                                                                  v4rhotau3,
                                                                  v4sigma4,
                                                                  v4sigma3lapl,
                                                                  v4sigma3tau,
                                                                  v4sigma2lapl2,
                                                                  v4sigma2lapltau,
                                                                  v4sigma2tau2,
                                                                  v4sigmalapl3,
                                                                  v4sigmalapl2tau,
                                                                  v4sigmalapltau2,
                                                                  v4sigmatau3,
                                                                  v4lapl4,
                                                                  v4lapl3tau,
                                                                  v4lapl2tau2,
                                                                  v4lapltau3,
                                                                  v4tau4,
                                                                  rwdengridcube,
                                                                  rw3dengrid,
                                                                  idensity,
                                                                  timer);

            // distribute partial Lxc to full Fock matrix

            timer.start("Lxc matrix dist.");

            dftsubmat::distributeSubMatrixToFock(aoFockPointers, (idensity + numdens_rw2), partial_mat_Lxc, aoinds, naos);

            timer.stop("Lxc matrix dist.");
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
integratePartialKxcFockForMGGA2(const CXCFunctional&     xcFunctional,
                                const double*            weights,
                                const CDenseMatrix&      gtoValues,
                                const CDenseMatrix&      gtoValuesX,
                                const CDenseMatrix&      gtoValuesY,
                                const CDenseMatrix&      gtoValuesZ,
                                const double*            rhograd,
                                const double*            vsigma,
                                const double*            v2rho2,
                                const double*            v2rhosigma,
                                const double*            v2rholapl,
                                const double*            v2rhotau,
                                const double*            v2sigma2,
                                const double*            v2sigmalapl,
                                const double*            v2sigmatau,
                                const double*            v2lapl2,
                                const double*            v2lapltau,
                                const double*            v2tau2,
                                const double*            v3rho3,
                                const double*            v3rho2sigma,
                                const double*            v3rho2lapl,
                                const double*            v3rho2tau,
                                const double*            v3rhosigma2,
                                const double*            v3rhosigmalapl,
                                const double*            v3rhosigmatau,
                                const double*            v3rholapl2,
                                const double*            v3rholapltau,
                                const double*            v3rhotau2,
                                const double*            v3sigma3,
                                const double*            v3sigma2lapl,
                                const double*            v3sigma2tau,
                                const double*            v3sigmalapl2,
                                const double*            v3sigmalapltau,
                                const double*            v3sigmatau2,
                                const double*            v3lapl3,
                                const double*            v3lapl2tau,
                                const double*            v3lapltau2,
                                const double*            v3tau3,
                                const CDensityGridCubic& rwDensityGridCubic,
                                const CDensityGrid&      rw2DensityGrid,
                                const int                iFock,
                                CMultiTimer&             timer) -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    timer.start("Kxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    // pointers to perturbed densities

    auto gam    = rwDensityGridCubic.gam2(iFock);
    auto rt_gam = rwDensityGridCubic.rt_gam2(iFock);
    auto rl_gam = rwDensityGridCubic.rl_gam2(iFock);
    auto tt_gam = rwDensityGridCubic.tt_gam2(iFock);
    auto tl_gam = rwDensityGridCubic.tl_gam2(iFock);
    auto ll_gam = rwDensityGridCubic.ll_gam2(iFock);

    auto gamx    = rwDensityGridCubic.gam2X(iFock);
    auto gamy    = rwDensityGridCubic.gam2Y(iFock);
    auto gamz    = rwDensityGridCubic.gam2Z(iFock);
    auto st_gamx = rwDensityGridCubic.st_gam2X(iFock);
    auto st_gamy = rwDensityGridCubic.st_gam2Y(iFock);
    auto st_gamz = rwDensityGridCubic.st_gam2Z(iFock);
    auto sl_gamx = rwDensityGridCubic.sl_gam2X(iFock);
    auto sl_gamy = rwDensityGridCubic.sl_gam2Y(iFock);
    auto sl_gamz = rwDensityGridCubic.sl_gam2Z(iFock);

    auto gamxx = rwDensityGridCubic.gam2XX(iFock);
    auto gamxy = rwDensityGridCubic.gam2XY(iFock);
    auto gamxz = rwDensityGridCubic.gam2XZ(iFock);
    auto gamyx = rwDensityGridCubic.gam2YX(iFock);
    auto gamyy = rwDensityGridCubic.gam2YY(iFock);
    auto gamyz = rwDensityGridCubic.gam2YZ(iFock);
    auto gamzx = rwDensityGridCubic.gam2ZX(iFock);
    auto gamzy = rwDensityGridCubic.gam2ZY(iFock);
    auto gamzz = rwDensityGridCubic.gam2ZZ(iFock);

    auto rhow    = rw2DensityGrid.alphaDensity(iFock);
    auto tauw    = rw2DensityGrid.alphaDensitytau(iFock);
    auto laplw   = rw2DensityGrid.alphaDensitylapl(iFock);
    auto gradw_x = rw2DensityGrid.alphaDensityGradientX(iFock);
    auto gradw_y = rw2DensityGrid.alphaDensityGradientY(iFock);
    auto gradw_z = rw2DensityGrid.alphaDensityGradientZ(iFock);

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

                // sigma and gamma
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

                double grada_x_g = rhograd[6 * g + 0];
                double grada_y_g = rhograd[6 * g + 1];
                double grada_z_g = rhograd[6 * g + 2];

                double l2contract = grada_x_g * gradw_x[g] + grada_y_g * gradw_y[g] + grada_z_g * gradw_z[g];

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

                double rho_0 = rr * rhow[g] + rx * l2contract + rt * tauw[g];
                //+ rl * laplw[g];

                // double lap_0 =    lr * rhow[g]
                //                 + lx * l2contract
                //                 + lt * tauw[g]
                //                 + ll * laplw[g];

                double tau_0 = tr * rhow[g] + tx * l2contract + tt * tauw[g];
                //+ tl * laplw[g];

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

                // xcomp +=  grada_x_g * ( xr * rhow[g] + xt * tauw[g] + xl * laplw[g])
                xcomp += grada_x_g * (xr * rhow[g] + xt * tauw[g]) + x * gradw_x[g] + xx * l5contract_x;

                // ycomp += grada_y_g * ( xr * rhow[g] + xt * tauw[g] + xl * laplw[g])
                ycomp += grada_y_g * (xr * rhow[g] + xt * tauw[g]) + x * gradw_y[g] + xx * l5contract_y;

                // zcomp += grada_z_g * ( xr * rhow[g] + xt * tauw[g] + xl * laplw[g])
                zcomp += grada_z_g * (xr * rhow[g] + xt * tauw[g]) + x * gradw_z[g] + xx * l5contract_z;

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

auto
integratePartialLxcFockForMGGA(const CXCFunctional&     xcFunctional,
                               const double*            weights,
                               const CDenseMatrix&      gtoValues,
                               const CDenseMatrix&      gtoValuesX,
                               const CDenseMatrix&      gtoValuesY,
                               const CDenseMatrix&      gtoValuesZ,
                               const double*            rhograd,
                               const double*            vsigma,
                               const double*            v2rho2,
                               const double*            v2rhosigma,
                               const double*            v2rholapl,
                               const double*            v2rhotau,
                               const double*            v2sigma2,
                               const double*            v2sigmalapl,
                               const double*            v2sigmatau,
                               const double*            v2lapl2,
                               const double*            v2lapltau,
                               const double*            v2tau2,
                               const double*            v3rho3,
                               const double*            v3rho2sigma,
                               const double*            v3rho2lapl,
                               const double*            v3rho2tau,
                               const double*            v3rhosigma2,
                               const double*            v3rhosigmalapl,
                               const double*            v3rhosigmatau,
                               const double*            v3rholapl2,
                               const double*            v3rholapltau,
                               const double*            v3rhotau2,
                               const double*            v3sigma3,
                               const double*            v3sigma2lapl,
                               const double*            v3sigma2tau,
                               const double*            v3sigmalapl2,
                               const double*            v3sigmalapltau,
                               const double*            v3sigmatau2,
                               const double*            v3lapl3,
                               const double*            v3lapl2tau,
                               const double*            v3lapltau2,
                               const double*            v3tau3,
                               const double*            v4rho4,
                               const double*            v4rho3sigma,
                               const double*            v4rho3lapl,
                               const double*            v4rho3tau,
                               const double*            v4rho2sigma2,
                               const double*            v4rho2sigmalapl,
                               const double*            v4rho2sigmatau,
                               const double*            v4rho2lapl2,
                               const double*            v4rho2lapltau,
                               const double*            v4rho2tau2,
                               const double*            v4rhosigma3,
                               const double*            v4rhosigma2lapl,
                               const double*            v4rhosigma2tau,
                               const double*            v4rhosigmalapl2,
                               const double*            v4rhosigmalapltau,
                               const double*            v4rhosigmatau2,
                               const double*            v4rholapl3,
                               const double*            v4rholapl2tau,
                               const double*            v4rholapltau2,
                               const double*            v4rhotau3,
                               const double*            v4sigma4,
                               const double*            v4sigma3lapl,
                               const double*            v4sigma3tau,
                               const double*            v4sigma2lapl2,
                               const double*            v4sigma2lapltau,
                               const double*            v4sigma2tau2,
                               const double*            v4sigmalapl3,
                               const double*            v4sigmalapl2tau,
                               const double*            v4sigmalapltau2,
                               const double*            v4sigmatau3,
                               const double*            v4lapl4,
                               const double*            v4lapl3tau,
                               const double*            v4lapl2tau2,
                               const double*            v4lapltau3,
                               const double*            v4tau4,
                               const CDensityGridCubic& rwDensityGridCubic,
                               const CDensityGrid&      rw3DensityGrid,
                               const int                iFock,
                               CMultiTimer&             timer) -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    timer.start("Lxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    // pointers to perturbed density gradient norms

    auto gam    = rwDensityGridCubic.gam(iFock);
    auto rt_gam = rwDensityGridCubic.rt_gam(iFock);
    auto rl_gam = rwDensityGridCubic.rl_gam(iFock);
    auto tt_gam = rwDensityGridCubic.tt_gam(iFock);
    auto tl_gam = rwDensityGridCubic.tl_gam(iFock);
    auto ll_gam = rwDensityGridCubic.ll_gam(iFock);

    auto gamx    = rwDensityGridCubic.gamX(iFock);
    auto gamy    = rwDensityGridCubic.gamY(iFock);
    auto gamz    = rwDensityGridCubic.gamZ(iFock);
    auto st_gamx = rwDensityGridCubic.st_gamX(iFock);
    auto st_gamy = rwDensityGridCubic.st_gamY(iFock);
    auto st_gamz = rwDensityGridCubic.st_gamZ(iFock);
    auto sl_gamx = rwDensityGridCubic.sl_gamX(iFock);
    auto sl_gamy = rwDensityGridCubic.sl_gamY(iFock);
    auto sl_gamz = rwDensityGridCubic.sl_gamZ(iFock);

    auto gamxx = rwDensityGridCubic.gamXX(iFock);
    auto gamxy = rwDensityGridCubic.gamXY(iFock);
    auto gamxz = rwDensityGridCubic.gamXZ(iFock);
    auto gamyx = rwDensityGridCubic.gamYX(iFock);
    auto gamyy = rwDensityGridCubic.gamYY(iFock);
    auto gamyz = rwDensityGridCubic.gamYZ(iFock);
    auto gamzx = rwDensityGridCubic.gamZX(iFock);
    auto gamzy = rwDensityGridCubic.gamZY(iFock);
    auto gamzz = rwDensityGridCubic.gamZZ(iFock);

    auto pi     = rwDensityGridCubic.pi(iFock);
    auto rrt_pi = rwDensityGridCubic.rrt_pi(iFock);
    auto rrl_pi = rwDensityGridCubic.rrl_pi(iFock);
    auto rtt_pi = rwDensityGridCubic.rtt_pi(iFock);
    auto rtl_pi = rwDensityGridCubic.rtl_pi(iFock);
    auto rll_pi = rwDensityGridCubic.rll_pi(iFock);
    auto ttt_pi = rwDensityGridCubic.ttt_pi(iFock);
    auto ttl_pi = rwDensityGridCubic.ttl_pi(iFock);
    auto tll_pi = rwDensityGridCubic.tll_pi(iFock);
    auto lll_pi = rwDensityGridCubic.lll_pi(iFock);

    auto pix    = rwDensityGridCubic.piX(iFock);
    auto rt_pix = rwDensityGridCubic.rt_piX(iFock);
    auto rl_pix = rwDensityGridCubic.rl_piX(iFock);
    auto ll_pix = rwDensityGridCubic.ll_piX(iFock);
    auto tt_pix = rwDensityGridCubic.tt_piX(iFock);
    auto tl_pix = rwDensityGridCubic.tl_piX(iFock);

    auto piy    = rwDensityGridCubic.piY(iFock);
    auto rt_piy = rwDensityGridCubic.rt_piY(iFock);
    auto rl_piy = rwDensityGridCubic.rl_piY(iFock);
    auto ll_piy = rwDensityGridCubic.ll_piY(iFock);
    auto tt_piy = rwDensityGridCubic.tt_piY(iFock);
    auto tl_piy = rwDensityGridCubic.tl_piY(iFock);

    auto piz    = rwDensityGridCubic.piZ(iFock);
    auto rt_piz = rwDensityGridCubic.rt_piZ(iFock);
    auto rl_piz = rwDensityGridCubic.rl_piZ(iFock);
    auto ll_piz = rwDensityGridCubic.ll_piZ(iFock);
    auto tt_piz = rwDensityGridCubic.tt_piZ(iFock);
    auto tl_piz = rwDensityGridCubic.tl_piZ(iFock);

    auto pixx = rwDensityGridCubic.piXX(iFock);
    auto pixy = rwDensityGridCubic.piXY(iFock);
    auto pixz = rwDensityGridCubic.piXZ(iFock);
    auto piyx = rwDensityGridCubic.piYX(iFock);
    auto piyy = rwDensityGridCubic.piYY(iFock);
    auto piyz = rwDensityGridCubic.piYZ(iFock);
    auto pizx = rwDensityGridCubic.piZX(iFock);
    auto pizy = rwDensityGridCubic.piZY(iFock);
    auto pizz = rwDensityGridCubic.piZZ(iFock);

    auto l_pixx = rwDensityGridCubic.l_piXX(iFock);
    auto l_pixy = rwDensityGridCubic.l_piXY(iFock);
    auto l_pixz = rwDensityGridCubic.l_piXZ(iFock);
    auto l_piyx = rwDensityGridCubic.l_piYX(iFock);
    auto l_piyy = rwDensityGridCubic.l_piYY(iFock);
    auto l_piyz = rwDensityGridCubic.l_piYZ(iFock);
    auto l_pizx = rwDensityGridCubic.l_piZX(iFock);
    auto l_pizy = rwDensityGridCubic.l_piZY(iFock);
    auto l_pizz = rwDensityGridCubic.l_piZZ(iFock);

    auto t_pixx = rwDensityGridCubic.t_piXX(iFock);
    auto t_pixy = rwDensityGridCubic.t_piXY(iFock);
    auto t_pixz = rwDensityGridCubic.t_piXZ(iFock);
    auto t_piyx = rwDensityGridCubic.t_piYX(iFock);
    auto t_piyy = rwDensityGridCubic.t_piYY(iFock);
    auto t_piyz = rwDensityGridCubic.t_piYZ(iFock);
    auto t_pizx = rwDensityGridCubic.t_piZX(iFock);
    auto t_pizy = rwDensityGridCubic.t_piZY(iFock);
    auto t_pizz = rwDensityGridCubic.t_piZZ(iFock);

    auto pixxx = rwDensityGridCubic.piXXX(iFock);
    auto pixxy = rwDensityGridCubic.piXXY(iFock);
    auto pixxz = rwDensityGridCubic.piXXZ(iFock);
    auto pixyx = rwDensityGridCubic.piXYX(iFock);
    auto pixyy = rwDensityGridCubic.piXYY(iFock);
    auto pixyz = rwDensityGridCubic.piXYZ(iFock);
    auto pixzx = rwDensityGridCubic.piXZX(iFock);
    auto pixzy = rwDensityGridCubic.piXZY(iFock);
    auto pixzz = rwDensityGridCubic.piXZZ(iFock);
    auto piyxx = rwDensityGridCubic.piYXX(iFock);
    auto piyxy = rwDensityGridCubic.piYXY(iFock);
    auto piyxz = rwDensityGridCubic.piYXZ(iFock);
    auto piyyx = rwDensityGridCubic.piYYX(iFock);
    auto piyyy = rwDensityGridCubic.piYYY(iFock);
    auto piyyz = rwDensityGridCubic.piYYZ(iFock);
    auto piyzx = rwDensityGridCubic.piYZX(iFock);
    auto piyzy = rwDensityGridCubic.piYZY(iFock);
    auto piyzz = rwDensityGridCubic.piYZZ(iFock);
    auto pizxx = rwDensityGridCubic.piZXX(iFock);
    auto pizxy = rwDensityGridCubic.piZXY(iFock);
    auto pizxz = rwDensityGridCubic.piZXZ(iFock);
    auto pizyx = rwDensityGridCubic.piZYX(iFock);
    auto pizyy = rwDensityGridCubic.piZYY(iFock);
    auto pizyz = rwDensityGridCubic.piZYZ(iFock);
    auto pizzx = rwDensityGridCubic.piZZX(iFock);
    auto pizzy = rwDensityGridCubic.piZZY(iFock);
    auto pizzz = rwDensityGridCubic.piZZZ(iFock);

    auto rhow    = rw3DensityGrid.alphaDensity(iFock);
    auto tauw    = rw3DensityGrid.alphaDensitytau(iFock);
    auto laplw   = rw3DensityGrid.alphaDensitylapl(iFock);
    auto gradw_x = rw3DensityGrid.alphaDensityGradientX(iFock);
    auto gradw_y = rw3DensityGrid.alphaDensityGradientY(iFock);
    auto gradw_z = rw3DensityGrid.alphaDensityGradientZ(iFock);

    timer.stop("Lxc matrix prep.");

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Lxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    // LDA contribution
    CDenseMatrix mat_G(naos, npoints);

    // GGA contribution
    CDenseMatrix mat_G_gga(naos, npoints);

    // tau contribution
    CDenseMatrix mat_G_gga_x(naos, npoints);
    CDenseMatrix mat_G_gga_y(naos, npoints);
    CDenseMatrix mat_G_gga_z(naos, npoints);

    auto G_val     = mat_G.values();
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
                double w = weights[g];

                double grada_x_g = rhograd[6 * g + 0];
                double grada_y_g = rhograd[6 * g + 1];
                double grada_z_g = rhograd[6 * g + 2];

                // vx1 contractions

                double l2contract = grada_x_g * gradw_x[g] + grada_y_g * gradw_y[g] + grada_z_g * gradw_z[g];

                double l5contract_x = grada_x_g * l2contract;
                double l5contract_y = grada_y_g * l2contract;
                double l5contract_z = grada_z_g * l2contract;

                // vx2 contractions

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

                // vx3 contractions

                double c2    = grada_x_g * pix[g] + grada_y_g * piy[g] + grada_z_g * piz[g];
                double rt_c2 = grada_x_g * rt_pix[g] + grada_y_g * rt_piy[g] + grada_z_g * rt_piz[g];
                // double rl_c2 = grada_x_g * rl_pix[g] + grada_y_g * rl_piy[g] + grada_z_g * rl_piz[g];
                // double ll_c2 = grada_x_g * ll_pix[g] + grada_y_g * ll_piy[g] + grada_z_g * ll_piz[g];
                double tt_c2 = grada_x_g * tt_pix[g] + grada_y_g * tt_piy[g] + grada_z_g * tt_piz[g];
                // double tl_c2 = grada_x_g * tl_pix[g] + grada_y_g * tl_piy[g] + grada_z_g * tl_piz[g];

                double c3 = grada_x_g * grada_x_g * pixx[g] + grada_x_g * grada_y_g * pixy[g] + grada_x_g * grada_z_g * pixz[g] +
                            grada_y_g * grada_x_g * piyx[g] + grada_y_g * grada_y_g * piyy[g] + grada_y_g * grada_z_g * piyz[g] +
                            grada_z_g * grada_x_g * pizx[g] + grada_z_g * grada_y_g * pizy[g] + grada_z_g * grada_z_g * pizz[g];

                // double l_c3 =  grada_x_g * grada_x_g * l_pixx[g]
                //              + grada_x_g * grada_y_g * l_pixy[g]
                //              + grada_x_g * grada_z_g * l_pixz[g]
                //              + grada_y_g * grada_x_g * l_piyx[g]
                //              + grada_y_g * grada_y_g * l_piyy[g]
                //              + grada_y_g * grada_z_g * l_piyz[g]
                //              + grada_z_g * grada_x_g * l_pizx[g]
                //              + grada_z_g * grada_y_g * l_pizy[g]
                //              + grada_z_g * grada_z_g * l_pizz[g];

                double t_c3 = grada_x_g * grada_x_g * t_pixx[g] + grada_x_g * grada_y_g * t_pixy[g] + grada_x_g * grada_z_g * t_pixz[g] +
                              grada_y_g * grada_x_g * t_piyx[g] + grada_y_g * grada_y_g * t_piyy[g] + grada_y_g * grada_z_g * t_piyz[g] +
                              grada_z_g * grada_x_g * t_pizx[g] + grada_z_g * grada_y_g * t_pizy[g] + grada_z_g * grada_z_g * t_pizz[g];

                double c4 = pixx[g] + piyy[g] + pizz[g];
                // double l_c4 = l_pixx[g] + l_piyy[g] + l_pizz[g];
                double t_c4 = t_pixx[g] + t_piyy[g] + t_pizz[g];

                double c5_6 = grada_x_g * (pixxx[g] + pixxx[g]) + grada_x_g * (piyxy[g] + pixyy[g]) + grada_x_g * (pizxz[g] + pixzz[g]) +
                              grada_y_g * (pixyx[g] + piyxx[g]) + grada_y_g * (piyyy[g] + piyyy[g]) + grada_y_g * (pizyz[g] + piyzz[g]) +
                              grada_z_g * (pixzx[g] + pizxx[g]) + grada_z_g * (piyzy[g] + pizyy[g]) + grada_z_g * (pizzz[g] + pizzz[g]);

                double c7 = grada_x_g * grada_x_g * grada_x_g * pixxx[g] + grada_x_g * grada_x_g * grada_y_g * pixxy[g] +
                            grada_x_g * grada_x_g * grada_z_g * pixxz[g] + grada_x_g * grada_y_g * grada_x_g * pixyx[g] +
                            grada_x_g * grada_y_g * grada_y_g * pixyy[g] + grada_x_g * grada_y_g * grada_z_g * pixyz[g] +
                            grada_x_g * grada_z_g * grada_x_g * pixzx[g] + grada_x_g * grada_z_g * grada_y_g * pixzy[g] +
                            grada_x_g * grada_z_g * grada_z_g * pixzz[g] + grada_y_g * grada_x_g * grada_x_g * piyxx[g] +
                            grada_y_g * grada_x_g * grada_y_g * piyxy[g] + grada_y_g * grada_x_g * grada_z_g * piyxz[g] +
                            grada_y_g * grada_y_g * grada_x_g * piyyx[g] + grada_y_g * grada_y_g * grada_y_g * piyyy[g] +
                            grada_y_g * grada_y_g * grada_z_g * piyyz[g] + grada_y_g * grada_z_g * grada_x_g * piyzx[g] +
                            grada_y_g * grada_z_g * grada_y_g * piyzy[g] + grada_y_g * grada_z_g * grada_z_g * piyzz[g] +
                            grada_z_g * grada_x_g * grada_x_g * pizxx[g] + grada_z_g * grada_x_g * grada_y_g * pizxy[g] +
                            grada_z_g * grada_x_g * grada_z_g * pizxz[g] + grada_z_g * grada_y_g * grada_x_g * pizyx[g] +
                            grada_z_g * grada_y_g * grada_y_g * pizyy[g] + grada_z_g * grada_y_g * grada_z_g * pizyz[g] +
                            grada_z_g * grada_z_g * grada_x_g * pizzx[g] + grada_z_g * grada_z_g * grada_y_g * pizzy[g] +
                            grada_z_g * grada_z_g * grada_z_g * pizzz[g];

                double c8 = grada_x_g * pixxx[g] + grada_y_g * pixxy[g] + grada_z_g * pixxz[g] + grada_x_g * piyyx[g] + grada_y_g * piyyy[g] +
                            grada_z_g * piyyz[g] + grada_x_g * pizzx[g] + grada_y_g * pizzy[g] + grada_z_g * pizzz[g];

                // double c9_x = grada_x_g * pi[g];
                // double c9_y = grada_y_g * pi[g];
                // double c9_z = grada_z_g * pi[g];

                // double c10_x = pix[g];
                // double c10_y = piy[g];
                // double c10_z = piz[g];

                // double c11_x = c2 * grada_x_g;
                // double c11_y = c2 * grada_y_g;
                // double c11_z = c2 * grada_z_g;

                double c12_c14_x = grada_x_g * (pixx[g] + pixx[g]) + grada_y_g * (pixy[g] + piyx[g]) + grada_z_g * (pixz[g] + pizx[g]);

                // double l_c12_c14_x =   grada_x_g * (l_pixx[g] + l_pixx[g])
                //                      + grada_y_g * (l_pixy[g] + l_piyx[g])
                //                      + grada_z_g * (l_pixz[g] + l_pizx[g]);

                double t_c12_c14_x = grada_x_g * (t_pixx[g] + t_pixx[g]) + grada_y_g * (t_pixy[g] + t_piyx[g]) + grada_z_g * (t_pixz[g] + t_pizx[g]);

                double c12_c14_y = grada_x_g * (piyx[g] + pixy[g]) + grada_y_g * (piyy[g] + piyy[g]) + grada_z_g * (piyz[g] + pizy[g]);

                // double l_c12_c14_y=  grada_x_g * (l_piyx[g] + l_pixy[g])
                //                    + grada_y_g * (l_piyy[g] + l_piyy[g])
                //                    + grada_z_g * (l_piyz[g] + l_pizy[g]);

                double t_c12_c14_y = grada_x_g * (t_piyx[g] + t_pixy[g]) + grada_y_g * (t_piyy[g] + t_piyy[g]) + grada_z_g * (t_piyz[g] + t_pizy[g]);

                double c12_c14_z = grada_x_g * (pizx[g] + pixz[g]) + grada_y_g * (pizy[g] + piyz[g]) + grada_z_g * (pizz[g] + pizz[g]);

                // double l_c12_c14_z=  grada_x_g * (l_pizx[g] + l_pixz[g])
                //                    + grada_y_g * (l_pizy[g] + l_piyz[g])
                //                    + grada_z_g * (l_pizz[g] + l_pizz[g]);

                double t_c12_c14_z = grada_x_g * (t_pizx[g] + t_pixz[g]) + grada_y_g * (t_pizy[g] + t_piyz[g]) + grada_z_g * (t_pizz[g] + t_pizz[g]);

                double c13 = grada_x_g * grada_x_g * pixx[g] + grada_x_g * grada_y_g * pixy[g] + grada_x_g * grada_z_g * pixz[g] +
                             grada_y_g * grada_x_g * piyx[g] + grada_y_g * grada_y_g * piyy[g] + grada_y_g * grada_z_g * piyz[g] +
                             grada_z_g * grada_x_g * pizx[g] + grada_z_g * grada_y_g * pizy[g] + grada_z_g * grada_z_g * pizz[g];

                // double l_c13 = grada_x_g * grada_x_g * l_pixx[g]
                //              + grada_x_g * grada_y_g * l_pixy[g]
                //              + grada_x_g * grada_z_g * l_pixz[g]
                //              + grada_y_g * grada_x_g * l_piyx[g]
                //              + grada_y_g * grada_y_g * l_piyy[g]
                //              + grada_y_g * grada_z_g * l_piyz[g]
                //              + grada_z_g * grada_x_g * l_pizx[g]
                //              + grada_z_g * grada_y_g * l_pizy[g]
                //              + grada_z_g * grada_z_g * l_pizz[g];

                double t_c13 = grada_x_g * grada_x_g * t_pixx[g] + grada_x_g * grada_y_g * t_pixy[g] + grada_x_g * grada_z_g * t_pixz[g] +
                               grada_y_g * grada_x_g * t_piyx[g] + grada_y_g * grada_y_g * t_piyy[g] + grada_y_g * grada_z_g * t_piyz[g] +
                               grada_z_g * grada_x_g * t_pizx[g] + grada_z_g * grada_y_g * t_pizy[g] + grada_z_g * grada_z_g * t_pizz[g];

                double c13_x = c13 * grada_x_g;
                // double l_c13_x = l_c13 * grada_x_g;
                double t_c13_x = t_c13 * grada_x_g;

                double c13_y = c13 * grada_y_g;
                // double l_c13_y = l_c13 * grada_y_g;
                double t_c13_y = t_c13 * grada_y_g;

                double c13_z = c13 * grada_z_g;
                // double l_c13_z = l_c13 * grada_z_g;
                double t_c13_z = t_c13 * grada_z_g;

                double c15_x = grada_x_g * c4;
                // double l_c15_x = grada_x_g * l_c4;
                double t_c15_x = grada_x_g * t_c4;

                double c15_y = grada_y_g * c4;
                // double l_c15_y = grada_y_g * l_c4;
                double t_c15_y = grada_y_g * t_c4;

                double c15_z = grada_z_g * c4;
                // double l_c15_z = grada_z_g * l_c4;
                double t_c15_z = grada_z_g * t_c4;

                double c16_19_22_x = grada_x_g * grada_x_g * pixxx[g] + grada_x_g * grada_x_g * pixxx[g] + grada_x_g * grada_x_g * pixxx[g] +
                                     grada_x_g * grada_y_g * pixxy[g] + grada_x_g * grada_y_g * pixyx[g] + grada_x_g * grada_y_g * pixxy[g] +
                                     grada_x_g * grada_z_g * pixxz[g] + grada_x_g * grada_z_g * pixzx[g] + grada_x_g * grada_z_g * pixxz[g] +
                                     grada_y_g * grada_x_g * pixyx[g] + grada_y_g * grada_x_g * piyxx[g] + grada_y_g * grada_x_g * piyxx[g] +
                                     grada_y_g * grada_y_g * pixyy[g] + grada_y_g * grada_y_g * piyyx[g] + grada_y_g * grada_y_g * piyxy[g] +
                                     grada_y_g * grada_z_g * pixyz[g] + grada_y_g * grada_z_g * piyzx[g] + grada_y_g * grada_z_g * piyxz[g] +
                                     grada_z_g * grada_x_g * pixzx[g] + grada_z_g * grada_x_g * pizxx[g] + grada_z_g * grada_x_g * pizxx[g] +
                                     grada_z_g * grada_y_g * pixzy[g] + grada_z_g * grada_y_g * pizyx[g] + grada_z_g * grada_y_g * pizxy[g] +
                                     grada_z_g * grada_z_g * pixzz[g] + grada_z_g * grada_z_g * pizzx[g] + grada_z_g * grada_z_g * pizxz[g];

                double c16_19_22_y = grada_x_g * grada_x_g * piyxx[g] + grada_x_g * grada_x_g * pixxy[g] + grada_x_g * grada_x_g * pixyx[g] +
                                     grada_x_g * grada_y_g * piyxy[g] + grada_x_g * grada_y_g * pixyy[g] + grada_x_g * grada_y_g * pixyy[g] +
                                     grada_x_g * grada_z_g * piyxz[g] + grada_x_g * grada_z_g * pixzy[g] + grada_x_g * grada_z_g * pixyz[g] +
                                     grada_y_g * grada_x_g * piyyx[g] + grada_y_g * grada_x_g * piyxy[g] + grada_y_g * grada_x_g * piyyx[g] +
                                     grada_y_g * grada_y_g * piyyy[g] + grada_y_g * grada_y_g * piyyy[g] + grada_y_g * grada_y_g * piyyy[g] +
                                     grada_y_g * grada_z_g * piyyz[g] + grada_y_g * grada_z_g * piyzy[g] + grada_y_g * grada_z_g * piyyz[g] +
                                     grada_z_g * grada_x_g * piyzx[g] + grada_z_g * grada_x_g * pizxy[g] + grada_z_g * grada_x_g * pizyx[g] +
                                     grada_z_g * grada_y_g * piyzy[g] + grada_z_g * grada_y_g * pizyy[g] + grada_z_g * grada_y_g * pizyy[g] +
                                     grada_z_g * grada_z_g * piyzz[g] + grada_z_g * grada_z_g * pizzy[g] + grada_z_g * grada_z_g * pizyz[g];

                double c16_19_22_z = grada_x_g * grada_x_g * pizxx[g] + grada_x_g * grada_x_g * pixxz[g] + grada_x_g * grada_x_g * pixzx[g] +
                                     grada_x_g * grada_y_g * pizxy[g] + grada_x_g * grada_y_g * pixyz[g] + grada_x_g * grada_y_g * pixzy[g] +
                                     grada_x_g * grada_z_g * pizxz[g] + grada_x_g * grada_z_g * pixzz[g] + grada_x_g * grada_z_g * pixzz[g] +
                                     grada_y_g * grada_x_g * pizyx[g] + grada_y_g * grada_x_g * piyxz[g] + grada_y_g * grada_x_g * piyzx[g] +
                                     grada_y_g * grada_y_g * pizyy[g] + grada_y_g * grada_y_g * piyyz[g] + grada_y_g * grada_y_g * piyzy[g] +
                                     grada_y_g * grada_z_g * pizyz[g] + grada_y_g * grada_z_g * piyzz[g] + grada_y_g * grada_z_g * piyzz[g] +
                                     grada_z_g * grada_x_g * pizzx[g] + grada_z_g * grada_x_g * pizxz[g] + grada_z_g * grada_x_g * pizzx[g] +
                                     grada_z_g * grada_y_g * pizzy[g] + grada_z_g * grada_y_g * pizyz[g] + grada_z_g * grada_y_g * pizzy[g] +
                                     grada_z_g * grada_z_g * pizzz[g] + grada_z_g * grada_z_g * pizzz[g] + grada_z_g * grada_z_g * pizzz[g];

                double c17_24_25_x = pixxx[g] + pixxx[g] + pixxx[g] + pixyy[g] + piyxy[g] + piyyx[g] + pixzz[g] + pizxz[g] + pizzx[g];

                double c17_24_25_y = piyxx[g] + pixyx[g] + pixxy[g] + piyyy[g] + piyyy[g] + piyyy[g] + piyzz[g] + pizyz[g] + pizzy[g];

                double c17_24_25_z = pizxx[g] + pixzx[g] + pixxz[g] + pizyy[g] + piyzy[g] + piyyz[g] + pizzz[g] + pizzz[g] + pizzz[g];

                double c18_x = c7 * grada_x_g;
                double c18_y = c7 * grada_y_g;
                double c18_z = c7 * grada_z_g;

                double c20_21_23 = grada_x_g * (pixxx[g] + pixxx[g] + pixxx[g]) + grada_x_g * (piyxy[g] + pixyy[g] + piyyx[g]) +
                                   grada_x_g * (pizxz[g] + pixzz[g] + pizzx[g]) + grada_y_g * (pixyx[g] + piyxx[g] + pixxy[g]) +
                                   grada_y_g * (piyyy[g] + piyyy[g] + piyyy[g]) + grada_y_g * (pizyz[g] + piyzz[g] + pizzy[g]) +
                                   grada_z_g * (pixzx[g] + pizxx[g] + pixxz[g]) + grada_z_g * (piyzy[g] + pizyy[g] + piyyz[g]) +
                                   grada_z_g * (pizzz[g] + pizzz[g] + pizzz[g]);

                double c20_21_23_x = grada_x_g * c20_21_23;
                double c20_21_23_y = grada_y_g * c20_21_23;
                double c20_21_23_z = grada_z_g * c20_21_23;

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

                // fourth-order

                auto v4rho4_aaaa = v4rho4[dim->v4rho4 * g + 0];
                auto v4rho4_aaab = v4rho4[dim->v4rho4 * g + 1];
                auto v4rho4_aabb = v4rho4[dim->v4rho4 * g + 2];
                auto v4rho4_abbb = v4rho4[dim->v4rho4 * g + 3];

                auto v4rho3sigma_aaaa = v4rho3sigma[dim->v4rho3sigma * g + 0];
                auto v4rho3sigma_aaac = v4rho3sigma[dim->v4rho3sigma * g + 1];
                auto v4rho3sigma_aaab = v4rho3sigma[dim->v4rho3sigma * g + 2];
                auto v4rho3sigma_aaba = v4rho3sigma[dim->v4rho3sigma * g + 3];
                auto v4rho3sigma_aabc = v4rho3sigma[dim->v4rho3sigma * g + 4];
                auto v4rho3sigma_aabb = v4rho3sigma[dim->v4rho3sigma * g + 5];
                auto v4rho3sigma_abba = v4rho3sigma[dim->v4rho3sigma * g + 6];
                auto v4rho3sigma_abbc = v4rho3sigma[dim->v4rho3sigma * g + 7];
                auto v4rho3sigma_abbb = v4rho3sigma[dim->v4rho3sigma * g + 8];
                auto v4rho3sigma_bbba = v4rho3sigma[dim->v4rho3sigma * g + 9];
                auto v4rho3sigma_bbbc = v4rho3sigma[dim->v4rho3sigma * g + 10];

                // auto v4rho3lapl_aaaa = v4rho3lapl[dim->v4rho3lapl * g + 0];
                // auto v4rho3lapl_aaab = v4rho3lapl[dim->v4rho3lapl * g + 1];
                // auto v4rho3lapl_aaba = v4rho3lapl[dim->v4rho3lapl * g + 2];
                // auto v4rho3lapl_aabb = v4rho3lapl[dim->v4rho3lapl * g + 3];
                // auto v4rho3lapl_abba = v4rho3lapl[dim->v4rho3lapl * g + 4];
                // auto v4rho3lapl_abbb = v4rho3lapl[dim->v4rho3lapl * g + 5];
                // auto v4rho3lapl_bbba = v4rho3lapl[dim->v4rho3lapl * g + 6];

                auto v4rho3tau_aaaa = v4rho3tau[dim->v4rho3tau * g + 0];
                auto v4rho3tau_aaab = v4rho3tau[dim->v4rho3tau * g + 1];
                auto v4rho3tau_aaba = v4rho3tau[dim->v4rho3tau * g + 2];
                auto v4rho3tau_aabb = v4rho3tau[dim->v4rho3tau * g + 3];
                auto v4rho3tau_abba = v4rho3tau[dim->v4rho3tau * g + 4];
                auto v4rho3tau_abbb = v4rho3tau[dim->v4rho3tau * g + 5];
                auto v4rho3tau_bbba = v4rho3tau[dim->v4rho3tau * g + 6];

                auto v4rho2sigma2_aaaa = v4rho2sigma2[dim->v4rho2sigma2 * g + 0];
                auto v4rho2sigma2_aaac = v4rho2sigma2[dim->v4rho2sigma2 * g + 1];
                auto v4rho2sigma2_aaab = v4rho2sigma2[dim->v4rho2sigma2 * g + 2];
                auto v4rho2sigma2_aacc = v4rho2sigma2[dim->v4rho2sigma2 * g + 3];
                auto v4rho2sigma2_aacb = v4rho2sigma2[dim->v4rho2sigma2 * g + 4];
                auto v4rho2sigma2_aabb = v4rho2sigma2[dim->v4rho2sigma2 * g + 5];
                auto v4rho2sigma2_abaa = v4rho2sigma2[dim->v4rho2sigma2 * g + 6];
                auto v4rho2sigma2_abac = v4rho2sigma2[dim->v4rho2sigma2 * g + 7];
                auto v4rho2sigma2_abab = v4rho2sigma2[dim->v4rho2sigma2 * g + 8];
                auto v4rho2sigma2_abcc = v4rho2sigma2[dim->v4rho2sigma2 * g + 9];
                auto v4rho2sigma2_abcb = v4rho2sigma2[dim->v4rho2sigma2 * g + 10];
                auto v4rho2sigma2_abbb = v4rho2sigma2[dim->v4rho2sigma2 * g + 11];
                auto v4rho2sigma2_bbaa = v4rho2sigma2[dim->v4rho2sigma2 * g + 12];
                auto v4rho2sigma2_bbac = v4rho2sigma2[dim->v4rho2sigma2 * g + 13];
                auto v4rho2sigma2_bbab = v4rho2sigma2[dim->v4rho2sigma2 * g + 14];
                auto v4rho2sigma2_bbcc = v4rho2sigma2[dim->v4rho2sigma2 * g + 15];
                auto v4rho2sigma2_bbcb = v4rho2sigma2[dim->v4rho2sigma2 * g + 16];

                // auto v4rho2sigmalapl_aaaa = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 0];
                // auto v4rho2sigmalapl_aaab = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 1];
                // auto v4rho2sigmalapl_aaca = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 2];
                // auto v4rho2sigmalapl_aacb = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 3];
                // auto v4rho2sigmalapl_aaba = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 4];
                // auto v4rho2sigmalapl_aabb = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 5];
                // auto v4rho2sigmalapl_abaa = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 6];
                // auto v4rho2sigmalapl_abab = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 7];
                // auto v4rho2sigmalapl_abca = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 8];
                // auto v4rho2sigmalapl_abcb = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 9];
                // auto v4rho2sigmalapl_abba = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 10];
                // auto v4rho2sigmalapl_abbb = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 11];
                // auto v4rho2sigmalapl_bbaa = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 12];
                // auto v4rho2sigmalapl_bbab = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 13];
                // auto v4rho2sigmalapl_bbca = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 14];
                // auto v4rho2sigmalapl_bbcb = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 15];
                // auto v4rho2sigmalapl_bbba = v4rho2sigmalapl[dim->v4rho2sigmalapl * g + 16];

                auto v4rho2sigmatau_aaaa = v4rho2sigmatau[dim->v4rho2sigmatau * g + 0];
                auto v4rho2sigmatau_aaab = v4rho2sigmatau[dim->v4rho2sigmatau * g + 1];
                auto v4rho2sigmatau_aaca = v4rho2sigmatau[dim->v4rho2sigmatau * g + 2];
                auto v4rho2sigmatau_aacb = v4rho2sigmatau[dim->v4rho2sigmatau * g + 3];
                auto v4rho2sigmatau_aaba = v4rho2sigmatau[dim->v4rho2sigmatau * g + 4];
                auto v4rho2sigmatau_aabb = v4rho2sigmatau[dim->v4rho2sigmatau * g + 5];
                auto v4rho2sigmatau_abaa = v4rho2sigmatau[dim->v4rho2sigmatau * g + 6];
                auto v4rho2sigmatau_abab = v4rho2sigmatau[dim->v4rho2sigmatau * g + 7];
                auto v4rho2sigmatau_abca = v4rho2sigmatau[dim->v4rho2sigmatau * g + 8];
                auto v4rho2sigmatau_abcb = v4rho2sigmatau[dim->v4rho2sigmatau * g + 9];
                auto v4rho2sigmatau_abba = v4rho2sigmatau[dim->v4rho2sigmatau * g + 10];
                auto v4rho2sigmatau_abbb = v4rho2sigmatau[dim->v4rho2sigmatau * g + 11];
                auto v4rho2sigmatau_bbaa = v4rho2sigmatau[dim->v4rho2sigmatau * g + 12];
                auto v4rho2sigmatau_bbab = v4rho2sigmatau[dim->v4rho2sigmatau * g + 13];
                auto v4rho2sigmatau_bbca = v4rho2sigmatau[dim->v4rho2sigmatau * g + 14];
                auto v4rho2sigmatau_bbcb = v4rho2sigmatau[dim->v4rho2sigmatau * g + 15];
                auto v4rho2sigmatau_bbba = v4rho2sigmatau[dim->v4rho2sigmatau * g + 16];

                // auto v4rho2lapl2_aaaa = v4rho2lapl2[dim->v4rho2lapl2 * g + 0];
                // auto v4rho2lapl2_aaab = v4rho2lapl2[dim->v4rho2lapl2 * g + 1];
                // auto v4rho2lapl2_aabb = v4rho2lapl2[dim->v4rho2lapl2 * g + 2];
                // auto v4rho2lapl2_abaa = v4rho2lapl2[dim->v4rho2lapl2 * g + 3];
                // auto v4rho2lapl2_abab = v4rho2lapl2[dim->v4rho2lapl2 * g + 4];
                // auto v4rho2lapl2_abbb = v4rho2lapl2[dim->v4rho2lapl2 * g + 5];
                // auto v4rho2lapl2_bbaa = v4rho2lapl2[dim->v4rho2lapl2 * g + 6];
                // auto v4rho2lapl2_bbab = v4rho2lapl2[dim->v4rho2lapl2 * g + 7];

                // auto v4rho2lapltau_aaaa = v4rho2lapltau[dim->v4rho2lapltau * g + 0];
                // auto v4rho2lapltau_aaab = v4rho2lapltau[dim->v4rho2lapltau * g + 1];
                // auto v4rho2lapltau_aaba = v4rho2lapltau[dim->v4rho2lapltau * g + 2];
                // auto v4rho2lapltau_aabb = v4rho2lapltau[dim->v4rho2lapltau * g + 3];
                // auto v4rho2lapltau_abaa = v4rho2lapltau[dim->v4rho2lapltau * g + 4];
                // auto v4rho2lapltau_abab = v4rho2lapltau[dim->v4rho2lapltau * g + 5];
                // auto v4rho2lapltau_abba = v4rho2lapltau[dim->v4rho2lapltau * g + 6];
                // auto v4rho2lapltau_abbb = v4rho2lapltau[dim->v4rho2lapltau * g + 7];
                // auto v4rho2lapltau_bbaa = v4rho2lapltau[dim->v4rho2lapltau * g + 8];
                // auto v4rho2lapltau_bbab = v4rho2lapltau[dim->v4rho2lapltau * g + 9];
                // auto v4rho2lapltau_bbba = v4rho2lapltau[dim->v4rho2lapltau * g + 10];

                auto v4rho2tau2_aaaa = v4rho2tau2[dim->v4rho2tau2 * g + 0];
                auto v4rho2tau2_aaab = v4rho2tau2[dim->v4rho2tau2 * g + 1];
                auto v4rho2tau2_aabb = v4rho2tau2[dim->v4rho2tau2 * g + 2];
                auto v4rho2tau2_abaa = v4rho2tau2[dim->v4rho2tau2 * g + 3];
                auto v4rho2tau2_abab = v4rho2tau2[dim->v4rho2tau2 * g + 4];
                auto v4rho2tau2_abbb = v4rho2tau2[dim->v4rho2tau2 * g + 5];
                auto v4rho2tau2_bbaa = v4rho2tau2[dim->v4rho2tau2 * g + 6];
                auto v4rho2tau2_bbab = v4rho2tau2[dim->v4rho2tau2 * g + 7];

                auto v4rhosigma3_aaaa = v4rhosigma3[dim->v4rhosigma3 * g + 0];
                auto v4rhosigma3_aaac = v4rhosigma3[dim->v4rhosigma3 * g + 1];
                auto v4rhosigma3_aaab = v4rhosigma3[dim->v4rhosigma3 * g + 2];
                auto v4rhosigma3_aacc = v4rhosigma3[dim->v4rhosigma3 * g + 3];
                auto v4rhosigma3_aacb = v4rhosigma3[dim->v4rhosigma3 * g + 4];
                auto v4rhosigma3_aabb = v4rhosigma3[dim->v4rhosigma3 * g + 5];
                auto v4rhosigma3_accc = v4rhosigma3[dim->v4rhosigma3 * g + 6];
                auto v4rhosigma3_accb = v4rhosigma3[dim->v4rhosigma3 * g + 7];
                auto v4rhosigma3_acbb = v4rhosigma3[dim->v4rhosigma3 * g + 8];
                auto v4rhosigma3_abbb = v4rhosigma3[dim->v4rhosigma3 * g + 9];
                auto v4rhosigma3_baaa = v4rhosigma3[dim->v4rhosigma3 * g + 10];
                auto v4rhosigma3_baac = v4rhosigma3[dim->v4rhosigma3 * g + 11];
                auto v4rhosigma3_baab = v4rhosigma3[dim->v4rhosigma3 * g + 12];
                auto v4rhosigma3_bacc = v4rhosigma3[dim->v4rhosigma3 * g + 13];
                auto v4rhosigma3_bacb = v4rhosigma3[dim->v4rhosigma3 * g + 14];
                auto v4rhosigma3_babb = v4rhosigma3[dim->v4rhosigma3 * g + 15];
                auto v4rhosigma3_bccc = v4rhosigma3[dim->v4rhosigma3 * g + 16];
                auto v4rhosigma3_bccb = v4rhosigma3[dim->v4rhosigma3 * g + 17];
                auto v4rhosigma3_bcbb = v4rhosigma3[dim->v4rhosigma3 * g + 18];

                // auto v4rhosigma2lapl_aaaa = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 0];
                // auto v4rhosigma2lapl_aaab = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 1];
                // auto v4rhosigma2lapl_aaca = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 2];
                // auto v4rhosigma2lapl_aacb = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 3];
                // auto v4rhosigma2lapl_aaba = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 4];
                // auto v4rhosigma2lapl_aabb = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 5];
                // auto v4rhosigma2lapl_acca = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 6];
                // auto v4rhosigma2lapl_accb = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 7];
                // auto v4rhosigma2lapl_acba = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 8];
                // auto v4rhosigma2lapl_acbb = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 9];
                // auto v4rhosigma2lapl_abba = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 10];
                // auto v4rhosigma2lapl_abbb = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 11];
                // auto v4rhosigma2lapl_baaa = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 12];
                // auto v4rhosigma2lapl_baab = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 13];
                // auto v4rhosigma2lapl_baca = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 14];
                // auto v4rhosigma2lapl_bacb = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 15];
                // auto v4rhosigma2lapl_baba = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 16];
                // auto v4rhosigma2lapl_babb = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 17];
                // auto v4rhosigma2lapl_bcca = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 18];
                // auto v4rhosigma2lapl_bccb = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 19];
                // auto v4rhosigma2lapl_bcba = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 20];
                // auto v4rhosigma2lapl_bcbb = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 21];
                // auto v4rhosigma2lapl_bbba = v4rhosigma2lapl[dim->v4rhosigma2lapl * g + 22];

                auto v4rhosigma2tau_aaaa = v4rhosigma2tau[dim->v4rhosigma2tau * g + 0];
                auto v4rhosigma2tau_aaab = v4rhosigma2tau[dim->v4rhosigma2tau * g + 1];
                auto v4rhosigma2tau_aaca = v4rhosigma2tau[dim->v4rhosigma2tau * g + 2];
                auto v4rhosigma2tau_aacb = v4rhosigma2tau[dim->v4rhosigma2tau * g + 3];
                auto v4rhosigma2tau_aaba = v4rhosigma2tau[dim->v4rhosigma2tau * g + 4];
                auto v4rhosigma2tau_aabb = v4rhosigma2tau[dim->v4rhosigma2tau * g + 5];
                auto v4rhosigma2tau_acca = v4rhosigma2tau[dim->v4rhosigma2tau * g + 6];
                auto v4rhosigma2tau_accb = v4rhosigma2tau[dim->v4rhosigma2tau * g + 7];
                auto v4rhosigma2tau_acba = v4rhosigma2tau[dim->v4rhosigma2tau * g + 8];
                auto v4rhosigma2tau_acbb = v4rhosigma2tau[dim->v4rhosigma2tau * g + 9];
                auto v4rhosigma2tau_abba = v4rhosigma2tau[dim->v4rhosigma2tau * g + 10];
                auto v4rhosigma2tau_abbb = v4rhosigma2tau[dim->v4rhosigma2tau * g + 11];
                auto v4rhosigma2tau_baaa = v4rhosigma2tau[dim->v4rhosigma2tau * g + 12];
                auto v4rhosigma2tau_baab = v4rhosigma2tau[dim->v4rhosigma2tau * g + 13];
                auto v4rhosigma2tau_baca = v4rhosigma2tau[dim->v4rhosigma2tau * g + 14];
                auto v4rhosigma2tau_bacb = v4rhosigma2tau[dim->v4rhosigma2tau * g + 15];
                auto v4rhosigma2tau_baba = v4rhosigma2tau[dim->v4rhosigma2tau * g + 16];
                auto v4rhosigma2tau_babb = v4rhosigma2tau[dim->v4rhosigma2tau * g + 17];
                auto v4rhosigma2tau_bcca = v4rhosigma2tau[dim->v4rhosigma2tau * g + 18];
                auto v4rhosigma2tau_bccb = v4rhosigma2tau[dim->v4rhosigma2tau * g + 19];
                auto v4rhosigma2tau_bcba = v4rhosigma2tau[dim->v4rhosigma2tau * g + 20];
                auto v4rhosigma2tau_bcbb = v4rhosigma2tau[dim->v4rhosigma2tau * g + 21];
                auto v4rhosigma2tau_bbba = v4rhosigma2tau[dim->v4rhosigma2tau * g + 22];

                // auto v4rhosigmalapl2_aaaa = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 0];
                // auto v4rhosigmalapl2_aaab = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 1];
                // auto v4rhosigmalapl2_aabb = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 2];
                // auto v4rhosigmalapl2_acaa = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 3];
                // auto v4rhosigmalapl2_acab = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 4];
                // auto v4rhosigmalapl2_acbb = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 5];
                // auto v4rhosigmalapl2_abaa = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 6];
                // auto v4rhosigmalapl2_abab = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 7];
                // auto v4rhosigmalapl2_abbb = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 8];
                // auto v4rhosigmalapl2_baaa = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 9];
                // auto v4rhosigmalapl2_baab = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 10];
                // auto v4rhosigmalapl2_babb = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 11];
                // auto v4rhosigmalapl2_bcaa = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 12];
                // auto v4rhosigmalapl2_bcab = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 13];
                // auto v4rhosigmalapl2_bcbb = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 14];
                // auto v4rhosigmalapl2_bbaa = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 15];
                // auto v4rhosigmalapl2_bbab = v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + 16];

                // auto v4rhosigmalapltau_aaaa = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 0];
                // auto v4rhosigmalapltau_aaab = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 1];
                // auto v4rhosigmalapltau_aaba = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 2];
                // auto v4rhosigmalapltau_aabb = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 3];
                // auto v4rhosigmalapltau_acaa = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 4];
                // auto v4rhosigmalapltau_acab = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 5];
                // auto v4rhosigmalapltau_acba = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 6];
                // auto v4rhosigmalapltau_acbb = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 7];
                // auto v4rhosigmalapltau_abaa = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 8];
                // auto v4rhosigmalapltau_abab = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 9];
                // auto v4rhosigmalapltau_abba = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 10];
                // auto v4rhosigmalapltau_abbb = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 11];
                // auto v4rhosigmalapltau_baaa = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 12];
                // auto v4rhosigmalapltau_baab = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 13];
                // auto v4rhosigmalapltau_baba = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 14];
                // auto v4rhosigmalapltau_babb = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 15];
                // auto v4rhosigmalapltau_bcaa = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 16];
                // auto v4rhosigmalapltau_bcab = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 17];
                // auto v4rhosigmalapltau_bcba = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 18];
                // auto v4rhosigmalapltau_bcbb = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 19];
                // auto v4rhosigmalapltau_bbaa = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 20];
                // auto v4rhosigmalapltau_bbab = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 21];
                // auto v4rhosigmalapltau_bbba = v4rhosigmalapltau[dim->v4rhosigmalapltau * g + 22];

                auto v4rhosigmatau2_aaaa = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 0];
                auto v4rhosigmatau2_aaab = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 1];
                auto v4rhosigmatau2_aabb = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 2];
                auto v4rhosigmatau2_acaa = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 3];
                auto v4rhosigmatau2_acab = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 4];
                auto v4rhosigmatau2_acbb = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 5];
                auto v4rhosigmatau2_abaa = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 6];
                auto v4rhosigmatau2_abab = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 7];
                auto v4rhosigmatau2_abbb = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 8];
                auto v4rhosigmatau2_baaa = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 9];
                auto v4rhosigmatau2_baab = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 10];
                auto v4rhosigmatau2_babb = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 11];
                auto v4rhosigmatau2_bcaa = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 12];
                auto v4rhosigmatau2_bcab = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 13];
                auto v4rhosigmatau2_bcbb = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 14];
                auto v4rhosigmatau2_bbaa = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 15];
                auto v4rhosigmatau2_bbab = v4rhosigmatau2[dim->v4rhosigmatau2 * g + 16];

                // auto v4rholapl3_aaaa = v4rholapl3[dim->v4rholapl3 * g + 0];
                // auto v4rholapl3_aaab = v4rholapl3[dim->v4rholapl3 * g + 1];
                // auto v4rholapl3_aabb = v4rholapl3[dim->v4rholapl3 * g + 2];
                // auto v4rholapl3_abbb = v4rholapl3[dim->v4rholapl3 * g + 3];
                // auto v4rholapl3_baaa = v4rholapl3[dim->v4rholapl3 * g + 4];
                // auto v4rholapl3_baab = v4rholapl3[dim->v4rholapl3 * g + 5];
                // auto v4rholapl3_babb = v4rholapl3[dim->v4rholapl3 * g + 6];

                // auto v4rholapl2tau_aaaa = v4rholapl2tau[dim->v4rholapl2tau * g + 0];
                // auto v4rholapl2tau_aaab = v4rholapl2tau[dim->v4rholapl2tau * g + 1];
                // auto v4rholapl2tau_aaba = v4rholapl2tau[dim->v4rholapl2tau * g + 2];
                // auto v4rholapl2tau_aabb = v4rholapl2tau[dim->v4rholapl2tau * g + 3];
                // auto v4rholapl2tau_abba = v4rholapl2tau[dim->v4rholapl2tau * g + 4];
                // auto v4rholapl2tau_abbb = v4rholapl2tau[dim->v4rholapl2tau * g + 5];
                // auto v4rholapl2tau_baaa = v4rholapl2tau[dim->v4rholapl2tau * g + 6];
                // auto v4rholapl2tau_baab = v4rholapl2tau[dim->v4rholapl2tau * g + 7];
                // auto v4rholapl2tau_baba = v4rholapl2tau[dim->v4rholapl2tau * g + 8];
                // auto v4rholapl2tau_babb = v4rholapl2tau[dim->v4rholapl2tau * g + 9];
                // auto v4rholapl2tau_bbba = v4rholapl2tau[dim->v4rholapl2tau * g + 10];

                // auto v4rholapltau2_aaaa = v4rholapltau2[dim->v4rholapltau2 * g + 0];
                // auto v4rholapltau2_aaab = v4rholapltau2[dim->v4rholapltau2 * g + 1];
                // auto v4rholapltau2_aabb = v4rholapltau2[dim->v4rholapltau2 * g + 2];
                // auto v4rholapltau2_abaa = v4rholapltau2[dim->v4rholapltau2 * g + 3];
                // auto v4rholapltau2_abab = v4rholapltau2[dim->v4rholapltau2 * g + 4];
                // auto v4rholapltau2_abbb = v4rholapltau2[dim->v4rholapltau2 * g + 5];
                // auto v4rholapltau2_baaa = v4rholapltau2[dim->v4rholapltau2 * g + 6];
                // auto v4rholapltau2_baab = v4rholapltau2[dim->v4rholapltau2 * g + 7];
                // auto v4rholapltau2_babb = v4rholapltau2[dim->v4rholapltau2 * g + 8];
                // auto v4rholapltau2_bbaa = v4rholapltau2[dim->v4rholapltau2 * g + 9];
                // auto v4rholapltau2_bbab = v4rholapltau2[dim->v4rholapltau2 * g + 10];

                auto v4rhotau3_aaaa = v4rhotau3[dim->v4rhotau3 * g + 0];
                auto v4rhotau3_aaab = v4rhotau3[dim->v4rhotau3 * g + 1];
                auto v4rhotau3_aabb = v4rhotau3[dim->v4rhotau3 * g + 2];
                auto v4rhotau3_abbb = v4rhotau3[dim->v4rhotau3 * g + 3];
                auto v4rhotau3_baaa = v4rhotau3[dim->v4rhotau3 * g + 4];
                auto v4rhotau3_baab = v4rhotau3[dim->v4rhotau3 * g + 5];
                auto v4rhotau3_babb = v4rhotau3[dim->v4rhotau3 * g + 6];

                auto v4sigma4_aaaa = v4sigma4[dim->v4sigma4 * g + 0];
                auto v4sigma4_aaac = v4sigma4[dim->v4sigma4 * g + 1];
                auto v4sigma4_aaab = v4sigma4[dim->v4sigma4 * g + 2];
                auto v4sigma4_aacc = v4sigma4[dim->v4sigma4 * g + 3];
                auto v4sigma4_aacb = v4sigma4[dim->v4sigma4 * g + 4];
                auto v4sigma4_aabb = v4sigma4[dim->v4sigma4 * g + 5];
                auto v4sigma4_accc = v4sigma4[dim->v4sigma4 * g + 6];
                auto v4sigma4_accb = v4sigma4[dim->v4sigma4 * g + 7];
                auto v4sigma4_acbb = v4sigma4[dim->v4sigma4 * g + 8];
                auto v4sigma4_abbb = v4sigma4[dim->v4sigma4 * g + 9];
                auto v4sigma4_cccc = v4sigma4[dim->v4sigma4 * g + 10];
                auto v4sigma4_cccb = v4sigma4[dim->v4sigma4 * g + 11];
                auto v4sigma4_ccbb = v4sigma4[dim->v4sigma4 * g + 12];
                auto v4sigma4_cbbb = v4sigma4[dim->v4sigma4 * g + 13];

                // auto v4sigma3lapl_aaaa = v4sigma3lapl[dim->v4sigma3lapl * g + 0];
                // auto v4sigma3lapl_aaab = v4sigma3lapl[dim->v4sigma3lapl * g + 1];
                // auto v4sigma3lapl_aaca = v4sigma3lapl[dim->v4sigma3lapl * g + 2];
                // auto v4sigma3lapl_aacb = v4sigma3lapl[dim->v4sigma3lapl * g + 3];
                // auto v4sigma3lapl_aaba = v4sigma3lapl[dim->v4sigma3lapl * g + 4];
                // auto v4sigma3lapl_aabb = v4sigma3lapl[dim->v4sigma3lapl * g + 5];
                // auto v4sigma3lapl_acca = v4sigma3lapl[dim->v4sigma3lapl * g + 6];
                // auto v4sigma3lapl_accb = v4sigma3lapl[dim->v4sigma3lapl * g + 7];
                // auto v4sigma3lapl_acba = v4sigma3lapl[dim->v4sigma3lapl * g + 8];
                // auto v4sigma3lapl_acbb = v4sigma3lapl[dim->v4sigma3lapl * g + 9];
                // auto v4sigma3lapl_abba = v4sigma3lapl[dim->v4sigma3lapl * g + 10];
                // auto v4sigma3lapl_abbb = v4sigma3lapl[dim->v4sigma3lapl * g + 11];
                // auto v4sigma3lapl_ccca = v4sigma3lapl[dim->v4sigma3lapl * g + 12];
                // auto v4sigma3lapl_cccb = v4sigma3lapl[dim->v4sigma3lapl * g + 13];
                // auto v4sigma3lapl_ccba = v4sigma3lapl[dim->v4sigma3lapl * g + 14];
                // auto v4sigma3lapl_ccbb = v4sigma3lapl[dim->v4sigma3lapl * g + 15];
                // auto v4sigma3lapl_cbba = v4sigma3lapl[dim->v4sigma3lapl * g + 16];
                // auto v4sigma3lapl_cbbb = v4sigma3lapl[dim->v4sigma3lapl * g + 17];
                // auto v4sigma3lapl_bbba = v4sigma3lapl[dim->v4sigma3lapl * g + 18];

                auto v4sigma3tau_aaaa = v4sigma3tau[dim->v4sigma3tau * g + 0];
                auto v4sigma3tau_aaab = v4sigma3tau[dim->v4sigma3tau * g + 1];
                auto v4sigma3tau_aaca = v4sigma3tau[dim->v4sigma3tau * g + 2];
                auto v4sigma3tau_aacb = v4sigma3tau[dim->v4sigma3tau * g + 3];
                auto v4sigma3tau_aaba = v4sigma3tau[dim->v4sigma3tau * g + 4];
                auto v4sigma3tau_aabb = v4sigma3tau[dim->v4sigma3tau * g + 5];
                auto v4sigma3tau_acca = v4sigma3tau[dim->v4sigma3tau * g + 6];
                auto v4sigma3tau_accb = v4sigma3tau[dim->v4sigma3tau * g + 7];
                auto v4sigma3tau_acba = v4sigma3tau[dim->v4sigma3tau * g + 8];
                auto v4sigma3tau_acbb = v4sigma3tau[dim->v4sigma3tau * g + 9];
                auto v4sigma3tau_abba = v4sigma3tau[dim->v4sigma3tau * g + 10];
                auto v4sigma3tau_abbb = v4sigma3tau[dim->v4sigma3tau * g + 11];
                auto v4sigma3tau_ccca = v4sigma3tau[dim->v4sigma3tau * g + 12];
                auto v4sigma3tau_cccb = v4sigma3tau[dim->v4sigma3tau * g + 13];
                auto v4sigma3tau_ccba = v4sigma3tau[dim->v4sigma3tau * g + 14];
                auto v4sigma3tau_ccbb = v4sigma3tau[dim->v4sigma3tau * g + 15];
                auto v4sigma3tau_cbba = v4sigma3tau[dim->v4sigma3tau * g + 16];
                auto v4sigma3tau_cbbb = v4sigma3tau[dim->v4sigma3tau * g + 17];
                auto v4sigma3tau_bbba = v4sigma3tau[dim->v4sigma3tau * g + 18];

                // auto v4sigma2lapl2_aaaa = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 0];
                // auto v4sigma2lapl2_aaab = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 1];
                // auto v4sigma2lapl2_aabb = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 2];
                // auto v4sigma2lapl2_acaa = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 3];
                // auto v4sigma2lapl2_acab = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 4];
                // auto v4sigma2lapl2_acbb = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 5];
                // auto v4sigma2lapl2_abaa = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 6];
                // auto v4sigma2lapl2_abab = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 7];
                // auto v4sigma2lapl2_abbb = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 8];
                // auto v4sigma2lapl2_ccaa = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 9];
                // auto v4sigma2lapl2_ccab = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 10];
                // auto v4sigma2lapl2_ccbb = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 11];
                // auto v4sigma2lapl2_cbaa = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 12];
                // auto v4sigma2lapl2_cbab = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 13];
                // auto v4sigma2lapl2_cbbb = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 14];
                // auto v4sigma2lapl2_bbaa = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 15];
                // auto v4sigma2lapl2_bbab = v4sigma2lapl2[dim->v4sigma2lapl2 * g + 16];

                // auto v4sigma2lapltau_aaaa = v4sigma2lapltau[dim->v4sigma2lapltau * g + 0];
                // auto v4sigma2lapltau_aaab = v4sigma2lapltau[dim->v4sigma2lapltau * g + 1];
                // auto v4sigma2lapltau_aaba = v4sigma2lapltau[dim->v4sigma2lapltau * g + 2];
                // auto v4sigma2lapltau_aabb = v4sigma2lapltau[dim->v4sigma2lapltau * g + 3];
                // auto v4sigma2lapltau_acaa = v4sigma2lapltau[dim->v4sigma2lapltau * g + 4];
                // auto v4sigma2lapltau_acab = v4sigma2lapltau[dim->v4sigma2lapltau * g + 5];
                // auto v4sigma2lapltau_acba = v4sigma2lapltau[dim->v4sigma2lapltau * g + 6];
                // auto v4sigma2lapltau_acbb = v4sigma2lapltau[dim->v4sigma2lapltau * g + 7];
                // auto v4sigma2lapltau_abaa = v4sigma2lapltau[dim->v4sigma2lapltau * g + 8];
                // auto v4sigma2lapltau_abab = v4sigma2lapltau[dim->v4sigma2lapltau * g + 9];
                // auto v4sigma2lapltau_abba = v4sigma2lapltau[dim->v4sigma2lapltau * g + 10];
                // auto v4sigma2lapltau_abbb = v4sigma2lapltau[dim->v4sigma2lapltau * g + 11];
                // auto v4sigma2lapltau_ccaa = v4sigma2lapltau[dim->v4sigma2lapltau * g + 12];
                // auto v4sigma2lapltau_ccab = v4sigma2lapltau[dim->v4sigma2lapltau * g + 13];
                // auto v4sigma2lapltau_ccba = v4sigma2lapltau[dim->v4sigma2lapltau * g + 14];
                // auto v4sigma2lapltau_ccbb = v4sigma2lapltau[dim->v4sigma2lapltau * g + 15];
                // auto v4sigma2lapltau_cbaa = v4sigma2lapltau[dim->v4sigma2lapltau * g + 16];
                // auto v4sigma2lapltau_cbab = v4sigma2lapltau[dim->v4sigma2lapltau * g + 17];
                // auto v4sigma2lapltau_cbba = v4sigma2lapltau[dim->v4sigma2lapltau * g + 18];
                // auto v4sigma2lapltau_cbbb = v4sigma2lapltau[dim->v4sigma2lapltau * g + 19];
                // auto v4sigma2lapltau_bbaa = v4sigma2lapltau[dim->v4sigma2lapltau * g + 20];
                // auto v4sigma2lapltau_bbab = v4sigma2lapltau[dim->v4sigma2lapltau * g + 21];
                // auto v4sigma2lapltau_bbba = v4sigma2lapltau[dim->v4sigma2lapltau * g + 22];

                auto v4sigma2tau2_aaaa = v4sigma2tau2[dim->v4sigma2tau2 * g + 0];
                auto v4sigma2tau2_aaab = v4sigma2tau2[dim->v4sigma2tau2 * g + 1];
                auto v4sigma2tau2_aabb = v4sigma2tau2[dim->v4sigma2tau2 * g + 2];
                auto v4sigma2tau2_acaa = v4sigma2tau2[dim->v4sigma2tau2 * g + 3];
                auto v4sigma2tau2_acab = v4sigma2tau2[dim->v4sigma2tau2 * g + 4];
                auto v4sigma2tau2_acbb = v4sigma2tau2[dim->v4sigma2tau2 * g + 5];
                auto v4sigma2tau2_abaa = v4sigma2tau2[dim->v4sigma2tau2 * g + 6];
                auto v4sigma2tau2_abab = v4sigma2tau2[dim->v4sigma2tau2 * g + 7];
                auto v4sigma2tau2_abbb = v4sigma2tau2[dim->v4sigma2tau2 * g + 8];
                auto v4sigma2tau2_ccaa = v4sigma2tau2[dim->v4sigma2tau2 * g + 9];
                auto v4sigma2tau2_ccab = v4sigma2tau2[dim->v4sigma2tau2 * g + 10];
                auto v4sigma2tau2_ccbb = v4sigma2tau2[dim->v4sigma2tau2 * g + 11];
                auto v4sigma2tau2_cbaa = v4sigma2tau2[dim->v4sigma2tau2 * g + 12];
                auto v4sigma2tau2_cbab = v4sigma2tau2[dim->v4sigma2tau2 * g + 13];
                auto v4sigma2tau2_cbbb = v4sigma2tau2[dim->v4sigma2tau2 * g + 14];
                auto v4sigma2tau2_bbaa = v4sigma2tau2[dim->v4sigma2tau2 * g + 15];
                auto v4sigma2tau2_bbab = v4sigma2tau2[dim->v4sigma2tau2 * g + 16];

                // auto v4sigmalapl3_aaaa = v4sigmalapl3[dim->v4sigmalapl3 * g + 0];
                // auto v4sigmalapl3_aaab = v4sigmalapl3[dim->v4sigmalapl3 * g + 1];
                // auto v4sigmalapl3_aabb = v4sigmalapl3[dim->v4sigmalapl3 * g + 2];
                // auto v4sigmalapl3_abbb = v4sigmalapl3[dim->v4sigmalapl3 * g + 3];
                // auto v4sigmalapl3_caaa = v4sigmalapl3[dim->v4sigmalapl3 * g + 4];
                // auto v4sigmalapl3_caab = v4sigmalapl3[dim->v4sigmalapl3 * g + 5];
                // auto v4sigmalapl3_cabb = v4sigmalapl3[dim->v4sigmalapl3 * g + 6];
                // auto v4sigmalapl3_cbbb = v4sigmalapl3[dim->v4sigmalapl3 * g + 7];
                // auto v4sigmalapl3_baaa = v4sigmalapl3[dim->v4sigmalapl3 * g + 8];
                // auto v4sigmalapl3_baab = v4sigmalapl3[dim->v4sigmalapl3 * g + 9];
                // auto v4sigmalapl3_babb = v4sigmalapl3[dim->v4sigmalapl3 * g + 10];

                // auto v4sigmalapl2tau_aaaa = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 0];
                // auto v4sigmalapl2tau_aaab = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 1];
                // auto v4sigmalapl2tau_aaba = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 2];
                // auto v4sigmalapl2tau_aabb = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 3];
                // auto v4sigmalapl2tau_abba = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 4];
                // auto v4sigmalapl2tau_abbb = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 5];
                // auto v4sigmalapl2tau_caaa = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 6];
                // auto v4sigmalapl2tau_caab = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 7];
                // auto v4sigmalapl2tau_caba = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 8];
                // auto v4sigmalapl2tau_cabb = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 9];
                // auto v4sigmalapl2tau_cbba = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 10];
                // auto v4sigmalapl2tau_cbbb = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 11];
                // auto v4sigmalapl2tau_baaa = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 12];
                // auto v4sigmalapl2tau_baab = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 13];
                // auto v4sigmalapl2tau_baba = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 14];
                // auto v4sigmalapl2tau_babb = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 15];
                // auto v4sigmalapl2tau_bbba = v4sigmalapl2tau[dim->v4sigmalapl2tau * g + 16];

                // auto v4sigmalapltau2_aaaa = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 0];
                // auto v4sigmalapltau2_aaab = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 1];
                // auto v4sigmalapltau2_aabb = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 2];
                // auto v4sigmalapltau2_abaa = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 3];
                // auto v4sigmalapltau2_abab = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 4];
                // auto v4sigmalapltau2_abbb = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 5];
                // auto v4sigmalapltau2_caaa = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 6];
                // auto v4sigmalapltau2_caab = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 7];
                // auto v4sigmalapltau2_cabb = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 8];
                // auto v4sigmalapltau2_cbaa = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 9];
                // auto v4sigmalapltau2_cbab = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 10];
                // auto v4sigmalapltau2_cbbb = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 11];
                // auto v4sigmalapltau2_baaa = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 12];
                // auto v4sigmalapltau2_baab = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 13];
                // auto v4sigmalapltau2_babb = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 14];
                // auto v4sigmalapltau2_bbaa = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 15];
                // auto v4sigmalapltau2_bbab = v4sigmalapltau2[dim->v4sigmalapltau2 * g + 16];

                auto v4sigmatau3_aaaa = v4sigmatau3[dim->v4sigmatau3 * g + 0];
                auto v4sigmatau3_aaab = v4sigmatau3[dim->v4sigmatau3 * g + 1];
                auto v4sigmatau3_aabb = v4sigmatau3[dim->v4sigmatau3 * g + 2];
                auto v4sigmatau3_abbb = v4sigmatau3[dim->v4sigmatau3 * g + 3];
                auto v4sigmatau3_caaa = v4sigmatau3[dim->v4sigmatau3 * g + 4];
                auto v4sigmatau3_caab = v4sigmatau3[dim->v4sigmatau3 * g + 5];
                auto v4sigmatau3_cabb = v4sigmatau3[dim->v4sigmatau3 * g + 6];
                auto v4sigmatau3_cbbb = v4sigmatau3[dim->v4sigmatau3 * g + 7];
                auto v4sigmatau3_baaa = v4sigmatau3[dim->v4sigmatau3 * g + 8];
                auto v4sigmatau3_baab = v4sigmatau3[dim->v4sigmatau3 * g + 9];
                auto v4sigmatau3_babb = v4sigmatau3[dim->v4sigmatau3 * g + 10];

                // auto v4lapl4_aaaa = v4lapl4[dim->v4lapl4 * g + 0];
                // auto v4lapl4_aaab = v4lapl4[dim->v4lapl4 * g + 1];
                // auto v4lapl4_aabb = v4lapl4[dim->v4lapl4 * g + 2];
                // auto v4lapl4_abbb = v4lapl4[dim->v4lapl4 * g + 3];

                // auto v4lapl3tau_aaaa = v4lapl3tau[dim->v4lapl3tau * g + 0];
                // auto v4lapl3tau_aaab = v4lapl3tau[dim->v4lapl3tau * g + 1];
                // auto v4lapl3tau_aaba = v4lapl3tau[dim->v4lapl3tau * g + 2];
                // auto v4lapl3tau_aabb = v4lapl3tau[dim->v4lapl3tau * g + 3];
                // auto v4lapl3tau_abba = v4lapl3tau[dim->v4lapl3tau * g + 4];
                // auto v4lapl3tau_abbb = v4lapl3tau[dim->v4lapl3tau * g + 5];
                // auto v4lapl3tau_bbba = v4lapl3tau[dim->v4lapl3tau * g + 6];

                // auto v4lapl2tau2_aaaa = v4lapl2tau2[dim->v4lapl2tau2 * g + 0];
                // auto v4lapl2tau2_aaab = v4lapl2tau2[dim->v4lapl2tau2 * g + 1];
                // auto v4lapl2tau2_aabb = v4lapl2tau2[dim->v4lapl2tau2 * g + 2];
                // auto v4lapl2tau2_abaa = v4lapl2tau2[dim->v4lapl2tau2 * g + 3];
                // auto v4lapl2tau2_abab = v4lapl2tau2[dim->v4lapl2tau2 * g + 4];
                // auto v4lapl2tau2_abbb = v4lapl2tau2[dim->v4lapl2tau2 * g + 5];
                // auto v4lapl2tau2_bbaa = v4lapl2tau2[dim->v4lapl2tau2 * g + 6];
                // auto v4lapl2tau2_bbab = v4lapl2tau2[dim->v4lapl2tau2 * g + 7];

                // auto v4lapltau3_aaaa = v4lapltau3[dim->v4lapltau3 * g + 0];
                // auto v4lapltau3_aaab = v4lapltau3[dim->v4lapltau3 * g + 1];
                // auto v4lapltau3_aabb = v4lapltau3[dim->v4lapltau3 * g + 2];
                // auto v4lapltau3_abbb = v4lapltau3[dim->v4lapltau3 * g + 3];
                // auto v4lapltau3_baaa = v4lapltau3[dim->v4lapltau3 * g + 4];
                // auto v4lapltau3_baab = v4lapltau3[dim->v4lapltau3 * g + 5];
                // auto v4lapltau3_babb = v4lapltau3[dim->v4lapltau3 * g + 6];

                auto v4tau4_aaaa = v4tau4[dim->v4tau4 * g + 0];
                auto v4tau4_aaab = v4tau4[dim->v4tau4 * g + 1];
                auto v4tau4_aabb = v4tau4[dim->v4tau4 * g + 2];
                auto v4tau4_abbb = v4tau4[dim->v4tau4 * g + 3];

                // sums of functional derivatives

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

                // sigma and gamma
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

                // fourth-order

                double rrrr = v4rho4_abbb + 3.0 * v4rho4_aabb + 3.0 * v4rho4_aaab + v4rho4_aaaa;
                double rrrx = 2.0 * v4rho3sigma_abbc + 2.0 * v4rho3sigma_abbb + 2.0 * v4rho3sigma_abba + 4.0 * v4rho3sigma_aabc +
                              4.0 * v4rho3sigma_aabb + 4.0 * v4rho3sigma_aaba + 2.0 * v4rho3sigma_aaac + 2.0 * v4rho3sigma_aaab +
                              2.0 * v4rho3sigma_aaaa;
                // double rrrl  = v4rho3lapl_abbb + v4rho3lapl_abba + 2.0 * v4rho3lapl_aabb + 2.0 * v4rho3lapl_aaba + v4rho3lapl_aaab +
                // v4rho3lapl_aaaa;
                double rrrt = v4rho3tau_abbb + v4rho3tau_abba + 2.0 * v4rho3tau_aabb + 2.0 * v4rho3tau_aaba + v4rho3tau_aaab + v4rho3tau_aaaa;
                double rrxx = 4.0 * v4rho2sigma2_abcc + 8.0 * v4rho2sigma2_abcb + 4.0 * v4rho2sigma2_abbb + 8.0 * v4rho2sigma2_abac +
                              8.0 * v4rho2sigma2_abab + 4.0 * v4rho2sigma2_abaa + 4.0 * v4rho2sigma2_aacc + 8.0 * v4rho2sigma2_aacb +
                              4.0 * v4rho2sigma2_aabb + 8.0 * v4rho2sigma2_aaac + 8.0 * v4rho2sigma2_aaab + 4.0 * v4rho2sigma2_aaaa;
                // double rrxl  = 2.0 * v4rho2sigmalapl_abcb + 2.0 * v4rho2sigmalapl_abca + 2.0 * v4rho2sigmalapl_abbb + 2.0 * v4rho2sigmalapl_abba
                // + 2.0 * v4rho2sigmalapl_abab + 2.0 * v4rho2sigmalapl_abaa + 2.0 * v4rho2sigmalapl_aacb + 2.0 * v4rho2sigmalapl_aaca + 2.0 *
                // v4rho2sigmalapl_aabb + 2.0 * v4rho2sigmalapl_aaba + 2.0 * v4rho2sigmalapl_aaab + 2.0 * v4rho2sigmalapl_aaaa;
                double rrxt = 2.0 * v4rho2sigmatau_abcb + 2.0 * v4rho2sigmatau_abca + 2.0 * v4rho2sigmatau_abbb + 2.0 * v4rho2sigmatau_abba +
                              2.0 * v4rho2sigmatau_abab + 2.0 * v4rho2sigmatau_abaa + 2.0 * v4rho2sigmatau_aacb + 2.0 * v4rho2sigmatau_aaca +
                              2.0 * v4rho2sigmatau_aabb + 2.0 * v4rho2sigmatau_aaba + 2.0 * v4rho2sigmatau_aaab + 2.0 * v4rho2sigmatau_aaaa;
                // double rrll  = v4rho2lapl2_abbb + 2.0 * v4rho2lapl2_abab + v4rho2lapl2_abaa + v4rho2lapl2_aabb + 2.0 * v4rho2lapl2_aaab +
                // v4rho2lapl2_aaaa; double rrlt  = v4rho2lapltau_abbb + v4rho2lapltau_abba + v4rho2lapltau_abab + v4rho2lapltau_abaa +
                // v4rho2lapltau_aabb + v4rho2lapltau_aaba + v4rho2lapltau_aaab + v4rho2lapltau_aaaa;
                double rrtt = v4rho2tau2_abbb + 2.0 * v4rho2tau2_abab + v4rho2tau2_abaa + v4rho2tau2_aabb + 2.0 * v4rho2tau2_aaab + v4rho2tau2_aaaa;
                double rxxx = 8.0 * v4rhosigma3_accc + 24.0 * v4rhosigma3_accb + 24.0 * v4rhosigma3_acbb + 8.0 * v4rhosigma3_abbb +
                              24.0 * v4rhosigma3_aacc + 48.0 * v4rhosigma3_aacb + 24.0 * v4rhosigma3_aabb + 24.0 * v4rhosigma3_aaac +
                              24.0 * v4rhosigma3_aaab + 8.0 * v4rhosigma3_aaaa;
                // double rxxl  = 4.0 * v4rhosigma2lapl_accb + 4.0 * v4rhosigma2lapl_acca + 8.0 * v4rhosigma2lapl_acbb
                //              + 8.0 * v4rhosigma2lapl_acba + 4.0 * v4rhosigma2lapl_abbb + 4.0 * v4rhosigma2lapl_abba
                //              + 8.0 * v4rhosigma2lapl_aacb + 8.0 * v4rhosigma2lapl_aaca + 8.0 * v4rhosigma2lapl_aabb
                //              + 8.0 * v4rhosigma2lapl_aaba + 4.0 * v4rhosigma2lapl_aaab + 4.0 * v4rhosigma2lapl_aaaa;
                double rxxt = 4.0 * v4rhosigma2tau_accb + 4.0 * v4rhosigma2tau_acca + 8.0 * v4rhosigma2tau_acbb + 8.0 * v4rhosigma2tau_acba +
                              4.0 * v4rhosigma2tau_abbb + 4.0 * v4rhosigma2tau_abba + 8.0 * v4rhosigma2tau_aacb + 8.0 * v4rhosigma2tau_aaca +
                              8.0 * v4rhosigma2tau_aabb + 8.0 * v4rhosigma2tau_aaba + 4.0 * v4rhosigma2tau_aaab + 4.0 * v4rhosigma2tau_aaaa;
                // double rxll  = 2.0 * v4rhosigmalapl2_acbb + 4.0 * v4rhosigmalapl2_acab + 2.0 * v4rhosigmalapl2_acaa
                //              + 2.0 * v4rhosigmalapl2_abbb + 4.0 * v4rhosigmalapl2_abab + 2.0 * v4rhosigmalapl2_abaa
                //              + 2.0 * v4rhosigmalapl2_aabb + 4.0 * v4rhosigmalapl2_aaab + 2.0 * v4rhosigmalapl2_aaaa;
                // double rxlt  = 2.0 * v4rhosigmalapltau_acbb + 2.0 * v4rhosigmalapltau_acba + 2.0 * v4rhosigmalapltau_acab
                //              + 2.0 * v4rhosigmalapltau_acaa + 2.0 * v4rhosigmalapltau_abbb + 2.0 * v4rhosigmalapltau_abba
                //              + 2.0 * v4rhosigmalapltau_abab + 2.0 * v4rhosigmalapltau_abaa + 2.0 * v4rhosigmalapltau_aabb
                //              + 2.0 * v4rhosigmalapltau_aaba + 2.0 * v4rhosigmalapltau_aaab + 2.0 * v4rhosigmalapltau_aaaa;
                double rxtt = 2.0 * v4rhosigmatau2_acbb + 4.0 * v4rhosigmatau2_acab + 2.0 * v4rhosigmatau2_acaa + 2.0 * v4rhosigmatau2_abbb +
                              4.0 * v4rhosigmatau2_abab + 2.0 * v4rhosigmatau2_abaa + 2.0 * v4rhosigmatau2_aabb + 4.0 * v4rhosigmatau2_aaab +
                              2.0 * v4rhosigmatau2_aaaa;
                // double rlll  = v4rholapl3_abbb + 3.0 * v4rholapl3_aabb + 3.0 * v4rholapl3_aaab + v4rholapl3_aaaa;
                // double rllt  = v4rholapl2tau_abbb + v4rholapl2tau_abba
                //              + 2.0 * v4rholapl2tau_aabb + 2.0 * v4rholapl2tau_aaba
                //             + v4rholapl2tau_aaab + v4rholapl2tau_aaaa;
                // double rltt  = v4rholapltau2_abbb + 2.0 * v4rholapltau2_abab + v4rholapltau2_abaa
                //              + v4rholapltau2_aabb + 2.0 * v4rholapltau2_aaab + v4rholapltau2_aaaa;
                double rttt = v4rhotau3_abbb + 3.0 * v4rhotau3_aabb + 3.0 * v4rhotau3_aaab + v4rhotau3_aaaa;
                double xxxx = 8.0 * v4sigma4_cccc + 24.0 * v4sigma4_cccb + 24.0 * v4sigma4_ccbb + 8.0 * v4sigma4_cbbb + 40.0 * v4sigma4_accc +
                              96.0 * v4sigma4_accb + 72.0 * v4sigma4_acbb + 16.0 * v4sigma4_abbb + 72.0 * v4sigma4_aacc + 120.0 * v4sigma4_aacb +
                              48.0 * v4sigma4_aabb + 56.0 * v4sigma4_aaac + 48.0 * v4sigma4_aaab + 16.0 * v4sigma4_aaaa;
                double xxxr = 4.0 * v4rhosigma3_bccc + 8.0 * v4rhosigma3_bccb + 4.0 * v4rhosigma3_bcbb + 16.0 * v4rhosigma3_bacc +
                              24.0 * v4rhosigma3_bacb + 8.0 * v4rhosigma3_babb + 20.0 * v4rhosigma3_baac + 16.0 * v4rhosigma3_baab +
                              8.0 * v4rhosigma3_baaa + 4.0 * v4rhosigma3_accc + 8.0 * v4rhosigma3_accb + 4.0 * v4rhosigma3_acbb +
                              16.0 * v4rhosigma3_aacc + 24.0 * v4rhosigma3_aacb + 8.0 * v4rhosigma3_aabb + 20.0 * v4rhosigma3_aaac +
                              16.0 * v4rhosigma3_aaab + 8.0 * v4rhosigma3_aaaa;
                // double xxxl  = 4.0 * v4sigma3lapl_cccb + 4.0 * v4sigma3lapl_ccca + 8.0 * v4sigma3lapl_ccbb + 8.0 * v4sigma3lapl_ccba
                //              + 4.0 * v4sigma3lapl_cbbb + 4.0 * v4sigma3lapl_cbba + 16.0 * v4sigma3lapl_accb + 16.0 * v4sigma3lapl_acca
                //              + 24.0 * v4sigma3lapl_acbb + 24.0 * v4sigma3lapl_acba + 8.0 * v4sigma3lapl_abbb + 8.0 * v4sigma3lapl_abba
                //              + 20.0 * v4sigma3lapl_aacb + 20.0 * v4sigma3lapl_aaca + 16.0 * v4sigma3lapl_aabb + 16.0 * v4sigma3lapl_aaba
                //              + 8.0 * v4sigma3lapl_aaab + 8.0 * v4sigma3lapl_aaaa;
                double xxxt = 4.0 * v4sigma3tau_cccb + 4.0 * v4sigma3tau_ccca + 8.0 * v4sigma3tau_ccbb + 8.0 * v4sigma3tau_ccba +
                              4.0 * v4sigma3tau_cbbb + 4.0 * v4sigma3tau_cbba + 16.0 * v4sigma3tau_accb + 16.0 * v4sigma3tau_acca +
                              24.0 * v4sigma3tau_acbb + 24.0 * v4sigma3tau_acba + 8.0 * v4sigma3tau_abbb + 8.0 * v4sigma3tau_abba +
                              20.0 * v4sigma3tau_aacb + 20.0 * v4sigma3tau_aaca + 16.0 * v4sigma3tau_aabb + 16.0 * v4sigma3tau_aaba +
                              8.0 * v4sigma3tau_aaab + 8.0 * v4sigma3tau_aaaa;
                double xxrr = 2.0 * v4rho2sigma2_bbcc + 2.0 * v4rho2sigma2_bbcb + 6.0 * v4rho2sigma2_bbac + 4.0 * v4rho2sigma2_bbab +
                              4.0 * v4rho2sigma2_bbaa + 4.0 * v4rho2sigma2_abcc + 4.0 * v4rho2sigma2_abcb + 12.0 * v4rho2sigma2_abac +
                              8.0 * v4rho2sigma2_abab + 8.0 * v4rho2sigma2_abaa + 2.0 * v4rho2sigma2_aacc + 2.0 * v4rho2sigma2_aacb +
                              6.0 * v4rho2sigma2_aaac + 4.0 * v4rho2sigma2_aaab + 4.0 * v4rho2sigma2_aaaa;
                // double xxrl  = 2.0 * v4rhosigma2lapl_bccb + 2.0 * v4rhosigma2lapl_bcca + 2.0 * v4rhosigma2lapl_bcbb
                //              + 2.0 * v4rhosigma2lapl_bcba + 6.0 * v4rhosigma2lapl_bacb + 6.0 * v4rhosigma2lapl_baca
                //              + 4.0 * v4rhosigma2lapl_babb + 4.0 * v4rhosigma2lapl_baba + 4.0 * v4rhosigma2lapl_baab
                //              + 4.0 * v4rhosigma2lapl_baaa + 2.0 * v4rhosigma2lapl_accb + 2.0 * v4rhosigma2lapl_acca
                //              + 2.0 * v4rhosigma2lapl_acbb + 2.0 * v4rhosigma2lapl_acba + 6.0 * v4rhosigma2lapl_aacb
                //              + 6.0 * v4rhosigma2lapl_aaca + 4.0 * v4rhosigma2lapl_aabb + 4.0 * v4rhosigma2lapl_aaba
                //              + 4.0 * v4rhosigma2lapl_aaab + 4.0 * v4rhosigma2lapl_aaaa;
                double xxrt = 2.0 * v4rhosigma2tau_bccb + 2.0 * v4rhosigma2tau_bcca + 2.0 * v4rhosigma2tau_bcbb + 2.0 * v4rhosigma2tau_bcba +
                              6.0 * v4rhosigma2tau_bacb + 6.0 * v4rhosigma2tau_baca + 4.0 * v4rhosigma2tau_babb + 4.0 * v4rhosigma2tau_baba +
                              4.0 * v4rhosigma2tau_baab + 4.0 * v4rhosigma2tau_baaa + 2.0 * v4rhosigma2tau_accb + 2.0 * v4rhosigma2tau_acca +
                              2.0 * v4rhosigma2tau_acbb + 2.0 * v4rhosigma2tau_acba + 6.0 * v4rhosigma2tau_aacb + 6.0 * v4rhosigma2tau_aaca +
                              4.0 * v4rhosigma2tau_aabb + 4.0 * v4rhosigma2tau_aaba + 4.0 * v4rhosigma2tau_aaab + 4.0 * v4rhosigma2tau_aaaa;
                // double xxll  = 2.0 * v4sigma2lapl2_ccbb + 4.0 * v4sigma2lapl2_ccab + 2.0 * v4sigma2lapl2_ccaa + 2.0 * v4sigma2lapl2_cbbb + 4.0 *
                // v4sigma2lapl2_cbab + 2.0 * v4sigma2lapl2_cbaa + 6.0 * v4sigma2lapl2_acbb + 12.0 * v4sigma2lapl2_acab + 6.0 * v4sigma2lapl2_acaa
                // + 4.0 * v4sigma2lapl2_abbb + 8.0 * v4sigma2lapl2_abab + 4.0 * v4sigma2lapl2_abaa + 4.0 * v4sigma2lapl2_aabb + 8.0 *
                // v4sigma2lapl2_aaab + 4.0 * v4sigma2lapl2_aaaa; double xxlt  = 2.0 * v4sigma2lapltau_ccbb + 2.0 * v4sigma2lapltau_ccba + 2.0 *
                // v4sigma2lapltau_ccab + 2.0 * v4sigma2lapltau_ccaa + 2.0 * v4sigma2lapltau_cbbb + 2.0 * v4sigma2lapltau_cbba + 2.0 *
                // v4sigma2lapltau_cbab + 2.0 * v4sigma2lapltau_cbaa + 6.0 * v4sigma2lapltau_acbb + 6.0 * v4sigma2lapltau_acba + 6.0 *
                // v4sigma2lapltau_acab + 6.0 * v4sigma2lapltau_acaa + 4.0 * v4sigma2lapltau_abbb + 4.0 * v4sigma2lapltau_abba + 4.0 *
                // v4sigma2lapltau_abab + 4.0 * v4sigma2lapltau_abaa + 4.0 * v4sigma2lapltau_aabb + 4.0 * v4sigma2lapltau_aaba + 4.0 *
                // v4sigma2lapltau_aaab + 4.0 * v4sigma2lapltau_aaaa;
                double xxtt = 2.0 * v4sigma2tau2_ccbb + 4.0 * v4sigma2tau2_ccab + 2.0 * v4sigma2tau2_ccaa + 2.0 * v4sigma2tau2_cbbb +
                              4.0 * v4sigma2tau2_cbab + 2.0 * v4sigma2tau2_cbaa + 6.0 * v4sigma2tau2_acbb + 12.0 * v4sigma2tau2_acab +
                              6.0 * v4sigma2tau2_acaa + 4.0 * v4sigma2tau2_abbb + 8.0 * v4sigma2tau2_abab + 4.0 * v4sigma2tau2_abaa +
                              4.0 * v4sigma2tau2_aabb + 8.0 * v4sigma2tau2_aaab + 4.0 * v4sigma2tau2_aaaa;
                double xrrr = v4rho3sigma_bbbc + 2.0 * v4rho3sigma_bbba + 3.0 * v4rho3sigma_abbc + 6.0 * v4rho3sigma_abba + 3.0 * v4rho3sigma_aabc +
                              6.0 * v4rho3sigma_aaba + v4rho3sigma_aaac + 2.0 * v4rho3sigma_aaaa;
                // double xrrl  = v4rho2sigmalapl_bbcb + v4rho2sigmalapl_bbca + 2.0 * v4rho2sigmalapl_bbab + 2.0 * v4rho2sigmalapl_bbaa + 2.0 *
                // v4rho2sigmalapl_abcb + 2.0 * v4rho2sigmalapl_abca + 4.0 * v4rho2sigmalapl_abab + 4.0 * v4rho2sigmalapl_abaa + v4rho2sigmalapl_aacb
                // + v4rho2sigmalapl_aaca + 2.0 * v4rho2sigmalapl_aaab + 2.0 * v4rho2sigmalapl_aaaa;
                double xrrt = v4rho2sigmatau_bbcb + v4rho2sigmatau_bbca + 2.0 * v4rho2sigmatau_bbab + 2.0 * v4rho2sigmatau_bbaa +
                              2.0 * v4rho2sigmatau_abcb + 2.0 * v4rho2sigmatau_abca + 4.0 * v4rho2sigmatau_abab + 4.0 * v4rho2sigmatau_abaa +
                              v4rho2sigmatau_aacb + v4rho2sigmatau_aaca + 2.0 * v4rho2sigmatau_aaab + 2.0 * v4rho2sigmatau_aaaa;
                // double xrll  = v4rhosigmalapl2_bcbb + 2.0 * v4rhosigmalapl2_bcab + v4rhosigmalapl2_bcaa + 2.0 * v4rhosigmalapl2_babb + 4.0 *
                // v4rhosigmalapl2_baab + 2.0 * v4rhosigmalapl2_baaa + v4rhosigmalapl2_acbb + 2.0 * v4rhosigmalapl2_acab + v4rhosigmalapl2_acaa + 2.0
                // * v4rhosigmalapl2_aabb + 4.0 * v4rhosigmalapl2_aaab + 2.0 * v4rhosigmalapl2_aaaa; double xrlt  = v4rhosigmalapltau_bcbb +
                // v4rhosigmalapltau_bcba + v4rhosigmalapltau_bcab + v4rhosigmalapltau_bcaa + 2.0 * v4rhosigmalapltau_babb + 2.0 *
                // v4rhosigmalapltau_baba + 2.0 * v4rhosigmalapltau_baab + 2.0 * v4rhosigmalapltau_baaa + v4rhosigmalapltau_acbb +
                // v4rhosigmalapltau_acba + v4rhosigmalapltau_acab + v4rhosigmalapltau_acaa + 2.0 * v4rhosigmalapltau_aabb + 2.0 *
                // v4rhosigmalapltau_aaba + 2.0 * v4rhosigmalapltau_aaab + 2.0 * v4rhosigmalapltau_aaaa;
                double xrtt = v4rhosigmatau2_bcbb + 2.0 * v4rhosigmatau2_bcab + v4rhosigmatau2_bcaa + 2.0 * v4rhosigmatau2_babb +
                              4.0 * v4rhosigmatau2_baab + 2.0 * v4rhosigmatau2_baaa + v4rhosigmatau2_acbb + 2.0 * v4rhosigmatau2_acab +
                              v4rhosigmatau2_acaa + 2.0 * v4rhosigmatau2_aabb + 4.0 * v4rhosigmatau2_aaab + 2.0 * v4rhosigmatau2_aaaa;
                // double xlll  = v4sigmalapl3_cbbb + 3.0 * v4sigmalapl3_cabb + 3.0 * v4sigmalapl3_caab + v4sigmalapl3_caaa + 2.0 * v4sigmalapl3_abbb
                // + 6.0 * v4sigmalapl3_aabb + 6.0 * v4sigmalapl3_aaab + 2.0 * v4sigmalapl3_aaaa; double xllt  = v4sigmalapl2tau_cbbb +
                // v4sigmalapl2tau_cbba + 2.0 * v4sigmalapl2tau_cabb + 2.0 * v4sigmalapl2tau_caba
                //              + v4sigmalapl2tau_caab + v4sigmalapl2tau_caaa + 2.0 * v4sigmalapl2tau_abbb + 2.0 * v4sigmalapl2tau_abba
                //              + 4.0 * v4sigmalapl2tau_aabb + 4.0 * v4sigmalapl2tau_aaba + 2.0 * v4sigmalapl2tau_aaab + 2.0 * v4sigmalapl2tau_aaaa;
                // double xltt  = v4sigmalapltau2_cbbb + 2.0 * v4sigmalapltau2_cbab + v4sigmalapltau2_cbaa + v4sigmalapltau2_cabb + 2.0 *
                // v4sigmalapltau2_caab + v4sigmalapltau2_caaa + 2.0 * v4sigmalapltau2_abbb + 4.0 * v4sigmalapltau2_abab + 2.0 * v4sigmalapltau2_abaa
                // + 2.0 * v4sigmalapltau2_aabb + 4.0 * v4sigmalapltau2_aaab + 2.0 * v4sigmalapltau2_aaaa;
                double xttt = v4sigmatau3_cbbb + 3.0 * v4sigmatau3_cabb + 3.0 * v4sigmatau3_caab + v4sigmatau3_caaa + 2.0 * v4sigmatau3_abbb +
                              6.0 * v4sigmatau3_aabb + 6.0 * v4sigmatau3_aaab + 2.0 * v4sigmatau3_aaaa;
                // double llll  = v4lapl4_abbb + 3.0 * v4lapl4_aabb + 3.0 * v4lapl4_aaab + v4lapl4_aaaa;
                // double lllr  = v4rholapl3_babb + 2.0 * v4rholapl3_baab + v4rholapl3_baaa + v4rholapl3_aabb + 2.0 * v4rholapl3_aaab +
                // v4rholapl3_aaaa; double lllx  = 2.0 * v4sigmalapl3_cabb + 4.0 * v4sigmalapl3_caab + 2.0 * v4sigmalapl3_caaa + 2.0 *
                // v4sigmalapl3_babb + 4.0 * v4sigmalapl3_baab + 2.0 * v4sigmalapl3_baaa + 2.0 * v4sigmalapl3_aabb + 4.0 * v4sigmalapl3_aaab + 2.0 *
                // v4sigmalapl3_aaaa; double lllt  = v4lapl3tau_abbb + v4lapl3tau_abba + 2.0 * v4lapl3tau_aabb + 2.0 * v4lapl3tau_aaba +
                // v4lapl3tau_aaab + v4lapl3tau_aaaa; double llrr  = v4rho2lapl2_bbab + v4rho2lapl2_bbaa + 2.0 * v4rho2lapl2_abab + 2.0 *
                // v4rho2lapl2_abaa + v4rho2lapl2_aaab + v4rho2lapl2_aaaa; double llrx  = 2.0 * v4rhosigmalapl2_bcab + 2.0 * v4rhosigmalapl2_bcaa
                // + 2.0 * v4rhosigmalapl2_bbab + 2.0 * v4rhosigmalapl2_bbaa + 2.0 * v4rhosigmalapl2_baab + 2.0 * v4rhosigmalapl2_baaa + 2.0 *
                // v4rhosigmalapl2_acab + 2.0 * v4rhosigmalapl2_acaa + 2.0 * v4rhosigmalapl2_abab + 2.0 * v4rhosigmalapl2_abaa + 2.0 *
                // v4rhosigmalapl2_aaab + 2.0 * v4rhosigmalapl2_aaaa; double llrt  = v4rholapl2tau_babb + v4rholapl2tau_baba + v4rholapl2tau_baab +
                // v4rholapl2tau_baaa + v4rholapl2tau_aabb + v4rholapl2tau_aaba + v4rholapl2tau_aaab + v4rholapl2tau_aaaa; double llxx  = 4.0 *
                // v4sigma2lapl2_ccab + 4.0 * v4sigma2lapl2_ccaa + 8.0 * v4sigma2lapl2_cbab + 8.0 * v4sigma2lapl2_cbaa + 4.0 * v4sigma2lapl2_bbab
                // + 4.0 * v4sigma2lapl2_bbaa + 8.0 * v4sigma2lapl2_acab + 8.0 * v4sigma2lapl2_acaa + 8.0 * v4sigma2lapl2_abab + 8.0 *
                // v4sigma2lapl2_abaa + 4.0 * v4sigma2lapl2_aaab + 4.0 * v4sigma2lapl2_aaaa; double llxt  = 2.0 * v4sigmalapl2tau_cabb + 2.0 *
                // v4sigmalapl2tau_caba + 2.0 * v4sigmalapl2tau_caab + 2.0 * v4sigmalapl2tau_caaa + 2.0 * v4sigmalapl2tau_babb + 2.0 *
                // v4sigmalapl2tau_baba + 2.0 * v4sigmalapl2tau_baab + 2.0 * v4sigmalapl2tau_baaa + 2.0 * v4sigmalapl2tau_aabb + 2.0 *
                // v4sigmalapl2tau_aaba + 2.0 * v4sigmalapl2tau_aaab + 2.0 * v4sigmalapl2tau_aaaa; double lltt  = v4lapl2tau2_abbb + 2.0 *
                // v4lapl2tau2_abab + v4lapl2tau2_abaa + v4lapl2tau2_aabb + 2.0 * v4lapl2tau2_aaab + v4lapl2tau2_aaaa; double lrrr  = v4rho3lapl_bbba
                // + 3.0 * v4rho3lapl_abba + 3.0 * v4rho3lapl_aaba + v4rho3lapl_aaaa; double lrrx  = 2.0 * v4rho2sigmalapl_bbca + 2.0 *
                // v4rho2sigmalapl_bbba + 2.0 * v4rho2sigmalapl_bbaa + 4.0 * v4rho2sigmalapl_abca + 4.0 * v4rho2sigmalapl_abba + 4.0 *
                // v4rho2sigmalapl_abaa + 2.0 * v4rho2sigmalapl_aaca + 2.0 * v4rho2sigmalapl_aaba + 2.0 * v4rho2sigmalapl_aaaa; double lrrt  =
                // v4rho2lapltau_bbab + v4rho2lapltau_bbaa + 2.0 * v4rho2lapltau_abab + 2.0 * v4rho2lapltau_abaa + v4rho2lapltau_aaab +
                // v4rho2lapltau_aaaa; double lrxx  = 4.0 * v4rhosigma2lapl_bcca + 8.0 * v4rhosigma2lapl_bcba + 4.0 * v4rhosigma2lapl_bbba + 8.0 *
                // v4rhosigma2lapl_baca
                //              + 8.0 * v4rhosigma2lapl_baba + 4.0 * v4rhosigma2lapl_baaa + 4.0 * v4rhosigma2lapl_acca + 8.0 * v4rhosigma2lapl_acba
                //              + 4.0 * v4rhosigma2lapl_abba + 8.0 * v4rhosigma2lapl_aaca + 8.0 * v4rhosigma2lapl_aaba + 4.0 * v4rhosigma2lapl_aaaa;
                // double lrxt  = 2.0 * v4rhosigmalapltau_bcab + 2.0 * v4rhosigmalapltau_bcaa + 2.0 * v4rhosigmalapltau_bbab + 2.0 *
                // v4rhosigmalapltau_bbaa + 2.0 * v4rhosigmalapltau_baab + 2.0 * v4rhosigmalapltau_baaa + 2.0 * v4rhosigmalapltau_acab + 2.0 *
                // v4rhosigmalapltau_acaa + 2.0 * v4rhosigmalapltau_abab + 2.0 * v4rhosigmalapltau_abaa + 2.0 * v4rhosigmalapltau_aaab + 2.0 *
                // v4rhosigmalapltau_aaaa; double lrtt  = v4rholapltau2_babb + 2.0 * v4rholapltau2_baab + v4rholapltau2_baaa + v4rholapltau2_aabb
                // + 2.0 * v4rholapltau2_aaab + v4rholapltau2_aaaa; double lxxx  = 8.0 * v4sigma3lapl_ccca + 24.0 * v4sigma3lapl_ccba + 24.0 *
                // v4sigma3lapl_cbba + 8.0 * v4sigma3lapl_bbba + 24.0 * v4sigma3lapl_acca + 48.0 * v4sigma3lapl_acba + 24.0 * v4sigma3lapl_abba + 24.0
                // * v4sigma3lapl_aaca + 24.0 * v4sigma3lapl_aaba + 8.0 * v4sigma3lapl_aaaa; double lxxt  = 4.0 * v4sigma2lapltau_ccab + 4.0 *
                // v4sigma2lapltau_ccaa + 8.0 * v4sigma2lapltau_cbab + 8.0 * v4sigma2lapltau_cbaa + 4.0 * v4sigma2lapltau_bbab + 4.0 *
                // v4sigma2lapltau_bbaa + 8.0 * v4sigma2lapltau_acab + 8.0 * v4sigma2lapltau_acaa + 8.0 * v4sigma2lapltau_abab + 8.0 *
                // v4sigma2lapltau_abaa + 4.0 * v4sigma2lapltau_aaab + 4.0 * v4sigma2lapltau_aaaa; double lxtt  = 2.0 * v4sigmalapltau2_cabb + 4.0 *
                // v4sigmalapltau2_caab + 2.0 * v4sigmalapltau2_caaa + 2.0 * v4sigmalapltau2_babb
                //              + 4.0 * v4sigmalapltau2_baab + 2.0 * v4sigmalapltau2_baaa + 2.0 * v4sigmalapltau2_aabb + 4.0 * v4sigmalapltau2_aaab
                //              + 2.0 * v4sigmalapltau2_aaaa;
                // double lttt  = v4lapltau3_abbb + 3.0 * v4lapltau3_aabb + 3.0 * v4lapltau3_aaab + v4lapltau3_aaaa;
                double tttt = v4tau4_abbb + 3.0 * v4tau4_aabb + 3.0 * v4tau4_aaab + v4tau4_aaaa;
                double tttr = v4rhotau3_babb + 2.0 * v4rhotau3_baab + v4rhotau3_baaa + v4rhotau3_aabb + 2.0 * v4rhotau3_aaab + v4rhotau3_aaaa;
                double tttx = 2.0 * v4sigmatau3_cabb + 4.0 * v4sigmatau3_caab + 2.0 * v4sigmatau3_caaa + 2.0 * v4sigmatau3_babb +
                              4.0 * v4sigmatau3_baab + 2.0 * v4sigmatau3_baaa + 2.0 * v4sigmatau3_aabb + 4.0 * v4sigmatau3_aaab +
                              2.0 * v4sigmatau3_aaaa;
                // double tttl  = v4lapltau3_babb + 2.0 * v4lapltau3_baab + v4lapltau3_baaa + v4lapltau3_aabb + 2.0 * v4lapltau3_aaab +
                // v4lapltau3_aaaa;
                double ttrr = v4rho2tau2_bbab + v4rho2tau2_bbaa + 2.0 * v4rho2tau2_abab + 2.0 * v4rho2tau2_abaa + v4rho2tau2_aaab + v4rho2tau2_aaaa;
                double ttrx = 2.0 * v4rhosigmatau2_bcab + 2.0 * v4rhosigmatau2_bcaa + 2.0 * v4rhosigmatau2_bbab + 2.0 * v4rhosigmatau2_bbaa +
                              2.0 * v4rhosigmatau2_baab + 2.0 * v4rhosigmatau2_baaa + 2.0 * v4rhosigmatau2_acab + 2.0 * v4rhosigmatau2_acaa +
                              2.0 * v4rhosigmatau2_abab + 2.0 * v4rhosigmatau2_abaa + 2.0 * v4rhosigmatau2_aaab + 2.0 * v4rhosigmatau2_aaaa;
                // double ttrl  = v4rholapltau2_bbab + v4rholapltau2_bbaa + v4rholapltau2_baab + v4rholapltau2_baaa + v4rholapltau2_abab +
                // v4rholapltau2_abaa + v4rholapltau2_aaab + v4rholapltau2_aaaa;
                double ttxx = 4.0 * v4sigma2tau2_ccab + 4.0 * v4sigma2tau2_ccaa + 8.0 * v4sigma2tau2_cbab + 8.0 * v4sigma2tau2_cbaa +
                              4.0 * v4sigma2tau2_bbab + 4.0 * v4sigma2tau2_bbaa + 8.0 * v4sigma2tau2_acab + 8.0 * v4sigma2tau2_acaa +
                              8.0 * v4sigma2tau2_abab + 8.0 * v4sigma2tau2_abaa + 4.0 * v4sigma2tau2_aaab + 4.0 * v4sigma2tau2_aaaa;
                // double ttxl  = 2.0 * v4sigmalapltau2_cbab + 2.0 * v4sigmalapltau2_cbaa + 2.0 * v4sigmalapltau2_caab + 2.0 * v4sigmalapltau2_caaa
                // + 2.0 * v4sigmalapltau2_bbab + 2.0 * v4sigmalapltau2_bbaa + 2.0 * v4sigmalapltau2_baab + 2.0 * v4sigmalapltau2_baaa + 2.0 *
                // v4sigmalapltau2_abab + 2.0 * v4sigmalapltau2_abaa + 2.0 * v4sigmalapltau2_aaab + 2.0 * v4sigmalapltau2_aaaa; double ttll  =
                // v4lapl2tau2_bbab + v4lapl2tau2_bbaa + 2.0 * v4lapl2tau2_abab + 2.0 * v4lapl2tau2_abaa + v4lapl2tau2_aaab + v4lapl2tau2_aaaa;
                double trrr = v4rho3tau_bbba + 3.0 * v4rho3tau_abba + 3.0 * v4rho3tau_aaba + v4rho3tau_aaaa;
                double trrx = 2.0 * v4rho2sigmatau_bbca + 2.0 * v4rho2sigmatau_bbba + 2.0 * v4rho2sigmatau_bbaa + 4.0 * v4rho2sigmatau_abca +
                              4.0 * v4rho2sigmatau_abba + 4.0 * v4rho2sigmatau_abaa + 2.0 * v4rho2sigmatau_aaca + 2.0 * v4rho2sigmatau_aaba +
                              2.0 * v4rho2sigmatau_aaaa;
                // double trrl  = v4rho2lapltau_bbba + v4rho2lapltau_bbaa + 2.0 * v4rho2lapltau_abba + 2.0 * v4rho2lapltau_abaa + v4rho2lapltau_aaba +
                // v4rho2lapltau_aaaa;
                double trxx = 4.0 * v4rhosigma2tau_bcca + 8.0 * v4rhosigma2tau_bcba + 4.0 * v4rhosigma2tau_bbba + 8.0 * v4rhosigma2tau_baca +
                              8.0 * v4rhosigma2tau_baba + 4.0 * v4rhosigma2tau_baaa + 4.0 * v4rhosigma2tau_acca + 8.0 * v4rhosigma2tau_acba +
                              4.0 * v4rhosigma2tau_abba + 8.0 * v4rhosigma2tau_aaca + 8.0 * v4rhosigma2tau_aaba + 4.0 * v4rhosigma2tau_aaaa;
                // double trxl  = 2.0 * v4rhosigmalapltau_bcba + 2.0 * v4rhosigmalapltau_bcaa + 2.0 * v4rhosigmalapltau_bbba + 2.0 *
                // v4rhosigmalapltau_bbaa + 2.0 * v4rhosigmalapltau_baba + 2.0 * v4rhosigmalapltau_baaa + 2.0 * v4rhosigmalapltau_acba + 2.0 *
                // v4rhosigmalapltau_acaa + 2.0 * v4rhosigmalapltau_abba + 2.0 * v4rhosigmalapltau_abaa + 2.0 * v4rhosigmalapltau_aaba + 2.0 *
                // v4rhosigmalapltau_aaaa; double trll  = v4rholapl2tau_bbba + 2.0 * v4rholapl2tau_baba + v4rholapl2tau_baaa + v4rholapl2tau_abba
                // + 2.0 * v4rholapl2tau_aaba + v4rholapl2tau_aaaa;
                double txxx = 8.0 * v4sigma3tau_ccca + 24.0 * v4sigma3tau_ccba + 24.0 * v4sigma3tau_cbba + 8.0 * v4sigma3tau_bbba +
                              24.0 * v4sigma3tau_acca + 48.0 * v4sigma3tau_acba + 24.0 * v4sigma3tau_abba + 24.0 * v4sigma3tau_aaca +
                              24.0 * v4sigma3tau_aaba + 8.0 * v4sigma3tau_aaaa;
                // double txxl  = 4.0 * v4sigma2lapltau_ccba + 4.0 * v4sigma2lapltau_ccaa + 8.0 * v4sigma2lapltau_cbba + 8.0 * v4sigma2lapltau_cbaa
                // + 4.0 * v4sigma2lapltau_bbba + 4.0 * v4sigma2lapltau_bbaa + 8.0 * v4sigma2lapltau_acba + 8.0 * v4sigma2lapltau_acaa + 8.0 *
                // v4sigma2lapltau_abba + 8.0 * v4sigma2lapltau_abaa + 4.0 * v4sigma2lapltau_aaba + 4.0 * v4sigma2lapltau_aaaa; double txll  = 2.0 *
                // v4sigmalapl2tau_cbba + 4.0 * v4sigmalapl2tau_caba + 2.0 * v4sigmalapl2tau_caaa + 2.0 * v4sigmalapl2tau_bbba + 4.0 *
                // v4sigmalapl2tau_baba + 2.0 * v4sigmalapl2tau_baaa + 2.0 * v4sigmalapl2tau_abba + 4.0 * v4sigmalapl2tau_aaba + 2.0 *
                // v4sigmalapl2tau_aaaa; double tlll  = v4lapl3tau_bbba + 3.0 * v4lapl3tau_abba + 3.0 * v4lapl3tau_aaba + v4lapl3tau_aaaa;

                // Scalar contribution

                double tau_0 = 0.0;
                double rho_0 = 0.0;
                // double lap_0 = 0.0;

                // vxc 1 contributions

                rho_0 += rr * rhow[g] + rx * l2contract + rt * tauw[g];
                //+ rl * laplw[g];

                // lap_0 +=   lr * rhow[g]
                //          + lx * l2contract
                //          + lt * tauw[g]
                //          + ll * laplw[g];

                tau_0 += tr * rhow[g] + tx * l2contract + tt * tauw[g];
                //+ tl * laplw[g];

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

                // vxc 3 contributions

                rho_0 += rrrr * pi[g] +
                         rrrt * rrt_pi[g]
                         //+ rrrl * rrl_pi[g]
                         + rrtt * rtt_pi[g]
                         //+ rrlt * rtl_pi[g]
                         //+ rrll * rll_pi[g]
                         + rttt * ttt_pi[g]
                         //+ rltt * ttl_pi[g]
                         //+ rllt * tll_pi[g]
                         //+ rlll * lll_pi[g]

                         + rrrx * c2 +
                         rrxt * rt_c2
                         //+ rrxl * rl_c2
                         //+ rxll * ll_c2
                         + rxtt * tt_c2
                         //+ rxlt * tl_c2

                         + rrxx * c3
                         //+ rxxl * l_c3
                         + rxxt * t_c3

                         + rrx * c4
                         //+ rlx * l_c4
                         + rtx * t_c4

                         + rxx * (c5_6 + c8) + rxxx * c7;

                // lap_0 +=   lrrr * pi[g]
                //          + lrrt * rrt_pi[g]
                //          + llrr * rrl_pi[g]
                //          + lrtt * rtt_pi[g]
                //          + llrt * rtl_pi[g]
                //          + lllr * rll_pi[g]
                //          + lttt * ttt_pi[g]
                //          + lltt * ttl_pi[g]
                //          + lllt * tll_pi[g]
                //          + llll * lll_pi[g]
                //
                //          + lrrx * c2
                //          + lrxt * rt_c2
                //          + llrx * rl_c2
                //          + lllx * ll_c2
                //          + lxtt * tt_c2
                //          + llxt * tl_c2
                //
                //          + lrxx * c3
                //          + llxx * l_c3
                //          + lxxt * t_c3
                //
                //          + lrx * c4
                //          + llx * l_c4
                //          + ltx * t_c4
                //
                //          + lxx * (c5_6 + c8)
                //          + lxxx * c7;

                tau_0 += tttt * ttt_pi[g] +
                         tttr * rtt_pi[g]
                         //+ tttl * ttl_pi[g]
                         + ttrr * rrt_pi[g]
                         //+ ttrl * rtl_pi[g]
                         //+ ttll * tll_pi[g]
                         + trrr * pi[g]
                         //+ trrl * rrl_pi[g]
                         //+ trll * rll_pi[g]

                         + trrx * c2 +
                         ttrx * rt_c2
                         //+ trxl *rl_c2
                         //+ txll *ll_c2
                         + tttx * tt_c2
                         //+ ttxl *tl_c2

                         + trxx * c3
                         //+ txxl * l_c3
                         + ttxx * t_c3

                         + trx * c4
                         //+ tlx * l_c4
                         + ttx * t_c4

                         + txx * (c5_6 + c8) + txxx * c7;

                // vector contribution

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                // vxc 1 contribution

                xcomp += xr * grada_x_g * rhow[g] +
                         xt * grada_x_g * tauw[g]
                         //+ xl * grada_x_g * laplw[g]
                         + x * gradw_x[g] + xx * l5contract_x;

                ycomp += xr * grada_y_g * rhow[g] +
                         xt * grada_y_g * tauw[g]
                         //+ xl * grada_y_g * laplw[g]
                         + x * gradw_y[g] + xx * l5contract_y;

                zcomp += xr * grada_z_g * rhow[g] +
                         xt * grada_z_g * tauw[g]
                         //+ xl * grada_z_g * laplw[g]
                         + x * gradw_z[g] + xx * l5contract_z;

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

                // vxc 3 contributions

                xcomp += xrrr * grada_x_g * pi[g]  // c9 terms
                         + xrrt * grada_x_g * rrt_pi[g]
                         //+ xrrl * grada_x_g * rrl_pi[g]
                         + xrtt * grada_x_g * rtt_pi[g]
                         //+ xrlt * grada_x_g * rtl_pi[g]
                         //+ xrll * grada_x_g * rll_pi[g]
                         + xttt * grada_x_g * ttt_pi[g]
                         //+ xltt * grada_x_g * ttl_pi[g]
                         //+ xllt * grada_x_g * tll_pi[g]
                         //+ xlll * grada_x_g * lll_pi[g]

                         + xrr * pix[g]  // c10 terms
                         + xrt * rt_pix[g]
                         //+ xrl * rl_pix[g]
                         //+ xll * ll_pix[g]
                         + xtt * tt_pix[g]
                         //+ xtl * tl_pix[g]

                         + xxrr * c2 * grada_x_g  // c11 terms
                         + xxrt * rt_c2 * grada_x_g
                         //+ xxrl * rl_c2 * grada_x_g
                         //+ xxll * ll_c2 * grada_x_g
                         + xxtt * tt_c2 * grada_x_g
                         //+ xxlt * tl_c2 * grada_x_g

                         + xxr * (c12_c14_x + c15_x)
                         //+ xxl * (l_c12_c14_x + l_c15_x)
                         + xxt * (t_c12_c14_x + t_c15_x)

                         + xxxr * c13_x
                         //+ xxxl * l_c13_x
                         + xxxt * t_c13_x

                         + xx * c17_24_25_x + xxxx * c18_x + xxx * (c16_19_22_x + c20_21_23_x);

                ycomp += xrrr * grada_y_g * pi[g]  // c9 terms
                         + xrrt * grada_y_g * rrt_pi[g]
                         //+ xrrl * grada_y_g * rrl_pi[g]
                         + xrtt * grada_y_g * rtt_pi[g]
                         //+ xrlt * grada_y_g * rtl_pi[g]
                         //+ xrll * grada_y_g * rll_pi[g]
                         + xttt * grada_y_g * ttt_pi[g]
                         //+ xltt * grada_y_g * ttl_pi[g]
                         //+ xllt * grada_y_g * tll_pi[g]
                         //+ xlll * grada_y_g * lll_pi[g]

                         + xrr * piy[g]  // c10 terms
                         + xrt * rt_piy[g]
                         //+ xrl * rl_piy[g]
                         //+ xll * ll_piy[g]
                         + xtt * tt_piy[g]
                         //+ xtl * tl_piy[g]

                         + xxrr * c2 * grada_y_g  // c11 terms
                         + xxrt * rt_c2 * grada_y_g
                         //+ xxrl * rl_c2 * grada_y_g
                         //+ xxll * ll_c2 * grada_y_g
                         + xxtt * tt_c2 * grada_y_g
                         //+ xxlt * tl_c2 * grada_y_g

                         + xxr * (c12_c14_y + c15_y)
                         //+ xxl * (l_c12_c14_y + l_c15_y)
                         + xxt * (t_c12_c14_y + t_c15_y)

                         + xxxr * c13_y
                         //+ xxxl * l_c13_y
                         + xxxt * t_c13_y

                         + xx * c17_24_25_y + xxxx * c18_y + xxx * (c16_19_22_y + c20_21_23_y);

                zcomp += xrrr * grada_z_g * pi[g]  // c9 terms
                         + xrrt * grada_z_g * rrt_pi[g]
                         //+ xrrl * grada_z_g * rrl_pi[g]
                         + xrtt * grada_z_g * rtt_pi[g]
                         //+ xrlt * grada_z_g * rtl_pi[g]
                         //+ xrll * grada_z_g * rll_pi[g]
                         + xttt * grada_z_g * ttt_pi[g]
                         //+ xltt * grada_z_g * ttl_pi[g]
                         //+ xllt * grada_z_g * tll_pi[g]
                         //+ xlll * grada_z_g * lll_pi[g]

                         + xrr * piz[g]  // c10 terms
                         + xrt * rt_piz[g]
                         //+ xrl * rl_piz[g]
                         //+ xll * ll_piz[g]
                         + xtt * tt_piz[g]
                         //+ xtl * tl_piz[g]

                         + xxrr * c2 * grada_z_g  // c11 terms
                         + xxrt * rt_c2 * grada_z_g
                         //+ xxrl * rl_c2 * grada_z_g
                         //+ xxll * ll_c2 * grada_z_g
                         + xxtt * tt_c2 * grada_z_g
                         //+ xxlt * tl_c2 * grada_z_g

                         + xxr * (c12_c14_z + c15_z)
                         //+ xxl * (l_c12_c14_z + l_c15_z)
                         + xxt * (t_c12_c14_z + t_c15_z)

                         + xxxr * c13_z
                         //+ xxxl * l_c13_z
                         + xxxt * t_c13_z

                         + xx * c17_24_25_z + xxxx * c18_z + xxx * (c16_19_22_z + c20_21_23_z);

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

    timer.stop("Lxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Lxc matrix matmul");

    // LDA and GGA contribution
    auto mat_Lxc = denblas::multABt(gtoValues, mat_G);

    auto mat_Lxc_gga = denblas::multABt(gtoValues, mat_G_gga);

    mat_Lxc_gga.symmetrize();  // (matrix + matrix.T)

    mat_Lxc = denblas::addAB(mat_Lxc, mat_Lxc_gga, 1.0);

    // tau contribution
    auto mat_Lxc_x = denblas::multABt(gtoValuesX, mat_G_gga_x);
    auto mat_Lxc_y = denblas::multABt(gtoValuesY, mat_G_gga_y);
    auto mat_Lxc_z = denblas::multABt(gtoValuesZ, mat_G_gga_z);

    mat_Lxc = denblas::addAB(mat_Lxc, mat_Lxc_x, 0.5);
    mat_Lxc = denblas::addAB(mat_Lxc, mat_Lxc_y, 0.5);
    mat_Lxc = denblas::addAB(mat_Lxc, mat_Lxc_z, 0.5);

    timer.stop("Lxc matrix matmul");

    return mat_Lxc;
}

}  // namespace xcintmgga
