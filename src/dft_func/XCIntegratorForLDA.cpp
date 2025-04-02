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

#include "XCIntegratorForLDA.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <sstream>

#include "DenseMatrix.hpp"
#include "DftSubMatrix.hpp"
#include "ErrorHandler.hpp"
#include "GtoFunc.hpp"
#include "GtoValues.hpp"
#include "MathFunc.hpp"
#include "MultiTimer.hpp"
#include "Prescreener.hpp"
#include "SerialDenseLinearAlgebra.hpp"
#include "SerialDensityGridGenerator.hpp"
#include "StringFormat.hpp"

namespace xcintlda {  // xcintlda namespace

auto
newIntegrateVxcFockForLdaClosedShell(const CMolecule&                  molecule,
                                     const CMolecularBasis&            basis,
                                     const std::vector<const double*>& gsDensityPointers,
                                     const CMolecularGrid&             molecularGrid,
                                     const double                      screeningThresholdForGTOValues,
                                     const CXCFunctional&              xcFunctional) -> CAOKohnShamMatrix
{
    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // Kohn-Sham matrix

    CAOKohnShamMatrix mat_Vxc(naos, naos, std::string("closedshell"));

    mat_Vxc.zero();

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

    const auto nboxes = counts.size();

    const auto ngblocks = gto_blocks.size();

    // set up pointers to OMP data

    auto ptr_counts = counts.data();

    auto ptr_displacements = displacements.data();

    auto ptr_gto_blocks = gto_blocks.data();

    auto ptr_gsDensityPointers = gsDensityPointers.data();

    auto ptr_xcFunctional = &xcFunctional;

#pragma omp parallel shared(ptr_counts, ptr_displacements, xcoords, ycoords, zcoords, ptr_gto_blocks, ptr_gsDensityPointers, ptr_xcFunctional, nboxes, ngblocks, naos, nele, xcene, mat_Vxc)
    {
#pragma omp single nowait
        {
            for (size_t box_id = 0; box_id < nboxes; box_id++)
            {
                #pragma omp task firstprivate(box_id)
                {
                    // grid points in box

                    auto npoints = ptr_counts[box_id];

                    auto gridblockpos = ptr_displacements[box_id];

                    // dimension of grid box

                    auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

                    // prescreening

                    std::vector<std::vector<int>> cgto_mask_blocks, pre_ao_inds_blocks;

                    std::vector<int> aoinds;

                    // reserve prescreening vectors

                    cgto_mask_blocks.reserve(ngblocks);

                    pre_ao_inds_blocks.reserve(ngblocks);

                    aoinds.reserve(naos); 

                    for (size_t i = 0; i < ngblocks; i++)
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

                    // compute VXC contributions

                    if (aocount > 0)
                    {
                        // density matrix slicing

                        auto sub_dens_mat = dftsubmat::getSubDensityMatrix(ptr_gsDensityPointers[0], aoinds, naos);

                        // GTO values on grid points

                        CDenseMatrix mat_chi(aocount, npoints);

                        const auto grid_x_ptr = xcoords + gridblockpos;

                        const auto grid_y_ptr = ycoords + gridblockpos;

                        const auto grid_z_ptr = zcoords + gridblockpos;

                        std::vector<double> grid_x(grid_x_ptr, grid_x_ptr + npoints);

                        std::vector<double> grid_y(grid_y_ptr, grid_y_ptr + npoints);

                        std::vector<double> grid_z(grid_z_ptr, grid_z_ptr + npoints);

                        // go through GTO blocks

                        for (size_t i_block = 0, idx = 0; i_block < ngblocks; i_block++)
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

                        // set up XC functional data

                        auto ldafunc = ptr_xcFunctional->getFunctionalPointerToLdaComponent();

                        const auto dim = &(ldafunc->dim);

                        auto rho = std::vector<double>(dim->rho * npoints);

                        auto exc = std::vector<double>(dim->zk * npoints);

                        auto vrho = std::vector<double>(dim->vrho * npoints);

                        // compute density on grid

                        sdengridgen::serialGenerateDensityForLDA(rho.data(), mat_chi, sub_dens_mat);

                        // compute XC functional values

                        auto loc_xcfunc = CXCFunctional(*ptr_xcFunctional);

                        loc_xcfunc.compute_exc_vxc_for_lda(npoints, rho.data(), exc.data(), vrho.data());

                        // copy weights to local vector

                        std::vector<double> local_weights(weights + gridblockpos, weights + gridblockpos + npoints);

                        // compute Vxc matrix

                        CDenseMatrix mat_G(aocount, npoints);

                        auto G_val = mat_G.values();

                        auto chi_val = mat_chi.values();

                        for (int nu = 0; nu < aocount; nu++)
                        {
                            auto nu_offset = nu * npoints;

                            #pragma omp simd
                            for (int g = 0; g < npoints; g++)
                            {
                                G_val[nu_offset + g] = local_weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];
                            }
                        }

                        auto local_mat_Vxc = sdenblas::serialMultABt(mat_chi, mat_G);

                        // compute local number of electrons and XC energy

                        double local_nele = 0.0, local_xcene = 0.0;

                        for (int g = 0; g < npoints; g++)
                        {
                            auto rho_total = rho[2 * g + 0] + rho[2 * g + 1];

                            local_nele += local_weights[g] * rho_total;

                            local_xcene += local_weights[g] * exc[g] * rho_total;
                        }

                        // update Vxc matrix, number of electrons, and XC energy

                        #pragma omp critical
                        {
                            nele += local_nele;

                            xcene += local_xcene;

                            dftsubmat::distributeSubMatrixToKohnSham(mat_Vxc, local_mat_Vxc, aoinds);
                        }
                    }
                }
            }
        }
    }

    mat_Vxc.setNumberOfElectrons(nele);

    mat_Vxc.setExchangeCorrelationEnergy(xcene);

    return mat_Vxc;
}

auto
integrateVxcFockForLdaClosedShell(const CMolecule&                  molecule,
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

    auto       ldafunc = xcFunctional.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<CXCFunctional> omp_xcfuncs(nthreads, CXCFunctional(xcFunctional));

    std::vector<std::vector<double>> omp_local_weights_data(nthreads, std::vector<double>(omp_max_npoints));

    std::vector<std::vector<double>> omp_rho_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));

    std::vector<std::vector<double>> omp_exc_data(nthreads, std::vector<double>(dim->zk * omp_max_npoints));
    std::vector<std::vector<double>> omp_vrho_data(nthreads, std::vector<double>(dim->vrho * omp_max_npoints));

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
                    std::memcpy(mat_chi.row(idx), submat_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");

            omptimers[thread_id].start("Generate density grid");

            auto local_weights = omp_local_weights_data[thread_id].data();

            auto rho     = omp_rho_data[thread_id].data();

            auto exc    = omp_exc_data[thread_id].data();
            auto vrho   = omp_vrho_data[thread_id].data();

            sdengridgen::serialGenerateDensityForLDA(rho, mat_chi, sub_dens_mat);

            omptimers[thread_id].stop("Generate density grid");

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_exc_vxc_for_lda(grid_batch_size, rho, exc, vrho);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            omptimers[thread_id].start("Vxc matrix G");

            CDenseMatrix mat_G(aocount, grid_batch_size);

            auto G_val = mat_G.values();

            auto chi_val = mat_chi.values();

            for (int nu = 0; nu < aocount; nu++)
            {
                auto nu_offset = nu * grid_batch_size;

#pragma omp simd
                for (int g = 0; g < grid_batch_size; g++)
                {
                    G_val[nu_offset + g] = local_weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];
                }
            }

            omptimers[thread_id].stop("Vxc matrix G");

            omptimers[thread_id].start("Vxc matmul");

            auto partial_mat_Vxc = sdenblas::serialMultABt(mat_chi, mat_G);

            omptimers[thread_id].stop("Vxc matmul");

            omptimers[thread_id].start("Vxc local matrix dist.");

#pragma omp critical
            sdenblas::serialInPlaceAddAB(sum_partial_mat_Vxc, partial_mat_Vxc);

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
integrateVxcFockForLdaOpenShell(const CMolecule&                  molecule,
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

    auto       ldafunc = xcFunctional.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<CXCFunctional> omp_xcfuncs(nthreads, CXCFunctional(xcFunctional));

    std::vector<std::vector<double>> omp_local_weights_data(nthreads, std::vector<double>(omp_max_npoints));

    std::vector<std::vector<double>> omp_rho_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));

    std::vector<std::vector<double>> omp_exc_data(nthreads, std::vector<double>(dim->zk * omp_max_npoints));
    std::vector<std::vector<double>> omp_vrho_data(nthreads, std::vector<double>(dim->vrho * omp_max_npoints));

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
                    std::memcpy(mat_chi.row(idx), submat_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");

            omptimers[thread_id].start("Generate density grid");

            auto local_weights = omp_local_weights_data[thread_id].data();

            auto rho     = omp_rho_data[thread_id].data();

            auto exc    = omp_exc_data[thread_id].data();
            auto vrho   = omp_vrho_data[thread_id].data();

            sdengridgen::serialGenerateDensityForLDA(rho, mat_chi, sub_dens_mat_a, sub_dens_mat_b);

            omptimers[thread_id].stop("Generate density grid");

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_exc_vxc_for_lda(grid_batch_size, rho, exc, vrho);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            omptimers[thread_id].start("Vxc matrix G");

            CDenseMatrix mat_G_a(aocount, grid_batch_size);

            CDenseMatrix mat_G_b(aocount, grid_batch_size);

            auto G_a_val = mat_G_a.values();

            auto G_b_val = mat_G_b.values();

            auto chi_val = mat_chi.values();

            for (int nu = 0; nu < aocount; nu++)
            {
                auto nu_offset = nu * grid_batch_size;

#pragma omp simd
                for (int g = 0; g < grid_batch_size; g++)
                {
                    G_a_val[nu_offset + g] = local_weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];

                    G_b_val[nu_offset + g] = local_weights[g] * vrho[2 * g + 1] * chi_val[nu_offset + g];
                }
            }

            omptimers[thread_id].stop("Vxc matrix G");

            omptimers[thread_id].start("Vxc matmul");

            auto partial_mat_Vxc_a = sdenblas::serialMultABt(mat_chi, mat_G_a);
            auto partial_mat_Vxc_b = sdenblas::serialMultABt(mat_chi, mat_G_b);

            omptimers[thread_id].stop("Vxc matmul");

            omptimers[thread_id].start("Vxc local matrix dist.");

#pragma omp critical
            {
                sdenblas::serialInPlaceAddAB(sum_partial_mat_Vxc_a, partial_mat_Vxc_a);
                sdenblas::serialInPlaceAddAB(sum_partial_mat_Vxc_b, partial_mat_Vxc_b);
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
integrateFxcFockForLdaClosedShell(const std::vector<double*>&       aoFockPointers,
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

    auto omp_max_npoints = max_npoints_per_box / nthreads;
    if (max_npoints_per_box % nthreads != 0) omp_max_npoints++;

    // density and functional derivatives

    auto       ldafunc = xcFunctional.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<CXCFunctional> omp_xcfuncs(nthreads, CXCFunctional(xcFunctional));

    std::vector<std::vector<double>> omp_local_weights_data(nthreads, std::vector<double>(omp_max_npoints));

    std::vector<std::vector<double>> omp_rho_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));
    std::vector<std::vector<double>> omp_rhow_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));

    std::vector<std::vector<double>> omp_v2rho2_data(nthreads, std::vector<double>(dim->v2rho2 * omp_max_npoints));

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

        // generate sub density matrix and density grid

        timer.start("Density matrix slicing");

        auto sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

        std::vector<CDenseMatrix> rw_sub_dens_mat_vec(rwDensityPointers.size());

        for (size_t idensity = 0; idensity < rwDensityPointers.size(); idensity++)
        {
            rw_sub_dens_mat_vec[idensity] = dftsubmat::getSubDensityMatrix(rwDensityPointers[idensity], aoinds, naos);
        }

        timer.stop("Density matrix slicing");

        // GTO values on grid points

        timer.start("OMP Fxc calc.");

        std::vector<CDenseMatrix> sum_partial_mat_Fxc(rwDensityPointers.size(), CDenseMatrix(aocount, aocount));

#pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            omptimers[thread_id].start("gtoeval");

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            CDenseMatrix mat_chi(aocount, grid_batch_size);

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
                    std::memcpy(mat_chi.row(idx), submat_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");

            omptimers[thread_id].start("Generate density grid");

            auto local_weights = omp_local_weights_data[thread_id].data();

            auto rho     = omp_rho_data[thread_id].data();
            auto rhow     = omp_rhow_data[thread_id].data();

            auto v2rho2     = omp_v2rho2_data[thread_id].data();

            sdengridgen::serialGenerateDensityForLDA(rho, mat_chi, sub_dens_mat);

            omptimers[thread_id].stop("Generate density grid");

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_fxc_for_lda(grid_batch_size, rho, v2rho2);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            // go through rhow density matrices

            for (size_t idensity = 0; idensity < rwDensityPointers.size(); idensity++)
            {
                omptimers[thread_id].start("Generate density grid");

                sdengridgen::serialGenerateDensityForLDA(rhow, mat_chi, rw_sub_dens_mat_vec[idensity]);

                omptimers[thread_id].stop("Generate density grid");

                omptimers[thread_id].start("Fxc matrix G");

                CDenseMatrix mat_G(aocount, grid_batch_size);

                auto G_val = mat_G.values();

                auto chi_val = mat_chi.values();

                auto       ldafunc = omp_xcfuncs[thread_id].getFunctionalPointerToLdaComponent();
                const auto dim     = &(ldafunc->dim);

                for (int nu = 0; nu < aocount; nu++)
                {
                    auto nu_offset = nu * grid_batch_size;

#pragma omp simd
                    for (int g = 0; g < grid_batch_size; g++)
                    {
                        auto rhow_a = rhow[2 * g + 0];
                        auto rhow_b = rhow[2 * g + 1];

                        // functional derivatives

                        // first-order

                        // second-order

                        auto v2rho2_aa = v2rho2[dim->v2rho2 * g + 0];
                        auto v2rho2_ab = v2rho2[dim->v2rho2 * g + 1];

                        G_val[nu_offset + g] = local_weights[g] * (v2rho2_aa * rhow_a + v2rho2_ab * rhow_b) * chi_val[nu_offset + g];
                    }
                }

                omptimers[thread_id].stop("Fxc matrix G");

                omptimers[thread_id].start("Fxc matrix matmul");

                auto partial_mat_Fxc = sdenblas::serialMultABt(mat_chi, mat_G);

                omptimers[thread_id].stop("Fxc matrix matmul");

                omptimers[thread_id].start("Fxc local matrix dist.");

#pragma omp critical
                sdenblas::serialInPlaceAddAB(sum_partial_mat_Fxc[idensity], partial_mat_Fxc);

                omptimers[thread_id].stop("Fxc local matrix dist.");
            }
        }

        timer.stop("OMP Fxc calc.");

        timer.start("Fxc matrix dist.");

        for (size_t idensity = 0; idensity < rwDensityPointers.size(); idensity++)
        {
            dftsubmat::distributeSubMatrixToFock(aoFockPointers, idensity, sum_partial_mat_Fxc[idensity], aoinds, naos);
        }

        timer.stop("Fxc matrix dist.");
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
integrateKxcFockForLdaClosedShell(const std::vector<double*>& aoFockPointers,
                                  const CMolecule&        molecule,
                                  const CMolecularBasis&  basis,
                                  const std::vector<const double*>& rwDensityPointers,
                                  const std::vector<const double*>& rw2DensityPointers,
                                  const std::vector<const double*>& gsDensityPointers,
                                  const CMolecularGrid&   molecularGrid,
                                  const double            screeningThresholdForGTOValues,
                                  const CXCFunctional&    xcFunctional,
                                  const std::string&      quadMode) -> void
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

    auto omp_max_npoints = max_npoints_per_box / nthreads;
    if (max_npoints_per_box % nthreads != 0) omp_max_npoints++;

    // density and functional derivatives

    auto       ldafunc = xcFunctional.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<CXCFunctional> omp_xcfuncs(nthreads, CXCFunctional(xcFunctional));

    std::vector<std::vector<double>> omp_local_weights_data(nthreads, std::vector<double>(omp_max_npoints));

    std::vector<std::vector<double>> omp_rho_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));

    std::vector<std::vector<double>> omp_v2rho2_data(nthreads, std::vector<double>(dim->v2rho2 * omp_max_npoints));
    std::vector<std::vector<double>> omp_v3rho3_data(nthreads, std::vector<double>(dim->v3rho3 * omp_max_npoints));

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

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

        auto rw_sub_dens_mat = dftsubmat::getSubAODensityMatrix(rwDensityPointers, aoinds, naos);

        auto rw2_sub_dens_mat = dftsubmat::getSubAODensityMatrix(rw2DensityPointers, aoinds, naos);

        timer.stop("Density matrix slicing");

        // GTO values on grid points

        timer.start("OMP Kxc calc.");

        std::vector<CDenseMatrix> sum_partial_mat_Kxc(rw2DensityPointers.size(), CDenseMatrix(aocount, aocount));

#pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            omptimers[thread_id].start("gtoeval");

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            CDenseMatrix mat_chi(aocount, grid_batch_size);

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
                    std::memcpy(mat_chi.row(idx), submat_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");

            omptimers[thread_id].start("Generate density grid");

            auto local_weights = omp_local_weights_data[thread_id].data();

            auto rho        = omp_rho_data[thread_id].data();

            auto v2rho2     = omp_v2rho2_data[thread_id].data();
            auto v3rho3     = omp_v3rho3_data[thread_id].data();

            sdengridgen::serialGenerateDensityForLDA(rho, mat_chi, gs_sub_dens_mat);

            auto xcfuntype = omp_xcfuncs[thread_id].getFunctionalType();

            auto rwdengrid = sdengridgen::serialGenerateDensityGridForLDA(mat_chi, rw_sub_dens_mat, xcfuntype);

            auto rw2dengrid = sdengridgen::serialGenerateDensityGridForLDA(mat_chi, rw2_sub_dens_mat, xcfuntype);

            omptimers[thread_id].stop("Generate density grid");

            // compute perturbed density

            omptimers[thread_id].start("Density grid quad");

            auto numdens_rw2 = static_cast<int>(rw2DensityPointers.size());

            CDensityGridQuad rwdengridquad(grid_batch_size, numdens_rw2, xcfuntype, dengrid::ab);

            rwdengridquad.DensityProd(rwdengrid, xcfuntype, numdens_rw2, quadMode);

            omptimers[thread_id].stop("Density grid quad");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_fxc_for_lda(grid_batch_size, rho, v2rho2);

            omp_xcfuncs[thread_id].compute_kxc_for_lda(grid_batch_size, rho, v3rho3);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            // go through density matrices

            for (int idensity = 0; idensity < numdens_rw2; idensity++)
            {
                // compute partial contribution to Kxc matrix

                auto rhow1a = rwdengridquad.gam(idensity);

                auto rhow12a = rw2dengrid.alphaDensity(idensity);

                auto rhow12b = rw2dengrid.betaDensity(idensity);

                std::vector<const double*> rwdengrid_pointers({rhow1a});
                std::vector<const double*> rw2dengrid_pointers({rhow12a, rhow12b});

                auto partial_mat_Kxc = integratePartialKxcFockForLdaClosedShell(
                    omp_xcfuncs[thread_id], local_weights, mat_chi, v2rho2, v3rho3,
                    rwdengrid_pointers, rw2dengrid_pointers, omptimers[thread_id]);

                // accumulate partial Kxc

                omptimers[thread_id].start("Kxc local matrix dist.");

                #pragma omp critical
                sdenblas::serialInPlaceAddAB(sum_partial_mat_Kxc[idensity], partial_mat_Kxc);

                omptimers[thread_id].stop("Kxc local matrix dist.");
            }
        }

        timer.stop("OMP Kxc calc.");

        timer.start("Kxc matrix dist.");

        for (size_t idensity = 0; idensity < rw2DensityPointers.size(); idensity++)
        {
            dftsubmat::distributeSubMatrixToFock(aoFockPointers, idensity, sum_partial_mat_Kxc[idensity], aoinds, naos);
        }

        timer.stop("Kxc matrix dist.");
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
integrateKxcLxcFockForLdaClosedShell(const std::vector<double*>& aoFockPointers,
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

    auto omp_max_npoints = max_npoints_per_box / nthreads;
    if (max_npoints_per_box % nthreads != 0) omp_max_npoints++;

    // density and functional derivatives

    auto       ldafunc = xcFunctional.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<CXCFunctional> omp_xcfuncs(nthreads, CXCFunctional(xcFunctional));

    std::vector<std::vector<double>> omp_local_weights_data(nthreads, std::vector<double>(omp_max_npoints));

    std::vector<std::vector<double>> omp_rho_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));

    std::vector<std::vector<double>> omp_v2rho2_data(nthreads, std::vector<double>(dim->v2rho2 * omp_max_npoints));
    std::vector<std::vector<double>> omp_v3rho3_data(nthreads, std::vector<double>(dim->v3rho3 * omp_max_npoints));
    std::vector<std::vector<double>> omp_v4rho4_data(nthreads, std::vector<double>(dim->v4rho4 * omp_max_npoints));

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

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

        auto rw_sub_dens_mat = dftsubmat::getSubAODensityMatrix(rwDensityPointers, aoinds, naos);

        auto rw2_sub_dens_mat = dftsubmat::getSubAODensityMatrix(rw2DensityPointers, aoinds, naos);

        auto rw3_sub_dens_mat = dftsubmat::getSubAODensityMatrix(rw3DensityPointers, aoinds, naos);

        timer.stop("Density matrix slicing");

        // GTO values on grid points

        timer.start("OMP KxcLxc calc.");

        std::vector<CDenseMatrix> sum_partial_mat_Kxc(rw2DensityPointers.size(), CDenseMatrix(aocount, aocount));

        std::vector<CDenseMatrix> sum_partial_mat_Lxc(rw3DensityPointers.size(), CDenseMatrix(aocount, aocount));

#pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            omptimers[thread_id].start("gtoeval");

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            CDenseMatrix mat_chi(aocount, grid_batch_size);

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
                    std::memcpy(mat_chi.row(idx), submat_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");

            // generate density grid

            omptimers[thread_id].start("Generate density grid");

            auto local_weights = omp_local_weights_data[thread_id].data();

            auto rho = omp_rho_data[thread_id].data();

            auto v2rho2 = omp_v2rho2_data[thread_id].data();
            auto v3rho3 = omp_v3rho3_data[thread_id].data();
            auto v4rho4 = omp_v4rho4_data[thread_id].data();

            sdengridgen::serialGenerateDensityForLDA(rho, mat_chi, gs_sub_dens_mat);

            auto xcfuntype = omp_xcfuncs[thread_id].getFunctionalType();

            auto rwdengrid = sdengridgen::serialGenerateDensityGridForLDA(mat_chi, rw_sub_dens_mat, xcfuntype);

            auto rw2dengrid = sdengridgen::serialGenerateDensityGridForLDA(mat_chi, rw2_sub_dens_mat, xcfuntype);

            auto rw3dengrid = sdengridgen::serialGenerateDensityGridForLDA(mat_chi, rw3_sub_dens_mat, xcfuntype);

            // compute perturbed density

            omptimers[thread_id].start("Density grid cube");

            auto numdens_rw3 = static_cast<int>(rw3DensityPointers.size());

            auto numdens_rw2 = static_cast<int>(rw2DensityPointers.size());

            CDensityGridCubic rwdengridcube(grid_batch_size, numdens_rw2 + numdens_rw3, xcfuntype, dengrid::ab);

            rwdengridcube.DensityProd(rwdengrid, rw2dengrid, xcfuntype, (numdens_rw2 + numdens_rw3), cubeMode);

            omptimers[thread_id].stop("Density grid cubic");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_fxc_for_lda(grid_batch_size, rho, v2rho2);

            omp_xcfuncs[thread_id].compute_kxc_for_lda(grid_batch_size, rho, v3rho3);

            omp_xcfuncs[thread_id].compute_lxc_for_lda(grid_batch_size, rho, v4rho4);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            // go through density matrices

            for (int idensity = 0; idensity < numdens_rw2; idensity++)
            {
                // compute partial contribution to Kxc matrix

                auto rhow1a = rwdengridcube.gam2(idensity);

                auto rhow12a = rw2dengrid.alphaDensity(idensity);

                auto rhow12b = rw2dengrid.betaDensity(idensity);

                std::vector<const double*> rwdengrid_pointers({rhow1a});
                std::vector<const double*> rw2dengrid_pointers({rhow12a, rhow12b});

                auto partial_mat_Kxc = integratePartialKxcFockForLdaClosedShell(
                    omp_xcfuncs[thread_id], local_weights, mat_chi, v2rho2, v3rho3,
                    rwdengrid_pointers, rw2dengrid_pointers, omptimers[thread_id]);

                // accumulate partial Kxc

                omptimers[thread_id].start("Kxc local matrix dist.");

                #pragma omp critical
                sdenblas::serialInPlaceAddAB(sum_partial_mat_Kxc[idensity], partial_mat_Kxc);

                omptimers[thread_id].stop("Kxc local matrix dist.");
            }

            for (int idensity = 0; idensity < numdens_rw3; idensity++)
            {
                // compute partial contribution to Lxc matrix

                auto partial_mat_Lxc = integratePartialLxcFockForLdaClosedShell(
                    omp_xcfuncs[thread_id], local_weights, mat_chi, v2rho2, v3rho3, v4rho4, rwdengridcube, rw3dengrid, idensity, omptimers[thread_id]);

                // accumulate partial Lxc

                omptimers[thread_id].start("Lxc local matrix dist.");

                #pragma omp critical
                sdenblas::serialInPlaceAddAB(sum_partial_mat_Lxc[idensity], partial_mat_Lxc);

                omptimers[thread_id].stop("Lxc local matrix dist.");
            }
        }

        timer.stop("OMP KxcLxc calc.");

        timer.start("Kxc matrix dist.");

        for (size_t idensity = 0; idensity < rw2DensityPointers.size(); idensity++)
        {
            dftsubmat::distributeSubMatrixToFock(aoFockPointers, idensity, sum_partial_mat_Kxc[idensity], aoinds, naos);
        }

        timer.stop("Kxc matrix dist.");

        timer.start("Lxc matrix dist.");

        for (size_t idensity = 0; idensity < rw3DensityPointers.size(); idensity++)
        {
            dftsubmat::distributeSubMatrixToFock(aoFockPointers, idensity + rw2DensityPointers.size(), sum_partial_mat_Lxc[idensity], aoinds, naos);
        }

        timer.stop("Lxc matrix dist.");
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
integratePartialKxcFockForLdaClosedShell(const CXCFunctional&     xcFunctional,
                                         const double*            weights,
                                         const CDenseMatrix&      gtoValues,
                                         const double*            v2rho2,
                                         const double*            v3rho3,
                                         const std::vector<const double*>& rwDensityGridPointers,
                                         const std::vector<const double*>& rw2DensityGridPointers,
                                         CMultiTimer&             timer) -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto rhow1a = rwDensityGridPointers[0];

    auto rhow12a = rw2DensityGridPointers[0];
    auto rhow12b = rw2DensityGridPointers[1];

    timer.start("Kxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);

    auto G_val = mat_G.values();

    auto       ldafunc = xcFunctional.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    {
        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd 
            for (int g = 0; g < npoints; g++)
            {
                // functional derivatives

                // first-order

                // second-order

                auto v2rho2_aa = v2rho2[dim->v2rho2 * g + 0];
                auto v2rho2_ab = v2rho2[dim->v2rho2 * g + 1];

                // third-order

                auto v3rho3_aaa = v3rho3[dim->v3rho3 * g + 0];
                auto v3rho3_aab = v3rho3[dim->v3rho3 * g + 1];
                auto v3rho3_abb = v3rho3[dim->v3rho3 * g + 2];

                G_val[nu_offset + g] = weights[g] *
                                       ((v3rho3_aaa + 2.0 * v3rho3_aab + v3rho3_abb) * rhow1a[g] + v2rho2_aa * rhow12a[g] + v2rho2_ab * rhow12b[g]) *
                                       chi_val[nu_offset + g];
            }
        }
    }

    timer.stop("Kxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Kxc matrix matmul");

    auto mat_Kxc = sdenblas::serialMultABt(gtoValues, mat_G);

    timer.stop("Kxc matrix matmul");

    return mat_Kxc;
}

auto
integratePartialLxcFockForLdaClosedShell(const CXCFunctional&     xcFunctional,
                                         const double*            weights,
                                         const CDenseMatrix&      gtoValues,
                                         const double*            v2rho2,
                                         const double*            v3rho3,
                                         const double*            v4rho4,
                                         const CDensityGridCubic& rwDensityGridCubic,
                                         const CDensityGrid&      rw3DensityGrid,
                                         const int                iFock,
                                         CMultiTimer&             timer) -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    timer.start("Lxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    // pointers to perturbed density

    auto rx_ry_rz = rwDensityGridCubic.pi(iFock);
    auto rxy_rz   = rwDensityGridCubic.gam(iFock);
    auto r_xyz    = rw3DensityGrid.alphaDensity(iFock);

    timer.stop("Lxc matrix prep.");

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Lxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);

    auto G_val = mat_G.values();

    auto       ldafunc = xcFunctional.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    {
        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd 
            for (int g = 0; g < npoints; g++)
            {
                // functional derivatives

                // first-order

                // second-order

                auto v2rho2_aa = v2rho2[dim->v2rho2 * g + 0];
                auto v2rho2_ab = v2rho2[dim->v2rho2 * g + 1];

                // third-order

                auto v3rho3_aaa = v3rho3[dim->v3rho3 * g + 0];
                auto v3rho3_aab = v3rho3[dim->v3rho3 * g + 1];
                auto v3rho3_abb = v3rho3[dim->v3rho3 * g + 2];

                // fourth-order

                auto v4rho4_aaaa = v4rho4[dim->v4rho4 * g + 0];
                auto v4rho4_aaab = v4rho4[dim->v4rho4 * g + 1];
                auto v4rho4_aabb = v4rho4[dim->v4rho4 * g + 2];
                auto v4rho4_abbb = v4rho4[dim->v4rho4 * g + 3];

                double rr   = (v2rho2_aa + v2rho2_ab);
                double rrr  = (v3rho3_aaa + 2.0 * v3rho3_aab + v3rho3_abb);
                double rrrr = (v4rho4_aaaa + 3.0 * v4rho4_aaab + 3.0 * v4rho4_aabb + v4rho4_abbb);

                G_val[nu_offset + g] = weights[g] * (rr * r_xyz[g] + rrr * rxy_rz[g] + rrrr * rx_ry_rz[g]) * chi_val[nu_offset + g];
            }
        }
    }

    timer.stop("Lxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Lxc matrix matmul");

    auto mat_Lxc = sdenblas::serialMultABt(gtoValues, mat_G);

    timer.stop("Lxc matrix matmul");

    return mat_Lxc;
}

}  // namespace xcintlda
