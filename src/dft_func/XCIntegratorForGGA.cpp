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

#include "XCIntegratorForGGA.hpp"

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

namespace xcintgga {  // xcintgga namespace

auto
integrateVxcFockForGgaClosedShell(const CMolecule&                  molecule,
                                  const CMolecularBasis&            basis,
                                  const std::vector<const double*>& gsDensityPointers,
                                  const CMolecularGrid&             molecularGrid,
                                  const double                      screeningThresholdForGTOValues,
                                  const CXCFunctional&              xcFunctional) -> CAOKohnShamMatrix
{
    CMultiTimer timer;

    timer.start("Total timing");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

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

    const auto n_boxes = counts.size();

    const auto n_gto_blocks = gto_blocks.size();

    // set up pointers to OMP data

    auto ptr_counts = counts.data();

    auto ptr_displacements = displacements.data();

    auto ptr_gto_blocks = gto_blocks.data();

    auto ptr_gsDensityPointers = gsDensityPointers.data();

    auto ptr_xcFunctional = &xcFunctional;

#pragma omp parallel shared(ptr_counts, ptr_displacements, xcoords, ycoords, zcoords, \
                            ptr_gto_blocks, ptr_gsDensityPointers, ptr_xcFunctional, \
                            n_boxes, n_gto_blocks, naos, nele, xcene, mat_Vxc)
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

            auto sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

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

            auto local_xcfunc = CXCFunctional(*ptr_xcFunctional);

            auto       ggafunc = local_xcfunc.getFunctionalPointerToGgaComponent();
            const auto dim     = &(ggafunc->dim);

            std::vector<double> local_weights_data(weights + gridblockpos, weights + gridblockpos + npoints);

            auto rho_data     = std::vector<double>(dim->rho * npoints);
            auto rhograd_data = std::vector<double>(dim->rho * 3 * npoints);
            auto sigma_data   = std::vector<double>(dim->sigma * npoints);

            auto exc_data    = std::vector<double>(dim->zk * npoints);
            auto vrho_data   = std::vector<double>(dim->vrho * npoints);
            auto vsigma_data = std::vector<double>(dim->vsigma * npoints);

            auto local_weights = local_weights_data.data();

            auto rho     = rho_data.data();
            auto rhograd = rhograd_data.data();
            auto sigma   = sigma_data.data();

            auto exc    = exc_data.data();
            auto vrho   = vrho_data.data();
            auto vsigma = vsigma_data.data();

            sdengridgen::serialGenerateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat);

            omptimers[thread_id].stop("Generate density grid");

            omptimers[thread_id].start("XC functional eval.");

            local_xcfunc.compute_exc_vxc_for_gga(npoints, rho, sigma, exc, vrho, vsigma);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Vxc matrix G");

            CDenseMatrix mat_G(aocount, npoints);
            CDenseMatrix mat_G_gga(aocount, npoints);

            auto G_val     = mat_G.values();
            auto G_gga_val = mat_G_gga.values();

            auto chi_val   = mat_chi.values();
            auto chi_x_val = mat_chi_x.values();
            auto chi_y_val = mat_chi_y.values();
            auto chi_z_val = mat_chi_z.values();

            for (int nu = 0; nu < aocount; nu++)
            {
                auto nu_offset = nu * npoints;

                #pragma omp simd
                for (int g = 0; g < npoints; g++)
                {
                    auto vx = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                    auto vy = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                    auto vz = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                    G_val[nu_offset + g] = local_weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];

                    G_gga_val[nu_offset + g] =
                        local_weights[g] * (vx * chi_x_val[nu_offset + g] + vy * chi_y_val[nu_offset + g] + vz * chi_z_val[nu_offset + g]);
                }
            }

            omptimers[thread_id].stop("Vxc matrix G");

            omptimers[thread_id].start("Vxc matmul and symm.");

            auto partial_mat_Vxc = sdenblas::serialMultABt(mat_chi, sdenblas::serialAddAB(mat_G, mat_G_gga, 2.0));

            partial_mat_Vxc.symmetrizeAndScale(0.5);

            omptimers[thread_id].stop("Vxc matmul and symm.");

            omptimers[thread_id].start("Vxc dist.");

            double local_nele = 0.0, local_xcene = 0.0;

            for (int g = 0; g < npoints; g++)
            {
                auto rho_total = rho[2 * g + 0] + rho[2 * g + 1];

                local_nele += local_weights[g] * rho_total;

                local_xcene += local_weights[g] * exc[g] * rho_total;
            }

            #pragma omp critical
            {
                nele += local_nele;

                xcene += local_xcene;

                dftsubmat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc, aoinds);
            }

            omptimers[thread_id].stop("Vxc dist.");
        }
    }
    }
    }
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
integrateVxcFockForGgaOpenShell(const CMolecule&                  molecule,
                                const CMolecularBasis&            basis,
                                const std::vector<const double*>& gsDensityPointers,
                                const CMolecularGrid&             molecularGrid,
                                const double                      screeningThresholdForGTOValues,
                                const CXCFunctional&              xcFunctional) -> CAOKohnShamMatrix
{
    CMultiTimer timer;

    timer.start("Total timing");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // Kohn-Sham matrix

    CAOKohnShamMatrix mat_Vxc(naos, naos, std::string("openshell"));

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

    const auto n_boxes = counts.size();

    const auto n_gto_blocks = gto_blocks.size();

    // set up pointers to OMP data

    auto ptr_counts = counts.data();

    auto ptr_displacements = displacements.data();

    auto ptr_gto_blocks = gto_blocks.data();

    auto ptr_gsDensityPointers = gsDensityPointers.data();

    auto ptr_xcFunctional = &xcFunctional;

#pragma omp parallel shared(ptr_counts, ptr_displacements, xcoords, ycoords, zcoords, \
                            ptr_gto_blocks, ptr_gsDensityPointers, ptr_xcFunctional, \
                            n_boxes, n_gto_blocks, naos, nele, xcene, mat_Vxc)
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

            auto sub_dens_mat_a = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);
            auto sub_dens_mat_b = dftsubmat::getSubDensityMatrix(gsDensityPointers[1], aoinds, naos);

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

            auto local_xcfunc = CXCFunctional(*ptr_xcFunctional);

            auto       ggafunc = local_xcfunc.getFunctionalPointerToGgaComponent();
            const auto dim     = &(ggafunc->dim);

            std::vector<double> local_weights_data(weights + gridblockpos, weights + gridblockpos + npoints);

            std::vector<double> rho_data(dim->rho * npoints);
            std::vector<double> rhograd_data(dim->rho * 3 * npoints);
            std::vector<double> sigma_data(dim->sigma * npoints);

            std::vector<double> exc_data(dim->zk * npoints);
            std::vector<double> vrho_data(dim->vrho * npoints);
            std::vector<double> vsigma_data(dim->vsigma * npoints);

            auto local_weights = local_weights_data.data();

            auto rho     = rho_data.data();
            auto rhograd = rhograd_data.data();
            auto sigma   = sigma_data.data();

            auto exc    = exc_data.data();
            auto vrho   = vrho_data.data();
            auto vsigma = vsigma_data.data();

            sdengridgen::serialGenerateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat_a, sub_dens_mat_b);

            omptimers[thread_id].stop("Generate density grid");

            omptimers[thread_id].start("XC functional eval.");

            local_xcfunc.compute_exc_vxc_for_gga(npoints, rho, sigma, exc, vrho, vsigma);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Vxc matrix G");

            CDenseMatrix mat_G_a(aocount, npoints);
            CDenseMatrix mat_G_a_gga(aocount, npoints);

            CDenseMatrix mat_G_b(aocount, npoints);
            CDenseMatrix mat_G_b_gga(aocount, npoints);

            auto G_a_val     = mat_G_a.values();
            auto G_a_gga_val = mat_G_a_gga.values();

            auto G_b_val     = mat_G_b.values();
            auto G_b_gga_val = mat_G_b_gga.values();

            auto chi_val   = mat_chi.values();
            auto chi_x_val = mat_chi_x.values();
            auto chi_y_val = mat_chi_y.values();
            auto chi_z_val = mat_chi_z.values();

            for (int nu = 0; nu < aocount; nu++)
            {
                auto nu_offset = nu * npoints;

                #pragma omp simd
                for (int g = 0; g < npoints; g++)
                {
                    auto vxa = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                    auto vya = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                    auto vza = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                    auto vxb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 3] + vsigma[3 * g + 1] * rhograd[6 * g + 0];
                    auto vyb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 4] + vsigma[3 * g + 1] * rhograd[6 * g + 1];
                    auto vzb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 5] + vsigma[3 * g + 1] * rhograd[6 * g + 2];

                    G_a_val[nu_offset + g] = local_weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];
                    G_b_val[nu_offset + g] = local_weights[g] * vrho[2 * g + 1] * chi_val[nu_offset + g];

                    G_a_gga_val[nu_offset + g] =
                        local_weights[g] * (vxa * chi_x_val[nu_offset + g] + vya * chi_y_val[nu_offset + g] + vza * chi_z_val[nu_offset + g]);
                    G_b_gga_val[nu_offset + g] =
                        local_weights[g] * (vxb * chi_x_val[nu_offset + g] + vyb * chi_y_val[nu_offset + g] + vzb * chi_z_val[nu_offset + g]);
                }
            }

            omptimers[thread_id].stop("Vxc matrix G");

            omptimers[thread_id].start("Vxc matmul and symm.");

            auto partial_mat_Vxc_a = sdenblas::serialMultABt(mat_chi, sdenblas::serialAddAB(mat_G_a, mat_G_a_gga, 2.0));
            auto partial_mat_Vxc_b = sdenblas::serialMultABt(mat_chi, sdenblas::serialAddAB(mat_G_b, mat_G_b_gga, 2.0));

            partial_mat_Vxc_a.symmetrizeAndScale(0.5);
            partial_mat_Vxc_b.symmetrizeAndScale(0.5);

            omptimers[thread_id].stop("Vxc matmul and symm.");

            omptimers[thread_id].start("Vxc dist.");

            double local_nele = 0.0, local_xcene = 0.0;

            for (int g = 0; g < npoints; g++)
            {
                auto rho_total = rho[2 * g + 0] + rho[2 * g + 1];

                local_nele += local_weights[g] * rho_total;

                local_xcene += local_weights[g] * exc[g] * rho_total;
            }

            #pragma omp critical
            {
                nele += local_nele;

                xcene += local_xcene;

                dftsubmat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc_a, partial_mat_Vxc_b, aoinds);
            }

            omptimers[thread_id].stop("Vxc dist.");
        }
    }
    }
    }
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
integrateFxcFockForGgaClosedShell(const std::vector<double*>&       aoFockPointers,
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

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

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

    const auto n_rw_densities = rwDensityPointers.size();

    // set up pointers to OMP data

    auto ptr_counts = counts.data();

    auto ptr_displacements = displacements.data();

    auto ptr_gto_blocks = gto_blocks.data();

    auto ptr_gsDensityPointers = gsDensityPointers.data();

    auto ptr_xcFunctional = &xcFunctional;

#pragma omp parallel shared(ptr_counts, ptr_displacements, xcoords, ycoords, zcoords, \
                            ptr_gto_blocks, ptr_gsDensityPointers, ptr_xcFunctional, \
                            n_boxes, n_gto_blocks, n_rw_densities, naos, \
                            aoFockPointers, rwDensityPointers)
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

            auto sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

            std::vector<CDenseMatrix> rw_sub_dens_mat_vec(n_rw_densities);

            for (size_t idensity = 0; idensity < n_rw_densities; idensity++)
            {
                rw_sub_dens_mat_vec[idensity] = dftsubmat::getSubDensityMatrix(rwDensityPointers[idensity], aoinds, naos);
            }

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

            auto local_xcfunc = CXCFunctional(*ptr_xcFunctional);

            auto       ggafunc = local_xcfunc.getFunctionalPointerToGgaComponent();
            const auto dim     = &(ggafunc->dim);

            std::vector<double> local_weights_data(weights + gridblockpos, weights + gridblockpos + npoints);

            std::vector<double> rho_data(dim->rho * npoints);
            std::vector<double> rhograd_data(dim->rho * 3 * npoints);
            std::vector<double> sigma_data(dim->sigma * npoints);

            std::vector<double> rhow_data(dim->rho * npoints);
            std::vector<double> rhowgrad_data(dim->rho * 3 * npoints);

            std::vector<double> vrho_data(dim->vrho * npoints);
            std::vector<double> vsigma_data(dim->vsigma * npoints);

            std::vector<double> v2rho2_data(dim->v2rho2 * npoints);
            std::vector<double> v2rhosigma_data(dim->v2rhosigma * npoints);
            std::vector<double> v2sigma2_data(dim->v2sigma2 * npoints);

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

            sdengridgen::serialGenerateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat);

            omptimers[thread_id].stop("Generate density grid");

            omptimers[thread_id].start("XC functional eval.");

            local_xcfunc.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

            local_xcfunc.compute_fxc_for_gga(npoints, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

            omptimers[thread_id].stop("XC functional eval.");

            // go through rhow density matrices

            for (size_t idensity = 0; idensity < n_rw_densities; idensity++)
            {
                omptimers[thread_id].start("Generate density grid");

                sdengridgen::serialGenerateDensityForGGA(rhow, rhowgrad, nullptr, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw_sub_dens_mat_vec[idensity]);

                omptimers[thread_id].stop("Generate density grid");

                omptimers[thread_id].start("Fxc matrix G");

                CDenseMatrix mat_G(aocount, npoints);
                CDenseMatrix mat_G_gga(aocount, npoints);

                auto G_val     = mat_G.values();
                auto G_gga_val = mat_G_gga.values();

                auto chi_val   = mat_chi.values();
                auto chi_x_val = mat_chi_x.values();
                auto chi_y_val = mat_chi_y.values();
                auto chi_z_val = mat_chi_z.values();

                for (int nu = 0; nu < aocount; nu++)
                {
                    auto nu_offset = nu * npoints;

                    #pragma omp simd
                    for (int g = 0; g < npoints; g++)
                    {
                        double w = local_weights[g];

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

                omptimers[thread_id].stop("Fxc matrix G");

                omptimers[thread_id].start("Fxc matmul and symm.");

                auto partial_mat_Fxc = sdenblas::serialMultABt(mat_chi, mat_G);

                auto partial_mat_Fxc_gga = sdenblas::serialMultABt(mat_chi, mat_G_gga);

                partial_mat_Fxc_gga.symmetrize();  // matrix + matrix.T

                sdenblas::serialInPlaceAddAB(partial_mat_Fxc, partial_mat_Fxc_gga);

                omptimers[thread_id].stop("Fxc matmul and symm.");

                omptimers[thread_id].start("Fxc dist.");

                #pragma omp critical
                dftsubmat::distributeSubMatrixToFock(aoFockPointers, idensity, partial_mat_Fxc, aoinds, naos);

                omptimers[thread_id].stop("Fxc dist.");
            }
        }
    }
    }
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
integrateKxcFockForGgaClosedShell(const std::vector<double*>& aoFockPointers,
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

    auto omp_max_npoints = max_npoints_per_box / nthreads;
    if (max_npoints_per_box % nthreads != 0) omp_max_npoints++;

    // density and functional derivatives

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<CXCFunctional> omp_xcfuncs(nthreads, CXCFunctional(xcFunctional));

    std::vector<std::vector<double>> omp_local_weights_data(nthreads, std::vector<double>(omp_max_npoints));

    std::vector<std::vector<double>> omp_rho_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));
    std::vector<std::vector<double>> omp_rhograd_data(nthreads, std::vector<double>(dim->rho * 3 * omp_max_npoints));
    std::vector<std::vector<double>> omp_sigma_data(nthreads, std::vector<double>(dim->sigma * omp_max_npoints));

    std::vector<std::vector<double>> omp_vrho_data(nthreads, std::vector<double>(dim->vrho * omp_max_npoints));
    std::vector<std::vector<double>> omp_vsigma_data(nthreads, std::vector<double>(dim->vsigma * omp_max_npoints));

    std::vector<std::vector<double>> omp_v2rho2_data(nthreads, std::vector<double>(dim->v2rho2 * omp_max_npoints));
    std::vector<std::vector<double>> omp_v2rhosigma_data(nthreads, std::vector<double>(dim->v2rhosigma * omp_max_npoints));
    std::vector<std::vector<double>> omp_v2sigma2_data(nthreads, std::vector<double>(dim->v2sigma2 * omp_max_npoints));

    std::vector<std::vector<double>> omp_v3rho3_data(nthreads, std::vector<double>(dim->v3rho3 * omp_max_npoints));
    std::vector<std::vector<double>> omp_v3rho2sigma_data(nthreads, std::vector<double>(dim->v3rho2sigma * omp_max_npoints));
    std::vector<std::vector<double>> omp_v3rhosigma2_data(nthreads, std::vector<double>(dim->v3rhosigma2 * omp_max_npoints));
    std::vector<std::vector<double>> omp_v3sigma3_data(nthreads, std::vector<double>(dim->v3sigma3 * omp_max_npoints));

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

            auto vrho   = omp_vrho_data[thread_id].data();
            auto vsigma = omp_vsigma_data[thread_id].data();

            auto v2rho2     = omp_v2rho2_data[thread_id].data();
            auto v2rhosigma = omp_v2rhosigma_data[thread_id].data();
            auto v2sigma2   = omp_v2sigma2_data[thread_id].data();

            auto v3rho3      = omp_v3rho3_data[thread_id].data();
            auto v3rho2sigma = omp_v3rho2sigma_data[thread_id].data();
            auto v3rhosigma2 = omp_v3rhosigma2_data[thread_id].data();
            auto v3sigma3    = omp_v3sigma3_data[thread_id].data();

            // generate density grid

            sdengridgen::serialGenerateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat);

            auto xcfuntype = omp_xcfuncs[thread_id].getFunctionalType();

            auto rwdengrid = sdengridgen::serialGenerateDensityGridForGGA(mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw_sub_dens_mat, xcfuntype);

            auto rw2dengrid =
                sdengridgen::serialGenerateDensityGridForGGA(mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw2_sub_dens_mat, xcfuntype);

            omptimers[thread_id].stop("Generate density grid");

            // compute perturbed density

            omptimers[thread_id].start("Density grid quad");

            auto numdens_rw2 = static_cast<int>(rw2DensityPointers.size());

            CDensityGridQuad rwdengridquad(grid_batch_size, numdens_rw2, xcfuntype, dengrid::ab);

            rwdengridquad.DensityProd(rwdengrid, xcfuntype, numdens_rw2, quadMode);

            omptimers[thread_id].stop("Density grid quad");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_vxc_for_gga(grid_batch_size, rho, sigma, vrho, vsigma);

            omp_xcfuncs[thread_id].compute_fxc_for_gga(grid_batch_size, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

            omp_xcfuncs[thread_id].compute_kxc_for_gga(grid_batch_size, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            // go through density matrices

            for (int idensity = 0; idensity < numdens_rw2; idensity++)
            {
                // compute partial contribution to Kxc matrix

                auto gam = rwdengridquad.gam(idensity);

                auto gamx = rwdengridquad.gamX(idensity);
                auto gamy = rwdengridquad.gamY(idensity);
                auto gamz = rwdengridquad.gamZ(idensity);

                auto gamxx = rwdengridquad.gamXX(idensity);
                auto gamxy = rwdengridquad.gamXY(idensity);
                auto gamxz = rwdengridquad.gamXZ(idensity);

                auto gamyx = rwdengridquad.gamYX(idensity);
                auto gamyy = rwdengridquad.gamYY(idensity);
                auto gamyz = rwdengridquad.gamYZ(idensity);

                auto gamzx = rwdengridquad.gamZX(idensity);
                auto gamzy = rwdengridquad.gamZY(idensity);
                auto gamzz = rwdengridquad.gamZZ(idensity);

                auto rhow12a = rw2dengrid.alphaDensity(idensity);

                auto gradw12a_x = rw2dengrid.alphaDensityGradientX(idensity);
                auto gradw12a_y = rw2dengrid.alphaDensityGradientY(idensity);
                auto gradw12a_z = rw2dengrid.alphaDensityGradientZ(idensity);

                std::vector<const double*> rwdengrid_pointers({gam,
                                                               gamx, gamy, gamz,
                                                               gamxx, gamxy, gamxz,
                                                               gamyx, gamyy, gamyz,
                                                               gamzx, gamzy, gamzz});

                std::vector<const double*> rw2dengrid_pointers({rhow12a,
                                                                gradw12a_x, gradw12a_y, gradw12a_z});

                auto partial_mat_Kxc = integratePartialKxcFockForGgaClosedShell(omp_xcfuncs[thread_id],
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
                                                                  rwdengrid_pointers,
                                                                  rw2dengrid_pointers,
                                                                  omptimers[thread_id]);

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
integrateKxcLxcFockForGgaClosedShell(const std::vector<double*>& aoFockPointers,
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

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<CXCFunctional> omp_xcfuncs(nthreads, CXCFunctional(xcFunctional));

    std::vector<std::vector<double>> omp_local_weights_data(nthreads, std::vector<double>(omp_max_npoints));

    std::vector<std::vector<double>> omp_rho_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));
    std::vector<std::vector<double>> omp_rhograd_data(nthreads, std::vector<double>(dim->rho * 3 * omp_max_npoints));
    std::vector<std::vector<double>> omp_sigma_data(nthreads, std::vector<double>(dim->sigma * omp_max_npoints));

    std::vector<std::vector<double>> omp_vrho_data(nthreads, std::vector<double>(dim->vrho * omp_max_npoints));
    std::vector<std::vector<double>> omp_vsigma_data(nthreads, std::vector<double>(dim->vsigma * omp_max_npoints));

    std::vector<std::vector<double>> omp_v2rho2_data(nthreads, std::vector<double>(dim->v2rho2 * omp_max_npoints));
    std::vector<std::vector<double>> omp_v2rhosigma_data(nthreads, std::vector<double>(dim->v2rhosigma * omp_max_npoints));
    std::vector<std::vector<double>> omp_v2sigma2_data(nthreads, std::vector<double>(dim->v2sigma2 * omp_max_npoints));

    std::vector<std::vector<double>> omp_v3rho3_data(nthreads, std::vector<double>(dim->v3rho3 * omp_max_npoints));
    std::vector<std::vector<double>> omp_v3rho2sigma_data(nthreads, std::vector<double>(dim->v3rho2sigma * omp_max_npoints));
    std::vector<std::vector<double>> omp_v3rhosigma2_data(nthreads, std::vector<double>(dim->v3rhosigma2 * omp_max_npoints));
    std::vector<std::vector<double>> omp_v3sigma3_data(nthreads, std::vector<double>(dim->v3sigma3 * omp_max_npoints));

    std::vector<std::vector<double>> omp_v4rho4_data(nthreads, std::vector<double>(dim->v4rho4 * omp_max_npoints));
    std::vector<std::vector<double>> omp_v4rho3sigma_data(nthreads, std::vector<double>(dim->v4rho3sigma * omp_max_npoints));
    std::vector<std::vector<double>> omp_v4rho2sigma2_data(nthreads, std::vector<double>(dim->v4rho2sigma2 * omp_max_npoints));
    std::vector<std::vector<double>> omp_v4rhosigma3_data(nthreads, std::vector<double>(dim->v4rhosigma3 * omp_max_npoints));
    std::vector<std::vector<double>> omp_v4sigma4_data(nthreads, std::vector<double>(dim->v4sigma4 * omp_max_npoints));

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

            // generate density grid

            omptimers[thread_id].start("Generate density grid");

            auto local_weights = omp_local_weights_data[thread_id].data();

            auto rho     = omp_rho_data[thread_id].data();
            auto rhograd = omp_rhograd_data[thread_id].data();
            auto sigma   = omp_sigma_data[thread_id].data();

            auto vrho   = omp_vrho_data[thread_id].data();
            auto vsigma = omp_vsigma_data[thread_id].data();

            auto v2rho2     = omp_v2rho2_data[thread_id].data();
            auto v2rhosigma = omp_v2rhosigma_data[thread_id].data();
            auto v2sigma2   = omp_v2sigma2_data[thread_id].data();

            auto v3rho3      = omp_v3rho3_data[thread_id].data();
            auto v3rho2sigma = omp_v3rho2sigma_data[thread_id].data();
            auto v3rhosigma2 = omp_v3rhosigma2_data[thread_id].data();
            auto v3sigma3    = omp_v3sigma3_data[thread_id].data();

            auto v4rho4       = omp_v4rho4_data[thread_id].data();
            auto v4rho3sigma  = omp_v4rho3sigma_data[thread_id].data();
            auto v4rho2sigma2 = omp_v4rho2sigma2_data[thread_id].data();
            auto v4rhosigma3  = omp_v4rhosigma3_data[thread_id].data();
            auto v4sigma4     = omp_v4sigma4_data[thread_id].data();

            sdengridgen::serialGenerateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat);

            auto xcfuntype = omp_xcfuncs[thread_id].getFunctionalType();

            auto rwdengrid =
                sdengridgen::serialGenerateDensityGridForGGA(mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw_sub_dens_mat, xcfuntype);

            auto rw2dengrid =
                sdengridgen::serialGenerateDensityGridForGGA(mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw2_sub_dens_mat, xcfuntype);

            auto rw3dengrid =
                sdengridgen::serialGenerateDensityGridForGGA(mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw3_sub_dens_mat, xcfuntype);

            omptimers[thread_id].stop("Generate density grid");

            // compute perturbed density

            omptimers[thread_id].start("Density grid cubic");

            auto numdens_rw3 = static_cast<int>(rw3DensityPointers.size());

            auto numdens_rw2 = static_cast<int>(rw2DensityPointers.size());

            CDensityGridCubic rwdengridcube(grid_batch_size, numdens_rw2 + numdens_rw3, xcfuntype, dengrid::ab);

            rwdengridcube.DensityProd(rwdengrid, rw2dengrid, xcfuntype, (numdens_rw2 + numdens_rw3), cubeMode);

            omptimers[thread_id].stop("Density grid cubic");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_vxc_for_gga(grid_batch_size, rho, sigma, vrho, vsigma);

            omp_xcfuncs[thread_id].compute_fxc_for_gga(grid_batch_size, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

            omp_xcfuncs[thread_id].compute_kxc_for_gga(grid_batch_size, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);

            omp_xcfuncs[thread_id].compute_lxc_for_gga(grid_batch_size, rho, sigma, v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            // go through density matrices

            for (int idensity = 0; idensity < numdens_rw2; idensity++)
            {
                // compute partial contribution to Kxc matrix

                auto rhow1rhow2 = rwdengridcube.gam2(idensity);

                auto rxw1rhow2 = rwdengridcube.gam2X(idensity);
                auto ryw1rhow2 = rwdengridcube.gam2Y(idensity);
                auto rzw1rhow2 = rwdengridcube.gam2Z(idensity);

                auto rxw1rxw2 = rwdengridcube.gam2XX(idensity);
                auto rxw1ryw2 = rwdengridcube.gam2XY(idensity);
                auto rxw1rzw2 = rwdengridcube.gam2XZ(idensity);

                auto ryw1rxw2 = rwdengridcube.gam2YX(idensity);
                auto ryw1ryw2 = rwdengridcube.gam2YY(idensity);
                auto ryw1rzw2 = rwdengridcube.gam2YZ(idensity);

                auto rzw1rxw2 = rwdengridcube.gam2ZX(idensity);
                auto rzw1ryw2 = rwdengridcube.gam2ZY(idensity);
                auto rzw1rzw2 = rwdengridcube.gam2ZZ(idensity);

                auto rhow12a = rw2dengrid.alphaDensity(idensity);

                auto gradw12a_x = rw2dengrid.alphaDensityGradientX(idensity);
                auto gradw12a_y = rw2dengrid.alphaDensityGradientY(idensity);
                auto gradw12a_z = rw2dengrid.alphaDensityGradientZ(idensity);

                std::vector<const double*> rwdengrid_pointers({rhow1rhow2,
                                                               rxw1rhow2, ryw1rhow2, rzw1rhow2,
                                                               rxw1rxw2, rxw1ryw2, rxw1rzw2, 
                                                               ryw1rxw2, ryw1ryw2, ryw1rzw2, 
                                                               rzw1rxw2, rzw1ryw2, rzw1rzw2});

                std::vector<const double*> rw2dengrid_pointers({rhow12a,
                                                                gradw12a_x, gradw12a_y, gradw12a_z});

                auto partial_mat_Kxc = integratePartialKxcFockForGgaClosedShell(omp_xcfuncs[thread_id],
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
                                                                  rwdengrid_pointers,
                                                                  rw2dengrid_pointers,
                                                                  omptimers[thread_id]);

                // accumulate partial Kxc

                omptimers[thread_id].start("Kxc local matrix dist.");

                #pragma omp critical
                sdenblas::serialInPlaceAddAB(sum_partial_mat_Kxc[idensity], partial_mat_Kxc);

                omptimers[thread_id].stop("Kxc local matrix dist.");
            }

            for (int idensity = 0; idensity < numdens_rw3; idensity++)
            {
                // compute partial contribution to Lxc matrix

                auto partial_mat_Lxc = integratePartialLxcFockForGgaClosedShell(omp_xcfuncs[thread_id],
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
                                                                 v4rho4,
                                                                 v4rho3sigma,
                                                                 v4rho2sigma2,
                                                                 v4rhosigma3,
                                                                 v4sigma4,
                                                                 rwdengridcube,
                                                                 rw3dengrid,
                                                                 idensity,
                                                                 omptimers[thread_id]);

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
integratePartialKxcFockForGgaClosedShell(const CXCFunctional&     xcFunctional,
                                         const double*            weights,
                                         const CDenseMatrix&      gtoValues,
                                         const CDenseMatrix&      gtoValuesX,
                                         const CDenseMatrix&      gtoValuesY,
                                         const CDenseMatrix&      gtoValuesZ,
                                         const double*            rhograd,
                                         const double*            vsigma,
                                         const double*            v2rho2,
                                         const double*            v2rhosigma,
                                         const double*            v2sigma2,
                                         const double*            v3rho3,
                                         const double*            v3rho2sigma,
                                         const double*            v3rhosigma2,
                                         const double*            v3sigma3,
                                         const std::vector<const double*>& rwDensityGridPointers,
                                         const std::vector<const double*>& rw2DensityGridPointers,
                                         CMultiTimer&             timer) -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    // pointers to perturbed density gradient norms

    auto rhow1rhow2 = rwDensityGridPointers[0];

    auto rxw1rhow2 = rwDensityGridPointers[1];
    auto ryw1rhow2 = rwDensityGridPointers[2];
    auto rzw1rhow2 = rwDensityGridPointers[3];

    auto rxw1rxw2 = rwDensityGridPointers[4];
    auto rxw1ryw2 = rwDensityGridPointers[5];
    auto rxw1rzw2 = rwDensityGridPointers[6];

    auto ryw1rxw2 = rwDensityGridPointers[7];
    auto ryw1ryw2 = rwDensityGridPointers[8];
    auto ryw1rzw2 = rwDensityGridPointers[9];

    auto rzw1rxw2 = rwDensityGridPointers[10];
    auto rzw1ryw2 = rwDensityGridPointers[11];
    auto rzw1rzw2 = rwDensityGridPointers[12];

    auto rhow12a = rw2DensityGridPointers[0];

    auto gradw12a_x = rw2DensityGridPointers[1];
    auto gradw12a_y = rw2DensityGridPointers[2];
    auto gradw12a_z = rw2DensityGridPointers[3];

    timer.start("Kxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);
    CDenseMatrix mat_G_gga(naos, npoints);

    auto G_val     = mat_G.values();
    auto G_gga_val = mat_G_gga.values();

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    {
        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd 
            for (int g = 0; g < npoints; g++)
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
                double q2contract   = grada_x_g * rxw1rhow2[g] + grada_y_g * ryw1rhow2[g] + grada_z_g * rzw1rhow2[g];
                double q3contract = grada_x_g * grada_x_g * rxw1rxw2[g] + grada_x_g * grada_y_g * rxw1ryw2[g] + grada_x_g * grada_z_g * rxw1rzw2[g] +
                                    grada_y_g * grada_x_g * ryw1rxw2[g] + grada_y_g * grada_y_g * ryw1ryw2[g] + grada_y_g * grada_z_g * ryw1rzw2[g] +
                                    grada_z_g * grada_x_g * rzw1rxw2[g] + grada_z_g * grada_y_g * rzw1ryw2[g] + grada_z_g * grada_z_g * rzw1rzw2[g];

                double q4contract = rxw1rxw2[g] + ryw1ryw2[g] + rzw1rzw2[g];
                double q7contract_x =
                    grada_x_g * grada_x_g * rxw1rhow2[g] + grada_x_g * grada_y_g * ryw1rhow2[g] + grada_x_g * grada_z_g * rzw1rhow2[g];
                double q7contract_y =
                    grada_y_g * grada_x_g * rxw1rhow2[g] + grada_y_g * grada_y_g * ryw1rhow2[g] + grada_y_g * grada_z_g * rzw1rhow2[g];
                double q7contract_z =
                    grada_z_g * grada_x_g * rxw1rhow2[g] + grada_z_g * grada_y_g * ryw1rhow2[g] + grada_z_g * grada_z_g * rzw1rhow2[g];

                double q8contract_x = grada_x_g * rxw1rxw2[g] + grada_y_g * rxw1ryw2[g] + grada_z_g * rxw1rzw2[g];
                double q8contract_y = grada_x_g * ryw1rxw2[g] + grada_y_g * ryw1ryw2[g] + grada_z_g * ryw1rzw2[g];
                double q8contract_z = grada_x_g * rzw1rxw2[g] + grada_y_g * rzw1ryw2[g] + grada_z_g * rzw1rzw2[g];

                double q9contract_x = grada_x_g * q3contract;
                double q9contract_y = grada_y_g * q3contract;
                double q9contract_z = grada_z_g * q3contract;

                double q10contract_x = grada_x_g * rxw1rxw2[g] + grada_y_g * ryw1rxw2[g] + grada_z_g * rzw1rxw2[g];
                double q10contract_y = grada_x_g * rxw1ryw2[g] + grada_y_g * ryw1ryw2[g] + grada_z_g * rzw1ryw2[g];
                double q10contract_z = grada_x_g * rxw1rzw2[g] + grada_y_g * ryw1rzw2[g] + grada_z_g * rzw1rzw2[g];

                double q11contract_x = grada_x_g * rxw1rxw2[g] + grada_x_g * ryw1ryw2[g] + grada_x_g * rzw1rzw2[g];
                double q11contract_y = grada_y_g * rxw1rxw2[g] + grada_y_g * ryw1ryw2[g] + grada_y_g * rzw1rzw2[g];
                double q11contract_z = grada_z_g * rxw1rxw2[g] + grada_z_g * ryw1ryw2[g] + grada_z_g * rzw1rzw2[g];

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

                prefac += rrr * rhow1rhow2[g]  // q1
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

                xcomp += xrr * grada_x_g * rhow1rhow2[g]  // q5
                         + xr * rxw1rhow2[g]              // q6
                         + xxr * q7contract_x + xx * (q8contract_x + q10contract_x + q11contract_x) + xxx * q9contract_x;

                ycomp += xrr * grada_y_g * rhow1rhow2[g]  // q5
                         + xr * ryw1rhow2[g]              // q6
                         + xxr * q7contract_y + xx * (q8contract_y + q10contract_y + q11contract_y) + xxx * q9contract_y;

                zcomp += xrr * grada_z_g * rhow1rhow2[g]  // q5
                         + xr * rzw1rhow2[g]              // q6
                         + xxr * q7contract_z + xx * (q8contract_z + q10contract_z + q11contract_z) + xxx * q9contract_z;

                G_gga_val[nu_offset + g] =
                    w * (xcomp * chi_x_val[nu_offset + g] + ycomp * chi_y_val[nu_offset + g] + zcomp * chi_z_val[nu_offset + g]);
            }
        }
    }

    timer.stop("Kxc matrix G");

    timer.start("Kxc matrix matmul");

    auto mat_Kxc = sdenblas::serialMultABt(gtoValues, mat_G);

    auto mat_Kxc_gga = sdenblas::serialMultABt(gtoValues, mat_G_gga);

    mat_Kxc_gga.symmetrize();  // matrix + matrix.T

    sdenblas::serialInPlaceAddAB(mat_Kxc, mat_Kxc_gga);

    timer.stop("Kxc matrix matmul");

    return mat_Kxc;
}

auto
integratePartialLxcFockForGgaClosedShell(const CXCFunctional&     xcFunctional,
                              const double*            weights,
                              const CDenseMatrix&      gtoValues,
                              const CDenseMatrix&      gtoValuesX,
                              const CDenseMatrix&      gtoValuesY,
                              const CDenseMatrix&      gtoValuesZ,
                              const double*            rhograd,
                              const double*            vsigma,
                              const double*            v2rho2,
                              const double*            v2rhosigma,
                              const double*            v2sigma2,
                              const double*            v3rho3,
                              const double*            v3rho2sigma,
                              const double*            v3rhosigma2,
                              const double*            v3sigma3,
                              const double*            v4rho4,
                              const double*            v4rho3sigma,
                              const double*            v4rho2sigma2,
                              const double*            v4rhosigma3,
                              const double*            v4sigma4,
                              const CDensityGridCubic& rwDensityGridCubic,
                              const CDensityGrid&      rw3DensityGrid,
                              const int            iFock,
                              CMultiTimer&             timer) -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    // pointers to perturbed density gradient norms

    auto gam = rwDensityGridCubic.gam(iFock);

    auto gamx = rwDensityGridCubic.gamX(iFock);
    auto gamy = rwDensityGridCubic.gamY(iFock);
    auto gamz = rwDensityGridCubic.gamZ(iFock);

    auto gamxx = rwDensityGridCubic.gamXX(iFock);
    auto gamxy = rwDensityGridCubic.gamXY(iFock);
    auto gamxz = rwDensityGridCubic.gamXZ(iFock);

    auto gamyx = rwDensityGridCubic.gamYX(iFock);
    auto gamyy = rwDensityGridCubic.gamYY(iFock);
    auto gamyz = rwDensityGridCubic.gamYZ(iFock);

    auto gamzx = rwDensityGridCubic.gamZX(iFock);
    auto gamzy = rwDensityGridCubic.gamZY(iFock);
    auto gamzz = rwDensityGridCubic.gamZZ(iFock);

    auto pi = rwDensityGridCubic.pi(iFock);

    auto pix = rwDensityGridCubic.piX(iFock);
    auto piy = rwDensityGridCubic.piY(iFock);
    auto piz = rwDensityGridCubic.piZ(iFock);

    auto pixx = rwDensityGridCubic.piXX(iFock);
    auto pixy = rwDensityGridCubic.piXY(iFock);
    auto pixz = rwDensityGridCubic.piXZ(iFock);

    auto piyx = rwDensityGridCubic.piYX(iFock);
    auto piyy = rwDensityGridCubic.piYY(iFock);
    auto piyz = rwDensityGridCubic.piYZ(iFock);

    auto pizx = rwDensityGridCubic.piZX(iFock);
    auto pizy = rwDensityGridCubic.piZY(iFock);
    auto pizz = rwDensityGridCubic.piZZ(iFock);

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

    auto rho3 = rw3DensityGrid.alphaDensity(iFock);

    auto grad3_x = rw3DensityGrid.alphaDensityGradientX(iFock);
    auto grad3_y = rw3DensityGrid.alphaDensityGradientY(iFock);
    auto grad3_z = rw3DensityGrid.alphaDensityGradientZ(iFock);

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Lxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);
    CDenseMatrix mat_G_gga(naos, npoints);

    auto G_val     = mat_G.values();
    auto G_gga_val = mat_G_gga.values();

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    {
        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd 
            for (int g = 0; g < npoints; g++)
            {
                double w = weights[g];

                double rxw123a = grad3_x[g];
                double ryw123a = grad3_y[g];
                double rzw123a = grad3_z[g];

                double grada_x_g = rhograd[6 * g + 0];
                double grada_y_g = rhograd[6 * g + 1];
                double grada_z_g = rhograd[6 * g + 2];

                double l2contract   = grada_x_g * rxw123a + grada_y_g * ryw123a + grada_z_g * rzw123a;
                double l5contract_x = grada_x_g * l2contract;
                double l5contract_y = grada_y_g * l2contract;
                double l5contract_z = grada_z_g * l2contract;

                // vx2 terms
                double q2contract = grada_x_g * gamx[g] + grada_y_g * gamy[g] + grada_z_g * gamz[g];

                double q3contract = grada_x_g * grada_x_g * gamxx[g] + grada_x_g * grada_y_g * gamxy[g] + grada_x_g * grada_z_g * gamxz[g] +
                                    grada_y_g * grada_x_g * gamyx[g] + grada_y_g * grada_y_g * gamyy[g] + grada_y_g * grada_z_g * gamyz[g] +
                                    grada_z_g * grada_x_g * gamzx[g] + grada_z_g * grada_y_g * gamzy[g] + grada_z_g * grada_z_g * gamzz[g];

                double q4contract = gamxx[g] + gamyy[g] + gamzz[g];

                double q7contract_x = grada_x_g * q2contract;
                double q7contract_y = grada_y_g * q2contract;
                double q7contract_z = grada_z_g * q2contract;

                double q8contract_x = grada_x_g * gamxx[g] + grada_y_g * gamxy[g] + grada_z_g * gamxz[g];
                double q8contract_y = grada_x_g * gamyx[g] + grada_y_g * gamyy[g] + grada_z_g * gamyz[g];
                double q8contract_z = grada_x_g * gamzx[g] + grada_y_g * gamzy[g] + grada_z_g * gamzz[g];

                double q10contract_x = grada_x_g * gamxx[g] + grada_y_g * gamyx[g] + grada_z_g * gamzx[g];
                double q10contract_y = grada_x_g * gamxy[g] + grada_y_g * gamyy[g] + grada_z_g * gamzy[g];
                double q10contract_z = grada_x_g * gamxz[g] + grada_y_g * gamyz[g] + grada_z_g * gamzz[g];

                double q11contract_x = grada_x_g * gamxx[g] + grada_x_g * gamyy[g] + grada_x_g * gamzz[g];
                double q11contract_y = grada_y_g * gamxx[g] + grada_y_g * gamyy[g] + grada_y_g * gamzz[g];
                double q11contract_z = grada_z_g * gamxx[g] + grada_z_g * gamyy[g] + grada_z_g * gamzz[g];

                double q9contract_x = grada_x_g * q3contract;
                double q9contract_y = grada_y_g * q3contract;
                double q9contract_z = grada_z_g * q3contract;

                // vx3 terms
                double c1 = pi[g];

                double c2 = grada_x_g * pix[g] + grada_y_g * piy[g] + grada_z_g * piz[g];

                double c3 = grada_x_g * grada_x_g * pixx[g] + grada_x_g * grada_y_g * pixy[g] + grada_x_g * grada_z_g * pixz[g] +
                            grada_y_g * grada_x_g * piyx[g] + grada_y_g * grada_y_g * piyy[g] + grada_y_g * grada_z_g * piyz[g] +
                            grada_z_g * grada_x_g * pizx[g] + grada_z_g * grada_y_g * pizy[g] + grada_z_g * grada_z_g * pizz[g];

                double c4 = pixx[g] + piyy[g] + pizz[g];

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

                double c9_x = grada_x_g * pi[g];
                double c9_y = grada_y_g * pi[g];
                double c9_z = grada_z_g * pi[g];

                double c10_x = pix[g];
                double c10_y = piy[g];
                double c10_z = piz[g];

                double c11_x = c2 * grada_x_g;
                double c11_y = c2 * grada_y_g;
                double c11_z = c2 * grada_z_g;

                double c12_c14_x = grada_x_g * (pixx[g] + pixx[g]) + grada_y_g * (pixy[g] + piyx[g]) + grada_z_g * (pixz[g] + pizx[g]);

                double c12_c14_y = grada_x_g * (piyx[g] + pixy[g]) + grada_y_g * (piyy[g] + piyy[g]) + grada_z_g * (piyz[g] + pizy[g]);

                double c12_c14_z = grada_x_g * (pizx[g] + pixz[g]) + grada_y_g * (pizy[g] + piyz[g]) + grada_z_g * (pizz[g] + pizz[g]);

                double c13 = grada_x_g * grada_x_g * pixx[g] + grada_x_g * grada_y_g * pixy[g] + grada_x_g * grada_z_g * pixz[g] +
                             grada_y_g * grada_x_g * piyx[g] + grada_y_g * grada_y_g * piyy[g] + grada_y_g * grada_z_g * piyz[g] +
                             grada_z_g * grada_x_g * pizx[g] + grada_z_g * grada_y_g * pizy[g] + grada_z_g * grada_z_g * pizz[g];

                double c13_x = c13 * grada_x_g;
                double c13_y = c13 * grada_y_g;
                double c13_z = c13 * grada_z_g;

                double c15_x = grada_x_g * c4;
                double c15_y = grada_y_g * c4;
                double c15_z = grada_z_g * c4;

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

                // Transformation of derivatives

                // first-order
                double x = vsigma_c + 2.0 * vsigma_a;

                // second-order
                double rr = (v2rho2_aa + v2rho2_ab);
                double rx = (2.0 * v2rhosigma_ac + 2.0 * v2rhosigma_ab + 2.0 * v2rhosigma_aa);
                double xr = v2rhosigma_bc + 2.0 * v2rhosigma_ba + v2rhosigma_ac + 2.0 * v2rhosigma_aa;
                double xx = 2.0 * v2sigma2_cc + 2.0 * v2sigma2_cb + 6.0 * v2sigma2_ac + 4.0 * v2sigma2_ab + 4.0 * v2sigma2_aa;

                // third-order
                double rrr = (v3rho3_aaa + 2.0 * v3rho3_aab + v3rho3_abb);
                double rxr = (2.0 * v3rho2sigma_abc + 2.0 * v3rho2sigma_abb + 2.0 * v3rho2sigma_aba + 2.0 * v3rho2sigma_aac + 2.0 * v3rho2sigma_aab +
                              2.0 * v3rho2sigma_aaa);
                double rxx = (4.0 * v3rhosigma2_acc + 8.0 * v3rhosigma2_acb + 4.0 * v3rhosigma2_abb + 8.0 * v3rhosigma2_aac + 8.0 * v3rhosigma2_aab +
                              4.0 * v3rhosigma2_aaa);
                double xrr =
                    v3rho2sigma_bbc + 2.0 * v3rho2sigma_bba + 2.0 * v3rho2sigma_abc + 4.0 * v3rho2sigma_aba + v3rho2sigma_aac + 2.0 * v3rho2sigma_aaa;
                double xxr = 2.0 * v3rhosigma2_bcc + 2.0 * v3rhosigma2_bcb + 6.0 * v3rhosigma2_bac + 4.0 * v3rhosigma2_bab + 4.0 * v3rhosigma2_baa +
                             2.0 * v3rhosigma2_acc + 2.0 * v3rhosigma2_acb + 6.0 * v3rhosigma2_aac + 4.0 * v3rhosigma2_aab + 4.0 * v3rhosigma2_aaa;
                double xxx = 4.0 * v3sigma3_ccc + 8.0 * v3sigma3_ccb + 4.0 * v3sigma3_cbb + 16.0 * v3sigma3_acc + 24.0 * v3sigma3_acb +
                             8.0 * v3sigma3_abb + 20.0 * v3sigma3_aac + 16.0 * v3sigma3_aab + 8.0 * v3sigma3_aaa;

                // fourth-order

                double xrrr = v4rho3sigma_bbbc + 2.0 * v4rho3sigma_bbba + 3.0 * v4rho3sigma_abbc + 6.0 * v4rho3sigma_abba + 3.0 * v4rho3sigma_aabc +
                              6.0 * v4rho3sigma_aaba + v4rho3sigma_aaac + 2.0 * v4rho3sigma_aaaa;

                double xxrr = 2.0 * v4rho2sigma2_bbcc + 2.0 * v4rho2sigma2_bbcb + 6.0 * v4rho2sigma2_bbac + 4.0 * v4rho2sigma2_bbab +
                              4.0 * v4rho2sigma2_bbaa + 4.0 * v4rho2sigma2_abcc + 4.0 * v4rho2sigma2_abcb + 12.0 * v4rho2sigma2_abac +
                              8.0 * v4rho2sigma2_abab + 8.0 * v4rho2sigma2_abaa + 2.0 * v4rho2sigma2_aacc + 2.0 * v4rho2sigma2_aacb +
                              6.0 * v4rho2sigma2_aaac + 4.0 * v4rho2sigma2_aaab + 4.0 * v4rho2sigma2_aaaa;

                double xxxr = 4.0 * v4rhosigma3_bccc + 8.0 * v4rhosigma3_bccb + 4.0 * v4rhosigma3_bcbb + 16.0 * v4rhosigma3_bacc +
                              24.0 * v4rhosigma3_bacb + 8.0 * v4rhosigma3_babb + 20.0 * v4rhosigma3_baac + 16.0 * v4rhosigma3_baab +
                              8.0 * v4rhosigma3_baaa + 4.0 * v4rhosigma3_accc + 8.0 * v4rhosigma3_accb + 4.0 * v4rhosigma3_acbb +
                              16.0 * v4rhosigma3_aacc + 24.0 * v4rhosigma3_aacb + 8.0 * v4rhosigma3_aabb + 20.0 * v4rhosigma3_aaac +
                              16.0 * v4rhosigma3_aaab + 8.0 * v4rhosigma3_aaaa;

                double xxxx = 8.0 * v4sigma4_cccc + 24.0 * v4sigma4_cccb + 24.0 * v4sigma4_ccbb + 8.0 * v4sigma4_cbbb + 40.0 * v4sigma4_accc +
                              96.0 * v4sigma4_accb + 72.0 * v4sigma4_acbb + 16.0 * v4sigma4_abbb + 72.0 * v4sigma4_aacc + 120.0 * v4sigma4_aacb +
                              48.0 * v4sigma4_aabb + 56.0 * v4sigma4_aaac + 48.0 * v4sigma4_aaab + 16.0 * v4sigma4_aaaa;

                double rrrr = v4rho4_aaaa + 3.0 * v4rho4_aaab + 3.0 * v4rho4_aabb + v4rho4_abbb;

                double rxrr = 2.0 * v4rho3sigma_abbc + 2.0 * v4rho3sigma_abbb + 2.0 * v4rho3sigma_abba + 4.0 * v4rho3sigma_aabc +
                              4.0 * v4rho3sigma_aabb + 4.0 * v4rho3sigma_aaba + 2.0 * v4rho3sigma_aaac + 2.0 * v4rho3sigma_aaab +
                              2.0 * v4rho3sigma_aaaa;

                double rxxr = 4.0 * v4rho2sigma2_abcc + 8.0 * v4rho2sigma2_abcb + 4.0 * v4rho2sigma2_abbb + 8.0 * v4rho2sigma2_abac +
                              8.0 * v4rho2sigma2_abab + 4.0 * v4rho2sigma2_abaa + 4.0 * v4rho2sigma2_aacc + 8.0 * v4rho2sigma2_aacb +
                              4.0 * v4rho2sigma2_aabb + 8.0 * v4rho2sigma2_aaac + 8.0 * v4rho2sigma2_aaab + 4.0 * v4rho2sigma2_aaaa;

                double rxxx = 8.0 * v4rhosigma3_accc + 24.0 * v4rhosigma3_accb + 24.0 * v4rhosigma3_acbb + 8.0 * v4rhosigma3_abbb +
                              24.0 * v4rhosigma3_aacc + 48.0 * v4rhosigma3_aacb + 24.0 * v4rhosigma3_aabb + 24.0 * v4rhosigma3_aaac +
                              24.0 * v4rhosigma3_aaab + 8.0 * v4rhosigma3_aaaa;

                // Scalar contribution

                double prefac = 0.0;

                // vxc 1 contributions

                prefac += rr * rho3[g]  // l1
                          + rx * l2contract;

                // // vxc 2 contributions

                prefac += rrr * gam[g]  // q1
                          + rxr * q2contract + rxx * q3contract + rx * q4contract;

                // // // vxc 3 contributions
                prefac += rrrr * c1 + rxrr * c2 + rxxr * c3 + rxr * c4 + rxx * (c5_6 + c8) + rxxx * c7;

                G_val[nu_offset + g] = w * prefac * chi_val[nu_offset + g];

                // vector contribution

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                // vxc 1 contributions

                xcomp += xr * grada_x_g * rho3[g] + x * rxw123a + xx * l5contract_x;

                ycomp += xr * grada_y_g * rho3[g] + x * ryw123a + xx * l5contract_y;

                zcomp += xr * grada_z_g * rho3[g] + x * rzw123a + xx * l5contract_z;

                // // vxc 2 contributions

                xcomp += xrr * grada_x_g * gam[g] + xr * gamx[g] + xxr * q7contract_x + xx * (q8contract_x + q10contract_x + q11contract_x) +
                         xxx * q9contract_x;

                ycomp += xrr * grada_y_g * gam[g] + xr * gamy[g] + xxr * q7contract_y + xx * (q8contract_y + q10contract_y + q11contract_y) +
                         xxx * q9contract_y;

                zcomp += xrr * grada_z_g * gam[g] + xr * gamz[g] + xxr * q7contract_z + xx * (q8contract_z + q10contract_z + q11contract_z) +
                         xxx * q9contract_z;

                // vxc 3 contributions
                xcomp += xrrr * c9_x + xrr * c10_x + xxrr * c11_x + xxr * (c12_c14_x + c15_x) + xxxr * c13_x + xx * c17_24_25_x + xxxx * c18_x +
                         xxx * (c16_19_22_x + c20_21_23_x);

                ycomp += xrrr * c9_y + xrr * c10_y + xxrr * c11_y + xxr * (c12_c14_y + c15_y) + xxxr * c13_y + xx * c17_24_25_y + xxxx * c18_y +
                         xxx * (c16_19_22_y + c20_21_23_y);

                zcomp += xrrr * c9_z + xrr * c10_z + xxrr * c11_z + xxr * (c12_c14_z + c15_z) + xxxr * c13_z + xx * c17_24_25_z + xxxx * c18_z +
                         xxx * (c16_19_22_z + c20_21_23_z);

                G_gga_val[nu_offset + g] =
                    w * (xcomp * chi_x_val[nu_offset + g] + ycomp * chi_y_val[nu_offset + g] + zcomp * chi_z_val[nu_offset + g]);
            }
        }
    }

    timer.stop("Lxc matrix G");

    timer.start("Lxc matrix matmul");

    auto mat_Lxc = sdenblas::serialMultABt(gtoValues, mat_G);

    auto mat_Lxc_gga = sdenblas::serialMultABt(gtoValues, mat_G_gga);

    mat_Lxc_gga.symmetrize();  // matrix + matrix.T

    sdenblas::serialInPlaceAddAB(mat_Lxc, mat_Lxc_gga);

    timer.stop("Lxc matrix matmul");

    return mat_Lxc;
}

}  // namespace xcintgga
