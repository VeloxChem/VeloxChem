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

#include "XCMolecularGradientForLDA.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "AODensityMatrix.hpp"
#include "AOIndices.hpp"
#include "DensityGridQuad.hpp"
#include "DftSubMatrix.hpp"
#include "ErrorHandler.hpp"
#include "FunctionalParser.hpp"
#include "GridScreener.hpp"
#include "GtoFunc.hpp"
#include "GtoValues.hpp"
#include "Prescreener.hpp"
#include "SerialDenseLinearAlgebra.hpp"
#include "SerialDensityGridGenerator.hpp"
#include "XCFunctional.hpp"
#include "XCMolecularGradientForLDA.hpp"

namespace xcgradlda {  // xcgradlda namespace

auto
integrateVxcGradientForLdaClosedShell(const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const std::vector<const double*>& rwDensityPointers,
                                      const std::vector<const double*>& gsDensityPointers,
                                      const CMolecularGrid&   molecularGrid,
                                      const double            screeningThresholdForGTOValues,
                                      const CXCFunctional&    xcFunctional) -> CDenseMatrix
{
    CMultiTimer timer;

    timer.start("Total timing");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // AO-to-atom mapping

    std::vector<int> ao_to_atom_ids(naos);

    aoindices::computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.number_of_atoms();

    CDenseMatrix molgrad_threads(nthreads, natoms * 3);

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
    auto ptr_rwDensityPointers = rwDensityPointers.data();

    auto ptr_xcFunctional = &xcFunctional;

    auto ptr_molgrad_threads = &molgrad_threads;

#pragma omp parallel shared(ptr_counts, ptr_displacements, xcoords, ycoords, zcoords, \
                            ptr_gto_blocks, ptr_gsDensityPointers, ptr_rwDensityPointers, \
                            ptr_xcFunctional, ptr_molgrad_threads, \
                            n_boxes, n_gto_blocks, naos)
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

            auto gs_sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);
            auto rw_sub_dens_mat = dftsubmat::getSubDensityMatrix(rwDensityPointers[0], aoinds, naos);

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

            auto       ldafunc = local_xcfunc.getFunctionalPointerToLdaComponent();
            const auto dim     = &(ldafunc->dim);

            std::vector<double> local_weights_data(weights + gridblockpos, weights + gridblockpos + npoints);

            std::vector<double> rho_data(dim->rho * npoints);

            std::vector<double> vrho_data(dim->vrho * npoints);

            auto local_weights = local_weights_data.data();

            auto rho  = rho_data.data();

            auto vrho = vrho_data.data();

            sdengridgen::serialGenerateDensityForLDA(rho, mat_chi, gs_sub_dens_mat);

            omptimers[thread_id].stop("Generate density grid");

            // generate density gradient grid

            omptimers[thread_id].start("Density grad. grid prep.");

            CDenseMatrix dengradx(natoms, npoints);
            CDenseMatrix dengrady(natoms, npoints);
            CDenseMatrix dengradz(natoms, npoints);

            omptimers[thread_id].stop("Density grad. grid prep.");

            // eq.(26), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid matmul");

            auto mat_F = sdenblas::serialMultAB(rw_sub_dens_mat, mat_chi);

            omptimers[thread_id].stop("Density grad. grid matmul");

            // eq.(34), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid rho");

            auto F_val = mat_F.values();

            auto chi_x_val = mat_chi_x.values();
            auto chi_y_val = mat_chi_y.values();
            auto chi_z_val = mat_chi_z.values();

            auto gdenx = dengradx.values();
            auto gdeny = dengrady.values();
            auto gdenz = dengradz.values();

            for (int nu = 0; nu < aocount; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * npoints;

                auto nu_offset = nu * npoints;

                #pragma omp simd
                for (int g = 0; g < npoints; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];
                    gdeny[atom_g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];
                    gdenz[atom_g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];
                }
            }

            omptimers[thread_id].stop("Density grad. grid rho");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            local_xcfunc.compute_vxc_for_lda(npoints, rho, vrho);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Accumulate gradient");

            auto gatm = ptr_molgrad_threads->row(thread_id);

            for (int iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz)
                for (int g = 0; g < npoints; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac = local_weights[g] * vrho[2 * g + 0];

                    gatmx += prefac * gdenx[atom_g];
                    gatmy += prefac * gdeny[atom_g];
                    gatmz += prefac * gdenz[atom_g];
                }

                // factor of 2 from sum of alpha and beta contributions

                gatm[iatom * 3 + 0] += 2.0 * gatmx;
                gatm[iatom * 3 + 1] += 2.0 * gatmy;
                gatm[iatom * 3 + 2] += 2.0 * gatmz;
            }

            omptimers[thread_id].stop("Accumulate gradient");
        }
    }
    }
    }
    }

    timer.stop("Total timing");

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;

    CDenseMatrix molgrad(natoms, 3);

    for (int iatom = 0; iatom < natoms; iatom++)
    {
        for (int thread_id = 0; thread_id < nthreads; thread_id++)
        {
            molgrad.row(iatom)[0] += molgrad_threads.row(thread_id)[iatom * 3 + 0];
            molgrad.row(iatom)[1] += molgrad_threads.row(thread_id)[iatom * 3 + 1];
            molgrad.row(iatom)[2] += molgrad_threads.row(thread_id)[iatom * 3 + 2];
        }
    }

    return molgrad;
}

auto
integrateVxcGradientForLdaOpenShell(const CMolecule&        molecule,
                                    const CMolecularBasis&  basis,
                                    const std::vector<const double*>& rwDensityPointers,
                                    const std::vector<const double*>& gsDensityPointers,
                                    const CMolecularGrid&   molecularGrid,
                                    const double            screeningThresholdForGTOValues,
                                    const CXCFunctional&    xcFunctional) -> CDenseMatrix
{
    CMultiTimer timer;

    timer.start("Total timing");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // AO-to-atom mapping

    std::vector<int> ao_to_atom_ids(naos);

    aoindices::computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.number_of_atoms();

    CDenseMatrix molgrad_threads(nthreads, natoms * 3);

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
    auto ptr_rwDensityPointers = rwDensityPointers.data();

    auto ptr_xcFunctional = &xcFunctional;

    auto ptr_molgrad_threads = &molgrad_threads;

#pragma omp parallel shared(ptr_counts, ptr_displacements, xcoords, ycoords, zcoords, \
                            ptr_gto_blocks, ptr_gsDensityPointers, ptr_rwDensityPointers, \
                            ptr_xcFunctional, ptr_molgrad_threads, \
                            n_boxes, n_gto_blocks, naos)
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

            auto gs_sub_dens_mat_a = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);
            auto gs_sub_dens_mat_b = dftsubmat::getSubDensityMatrix(gsDensityPointers[1], aoinds, naos);

            auto rw_sub_dens_mat_a = dftsubmat::getSubDensityMatrix(rwDensityPointers[0], aoinds, naos);
            auto rw_sub_dens_mat_b = dftsubmat::getSubDensityMatrix(rwDensityPointers[1], aoinds, naos);

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

            // generate density grid

            omptimers[thread_id].start("Generate density grid");

            auto local_xcfunc = CXCFunctional(*ptr_xcFunctional);

            auto       ldafunc = local_xcfunc.getFunctionalPointerToLdaComponent();
            const auto dim     = &(ldafunc->dim);

            std::vector<double> local_weights_data(weights + gridblockpos, weights + gridblockpos + npoints);

            std::vector<double> rho_data(dim->rho * npoints);

            std::vector<double> vrho_data(dim->vrho * npoints);

            auto local_weights = local_weights_data.data();

            auto rho  = rho_data.data();

            auto vrho = vrho_data.data();

            sdengridgen::serialGenerateDensityForLDA(rho, mat_chi, gs_sub_dens_mat_a, gs_sub_dens_mat_b);

            omptimers[thread_id].stop("Generate density grid");

            // generate density gradient grid

            omptimers[thread_id].start("Density grad. grid prep.");

            CDenseMatrix dengrad_a_x(natoms, npoints);
            CDenseMatrix dengrad_a_y(natoms, npoints);
            CDenseMatrix dengrad_a_z(natoms, npoints);

            CDenseMatrix dengrad_b_x(natoms, npoints);
            CDenseMatrix dengrad_b_y(natoms, npoints);
            CDenseMatrix dengrad_b_z(natoms, npoints);

            omptimers[thread_id].stop("Density grad. grid prep.");

            // eq.(26), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid matmul");

            auto mat_F_a = sdenblas::serialMultAB(rw_sub_dens_mat_a, mat_chi);
            auto mat_F_b = sdenblas::serialMultAB(rw_sub_dens_mat_b, mat_chi);

            omptimers[thread_id].stop("Density grad. grid matmul");

            // eq.(34), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid rho");

            auto F_a_val = mat_F_a.values();
            auto F_b_val = mat_F_b.values();

            auto chi_x_val = mat_chi_x.values();
            auto chi_y_val = mat_chi_y.values();
            auto chi_z_val = mat_chi_z.values();

            auto gden_a_x = dengrad_a_x.values();
            auto gden_a_y = dengrad_a_y.values();
            auto gden_a_z = dengrad_a_z.values();

            auto gden_b_x = dengrad_b_x.values();
            auto gden_b_y = dengrad_b_y.values();
            auto gden_b_z = dengrad_b_z.values();

            for (int nu = 0; nu < aocount; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * npoints;

                auto nu_offset = nu * npoints;

                #pragma omp simd
                for (int g = 0; g < npoints; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gden_a_x[atom_g] -= 2.0 * F_a_val[nu_g] * chi_x_val[nu_g];
                    gden_a_y[atom_g] -= 2.0 * F_a_val[nu_g] * chi_y_val[nu_g];
                    gden_a_z[atom_g] -= 2.0 * F_a_val[nu_g] * chi_z_val[nu_g];

                    gden_b_x[atom_g] -= 2.0 * F_b_val[nu_g] * chi_x_val[nu_g];
                    gden_b_y[atom_g] -= 2.0 * F_b_val[nu_g] * chi_y_val[nu_g];
                    gden_b_z[atom_g] -= 2.0 * F_b_val[nu_g] * chi_z_val[nu_g];
                }
            }

            omptimers[thread_id].stop("Density grad. grid rho");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            local_xcfunc.compute_vxc_for_lda(npoints, rho, vrho);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Accumulate gradient");

            auto gatm = ptr_molgrad_threads->row(thread_id);

            for (int iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz)
                for (int g = 0; g < npoints; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac_a = local_weights[g] * vrho[2 * g + 0];
                    double prefac_b = local_weights[g] * vrho[2 * g + 1];

                    gatmx += prefac_a * gden_a_x[atom_g];
                    gatmy += prefac_a * gden_a_y[atom_g];
                    gatmz += prefac_a * gden_a_z[atom_g];

                    gatmx += prefac_b * gden_b_x[atom_g];
                    gatmy += prefac_b * gden_b_y[atom_g];
                    gatmz += prefac_b * gden_b_z[atom_g];
                }

                gatm[iatom * 3 + 0] += gatmx;
                gatm[iatom * 3 + 1] += gatmy;
                gatm[iatom * 3 + 2] += gatmz;
            }

            omptimers[thread_id].stop("Accumulate gradient");
        }
    }
    }
    }
    }

    timer.stop("Total timing");

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;

    CDenseMatrix molgrad(natoms, 3);

    for (int iatom = 0; iatom < natoms; iatom++)
    {
        for (int thread_id = 0; thread_id < nthreads; thread_id++)
        {
            molgrad.row(iatom)[0] += molgrad_threads.row(thread_id)[iatom * 3 + 0];
            molgrad.row(iatom)[1] += molgrad_threads.row(thread_id)[iatom * 3 + 1];
            molgrad.row(iatom)[2] += molgrad_threads.row(thread_id)[iatom * 3 + 2];
        }
    }

    return molgrad;
}

auto
integrateFxcGradientForLdaClosedShell(const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const std::vector<const double*>& rwDensityPointersOne,
                                      const std::vector<const double*>& rwDensityPointersTwo,
                                      const std::vector<const double*>& gsDensityPointers,
                                      const CMolecularGrid&   molecularGrid,
                                      const double            screeningThresholdForGTOValues,
                                      const CXCFunctional&    xcFunctional) -> CDenseMatrix
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // AO-to-atom mapping

    std::vector<int> ao_to_atom_ids(naos);

    aoindices::computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.number_of_atoms();

    CDenseMatrix molgrad_threads(nthreads, natoms * 3);

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

        // generate sub density matrix and density grid

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

        auto rw_sub_dens_mat_one = dftsubmat::getSubDensityMatrix(rwDensityPointersOne[0], aoinds, naos);
        auto rw_sub_dens_mat_two = dftsubmat::getSubDensityMatrix(rwDensityPointersTwo[0], aoinds, naos);

        timer.stop("Density matrix slicing");

        // GTO values on grid points

        timer.start("OMP FxcGrad calc.");

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

            auto rho    = omp_rho_data[thread_id].data();
            auto rhow   = omp_rhow_data[thread_id].data();
            auto v2rho2 = omp_v2rho2_data[thread_id].data();

            sdengridgen::serialGenerateDensityForLDA(rho, mat_chi, gs_sub_dens_mat);

            sdengridgen::serialGenerateDensityForLDA(rhow, mat_chi, rw_sub_dens_mat_one);

            omptimers[thread_id].stop("Generate density grid");

            // generate density gradient grid

            omptimers[thread_id].start("Density grad. grid prep.");

            CDenseMatrix dengradx(natoms, grid_batch_size);
            CDenseMatrix dengrady(natoms, grid_batch_size);
            CDenseMatrix dengradz(natoms, grid_batch_size);

            omptimers[thread_id].stop("Density grad. grid prep.");

            // eq.(26), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid matmul");

            auto mat_F = sdenblas::serialMultAB(rw_sub_dens_mat_two, mat_chi);

            omptimers[thread_id].stop("Density grad. grid matmul");

            // eq.(34), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid rho");

            auto F_val = mat_F.values();

            auto chi_x_val = mat_chi_x.values();
            auto chi_y_val = mat_chi_y.values();
            auto chi_z_val = mat_chi_z.values();

            auto gdenx = dengradx.values();
            auto gdeny = dengrady.values();
            auto gdenz = dengradz.values();

            for (int nu = 0; nu < aocount; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * grid_batch_size;

                auto nu_offset = nu * grid_batch_size;

                #pragma omp simd
                for (int g = 0; g < grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];
                    gdeny[atom_g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];
                    gdenz[atom_g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];
                }
            }

            omptimers[thread_id].stop("Density grad. grid rho");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_fxc_for_lda(grid_batch_size, rho, v2rho2);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            // eq.(32), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Accumulate gradient");

            auto gatm = molgrad_threads.row(thread_id);

            for (int iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * grid_batch_size;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz)
                for (int g = 0; g < grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac = local_weights[g] * (v2rho2[3 * g + 0] * rhow[2 * g + 0] + v2rho2[3 * g + 1] * rhow[2 * g + 1]);

                    gatmx += prefac * gdenx[atom_g];
                    gatmy += prefac * gdeny[atom_g];
                    gatmz += prefac * gdenz[atom_g];
                }

                // factor of 2 from sum of alpha and beta contributions

                gatm[iatom * 3 + 0] += 2.0 * gatmx;
                gatm[iatom * 3 + 1] += 2.0 * gatmy;
                gatm[iatom * 3 + 2] += 2.0 * gatmz;
            }

            omptimers[thread_id].stop("Accumulate gradient");
        }

        timer.stop("OMP FxcGrad calc.");
    }

    timer.stop("Total timing");

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;

    CDenseMatrix molgrad(natoms, 3);

    for (int iatom = 0; iatom < natoms; iatom++)
    {
        for (int thread_id = 0; thread_id < nthreads; thread_id++)
        {
            molgrad.row(iatom)[0] += molgrad_threads.row(thread_id)[iatom * 3 + 0];
            molgrad.row(iatom)[1] += molgrad_threads.row(thread_id)[iatom * 3 + 1];
            molgrad.row(iatom)[2] += molgrad_threads.row(thread_id)[iatom * 3 + 2];
        }
    }

    return molgrad;
}

auto
integrateKxcGradientForLdaClosedShell(const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const std::vector<const double*>& rwDensityPointersOne,
                                      const std::vector<const double*>& rwDensityPointersTwo,
                                      const std::vector<const double*>& gsDensityPointers,
                                      const CMolecularGrid&   molecularGrid,
                                      const double            screeningThresholdForGTOValues,
                                      const CXCFunctional&    xcFunctional) -> CDenseMatrix
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // AO-to-atom mapping

    std::vector<int> ao_to_atom_ids(naos);

    aoindices::computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.number_of_atoms();

    CDenseMatrix molgrad_threads(nthreads, natoms * 3);

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

        // generate sub density matrix and density grid

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

        auto rw_sub_dens_mat_one = dftsubmat::getSubDensityMatrix(rwDensityPointersOne[0], aoinds, naos);
        auto rw_sub_dens_mat_two = dftsubmat::getSubDensityMatrix(rwDensityPointersTwo[0], aoinds, naos);

        timer.stop("Density matrix slicing");

        // GTO values on grid points

        timer.start("OMP KxcGrad calc.");

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

            auto rho = omp_rho_data[thread_id].data();

            auto v2rho2 = omp_v2rho2_data[thread_id].data();
            auto v3rho3 = omp_v3rho3_data[thread_id].data();

            sdengridgen::serialGenerateDensityForLDA(rho, mat_chi, gs_sub_dens_mat);

            // compute perturbed density

            // prepare rwdenmat for quadratic response

            omptimers[thread_id].start("Density grid quad");

            CDenseMatrix zero_sub_den_mat_one(rw_sub_dens_mat_one);
            CDenseMatrix zero_sub_den_mat_two(rw_sub_dens_mat_two);

            zero_sub_den_mat_one.zero();
            zero_sub_den_mat_two.zero();

            CAODensityMatrix rwdenmat(std::vector<CDenseMatrix>({rw_sub_dens_mat_one, zero_sub_den_mat_one, rw_sub_dens_mat_two, zero_sub_den_mat_two}),
                                      denmat::rest);

            // Note: We use quadratic response (quadMode == "QRF") to calculate
            // third-order functional derivative contribution. The rw2DensityMatrix
            // contains zero matrices and is therefore removed from the following code.
            // Same for rw2dengrid.

            // For "QRF" we have rwDensityMatrix.getNumberOfDensityMatrices() ==
            // 2 * rw2DensityMatrix.getNumberOfDensityMatrices()

            std::string quadMode("QRF");

            auto numdens_rw2 = rwdenmat.getNumberOfDensityMatrices() / 2;

            auto xcfuntype = omp_xcfuncs[thread_id].getFunctionalType();

            auto rwdengrid = sdengridgen::serialGenerateDensityGridForLDA(mat_chi, rwdenmat, xcfuntype);

            CDensityGridQuad rwdengridquad(grid_batch_size, numdens_rw2, xcfuntype, dengrid::ab);

            rwdengridquad.DensityProd(rwdengrid, xcfuntype, numdens_rw2, quadMode);

            omptimers[thread_id].stop("Density grid quad");

            // generate density gradient grid

            omptimers[thread_id].start("Density grad. grid prep.");

            CDenseMatrix dengradx(natoms, grid_batch_size);
            CDenseMatrix dengrady(natoms, grid_batch_size);
            CDenseMatrix dengradz(natoms, grid_batch_size);

            omptimers[thread_id].stop("Density grad. grid prep.");

            // eq.(26), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid matmul");

            auto mat_F = sdenblas::serialMultAB(gs_sub_dens_mat, mat_chi);

            omptimers[thread_id].stop("Density grad. grid matmul");

            // eq.(34), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid rho");

            auto F_val = mat_F.values();

            auto chi_x_val = mat_chi_x.values();
            auto chi_y_val = mat_chi_y.values();
            auto chi_z_val = mat_chi_z.values();

            auto gdenx = dengradx.values();
            auto gdeny = dengrady.values();
            auto gdenz = dengradz.values();

            for (int nu = 0; nu < aocount; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * grid_batch_size;

                auto nu_offset = nu * grid_batch_size;

                #pragma omp simd
                for (int g = 0; g < grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];
                    gdeny[atom_g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];
                    gdenz[atom_g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];
                }
            }

            omptimers[thread_id].stop("Density grad. grid rho");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_fxc_for_lda(grid_batch_size, rho, v2rho2);

            omp_xcfuncs[thread_id].compute_kxc_for_lda(grid_batch_size, rho, v3rho3);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            // pointers to perturbed density gradient norms

            auto rhow1a = rwdengridquad.gam(0);

            // Note: rw2DensityMatrix is zero in KxcGradientForLDA
            // auto rhow12a = rw2DensityGrid.alphaDensity(iFock);
            // auto rhow12b = rw2DensityGrid.betaDensity(iFock);

            // eq.(32), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Accumulate gradient");

            auto gatm = molgrad_threads.row(thread_id);

            for (int iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * grid_batch_size;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz)
                for (int g = 0; g < grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac = local_weights[g] * (v3rho3[4 * g + 0] + 2.0 * v3rho3[4 * g + 1] + v3rho3[4 * g + 2]) * rhow1a[g];

                    gatmx += prefac * gdenx[atom_g];
                    gatmy += prefac * gdeny[atom_g];
                    gatmz += prefac * gdenz[atom_g];
                }

                // factor of 2 from sum of alpha and beta contributions
                // factor of 0.25 from quadratic response

                gatm[iatom * 3 + 0] += 0.25 * (2.0 * gatmx);
                gatm[iatom * 3 + 1] += 0.25 * (2.0 * gatmy);
                gatm[iatom * 3 + 2] += 0.25 * (2.0 * gatmz);
            }

            omptimers[thread_id].stop("Accumulate gradient");
        }

        timer.stop("OMP KxcGrad calc.");
    }

    timer.stop("Total timing");

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;

    CDenseMatrix molgrad(natoms, 3);

    for (int iatom = 0; iatom < natoms; iatom++)
    {
        for (int thread_id = 0; thread_id < nthreads; thread_id++)
        {
            molgrad.row(iatom)[0] += molgrad_threads.row(thread_id)[iatom * 3 + 0];
            molgrad.row(iatom)[1] += molgrad_threads.row(thread_id)[iatom * 3 + 1];
            molgrad.row(iatom)[2] += molgrad_threads.row(thread_id)[iatom * 3 + 2];
        }
    }

    return molgrad;
}

}  // namespace xcgradlda
