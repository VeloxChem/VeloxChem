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

#include "XCMolecularGradientForPDFT.hpp"

#include <omp.h>
#include <cstring>

#include "AOIndices.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DftSubMatrix.hpp"
#include "GtoFunc.hpp"
#include "GtoValues.hpp"
#include "MultiTimer.hpp"
#include "PairDensityGridGenerator.hpp"
#include "Prescreener.hpp"
#include "SerialDenseLinearAlgebra.hpp"

namespace xcgradpdft {  // xcgradpdft namespace

CDenseMatrix
integrateVxcPDFTGradientForLDA(const CMolecule&                molecule,
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

    auto ptr_twoBodyDensityMatrix = &twoBodyDensityMatrix;
    auto ptr_activeMOs = &activeMOs;

    auto ptr_xcFunctional = &xcFunctional;

#pragma omp parallel shared(ptr_counts, ptr_displacements, xcoords, ycoords, zcoords, \
                            ptr_gto_blocks, densityMatrixPointer, \
                            ptr_twoBodyDensityMatrix, ptr_activeMOs, ptr_xcFunctional, \
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

            std::vector<double> exc_data(1 * npoints); //Not needed but always provided for now
            std::vector<double> vrho_data(2 * npoints);

            auto local_weights = local_weights_data.data();

            auto rho  = rho_data.data();

            auto exc  = exc_data.data();
            auto vrho = vrho_data.data();

            // generate density and on-top pair density on the grid

            pairdengridgen::serialGeneratePairDensityForLDA(rho, mat_chi, sub_dens_mat_a, sub_active_mos, twoBodyDensityMatrix);

            // generate density gradient grid

            omptimers[thread_id].start("Density grad. grid prep.");

            CDenseMatrix dengradx(natoms, npoints);
            CDenseMatrix dengrady(natoms, npoints);
            CDenseMatrix dengradz(natoms, npoints);

            omptimers[thread_id].stop("Density grad. grid prep.");

            // eq.(26), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid matmul");

            auto mat_F = sdenblas::serialMultAB(sub_dens_mat_a, mat_chi);

            omptimers[thread_id].stop("Density grad. grid matmul");

            // Pair-density parts (this recomputes a lot of things)

            omptimers[thread_id].start("Density grad mo pair");

            auto n_active = activeMOs.getNumberOfRows();

            //1) \phi_t(r) = C_mu^t \phi_\mu(r)
            CDenseMatrix mos_on_grid;

            if (n_active > 0)
            {
                mos_on_grid = sdenblas::serialMultAB(sub_active_mos, mat_chi);
            }

            auto n_active2 = n_active * n_active;

            //2) \phi_tu(r) = \phi_t(r) \phi_u(r)
            CDenseMatrix mo_pair(n_active2, npoints);

            auto mo_pair_val = mo_pair.values();

            for (int t = 0; t < n_active; t++)
            {
                auto MOt = mos_on_grid.row(t);

                auto t_offset = t * n_active * npoints;

                for (int u = 0; u < n_active; u++)
                {
                    auto MOu = mos_on_grid.row(u);

                    auto tu_offset = t_offset + u * npoints;

                    #pragma omp simd 
                    for (int g = 0; g < npoints; g++)
                    {
                        mo_pair_val[tu_offset + g] += MOt[g] * MOu[g];
                    }
                }
            }

            omptimers[thread_id].stop("Density grad mo pair");

            omptimers[thread_id].start("Density grad pi matmul");

            //3) d_tu(r) = d_tuvw \phi_v(r) \phi_w(r)
            auto mat_d = sdenblas::serialMultAB(twoBodyDensityMatrix, mo_pair);

            omptimers[thread_id].stop("Density grad pi matmul");

            //4) g_t(r) = d_tu(r) phi_u(r)
            CDenseMatrix mat_g(n_active, npoints);

            auto g_val = mat_g.values();

            auto d_val = mat_d.values();

            for (int v = 0; v < n_active; v++)
            {
                auto v_offset = v * npoints;

                for (int w = 0; w < n_active; w++)
                {
                    auto vw_offset = (v*n_active + w) * npoints;

                    auto MOw = mos_on_grid.row(w);

                    #pragma omp simd 
                    for (int g = 0; g < npoints; g++)
                    {
                        g_val[v_offset + g] += d_val[vw_offset + g] * MOw[g];
                    }
                }

            }

            //5) k_mu(r) = C_mu^t g_t(r)
            CDenseMatrix mat_k(aocount, npoints);

            if (n_active > 0)
            {
                mat_k = sdenblas::serialMultAtB(sub_active_mos, mat_g);
            }
            else
            {
                mat_k.zero();
            }

            omptimers[thread_id].start("Density grad. grid rho");

            auto F_val = mat_F.values();
            auto k_val = mat_k.values();

            auto chi_x_val = mat_chi_x.values();
            auto chi_y_val = mat_chi_y.values();
            auto chi_z_val = mat_chi_z.values();

            auto gdenx = dengradx.values();
            auto gdeny = dengrady.values();
            auto gdenz = dengradz.values();

            CDenseMatrix dengradpi_x(natoms, npoints);
            CDenseMatrix dengradpi_y(natoms, npoints);
            CDenseMatrix dengradpi_z(natoms, npoints);

            auto gdenpi_x = dengradpi_x.values();
            auto gdenpi_y = dengradpi_y.values();
            auto gdenpi_z = dengradpi_z.values();

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

                    gdenpi_x[atom_g] -= 4.0 * k_val[nu_g] * chi_x_val[nu_g];
                    gdenpi_y[atom_g] -= 4.0 * k_val[nu_g] * chi_y_val[nu_g];
                    gdenpi_z[atom_g] -= 4.0 * k_val[nu_g] * chi_z_val[nu_g];
                }
            }

            omptimers[thread_id].stop("Density grad. grid rho");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            local_xcfunc.compute_exc_vxc_for_plda(npoints, rho, exc, vrho, rs_omega);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Accumulate gradient");

            auto gatm = molgrad_threads.row(thread_id);
            
            for (int iatom = 0; iatom < natoms; iatom++)
            {   
                auto atom_offset = iatom * npoints;
                
                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;
                
                #pragma omp simd
                for (int g = 0; g < npoints; g++)
                {   
                    auto atom_g = atom_offset + g;
                    
                    double prefac = local_weights[g] * vrho[2 * g + 0];
                    double prefacpi = local_weights[g] * vrho[2 * g + 1];
                    
                    gatmx += prefac * gdenx[atom_g] + prefacpi * gdenpi_x[atom_g];
                    gatmy += prefac * gdeny[atom_g] + prefacpi * gdenpi_y[atom_g];
                    gatmz += prefac * gdenz[atom_g] + prefacpi * gdenpi_z[atom_g];
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

CDenseMatrix
integrateVxcPDFTGradientForGGA(const CMolecule&                molecule,
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

    // density and functional derivatives

    std::vector<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    std::vector<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    std::vector<double> rhograd_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    std::vector<double> sigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    std::vector<double> exc_data(1 * max_npoints_per_box); //Not needed but always provided for now
    std::vector<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    std::vector<double> vsigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();

    auto exc  = exc_data.data();
    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();

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
            // 2nd order GTO derivative
            auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(gto_block, 2, screeningThresholdForGTOValues, boxdim);

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

        CDenseMatrix mat_chi_xx(aocount, npoints);
        CDenseMatrix mat_chi_xy(aocount, npoints);
        CDenseMatrix mat_chi_xz(aocount, npoints);
        CDenseMatrix mat_chi_yy(aocount, npoints);
        CDenseMatrix mat_chi_yz(aocount, npoints);
        CDenseMatrix mat_chi_zz(aocount, npoints);

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

                auto cmat = gtoval::get_gto_values_for_mgga(gto_block, grid_x, grid_y, grid_z, cgto_mask);

                if (cmat.is_empty()) continue;

                auto submat_0_ptr = cmat.sub_matrix({0, 0});

                auto submat_x_ptr = cmat.sub_matrix({1, 0});
                auto submat_y_ptr = cmat.sub_matrix({1, 1});
                auto submat_z_ptr = cmat.sub_matrix({1, 2});

                auto submat_xx_ptr = cmat.sub_matrix({2, 0});
                auto submat_xy_ptr = cmat.sub_matrix({2, 1});
                auto submat_xz_ptr = cmat.sub_matrix({2, 2});
                auto submat_yy_ptr = cmat.sub_matrix({2, 3});
                auto submat_yz_ptr = cmat.sub_matrix({2, 4});
                auto submat_zz_ptr = cmat.sub_matrix({2, 5});

                auto submat_0_data = submat_0_ptr->data();

                auto submat_x_data = submat_x_ptr->data();
                auto submat_y_data = submat_y_ptr->data();
                auto submat_z_data = submat_z_ptr->data();

                auto submat_xx_data = submat_xx_ptr->data();
                auto submat_xy_data = submat_xy_ptr->data();
                auto submat_xz_data = submat_xz_ptr->data();
                auto submat_yy_data = submat_yy_ptr->data();
                auto submat_yz_data = submat_yz_ptr->data();
                auto submat_zz_data = submat_zz_ptr->data();

                for (int nu = 0; nu < static_cast<int>(pre_ao_inds.size()); nu++, idx++)
                {
                    std::memcpy(mat_chi.row(idx) + grid_batch_offset, submat_0_data + nu * grid_batch_size, grid_batch_size * sizeof(double));

                    std::memcpy(mat_chi_x.row(idx) + grid_batch_offset, submat_x_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_y.row(idx) + grid_batch_offset, submat_y_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_z.row(idx) + grid_batch_offset, submat_z_data + nu * grid_batch_size, grid_batch_size * sizeof(double));

                    std::memcpy(mat_chi_xx.row(idx) + grid_batch_offset, submat_xx_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xy.row(idx) + grid_batch_offset, submat_xy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xz.row(idx) + grid_batch_offset, submat_xz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yy.row(idx) + grid_batch_offset, submat_yy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yz.row(idx) + grid_batch_offset, submat_yz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_zz.row(idx) + grid_batch_offset, submat_zz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");
        }

        timer.stop("OMP GTO evaluation");

        // generate sub density matrix and density grid

        timer.start("Density matrix slicing");

        auto sub_dens_mat_a = dftsubmat::getSubDensityMatrix(densityMatrixPointer, aoinds, naos);

        auto sub_active_mos = dftsubmat::getSubMatrixByColumnSlicing(activeMOs, aoinds, naos);

        timer.stop("Density matrix slicing");

        // generate density and on-top pair density on the grid

        pairdengridgen::serialGeneratePairDensityForGGA(
            rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat_a, sub_active_mos, twoBodyDensityMatrix);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, npoints);
        CDenseMatrix dengrady(natoms, npoints);
        CDenseMatrix dengradz(natoms, npoints);

        CDenseMatrix dengradxx(natoms, npoints);
        CDenseMatrix dengradxy(natoms, npoints);
        CDenseMatrix dengradxz(natoms, npoints);

        CDenseMatrix dengradyx(natoms, npoints);
        CDenseMatrix dengradyy(natoms, npoints);
        CDenseMatrix dengradyz(natoms, npoints);

        CDenseMatrix dengradzx(natoms, npoints);
        CDenseMatrix dengradzy(natoms, npoints);
        CDenseMatrix dengradzz(natoms, npoints);

        timer.stop("Density grad. grid prep.");

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(sub_dens_mat_a, mat_chi);

        auto mat_F_x = denblas::multAB(sub_dens_mat_a, mat_chi_x);
        auto mat_F_y = denblas::multAB(sub_dens_mat_a, mat_chi_y);
        auto mat_F_z = denblas::multAB(sub_dens_mat_a, mat_chi_z);

        timer.stop("Density grad. grid matmul");

        //PDFT parts (this recomputes a lot of things)
        auto n_active = activeMOs.getNumberOfRows();

        //1) \phi_t(r) = C_mu^t \phi_\mu(r)
        CDenseMatrix mos_on_grid;
        //CDenseMatrix mos_x_on_grid; // For now we do not have dependence on pi'
        //CDenseMatrix mos_y_on_grid;
        //CDenseMatrix mos_z_on_grid;
        if (n_active > 0)
        {
            mos_on_grid = denblas::multAB(sub_active_mos, mat_chi);
            //mos_x_on_grid = denblas::multAB(sub_active_mos, mat_chi_x);
            //mos_y_on_grid = denblas::multAB(sub_active_mos, mat_chi_y);
            //mos_z_on_grid = denblas::multAB(sub_active_mos, mat_chi_z);
        }

        timer.start("Density grad mo pair");

        auto n_active2 = n_active * n_active;

        //2) \phi_tu(r) = \phi_t(r) \phi_u(r)
        CDenseMatrix mo_pair(n_active2, npoints);

        auto mo_pair_val = mo_pair.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            for (int t = 0; t < n_active; t++)
            {
                auto MOt = mos_on_grid.row(t);

                auto t_offset = t * n_active * npoints;

                for (int u = 0; u < n_active; u++)
                {
                    auto MOu = mos_on_grid.row(u);

                    auto tu_offset = t_offset + u * npoints;

                    #pragma omp simd 
                    for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                    {
                        mo_pair_val[tu_offset + g] += MOt[g] * MOu[g];
                    }
                }
            }
        }

        timer.stop("Density grad mo pair");

        timer.start("Density grad pi matmul");

        //3) d_tu(r) = d_tuvw \phi_v(r) \phi_w(r)
        auto mat_d = denblas::multAB(twoBodyDensityMatrix, mo_pair);

        timer.stop("Density grad pi matmul");

        //4) g_t(r) = d_tu(r) phi_u(r)
        CDenseMatrix mat_g(n_active, npoints);
        //CDenseMatrix mat_g_x(n_active, npoints);
        //CDenseMatrix mat_g_y(n_active, npoints);
        //CDenseMatrix mat_g_z(n_active, npoints);

        auto g_val = mat_g.values();
        //auto g_x_val = mat_g_x.values();
        //auto g_y_val = mat_g_y.values();
        //auto g_x_val = mat_g_z.values();

        auto d_val = mat_d.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            for (int v = 0; v < n_active; v++)
            {
                auto v_offset = v * npoints;

                for (int w = 0; w < n_active; w++)
                {
                    auto vw_offset = (v*n_active + w) * npoints;

                    auto MOw = mos_on_grid.row(w);
                    //auto MOw_x = mos_x_on_grid.row(w);
                    //auto MOw_y = mos_y_on_grid.row(w);
                    //auto MOw_z = mos_z_on_grid.row(w);

                    #pragma omp simd 
                    for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                    {
                        g_val[v_offset + g] += d_val[vw_offset + g] * MOw[g];
                        //g_x_val[v_offset + g] += d_val[vw_offset + g] * MOw_x[g];
                        //g_y_val[v_offset + g] += d_val[vw_offset + g] * MOw_y[g];
                        //g_z_val[v_offset + g] += d_val[vw_offset + g] * MOw_z[g];
                    }
                }

            }
        }

        //5) k_mu(r) = C_mu^t g_t(r)
        CDenseMatrix mat_k(aocount, npoints);
        //CDenseMatrix mat_k_x(aocount, npoints);
        //CDenseMatrix mat_k_y(aocount, npoints);
        //CDenseMatrix mat_k_z(aocount, npoints);
        if (n_active > 0)
        {
            mat_k = denblas::multAtB(sub_active_mos, mat_g);
            //mat_k_x = denblas::multAtB(sub_active_mos, mat_g_x);
            //mat_k_y = denblas::multAtB(sub_active_mos, mat_g_y);
            //mat_k_z = denblas::multAtB(sub_active_mos, mat_g_z);
        }
        else
        {
            mat_k.zero();
            //mat_k_x.zero();
            //mat_k_y.zero();
            //mat_k_z.zero();
        }

        timer.start("Density grad. grid rho");

        auto F_val = mat_F.values();
        auto F_x_val = mat_F_x.values();
        auto F_y_val = mat_F_y.values();
        auto F_z_val = mat_F_z.values();

        auto k_val = mat_k.values();
        //auto k_x_val = mat_k_x.values();
        //auto k_y_val = mat_k_y.values();
        //auto k_z_val = mat_k_z.values();

        auto chi_x_val = mat_chi_x.values();
        auto chi_y_val = mat_chi_y.values();
        auto chi_z_val = mat_chi_z.values();

        auto chi_xx_val = mat_chi_xx.values();
        auto chi_xy_val = mat_chi_xy.values();
        auto chi_xz_val = mat_chi_xz.values();
        auto chi_yy_val = mat_chi_yy.values();
        auto chi_yz_val = mat_chi_yz.values();
        auto chi_zz_val = mat_chi_zz.values();

        auto gdenx = dengradx.values();
        auto gdeny = dengrady.values();
        auto gdenz = dengradz.values();

        auto gdenxx = dengradxx.values();
        auto gdenxy = dengradxy.values();
        auto gdenxz = dengradxz.values();

        auto gdenyx = dengradyx.values();
        auto gdenyy = dengradyy.values();
        auto gdenyz = dengradyz.values();

        auto gdenzx = dengradzx.values();
        auto gdenzy = dengradzy.values();
        auto gdenzz = dengradzz.values();

        CDenseMatrix dengradpi_x(natoms, npoints);
        CDenseMatrix dengradpi_y(natoms, npoints);
        CDenseMatrix dengradpi_z(natoms, npoints);
        auto gdenpi_x = dengradpi_x.values();
        auto gdenpi_y = dengradpi_y.values();
        auto gdenpi_z = dengradpi_z.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            for (int nu = 0; nu < aocount; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * npoints;

                auto nu_offset = nu * npoints;

                #pragma omp simd 
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];
                    gdeny[atom_g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];
                    gdenz[atom_g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];

                    gdenxx[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xx_val[nu_g]);
                    gdenxy[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_xy_val[nu_g]);
                    gdenxz[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_xz_val[nu_g]);

                    gdenyx[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xy_val[nu_g]);
                    gdenyy[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_yy_val[nu_g]);
                    gdenyz[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_yz_val[nu_g]);

                    gdenzx[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xz_val[nu_g]);
                    gdenzy[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_yz_val[nu_g]);
                    gdenzz[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_zz_val[nu_g]);

                    gdenpi_x[atom_g] -= 4.0 * k_val[nu_g] * chi_x_val[nu_g];
                    gdenpi_y[atom_g] -= 4.0 * k_val[nu_g] * chi_y_val[nu_g];
                    gdenpi_z[atom_g] -= 4.0 * k_val[nu_g] * chi_z_val[nu_g];

/*                    gdenpi_xx[atom_g] -= 3.0 * k_x_val[nu_g] * chi_x_val[nu_g] + k_val[nu_g] * chi_xx_val[nu_g];
                    gdenpi_xy[atom_g] -= 3.0 * k_x_val[nu_g] * chi_y_val[nu_g] + k_val[nu_g] * chi_xy_val[nu_g];
                    gdenpi_xz[atom_g] -= 3.0 * k_x_val[nu_g] * chi_z_val[nu_g] + k_val[nu_g] * chi_xz_val[nu_g];

                    gdenpi_yx[atom_g] -= 3.0 * k_y_val[nu_g] * chi_x_val[nu_g] + k_val[nu_g] * chi_xy_val[nu_g];
                    gdenpi_yy[atom_g] -= 3.0 * k_y_val[nu_g] * chi_y_val[nu_g] + k_val[nu_g] * chi_yy_val[nu_g];
                    gdenpi_yz[atom_g] -= 3.0 * k_y_val[nu_g] * chi_z_val[nu_g] + k_val[nu_g] * chi_yz_val[nu_g];

                    gdenpi_zx[atom_g] -= 3.0 * k_z_val[nu_g] * chi_x_val[nu_g] + k_val[nu_g] * chi_xz_val[nu_g];
                    gdenpi_zy[atom_g] -= 3.0 * k_z_val[nu_g] * chi_y_val[nu_g] + k_val[nu_g] * chi_yz_val[nu_g];
                    gdenpi_zz[atom_g] -= 3.0 * k_z_val[nu_g] * chi_z_val[nu_g] + k_val[nu_g] * chi_zz_val[nu_g];*/
                }
            }
        }

        timer.stop("Density grad. grid rho");

        // compute exchange-correlation functional derivative


        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_pgga(npoints, rho, sigma, exc, vrho, vsigma, rs_omega);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.row(thread_id);

            for (int iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) 
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac = local_weights[g] * vrho[2 * g + 0];
                    double prefac_x = local_weights[g] * (vsigma[3 * g + 0] * 2.0 * rhograd[6 * g + 0]); //+ vsigma[3 * g + 1] * rhograd[6 * g + 3]);
                    double prefac_y = local_weights[g] * (vsigma[3 * g + 0] * 2.0 * rhograd[6 * g + 1]); //+ vsigma[3 * g + 1] * rhograd[6 * g + 4]);
                    double prefac_z = local_weights[g] * (vsigma[3 * g + 0] * 2.0 * rhograd[6 * g + 2]); //+ vsigma[3 * g + 1] * rhograd[6 * g + 5]);

                    double prefacpi = local_weights[g] * vrho[2 * g + 1];
                    //double prefacpi_x = local_weights[g] * (vsigma[3 * g + 2] * 2.0 * rhograd[6 * g + 3] + vsigma[3 * g + 1] * rhograd[6 * g + 0]);
                    //double prefacpi_y = local_weights[g] * (vsigma[3 * g + 2] * 2.0 * rhograd[6 * g + 4] + vsigma[3 * g + 1] * rhograd[6 * g + 1]);
                    //double prefacpi_z = local_weights[g] * (vsigma[3 * g + 2] * 2.0 * rhograd[6 * g + 5] + vsigma[3 * g + 1] * rhograd[6 * g + 2]);

                    gatmx += prefac * gdenx[atom_g] + prefacpi * gdenpi_x[atom_g];
                    gatmy += prefac * gdeny[atom_g] + prefacpi * gdenpi_y[atom_g];
                    gatmz += prefac * gdenz[atom_g] + prefacpi * gdenpi_z[atom_g];

                    gatmx += prefac_x * gdenxx[atom_g] + prefac_y * gdenyx[atom_g] + prefac_z * gdenzx[atom_g];
                    gatmy += prefac_x * gdenxy[atom_g] + prefac_y * gdenyy[atom_g] + prefac_z * gdenzy[atom_g];
                    gatmz += prefac_x * gdenxz[atom_g] + prefac_y * gdenyz[atom_g] + prefac_z * gdenzz[atom_g];

                    //gatmx += prefacpi_x * gdenpi_xx[atom_g] + prefacpi_y * gdenpi_yx[atom_g] + prefacpi_z * gdenpi_zx[atom_g];
                    //gatmy += prefacpi_x * gdenpi_xy[atom_g] + prefacpi_y * gdenpi_yy[atom_g] + prefacpi_z * gdenpi_zy[atom_g];
                    //gatmz += prefacpi_x * gdenpi_xz[atom_g] + prefacpi_y * gdenpi_yz[atom_g] + prefacpi_z * gdenpi_zz[atom_g];
                }

                gatm[iatom * 3 + 0] += gatmx;
                gatm[iatom * 3 + 1] += gatmy;
                gatm[iatom * 3 + 2] += gatmz;
            }
        }

        timer.stop("Accumulate gradient");
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

}   // namespace xcgradpdft
