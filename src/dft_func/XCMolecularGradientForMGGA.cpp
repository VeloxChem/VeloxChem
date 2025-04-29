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

#include "XCMolecularGradientForMGGA.hpp"

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
#include "XCMolecularGradientForGGA.hpp"

namespace xcgradmgga {  // xcgradmgga namespace

auto
integrateVxcGradientForMetaGgaClosedShell(const CMolecule&        molecule,
                                          const CMolecularBasis&  basis,
                                          const std::vector<const double*>& rwDensityPointers,
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

    auto       mggafunc = xcFunctional.getFunctionalPointerToMetaGgaComponent();
    const auto dim      = &(mggafunc->dim);

    std::vector<CXCFunctional> omp_xcfuncs(nthreads, CXCFunctional(xcFunctional));

    std::vector<std::vector<double>> omp_local_weights_data(nthreads, std::vector<double>(omp_max_npoints));

    std::vector<std::vector<double>> omp_rho_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));
    std::vector<std::vector<double>> omp_rhograd_data(nthreads, std::vector<double>(dim->rho * 3 * omp_max_npoints));
    std::vector<std::vector<double>> omp_sigma_data(nthreads, std::vector<double>(dim->sigma * omp_max_npoints));
    std::vector<std::vector<double>> omp_lapl_data(nthreads, std::vector<double>(dim->lapl * omp_max_npoints));
    std::vector<std::vector<double>> omp_tau_data(nthreads, std::vector<double>(dim->tau * omp_max_npoints));

    std::vector<std::vector<double>> omp_vrho_data(nthreads, std::vector<double>(dim->vrho * omp_max_npoints));
    std::vector<std::vector<double>> omp_vsigma_data(nthreads, std::vector<double>(dim->vsigma * omp_max_npoints));
    std::vector<std::vector<double>> omp_vlapl_data(nthreads, std::vector<double>(dim->vlapl * omp_max_npoints));
    std::vector<std::vector<double>> omp_vtau_data(nthreads, std::vector<double>(dim->vtau * omp_max_npoints));

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

        // generate sub density matrix and density grid

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);
        auto rw_sub_dens_mat = dftsubmat::getSubDensityMatrix(rwDensityPointers[0], aoinds, naos);

        timer.stop("Density matrix slicing");

        // GTO values on grid points

        timer.start("OMP VxcGrad calc.");

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

            CDenseMatrix mat_chi_xx(aocount, grid_batch_size);
            CDenseMatrix mat_chi_xy(aocount, grid_batch_size);
            CDenseMatrix mat_chi_xz(aocount, grid_batch_size);
            CDenseMatrix mat_chi_yy(aocount, grid_batch_size);
            CDenseMatrix mat_chi_yz(aocount, grid_batch_size);
            CDenseMatrix mat_chi_zz(aocount, grid_batch_size);

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
                    std::memcpy(mat_chi.row(idx), submat_0_data + nu * grid_batch_size, grid_batch_size * sizeof(double));

                    std::memcpy(mat_chi_x.row(idx), submat_x_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_y.row(idx), submat_y_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_z.row(idx), submat_z_data + nu * grid_batch_size, grid_batch_size * sizeof(double));

                    std::memcpy(mat_chi_xx.row(idx), submat_xx_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xy.row(idx), submat_xy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xz.row(idx), submat_xz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yy.row(idx), submat_yy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yz.row(idx), submat_yz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_zz.row(idx), submat_zz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
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

            auto vrho   = omp_vrho_data[thread_id].data();
            auto vsigma = omp_vsigma_data[thread_id].data();
            auto vlapl  = omp_vlapl_data[thread_id].data();
            auto vtau   = omp_vtau_data[thread_id].data();

            sdengridgen::serialGenerateDensityForMGGA(rho, rhograd, sigma, lapl, tau, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat);

            omptimers[thread_id].stop("Generate density grid");

            // generate density gradient grid

            omptimers[thread_id].start("Density grad. grid prep.");

            CDenseMatrix dengradx(natoms, grid_batch_size);
            CDenseMatrix dengrady(natoms, grid_batch_size);
            CDenseMatrix dengradz(natoms, grid_batch_size);

            CDenseMatrix dengradxx(natoms, grid_batch_size);
            CDenseMatrix dengradxy(natoms, grid_batch_size);
            CDenseMatrix dengradxz(natoms, grid_batch_size);

            CDenseMatrix dengradyx(natoms, grid_batch_size);
            CDenseMatrix dengradyy(natoms, grid_batch_size);
            CDenseMatrix dengradyz(natoms, grid_batch_size);

            CDenseMatrix dengradzx(natoms, grid_batch_size);
            CDenseMatrix dengradzy(natoms, grid_batch_size);
            CDenseMatrix dengradzz(natoms, grid_batch_size);

            CDenseMatrix taugradx(natoms, grid_batch_size);
            CDenseMatrix taugrady(natoms, grid_batch_size);
            CDenseMatrix taugradz(natoms, grid_batch_size);

            omptimers[thread_id].stop("Density grad. grid prep.");

            // eq.(26), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid matmul");

            auto mat_F = sdenblas::serialMultAB(rw_sub_dens_mat, mat_chi);

            auto mat_F_x = sdenblas::serialMultAB(rw_sub_dens_mat, mat_chi_x);
            auto mat_F_y = sdenblas::serialMultAB(rw_sub_dens_mat, mat_chi_y);
            auto mat_F_z = sdenblas::serialMultAB(rw_sub_dens_mat, mat_chi_z);

            omptimers[thread_id].stop("Density grad. grid matmul");

            // eq.(34), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid rho");

            auto F_val = mat_F.values();

            auto F_x_val = mat_F_x.values();
            auto F_y_val = mat_F_y.values();
            auto F_z_val = mat_F_z.values();

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

            auto gtaux = taugradx.values();
            auto gtauy = taugrady.values();
            auto gtauz = taugradz.values();

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

                    gdenxx[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xx_val[nu_g]);
                    gdenxy[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_xy_val[nu_g]);
                    gdenxz[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_xz_val[nu_g]);

                    gdenyx[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xy_val[nu_g]);
                    gdenyy[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_yy_val[nu_g]);
                    gdenyz[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_yz_val[nu_g]);

                    gdenzx[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xz_val[nu_g]);
                    gdenzy[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_yz_val[nu_g]);
                    gdenzz[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_zz_val[nu_g]);

                    gtaux[atom_g] -= (F_x_val[nu_g] * chi_xx_val[nu_g] +
                                      F_y_val[nu_g] * chi_xy_val[nu_g] +
                                      F_z_val[nu_g] * chi_xz_val[nu_g]);
                    gtauy[atom_g] -= (F_x_val[nu_g] * chi_xy_val[nu_g] +
                                      F_y_val[nu_g] * chi_yy_val[nu_g] +
                                      F_z_val[nu_g] * chi_yz_val[nu_g]);
                    gtauz[atom_g] -= (F_x_val[nu_g] * chi_xz_val[nu_g] +
                                      F_y_val[nu_g] * chi_yz_val[nu_g] +
                                      F_z_val[nu_g] * chi_zz_val[nu_g]);
                }
            }

            omptimers[thread_id].stop("Density grad. grid rho");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_vxc_for_mgga(grid_batch_size, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau);

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

                    auto vx = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                    auto vy = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                    auto vz = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                    gatmx += local_weights[g] * vrho[2 * g + 0] * gdenx[atom_g];
                    gatmy += local_weights[g] * vrho[2 * g + 0] * gdeny[atom_g];
                    gatmz += local_weights[g] * vrho[2 * g + 0] * gdenz[atom_g];

                    gatmx += local_weights[g] * (vx * gdenxx[atom_g] + vy * gdenyx[atom_g] + vz * gdenzx[atom_g]);
                    gatmy += local_weights[g] * (vx * gdenxy[atom_g] + vy * gdenyy[atom_g] + vz * gdenzy[atom_g]);
                    gatmz += local_weights[g] * (vx * gdenxz[atom_g] + vy * gdenyz[atom_g] + vz * gdenzz[atom_g]);

                    gatmx += local_weights[g] * vtau[2 * g + 0] * gtaux[atom_g];
                    gatmy += local_weights[g] * vtau[2 * g + 0] * gtauy[atom_g];
                    gatmz += local_weights[g] * vtau[2 * g + 0] * gtauz[atom_g];
                }

                // factor of 2 from sum of alpha and beta contributions

                gatm[iatom * 3 + 0] += 2.0 * gatmx;
                gatm[iatom * 3 + 1] += 2.0 * gatmy;
                gatm[iatom * 3 + 2] += 2.0 * gatmz;
            }

            omptimers[thread_id].stop("Accumulate gradient");
        }

        timer.stop("OMP VxcGrad calc.");
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
integrateVxcGradientForMetaGgaOpenShell(const CMolecule&        molecule,
                                        const CMolecularBasis&  basis,
                                        const std::vector<const double*>& rwDensityPointers,
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

    auto       mggafunc = xcFunctional.getFunctionalPointerToMetaGgaComponent();
    const auto dim      = &(mggafunc->dim);

    std::vector<CXCFunctional> omp_xcfuncs(nthreads, CXCFunctional(xcFunctional));

    std::vector<std::vector<double>> omp_local_weights_data(nthreads, std::vector<double>(omp_max_npoints));

    std::vector<std::vector<double>> omp_rho_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));
    std::vector<std::vector<double>> omp_rhograd_data(nthreads, std::vector<double>(dim->rho * 3 * omp_max_npoints));
    std::vector<std::vector<double>> omp_sigma_data(nthreads, std::vector<double>(dim->sigma * omp_max_npoints));
    std::vector<std::vector<double>> omp_lapl_data(nthreads, std::vector<double>(dim->lapl * omp_max_npoints));
    std::vector<std::vector<double>> omp_tau_data(nthreads, std::vector<double>(dim->tau * omp_max_npoints));

    std::vector<std::vector<double>> omp_vrho_data(nthreads, std::vector<double>(dim->vrho * omp_max_npoints));
    std::vector<std::vector<double>> omp_vsigma_data(nthreads, std::vector<double>(dim->vsigma * omp_max_npoints));
    std::vector<std::vector<double>> omp_vlapl_data(nthreads, std::vector<double>(dim->vlapl * omp_max_npoints));
    std::vector<std::vector<double>> omp_vtau_data(nthreads, std::vector<double>(dim->vtau * omp_max_npoints));

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

        // generate sub density matrix and density grid

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat_a = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);
        auto gs_sub_dens_mat_b = dftsubmat::getSubDensityMatrix(gsDensityPointers[1], aoinds, naos);

        auto rw_sub_dens_mat_a = dftsubmat::getSubDensityMatrix(rwDensityPointers[0], aoinds, naos);
        auto rw_sub_dens_mat_b = dftsubmat::getSubDensityMatrix(rwDensityPointers[1], aoinds, naos);

        timer.stop("Density matrix slicing");

        // GTO values on grid points

        timer.start("OMP VxcGrad calc.");

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

            CDenseMatrix mat_chi_xx(aocount, grid_batch_size);
            CDenseMatrix mat_chi_xy(aocount, grid_batch_size);
            CDenseMatrix mat_chi_xz(aocount, grid_batch_size);
            CDenseMatrix mat_chi_yy(aocount, grid_batch_size);
            CDenseMatrix mat_chi_yz(aocount, grid_batch_size);
            CDenseMatrix mat_chi_zz(aocount, grid_batch_size);


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
                    std::memcpy(mat_chi.row(idx), submat_0_data + nu * grid_batch_size, grid_batch_size * sizeof(double));

                    std::memcpy(mat_chi_x.row(idx), submat_x_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_y.row(idx), submat_y_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_z.row(idx), submat_z_data + nu * grid_batch_size, grid_batch_size * sizeof(double));

                    std::memcpy(mat_chi_xx.row(idx), submat_xx_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xy.row(idx), submat_xy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xz.row(idx), submat_xz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yy.row(idx), submat_yy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yz.row(idx), submat_yz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_zz.row(idx), submat_zz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
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

            auto vrho   = omp_vrho_data[thread_id].data();
            auto vsigma = omp_vsigma_data[thread_id].data();
            auto vlapl  = omp_vlapl_data[thread_id].data();
            auto vtau   = omp_vtau_data[thread_id].data();

            sdengridgen::serialGenerateDensityForMGGA(rho, rhograd, sigma, lapl, tau, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat_a, gs_sub_dens_mat_b);

            omptimers[thread_id].stop("Generate density grid");

            // generate density gradient grid

            omptimers[thread_id].start("Density grad. grid prep.");

            CDenseMatrix dengrad_a_x(natoms, grid_batch_size);
            CDenseMatrix dengrad_a_y(natoms, grid_batch_size);
            CDenseMatrix dengrad_a_z(natoms, grid_batch_size);

            CDenseMatrix dengrad_b_x(natoms, grid_batch_size);
            CDenseMatrix dengrad_b_y(natoms, grid_batch_size);
            CDenseMatrix dengrad_b_z(natoms, grid_batch_size);

            CDenseMatrix dengrad_a_xx(natoms, grid_batch_size);
            CDenseMatrix dengrad_a_xy(natoms, grid_batch_size);
            CDenseMatrix dengrad_a_xz(natoms, grid_batch_size);

            CDenseMatrix dengrad_b_xx(natoms, grid_batch_size);
            CDenseMatrix dengrad_b_xy(natoms, grid_batch_size);
            CDenseMatrix dengrad_b_xz(natoms, grid_batch_size);

            CDenseMatrix dengrad_a_yx(natoms, grid_batch_size);
            CDenseMatrix dengrad_a_yy(natoms, grid_batch_size);
            CDenseMatrix dengrad_a_yz(natoms, grid_batch_size);

            CDenseMatrix dengrad_b_yx(natoms, grid_batch_size);
            CDenseMatrix dengrad_b_yy(natoms, grid_batch_size);
            CDenseMatrix dengrad_b_yz(natoms, grid_batch_size);

            CDenseMatrix dengrad_a_zx(natoms, grid_batch_size);
            CDenseMatrix dengrad_a_zy(natoms, grid_batch_size);
            CDenseMatrix dengrad_a_zz(natoms, grid_batch_size);

            CDenseMatrix dengrad_b_zx(natoms, grid_batch_size);
            CDenseMatrix dengrad_b_zy(natoms, grid_batch_size);
            CDenseMatrix dengrad_b_zz(natoms, grid_batch_size);

            CDenseMatrix taugrad_a_x(natoms, grid_batch_size);
            CDenseMatrix taugrad_a_y(natoms, grid_batch_size);
            CDenseMatrix taugrad_a_z(natoms, grid_batch_size);

            CDenseMatrix taugrad_b_x(natoms, grid_batch_size);
            CDenseMatrix taugrad_b_y(natoms, grid_batch_size);
            CDenseMatrix taugrad_b_z(natoms, grid_batch_size);

            omptimers[thread_id].stop("Density grad. grid prep.");

            // eq.(26), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid matmul");

            auto mat_F_a = sdenblas::serialMultAB(rw_sub_dens_mat_a, mat_chi);
            auto mat_F_b = sdenblas::serialMultAB(rw_sub_dens_mat_b, mat_chi);

            auto mat_F_a_x = sdenblas::serialMultAB(rw_sub_dens_mat_a, mat_chi_x);
            auto mat_F_a_y = sdenblas::serialMultAB(rw_sub_dens_mat_a, mat_chi_y);
            auto mat_F_a_z = sdenblas::serialMultAB(rw_sub_dens_mat_a, mat_chi_z);

            auto mat_F_b_x = sdenblas::serialMultAB(rw_sub_dens_mat_b, mat_chi_x);
            auto mat_F_b_y = sdenblas::serialMultAB(rw_sub_dens_mat_b, mat_chi_y);
            auto mat_F_b_z = sdenblas::serialMultAB(rw_sub_dens_mat_b, mat_chi_z);

            omptimers[thread_id].stop("Density grad. grid matmul");

            // eq.(34), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid rho");

            auto F_a_val = mat_F_a.values();
            auto F_b_val = mat_F_b.values();

            auto F_a_x_val = mat_F_a_x.values();
            auto F_a_y_val = mat_F_a_y.values();
            auto F_a_z_val = mat_F_a_z.values();

            auto F_b_x_val = mat_F_b_x.values();
            auto F_b_y_val = mat_F_b_y.values();
            auto F_b_z_val = mat_F_b_z.values();

            auto chi_x_val = mat_chi_x.values();
            auto chi_y_val = mat_chi_y.values();
            auto chi_z_val = mat_chi_z.values();

            auto chi_xx_val = mat_chi_xx.values();
            auto chi_xy_val = mat_chi_xy.values();
            auto chi_xz_val = mat_chi_xz.values();

            auto chi_yy_val = mat_chi_yy.values();
            auto chi_yz_val = mat_chi_yz.values();
            auto chi_zz_val = mat_chi_zz.values();

            auto gden_a_x = dengrad_a_x.values();
            auto gden_a_y = dengrad_a_y.values();
            auto gden_a_z = dengrad_a_z.values();

            auto gden_b_x = dengrad_b_x.values();
            auto gden_b_y = dengrad_b_y.values();
            auto gden_b_z = dengrad_b_z.values();

            auto gden_a_xx = dengrad_a_xx.values();
            auto gden_a_xy = dengrad_a_xy.values();
            auto gden_a_xz = dengrad_a_xz.values();

            auto gden_b_xx = dengrad_b_xx.values();
            auto gden_b_xy = dengrad_b_xy.values();
            auto gden_b_xz = dengrad_b_xz.values();

            auto gden_a_yx = dengrad_a_yx.values();
            auto gden_a_yy = dengrad_a_yy.values();
            auto gden_a_yz = dengrad_a_yz.values();

            auto gden_b_yx = dengrad_b_yx.values();
            auto gden_b_yy = dengrad_b_yy.values();
            auto gden_b_yz = dengrad_b_yz.values();

            auto gden_a_zx = dengrad_a_zx.values();
            auto gden_a_zy = dengrad_a_zy.values();
            auto gden_a_zz = dengrad_a_zz.values();

            auto gden_b_zx = dengrad_b_zx.values();
            auto gden_b_zy = dengrad_b_zy.values();
            auto gden_b_zz = dengrad_b_zz.values();

            auto gtau_a_x = taugrad_a_x.values();
            auto gtau_a_y = taugrad_a_y.values();
            auto gtau_a_z = taugrad_a_z.values();

            auto gtau_b_x = taugrad_b_x.values();
            auto gtau_b_y = taugrad_b_y.values();
            auto gtau_b_z = taugrad_b_z.values();

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

                    gden_a_x[atom_g] -= 2.0 * F_a_val[nu_g] * chi_x_val[nu_g];
                    gden_a_y[atom_g] -= 2.0 * F_a_val[nu_g] * chi_y_val[nu_g];
                    gden_a_z[atom_g] -= 2.0 * F_a_val[nu_g] * chi_z_val[nu_g];

                    gden_b_x[atom_g] -= 2.0 * F_b_val[nu_g] * chi_x_val[nu_g];
                    gden_b_y[atom_g] -= 2.0 * F_b_val[nu_g] * chi_y_val[nu_g];
                    gden_b_z[atom_g] -= 2.0 * F_b_val[nu_g] * chi_z_val[nu_g];

                    gden_a_xx[atom_g] -= 2.0 * (F_a_x_val[nu_g] * chi_x_val[nu_g] + F_a_val[nu_g] * chi_xx_val[nu_g]);
                    gden_a_xy[atom_g] -= 2.0 * (F_a_x_val[nu_g] * chi_y_val[nu_g] + F_a_val[nu_g] * chi_xy_val[nu_g]);
                    gden_a_xz[atom_g] -= 2.0 * (F_a_x_val[nu_g] * chi_z_val[nu_g] + F_a_val[nu_g] * chi_xz_val[nu_g]);

                    gden_b_xx[atom_g] -= 2.0 * (F_b_x_val[nu_g] * chi_x_val[nu_g] + F_b_val[nu_g] * chi_xx_val[nu_g]);
                    gden_b_xy[atom_g] -= 2.0 * (F_b_x_val[nu_g] * chi_y_val[nu_g] + F_b_val[nu_g] * chi_xy_val[nu_g]);
                    gden_b_xz[atom_g] -= 2.0 * (F_b_x_val[nu_g] * chi_z_val[nu_g] + F_b_val[nu_g] * chi_xz_val[nu_g]);

                    gden_a_yx[atom_g] -= 2.0 * (F_a_y_val[nu_g] * chi_x_val[nu_g] + F_a_val[nu_g] * chi_xy_val[nu_g]);
                    gden_a_yy[atom_g] -= 2.0 * (F_a_y_val[nu_g] * chi_y_val[nu_g] + F_a_val[nu_g] * chi_yy_val[nu_g]);
                    gden_a_yz[atom_g] -= 2.0 * (F_a_y_val[nu_g] * chi_z_val[nu_g] + F_a_val[nu_g] * chi_yz_val[nu_g]);

                    gden_b_yx[atom_g] -= 2.0 * (F_b_y_val[nu_g] * chi_x_val[nu_g] + F_b_val[nu_g] * chi_xy_val[nu_g]);
                    gden_b_yy[atom_g] -= 2.0 * (F_b_y_val[nu_g] * chi_y_val[nu_g] + F_b_val[nu_g] * chi_yy_val[nu_g]);
                    gden_b_yz[atom_g] -= 2.0 * (F_b_y_val[nu_g] * chi_z_val[nu_g] + F_b_val[nu_g] * chi_yz_val[nu_g]);

                    gden_a_zx[atom_g] -= 2.0 * (F_a_z_val[nu_g] * chi_x_val[nu_g] + F_a_val[nu_g] * chi_xz_val[nu_g]);
                    gden_a_zy[atom_g] -= 2.0 * (F_a_z_val[nu_g] * chi_y_val[nu_g] + F_a_val[nu_g] * chi_yz_val[nu_g]);
                    gden_a_zz[atom_g] -= 2.0 * (F_a_z_val[nu_g] * chi_z_val[nu_g] + F_a_val[nu_g] * chi_zz_val[nu_g]);

                    gden_b_zx[atom_g] -= 2.0 * (F_b_z_val[nu_g] * chi_x_val[nu_g] + F_b_val[nu_g] * chi_xz_val[nu_g]);
                    gden_b_zy[atom_g] -= 2.0 * (F_b_z_val[nu_g] * chi_y_val[nu_g] + F_b_val[nu_g] * chi_yz_val[nu_g]);
                    gden_b_zz[atom_g] -= 2.0 * (F_b_z_val[nu_g] * chi_z_val[nu_g] + F_b_val[nu_g] * chi_zz_val[nu_g]);

                    gtau_a_x[atom_g] -= (F_a_x_val[nu_g] * chi_xx_val[nu_g] +
                                         F_a_y_val[nu_g] * chi_xy_val[nu_g] +
                                         F_a_z_val[nu_g] * chi_xz_val[nu_g]);
                    gtau_a_y[atom_g] -= (F_a_x_val[nu_g] * chi_xy_val[nu_g] +
                                         F_a_y_val[nu_g] * chi_yy_val[nu_g] +
                                         F_a_z_val[nu_g] * chi_yz_val[nu_g]);
                    gtau_a_z[atom_g] -= (F_a_x_val[nu_g] * chi_xz_val[nu_g] +
                                         F_a_y_val[nu_g] * chi_yz_val[nu_g] +
                                         F_a_z_val[nu_g] * chi_zz_val[nu_g]);

                    gtau_b_x[atom_g] -= (F_b_x_val[nu_g] * chi_xx_val[nu_g] +
                                         F_b_y_val[nu_g] * chi_xy_val[nu_g] +
                                         F_b_z_val[nu_g] * chi_xz_val[nu_g]);
                    gtau_b_y[atom_g] -= (F_b_x_val[nu_g] * chi_xy_val[nu_g] +
                                         F_b_y_val[nu_g] * chi_yy_val[nu_g] +
                                         F_b_z_val[nu_g] * chi_yz_val[nu_g]);
                    gtau_b_z[atom_g] -= (F_b_x_val[nu_g] * chi_xz_val[nu_g] +
                                         F_b_y_val[nu_g] * chi_yz_val[nu_g] +
                                         F_b_z_val[nu_g] * chi_zz_val[nu_g]);
                }
            }

            omptimers[thread_id].stop("Density grad. grid rho");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_vxc_for_mgga(grid_batch_size, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau);

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

                    auto vxa = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                    auto vya = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                    auto vza = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                    auto vxb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 3] + vsigma[3 * g + 1] * rhograd[6 * g + 0];
                    auto vyb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 4] + vsigma[3 * g + 1] * rhograd[6 * g + 1];
                    auto vzb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 5] + vsigma[3 * g + 1] * rhograd[6 * g + 2];

                    // alpha
                    gatmx += local_weights[g] * vrho[2 * g + 0] * gden_a_x[atom_g];
                    gatmy += local_weights[g] * vrho[2 * g + 0] * gden_a_y[atom_g];
                    gatmz += local_weights[g] * vrho[2 * g + 0] * gden_a_z[atom_g];

                    // beta
                    gatmx += local_weights[g] * vrho[2 * g + 1] * gden_b_x[atom_g];
                    gatmy += local_weights[g] * vrho[2 * g + 1] * gden_b_y[atom_g];
                    gatmz += local_weights[g] * vrho[2 * g + 1] * gden_b_z[atom_g];

                    // alpha
                    gatmx += local_weights[g] * (vxa * gden_a_xx[atom_g] + vya * gden_a_yx[atom_g] + vza * gden_a_zx[atom_g]);
                    gatmy += local_weights[g] * (vxa * gden_a_xy[atom_g] + vya * gden_a_yy[atom_g] + vza * gden_a_zy[atom_g]);
                    gatmz += local_weights[g] * (vxa * gden_a_xz[atom_g] + vya * gden_a_yz[atom_g] + vza * gden_a_zz[atom_g]);

                    // beta
                    gatmx += local_weights[g] * (vxb * gden_b_xx[atom_g] + vyb * gden_b_yx[atom_g] + vzb * gden_b_zx[atom_g]);
                    gatmy += local_weights[g] * (vxb * gden_b_xy[atom_g] + vyb * gden_b_yy[atom_g] + vzb * gden_b_zy[atom_g]);
                    gatmz += local_weights[g] * (vxb * gden_b_xz[atom_g] + vyb * gden_b_yz[atom_g] + vzb * gden_b_zz[atom_g]);

                    // alpha
                    gatmx += local_weights[g] * vtau[2 * g + 0] * gtau_a_x[atom_g];
                    gatmy += local_weights[g] * vtau[2 * g + 0] * gtau_a_y[atom_g];
                    gatmz += local_weights[g] * vtau[2 * g + 0] * gtau_a_z[atom_g];

                    // beta
                    gatmx += local_weights[g] * vtau[2 * g + 1] * gtau_b_x[atom_g];
                    gatmy += local_weights[g] * vtau[2 * g + 1] * gtau_b_y[atom_g];
                    gatmz += local_weights[g] * vtau[2 * g + 1] * gtau_b_z[atom_g];
                }

                gatm[iatom * 3 + 0] += gatmx;
                gatm[iatom * 3 + 1] += gatmy;
                gatm[iatom * 3 + 2] += gatmz;
            }

            omptimers[thread_id].stop("Accumulate gradient");
        }

        timer.stop("OMP VxcGrad calc.");
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

}  // namespace xcgradmgga
