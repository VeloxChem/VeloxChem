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

#include "XCMolecularGradientForGGA.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "AODensityMatrix.hpp"
#include "AOIndices.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DensityGridGenerator.hpp"
#include "DensityGridQuad.hpp"
#include "DftSubMatrix.hpp"
#include "ErrorHandler.hpp"
#include "FunctionalParser.hpp"
#include "GridScreener.hpp"
#include "GtoFunc.hpp"
#include "GtoValues.hpp"
#include "Prescreener.hpp"
#include "XCFunctional.hpp"
#include "XCMolecularGradientForGGA.hpp"

namespace xcgradgga {  // xcgradgga namespace

auto
integrateVxcGradientForGgaClosedShell(const CMolecule&        molecule,
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

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<CXCFunctional> omp_xcfuncs(nthreads, CXCFunctional(xcFunctional));

    std::vector<std::vector<double>> omp_local_weights_data(nthreads, std::vector<double>(omp_max_npoints));

    std::vector<std::vector<double>> omp_rho_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));
    std::vector<std::vector<double>> omp_rhograd_data(nthreads, std::vector<double>(dim->rho * 3 * omp_max_npoints));
    std::vector<std::vector<double>> omp_sigma_data(nthreads, std::vector<double>(dim->sigma * omp_max_npoints));

    std::vector<std::vector<double>> omp_vrho_data(nthreads, std::vector<double>(dim->vrho * omp_max_npoints));
    std::vector<std::vector<double>> omp_vsigma_data(nthreads, std::vector<double>(dim->vsigma * omp_max_npoints));

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

            auto vrho   = omp_vrho_data[thread_id].data();
            auto vsigma = omp_vsigma_data[thread_id].data();

            dengridgen::serialGenerateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat);

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

            omptimers[thread_id].stop("Density grad. grid prep.");

            // eq.(26), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid matmul");

            auto mat_F = denblas::serialMultAB(rw_sub_dens_mat, mat_chi);

            auto mat_F_x = denblas::serialMultAB(rw_sub_dens_mat, mat_chi_x);
            auto mat_F_y = denblas::serialMultAB(rw_sub_dens_mat, mat_chi_y);
            auto mat_F_z = denblas::serialMultAB(rw_sub_dens_mat, mat_chi_z);

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
                }
            }

            omptimers[thread_id].stop("Density grad. grid rho");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_vxc_for_gga(grid_batch_size, rho, sigma, vrho, vsigma);

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
integrateVxcGradientForGgaOpenShell(const CMolecule&        molecule,
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

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<CXCFunctional> omp_xcfuncs(nthreads, CXCFunctional(xcFunctional));

    std::vector<std::vector<double>> omp_local_weights_data(nthreads, std::vector<double>(omp_max_npoints));

    std::vector<std::vector<double>> omp_rho_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));
    std::vector<std::vector<double>> omp_rhograd_data(nthreads, std::vector<double>(dim->rho * 3 * omp_max_npoints));
    std::vector<std::vector<double>> omp_sigma_data(nthreads, std::vector<double>(dim->sigma * omp_max_npoints));

    std::vector<std::vector<double>> omp_vrho_data(nthreads, std::vector<double>(dim->vrho * omp_max_npoints));
    std::vector<std::vector<double>> omp_vsigma_data(nthreads, std::vector<double>(dim->vsigma * omp_max_npoints));

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

            auto vrho   = omp_vrho_data[thread_id].data();
            auto vsigma = omp_vsigma_data[thread_id].data();

            dengridgen::serialGenerateDensityForGGA(
                rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat_a, gs_sub_dens_mat_b);

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

            omptimers[thread_id].stop("Density grad. grid prep.");

            // eq.(26), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid matmul");

            auto mat_F_a = denblas::serialMultAB(rw_sub_dens_mat_a, mat_chi);
            auto mat_F_b = denblas::serialMultAB(rw_sub_dens_mat_b, mat_chi);

            auto mat_F_a_x = denblas::serialMultAB(rw_sub_dens_mat_a, mat_chi_x);
            auto mat_F_a_y = denblas::serialMultAB(rw_sub_dens_mat_a, mat_chi_y);
            auto mat_F_a_z = denblas::serialMultAB(rw_sub_dens_mat_a, mat_chi_z);

            auto mat_F_b_x = denblas::serialMultAB(rw_sub_dens_mat_b, mat_chi_x);
            auto mat_F_b_y = denblas::serialMultAB(rw_sub_dens_mat_b, mat_chi_y);
            auto mat_F_b_z = denblas::serialMultAB(rw_sub_dens_mat_b, mat_chi_z);

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
                }
            }

            omptimers[thread_id].stop("Density grad. grid rho");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_vxc_for_gga(grid_batch_size, rho, sigma, vrho, vsigma);

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

auto
integrateFxcGradientForGgaClosedShell(const CMolecule&        molecule,
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

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<CXCFunctional> omp_xcfuncs(nthreads, CXCFunctional(xcFunctional));

    std::vector<std::vector<double>> omp_local_weights_data(nthreads, std::vector<double>(omp_max_npoints));

    std::vector<std::vector<double>> omp_rho_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));
    std::vector<std::vector<double>> omp_rhow_data(nthreads, std::vector<double>(dim->rho * omp_max_npoints));

    std::vector<std::vector<double>> omp_rhograd_data(nthreads, std::vector<double>(dim->rho * 3 * omp_max_npoints));
    std::vector<std::vector<double>> omp_rhowgrad_data(nthreads, std::vector<double>(dim->rho * 3 * omp_max_npoints));

    std::vector<std::vector<double>> omp_sigma_data(nthreads, std::vector<double>(dim->sigma * omp_max_npoints));

    std::vector<std::vector<double>> omp_vrho_data(nthreads, std::vector<double>(dim->vrho * omp_max_npoints));
    std::vector<std::vector<double>> omp_vsigma_data(nthreads, std::vector<double>(dim->vsigma * omp_max_npoints));

    std::vector<std::vector<double>> omp_v2rho2_data(nthreads, std::vector<double>(dim->v2rho2 * omp_max_npoints));
    std::vector<std::vector<double>> omp_v2rhosigma_data(nthreads, std::vector<double>(dim->v2rhosigma * omp_max_npoints));
    std::vector<std::vector<double>> omp_v2sigma2_data(nthreads, std::vector<double>(dim->v2sigma2 * omp_max_npoints));

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

            // generate density grid

            omptimers[thread_id].start("Generate density grid");

            auto local_weights = omp_local_weights_data[thread_id].data();

            auto rho    = omp_rho_data[thread_id].data();
            auto rhow   = omp_rhow_data[thread_id].data();

            auto rhograd  = omp_rhograd_data[thread_id].data();
            auto rhowgrad = omp_rhowgrad_data[thread_id].data();

            auto sigma  = omp_sigma_data[thread_id].data();

            auto vrho   = omp_vrho_data[thread_id].data();
            auto vsigma = omp_vsigma_data[thread_id].data();

            auto v2rho2     = omp_v2rho2_data[thread_id].data();
            auto v2rhosigma = omp_v2rhosigma_data[thread_id].data();
            auto v2sigma2   = omp_v2sigma2_data[thread_id].data();

            dengridgen::serialGenerateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat);

            dengridgen::serialGenerateDensityForGGA(rhow, rhowgrad, nullptr, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rw_sub_dens_mat_one);

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

            omptimers[thread_id].stop("Density grad. grid prep.");

            // eq.(26), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid matmul");

            auto mat_F = denblas::serialMultAB(rw_sub_dens_mat_two, mat_chi);

            auto mat_F_x = denblas::serialMultAB(rw_sub_dens_mat_two, mat_chi_x);
            auto mat_F_y = denblas::serialMultAB(rw_sub_dens_mat_two, mat_chi_y);
            auto mat_F_z = denblas::serialMultAB(rw_sub_dens_mat_two, mat_chi_z);

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
                }
            }

            omptimers[thread_id].stop("Density grad. grid rho");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_vxc_for_gga(grid_batch_size, rho, sigma, vrho, vsigma);

            omp_xcfuncs[thread_id].compute_fxc_for_gga(grid_batch_size, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

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

                    double w = local_weights[g];

                    auto grhow_grho_aa = 2.0 * (rhowgrad[6 * g + 0] * rhograd[6 * g + 0] + rhowgrad[6 * g + 1] * rhograd[6 * g + 1] +
                                                rhowgrad[6 * g + 2] * rhograd[6 * g + 2]);

                    auto grhow_grho_bb = 2.0 * (rhowgrad[6 * g + 3] * rhograd[6 * g + 3] + rhowgrad[6 * g + 4] * rhograd[6 * g + 4] +
                                                rhowgrad[6 * g + 5] * rhograd[6 * g + 5]);

                    auto grhow_grho_ab = (rhowgrad[6 * g + 0] * rhograd[6 * g + 3] + rhowgrad[6 * g + 1] * rhograd[6 * g + 4] +
                                          rhowgrad[6 * g + 2] * rhograd[6 * g + 5] +

                                          rhowgrad[6 * g + 3] * rhograd[6 * g + 0] + rhowgrad[6 * g + 4] * rhograd[6 * g + 1] +
                                          rhowgrad[6 * g + 5] * rhograd[6 * g + 2]);

                    // scalar contribution, \nabla_A (\phi_mu \phi_nu)

                    double f_0 = v2rho2[3 * g + 0] * rhow[2 * g + 0] + v2rho2[3 * g + 1] * rhow[2 * g + 1] +

                                 v2rhosigma[6 * g + 0] * grhow_grho_aa + v2rhosigma[6 * g + 1] * grhow_grho_ab +
                                 v2rhosigma[6 * g + 2] * grhow_grho_bb;

                    gatmx += w * f_0 * gdenx[atom_g];
                    gatmy += w * f_0 * gdeny[atom_g];
                    gatmz += w * f_0 * gdenz[atom_g];

                    // vector contribution, \nabla_A (\nabla (\phi_mu \phi_nu))

                    double f_aa = v2rhosigma[6 * g + 0] * rhow[2 * g + 0] + v2rhosigma[6 * g + 3] * rhow[2 * g + 1] +

                                  v2sigma2[6 * g + 0] * grhow_grho_aa + v2sigma2[6 * g + 1] * grhow_grho_ab + v2sigma2[6 * g + 2] * grhow_grho_bb;

                    double f_ab = v2rhosigma[6 * g + 1] * rhow[2 * g + 0] + v2rhosigma[6 * g + 4] * rhow[2 * g + 1] +

                                  v2sigma2[6 * g + 1] * grhow_grho_aa + v2sigma2[6 * g + 3] * grhow_grho_ab + v2sigma2[6 * g + 4] * grhow_grho_bb;

                    double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                    xcomp += 2.0 * f_aa * rhograd[6 * g + 0] + f_ab * rhograd[6 * g + 3];
                    ycomp += 2.0 * f_aa * rhograd[6 * g + 1] + f_ab * rhograd[6 * g + 4];
                    zcomp += 2.0 * f_aa * rhograd[6 * g + 2] + f_ab * rhograd[6 * g + 5];

                    xcomp += 2.0 * vsigma[3 * g + 0] * rhowgrad[6 * g + 0] + vsigma[3 * g + 1] * rhowgrad[6 * g + 3];
                    ycomp += 2.0 * vsigma[3 * g + 0] * rhowgrad[6 * g + 1] + vsigma[3 * g + 1] * rhowgrad[6 * g + 4];
                    zcomp += 2.0 * vsigma[3 * g + 0] * rhowgrad[6 * g + 2] + vsigma[3 * g + 1] * rhowgrad[6 * g + 5];

                    gatmx += w * (xcomp * gdenxx[atom_g] + ycomp * gdenyx[atom_g] + zcomp * gdenzx[atom_g]);
                    gatmy += w * (xcomp * gdenxy[atom_g] + ycomp * gdenyy[atom_g] + zcomp * gdenzy[atom_g]);
                    gatmz += w * (xcomp * gdenxz[atom_g] + ycomp * gdenyz[atom_g] + zcomp * gdenzz[atom_g]);
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
integrateKxcGradientForGgaClosedShell(const CMolecule&        molecule,
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

            dengridgen::serialGenerateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat);

            omptimers[thread_id].stop("Generate density grid");

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

            auto rwdengrid = dengridgen::serialGenerateDensityGridForGGA(mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rwdenmat, xcfuntype);

            CDensityGridQuad rwdengridquad(grid_batch_size, numdens_rw2, xcfuntype, dengrid::ab);

            rwdengridquad.DensityProd(rwdengrid, xcfuntype, numdens_rw2, quadMode);

            omptimers[thread_id].stop("Density grid quad");

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

            omptimers[thread_id].stop("Density grad. grid prep.");

            // eq.(26), JCTC 2021, 17, 1512-1521

            omptimers[thread_id].start("Density grad. grid matmul");

            auto mat_F = denblas::serialMultAB(gs_sub_dens_mat, mat_chi);

            auto mat_F_x = denblas::serialMultAB(gs_sub_dens_mat, mat_chi_x);
            auto mat_F_y = denblas::serialMultAB(gs_sub_dens_mat, mat_chi_y);
            auto mat_F_z = denblas::serialMultAB(gs_sub_dens_mat, mat_chi_z);

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
                }
            }

            omptimers[thread_id].stop("Density grad. grid rho");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_vxc_for_gga(grid_batch_size, rho, sigma, vrho, vsigma);

            omp_xcfuncs[thread_id].compute_fxc_for_gga(grid_batch_size, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

            omp_xcfuncs[thread_id].compute_kxc_for_gga(grid_batch_size, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            // pointers to perturbed densities

            auto rhow1rhow2 = rwdengridquad.gam(0);

            auto rxw1rhow2 = rwdengridquad.gamX(0);
            auto ryw1rhow2 = rwdengridquad.gamY(0);
            auto rzw1rhow2 = rwdengridquad.gamZ(0);

            auto rxw1rxw2 = rwdengridquad.gamXX(0);
            auto rxw1ryw2 = rwdengridquad.gamXY(0);
            auto rxw1rzw2 = rwdengridquad.gamXZ(0);

            auto ryw1rxw2 = rwdengridquad.gamYX(0);
            auto ryw1ryw2 = rwdengridquad.gamYY(0);
            auto ryw1rzw2 = rwdengridquad.gamYZ(0);

            auto rzw1rxw2 = rwdengridquad.gamZX(0);
            auto rzw1ryw2 = rwdengridquad.gamZY(0);
            auto rzw1rzw2 = rwdengridquad.gamZZ(0);

            // Note: rw2DensityMatrix is zero in KxcGradientForGGA
            // auto rhow12a = rw2DensityGrid.alphaDensity(idensity);
            // auto gradw12a_x = rw2DensityGrid.alphaDensityGradientX(idensity);
            // auto gradw12a_y = rw2DensityGrid.alphaDensityGradientY(idensity);
            // auto gradw12a_z = rw2DensityGrid.alphaDensityGradientZ(idensity);

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

                    double w = local_weights[g];

                    // double rxw12a = gradw12a_x[g];
                    // double ryw12a = gradw12a_y[g];
                    // double rzw12a = gradw12a_z[g];

                    double grada_x_g = rhograd[6 * g + 0];
                    double grada_y_g = rhograd[6 * g + 1];
                    double grada_z_g = rhograd[6 * g + 2];

                    // double l2contract = grada_x_g * rxw12a + grada_y_g * ryw12a + grada_z_g * rzw12a;
                    // double l5contract_x = grada_x_g * l2contract;
                    // double l5contract_y = grada_y_g * l2contract;
                    // double l5contract_z = grada_z_g * l2contract;

                    double q2contract = grada_x_g * rxw1rhow2[g] + grada_y_g * ryw1rhow2[g] + grada_z_g * rzw1rhow2[g];

                    double q3contract =
                        grada_x_g * grada_x_g * rxw1rxw2[g] + grada_x_g * grada_y_g * rxw1ryw2[g] + grada_x_g * grada_z_g * rxw1rzw2[g] +
                        grada_y_g * grada_x_g * ryw1rxw2[g] + grada_y_g * grada_y_g * ryw1ryw2[g] + grada_y_g * grada_z_g * ryw1rzw2[g] +
                        grada_z_g * grada_x_g * rzw1rxw2[g] + grada_z_g * grada_y_g * rzw1ryw2[g] + grada_z_g * grada_z_g * rzw1rzw2[g];

                    double q4contract = rxw1rxw2[g] + ryw1ryw2[g] + rzw1rzw2[g];

                    double q7contract_x = grada_x_g * (grada_x_g * rxw1rhow2[g] + grada_y_g * ryw1rhow2[g] + grada_z_g * rzw1rhow2[g]);
                    double q7contract_y = grada_y_g * (grada_x_g * rxw1rhow2[g] + grada_y_g * ryw1rhow2[g] + grada_z_g * rzw1rhow2[g]);
                    double q7contract_z = grada_z_g * (grada_x_g * rxw1rhow2[g] + grada_y_g * ryw1rhow2[g] + grada_z_g * rzw1rhow2[g]);

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

                    // Scalar contribution

                    double prefac = 0.0;

                    // vxc 1 contributions

                    // L1
                    // v2rho2_aa = v2rho2[3 * g + 0];
                    // v2rho2_ab = v2rho2[3 * g + 1];
                    // prefac += (v2rho2[3 * g + 0] + v2rho2[3 * g + 1]) * rhow12a[g];

                    // L2
                    // v2rhosigma_aa = v2rhosigma[6 * g + 0];
                    // v2rhosigma_ac = v2rhosigma[6 * g + 1];
                    // v2rhosigma_ab = v2rhosigma[6 * g + 2];
                    // prefac += 2.0 * (v2rhosigma[6 * g + 0] + v2rhosigma[6 * g + 1] + v2rhosigma[6 * g + 2]) * l2contract;

                    // vxc 2 contributions

                    // Q1
                    // v3rho3_aaa = v3rho3[4 * g + 0];
                    // v3rho3_aab = v3rho3[4 * g + 1];
                    // v3rho3_abb = v3rho3[4 * g + 2];
                    prefac += (v3rho3[4 * g + 0] + 2.0 * v3rho3[4 * g + 1] + v3rho3[4 * g + 2]) * rhow1rhow2[g];

                    // Q2
                    // v3rho2sigma_aaa = v3rho2sigma[9 * g + 0];
                    // v3rho2sigma_aac = v3rho2sigma[9 * g + 1];
                    // v3rho2sigma_aab = v3rho2sigma[9 * g + 2];
                    // v3rho2sigma_aba = v3rho2sigma[9 * g + 3];
                    // v3rho2sigma_abc = v3rho2sigma[9 * g + 4];
                    // v3rho2sigma_abb = v3rho2sigma[9 * g + 5];
                    prefac += 2.0 *
                              (v3rho2sigma[9 * g + 0] + v3rho2sigma[9 * g + 1] + v3rho2sigma[9 * g + 2] + v3rho2sigma[9 * g + 3] +
                               v3rho2sigma[9 * g + 4] + v3rho2sigma[9 * g + 5]) *
                              q2contract;

                    // Q3
                    // v3rhosigma2_aaa = v3rhosigma2[12 * g + 0];
                    // v3rhosigma2_aac = v3rhosigma2[12 * g + 1];
                    // v3rhosigma2_aab = v3rhosigma2[12 * g + 2];
                    // v3rhosigma2_acc = v3rhosigma2[12 * g + 3];
                    // v3rhosigma2_acb = v3rhosigma2[12 * g + 4];
                    // v3rhosigma2_abb = v3rhosigma2[12 * g + 5];
                    prefac += 4.0 *
                              (v3rhosigma2[12 * g + 0] + 2.0 * v3rhosigma2[12 * g + 1] + 2.0 * v3rhosigma2[12 * g + 2] + v3rhosigma2[12 * g + 3] +
                               2.0 * v3rhosigma2[12 * g + 4] + v3rhosigma2[12 * g + 5]) *
                              q3contract;

                    // Q4
                    // v2rhosigma_aa = v2rhosigma[6 * g + 0];
                    // v2rhosigma_ac = v2rhosigma[6 * g + 1];
                    // v2rhosigma_ab = v2rhosigma[6 * g + 2];
                    prefac += 2.0 * (v2rhosigma[6 * g + 0] + v2rhosigma[6 * g + 1] + v2rhosigma[6 * g + 2]) * q4contract;

                    gatmx += w * prefac * gdenx[atom_g];
                    gatmy += w * prefac * gdeny[atom_g];
                    gatmz += w * prefac * gdenz[atom_g];

                    // vector contribution

                    double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                    // vxc 1 contributions

                    // L3
                    // v2rhosigma_aa = v2rhosigma[6 * g + 0];
                    // v2rhosigma_ac = v2rhosigma[6 * g + 1];
                    // v2rhosigma_ba = v2rhosigma[6 * g + 3];
                    // v2rhosigma_bc = v2rhosigma[6 * g + 4];
                    // double l3 = 2.0*v2rhosigma[6 * g + 0] + v2rhosigma[6 * g + 1] +
                    //             2.0*v2rhosigma[6 * g + 3] + v2rhosigma[6 * g + 4];
                    // xcomp += l3 * grada_x_g * rhow12a[g];
                    // ycomp += l3 * grada_y_g * rhow12a[g];
                    // zcomp += l3 * grada_z_g * rhow12a[g];

                    // L4
                    // vsigma_a = vsigma[3 * g + 0];
                    // vsigma_c = vsigma[3 * g + 1];
                    // double l4 = 2.0*vsigma[3 * g + 0] + vsigma[3 * g + 1];
                    // xcomp += l4 * rxw12a;
                    // ycomp += l4 * ryw12a;
                    // zcomp += l4 * rzw12a;

                    // L5
                    // v2sigma2_aa = v2sigma2[6 * g + 0];
                    // v2sigma2_ac = v2sigma2[6 * g + 1];
                    // v2sigma2_ab = v2sigma2[6 * g + 2];
                    // v2sigma2_cc = v2sigma2[6 * g + 3];
                    // v2sigma2_cb = v2sigma2[6 * g + 4];
                    // double l5 = 4.0*v2sigma2[6 * g + 0] + 6.0*v2sigma2[6 * g + 1] + 4.0*v2sigma2[6 * g + 2] +
                    //             2.0*v2sigma2[6 * g + 3] + 2.0*v2sigma2[6 * g + 4];
                    // xcomp += l5 * l5contract_x;
                    // ycomp += l5 * l5contract_y;
                    // zcomp += l5 * l5contract_z;

                    // vxc 2 contributions

                    // Q5
                    // v3rho2sigma_aaa = v3rho2sigma[9 * g + 0];
                    // v3rho2sigma_aac = v3rho2sigma[9 * g + 1];
                    // v3rho2sigma_aba = v3rho2sigma[9 * g + 3];
                    // v3rho2sigma_abc = v3rho2sigma[9 * g + 4];
                    // v3rho2sigma_bba = v3rho2sigma[9 * g + 6];
                    // v3rho2sigma_bbc = v3rho2sigma[9 * g + 7];
                    double q5 = 2.0 * v3rho2sigma[9 * g + 0] + v3rho2sigma[9 * g + 1] + 4.0 * v3rho2sigma[9 * g + 3] + 2.0 * v3rho2sigma[9 * g + 4] +
                                2.0 * v3rho2sigma[9 * g + 6] + v3rho2sigma[9 * g + 7];

                    xcomp += q5 * grada_x_g * rhow1rhow2[g];
                    ycomp += q5 * grada_y_g * rhow1rhow2[g];
                    zcomp += q5 * grada_z_g * rhow1rhow2[g];

                    // Q6
                    // v2rhosigma_aa = v2rhosigma[6 * g + 0];
                    // v2rhosigma_ac = v2rhosigma[6 * g + 1];
                    // v2rhosigma_ba = v2rhosigma[6 * g + 3];
                    // v2rhosigma_bc = v2rhosigma[6 * g + 4];
                    double q6 = 2.0 * v2rhosigma[6 * g + 0] + v2rhosigma[6 * g + 1] + 2.0 * v2rhosigma[6 * g + 3] + v2rhosigma[6 * g + 4];

                    xcomp += q6 * rxw1rhow2[g];
                    ycomp += q6 * ryw1rhow2[g];
                    zcomp += q6 * rzw1rhow2[g];

                    // Q7
                    // v3rhosigma2_aaa = v3rhosigma2[12 * g + 0];
                    // v3rhosigma2_aac = v3rhosigma2[12 * g + 1];
                    // v3rhosigma2_aab = v3rhosigma2[12 * g + 2];
                    // v3rhosigma2_acc = v3rhosigma2[12 * g + 3];
                    // v3rhosigma2_acb = v3rhosigma2[12 * g + 4];
                    // v3rhosigma2_baa = v3rhosigma2[12 * g + 6];
                    // v3rhosigma2_bac = v3rhosigma2[12 * g + 7];
                    // v3rhosigma2_bab = v3rhosigma2[12 * g + 8];
                    // v3rhosigma2_bcc = v3rhosigma2[12 * g + 9];
                    // v3rhosigma2_bcb = v3rhosigma2[12 * g + 10];
                    double q7 =
                        2.0 * (2.0 * v3rhosigma2[12 * g + 0] + 3.0 * v3rhosigma2[12 * g + 1] + 2.0 * v3rhosigma2[12 * g + 2] +
                               v3rhosigma2[12 * g + 3] + v3rhosigma2[12 * g + 4] + 2.0 * v3rhosigma2[12 * g + 6] + 3.0 * v3rhosigma2[12 * g + 7] +
                               2.0 * v3rhosigma2[12 * g + 8] + v3rhosigma2[12 * g + 9] + v3rhosigma2[12 * g + 10]);

                    xcomp += q7 * q7contract_x;
                    ycomp += q7 * q7contract_y;
                    zcomp += q7 * q7contract_z;

                    // Q8
                    // v2sigma2_aa = v2sigma2[6 * g + 0];
                    // v2sigma2_ac = v2sigma2[6 * g + 1];
                    // v2sigma2_ab = v2sigma2[6 * g + 2];
                    // v2sigma2_cc = v2sigma2[6 * g + 3];
                    // v2sigma2_cb = v2sigma2[6 * g + 4];
                    double q8 = 4.0 * v2sigma2[6 * g + 0] + 6.0 * v2sigma2[6 * g + 1] + 4.0 * v2sigma2[6 * g + 2] + 2.0 * v2sigma2[6 * g + 3] +
                                2.0 * v2sigma2[6 * g + 4];

                    xcomp += q8 * (q8contract_x + q10contract_x + q11contract_x);
                    ycomp += q8 * (q8contract_y + q10contract_y + q11contract_y);
                    zcomp += q8 * (q8contract_z + q10contract_z + q11contract_z);

                    // Q9
                    // v3sigma3_aaa = v3sigma3[10 * g + 0];
                    // v3sigma3_aac = v3sigma3[10 * g + 1];
                    // v3sigma3_aab = v3sigma3[10 * g + 2];
                    // v3sigma3_acc = v3sigma3[10 * g + 3];
                    // v3sigma3_acb = v3sigma3[10 * g + 4];
                    // v3sigma3_abb = v3sigma3[10 * g + 5];
                    // v3sigma3_ccc = v3sigma3[10 * g + 6];
                    // v3sigma3_ccb = v3sigma3[10 * g + 7];
                    // v3sigma3_cbb = v3sigma3[10 * g + 8];
                    double q9 = 8.0 * v3sigma3[10 * g + 0] + 20.0 * v3sigma3[10 * g + 1] + 16.0 * v3sigma3[10 * g + 2] + 16.0 * v3sigma3[10 * g + 3] +
                                24.0 * v3sigma3[10 * g + 4] + 8.0 * v3sigma3[10 * g + 5] + 4.0 * v3sigma3[10 * g + 6] + 8.0 * v3sigma3[10 * g + 7] +
                                4.0 * v3sigma3[10 * g + 8];

                    xcomp += q9 * q9contract_x;
                    ycomp += q9 * q9contract_y;
                    zcomp += q9 * q9contract_z;

                    gatmx += w * (xcomp * gdenxx[atom_g] + ycomp * gdenyx[atom_g] + zcomp * gdenzx[atom_g]);
                    gatmy += w * (xcomp * gdenxy[atom_g] + ycomp * gdenyy[atom_g] + zcomp * gdenzy[atom_g]);
                    gatmz += w * (xcomp * gdenxz[atom_g] + ycomp * gdenyz[atom_g] + zcomp * gdenzz[atom_g]);
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

}  // namespace xcgradgga
