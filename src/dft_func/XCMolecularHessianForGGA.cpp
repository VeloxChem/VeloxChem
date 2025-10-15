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

#include "XCMolecularHessianForGGA.hpp"

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
#include "MultiTimer.hpp"
#include "Prescreener.hpp"
#include "SerialDenseLinearAlgebra.hpp"
#include "SerialDensityGridGenerator.hpp"
#include "XCFunctional.hpp"

namespace xchessgga {  // xchessgga namespace

auto
integrateExcHessianForGgaClosedShell(const CMolecule&        molecule,
                                     const CMolecularBasis&  basis,
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

    // molecular Hessian

    auto natoms = molecule.number_of_atoms();

    CDenseMatrix molhess_threads(nthreads, (natoms * 3) * (natoms * 3));

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

    std::vector<std::vector<double>> omp_weighted_vrho(nthreads, std::vector<double>(omp_max_npoints));
    std::vector<std::vector<double>> omp_weighted_vnabla_x(nthreads, std::vector<double>(omp_max_npoints));
    std::vector<std::vector<double>> omp_weighted_vnabla_y(nthreads, std::vector<double>(omp_max_npoints));
    std::vector<std::vector<double>> omp_weighted_vnabla_z(nthreads, std::vector<double>(omp_max_npoints));

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
            // 3rd order GTO derivative
            auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(gto_block, 3, screeningThresholdForGTOValues, boxdim);

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

        timer.stop("Density matrix slicing");

        // GTO values on grid points

        timer.start("OMP Vxc Hessian evaluation");

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

            CDenseMatrix mat_chi_xxx(aocount, grid_batch_size);
            CDenseMatrix mat_chi_xxy(aocount, grid_batch_size);
            CDenseMatrix mat_chi_xxz(aocount, grid_batch_size);
            CDenseMatrix mat_chi_xyy(aocount, grid_batch_size);
            CDenseMatrix mat_chi_xyz(aocount, grid_batch_size);
            CDenseMatrix mat_chi_xzz(aocount, grid_batch_size);
            CDenseMatrix mat_chi_yyy(aocount, grid_batch_size);
            CDenseMatrix mat_chi_yyz(aocount, grid_batch_size);
            CDenseMatrix mat_chi_yzz(aocount, grid_batch_size);
            CDenseMatrix mat_chi_zzz(aocount, grid_batch_size);

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

                auto cmat = gtoval::get_gto_values_for_3rd_order(gto_block, grid_x, grid_y, grid_z, cgto_mask);

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

                auto submat_xxx_ptr = cmat.sub_matrix({3, 0});
                auto submat_xxy_ptr = cmat.sub_matrix({3, 1});
                auto submat_xxz_ptr = cmat.sub_matrix({3, 2});
                auto submat_xyy_ptr = cmat.sub_matrix({3, 3});
                auto submat_xyz_ptr = cmat.sub_matrix({3, 4});
                auto submat_xzz_ptr = cmat.sub_matrix({3, 5});
                auto submat_yyy_ptr = cmat.sub_matrix({3, 6});
                auto submat_yyz_ptr = cmat.sub_matrix({3, 7});
                auto submat_yzz_ptr = cmat.sub_matrix({3, 8});
                auto submat_zzz_ptr = cmat.sub_matrix({3, 9});

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

                auto submat_xxx_data = submat_xxx_ptr->data();
                auto submat_xxy_data = submat_xxy_ptr->data();
                auto submat_xxz_data = submat_xxz_ptr->data();
                auto submat_xyy_data = submat_xyy_ptr->data();
                auto submat_xyz_data = submat_xyz_ptr->data();
                auto submat_xzz_data = submat_xzz_ptr->data();
                auto submat_yyy_data = submat_yyy_ptr->data();
                auto submat_yyz_data = submat_yyz_ptr->data();
                auto submat_yzz_data = submat_yzz_ptr->data();
                auto submat_zzz_data = submat_zzz_ptr->data();

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

                    std::memcpy(mat_chi_xxx.row(idx), submat_xxx_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xxy.row(idx), submat_xxy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xxz.row(idx), submat_xxz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xyy.row(idx), submat_xyy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xyz.row(idx), submat_xyz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xzz.row(idx), submat_xzz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yyy.row(idx), submat_yyy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yyz.row(idx), submat_yyz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yzz.row(idx), submat_yzz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_zzz.row(idx), submat_zzz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
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

            auto w0 = omp_weighted_vrho[thread_id].data();
            auto wx = omp_weighted_vnabla_x[thread_id].data();
            auto wy = omp_weighted_vnabla_y[thread_id].data();
            auto wz = omp_weighted_vnabla_z[thread_id].data();

            sdengridgen::serialGenerateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat);

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

            omptimers[thread_id].stop("Density grad. grid prep.");

            omptimers[thread_id].start("Density grad. grid matmul");

            auto mat_F = sdenblas::serialMultAB(gs_sub_dens_mat, mat_chi);

            auto mat_F_x = sdenblas::serialMultAB(gs_sub_dens_mat, mat_chi_x);
            auto mat_F_y = sdenblas::serialMultAB(gs_sub_dens_mat, mat_chi_y);
            auto mat_F_z = sdenblas::serialMultAB(gs_sub_dens_mat, mat_chi_z);

            omptimers[thread_id].stop("Density grad. grid matmul");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_vxc_for_gga(grid_batch_size, rho, sigma, vrho, vsigma);

            omp_xcfuncs[thread_id].compute_fxc_for_gga(grid_batch_size, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            omptimers[thread_id].start("Accumulate Hessian");

            auto D_val = gs_sub_dens_mat.values();

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

            auto chi_xxx_val = mat_chi_xxx.values();
            auto chi_xxy_val = mat_chi_xxy.values();
            auto chi_xxz_val = mat_chi_xxz.values();
            auto chi_xyy_val = mat_chi_xyy.values();
            auto chi_xyz_val = mat_chi_xyz.values();
            auto chi_xzz_val = mat_chi_xzz.values();
            auto chi_yyy_val = mat_chi_yyy.values();
            auto chi_yyz_val = mat_chi_yyz.values();
            auto chi_yzz_val = mat_chi_yzz.values();
            auto chi_zzz_val = mat_chi_zzz.values();

            auto gatm = molhess_threads.row(thread_id);

            // prepare w0, wx, wy and wz

            #pragma omp simd 
            for (int g = 0; g < grid_batch_size; g++)
            {
                w0[g] = local_weights[g] * vrho[2 * g + 0];

                auto vx = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                auto vy = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                auto vz = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                wx[g] = local_weights[g] * vx;
                wy[g] = local_weights[g] * vy;
                wz[g] = local_weights[g] * vz;
            }

            // prepare gradient grid

            for (int mu = 0; mu < aocount; mu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[mu]];

                auto atom_offset = atomidx * grid_batch_size;

                auto mu_offset = mu * grid_batch_size;

                #pragma omp simd 
                for (int g = 0; g < grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto mu_g = mu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[mu_g] * chi_x_val[mu_g];
                    gdeny[atom_g] -= 2.0 * F_val[mu_g] * chi_y_val[mu_g];
                    gdenz[atom_g] -= 2.0 * F_val[mu_g] * chi_z_val[mu_g];

                    gdenxx[atom_g] -= 2.0 * (F_x_val[mu_g] * chi_x_val[mu_g] + F_val[mu_g] * chi_xx_val[mu_g]);
                    gdenxy[atom_g] -= 2.0 * (F_x_val[mu_g] * chi_y_val[mu_g] + F_val[mu_g] * chi_xy_val[mu_g]);
                    gdenxz[atom_g] -= 2.0 * (F_x_val[mu_g] * chi_z_val[mu_g] + F_val[mu_g] * chi_xz_val[mu_g]);

                    gdenyx[atom_g] -= 2.0 * (F_y_val[mu_g] * chi_x_val[mu_g] + F_val[mu_g] * chi_xy_val[mu_g]);
                    gdenyy[atom_g] -= 2.0 * (F_y_val[mu_g] * chi_y_val[mu_g] + F_val[mu_g] * chi_yy_val[mu_g]);
                    gdenyz[atom_g] -= 2.0 * (F_y_val[mu_g] * chi_z_val[mu_g] + F_val[mu_g] * chi_yz_val[mu_g]);

                    gdenzx[atom_g] -= 2.0 * (F_z_val[mu_g] * chi_x_val[mu_g] + F_val[mu_g] * chi_xz_val[mu_g]);
                    gdenzy[atom_g] -= 2.0 * (F_z_val[mu_g] * chi_y_val[mu_g] + F_val[mu_g] * chi_yz_val[mu_g]);
                    gdenzz[atom_g] -= 2.0 * (F_z_val[mu_g] * chi_z_val[mu_g] + F_val[mu_g] * chi_zz_val[mu_g]);
                }
            }

            // first contribution to
            // f_{\rho_{\alpha}} \rho_{\alpha}^{(\xi,\zeta)}
            // and
            // (2 f_{\sigma_{\alpha\alpha}} + f_{\sigma_{\alpha\beta}})
            // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\xi,\zeta)}
            // on the same atom

            for (int mu = 0; mu < aocount; mu++)
            {
                auto iatom = ao_to_atom_ids[aoinds[mu]];

                auto ix = iatom * 3 + 0;
                auto iy = iatom * 3 + 1;
                auto iz = iatom * 3 + 2;

                auto mu_offset = mu * grid_batch_size;

                double gatmxx = 0.0, gatmxy = 0.0, gatmxz = 0.0;
                double gatmyx = 0.0, gatmyy = 0.0, gatmyz = 0.0;
                double gatmzx = 0.0, gatmzy = 0.0, gatmzz = 0.0;

                #pragma omp simd reduction(+ : gatmxx, gatmxy, gatmxz, gatmyx, gatmyy, gatmyz, gatmzx, gatmzy, gatmzz) 
                for (int g = 0; g < grid_batch_size; g++)
                {
                    auto mu_g = mu_offset + g;

                    // f_{\rho_{\alpha}} \rho_{\alpha}^{(\xi,\zeta)}

                    // \rho_{\alpha}^{(\xi,\zeta)} (first contrib.)
                    // = 2 \sum_{\mu\nu} D^{\alpha}_{\mu\nu,\rm{sym}} \phi_{\mu}^{(\xi,\zeta)} \phi_{\nu}
                    // = 2 \sum_{\mu} F_{\mu} \phi_{\mu}^{(\xi,\zeta)}

                    // factor of 2 is added outside of the for loop

                    double gdenxx = F_val[mu_g] * chi_xx_val[mu_g];
                    double gdenxy = F_val[mu_g] * chi_xy_val[mu_g];
                    double gdenxz = F_val[mu_g] * chi_xz_val[mu_g];

                    double gdenyx = F_val[mu_g] * chi_xy_val[mu_g];
                    double gdenyy = F_val[mu_g] * chi_yy_val[mu_g];
                    double gdenyz = F_val[mu_g] * chi_yz_val[mu_g];

                    double gdenzx = F_val[mu_g] * chi_xz_val[mu_g];
                    double gdenzy = F_val[mu_g] * chi_yz_val[mu_g];
                    double gdenzz = F_val[mu_g] * chi_zz_val[mu_g];

                    // (2 f_{\sigma_{\alpha\alpha}} + f_{\sigma_{\alpha\beta}})
                    // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\xi,\zeta)}

                    // \nabla\rho_{\alpha}^{(\xi,\zeta)} (first contrib.)
                    // = 2 \sum_{\mu\nu} D^{\alpha}_{\mu\nu,\rm{sym}} (\nabla\phi_{\mu})^{(\xi,\zeta)} \phi_{\nu}
                    // + 2 \sum_{\mu\nu} D^{\alpha}_{\mu\nu,\rm{sym}} \phi_{\mu}^{(\xi,\zeta)} \nabla\phi_{\nu}
                    // = 2 \sum_{\mu} F_{\mu} (\nabla\phi_{\mu})^{(\xi,\zeta)}
                    // + 2 \sum_{\mu} (\nabla F_{\mu}) \phi_{\mu}^{(\xi,\zeta)}

                    // factor of 2 is added outside of the for loop

                    // ordering of components: nabla, xi, zeta

                    // === x ===

                    double gdenxxx = F_x_val[mu_g] * chi_xx_val[mu_g] + F_val[mu_g] * chi_xxx_val[mu_g];
                    double gdenxxy = F_x_val[mu_g] * chi_xy_val[mu_g] + F_val[mu_g] * chi_xxy_val[mu_g];
                    double gdenxxz = F_x_val[mu_g] * chi_xz_val[mu_g] + F_val[mu_g] * chi_xxz_val[mu_g];

                    double gdenxyx = F_x_val[mu_g] * chi_xy_val[mu_g] + F_val[mu_g] * chi_xxy_val[mu_g];
                    double gdenxyy = F_x_val[mu_g] * chi_yy_val[mu_g] + F_val[mu_g] * chi_xyy_val[mu_g];
                    double gdenxyz = F_x_val[mu_g] * chi_yz_val[mu_g] + F_val[mu_g] * chi_xyz_val[mu_g];

                    double gdenxzx = F_x_val[mu_g] * chi_xz_val[mu_g] + F_val[mu_g] * chi_xxz_val[mu_g];
                    double gdenxzy = F_x_val[mu_g] * chi_yz_val[mu_g] + F_val[mu_g] * chi_xyz_val[mu_g];
                    double gdenxzz = F_x_val[mu_g] * chi_zz_val[mu_g] + F_val[mu_g] * chi_xzz_val[mu_g];

                    // === y ===

                    double gdenyxx = F_y_val[mu_g] * chi_xx_val[mu_g] + F_val[mu_g] * chi_xxy_val[mu_g];
                    double gdenyxy = F_y_val[mu_g] * chi_xy_val[mu_g] + F_val[mu_g] * chi_xyy_val[mu_g];
                    double gdenyxz = F_y_val[mu_g] * chi_xz_val[mu_g] + F_val[mu_g] * chi_xyz_val[mu_g];

                    double gdenyyx = F_y_val[mu_g] * chi_xy_val[mu_g] + F_val[mu_g] * chi_xyy_val[mu_g];
                    double gdenyyy = F_y_val[mu_g] * chi_yy_val[mu_g] + F_val[mu_g] * chi_yyy_val[mu_g];
                    double gdenyyz = F_y_val[mu_g] * chi_yz_val[mu_g] + F_val[mu_g] * chi_yyz_val[mu_g];

                    double gdenyzx = F_y_val[mu_g] * chi_xz_val[mu_g] + F_val[mu_g] * chi_xyz_val[mu_g];
                    double gdenyzy = F_y_val[mu_g] * chi_yz_val[mu_g] + F_val[mu_g] * chi_yyz_val[mu_g];
                    double gdenyzz = F_y_val[mu_g] * chi_zz_val[mu_g] + F_val[mu_g] * chi_yzz_val[mu_g];

                    // === z ===

                    double gdenzxx = F_z_val[mu_g] * chi_xx_val[mu_g] + F_val[mu_g] * chi_xxz_val[mu_g];
                    double gdenzxy = F_z_val[mu_g] * chi_xy_val[mu_g] + F_val[mu_g] * chi_xyz_val[mu_g];
                    double gdenzxz = F_z_val[mu_g] * chi_xz_val[mu_g] + F_val[mu_g] * chi_xzz_val[mu_g];

                    double gdenzyx = F_z_val[mu_g] * chi_xy_val[mu_g] + F_val[mu_g] * chi_xyz_val[mu_g];
                    double gdenzyy = F_z_val[mu_g] * chi_yy_val[mu_g] + F_val[mu_g] * chi_yyz_val[mu_g];
                    double gdenzyz = F_z_val[mu_g] * chi_yz_val[mu_g] + F_val[mu_g] * chi_yzz_val[mu_g];

                    double gdenzzx = F_z_val[mu_g] * chi_xz_val[mu_g] + F_val[mu_g] * chi_xzz_val[mu_g];
                    double gdenzzy = F_z_val[mu_g] * chi_yz_val[mu_g] + F_val[mu_g] * chi_yzz_val[mu_g];
                    double gdenzzz = F_z_val[mu_g] * chi_zz_val[mu_g] + F_val[mu_g] * chi_zzz_val[mu_g];

                    // accumulate contribution

                    gatmxx += w0[g] * gdenxx + (wx[g] * gdenxxx + wy[g] * gdenyxx + wz[g] * gdenzxx);
                    gatmxy += w0[g] * gdenxy + (wx[g] * gdenxxy + wy[g] * gdenyxy + wz[g] * gdenzxy);
                    gatmxz += w0[g] * gdenxz + (wx[g] * gdenxxz + wy[g] * gdenyxz + wz[g] * gdenzxz);

                    gatmyx += w0[g] * gdenyx + (wx[g] * gdenxyx + wy[g] * gdenyyx + wz[g] * gdenzyx);
                    gatmyy += w0[g] * gdenyy + (wx[g] * gdenxyy + wy[g] * gdenyyy + wz[g] * gdenzyy);
                    gatmyz += w0[g] * gdenyz + (wx[g] * gdenxyz + wy[g] * gdenyyz + wz[g] * gdenzyz);

                    gatmzx += w0[g] * gdenzx + (wx[g] * gdenxzx + wy[g] * gdenyzx + wz[g] * gdenzzx);
                    gatmzy += w0[g] * gdenzy + (wx[g] * gdenxzy + wy[g] * gdenyzy + wz[g] * gdenzzy);
                    gatmzz += w0[g] * gdenzz + (wx[g] * gdenxzz + wy[g] * gdenyzz + wz[g] * gdenzzz);
                }

                // factor of 2 from differentiation
                // factor of 2 from sum of alpha and beta contributions

                gatm[ix * (natoms * 3) + ix] += 4.0 * gatmxx;
                gatm[ix * (natoms * 3) + iy] += 4.0 * gatmxy;
                gatm[ix * (natoms * 3) + iz] += 4.0 * gatmxz;

                gatm[iy * (natoms * 3) + ix] += 4.0 * gatmyx;
                gatm[iy * (natoms * 3) + iy] += 4.0 * gatmyy;
                gatm[iy * (natoms * 3) + iz] += 4.0 * gatmyz;

                gatm[iz * (natoms * 3) + ix] += 4.0 * gatmzx;
                gatm[iz * (natoms * 3) + iy] += 4.0 * gatmzy;
                gatm[iz * (natoms * 3) + iz] += 4.0 * gatmzz;
            }

            // second contribution to
            // f_{\rho_{\alpha}} \rho_{\alpha}^{(\xi,\zeta)}
            // and
            // (2 f_{\sigma_{\alpha\alpha}} + f_{\sigma_{\alpha\beta}})
            // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\xi,\zeta)}
            // on the same atom and on different atoms

            for (int mu = 0; mu < aocount; mu++)
            {
                auto iatom = ao_to_atom_ids[aoinds[mu]];

                auto ix = iatom * 3 + 0;
                auto iy = iatom * 3 + 1;
                auto iz = iatom * 3 + 2;

                auto mu_offset = mu * grid_batch_size;

                for (int nu = 0; nu < aocount; nu++)
                {
                    auto jatom = ao_to_atom_ids[aoinds[nu]];

                    // only consider the upper triangular part, i.e. iatom <= jatom

                    if (iatom > jatom) continue;

                    auto jx = jatom * 3 + 0;
                    auto jy = jatom * 3 + 1;
                    auto jz = jatom * 3 + 2;

                    auto nu_offset = nu * grid_batch_size;

                    double gatmxx = 0.0, gatmxy = 0.0, gatmxz = 0.0;
                    double gatmyx = 0.0, gatmyy = 0.0, gatmyz = 0.0;
                    double gatmzx = 0.0, gatmzy = 0.0, gatmzz = 0.0;

                    #pragma omp simd reduction(+ : gatmxx, gatmxy, gatmxz, gatmyx, gatmyy, gatmyz, gatmzx, gatmzy, gatmzz) 
                    for (int g = 0; g < grid_batch_size; g++)
                    {
                        auto mu_g = mu_offset + g;
                        auto nu_g = nu_offset + g;

                        // f_{\rho_{\alpha}} \rho_{\alpha}^{(\xi,\zeta)}

                        // \rho_{\alpha}^{(\xi,\zeta)} (second contrib.)
                        // = 2 \sum_{\mu\nu} D^{\alpha}_{\mu\nu,\rm{sym}} \phi_{\mu}^{(\xi)} \phi_{\nu}^{(\zeta)}

                        // factor of 2 and D_{\mu\nu,sym}^{\alpha} are added outside of the for loop

                        double gxx = chi_x_val[mu_g] * chi_x_val[nu_g];
                        double gxy = chi_x_val[mu_g] * chi_y_val[nu_g];
                        double gxz = chi_x_val[mu_g] * chi_z_val[nu_g];

                        double gyx = chi_y_val[mu_g] * chi_x_val[nu_g];
                        double gyy = chi_y_val[mu_g] * chi_y_val[nu_g];
                        double gyz = chi_y_val[mu_g] * chi_z_val[nu_g];

                        double gzx = chi_z_val[mu_g] * chi_x_val[nu_g];
                        double gzy = chi_z_val[mu_g] * chi_y_val[nu_g];
                        double gzz = chi_z_val[mu_g] * chi_z_val[nu_g];

                        // (2 f_{\sigma_{\alpha\alpha}} + f_{\sigma_{\alpha\beta}})
                        // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\xi,\zeta)}

                        // \nabla\rho_{\alpha}^{(\xi,\zeta)} (second contrib.)
                        // = 2 \sum_{\mu\nu} D^{\alpha}_{\mu\nu,\rm{sym}}
                        //   (\nabla\phi_{\mu}^{(\xi)} \phi_{\nu}^{(\zeta)} +
                        //    \nabla\phi_{\nu}^{(\zeta)} \phi_{\mu}^{(\xi)})

                        // factor of 2 and D^{\alpha}_{\mu\nu,\rm{sym}} are added outside of the for loop

                        // ordering of components: nabla, xi, zeta

                        // === x ===

                        double gxxx = chi_xx_val[mu_g] * chi_x_val[nu_g] + chi_xx_val[nu_g] * chi_x_val[mu_g];
                        double gxxy = chi_xx_val[mu_g] * chi_y_val[nu_g] + chi_xy_val[nu_g] * chi_x_val[mu_g];
                        double gxxz = chi_xx_val[mu_g] * chi_z_val[nu_g] + chi_xz_val[nu_g] * chi_x_val[mu_g];

                        double gxyx = chi_xy_val[mu_g] * chi_x_val[nu_g] + chi_xx_val[nu_g] * chi_y_val[mu_g];
                        double gxyy = chi_xy_val[mu_g] * chi_y_val[nu_g] + chi_xy_val[nu_g] * chi_y_val[mu_g];
                        double gxyz = chi_xy_val[mu_g] * chi_z_val[nu_g] + chi_xz_val[nu_g] * chi_y_val[mu_g];

                        double gxzx = chi_xz_val[mu_g] * chi_x_val[nu_g] + chi_xx_val[nu_g] * chi_z_val[mu_g];
                        double gxzy = chi_xz_val[mu_g] * chi_y_val[nu_g] + chi_xy_val[nu_g] * chi_z_val[mu_g];
                        double gxzz = chi_xz_val[mu_g] * chi_z_val[nu_g] + chi_xz_val[nu_g] * chi_z_val[mu_g];

                        // === y ===

                        double gyxx = chi_xy_val[mu_g] * chi_x_val[nu_g] + chi_xy_val[nu_g] * chi_x_val[mu_g];
                        double gyxy = chi_xy_val[mu_g] * chi_y_val[nu_g] + chi_yy_val[nu_g] * chi_x_val[mu_g];
                        double gyxz = chi_xy_val[mu_g] * chi_z_val[nu_g] + chi_yz_val[nu_g] * chi_x_val[mu_g];

                        double gyyx = chi_yy_val[mu_g] * chi_x_val[nu_g] + chi_xy_val[nu_g] * chi_y_val[mu_g];
                        double gyyy = chi_yy_val[mu_g] * chi_y_val[nu_g] + chi_yy_val[nu_g] * chi_y_val[mu_g];
                        double gyyz = chi_yy_val[mu_g] * chi_z_val[nu_g] + chi_yz_val[nu_g] * chi_y_val[mu_g];

                        double gyzx = chi_yz_val[mu_g] * chi_x_val[nu_g] + chi_xy_val[nu_g] * chi_z_val[mu_g];
                        double gyzy = chi_yz_val[mu_g] * chi_y_val[nu_g] + chi_yy_val[nu_g] * chi_z_val[mu_g];
                        double gyzz = chi_yz_val[mu_g] * chi_z_val[nu_g] + chi_yz_val[nu_g] * chi_z_val[mu_g];

                        // === z ===

                        double gzxx = chi_xz_val[mu_g] * chi_x_val[nu_g] + chi_xz_val[nu_g] * chi_x_val[mu_g];
                        double gzxy = chi_xz_val[mu_g] * chi_y_val[nu_g] + chi_yz_val[nu_g] * chi_x_val[mu_g];
                        double gzxz = chi_xz_val[mu_g] * chi_z_val[nu_g] + chi_zz_val[nu_g] * chi_x_val[mu_g];

                        double gzyx = chi_yz_val[mu_g] * chi_x_val[nu_g] + chi_xz_val[nu_g] * chi_y_val[mu_g];
                        double gzyy = chi_yz_val[mu_g] * chi_y_val[nu_g] + chi_yz_val[nu_g] * chi_y_val[mu_g];
                        double gzyz = chi_yz_val[mu_g] * chi_z_val[nu_g] + chi_zz_val[nu_g] * chi_y_val[mu_g];

                        double gzzx = chi_zz_val[mu_g] * chi_x_val[nu_g] + chi_xz_val[nu_g] * chi_z_val[mu_g];
                        double gzzy = chi_zz_val[mu_g] * chi_y_val[nu_g] + chi_yz_val[nu_g] * chi_z_val[mu_g];
                        double gzzz = chi_zz_val[mu_g] * chi_z_val[nu_g] + chi_zz_val[nu_g] * chi_z_val[mu_g];

                        // accumulate contributions

                        gatmxx += w0[g] * gxx + (wx[g] * gxxx + wy[g] * gyxx + wz[g] * gzxx);
                        gatmxy += w0[g] * gxy + (wx[g] * gxxy + wy[g] * gyxy + wz[g] * gzxy);
                        gatmxz += w0[g] * gxz + (wx[g] * gxxz + wy[g] * gyxz + wz[g] * gzxz);

                        gatmyx += w0[g] * gyx + (wx[g] * gxyx + wy[g] * gyyx + wz[g] * gzyx);
                        gatmyy += w0[g] * gyy + (wx[g] * gxyy + wy[g] * gyyy + wz[g] * gzyy);
                        gatmyz += w0[g] * gyz + (wx[g] * gxyz + wy[g] * gyyz + wz[g] * gzyz);

                        gatmzx += w0[g] * gzx + (wx[g] * gxzx + wy[g] * gyzx + wz[g] * gzzx);
                        gatmzy += w0[g] * gzy + (wx[g] * gxzy + wy[g] * gyzy + wz[g] * gzzy);
                        gatmzz += w0[g] * gzz + (wx[g] * gxzz + wy[g] * gyzz + wz[g] * gzzz);
                    }

                    auto D_mn = D_val[mu * aocount + nu];

                    // factor of 2 from differentiation
                    // factor of 2 from sum of alpha and beta contributions

                    gatm[ix * (natoms * 3) + jx] += 4.0 * gatmxx * D_mn;
                    gatm[ix * (natoms * 3) + jy] += 4.0 * gatmxy * D_mn;
                    gatm[ix * (natoms * 3) + jz] += 4.0 * gatmxz * D_mn;

                    gatm[iy * (natoms * 3) + jx] += 4.0 * gatmyx * D_mn;
                    gatm[iy * (natoms * 3) + jy] += 4.0 * gatmyy * D_mn;
                    gatm[iy * (natoms * 3) + jz] += 4.0 * gatmyz * D_mn;

                    gatm[iz * (natoms * 3) + jx] += 4.0 * gatmzx * D_mn;
                    gatm[iz * (natoms * 3) + jy] += 4.0 * gatmzy * D_mn;
                    gatm[iz * (natoms * 3) + jz] += 4.0 * gatmzz * D_mn;
                }
            }

            // other contributions

            for (int iatom = 0; iatom < natoms; iatom++)
            {
                auto ix = iatom * 3 + 0;
                auto iy = iatom * 3 + 1;
                auto iz = iatom * 3 + 2;

                auto i_offset = iatom * grid_batch_size;

                for (int jatom = iatom; jatom < natoms; jatom++)
                {
                    auto jx = jatom * 3 + 0;
                    auto jy = jatom * 3 + 1;
                    auto jz = jatom * 3 + 2;

                    auto j_offset = jatom * grid_batch_size;

                    double gatmxx = 0.0, gatmyx = 0.0, gatmzx = 0.0;
                    double gatmxy = 0.0, gatmyy = 0.0, gatmzy = 0.0;
                    double gatmxz = 0.0, gatmyz = 0.0, gatmzz = 0.0;

                    #pragma omp simd reduction(+ : gatmxx, gatmxy, gatmxz, gatmyx, gatmyy, gatmyz, gatmzx, gatmzy, gatmzz) 
                    for (int g = 0; g < grid_batch_size; g++)
                    {
                        auto ig = i_offset + g;
                        auto jg = j_offset + g;

                        double w = local_weights[g];

                        // (f_{\rho_{\alpha} \rho_{\alpha}} + f_{\rho_{\alpha} \rho_{\beta}})
                        // \rho_{\alpha}^{(\xi)} \rho_{\alpha}^{(\zeta)}

                        double prefac = w * (v2rho2[3 * g + 0] + v2rho2[3 * g + 1]);

                        gatmxx += prefac * gdenx[ig] * gdenx[jg];
                        gatmxy += prefac * gdenx[ig] * gdeny[jg];
                        gatmxz += prefac * gdenx[ig] * gdenz[jg];

                        gatmyx += prefac * gdeny[ig] * gdenx[jg];
                        gatmyy += prefac * gdeny[ig] * gdeny[jg];
                        gatmyz += prefac * gdeny[ig] * gdenz[jg];

                        gatmzx += prefac * gdenz[ig] * gdenx[jg];
                        gatmzy += prefac * gdenz[ig] * gdeny[jg];
                        gatmzz += prefac * gdenz[ig] * gdenz[jg];

                        // (2 f_{\rho_{\alpha} \sigma_{\alpha\alpha}} +
                        //  2 f_{\rho_{\alpha} \sigma_{\alpha\beta}} +
                        //  2 f_{\rho_{\alpha} \sigma_{\beta\beta}})
                        // \rho_{\alpha}^{(\xi)} \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\zeta)}

                        prefac = w * 2.0 * (v2rhosigma[6 * g + 0] + v2rhosigma[6 * g + 1] + v2rhosigma[6 * g + 2]);

                        auto gx = rhograd[6 * g + 0];
                        auto gy = rhograd[6 * g + 1];
                        auto gz = rhograd[6 * g + 2];

                        auto xcomp_j = (gx * gdenxx[jg] + gy * gdenyx[jg] + gz * gdenzx[jg]);
                        auto ycomp_j = (gx * gdenxy[jg] + gy * gdenyy[jg] + gz * gdenzy[jg]);
                        auto zcomp_j = (gx * gdenxz[jg] + gy * gdenyz[jg] + gz * gdenzz[jg]);

                        gatmxx += prefac * gdenx[ig] * xcomp_j;
                        gatmxy += prefac * gdenx[ig] * ycomp_j;
                        gatmxz += prefac * gdenx[ig] * zcomp_j;

                        gatmyx += prefac * gdeny[ig] * xcomp_j;
                        gatmyy += prefac * gdeny[ig] * ycomp_j;
                        gatmyz += prefac * gdeny[ig] * zcomp_j;

                        gatmzx += prefac * gdenz[ig] * xcomp_j;
                        gatmzy += prefac * gdenz[ig] * ycomp_j;
                        gatmzz += prefac * gdenz[ig] * zcomp_j;

                        // (2 f_{\rho_{\alpha} \sigma_{\alpha\alpha}} +
                        //  2 f_{\rho_{\beta} \sigma_{\alpha\alpha}} +
                        //  f_{\rho_{\alpha} \sigma_{\alpha\beta}} +
                        //  f_{\rho_{\beta} \sigma_{\alpha\beta}})
                        // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\xi)} \rho_{\alpha}^{(\zeta)}

                        auto f_aa = v2rhosigma[6 * g + 0] + v2rhosigma[6 * g + 3];
                        auto f_ab = v2rhosigma[6 * g + 1] + v2rhosigma[6 * g + 4];

                        prefac = w * (2.0 * f_aa + f_ab);

                        auto xcomp_i = (gx * gdenxx[ig] + gy * gdenyx[ig] + gz * gdenzx[ig]);
                        auto ycomp_i = (gx * gdenxy[ig] + gy * gdenyy[ig] + gz * gdenzy[ig]);
                        auto zcomp_i = (gx * gdenxz[ig] + gy * gdenyz[ig] + gz * gdenzz[ig]);

                        gatmxx += prefac * xcomp_i * gdenx[jg];
                        gatmxy += prefac * xcomp_i * gdeny[jg];
                        gatmxz += prefac * xcomp_i * gdenz[jg];

                        gatmyx += prefac * ycomp_i * gdenx[jg];
                        gatmyy += prefac * ycomp_i * gdeny[jg];
                        gatmyz += prefac * ycomp_i * gdenz[jg];

                        gatmzx += prefac * zcomp_i * gdenx[jg];
                        gatmzy += prefac * zcomp_i * gdeny[jg];
                        gatmzz += prefac * zcomp_i * gdenz[jg];

                        // (4 f_{\sigma_{\alpha\alpha} \sigma_{\alpha\alpha}} +
                        //  6 f_{\sigma_{\alpha\alpha} \sigma_{\alpha\beta}} +
                        //  4 f_{\sigma_{\alpha\alpha} \sigma_{\beta\beta}} +
                        //  2 f_{\sigma_{\alpha\beta} \sigma_{\alpha\beta}} +
                        //  2 f_{\sigma_{\alpha\beta} \sigma_{\beta\beta}})
                        // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\xi)}
                        // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\zeta)}

                        f_aa = v2sigma2[6 * g + 0] + v2sigma2[6 * g + 1] + v2sigma2[6 * g + 2];
                        f_ab = v2sigma2[6 * g + 1] + v2sigma2[6 * g + 3] + v2sigma2[6 * g + 4];

                        prefac = w * 2.0 * (2.0 * f_aa + f_ab);

                        gatmxx += prefac * xcomp_i * xcomp_j;
                        gatmxy += prefac * xcomp_i * ycomp_j;
                        gatmxz += prefac * xcomp_i * zcomp_j;

                        gatmyx += prefac * ycomp_i * xcomp_j;
                        gatmyy += prefac * ycomp_i * ycomp_j;
                        gatmyz += prefac * ycomp_i * zcomp_j;

                        gatmzx += prefac * zcomp_i * xcomp_j;
                        gatmzy += prefac * zcomp_i * ycomp_j;
                        gatmzz += prefac * zcomp_i * zcomp_j;

                        // (2 f_{\sigma_{\alpha\alpha}} + f_{\sigma_{\alpha\beta}})
                        // \nabla\rho_{\alpha}^{(\xi)} \cdot \nabla\rho_{\alpha}^{(\zeta)}

                        f_aa = vsigma[3 * g + 0];
                        f_ab = vsigma[3 * g + 1];

                        prefac = w * (2.0 * f_aa + f_ab);

                        gatmxx += prefac * (gdenxx[ig] * gdenxx[jg] + gdenyx[ig] * gdenyx[jg] + gdenzx[ig] * gdenzx[jg]);
                        gatmxy += prefac * (gdenxx[ig] * gdenxy[jg] + gdenyx[ig] * gdenyy[jg] + gdenzx[ig] * gdenzy[jg]);
                        gatmxz += prefac * (gdenxx[ig] * gdenxz[jg] + gdenyx[ig] * gdenyz[jg] + gdenzx[ig] * gdenzz[jg]);

                        gatmyx += prefac * (gdenxy[ig] * gdenxx[jg] + gdenyy[ig] * gdenyx[jg] + gdenzy[ig] * gdenzx[jg]);
                        gatmyy += prefac * (gdenxy[ig] * gdenxy[jg] + gdenyy[ig] * gdenyy[jg] + gdenzy[ig] * gdenzy[jg]);
                        gatmyz += prefac * (gdenxy[ig] * gdenxz[jg] + gdenyy[ig] * gdenyz[jg] + gdenzy[ig] * gdenzz[jg]);

                        gatmzx += prefac * (gdenxz[ig] * gdenxx[jg] + gdenyz[ig] * gdenyx[jg] + gdenzz[ig] * gdenzx[jg]);
                        gatmzy += prefac * (gdenxz[ig] * gdenxy[jg] + gdenyz[ig] * gdenyy[jg] + gdenzz[ig] * gdenzy[jg]);
                        gatmzz += prefac * (gdenxz[ig] * gdenxz[jg] + gdenyz[ig] * gdenyz[jg] + gdenzz[ig] * gdenzz[jg]);
                    }

                    // factor of 2 from sum of alpha and beta contributions

                    gatm[ix * (natoms * 3) + jx] += 2.0 * gatmxx;
                    gatm[ix * (natoms * 3) + jy] += 2.0 * gatmxy;
                    gatm[ix * (natoms * 3) + jz] += 2.0 * gatmxz;

                    gatm[iy * (natoms * 3) + jx] += 2.0 * gatmyx;
                    gatm[iy * (natoms * 3) + jy] += 2.0 * gatmyy;
                    gatm[iy * (natoms * 3) + jz] += 2.0 * gatmyz;

                    gatm[iz * (natoms * 3) + jx] += 2.0 * gatmzx;
                    gatm[iz * (natoms * 3) + jy] += 2.0 * gatmzy;
                    gatm[iz * (natoms * 3) + jz] += 2.0 * gatmzz;
                }
            }

            omptimers[thread_id].stop("Accumulate Hessian");
        }

        timer.stop("OMP Vxc Hessian evaluation");
    }

    timer.stop("Total timing");

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;

    CDenseMatrix molhess(natoms * 3, natoms * 3);

    for (int iatom = 0; iatom < natoms; iatom++)
    {
        auto ix = iatom * 3 + 0;
        auto iy = iatom * 3 + 1;
        auto iz = iatom * 3 + 2;

        for (int jatom = iatom; jatom < natoms; jatom++)
        {
            auto jx = jatom * 3 + 0;
            auto jy = jatom * 3 + 1;
            auto jz = jatom * 3 + 2;

            for (int thread_id = 0; thread_id < nthreads; thread_id++)
            {
                molhess.row(ix)[jx] += molhess_threads.row(thread_id)[ix * (natoms * 3) + jx];
                molhess.row(ix)[jy] += molhess_threads.row(thread_id)[ix * (natoms * 3) + jy];
                molhess.row(ix)[jz] += molhess_threads.row(thread_id)[ix * (natoms * 3) + jz];

                molhess.row(iy)[jx] += molhess_threads.row(thread_id)[iy * (natoms * 3) + jx];
                molhess.row(iy)[jy] += molhess_threads.row(thread_id)[iy * (natoms * 3) + jy];
                molhess.row(iy)[jz] += molhess_threads.row(thread_id)[iy * (natoms * 3) + jz];

                molhess.row(iz)[jx] += molhess_threads.row(thread_id)[iz * (natoms * 3) + jx];
                molhess.row(iz)[jy] += molhess_threads.row(thread_id)[iz * (natoms * 3) + jy];
                molhess.row(iz)[jz] += molhess_threads.row(thread_id)[iz * (natoms * 3) + jz];
            }

            if (jatom != iatom)
            {
                molhess.row(jx)[ix] = molhess.row(ix)[jx];
                molhess.row(jx)[iy] = molhess.row(iy)[jx];
                molhess.row(jx)[iz] = molhess.row(iz)[jx];

                molhess.row(jy)[ix] = molhess.row(ix)[jy];
                molhess.row(jy)[iy] = molhess.row(iy)[jy];
                molhess.row(jy)[iz] = molhess.row(iz)[jy];

                molhess.row(jz)[ix] = molhess.row(ix)[jz];
                molhess.row(jz)[iy] = molhess.row(iy)[jz];
                molhess.row(jz)[iz] = molhess.row(iz)[jz];
            }
        }
    }

    return molhess;
}

auto
integrateExcHessianForGgaOpenShell(const CMolecule&        molecule,
                                   const CMolecularBasis&  basis,
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

    // molecular Hessian

    auto natoms = molecule.number_of_atoms();

    CDenseMatrix molhess_threads(nthreads, (natoms * 3) * (natoms * 3));

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

    std::vector<std::vector<double>> omp_weighted_vrho(nthreads, std::vector<double>(omp_max_npoints));
    std::vector<std::vector<double>> omp_weighted_vnabla_x(nthreads, std::vector<double>(omp_max_npoints));
    std::vector<std::vector<double>> omp_weighted_vnabla_y(nthreads, std::vector<double>(omp_max_npoints));
    std::vector<std::vector<double>> omp_weighted_vnabla_z(nthreads, std::vector<double>(omp_max_npoints));

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
            // 3rd order GTO derivative
            auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(gto_block, 3, screeningThresholdForGTOValues, boxdim);

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

        timer.stop("Density matrix slicing");

        // GTO values on grid points

        timer.start("OMP Vxc Hessian evaluation");

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

            CDenseMatrix mat_chi_xxx(aocount, grid_batch_size);
            CDenseMatrix mat_chi_xxy(aocount, grid_batch_size);
            CDenseMatrix mat_chi_xxz(aocount, grid_batch_size);
            CDenseMatrix mat_chi_xyy(aocount, grid_batch_size);
            CDenseMatrix mat_chi_xyz(aocount, grid_batch_size);
            CDenseMatrix mat_chi_xzz(aocount, grid_batch_size);
            CDenseMatrix mat_chi_yyy(aocount, grid_batch_size);
            CDenseMatrix mat_chi_yyz(aocount, grid_batch_size);
            CDenseMatrix mat_chi_yzz(aocount, grid_batch_size);
            CDenseMatrix mat_chi_zzz(aocount, grid_batch_size);

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

                auto cmat = gtoval::get_gto_values_for_3rd_order(gto_block, grid_x, grid_y, grid_z, cgto_mask);

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

                auto submat_xxx_ptr = cmat.sub_matrix({3, 0});
                auto submat_xxy_ptr = cmat.sub_matrix({3, 1});
                auto submat_xxz_ptr = cmat.sub_matrix({3, 2});
                auto submat_xyy_ptr = cmat.sub_matrix({3, 3});
                auto submat_xyz_ptr = cmat.sub_matrix({3, 4});
                auto submat_xzz_ptr = cmat.sub_matrix({3, 5});
                auto submat_yyy_ptr = cmat.sub_matrix({3, 6});
                auto submat_yyz_ptr = cmat.sub_matrix({3, 7});
                auto submat_yzz_ptr = cmat.sub_matrix({3, 8});
                auto submat_zzz_ptr = cmat.sub_matrix({3, 9});

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

                auto submat_xxx_data = submat_xxx_ptr->data();
                auto submat_xxy_data = submat_xxy_ptr->data();
                auto submat_xxz_data = submat_xxz_ptr->data();
                auto submat_xyy_data = submat_xyy_ptr->data();
                auto submat_xyz_data = submat_xyz_ptr->data();
                auto submat_xzz_data = submat_xzz_ptr->data();
                auto submat_yyy_data = submat_yyy_ptr->data();
                auto submat_yyz_data = submat_yyz_ptr->data();
                auto submat_yzz_data = submat_yzz_ptr->data();
                auto submat_zzz_data = submat_zzz_ptr->data();

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

                    std::memcpy(mat_chi_xxx.row(idx), submat_xxx_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xxy.row(idx), submat_xxy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xxz.row(idx), submat_xxz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xyy.row(idx), submat_xyy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xyz.row(idx), submat_xyz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xzz.row(idx), submat_xzz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yyy.row(idx), submat_yyy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yyz.row(idx), submat_yyz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yzz.row(idx), submat_yzz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_zzz.row(idx), submat_zzz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
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

            auto w_a_0 = omp_weighted_vrho[thread_id].data();
            auto w_a_x = omp_weighted_vnabla_x[thread_id].data();
            auto w_a_y = omp_weighted_vnabla_y[thread_id].data();
            auto w_a_z = omp_weighted_vnabla_z[thread_id].data();

            auto w_b_0 = omp_weighted_vrho[thread_id].data();
            auto w_b_x = omp_weighted_vnabla_x[thread_id].data();
            auto w_b_y = omp_weighted_vnabla_y[thread_id].data();
            auto w_b_z = omp_weighted_vnabla_z[thread_id].data();

            sdengridgen::serialGenerateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat_a, gs_sub_dens_mat_b);

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

            CDenseMatrix dengrad_a_yx(natoms, grid_batch_size);
            CDenseMatrix dengrad_a_yy(natoms, grid_batch_size);
            CDenseMatrix dengrad_a_yz(natoms, grid_batch_size);

            CDenseMatrix dengrad_a_zx(natoms, grid_batch_size);
            CDenseMatrix dengrad_a_zy(natoms, grid_batch_size);
            CDenseMatrix dengrad_a_zz(natoms, grid_batch_size);

            CDenseMatrix dengrad_b_xx(natoms, grid_batch_size);
            CDenseMatrix dengrad_b_xy(natoms, grid_batch_size);
            CDenseMatrix dengrad_b_xz(natoms, grid_batch_size);

            CDenseMatrix dengrad_b_yx(natoms, grid_batch_size);
            CDenseMatrix dengrad_b_yy(natoms, grid_batch_size);
            CDenseMatrix dengrad_b_yz(natoms, grid_batch_size);

            CDenseMatrix dengrad_b_zx(natoms, grid_batch_size);
            CDenseMatrix dengrad_b_zy(natoms, grid_batch_size);
            CDenseMatrix dengrad_b_zz(natoms, grid_batch_size);

            auto gden_a_x = dengrad_a_x.values();
            auto gden_a_y = dengrad_a_y.values();
            auto gden_a_z = dengrad_a_z.values();

            auto gden_b_x = dengrad_b_x.values();
            auto gden_b_y = dengrad_b_y.values();
            auto gden_b_z = dengrad_b_z.values();

            auto gden_a_xx = dengrad_a_xx.values();
            auto gden_a_xy = dengrad_a_xy.values();
            auto gden_a_xz = dengrad_a_xz.values();

            auto gden_a_yx = dengrad_a_yx.values();
            auto gden_a_yy = dengrad_a_yy.values();
            auto gden_a_yz = dengrad_a_yz.values();

            auto gden_a_zx = dengrad_a_zx.values();
            auto gden_a_zy = dengrad_a_zy.values();
            auto gden_a_zz = dengrad_a_zz.values();

            auto gden_b_xx = dengrad_b_xx.values();
            auto gden_b_xy = dengrad_b_xy.values();
            auto gden_b_xz = dengrad_b_xz.values();

            auto gden_b_yx = dengrad_b_yx.values();
            auto gden_b_yy = dengrad_b_yy.values();
            auto gden_b_yz = dengrad_b_yz.values();

            auto gden_b_zx = dengrad_b_zx.values();
            auto gden_b_zy = dengrad_b_zy.values();
            auto gden_b_zz = dengrad_b_zz.values();

            omptimers[thread_id].stop("Density grad. grid prep.");

            omptimers[thread_id].start("Density grad. grid matmul");

            auto mat_F_a = sdenblas::serialMultAB(gs_sub_dens_mat_a, mat_chi);
            auto mat_F_b = sdenblas::serialMultAB(gs_sub_dens_mat_b, mat_chi);

            auto mat_F_a_x = sdenblas::serialMultAB(gs_sub_dens_mat_a, mat_chi_x);
            auto mat_F_a_y = sdenblas::serialMultAB(gs_sub_dens_mat_a, mat_chi_y);
            auto mat_F_a_z = sdenblas::serialMultAB(gs_sub_dens_mat_a, mat_chi_z);

            auto mat_F_b_x = sdenblas::serialMultAB(gs_sub_dens_mat_b, mat_chi_x);
            auto mat_F_b_y = sdenblas::serialMultAB(gs_sub_dens_mat_b, mat_chi_y);
            auto mat_F_b_z = sdenblas::serialMultAB(gs_sub_dens_mat_b, mat_chi_z);

            omptimers[thread_id].stop("Density grad. grid matmul");

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_vxc_for_gga(grid_batch_size, rho, sigma, vrho, vsigma);

            omp_xcfuncs[thread_id].compute_fxc_for_gga(grid_batch_size, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            omptimers[thread_id].start("Accumulate Hessian");

            auto D_a_val = gs_sub_dens_mat_a.values();
            auto D_b_val = gs_sub_dens_mat_b.values();

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

            auto chi_xxx_val = mat_chi_xxx.values();
            auto chi_xxy_val = mat_chi_xxy.values();
            auto chi_xxz_val = mat_chi_xxz.values();
            auto chi_xyy_val = mat_chi_xyy.values();
            auto chi_xyz_val = mat_chi_xyz.values();
            auto chi_xzz_val = mat_chi_xzz.values();
            auto chi_yyy_val = mat_chi_yyy.values();
            auto chi_yyz_val = mat_chi_yyz.values();
            auto chi_yzz_val = mat_chi_yzz.values();
            auto chi_zzz_val = mat_chi_zzz.values();

            auto gatm = molhess_threads.row(thread_id);

            // prepare w0, wx, wy and wz

            #pragma omp simd 
            for (int g = 0; g < grid_batch_size; g++)
            {
                w_a_0[g] = local_weights[g] * vrho[2 * g + 0];
                w_b_0[g] = local_weights[g] * vrho[2 * g + 1];

                auto vxa = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                auto vya = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                auto vza = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                w_a_x[g] = local_weights[g] * vxa;
                w_a_y[g] = local_weights[g] * vya;
                w_a_z[g] = local_weights[g] * vza;

                auto vxb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 3] + vsigma[3 * g + 1] * rhograd[6 * g + 0];
                auto vyb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 4] + vsigma[3 * g + 1] * rhograd[6 * g + 1];
                auto vzb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 5] + vsigma[3 * g + 1] * rhograd[6 * g + 2];

                w_b_x[g] = local_weights[g] * vxb;
                w_b_y[g] = local_weights[g] * vyb;
                w_b_z[g] = local_weights[g] * vzb;
            }

            // prepare gradient grid

            for (int mu = 0; mu < aocount; mu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[mu]];

                auto atom_offset = atomidx * grid_batch_size;

                auto mu_offset = mu * grid_batch_size;

                #pragma omp simd 
                for (int g = 0; g < grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto mu_g = mu_offset + g;

                    gden_a_x[atom_g] -= 2.0 * F_a_val[mu_g] * chi_x_val[mu_g];
                    gden_a_y[atom_g] -= 2.0 * F_a_val[mu_g] * chi_y_val[mu_g];
                    gden_a_z[atom_g] -= 2.0 * F_a_val[mu_g] * chi_z_val[mu_g];

                    gden_b_x[atom_g] -= 2.0 * F_b_val[mu_g] * chi_x_val[mu_g];
                    gden_b_y[atom_g] -= 2.0 * F_b_val[mu_g] * chi_y_val[mu_g];
                    gden_b_z[atom_g] -= 2.0 * F_b_val[mu_g] * chi_z_val[mu_g];

                    gden_a_xx[atom_g] -= 2.0 * (F_a_x_val[mu_g] * chi_x_val[mu_g] + F_a_val[mu_g] * chi_xx_val[mu_g]);
                    gden_a_xy[atom_g] -= 2.0 * (F_a_x_val[mu_g] * chi_y_val[mu_g] + F_a_val[mu_g] * chi_xy_val[mu_g]);
                    gden_a_xz[atom_g] -= 2.0 * (F_a_x_val[mu_g] * chi_z_val[mu_g] + F_a_val[mu_g] * chi_xz_val[mu_g]);

                    gden_a_yx[atom_g] -= 2.0 * (F_a_y_val[mu_g] * chi_x_val[mu_g] + F_a_val[mu_g] * chi_xy_val[mu_g]);
                    gden_a_yy[atom_g] -= 2.0 * (F_a_y_val[mu_g] * chi_y_val[mu_g] + F_a_val[mu_g] * chi_yy_val[mu_g]);
                    gden_a_yz[atom_g] -= 2.0 * (F_a_y_val[mu_g] * chi_z_val[mu_g] + F_a_val[mu_g] * chi_yz_val[mu_g]);

                    gden_a_zx[atom_g] -= 2.0 * (F_a_z_val[mu_g] * chi_x_val[mu_g] + F_a_val[mu_g] * chi_xz_val[mu_g]);
                    gden_a_zy[atom_g] -= 2.0 * (F_a_z_val[mu_g] * chi_y_val[mu_g] + F_a_val[mu_g] * chi_yz_val[mu_g]);
                    gden_a_zz[atom_g] -= 2.0 * (F_a_z_val[mu_g] * chi_z_val[mu_g] + F_a_val[mu_g] * chi_zz_val[mu_g]);

                    gden_b_xx[atom_g] -= 2.0 * (F_b_x_val[mu_g] * chi_x_val[mu_g] + F_b_val[mu_g] * chi_xx_val[mu_g]);
                    gden_b_xy[atom_g] -= 2.0 * (F_b_x_val[mu_g] * chi_y_val[mu_g] + F_b_val[mu_g] * chi_xy_val[mu_g]);
                    gden_b_xz[atom_g] -= 2.0 * (F_b_x_val[mu_g] * chi_z_val[mu_g] + F_b_val[mu_g] * chi_xz_val[mu_g]);

                    gden_b_yx[atom_g] -= 2.0 * (F_b_y_val[mu_g] * chi_x_val[mu_g] + F_b_val[mu_g] * chi_xy_val[mu_g]);
                    gden_b_yy[atom_g] -= 2.0 * (F_b_y_val[mu_g] * chi_y_val[mu_g] + F_b_val[mu_g] * chi_yy_val[mu_g]);
                    gden_b_yz[atom_g] -= 2.0 * (F_b_y_val[mu_g] * chi_z_val[mu_g] + F_b_val[mu_g] * chi_yz_val[mu_g]);

                    gden_b_zx[atom_g] -= 2.0 * (F_b_z_val[mu_g] * chi_x_val[mu_g] + F_b_val[mu_g] * chi_xz_val[mu_g]);
                    gden_b_zy[atom_g] -= 2.0 * (F_b_z_val[mu_g] * chi_y_val[mu_g] + F_b_val[mu_g] * chi_yz_val[mu_g]);
                    gden_b_zz[atom_g] -= 2.0 * (F_b_z_val[mu_g] * chi_z_val[mu_g] + F_b_val[mu_g] * chi_zz_val[mu_g]);
                }
            }

            // first contribution to
            // f_{\rho_{\alpha}} \rho_{\alpha}^{(\xi,\zeta)}
            // and
            // (2 f_{\sigma_{\alpha\alpha}} + f_{\sigma_{\alpha\beta}})
            // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\xi,\zeta)}
            // on the same atom

            // (exc_rho[a]) * rho[a]_xi_zeta
            // (exc_rho[b]) * rho[b]_xi_zeta
            // (2 exc_sigma[aa]) * nabla[0]_rho[a]         * nabla[0]_rho[a]_xi_zeta
            // (  exc_sigma[ab]) * nabla[0]_rho[a]_xi_zeta * nabla[0]_rho[b]
            // (  exc_sigma[ab]) * nabla[0]_rho[a]         * nabla[0]_rho[b]_xi_zeta
            // (2 exc_sigma[bb]) * nabla[0]_rho[b]         * nabla[0]_rho[b]_xi_zeta

            for (int mu = 0; mu < aocount; mu++)
            {
                auto iatom = ao_to_atom_ids[aoinds[mu]];

                auto ix = iatom * 3 + 0;
                auto iy = iatom * 3 + 1;
                auto iz = iatom * 3 + 2;

                auto mu_offset = mu * grid_batch_size;

                double gatmxx = 0.0, gatmxy = 0.0, gatmxz = 0.0;
                double gatmyx = 0.0, gatmyy = 0.0, gatmyz = 0.0;
                double gatmzx = 0.0, gatmzy = 0.0, gatmzz = 0.0;

                #pragma omp simd reduction(+ : gatmxx, gatmxy, gatmxz, gatmyx, gatmyy, gatmyz, gatmzx, gatmzy, gatmzz) 
                for (int g = 0; g < grid_batch_size; g++)
                {
                    auto mu_g = mu_offset + g;

                    double w = local_weights[g];

                    // f_{\rho_{\alpha}} \rho_{\alpha}^{(\xi,\zeta)}

                    // \rho_{\alpha}^{(\xi,\zeta)} (first contrib.)
                    // = 2 \sum_{\mu\nu} D^{\alpha}_{\mu\nu,\rm{sym}} \phi_{\mu}^{(\xi,\zeta)} \phi_{\nu}
                    // = 2 \sum_{\mu} F_{\mu} \phi_{\mu}^{(\xi,\zeta)}

                    double gdenxxa = F_a_val[mu_g] * chi_xx_val[mu_g];
                    double gdenxya = F_a_val[mu_g] * chi_xy_val[mu_g];
                    double gdenxza = F_a_val[mu_g] * chi_xz_val[mu_g];

                    double gdenyxa = F_a_val[mu_g] * chi_xy_val[mu_g];
                    double gdenyya = F_a_val[mu_g] * chi_yy_val[mu_g];
                    double gdenyza = F_a_val[mu_g] * chi_yz_val[mu_g];

                    double gdenzxa = F_a_val[mu_g] * chi_xz_val[mu_g];
                    double gdenzya = F_a_val[mu_g] * chi_yz_val[mu_g];
                    double gdenzza = F_a_val[mu_g] * chi_zz_val[mu_g];

                    double gdenxxb = F_b_val[mu_g] * chi_xx_val[mu_g];
                    double gdenxyb = F_b_val[mu_g] * chi_xy_val[mu_g];
                    double gdenxzb = F_b_val[mu_g] * chi_xz_val[mu_g];

                    double gdenyxb = F_b_val[mu_g] * chi_xy_val[mu_g];
                    double gdenyyb = F_b_val[mu_g] * chi_yy_val[mu_g];
                    double gdenyzb = F_b_val[mu_g] * chi_yz_val[mu_g];

                    double gdenzxb = F_b_val[mu_g] * chi_xz_val[mu_g];
                    double gdenzyb = F_b_val[mu_g] * chi_yz_val[mu_g];
                    double gdenzzb = F_b_val[mu_g] * chi_zz_val[mu_g];

                    // (2 f_{\sigma_{\alpha\alpha}} + f_{\sigma_{\alpha\beta}})
                    // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\xi,\zeta)}

                    // \nabla\rho_{\alpha}^{(\xi,\zeta)} (first contrib.)
                    // = 2 \sum_{\mu\nu} D^{\alpha}_{\mu\nu,\rm{sym}} (\nabla\phi_{\mu})^{(\xi,\zeta)} \phi_{\nu}
                    // + 2 \sum_{\mu\nu} D^{\alpha}_{\mu\nu,\rm{sym}} \phi_{\mu}^{(\xi,\zeta)} \nabla\phi_{\nu}
                    // = 2 \sum_{\mu} F_{\mu} (\nabla\phi_{\mu})^{(\xi,\zeta)}
                    // + 2 \sum_{\mu} (\nabla F_{\mu}) \phi_{\mu}^{(\xi,\zeta)}

                    // ordering of components: nabla, xi, zeta

                    // === x ===

                    double gdenxxxa = F_a_x_val[mu_g] * chi_xx_val[mu_g] + F_a_val[mu_g] * chi_xxx_val[mu_g];
                    double gdenxxya = F_a_x_val[mu_g] * chi_xy_val[mu_g] + F_a_val[mu_g] * chi_xxy_val[mu_g];
                    double gdenxxza = F_a_x_val[mu_g] * chi_xz_val[mu_g] + F_a_val[mu_g] * chi_xxz_val[mu_g];

                    double gdenxyxa = F_a_x_val[mu_g] * chi_xy_val[mu_g] + F_a_val[mu_g] * chi_xxy_val[mu_g];
                    double gdenxyya = F_a_x_val[mu_g] * chi_yy_val[mu_g] + F_a_val[mu_g] * chi_xyy_val[mu_g];
                    double gdenxyza = F_a_x_val[mu_g] * chi_yz_val[mu_g] + F_a_val[mu_g] * chi_xyz_val[mu_g];

                    double gdenxzxa = F_a_x_val[mu_g] * chi_xz_val[mu_g] + F_a_val[mu_g] * chi_xxz_val[mu_g];
                    double gdenxzya = F_a_x_val[mu_g] * chi_yz_val[mu_g] + F_a_val[mu_g] * chi_xyz_val[mu_g];
                    double gdenxzza = F_a_x_val[mu_g] * chi_zz_val[mu_g] + F_a_val[mu_g] * chi_xzz_val[mu_g];

                    double gdenxxxb = F_b_x_val[mu_g] * chi_xx_val[mu_g] + F_b_val[mu_g] * chi_xxx_val[mu_g];
                    double gdenxxyb = F_b_x_val[mu_g] * chi_xy_val[mu_g] + F_b_val[mu_g] * chi_xxy_val[mu_g];
                    double gdenxxzb = F_b_x_val[mu_g] * chi_xz_val[mu_g] + F_b_val[mu_g] * chi_xxz_val[mu_g];

                    double gdenxyxb = F_b_x_val[mu_g] * chi_xy_val[mu_g] + F_b_val[mu_g] * chi_xxy_val[mu_g];
                    double gdenxyyb = F_b_x_val[mu_g] * chi_yy_val[mu_g] + F_b_val[mu_g] * chi_xyy_val[mu_g];
                    double gdenxyzb = F_b_x_val[mu_g] * chi_yz_val[mu_g] + F_b_val[mu_g] * chi_xyz_val[mu_g];

                    double gdenxzxb = F_b_x_val[mu_g] * chi_xz_val[mu_g] + F_b_val[mu_g] * chi_xxz_val[mu_g];
                    double gdenxzyb = F_b_x_val[mu_g] * chi_yz_val[mu_g] + F_b_val[mu_g] * chi_xyz_val[mu_g];
                    double gdenxzzb = F_b_x_val[mu_g] * chi_zz_val[mu_g] + F_b_val[mu_g] * chi_xzz_val[mu_g];

                    // === y ===

                    double gdenyxxa = F_a_y_val[mu_g] * chi_xx_val[mu_g] + F_a_val[mu_g] * chi_xxy_val[mu_g];
                    double gdenyxya = F_a_y_val[mu_g] * chi_xy_val[mu_g] + F_a_val[mu_g] * chi_xyy_val[mu_g];
                    double gdenyxza = F_a_y_val[mu_g] * chi_xz_val[mu_g] + F_a_val[mu_g] * chi_xyz_val[mu_g];

                    double gdenyyxa = F_a_y_val[mu_g] * chi_xy_val[mu_g] + F_a_val[mu_g] * chi_xyy_val[mu_g];
                    double gdenyyya = F_a_y_val[mu_g] * chi_yy_val[mu_g] + F_a_val[mu_g] * chi_yyy_val[mu_g];
                    double gdenyyza = F_a_y_val[mu_g] * chi_yz_val[mu_g] + F_a_val[mu_g] * chi_yyz_val[mu_g];

                    double gdenyzxa = F_a_y_val[mu_g] * chi_xz_val[mu_g] + F_a_val[mu_g] * chi_xyz_val[mu_g];
                    double gdenyzya = F_a_y_val[mu_g] * chi_yz_val[mu_g] + F_a_val[mu_g] * chi_yyz_val[mu_g];
                    double gdenyzza = F_a_y_val[mu_g] * chi_zz_val[mu_g] + F_a_val[mu_g] * chi_yzz_val[mu_g];

                    double gdenyxxb = F_b_y_val[mu_g] * chi_xx_val[mu_g] + F_b_val[mu_g] * chi_xxy_val[mu_g];
                    double gdenyxyb = F_b_y_val[mu_g] * chi_xy_val[mu_g] + F_b_val[mu_g] * chi_xyy_val[mu_g];
                    double gdenyxzb = F_b_y_val[mu_g] * chi_xz_val[mu_g] + F_b_val[mu_g] * chi_xyz_val[mu_g];

                    double gdenyyxb = F_b_y_val[mu_g] * chi_xy_val[mu_g] + F_b_val[mu_g] * chi_xyy_val[mu_g];
                    double gdenyyyb = F_b_y_val[mu_g] * chi_yy_val[mu_g] + F_b_val[mu_g] * chi_yyy_val[mu_g];
                    double gdenyyzb = F_b_y_val[mu_g] * chi_yz_val[mu_g] + F_b_val[mu_g] * chi_yyz_val[mu_g];

                    double gdenyzxb = F_b_y_val[mu_g] * chi_xz_val[mu_g] + F_b_val[mu_g] * chi_xyz_val[mu_g];
                    double gdenyzyb = F_b_y_val[mu_g] * chi_yz_val[mu_g] + F_b_val[mu_g] * chi_yyz_val[mu_g];
                    double gdenyzzb = F_b_y_val[mu_g] * chi_zz_val[mu_g] + F_b_val[mu_g] * chi_yzz_val[mu_g];

                    // === z ===

                    double gdenzxxa = F_a_z_val[mu_g] * chi_xx_val[mu_g] + F_a_val[mu_g] * chi_xxz_val[mu_g];
                    double gdenzxya = F_a_z_val[mu_g] * chi_xy_val[mu_g] + F_a_val[mu_g] * chi_xyz_val[mu_g];
                    double gdenzxza = F_a_z_val[mu_g] * chi_xz_val[mu_g] + F_a_val[mu_g] * chi_xzz_val[mu_g];

                    double gdenzyxa = F_a_z_val[mu_g] * chi_xy_val[mu_g] + F_a_val[mu_g] * chi_xyz_val[mu_g];
                    double gdenzyya = F_a_z_val[mu_g] * chi_yy_val[mu_g] + F_a_val[mu_g] * chi_yyz_val[mu_g];
                    double gdenzyza = F_a_z_val[mu_g] * chi_yz_val[mu_g] + F_a_val[mu_g] * chi_yzz_val[mu_g];

                    double gdenzzxa = F_a_z_val[mu_g] * chi_xz_val[mu_g] + F_a_val[mu_g] * chi_xzz_val[mu_g];
                    double gdenzzya = F_a_z_val[mu_g] * chi_yz_val[mu_g] + F_a_val[mu_g] * chi_yzz_val[mu_g];
                    double gdenzzza = F_a_z_val[mu_g] * chi_zz_val[mu_g] + F_a_val[mu_g] * chi_zzz_val[mu_g];

                    double gdenzxxb = F_b_z_val[mu_g] * chi_xx_val[mu_g] + F_b_val[mu_g] * chi_xxz_val[mu_g];
                    double gdenzxyb = F_b_z_val[mu_g] * chi_xy_val[mu_g] + F_b_val[mu_g] * chi_xyz_val[mu_g];
                    double gdenzxzb = F_b_z_val[mu_g] * chi_xz_val[mu_g] + F_b_val[mu_g] * chi_xzz_val[mu_g];

                    double gdenzyxb = F_b_z_val[mu_g] * chi_xy_val[mu_g] + F_b_val[mu_g] * chi_xyz_val[mu_g];
                    double gdenzyyb = F_b_z_val[mu_g] * chi_yy_val[mu_g] + F_b_val[mu_g] * chi_yyz_val[mu_g];
                    double gdenzyzb = F_b_z_val[mu_g] * chi_yz_val[mu_g] + F_b_val[mu_g] * chi_yzz_val[mu_g];

                    double gdenzzxb = F_b_z_val[mu_g] * chi_xz_val[mu_g] + F_b_val[mu_g] * chi_xzz_val[mu_g];
                    double gdenzzyb = F_b_z_val[mu_g] * chi_yz_val[mu_g] + F_b_val[mu_g] * chi_yzz_val[mu_g];
                    double gdenzzzb = F_b_z_val[mu_g] * chi_zz_val[mu_g] + F_b_val[mu_g] * chi_zzz_val[mu_g];

                    // accumulate contribution

                    // (exc_rho[a]) * rho[a]_xi_zeta
                    // (exc_rho[b]) * rho[b]_xi_zeta

                    double f_a = vrho[2 * g + 0];
                    double f_b = vrho[2 * g + 1];

                    gatmxx += w * (f_a * gdenxxa + f_b * gdenxxb);
                    gatmxy += w * (f_a * gdenxya + f_b * gdenxyb);
                    gatmxz += w * (f_a * gdenxza + f_b * gdenxzb);

                    gatmyx += w * (f_a * gdenyxa + f_b * gdenyxb);
                    gatmyy += w * (f_a * gdenyya + f_b * gdenyyb);
                    gatmyz += w * (f_a * gdenyza + f_b * gdenyzb);

                    gatmzx += w * (f_a * gdenzxa + f_b * gdenzxb);
                    gatmzy += w * (f_a * gdenzya + f_b * gdenzyb);
                    gatmzz += w * (f_a * gdenzza + f_b * gdenzzb);

                    // (2 exc_sigma[aa]) * nabla[0]_rho[a] * nabla[0]_rho[a]_xi_zeta
                    // (  exc_sigma[ab]) * nabla[0]_rho[b] * nabla[0]_rho[a]_xi_zeta
                    // (  exc_sigma[ab]) * nabla[0]_rho[a] * nabla[0]_rho[b]_xi_zeta
                    // (2 exc_sigma[bb]) * nabla[0]_rho[b] * nabla[0]_rho[b]_xi_zeta

                    double f_aa = vsigma[3 * g + 0];
                    double f_ab = vsigma[3 * g + 1];
                    double f_bb = vsigma[3 * g + 2];

                    double g_a_x = rhograd[6 * g + 0];
                    double g_a_y = rhograd[6 * g + 1];
                    double g_a_z = rhograd[6 * g + 2];

                    double g_b_x = rhograd[6 * g + 3];
                    double g_b_y = rhograd[6 * g + 4];
                    double g_b_z = rhograd[6 * g + 5];

                    // (2 exc_sigma[aa]) * nabla[0]_rho[a] * nabla[0]_rho[a]_xi_zeta
                    // (  exc_sigma[ab]) * nabla[0]_rho[b] * nabla[0]_rho[a]_xi_zeta
                    // (  exc_sigma[ab]) * nabla[0]_rho[a] * nabla[0]_rho[b]_xi_zeta
                    // (2 exc_sigma[bb]) * nabla[0]_rho[b] * nabla[0]_rho[b]_xi_zeta

                    // xx,xy,xz

                    gatmxx += w * f_aa * 2.0 * (g_a_x * gdenxxxa + g_a_y * gdenyxxa + g_a_z * gdenzxxa);
                    gatmxy += w * f_aa * 2.0 * (g_a_x * gdenxxya + g_a_y * gdenyxya + g_a_z * gdenzxya);
                    gatmxz += w * f_aa * 2.0 * (g_a_x * gdenxxza + g_a_y * gdenyxza + g_a_z * gdenzxza);

                    gatmxx += w * f_ab * (g_b_x * gdenxxxa + g_b_y * gdenyxxa + g_b_z * gdenzxxa);
                    gatmxy += w * f_ab * (g_b_x * gdenxxya + g_b_y * gdenyxya + g_b_z * gdenzxya);
                    gatmxz += w * f_ab * (g_b_x * gdenxxza + g_b_y * gdenyxza + g_b_z * gdenzxza);

                    gatmxx += w * f_ab * (g_a_x * gdenxxxb + g_a_y * gdenyxxb + g_a_z * gdenzxxb);
                    gatmxy += w * f_ab * (g_a_x * gdenxxyb + g_a_y * gdenyxyb + g_a_z * gdenzxyb);
                    gatmxz += w * f_ab * (g_a_x * gdenxxzb + g_a_y * gdenyxzb + g_a_z * gdenzxzb);

                    gatmxx += w * f_bb * 2.0 * (g_b_x * gdenxxxb + g_b_y * gdenyxxb + g_b_z * gdenzxxb);
                    gatmxy += w * f_bb * 2.0 * (g_b_x * gdenxxyb + g_b_y * gdenyxyb + g_b_z * gdenzxyb);
                    gatmxz += w * f_bb * 2.0 * (g_b_x * gdenxxzb + g_b_y * gdenyxzb + g_b_z * gdenzxzb);

                    // yx,yy,yz

                    gatmyx += w * f_aa * 2.0 * (g_a_x * gdenxyxa + g_a_y * gdenyyxa + g_a_z * gdenzyxa);
                    gatmyy += w * f_aa * 2.0 * (g_a_x * gdenxyya + g_a_y * gdenyyya + g_a_z * gdenzyya);
                    gatmyz += w * f_aa * 2.0 * (g_a_x * gdenxyza + g_a_y * gdenyyza + g_a_z * gdenzyza);

                    gatmyx += w * f_ab * (g_b_x * gdenxyxa + g_b_y * gdenyyxa + g_b_z * gdenzyxa);
                    gatmyy += w * f_ab * (g_b_x * gdenxyya + g_b_y * gdenyyya + g_b_z * gdenzyya);
                    gatmyz += w * f_ab * (g_b_x * gdenxyza + g_b_y * gdenyyza + g_b_z * gdenzyza);

                    gatmyx += w * f_ab * (g_a_x * gdenxyxb + g_a_y * gdenyyxb + g_a_z * gdenzyxb);
                    gatmyy += w * f_ab * (g_a_x * gdenxyyb + g_a_y * gdenyyyb + g_a_z * gdenzyyb);
                    gatmyz += w * f_ab * (g_a_x * gdenxyzb + g_a_y * gdenyyzb + g_a_z * gdenzyzb);

                    gatmyx += w * f_bb * 2.0 * (g_b_x * gdenxyxb + g_b_y * gdenyyxb + g_b_z * gdenzyxb);
                    gatmyy += w * f_bb * 2.0 * (g_b_x * gdenxyyb + g_b_y * gdenyyyb + g_b_z * gdenzyyb);
                    gatmyz += w * f_bb * 2.0 * (g_b_x * gdenxyzb + g_b_y * gdenyyzb + g_b_z * gdenzyzb);

                    // zx,zy,zz

                    gatmzx += w * f_aa * 2.0 * (g_a_x * gdenxzxa + g_a_y * gdenyzxa + g_a_z * gdenzzxa);
                    gatmzy += w * f_aa * 2.0 * (g_a_x * gdenxzya + g_a_y * gdenyzya + g_a_z * gdenzzya);
                    gatmzz += w * f_aa * 2.0 * (g_a_x * gdenxzza + g_a_y * gdenyzza + g_a_z * gdenzzza);

                    gatmzx += w * f_ab * (g_b_x * gdenxzxa + g_b_y * gdenyzxa + g_b_z * gdenzzxa);
                    gatmzy += w * f_ab * (g_b_x * gdenxzya + g_b_y * gdenyzya + g_b_z * gdenzzya);
                    gatmzz += w * f_ab * (g_b_x * gdenxzza + g_b_y * gdenyzza + g_b_z * gdenzzza);

                    gatmzx += w * f_ab * (g_a_x * gdenxzxb + g_a_y * gdenyzxb + g_a_z * gdenzzxb);
                    gatmzy += w * f_ab * (g_a_x * gdenxzyb + g_a_y * gdenyzyb + g_a_z * gdenzzyb);
                    gatmzz += w * f_ab * (g_a_x * gdenxzzb + g_a_y * gdenyzzb + g_a_z * gdenzzzb);

                    gatmzx += w * f_bb * 2.0 * (g_b_x * gdenxzxb + g_b_y * gdenyzxb + g_b_z * gdenzzxb);
                    gatmzy += w * f_bb * 2.0 * (g_b_x * gdenxzyb + g_b_y * gdenyzyb + g_b_z * gdenzzyb);
                    gatmzz += w * f_bb * 2.0 * (g_b_x * gdenxzzb + g_b_y * gdenyzzb + g_b_z * gdenzzzb);
                }

                // factor of 2 from differentiation
                // no factor of 2 for open-shell

                gatm[ix * (natoms * 3) + ix] += 2.0 * gatmxx;
                gatm[ix * (natoms * 3) + iy] += 2.0 * gatmxy;
                gatm[ix * (natoms * 3) + iz] += 2.0 * gatmxz;

                gatm[iy * (natoms * 3) + ix] += 2.0 * gatmyx;
                gatm[iy * (natoms * 3) + iy] += 2.0 * gatmyy;
                gatm[iy * (natoms * 3) + iz] += 2.0 * gatmyz;

                gatm[iz * (natoms * 3) + ix] += 2.0 * gatmzx;
                gatm[iz * (natoms * 3) + iy] += 2.0 * gatmzy;
                gatm[iz * (natoms * 3) + iz] += 2.0 * gatmzz;
            }

            // second contribution to
            // f_{\rho_{\alpha}} \rho_{\alpha}^{(\xi,\zeta)}
            // and
            // (2 f_{\sigma_{\alpha\alpha}} + f_{\sigma_{\alpha\beta}})
            // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\xi,\zeta)}
            // on the same atom and on different atoms

            for (int mu = 0; mu < aocount; mu++)
            {
                auto iatom = ao_to_atom_ids[aoinds[mu]];

                auto ix = iatom * 3 + 0;
                auto iy = iatom * 3 + 1;
                auto iz = iatom * 3 + 2;

                auto mu_offset = mu * grid_batch_size;

                for (int nu = 0; nu < aocount; nu++)
                {
                    auto jatom = ao_to_atom_ids[aoinds[nu]];

                    // only consider the upper triangular part, i.e. iatom <= jatom

                    if (iatom > jatom) continue;

                    auto jx = jatom * 3 + 0;
                    auto jy = jatom * 3 + 1;
                    auto jz = jatom * 3 + 2;

                    auto nu_offset = nu * grid_batch_size;

                    auto D_a_mn = D_a_val[mu * aocount + nu];
                    auto D_b_mn = D_b_val[mu * aocount + nu];

                    double gatmxx = 0.0, gatmxy = 0.0, gatmxz = 0.0;
                    double gatmyx = 0.0, gatmyy = 0.0, gatmyz = 0.0;
                    double gatmzx = 0.0, gatmzy = 0.0, gatmzz = 0.0;

                    #pragma omp simd reduction(+ : gatmxx, gatmxy, gatmxz, gatmyx, gatmyy, gatmyz, gatmzx, gatmzy, gatmzz) 
                    for (int g = 0; g < grid_batch_size; g++)
                    {
                        auto mu_g = mu_offset + g;
                        auto nu_g = nu_offset + g;

                        double w = local_weights[g];

                        // f_{\rho_{\alpha}} \rho_{\alpha}^{(\xi,\zeta)}

                        // \rho_{\alpha}^{(\xi,\zeta)} (second contrib.)
                        // = 2 \sum_{\mu\nu} D^{\alpha}_{\mu\nu,\rm{sym}} \phi_{\mu}^{(\xi)} \phi_{\nu}^{(\zeta)}

                        // factor of 2 and D_{\mu\nu,sym}^{\alpha} are added outside of the for loop

                        double gxx = chi_x_val[mu_g] * chi_x_val[nu_g];
                        double gxy = chi_x_val[mu_g] * chi_y_val[nu_g];
                        double gxz = chi_x_val[mu_g] * chi_z_val[nu_g];

                        double gyx = chi_y_val[mu_g] * chi_x_val[nu_g];
                        double gyy = chi_y_val[mu_g] * chi_y_val[nu_g];
                        double gyz = chi_y_val[mu_g] * chi_z_val[nu_g];

                        double gzx = chi_z_val[mu_g] * chi_x_val[nu_g];
                        double gzy = chi_z_val[mu_g] * chi_y_val[nu_g];
                        double gzz = chi_z_val[mu_g] * chi_z_val[nu_g];

                        // (2 f_{\sigma_{\alpha\alpha}} + f_{\sigma_{\alpha\beta}})
                        // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\xi,\zeta)}

                        // \nabla\rho_{\alpha}^{(\xi,\zeta)} (second contrib.)
                        // = 2 \sum_{\mu\nu} D^{\alpha}_{\mu\nu,\rm{sym}}
                        //   (\nabla\phi_{\mu}^{(\xi)} \phi_{\nu}^{(\zeta)} +
                        //    \nabla\phi_{\nu}^{(\zeta)} \phi_{\mu}^{(\xi)})

                        // ordering of components: nabla, xi, zeta

                        // === x ===

                        double gxxx = chi_xx_val[mu_g] * chi_x_val[nu_g] + chi_xx_val[nu_g] * chi_x_val[mu_g];
                        double gxxy = chi_xx_val[mu_g] * chi_y_val[nu_g] + chi_xy_val[nu_g] * chi_x_val[mu_g];
                        double gxxz = chi_xx_val[mu_g] * chi_z_val[nu_g] + chi_xz_val[nu_g] * chi_x_val[mu_g];

                        double gxyx = chi_xy_val[mu_g] * chi_x_val[nu_g] + chi_xx_val[nu_g] * chi_y_val[mu_g];
                        double gxyy = chi_xy_val[mu_g] * chi_y_val[nu_g] + chi_xy_val[nu_g] * chi_y_val[mu_g];
                        double gxyz = chi_xy_val[mu_g] * chi_z_val[nu_g] + chi_xz_val[nu_g] * chi_y_val[mu_g];

                        double gxzx = chi_xz_val[mu_g] * chi_x_val[nu_g] + chi_xx_val[nu_g] * chi_z_val[mu_g];
                        double gxzy = chi_xz_val[mu_g] * chi_y_val[nu_g] + chi_xy_val[nu_g] * chi_z_val[mu_g];
                        double gxzz = chi_xz_val[mu_g] * chi_z_val[nu_g] + chi_xz_val[nu_g] * chi_z_val[mu_g];

                        // === y ===

                        double gyxx = chi_xy_val[mu_g] * chi_x_val[nu_g] + chi_xy_val[nu_g] * chi_x_val[mu_g];
                        double gyxy = chi_xy_val[mu_g] * chi_y_val[nu_g] + chi_yy_val[nu_g] * chi_x_val[mu_g];
                        double gyxz = chi_xy_val[mu_g] * chi_z_val[nu_g] + chi_yz_val[nu_g] * chi_x_val[mu_g];

                        double gyyx = chi_yy_val[mu_g] * chi_x_val[nu_g] + chi_xy_val[nu_g] * chi_y_val[mu_g];
                        double gyyy = chi_yy_val[mu_g] * chi_y_val[nu_g] + chi_yy_val[nu_g] * chi_y_val[mu_g];
                        double gyyz = chi_yy_val[mu_g] * chi_z_val[nu_g] + chi_yz_val[nu_g] * chi_y_val[mu_g];

                        double gyzx = chi_yz_val[mu_g] * chi_x_val[nu_g] + chi_xy_val[nu_g] * chi_z_val[mu_g];
                        double gyzy = chi_yz_val[mu_g] * chi_y_val[nu_g] + chi_yy_val[nu_g] * chi_z_val[mu_g];
                        double gyzz = chi_yz_val[mu_g] * chi_z_val[nu_g] + chi_yz_val[nu_g] * chi_z_val[mu_g];

                        // === z ===

                        double gzxx = chi_xz_val[mu_g] * chi_x_val[nu_g] + chi_xz_val[nu_g] * chi_x_val[mu_g];
                        double gzxy = chi_xz_val[mu_g] * chi_y_val[nu_g] + chi_yz_val[nu_g] * chi_x_val[mu_g];
                        double gzxz = chi_xz_val[mu_g] * chi_z_val[nu_g] + chi_zz_val[nu_g] * chi_x_val[mu_g];

                        double gzyx = chi_yz_val[mu_g] * chi_x_val[nu_g] + chi_xz_val[nu_g] * chi_y_val[mu_g];
                        double gzyy = chi_yz_val[mu_g] * chi_y_val[nu_g] + chi_yz_val[nu_g] * chi_y_val[mu_g];
                        double gzyz = chi_yz_val[mu_g] * chi_z_val[nu_g] + chi_zz_val[nu_g] * chi_y_val[mu_g];

                        double gzzx = chi_zz_val[mu_g] * chi_x_val[nu_g] + chi_xz_val[nu_g] * chi_z_val[mu_g];
                        double gzzy = chi_zz_val[mu_g] * chi_y_val[nu_g] + chi_yz_val[nu_g] * chi_z_val[mu_g];
                        double gzzz = chi_zz_val[mu_g] * chi_z_val[nu_g] + chi_zz_val[nu_g] * chi_z_val[mu_g];

                        // accumulate contributions

                        // (exc_rho[a]) * rho[a]_xi_zeta
                        // (exc_rho[b]) * rho[b]_xi_zeta

                        double f_a = vrho[2 * g + 0];
                        double f_b = vrho[2 * g + 1];

                        gatmxx += w * (f_a * D_a_mn + f_b * D_b_mn) * gxx;
                        gatmxy += w * (f_a * D_a_mn + f_b * D_b_mn) * gxy;
                        gatmxz += w * (f_a * D_a_mn + f_b * D_b_mn) * gxz;

                        gatmyx += w * (f_a * D_a_mn + f_b * D_b_mn) * gyx;
                        gatmyy += w * (f_a * D_a_mn + f_b * D_b_mn) * gyy;
                        gatmyz += w * (f_a * D_a_mn + f_b * D_b_mn) * gyz;

                        gatmzx += w * (f_a * D_a_mn + f_b * D_b_mn) * gzx;
                        gatmzy += w * (f_a * D_a_mn + f_b * D_b_mn) * gzy;
                        gatmzz += w * (f_a * D_a_mn + f_b * D_b_mn) * gzz;

                        // (2 exc_sigma[aa]) * nabla[0]_rho[a] * nabla[0]_rho[a]_xi_zeta
                        // (  exc_sigma[ab]) * nabla[0]_rho[b] * nabla[0]_rho[a]_xi_zeta
                        // (  exc_sigma[ab]) * nabla[0]_rho[a] * nabla[0]_rho[b]_xi_zeta
                        // (2 exc_sigma[bb]) * nabla[0]_rho[b] * nabla[0]_rho[b]_xi_zeta

                        double f_aa = vsigma[3 * g + 0];
                        double f_ab = vsigma[3 * g + 1];
                        double f_bb = vsigma[3 * g + 2];

                        double g_a_x = rhograd[6 * g + 0];
                        double g_a_y = rhograd[6 * g + 1];
                        double g_a_z = rhograd[6 * g + 2];

                        double g_b_x = rhograd[6 * g + 3];
                        double g_b_y = rhograd[6 * g + 4];
                        double g_b_z = rhograd[6 * g + 5];

                        // xx,xy,xz

                        gatmxx += w * f_aa * 2.0 * (g_a_x * gxxx + g_a_y * gyxx + g_a_z * gzxx) * D_a_mn;
                        gatmxy += w * f_aa * 2.0 * (g_a_x * gxxy + g_a_y * gyxy + g_a_z * gzxy) * D_a_mn;
                        gatmxz += w * f_aa * 2.0 * (g_a_x * gxxz + g_a_y * gyxz + g_a_z * gzxz) * D_a_mn;

                        gatmxx += w * f_ab * (g_b_x * gxxx + g_b_y * gyxx + g_b_z * gzxx) * D_a_mn;
                        gatmxy += w * f_ab * (g_b_x * gxxy + g_b_y * gyxy + g_b_z * gzxy) * D_a_mn;
                        gatmxz += w * f_ab * (g_b_x * gxxz + g_b_y * gyxz + g_b_z * gzxz) * D_a_mn;

                        gatmxx += w * f_ab * (g_a_x * gxxx + g_a_y * gyxx + g_a_z * gzxx) * D_b_mn;
                        gatmxy += w * f_ab * (g_a_x * gxxy + g_a_y * gyxy + g_a_z * gzxy) * D_b_mn;
                        gatmxz += w * f_ab * (g_a_x * gxxz + g_a_y * gyxz + g_a_z * gzxz) * D_b_mn;

                        gatmxx += w * f_bb * 2.0 * (g_b_x * gxxx + g_b_y * gyxx + g_b_z * gzxx) * D_b_mn;
                        gatmxy += w * f_bb * 2.0 * (g_b_x * gxxy + g_b_y * gyxy + g_b_z * gzxy) * D_b_mn;
                        gatmxz += w * f_bb * 2.0 * (g_b_x * gxxz + g_b_y * gyxz + g_b_z * gzxz) * D_b_mn;

                        // yx,yy,yz

                        gatmyx += w * f_aa * 2.0 * (g_a_x * gxyx + g_a_y * gyyx + g_a_z * gzyx) * D_a_mn;
                        gatmyy += w * f_aa * 2.0 * (g_a_x * gxyy + g_a_y * gyyy + g_a_z * gzyy) * D_a_mn;
                        gatmyz += w * f_aa * 2.0 * (g_a_x * gxyz + g_a_y * gyyz + g_a_z * gzyz) * D_a_mn;

                        gatmyx += w * f_ab * (g_b_x * gxyx + g_b_y * gyyx + g_b_z * gzyx) * D_a_mn;
                        gatmyy += w * f_ab * (g_b_x * gxyy + g_b_y * gyyy + g_b_z * gzyy) * D_a_mn;
                        gatmyz += w * f_ab * (g_b_x * gxyz + g_b_y * gyyz + g_b_z * gzyz) * D_a_mn;

                        gatmyx += w * f_ab * (g_a_x * gxyx + g_a_y * gyyx + g_a_z * gzyx) * D_b_mn;
                        gatmyy += w * f_ab * (g_a_x * gxyy + g_a_y * gyyy + g_a_z * gzyy) * D_b_mn;
                        gatmyz += w * f_ab * (g_a_x * gxyz + g_a_y * gyyz + g_a_z * gzyz) * D_b_mn;

                        gatmyx += w * f_bb * 2.0 * (g_b_x * gxyx + g_b_y * gyyx + g_b_z * gzyx) * D_b_mn;
                        gatmyy += w * f_bb * 2.0 * (g_b_x * gxyy + g_b_y * gyyy + g_b_z * gzyy) * D_b_mn;
                        gatmyz += w * f_bb * 2.0 * (g_b_x * gxyz + g_b_y * gyyz + g_b_z * gzyz) * D_b_mn;

                        // zx,zy,zz

                        gatmzx += w * f_aa * 2.0 * (g_a_x * gxzx + g_a_y * gyzx + g_a_z * gzzx) * D_a_mn;
                        gatmzy += w * f_aa * 2.0 * (g_a_x * gxzy + g_a_y * gyzy + g_a_z * gzzy) * D_a_mn;
                        gatmzz += w * f_aa * 2.0 * (g_a_x * gxzz + g_a_y * gyzz + g_a_z * gzzz) * D_a_mn;

                        gatmzx += w * f_ab * (g_b_x * gxzx + g_b_y * gyzx + g_b_z * gzzx) * D_a_mn;
                        gatmzy += w * f_ab * (g_b_x * gxzy + g_b_y * gyzy + g_b_z * gzzy) * D_a_mn;
                        gatmzz += w * f_ab * (g_b_x * gxzz + g_b_y * gyzz + g_b_z * gzzz) * D_a_mn;

                        gatmzx += w * f_ab * (g_a_x * gxzx + g_a_y * gyzx + g_a_z * gzzx) * D_b_mn;
                        gatmzy += w * f_ab * (g_a_x * gxzy + g_a_y * gyzy + g_a_z * gzzy) * D_b_mn;
                        gatmzz += w * f_ab * (g_a_x * gxzz + g_a_y * gyzz + g_a_z * gzzz) * D_b_mn;

                        gatmzx += w * f_bb * 2.0 * (g_b_x * gxzx + g_b_y * gyzx + g_b_z * gzzx) * D_b_mn;
                        gatmzy += w * f_bb * 2.0 * (g_b_x * gxzy + g_b_y * gyzy + g_b_z * gzzy) * D_b_mn;
                        gatmzz += w * f_bb * 2.0 * (g_b_x * gxzz + g_b_y * gyzz + g_b_z * gzzz) * D_b_mn;
                    }

                    // factor of 2 from differentiation
                    // no factor of 2 for open-shell

                    gatm[ix * (natoms * 3) + jx] += 2.0 * gatmxx;
                    gatm[ix * (natoms * 3) + jy] += 2.0 * gatmxy;
                    gatm[ix * (natoms * 3) + jz] += 2.0 * gatmxz;

                    gatm[iy * (natoms * 3) + jx] += 2.0 * gatmyx;
                    gatm[iy * (natoms * 3) + jy] += 2.0 * gatmyy;
                    gatm[iy * (natoms * 3) + jz] += 2.0 * gatmyz;

                    gatm[iz * (natoms * 3) + jx] += 2.0 * gatmzx;
                    gatm[iz * (natoms * 3) + jy] += 2.0 * gatmzy;
                    gatm[iz * (natoms * 3) + jz] += 2.0 * gatmzz;
                }
            }

            // other contributions

            for (int iatom = 0; iatom < natoms; iatom++)
            {
                auto ix = iatom * 3 + 0;
                auto iy = iatom * 3 + 1;
                auto iz = iatom * 3 + 2;

                auto i_offset = iatom * grid_batch_size;

                for (int jatom = iatom; jatom < natoms; jatom++)
                {
                    auto jx = jatom * 3 + 0;
                    auto jy = jatom * 3 + 1;
                    auto jz = jatom * 3 + 2;

                    auto j_offset = jatom * grid_batch_size;

                    double gatmxx = 0.0, gatmyx = 0.0, gatmzx = 0.0;
                    double gatmxy = 0.0, gatmyy = 0.0, gatmzy = 0.0;
                    double gatmxz = 0.0, gatmyz = 0.0, gatmzz = 0.0;

                    #pragma omp simd reduction(+ : gatmxx, gatmxy, gatmxz, gatmyx, gatmyy, gatmyz, gatmzx, gatmzy, gatmzz) 
                    for (int g = 0; g < grid_batch_size; g++)
                    {
                        auto ig = i_offset + g;
                        auto jg = j_offset + g;

                        double w = local_weights[g];

                        // (f_{\rho_{\alpha} \rho_{\alpha}} + f_{\rho_{\alpha} \rho_{\beta}})
                        // \rho_{\alpha}^{(\xi)} \rho_{\alpha}^{(\zeta)}

                        // (exc_rho[a]_rho[a]) * rho[a]_xi * rho[a]_zeta
                        // (exc_rho[a]_rho[b]) * rho[a]_xi * rho[b]_zeta
                        // (exc_rho[a]_rho[b]) * rho[b]_xi * rho[a]_zeta
                        // (exc_rho[b]_rho[b]) * rho[b]_xi * rho[b]_zeta

                        double prefac_aa = w * v2rho2[3 * g + 0];
                        double prefac_ab = w * v2rho2[3 * g + 1];
                        double prefac_bb = w * v2rho2[3 * g + 2];

                        // xx,xy,xz

                        gatmxx += prefac_aa * gden_a_x[ig] * gden_a_x[jg];
                        gatmxy += prefac_aa * gden_a_x[ig] * gden_a_y[jg];
                        gatmxz += prefac_aa * gden_a_x[ig] * gden_a_z[jg];

                        gatmxx += prefac_ab * gden_a_x[ig] * gden_b_x[jg];
                        gatmxy += prefac_ab * gden_a_x[ig] * gden_b_y[jg];
                        gatmxz += prefac_ab * gden_a_x[ig] * gden_b_z[jg];

                        gatmxx += prefac_ab * gden_b_x[ig] * gden_a_x[jg];
                        gatmxy += prefac_ab * gden_b_x[ig] * gden_a_y[jg];
                        gatmxz += prefac_ab * gden_b_x[ig] * gden_a_z[jg];

                        gatmxx += prefac_bb * gden_b_x[ig] * gden_b_x[jg];
                        gatmxy += prefac_bb * gden_b_x[ig] * gden_b_y[jg];
                        gatmxz += prefac_bb * gden_b_x[ig] * gden_b_z[jg];

                        // yx,yy,yz

                        gatmyx += prefac_aa * gden_a_y[ig] * gden_a_x[jg];
                        gatmyy += prefac_aa * gden_a_y[ig] * gden_a_y[jg];
                        gatmyz += prefac_aa * gden_a_y[ig] * gden_a_z[jg];

                        gatmyx += prefac_ab * gden_a_y[ig] * gden_b_x[jg];
                        gatmyy += prefac_ab * gden_a_y[ig] * gden_b_y[jg];
                        gatmyz += prefac_ab * gden_a_y[ig] * gden_b_z[jg];

                        gatmyx += prefac_ab * gden_b_y[ig] * gden_a_x[jg];
                        gatmyy += prefac_ab * gden_b_y[ig] * gden_a_y[jg];
                        gatmyz += prefac_ab * gden_b_y[ig] * gden_a_z[jg];

                        gatmyx += prefac_bb * gden_b_y[ig] * gden_b_x[jg];
                        gatmyy += prefac_bb * gden_b_y[ig] * gden_b_y[jg];
                        gatmyz += prefac_bb * gden_b_y[ig] * gden_b_z[jg];

                        // zx,zy,zz

                        gatmzx += prefac_aa * gden_a_z[ig] * gden_a_x[jg];
                        gatmzy += prefac_aa * gden_a_z[ig] * gden_a_y[jg];
                        gatmzz += prefac_aa * gden_a_z[ig] * gden_a_z[jg];

                        gatmzx += prefac_ab * gden_a_z[ig] * gden_b_x[jg];
                        gatmzy += prefac_ab * gden_a_z[ig] * gden_b_y[jg];
                        gatmzz += prefac_ab * gden_a_z[ig] * gden_b_z[jg];

                        gatmzx += prefac_ab * gden_b_z[ig] * gden_a_x[jg];
                        gatmzy += prefac_ab * gden_b_z[ig] * gden_a_y[jg];
                        gatmzz += prefac_ab * gden_b_z[ig] * gden_a_z[jg];

                        gatmzx += prefac_bb * gden_b_z[ig] * gden_b_x[jg];
                        gatmzy += prefac_bb * gden_b_z[ig] * gden_b_y[jg];
                        gatmzz += prefac_bb * gden_b_z[ig] * gden_b_z[jg];

                        // (2 f_{\rho_{\alpha} \sigma_{\alpha\alpha}} +
                        //  2 f_{\rho_{\alpha} \sigma_{\alpha\beta}} +
                        //  2 f_{\rho_{\alpha} \sigma_{\beta\beta}})
                        // \rho_{\alpha}^{(\xi)} \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\zeta)}

                        // (2 exc_rho[a]_sigma[aa]) * nabla[1]_rho[a] * nabla[1]_rho[a]_zeta * rho[a]_xi
                        // (  exc_rho[a]_sigma[ab]) * nabla[1]_rho[b] * nabla[1]_rho[a]_zeta * rho[a]_xi
                        // (  exc_rho[a]_sigma[ab]) * nabla[1]_rho[a] * nabla[1]_rho[b]_zeta * rho[a]_xi
                        // (2 exc_rho[a]_sigma[bb]) * nabla[1]_rho[b] * nabla[1]_rho[b]_zeta * rho[a]_xi
                        // (2 exc_rho[b]_sigma[aa]) * nabla[1]_rho[a] * nabla[1]_rho[a]_zeta * rho[b]_xi
                        // (  exc_rho[b]_sigma[ab]) * nabla[1]_rho[b] * nabla[1]_rho[a]_zeta * rho[b]_xi
                        // (  exc_rho[b]_sigma[ab]) * nabla[1]_rho[a] * nabla[1]_rho[b]_zeta * rho[b]_xi
                        // (2 exc_rho[b]_sigma[bb]) * nabla[1]_rho[b] * nabla[1]_rho[b]_zeta * rho[b]_xi

                        auto f_a_aa = v2rhosigma[6 * g + 0];
                        auto f_a_ab = v2rhosigma[6 * g + 1];
                        auto f_a_bb = v2rhosigma[6 * g + 2];
                        auto f_b_aa = v2rhosigma[6 * g + 3];
                        auto f_b_ab = v2rhosigma[6 * g + 4];
                        auto f_b_bb = v2rhosigma[6 * g + 5];

                        auto g_a_x = rhograd[6 * g + 0];
                        auto g_a_y = rhograd[6 * g + 1];
                        auto g_a_z = rhograd[6 * g + 2];

                        auto g_b_x = rhograd[6 * g + 3];
                        auto g_b_y = rhograd[6 * g + 4];
                        auto g_b_z = rhograd[6 * g + 5];

                        auto xcomp_j_aa = (g_a_x * gden_a_xx[jg] + g_a_y * gden_a_yx[jg] + g_a_z * gden_a_zx[jg]);
                        auto ycomp_j_aa = (g_a_x * gden_a_xy[jg] + g_a_y * gden_a_yy[jg] + g_a_z * gden_a_zy[jg]);
                        auto zcomp_j_aa = (g_a_x * gden_a_xz[jg] + g_a_y * gden_a_yz[jg] + g_a_z * gden_a_zz[jg]);

                        auto xcomp_j_ba = (g_b_x * gden_a_xx[jg] + g_b_y * gden_a_yx[jg] + g_b_z * gden_a_zx[jg]);
                        auto ycomp_j_ba = (g_b_x * gden_a_xy[jg] + g_b_y * gden_a_yy[jg] + g_b_z * gden_a_zy[jg]);
                        auto zcomp_j_ba = (g_b_x * gden_a_xz[jg] + g_b_y * gden_a_yz[jg] + g_b_z * gden_a_zz[jg]);

                        auto xcomp_j_ab = (g_a_x * gden_b_xx[jg] + g_a_y * gden_b_yx[jg] + g_a_z * gden_b_zx[jg]);
                        auto ycomp_j_ab = (g_a_x * gden_b_xy[jg] + g_a_y * gden_b_yy[jg] + g_a_z * gden_b_zy[jg]);
                        auto zcomp_j_ab = (g_a_x * gden_b_xz[jg] + g_a_y * gden_b_yz[jg] + g_a_z * gden_b_zz[jg]);

                        auto xcomp_j_bb = (g_b_x * gden_b_xx[jg] + g_b_y * gden_b_yx[jg] + g_b_z * gden_b_zx[jg]);
                        auto ycomp_j_bb = (g_b_x * gden_b_xy[jg] + g_b_y * gden_b_yy[jg] + g_b_z * gden_b_zy[jg]);
                        auto zcomp_j_bb = (g_b_x * gden_b_xz[jg] + g_b_y * gden_b_yz[jg] + g_b_z * gden_b_zz[jg]);

                        // (2 exc_rho[a]_sigma[aa]) * nabla[1]_rho[a] * nabla[1]_rho[a]_zeta * rho[a]_xi
                        // (  exc_rho[a]_sigma[ab]) * nabla[1]_rho[b] * nabla[1]_rho[a]_zeta * rho[a]_xi
                        // (  exc_rho[a]_sigma[ab]) * nabla[1]_rho[a] * nabla[1]_rho[b]_zeta * rho[a]_xi
                        // (2 exc_rho[a]_sigma[bb]) * nabla[1]_rho[b] * nabla[1]_rho[b]_zeta * rho[a]_xi
                        // (2 exc_rho[b]_sigma[aa]) * nabla[1]_rho[a] * nabla[1]_rho[a]_zeta * rho[b]_xi
                        // (  exc_rho[b]_sigma[ab]) * nabla[1]_rho[b] * nabla[1]_rho[a]_zeta * rho[b]_xi
                        // (  exc_rho[b]_sigma[ab]) * nabla[1]_rho[a] * nabla[1]_rho[b]_zeta * rho[b]_xi
                        // (2 exc_rho[b]_sigma[bb]) * nabla[1]_rho[b] * nabla[1]_rho[b]_zeta * rho[b]_xi

                        // xx,xy,xz

                        gatmxx += w * 2.0 * (f_a_aa * gden_a_x[ig] + f_b_aa * gden_b_x[ig]) * xcomp_j_aa;
                        gatmxy += w * 2.0 * (f_a_aa * gden_a_x[ig] + f_b_aa * gden_b_x[ig]) * ycomp_j_aa;
                        gatmxz += w * 2.0 * (f_a_aa * gden_a_x[ig] + f_b_aa * gden_b_x[ig]) * zcomp_j_aa;

                        gatmxx += w * 2.0 * (f_a_bb * gden_a_x[ig] + f_b_bb * gden_b_x[ig]) * xcomp_j_bb;
                        gatmxy += w * 2.0 * (f_a_bb * gden_a_x[ig] + f_b_bb * gden_b_x[ig]) * ycomp_j_bb;
                        gatmxz += w * 2.0 * (f_a_bb * gden_a_x[ig] + f_b_bb * gden_b_x[ig]) * zcomp_j_bb;
                                         
                        gatmxx += w * (f_a_ab * gden_a_x[ig] + f_b_ab * gden_b_x[ig]) * (xcomp_j_ba + xcomp_j_ab);
                        gatmxy += w * (f_a_ab * gden_a_x[ig] + f_b_ab * gden_b_x[ig]) * (ycomp_j_ba + ycomp_j_ab);
                        gatmxz += w * (f_a_ab * gden_a_x[ig] + f_b_ab * gden_b_x[ig]) * (zcomp_j_ba + zcomp_j_ab);

                        // yx,yy,yz

                        gatmyx += w * 2.0 * (f_a_aa * gden_a_y[ig] + f_b_aa * gden_b_y[ig]) * xcomp_j_aa;
                        gatmyy += w * 2.0 * (f_a_aa * gden_a_y[ig] + f_b_aa * gden_b_y[ig]) * ycomp_j_aa;
                        gatmyz += w * 2.0 * (f_a_aa * gden_a_y[ig] + f_b_aa * gden_b_y[ig]) * zcomp_j_aa;

                        gatmyx += w * 2.0 * (f_a_bb * gden_a_y[ig] + f_b_bb * gden_b_y[ig]) * xcomp_j_bb;
                        gatmyy += w * 2.0 * (f_a_bb * gden_a_y[ig] + f_b_bb * gden_b_y[ig]) * ycomp_j_bb;
                        gatmyz += w * 2.0 * (f_a_bb * gden_a_y[ig] + f_b_bb * gden_b_y[ig]) * zcomp_j_bb;
                                         
                        gatmyx += w * (f_a_ab * gden_a_y[ig] + f_b_ab * gden_b_y[ig]) * (xcomp_j_ba + xcomp_j_ab);
                        gatmyy += w * (f_a_ab * gden_a_y[ig] + f_b_ab * gden_b_y[ig]) * (ycomp_j_ba + ycomp_j_ab);
                        gatmyz += w * (f_a_ab * gden_a_y[ig] + f_b_ab * gden_b_y[ig]) * (zcomp_j_ba + zcomp_j_ab);

                        // zx,zy,zz

                        gatmzx += w * 2.0 * (f_a_aa * gden_a_z[ig] + f_b_aa * gden_b_z[ig]) * xcomp_j_aa;
                        gatmzy += w * 2.0 * (f_a_aa * gden_a_z[ig] + f_b_aa * gden_b_z[ig]) * ycomp_j_aa;
                        gatmzz += w * 2.0 * (f_a_aa * gden_a_z[ig] + f_b_aa * gden_b_z[ig]) * zcomp_j_aa;

                        gatmzx += w * 2.0 * (f_a_bb * gden_a_z[ig] + f_b_bb * gden_b_z[ig]) * xcomp_j_bb;
                        gatmzy += w * 2.0 * (f_a_bb * gden_a_z[ig] + f_b_bb * gden_b_z[ig]) * ycomp_j_bb;
                        gatmzz += w * 2.0 * (f_a_bb * gden_a_z[ig] + f_b_bb * gden_b_z[ig]) * zcomp_j_bb;
                                         
                        gatmzx += w * (f_a_ab * gden_a_z[ig] + f_b_ab * gden_b_z[ig]) * (xcomp_j_ba + xcomp_j_ab);
                        gatmzy += w * (f_a_ab * gden_a_z[ig] + f_b_ab * gden_b_z[ig]) * (ycomp_j_ba + ycomp_j_ab);
                        gatmzz += w * (f_a_ab * gden_a_z[ig] + f_b_ab * gden_b_z[ig]) * (zcomp_j_ba + zcomp_j_ab);

                        // (2 f_{\rho_{\alpha} \sigma_{\alpha\alpha}} +
                        //  2 f_{\rho_{\beta} \sigma_{\alpha\alpha}} +
                        //  f_{\rho_{\alpha} \sigma_{\alpha\beta}} +
                        //  f_{\rho_{\beta} \sigma_{\alpha\beta}})
                        // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\xi)} \rho_{\alpha}^{(\zeta)}

                        // (2 exc_rho[a]_sigma[aa]) * nabla[0]_rho[a] * nabla[0]_rho[a]_xi * rho[a]_zeta
                        // (2 exc_rho[b]_sigma[aa]) * nabla[0]_rho[a] * nabla[0]_rho[a]_xi * rho[b]_zeta
                        // (  exc_rho[a]_sigma[ab]) * nabla[0]_rho[b] * nabla[0]_rho[a]_xi * rho[a]_zeta
                        // (  exc_rho[b]_sigma[ab]) * nabla[0]_rho[b] * nabla[0]_rho[a]_xi * rho[b]_zeta
                        // (2 exc_rho[a]_sigma[bb]) * nabla[0]_rho[b] * nabla[0]_rho[b]_xi * rho[a]_zeta
                        // (2 exc_rho[b]_sigma[bb]) * nabla[0]_rho[b] * nabla[0]_rho[b]_xi * rho[b]_zeta
                        // (  exc_rho[a]_sigma[ab]) * nabla[0]_rho[a] * nabla[0]_rho[b]_xi * rho[a]_zeta
                        // (  exc_rho[b]_sigma[ab]) * nabla[0]_rho[a] * nabla[0]_rho[b]_xi * rho[b]_zeta

                        auto xcomp_i_aa = (g_a_x * gden_a_xx[ig] + g_a_y * gden_a_yx[ig] + g_a_z * gden_a_zx[ig]);
                        auto ycomp_i_aa = (g_a_x * gden_a_xy[ig] + g_a_y * gden_a_yy[ig] + g_a_z * gden_a_zy[ig]);
                        auto zcomp_i_aa = (g_a_x * gden_a_xz[ig] + g_a_y * gden_a_yz[ig] + g_a_z * gden_a_zz[ig]);

                        auto xcomp_i_ba = (g_b_x * gden_a_xx[ig] + g_b_y * gden_a_yx[ig] + g_b_z * gden_a_zx[ig]);
                        auto ycomp_i_ba = (g_b_x * gden_a_xy[ig] + g_b_y * gden_a_yy[ig] + g_b_z * gden_a_zy[ig]);
                        auto zcomp_i_ba = (g_b_x * gden_a_xz[ig] + g_b_y * gden_a_yz[ig] + g_b_z * gden_a_zz[ig]);

                        auto xcomp_i_ab = (g_a_x * gden_b_xx[ig] + g_a_y * gden_b_yx[ig] + g_a_z * gden_b_zx[ig]);
                        auto ycomp_i_ab = (g_a_x * gden_b_xy[ig] + g_a_y * gden_b_yy[ig] + g_a_z * gden_b_zy[ig]);
                        auto zcomp_i_ab = (g_a_x * gden_b_xz[ig] + g_a_y * gden_b_yz[ig] + g_a_z * gden_b_zz[ig]);

                        auto xcomp_i_bb = (g_b_x * gden_b_xx[ig] + g_b_y * gden_b_yx[ig] + g_b_z * gden_b_zx[ig]);
                        auto ycomp_i_bb = (g_b_x * gden_b_xy[ig] + g_b_y * gden_b_yy[ig] + g_b_z * gden_b_zy[ig]);
                        auto zcomp_i_bb = (g_b_x * gden_b_xz[ig] + g_b_y * gden_b_yz[ig] + g_b_z * gden_b_zz[ig]);

                        // (2 exc_rho[a]_sigma[aa]) * nabla[0]_rho[a] * nabla[0]_rho[a]_xi * rho[a]_zeta
                        // (2 exc_rho[b]_sigma[aa]) * nabla[0]_rho[a] * nabla[0]_rho[a]_xi * rho[b]_zeta
                        // (  exc_rho[a]_sigma[ab]) * nabla[0]_rho[b] * nabla[0]_rho[a]_xi * rho[a]_zeta
                        // (  exc_rho[b]_sigma[ab]) * nabla[0]_rho[b] * nabla[0]_rho[a]_xi * rho[b]_zeta
                        // (2 exc_rho[a]_sigma[bb]) * nabla[0]_rho[b] * nabla[0]_rho[b]_xi * rho[a]_zeta
                        // (2 exc_rho[b]_sigma[bb]) * nabla[0]_rho[b] * nabla[0]_rho[b]_xi * rho[b]_zeta
                        // (  exc_rho[a]_sigma[ab]) * nabla[0]_rho[a] * nabla[0]_rho[b]_xi * rho[a]_zeta
                        // (  exc_rho[b]_sigma[ab]) * nabla[0]_rho[a] * nabla[0]_rho[b]_xi * rho[b]_zeta

                        // xx,xy,xz

                        gatmxx += w * 2.0 * (f_a_aa * gden_a_x[jg] + f_b_aa * gden_b_x[jg]) * xcomp_i_aa;
                        gatmxy += w * 2.0 * (f_a_aa * gden_a_y[jg] + f_b_aa * gden_b_y[jg]) * xcomp_i_aa;
                        gatmxz += w * 2.0 * (f_a_aa * gden_a_z[jg] + f_b_aa * gden_b_z[jg]) * xcomp_i_aa;
                                                                              
                        gatmxx += w * (f_a_ab * gden_a_x[jg] + f_b_ab * gden_b_x[jg]) * (xcomp_i_ba + xcomp_i_ab);
                        gatmxy += w * (f_a_ab * gden_a_y[jg] + f_b_ab * gden_b_y[jg]) * (xcomp_i_ba + xcomp_i_ab);
                        gatmxz += w * (f_a_ab * gden_a_z[jg] + f_b_ab * gden_b_z[jg]) * (xcomp_i_ba + xcomp_i_ab);
                                                                        
                        gatmxx += w * 2.0 * (f_a_bb * gden_a_x[jg] + f_b_bb * gden_b_x[jg]) * xcomp_i_bb;
                        gatmxy += w * 2.0 * (f_a_bb * gden_a_y[jg] + f_b_bb * gden_b_y[jg]) * xcomp_i_bb;
                        gatmxz += w * 2.0 * (f_a_bb * gden_a_z[jg] + f_b_bb * gden_b_z[jg]) * xcomp_i_bb;

                        // yx,yy,yz

                        gatmyx += w * 2.0 * (f_a_aa * gden_a_x[jg] + f_b_aa * gden_b_x[jg]) * ycomp_i_aa;
                        gatmyy += w * 2.0 * (f_a_aa * gden_a_y[jg] + f_b_aa * gden_b_y[jg]) * ycomp_i_aa;
                        gatmyz += w * 2.0 * (f_a_aa * gden_a_z[jg] + f_b_aa * gden_b_z[jg]) * ycomp_i_aa;
                                                                              
                        gatmyx += w * (f_a_ab * gden_a_x[jg] + f_b_ab * gden_b_x[jg]) * (ycomp_i_ba + ycomp_i_ab);
                        gatmyy += w * (f_a_ab * gden_a_y[jg] + f_b_ab * gden_b_y[jg]) * (ycomp_i_ba + ycomp_i_ab);
                        gatmyz += w * (f_a_ab * gden_a_z[jg] + f_b_ab * gden_b_z[jg]) * (ycomp_i_ba + ycomp_i_ab);
                                                                        
                        gatmyx += w * 2.0 * (f_a_bb * gden_a_x[jg] + f_b_bb * gden_b_x[jg]) * ycomp_i_bb;
                        gatmyy += w * 2.0 * (f_a_bb * gden_a_y[jg] + f_b_bb * gden_b_y[jg]) * ycomp_i_bb;
                        gatmyz += w * 2.0 * (f_a_bb * gden_a_z[jg] + f_b_bb * gden_b_z[jg]) * ycomp_i_bb;

                        // zx,zy,zz

                        gatmzx += w * 2.0 * (f_a_aa * gden_a_x[jg] + f_b_aa * gden_b_x[jg]) * zcomp_i_aa;
                        gatmzy += w * 2.0 * (f_a_aa * gden_a_y[jg] + f_b_aa * gden_b_y[jg]) * zcomp_i_aa;
                        gatmzz += w * 2.0 * (f_a_aa * gden_a_z[jg] + f_b_aa * gden_b_z[jg]) * zcomp_i_aa;
                                                                              
                        gatmzx += w * (f_a_ab * gden_a_x[jg] + f_b_ab * gden_b_x[jg]) * (zcomp_i_ba + zcomp_i_ab);
                        gatmzy += w * (f_a_ab * gden_a_y[jg] + f_b_ab * gden_b_y[jg]) * (zcomp_i_ba + zcomp_i_ab);
                        gatmzz += w * (f_a_ab * gden_a_z[jg] + f_b_ab * gden_b_z[jg]) * (zcomp_i_ba + zcomp_i_ab);
                                                                        
                        gatmzx += w * 2.0 * (f_a_bb * gden_a_x[jg] + f_b_bb * gden_b_x[jg]) * zcomp_i_bb;
                        gatmzy += w * 2.0 * (f_a_bb * gden_a_y[jg] + f_b_bb * gden_b_y[jg]) * zcomp_i_bb;
                        gatmzz += w * 2.0 * (f_a_bb * gden_a_z[jg] + f_b_bb * gden_b_z[jg]) * zcomp_i_bb;

                        // (4 f_{\sigma_{\alpha\alpha} \sigma_{\alpha\alpha}} +
                        //  6 f_{\sigma_{\alpha\alpha} \sigma_{\alpha\beta}} +
                        //  4 f_{\sigma_{\alpha\alpha} \sigma_{\beta\beta}} +
                        //  2 f_{\sigma_{\alpha\beta} \sigma_{\alpha\beta}} +
                        //  2 f_{\sigma_{\alpha\beta} \sigma_{\beta\beta}})
                        // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\xi)}
                        // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\zeta)}

                        // (4 exc_sigma[aa]_sigma[aa]) * nabla[0]_rho[a] * nabla[0]_rho[a]_xi * nabla[1]_rho[a] * nabla[1]_rho[a]_zeta
                        // (2 exc_sigma[aa]_sigma[ab]) * nabla[0]_rho[a] * nabla[0]_rho[a]_xi * nabla[1]_rho[b] * nabla[1]_rho[a]_zeta
                        // (2 exc_sigma[aa]_sigma[ab]) * nabla[0]_rho[a] * nabla[0]_rho[a]_xi * nabla[1]_rho[a] * nabla[1]_rho[b]_zeta
                        // (4 exc_sigma[aa]_sigma[bb]) * nabla[0]_rho[a] * nabla[0]_rho[a]_xi * nabla[1]_rho[b] * nabla[1]_rho[b]_zeta
                        // (2 exc_sigma[aa]_sigma[ab]) * nabla[0]_rho[b] * nabla[0]_rho[a]_xi * nabla[1]_rho[a] * nabla[1]_rho[a]_zeta
                        // (  exc_sigma[ab]_sigma[ab]) * nabla[0]_rho[b] * nabla[0]_rho[a]_xi * nabla[1]_rho[b] * nabla[1]_rho[a]_zeta
                        // (  exc_sigma[ab]_sigma[ab]) * nabla[0]_rho[b] * nabla[0]_rho[a]_xi * nabla[1]_rho[a] * nabla[1]_rho[b]_zeta
                        // (2 exc_sigma[ab]_sigma[bb]) * nabla[0]_rho[b] * nabla[0]_rho[a]_xi * nabla[1]_rho[b] * nabla[1]_rho[b]_zeta
                        // (4 exc_sigma[aa]_sigma[bb]) * nabla[0]_rho[b] * nabla[0]_rho[b]_xi * nabla[1]_rho[a] * nabla[1]_rho[a]_zeta
                        // (2 exc_sigma[ab]_sigma[bb]) * nabla[0]_rho[b] * nabla[0]_rho[b]_xi * nabla[1]_rho[b] * nabla[1]_rho[a]_zeta
                        // (2 exc_sigma[ab]_sigma[bb]) * nabla[0]_rho[b] * nabla[0]_rho[b]_xi * nabla[1]_rho[a] * nabla[1]_rho[b]_zeta
                        // (4 exc_sigma[bb]_sigma[bb]) * nabla[0]_rho[b] * nabla[0]_rho[b]_xi * nabla[1]_rho[b] * nabla[1]_rho[b]_zeta
                        // (2 exc_sigma[aa]_sigma[ab]) * nabla[0]_rho[a] * nabla[0]_rho[b]_xi * nabla[1]_rho[a] * nabla[1]_rho[a]_zeta
                        // (  exc_sigma[ab]_sigma[ab]) * nabla[0]_rho[a] * nabla[0]_rho[b]_xi * nabla[1]_rho[b] * nabla[1]_rho[a]_zeta
                        // (  exc_sigma[ab]_sigma[ab]) * nabla[0]_rho[a] * nabla[0]_rho[b]_xi * nabla[1]_rho[a] * nabla[1]_rho[b]_zeta
                        // (2 exc_sigma[ab]_sigma[bb]) * nabla[0]_rho[a] * nabla[0]_rho[b]_xi * nabla[1]_rho[b] * nabla[1]_rho[b]_zeta

                        auto f_aa_aa = v2sigma2[6 * g + 0];
                        auto f_aa_ab = v2sigma2[6 * g + 1];
                        auto f_aa_bb = v2sigma2[6 * g + 2];
                        auto f_ab_ab = v2sigma2[6 * g + 3];
                        auto f_ab_bb = v2sigma2[6 * g + 4];
                        auto f_bb_bb = v2sigma2[6 * g + 5];

                        // 4 * f_aa_aa * xcomp_i_aa * xcomp_j_aa
                        // 2 * f_aa_ab * xcomp_i_aa * xcomp_j_ba
                        // 2 * f_aa_ab * xcomp_i_aa * xcomp_j_ab
                        // 4 * f_aa_bb * xcomp_i_aa * xcomp_j_bb
                        // 2 * f_aa_ab * xcomp_i_ba * xcomp_j_aa
                        //     f_ab_ab * xcomp_i_ba * xcomp_j_ba
                        //     f_ab_ab * xcomp_i_ba * xcomp_j_ab
                        // 2 * f_ab_bb * xcomp_i_ba * xcomp_j_bb
                        // 4 * f_aa_bb * xcomp_i_bb * xcomp_j_aa
                        // 2 * f_ab_bb * xcomp_i_bb * xcomp_j_ba
                        // 2 * f_ab_bb * xcomp_i_bb * xcomp_j_ab
                        // 4 * f_bb_bb * xcomp_i_bb * xcomp_j_bb
                        // 2 * f_aa_ab * xcomp_i_ab * xcomp_j_aa
                        //     f_ab_ab * xcomp_i_ab * xcomp_j_ba
                        //     f_ab_ab * xcomp_i_ab * xcomp_j_ab
                        // 2 * f_ab_bb * xcomp_i_ab * xcomp_j_bb

                        gatmxx += w * (4.0 * f_aa_aa * xcomp_i_aa * xcomp_j_aa +
                                  2.0 * f_aa_ab * xcomp_i_aa * xcomp_j_ba +
                                  2.0 * f_aa_ab * xcomp_i_aa * xcomp_j_ab +
                                  4.0 * f_aa_bb * xcomp_i_aa * xcomp_j_bb +
                                  2.0 * f_aa_ab * xcomp_i_ba * xcomp_j_aa +
                                        f_ab_ab * xcomp_i_ba * xcomp_j_ba +
                                        f_ab_ab * xcomp_i_ba * xcomp_j_ab +
                                  2.0 * f_ab_bb * xcomp_i_ba * xcomp_j_bb +
                                  4.0 * f_aa_bb * xcomp_i_bb * xcomp_j_aa +
                                  2.0 * f_ab_bb * xcomp_i_bb * xcomp_j_ba +
                                  2.0 * f_ab_bb * xcomp_i_bb * xcomp_j_ab +
                                  4.0 * f_bb_bb * xcomp_i_bb * xcomp_j_bb +
                                  2.0 * f_aa_ab * xcomp_i_ab * xcomp_j_aa +
                                        f_ab_ab * xcomp_i_ab * xcomp_j_ba +
                                        f_ab_ab * xcomp_i_ab * xcomp_j_ab +
                                  2.0 * f_ab_bb * xcomp_i_ab * xcomp_j_bb);

                        gatmxy += w * (4.0 * f_aa_aa * xcomp_i_aa * ycomp_j_aa +
                                  2.0 * f_aa_ab * xcomp_i_aa * ycomp_j_ba +
                                  2.0 * f_aa_ab * xcomp_i_aa * ycomp_j_ab +
                                  4.0 * f_aa_bb * xcomp_i_aa * ycomp_j_bb +
                                  2.0 * f_aa_ab * xcomp_i_ba * ycomp_j_aa +
                                        f_ab_ab * xcomp_i_ba * ycomp_j_ba +
                                        f_ab_ab * xcomp_i_ba * ycomp_j_ab +
                                  2.0 * f_ab_bb * xcomp_i_ba * ycomp_j_bb +
                                  4.0 * f_aa_bb * xcomp_i_bb * ycomp_j_aa +
                                  2.0 * f_ab_bb * xcomp_i_bb * ycomp_j_ba +
                                  2.0 * f_ab_bb * xcomp_i_bb * ycomp_j_ab +
                                  4.0 * f_bb_bb * xcomp_i_bb * ycomp_j_bb +
                                  2.0 * f_aa_ab * xcomp_i_ab * ycomp_j_aa +
                                        f_ab_ab * xcomp_i_ab * ycomp_j_ba +
                                        f_ab_ab * xcomp_i_ab * ycomp_j_ab +
                                  2.0 * f_ab_bb * xcomp_i_ab * ycomp_j_bb);

                        gatmxz += w * (4.0 * f_aa_aa * xcomp_i_aa * zcomp_j_aa +
                                  2.0 * f_aa_ab * xcomp_i_aa * zcomp_j_ba +
                                  2.0 * f_aa_ab * xcomp_i_aa * zcomp_j_ab +
                                  4.0 * f_aa_bb * xcomp_i_aa * zcomp_j_bb +
                                  2.0 * f_aa_ab * xcomp_i_ba * zcomp_j_aa +
                                        f_ab_ab * xcomp_i_ba * zcomp_j_ba +
                                        f_ab_ab * xcomp_i_ba * zcomp_j_ab +
                                  2.0 * f_ab_bb * xcomp_i_ba * zcomp_j_bb +
                                  4.0 * f_aa_bb * xcomp_i_bb * zcomp_j_aa +
                                  2.0 * f_ab_bb * xcomp_i_bb * zcomp_j_ba +
                                  2.0 * f_ab_bb * xcomp_i_bb * zcomp_j_ab +
                                  4.0 * f_bb_bb * xcomp_i_bb * zcomp_j_bb +
                                  2.0 * f_aa_ab * xcomp_i_ab * zcomp_j_aa +
                                        f_ab_ab * xcomp_i_ab * zcomp_j_ba +
                                        f_ab_ab * xcomp_i_ab * zcomp_j_ab +
                                  2.0 * f_ab_bb * xcomp_i_ab * zcomp_j_bb);

                        gatmyx += w * (4.0 * f_aa_aa * ycomp_i_aa * xcomp_j_aa +
                                  2.0 * f_aa_ab * ycomp_i_aa * xcomp_j_ba +
                                  2.0 * f_aa_ab * ycomp_i_aa * xcomp_j_ab +
                                  4.0 * f_aa_bb * ycomp_i_aa * xcomp_j_bb +
                                  2.0 * f_aa_ab * ycomp_i_ba * xcomp_j_aa +
                                        f_ab_ab * ycomp_i_ba * xcomp_j_ba +
                                        f_ab_ab * ycomp_i_ba * xcomp_j_ab +
                                  2.0 * f_ab_bb * ycomp_i_ba * xcomp_j_bb +
                                  4.0 * f_aa_bb * ycomp_i_bb * xcomp_j_aa +
                                  2.0 * f_ab_bb * ycomp_i_bb * xcomp_j_ba +
                                  2.0 * f_ab_bb * ycomp_i_bb * xcomp_j_ab +
                                  4.0 * f_bb_bb * ycomp_i_bb * xcomp_j_bb +
                                  2.0 * f_aa_ab * ycomp_i_ab * xcomp_j_aa +
                                        f_ab_ab * ycomp_i_ab * xcomp_j_ba +
                                        f_ab_ab * ycomp_i_ab * xcomp_j_ab +
                                  2.0 * f_ab_bb * ycomp_i_ab * xcomp_j_bb);

                        gatmyy += w * (4.0 * f_aa_aa * ycomp_i_aa * ycomp_j_aa +
                                  2.0 * f_aa_ab * ycomp_i_aa * ycomp_j_ba +
                                  2.0 * f_aa_ab * ycomp_i_aa * ycomp_j_ab +
                                  4.0 * f_aa_bb * ycomp_i_aa * ycomp_j_bb +
                                  2.0 * f_aa_ab * ycomp_i_ba * ycomp_j_aa +
                                        f_ab_ab * ycomp_i_ba * ycomp_j_ba +
                                        f_ab_ab * ycomp_i_ba * ycomp_j_ab +
                                  2.0 * f_ab_bb * ycomp_i_ba * ycomp_j_bb +
                                  4.0 * f_aa_bb * ycomp_i_bb * ycomp_j_aa +
                                  2.0 * f_ab_bb * ycomp_i_bb * ycomp_j_ba +
                                  2.0 * f_ab_bb * ycomp_i_bb * ycomp_j_ab +
                                  4.0 * f_bb_bb * ycomp_i_bb * ycomp_j_bb +
                                  2.0 * f_aa_ab * ycomp_i_ab * ycomp_j_aa +
                                        f_ab_ab * ycomp_i_ab * ycomp_j_ba +
                                        f_ab_ab * ycomp_i_ab * ycomp_j_ab +
                                  2.0 * f_ab_bb * ycomp_i_ab * ycomp_j_bb);

                        gatmyz += w * (4.0 * f_aa_aa * ycomp_i_aa * zcomp_j_aa +
                                  2.0 * f_aa_ab * ycomp_i_aa * zcomp_j_ba +
                                  2.0 * f_aa_ab * ycomp_i_aa * zcomp_j_ab +
                                  4.0 * f_aa_bb * ycomp_i_aa * zcomp_j_bb +
                                  2.0 * f_aa_ab * ycomp_i_ba * zcomp_j_aa +
                                        f_ab_ab * ycomp_i_ba * zcomp_j_ba +
                                        f_ab_ab * ycomp_i_ba * zcomp_j_ab +
                                  2.0 * f_ab_bb * ycomp_i_ba * zcomp_j_bb +
                                  4.0 * f_aa_bb * ycomp_i_bb * zcomp_j_aa +
                                  2.0 * f_ab_bb * ycomp_i_bb * zcomp_j_ba +
                                  2.0 * f_ab_bb * ycomp_i_bb * zcomp_j_ab +
                                  4.0 * f_bb_bb * ycomp_i_bb * zcomp_j_bb +
                                  2.0 * f_aa_ab * ycomp_i_ab * zcomp_j_aa +
                                        f_ab_ab * ycomp_i_ab * zcomp_j_ba +
                                        f_ab_ab * ycomp_i_ab * zcomp_j_ab +
                                  2.0 * f_ab_bb * ycomp_i_ab * zcomp_j_bb);

                        gatmzx += w * (4.0 * f_aa_aa * zcomp_i_aa * xcomp_j_aa +
                                  2.0 * f_aa_ab * zcomp_i_aa * xcomp_j_ba +
                                  2.0 * f_aa_ab * zcomp_i_aa * xcomp_j_ab +
                                  4.0 * f_aa_bb * zcomp_i_aa * xcomp_j_bb +
                                  2.0 * f_aa_ab * zcomp_i_ba * xcomp_j_aa +
                                        f_ab_ab * zcomp_i_ba * xcomp_j_ba +
                                        f_ab_ab * zcomp_i_ba * xcomp_j_ab +
                                  2.0 * f_ab_bb * zcomp_i_ba * xcomp_j_bb +
                                  4.0 * f_aa_bb * zcomp_i_bb * xcomp_j_aa +
                                  2.0 * f_ab_bb * zcomp_i_bb * xcomp_j_ba +
                                  2.0 * f_ab_bb * zcomp_i_bb * xcomp_j_ab +
                                  4.0 * f_bb_bb * zcomp_i_bb * xcomp_j_bb +
                                  2.0 * f_aa_ab * zcomp_i_ab * xcomp_j_aa +
                                        f_ab_ab * zcomp_i_ab * xcomp_j_ba +
                                        f_ab_ab * zcomp_i_ab * xcomp_j_ab +
                                  2.0 * f_ab_bb * zcomp_i_ab * xcomp_j_bb);

                        gatmzy += w * (4.0 * f_aa_aa * zcomp_i_aa * ycomp_j_aa +
                                  2.0 * f_aa_ab * zcomp_i_aa * ycomp_j_ba +
                                  2.0 * f_aa_ab * zcomp_i_aa * ycomp_j_ab +
                                  4.0 * f_aa_bb * zcomp_i_aa * ycomp_j_bb +
                                  2.0 * f_aa_ab * zcomp_i_ba * ycomp_j_aa +
                                        f_ab_ab * zcomp_i_ba * ycomp_j_ba +
                                        f_ab_ab * zcomp_i_ba * ycomp_j_ab +
                                  2.0 * f_ab_bb * zcomp_i_ba * ycomp_j_bb +
                                  4.0 * f_aa_bb * zcomp_i_bb * ycomp_j_aa +
                                  2.0 * f_ab_bb * zcomp_i_bb * ycomp_j_ba +
                                  2.0 * f_ab_bb * zcomp_i_bb * ycomp_j_ab +
                                  4.0 * f_bb_bb * zcomp_i_bb * ycomp_j_bb +
                                  2.0 * f_aa_ab * zcomp_i_ab * ycomp_j_aa +
                                        f_ab_ab * zcomp_i_ab * ycomp_j_ba +
                                        f_ab_ab * zcomp_i_ab * ycomp_j_ab +
                                  2.0 * f_ab_bb * zcomp_i_ab * ycomp_j_bb);

                        gatmzz += w * (4.0 * f_aa_aa * zcomp_i_aa * zcomp_j_aa +
                                  2.0 * f_aa_ab * zcomp_i_aa * zcomp_j_ba +
                                  2.0 * f_aa_ab * zcomp_i_aa * zcomp_j_ab +
                                  4.0 * f_aa_bb * zcomp_i_aa * zcomp_j_bb +
                                  2.0 * f_aa_ab * zcomp_i_ba * zcomp_j_aa +
                                        f_ab_ab * zcomp_i_ba * zcomp_j_ba +
                                        f_ab_ab * zcomp_i_ba * zcomp_j_ab +
                                  2.0 * f_ab_bb * zcomp_i_ba * zcomp_j_bb +
                                  4.0 * f_aa_bb * zcomp_i_bb * zcomp_j_aa +
                                  2.0 * f_ab_bb * zcomp_i_bb * zcomp_j_ba +
                                  2.0 * f_ab_bb * zcomp_i_bb * zcomp_j_ab +
                                  4.0 * f_bb_bb * zcomp_i_bb * zcomp_j_bb +
                                  2.0 * f_aa_ab * zcomp_i_ab * zcomp_j_aa +
                                        f_ab_ab * zcomp_i_ab * zcomp_j_ba +
                                        f_ab_ab * zcomp_i_ab * zcomp_j_ab +
                                  2.0 * f_ab_bb * zcomp_i_ab * zcomp_j_bb);

                        // (2 f_{\sigma_{\alpha\alpha}} + f_{\sigma_{\alpha\beta}})
                        // \nabla\rho_{\alpha}^{(\xi)} \cdot \nabla\rho_{\alpha}^{(\zeta)}

                        // (2 exc_sigma[aa]) * nabla[0]_rho[a]_xi * nabla[0]_rho[a]_zeta
                        // (  exc_sigma[ab]) * nabla[0]_rho[a]_xi * nabla[0]_rho[b]_zeta
                        // (2 exc_sigma[bb]) * nabla[0]_rho[b]_xi * nabla[0]_rho[b]_zeta
                        // (  exc_sigma[ab]) * nabla[0]_rho[b]_xi * nabla[0]_rho[a]_zeta

                        double f_aa = vsigma[3 * g + 0];
                        double f_ab = vsigma[3 * g + 1];
                        double f_bb = vsigma[3 * g + 2];

                        // xx,xy,xz

                        gatmxx += w * 2.0 * f_aa * (gden_a_xx[ig] * gden_a_xx[jg] +
                                                    gden_a_yx[ig] * gden_a_yx[jg] +
                                                    gden_a_zx[ig] * gden_a_zx[jg]);

                        gatmxx += w * f_ab * (gden_a_xx[ig] * gden_b_xx[jg] +
                                              gden_a_yx[ig] * gden_b_yx[jg] +
                                              gden_a_zx[ig] * gden_b_zx[jg]);

                        gatmxx += w * f_ab * (gden_b_xx[ig] * gden_a_xx[jg] +
                                              gden_b_yx[ig] * gden_a_yx[jg] +
                                              gden_b_zx[ig] * gden_a_zx[jg]);

                        gatmxx += w * 2.0 * f_bb * (gden_b_xx[ig] * gden_b_xx[jg] +
                                                    gden_b_yx[ig] * gden_b_yx[jg] +
                                                    gden_b_zx[ig] * gden_b_zx[jg]);

                        gatmxy += w * 2.0 * f_aa * (gden_a_xx[ig] * gden_a_xy[jg] +
                                                    gden_a_yx[ig] * gden_a_yy[jg] +
                                                    gden_a_zx[ig] * gden_a_zy[jg]);

                        gatmxy += w * f_ab * (gden_a_xx[ig] * gden_b_xy[jg] +
                                              gden_a_yx[ig] * gden_b_yy[jg] +
                                              gden_a_zx[ig] * gden_b_zy[jg]);

                        gatmxy += w * f_ab * (gden_b_xx[ig] * gden_a_xy[jg] +
                                              gden_b_yx[ig] * gden_a_yy[jg] +
                                              gden_b_zx[ig] * gden_a_zy[jg]);

                        gatmxy += w * 2.0 * f_bb * (gden_b_xx[ig] * gden_b_xy[jg] +
                                                    gden_b_yx[ig] * gden_b_yy[jg] +
                                                    gden_b_zx[ig] * gden_b_zy[jg]);

                        gatmxz += w * 2.0 * f_aa * (gden_a_xx[ig] * gden_a_xz[jg] +
                                                    gden_a_yx[ig] * gden_a_yz[jg] +
                                                    gden_a_zx[ig] * gden_a_zz[jg]);

                        gatmxz += w * f_ab * (gden_a_xx[ig] * gden_b_xz[jg] +
                                              gden_a_yx[ig] * gden_b_yz[jg] +
                                              gden_a_zx[ig] * gden_b_zz[jg]);

                        gatmxz += w * f_ab * (gden_b_xx[ig] * gden_a_xz[jg] +
                                              gden_b_yx[ig] * gden_a_yz[jg] +
                                              gden_b_zx[ig] * gden_a_zz[jg]);

                        gatmxz += w * 2.0 * f_bb * (gden_b_xx[ig] * gden_b_xz[jg] +
                                                    gden_b_yx[ig] * gden_b_yz[jg] +
                                                    gden_b_zx[ig] * gden_b_zz[jg]);

                        // yx,yy,yz

                        gatmyx += w * 2.0 * f_aa * (gden_a_xy[ig] * gden_a_xx[jg] +
                                                    gden_a_yy[ig] * gden_a_yx[jg] +
                                                    gden_a_zy[ig] * gden_a_zx[jg]);

                        gatmyx += w * f_ab * (gden_a_xy[ig] * gden_b_xx[jg] +
                                              gden_a_yy[ig] * gden_b_yx[jg] +
                                              gden_a_zy[ig] * gden_b_zx[jg]);

                        gatmyx += w * f_ab * (gden_b_xy[ig] * gden_a_xx[jg] +
                                              gden_b_yy[ig] * gden_a_yx[jg] +
                                              gden_b_zy[ig] * gden_a_zx[jg]);

                        gatmyx += w * 2.0 * f_bb * (gden_b_xy[ig] * gden_b_xx[jg] +
                                                    gden_b_yy[ig] * gden_b_yx[jg] +
                                                    gden_b_zy[ig] * gden_b_zx[jg]);

                        gatmyy += w * 2.0 * f_aa * (gden_a_xy[ig] * gden_a_xy[jg] +
                                                    gden_a_yy[ig] * gden_a_yy[jg] +
                                                    gden_a_zy[ig] * gden_a_zy[jg]);

                        gatmyy += w * f_ab * (gden_a_xy[ig] * gden_b_xy[jg] +
                                              gden_a_yy[ig] * gden_b_yy[jg] +
                                              gden_a_zy[ig] * gden_b_zy[jg]);

                        gatmyy += w * f_ab * (gden_b_xy[ig] * gden_a_xy[jg] +
                                              gden_b_yy[ig] * gden_a_yy[jg] +
                                              gden_b_zy[ig] * gden_a_zy[jg]);

                        gatmyy += w * 2.0 * f_bb * (gden_b_xy[ig] * gden_b_xy[jg] +
                                                    gden_b_yy[ig] * gden_b_yy[jg] +
                                                    gden_b_zy[ig] * gden_b_zy[jg]);

                        gatmyz += w * 2.0 * f_aa * (gden_a_xy[ig] * gden_a_xz[jg] +
                                                    gden_a_yy[ig] * gden_a_yz[jg] +
                                                    gden_a_zy[ig] * gden_a_zz[jg]);

                        gatmyz += w * f_ab * (gden_a_xy[ig] * gden_b_xz[jg] +
                                              gden_a_yy[ig] * gden_b_yz[jg] +
                                              gden_a_zy[ig] * gden_b_zz[jg]);

                        gatmyz += w * f_ab * (gden_b_xy[ig] * gden_a_xz[jg] +
                                              gden_b_yy[ig] * gden_a_yz[jg] +
                                              gden_b_zy[ig] * gden_a_zz[jg]);

                        gatmyz += w * 2.0 * f_bb * (gden_b_xy[ig] * gden_b_xz[jg] +
                                                    gden_b_yy[ig] * gden_b_yz[jg] +
                                                    gden_b_zy[ig] * gden_b_zz[jg]);

                        // zx,zy,zz

                        gatmzx += w * 2.0 * f_aa * (gden_a_xz[ig] * gden_a_xx[jg] +
                                                    gden_a_yz[ig] * gden_a_yx[jg] +
                                                    gden_a_zz[ig] * gden_a_zx[jg]);

                        gatmzx += w * f_ab * (gden_a_xz[ig] * gden_b_xx[jg] +
                                              gden_a_yz[ig] * gden_b_yx[jg] +
                                              gden_a_zz[ig] * gden_b_zx[jg]);

                        gatmzx += w * f_ab * (gden_b_xz[ig] * gden_a_xx[jg] +
                                              gden_b_yz[ig] * gden_a_yx[jg] +
                                              gden_b_zz[ig] * gden_a_zx[jg]);

                        gatmzx += w * 2.0 * f_bb * (gden_b_xz[ig] * gden_b_xx[jg] +
                                                    gden_b_yz[ig] * gden_b_yx[jg] +
                                                    gden_b_zz[ig] * gden_b_zx[jg]);

                        gatmzy += w * 2.0 * f_aa * (gden_a_xz[ig] * gden_a_xy[jg] +
                                                    gden_a_yz[ig] * gden_a_yy[jg] +
                                                    gden_a_zz[ig] * gden_a_zy[jg]);

                        gatmzy += w * f_ab * (gden_a_xz[ig] * gden_b_xy[jg] +
                                              gden_a_yz[ig] * gden_b_yy[jg] +
                                              gden_a_zz[ig] * gden_b_zy[jg]);

                        gatmzy += w * f_ab * (gden_b_xz[ig] * gden_a_xy[jg] +
                                              gden_b_yz[ig] * gden_a_yy[jg] +
                                              gden_b_zz[ig] * gden_a_zy[jg]);

                        gatmzy += w * 2.0 * f_bb * (gden_b_xz[ig] * gden_b_xy[jg] +
                                                    gden_b_yz[ig] * gden_b_yy[jg] +
                                                    gden_b_zz[ig] * gden_b_zy[jg]);

                        gatmzz += w * 2.0 * f_aa * (gden_a_xz[ig] * gden_a_xz[jg] +
                                                    gden_a_yz[ig] * gden_a_yz[jg] +
                                                    gden_a_zz[ig] * gden_a_zz[jg]);

                        gatmzz += w * f_ab * (gden_a_xz[ig] * gden_b_xz[jg] +
                                              gden_a_yz[ig] * gden_b_yz[jg] +
                                              gden_a_zz[ig] * gden_b_zz[jg]);

                        gatmzz += w * f_ab * (gden_b_xz[ig] * gden_a_xz[jg] +
                                              gden_b_yz[ig] * gden_a_yz[jg] +
                                              gden_b_zz[ig] * gden_a_zz[jg]);

                        gatmzz += w * 2.0 * f_bb * (gden_b_xz[ig] * gden_b_xz[jg] +
                                                    gden_b_yz[ig] * gden_b_yz[jg] +
                                                    gden_b_zz[ig] * gden_b_zz[jg]);
                    }

                    // no factor of 2 for open-shell

                    gatm[ix * (natoms * 3) + jx] += gatmxx;
                    gatm[ix * (natoms * 3) + jy] += gatmxy;
                    gatm[ix * (natoms * 3) + jz] += gatmxz;

                    gatm[iy * (natoms * 3) + jx] += gatmyx;
                    gatm[iy * (natoms * 3) + jy] += gatmyy;
                    gatm[iy * (natoms * 3) + jz] += gatmyz;

                    gatm[iz * (natoms * 3) + jx] += gatmzx;
                    gatm[iz * (natoms * 3) + jy] += gatmzy;
                    gatm[iz * (natoms * 3) + jz] += gatmzz;
                }
            }

            omptimers[thread_id].stop("Accumulate Hessian");
        }

        timer.stop("OMP Vxc Hessian evaluation");
    }

    timer.stop("Total timing");

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;

    CDenseMatrix molhess(natoms * 3, natoms * 3);

    for (int iatom = 0; iatom < natoms; iatom++)
    {
        auto ix = iatom * 3 + 0;
        auto iy = iatom * 3 + 1;
        auto iz = iatom * 3 + 2;

        for (int jatom = iatom; jatom < natoms; jatom++)
        {
            auto jx = jatom * 3 + 0;
            auto jy = jatom * 3 + 1;
            auto jz = jatom * 3 + 2;

            for (int thread_id = 0; thread_id < nthreads; thread_id++)
            {
                molhess.row(ix)[jx] += molhess_threads.row(thread_id)[ix * (natoms * 3) + jx];
                molhess.row(ix)[jy] += molhess_threads.row(thread_id)[ix * (natoms * 3) + jy];
                molhess.row(ix)[jz] += molhess_threads.row(thread_id)[ix * (natoms * 3) + jz];

                molhess.row(iy)[jx] += molhess_threads.row(thread_id)[iy * (natoms * 3) + jx];
                molhess.row(iy)[jy] += molhess_threads.row(thread_id)[iy * (natoms * 3) + jy];
                molhess.row(iy)[jz] += molhess_threads.row(thread_id)[iy * (natoms * 3) + jz];

                molhess.row(iz)[jx] += molhess_threads.row(thread_id)[iz * (natoms * 3) + jx];
                molhess.row(iz)[jy] += molhess_threads.row(thread_id)[iz * (natoms * 3) + jy];
                molhess.row(iz)[jz] += molhess_threads.row(thread_id)[iz * (natoms * 3) + jz];
            }

            if (jatom != iatom)
            {
                molhess.row(jx)[ix] = molhess.row(ix)[jx];
                molhess.row(jx)[iy] = molhess.row(iy)[jx];
                molhess.row(jx)[iz] = molhess.row(iz)[jx];

                molhess.row(jy)[ix] = molhess.row(ix)[jy];
                molhess.row(jy)[iy] = molhess.row(iy)[jy];
                molhess.row(jy)[iz] = molhess.row(iz)[jy];

                molhess.row(jz)[ix] = molhess.row(ix)[jz];
                molhess.row(jz)[iy] = molhess.row(iy)[jz];
                molhess.row(jz)[iz] = molhess.row(iz)[jz];
            }
        }
    }

    return molhess;
}

auto
integrateVxcFockGradientForGgaClosedShell(const CMolecule&        molecule,
                                          const CMolecularBasis&  basis,
                                          const std::vector<const double*>& gsDensityPointers,
                                          const CMolecularGrid&   molecularGrid,
                                          const double            screeningThresholdForGTOValues,
                                          const CXCFunctional&    xcFunctional,
                                          const std::vector<int>& atomIdxVec) -> std::vector<CDenseMatrix>
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

    // Vxc Fock gradeints (in x,y,z directions)

    auto natoms = static_cast<int>(atomIdxVec.size());

    std::vector<CDenseMatrix> vxcgrads(natoms * 3, CDenseMatrix(naos, naos));

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

        timer.stop("Density matrix slicing");

        // GTO values on grid points

        timer.start("OMP VxcFockGrad calc.");

        std::vector<CDenseMatrix> sum_partial_vxc_gx(natoms, CDenseMatrix(aocount, aocount));
        std::vector<CDenseMatrix> sum_partial_vxc_gy(natoms, CDenseMatrix(aocount, aocount));
        std::vector<CDenseMatrix> sum_partial_vxc_gz(natoms, CDenseMatrix(aocount, aocount));

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

            std::vector<CDenseMatrix> mat_atomvec_chi_x(natoms, CDenseMatrix(aocount, grid_batch_size));
            std::vector<CDenseMatrix> mat_atomvec_chi_y(natoms, CDenseMatrix(aocount, grid_batch_size));
            std::vector<CDenseMatrix> mat_atomvec_chi_z(natoms, CDenseMatrix(aocount, grid_batch_size));

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

                    auto iatom = ao_to_atom_ids[pre_ao_inds[nu]];

                    for (int vecind = 0; vecind < natoms; vecind++)
                    {
                        if (atomIdxVec[vecind] == iatom)
                        {
                            std::memcpy(mat_atomvec_chi_x[vecind].row(idx), submat_x_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                            std::memcpy(mat_atomvec_chi_y[vecind].row(idx), submat_y_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                            std::memcpy(mat_atomvec_chi_z[vecind].row(idx), submat_z_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                        }
                    }

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

            auto v2rho2     = omp_v2rho2_data[thread_id].data();
            auto v2rhosigma = omp_v2rhosigma_data[thread_id].data();
            auto v2sigma2   = omp_v2sigma2_data[thread_id].data();

            sdengridgen::serialGenerateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat);

            omptimers[thread_id].stop("Generate density grid");

            omptimers[thread_id].start("Density grad. grid matmul");

            auto mat_F = sdenblas::serialMultAB(gs_sub_dens_mat, mat_chi);

            auto mat_F_x = sdenblas::serialMultAB(gs_sub_dens_mat, mat_chi_x);
            auto mat_F_y = sdenblas::serialMultAB(gs_sub_dens_mat, mat_chi_y);
            auto mat_F_z = sdenblas::serialMultAB(gs_sub_dens_mat, mat_chi_z);

            omptimers[thread_id].stop("Density grad. grid matmul");

            auto F_val = mat_F.values();

            auto F_x_val = mat_F_x.values();
            auto F_y_val = mat_F_y.values();
            auto F_z_val = mat_F_z.values();

            auto chi_val = mat_chi.values();

            auto chi_x_val = mat_chi_x.values();
            auto chi_y_val = mat_chi_y.values();
            auto chi_z_val = mat_chi_z.values();

            auto chi_xx_val = mat_chi_xx.values();
            auto chi_xy_val = mat_chi_xy.values();
            auto chi_xz_val = mat_chi_xz.values();
            auto chi_yy_val = mat_chi_yy.values();
            auto chi_yz_val = mat_chi_yz.values();
            auto chi_zz_val = mat_chi_zz.values();

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_vxc_for_gga(grid_batch_size, rho, sigma, vrho, vsigma);

            omp_xcfuncs[thread_id].compute_fxc_for_gga(grid_batch_size, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            omptimers[thread_id].start("VxcFockGrad on atoms");

            CDenseMatrix dengradx(3, grid_batch_size);

            CDenseMatrix dengradxx(9, grid_batch_size);

            CDenseMatrix vxc_w(aocount, grid_batch_size);
            CDenseMatrix vxc_wx(aocount, grid_batch_size);
            CDenseMatrix vxc_wy(aocount, grid_batch_size);
            CDenseMatrix vxc_wz(aocount, grid_batch_size);

            for (int vecind = 0; vecind < natoms; vecind++)
            {
                auto atomIdx = atomIdxVec[vecind];

                dengradx.zero();

                auto gdenx = dengradx.row(0);
                auto gdeny = dengradx.row(1);
                auto gdenz = dengradx.row(2);

                dengradxx.zero();

                auto gdenxx = dengradxx.row(0);
                auto gdenxy = dengradxx.row(1);
                auto gdenxz = dengradxx.row(2);

                auto gdenyx = dengradxx.row(3);
                auto gdenyy = dengradxx.row(4);
                auto gdenyz = dengradxx.row(5);

                auto gdenzx = dengradxx.row(6);
                auto gdenzy = dengradxx.row(7);
                auto gdenzz = dengradxx.row(8);

                vxc_w.zero();
                vxc_wx.zero();
                vxc_wy.zero();
                vxc_wz.zero();

                auto vxc_w_val  = vxc_w.values();
                auto vxc_wx_val = vxc_wx.values();
                auto vxc_wy_val = vxc_wy.values();
                auto vxc_wz_val = vxc_wz.values();

                // prepare gradient density

                for (int nu = 0; nu < aocount; nu++)
                {
                    auto iatom = ao_to_atom_ids[aoinds[nu]];

                    if (iatom != atomIdx) continue;

                    auto nu_offset = nu * grid_batch_size;

                    #pragma omp simd 
                    for (int g = 0; g < grid_batch_size; g++)
                    {
                        auto nu_g = nu_offset + g;

                        gdenx[g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];
                        gdeny[g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];
                        gdenz[g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];

                        gdenxx[g] -= 2.0 * (F_x_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xx_val[nu_g]);
                        gdenxy[g] -= 2.0 * (F_x_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_xy_val[nu_g]);
                        gdenxz[g] -= 2.0 * (F_x_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_xz_val[nu_g]);

                        gdenyx[g] -= 2.0 * (F_y_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xy_val[nu_g]);
                        gdenyy[g] -= 2.0 * (F_y_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_yy_val[nu_g]);
                        gdenyz[g] -= 2.0 * (F_y_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_yz_val[nu_g]);

                        gdenzx[g] -= 2.0 * (F_z_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xz_val[nu_g]);
                        gdenzy[g] -= 2.0 * (F_z_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_yz_val[nu_g]);
                        gdenzz[g] -= 2.0 * (F_z_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_zz_val[nu_g]);
                    }
                }

                // Vxc matrix element gradient

                for (int nu = 0; nu < aocount; nu++)
                {
                    auto iatom = ao_to_atom_ids[aoinds[nu]];

                    bool nu_on_atom = (iatom == atomIdx);

                    auto nu_offset = nu * grid_batch_size;

                    #pragma omp simd 
                    for (int g = 0; g < grid_batch_size; g++)
                    {
                        auto nu_g = nu_offset + g;

                        auto w = local_weights[g];

                        // 2 f_{\rho_{\alpha}} \phi_{\mu}^{(\xi)} \phi_{\nu}

                        // note: \phi_{\mu}^{(\xi)} will be added later (from mat_atomvec_chi_{xyz})

                        auto prefac = w * vrho[2 * g + 0];

                        vxc_w_val[nu_g] += -2.0 * prefac * chi_val[nu_g];

                        // (f_{\rho_{\alpha} \rho_{\alpha}} + f_{\rho_{\alpha} \rho_{\beta}})
                        // \rho_{\alpha}^{(\xi)} \phi_{\mu} \phi_{\nu}

                        // note: \phi_{\mu} will be added later (from mat_chi)

                        prefac = w * (v2rho2[3 * g + 0] + v2rho2[3 * g + 1]);

                        vxc_wx_val[nu_g] += prefac * gdenx[g] * chi_val[nu_g];
                        vxc_wy_val[nu_g] += prefac * gdeny[g] * chi_val[nu_g];
                        vxc_wz_val[nu_g] += prefac * gdenz[g] * chi_val[nu_g];

                        // (2 f_{\rho_{\alpha} \sigma_{\alpha\alpha}} +
                        //  2 f_{\rho_{\alpha} \sigma_{\alpha\beta}} +
                        //  2 f_{\rho_{\alpha} \sigma_{\beta\beta}})
                        // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\xi)} \phi_{\mu} \phi_{\nu}

                        // note: \phi_{\mu} will be added later (from mat_chi)

                        prefac = w * 2.0 * (v2rhosigma[6 * g + 0] + v2rhosigma[6 * g + 1] + v2rhosigma[6 * g + 2]);

                        auto gx = rhograd[6 * g + 0];
                        auto gy = rhograd[6 * g + 1];
                        auto gz = rhograd[6 * g + 2];

                        auto xcomp = (gx * gdenxx[g] + gy * gdenyx[g] + gz * gdenzx[g]);
                        auto ycomp = (gx * gdenxy[g] + gy * gdenyy[g] + gz * gdenzy[g]);
                        auto zcomp = (gx * gdenxz[g] + gy * gdenyz[g] + gz * gdenzz[g]);

                        vxc_wx_val[nu_g] += prefac * xcomp * chi_val[nu_g];
                        vxc_wy_val[nu_g] += prefac * ycomp * chi_val[nu_g];
                        vxc_wz_val[nu_g] += prefac * zcomp * chi_val[nu_g];

                        // 2 (2 f_{\sigma_{\alpha\alpha}} + f_{\sigma_{\alpha\beta}})
                        // \nabla\rho_{\alpha}^{(\xi)} \cdot \nabla\phi_{\nu} \phi_{\mu}

                        // note: \phi_{\mu} will be added later (from mat_chi)

                        auto f_aa = vsigma[3 * g + 0];
                        auto f_ab = vsigma[3 * g + 1];

                        prefac = w * (2.0 * f_aa + f_ab);

                        xcomp = (chi_x_val[nu_g] * gdenxx[g] + chi_y_val[nu_g] * gdenyx[g] + chi_z_val[nu_g] * gdenzx[g]);
                        ycomp = (chi_x_val[nu_g] * gdenxy[g] + chi_y_val[nu_g] * gdenyy[g] + chi_z_val[nu_g] * gdenzy[g]);
                        zcomp = (chi_x_val[nu_g] * gdenxz[g] + chi_y_val[nu_g] * gdenyz[g] + chi_z_val[nu_g] * gdenzz[g]);

                        vxc_wx_val[nu_g] += 2.0 * prefac * xcomp;
                        vxc_wy_val[nu_g] += 2.0 * prefac * ycomp;
                        vxc_wz_val[nu_g] += 2.0 * prefac * zcomp;

                        // 2 (2 f_{\sigma_{\alpha\alpha}} + f_{\sigma_{\alpha\beta}})
                        // \nabla\rho_{\alpha} \cdot \nabla\phi_{\nu}^{(\xi)} \phi_{\mu}

                        // note: \phi_{\mu} will be added later (from mat_chi)

                        if (nu_on_atom)
                        {
                            prefac = w * (2.0 * f_aa + f_ab);

                            xcomp = (gx * chi_xx_val[nu_g] + gy * chi_xy_val[nu_g] + gz * chi_xz_val[nu_g]);
                            ycomp = (gx * chi_xy_val[nu_g] + gy * chi_yy_val[nu_g] + gz * chi_yz_val[nu_g]);
                            zcomp = (gx * chi_xz_val[nu_g] + gy * chi_yz_val[nu_g] + gz * chi_zz_val[nu_g]);

                            vxc_wx_val[nu_g] += -2.0 * prefac * xcomp;
                            vxc_wy_val[nu_g] += -2.0 * prefac * ycomp;
                            vxc_wz_val[nu_g] += -2.0 * prefac * zcomp;
                        }

                        // 2 (2 f_{\sigma_{\alpha\alpha}} + f_{\sigma_{\alpha\beta}})
                        // \nabla\rho_{\alpha} \cdot \nabla\phi_{\nu} \phi_{\mu}^{(\xi)}

                        // note: \phi_{\mu}^{(\xi)} will be added later (from mat_atomvec_chi_{xyz})

                        prefac = w * (2.0 * f_aa + f_ab);

                        auto dot_val = (gx * chi_x_val[nu_g] + gy * chi_y_val[nu_g] + gz * chi_z_val[nu_g]);

                        vxc_w_val[nu_g] += -2.0 * prefac * dot_val;

                        // 2 (2 f_{\rho_{\alpha} \sigma_{\alpha\alpha}} +
                        //    2 f_{\rho_{\beta} \sigma_{\alpha\alpha}} +
                        //    f_{\rho_{\alpha} \sigma_{\alpha\beta}} +
                        //    f_{\rho_{\beta} \sigma_{\alpha\beta}})
                        // \rho_{\alpha}^{(\xi)} \nabla\phi_{\nu} \cdot \nabla\rho_{\alpha} \phi_{\mu}

                        // note: \phi_{\mu} will be added later (from mat_chi)

                        f_aa = v2rhosigma[6 * g + 0] + v2rhosigma[6 * g + 3];
                        f_ab = v2rhosigma[6 * g + 1] + v2rhosigma[6 * g + 4];

                        prefac = w * (2.0 * f_aa + f_ab);

                        vxc_wx_val[nu_g] += 2.0 * prefac * dot_val * gdenx[g];
                        vxc_wy_val[nu_g] += 2.0 * prefac * dot_val * gdeny[g];
                        vxc_wz_val[nu_g] += 2.0 * prefac * dot_val * gdenz[g];

                        // 2 (4 f_{\sigma_{\alpha\alpha} \sigma_{\alpha\alpha}} +
                        //    6 f_{\sigma_{\alpha\alpha} \sigma_{\alpha\beta}} +
                        //    4 f_{\sigma_{\alpha\alpha} \sigma_{\beta\beta}} +
                        //    2 f_{\sigma_{\alpha\beta} \sigma_{\alpha\beta}} +
                        //    2 f_{\sigma_{\alpha\beta} \sigma_{\beta\beta}})
                        // \nabla\phi_{\nu} \cdot \nabla\rho_{\alpha}
                        // \nabla\rho_{\alpha} \cdot \nabla\rho_{\alpha}^{(\xi)} \phi_{\mu}

                        // note: \phi_{\mu} will be added later (from mat_chi)

                        f_aa = v2sigma2[6 * g + 0] + v2sigma2[6 * g + 1] + v2sigma2[6 * g + 2];
                        f_ab = v2sigma2[6 * g + 1] + v2sigma2[6 * g + 3] + v2sigma2[6 * g + 4];

                        prefac = w * 2.0 * (2.0 * f_aa + f_ab);

                        xcomp = (gx * gdenxx[g] + gy * gdenyx[g] + gz * gdenzx[g]);
                        ycomp = (gx * gdenxy[g] + gy * gdenyy[g] + gz * gdenzy[g]);
                        zcomp = (gx * gdenxz[g] + gy * gdenyz[g] + gz * gdenzz[g]);

                        vxc_wx_val[nu_g] += 2.0 * prefac * dot_val * xcomp;
                        vxc_wy_val[nu_g] += 2.0 * prefac * dot_val * ycomp;
                        vxc_wz_val[nu_g] += 2.0 * prefac * dot_val * zcomp;
                    }
                }

                auto vxc_gx = sdenblas::serialMultABt(mat_atomvec_chi_x[vecind], vxc_w);
                auto vxc_gy = sdenblas::serialMultABt(mat_atomvec_chi_y[vecind], vxc_w);
                auto vxc_gz = sdenblas::serialMultABt(mat_atomvec_chi_z[vecind], vxc_w);

                auto vxc_gx_2 = sdenblas::serialMultABt(mat_chi, vxc_wx);
                auto vxc_gy_2 = sdenblas::serialMultABt(mat_chi, vxc_wy);
                auto vxc_gz_2 = sdenblas::serialMultABt(mat_chi, vxc_wz);

                sdenblas::serialInPlaceAddAB(vxc_gx, vxc_gx_2);
                sdenblas::serialInPlaceAddAB(vxc_gy, vxc_gy_2);
                sdenblas::serialInPlaceAddAB(vxc_gz, vxc_gz_2);

                vxc_gx.symmetrizeAndScale(0.5);
                vxc_gy.symmetrizeAndScale(0.5);
                vxc_gz.symmetrizeAndScale(0.5);

                #pragma omp critical
                {
                    sdenblas::serialInPlaceAddAB(sum_partial_vxc_gx[vecind], vxc_gx);
                    sdenblas::serialInPlaceAddAB(sum_partial_vxc_gy[vecind], vxc_gy);
                    sdenblas::serialInPlaceAddAB(sum_partial_vxc_gz[vecind], vxc_gz);
                }
            }

            omptimers[thread_id].stop("VxcFockGrad on atoms");
        }

        timer.start("Vxc grad. dist.");

        for (int vecind = 0; vecind < natoms; vecind++)
        {
            dftsubmat::distributeSubMatrixToDenseMatrix(vxcgrads[vecind * 3 + 0], sum_partial_vxc_gx[vecind], aoinds, naos);
            dftsubmat::distributeSubMatrixToDenseMatrix(vxcgrads[vecind * 3 + 1], sum_partial_vxc_gy[vecind], aoinds, naos);
            dftsubmat::distributeSubMatrixToDenseMatrix(vxcgrads[vecind * 3 + 2], sum_partial_vxc_gz[vecind], aoinds, naos);
        }

        timer.stop("Vxc grad. dist.");
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

    return vxcgrads;
}

auto
integrateVxcFockGradientForGgaOpenShell(const CMolecule&        molecule,
                                        const CMolecularBasis&  basis,
                                        const std::vector<const double*>& gsDensityPointers,
                                        const CMolecularGrid&   molecularGrid,
                                        const double            screeningThresholdForGTOValues,
                                        const CXCFunctional&    xcFunctional,
                                        const std::vector<int>& atomIdxVec) -> std::vector<CDenseMatrix>
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

    // Vxc Fock gradeints (in x,y,z directions)

    auto natoms = static_cast<int>(atomIdxVec.size());

    std::vector<CDenseMatrix> vxcgrads(natoms * 6, CDenseMatrix(naos, naos));

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

        timer.stop("Density matrix slicing");

        // GTO values on grid points

        timer.start("OMP VxcFockGrad calc.");

        std::vector<CDenseMatrix> sum_partial_vxc_a_gx(natoms, CDenseMatrix(aocount, aocount));
        std::vector<CDenseMatrix> sum_partial_vxc_a_gy(natoms, CDenseMatrix(aocount, aocount));
        std::vector<CDenseMatrix> sum_partial_vxc_a_gz(natoms, CDenseMatrix(aocount, aocount));

        std::vector<CDenseMatrix> sum_partial_vxc_b_gx(natoms, CDenseMatrix(aocount, aocount));
        std::vector<CDenseMatrix> sum_partial_vxc_b_gy(natoms, CDenseMatrix(aocount, aocount));
        std::vector<CDenseMatrix> sum_partial_vxc_b_gz(natoms, CDenseMatrix(aocount, aocount));

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

            std::vector<CDenseMatrix> mat_atomvec_chi_x(natoms, CDenseMatrix(aocount, grid_batch_size));
            std::vector<CDenseMatrix> mat_atomvec_chi_y(natoms, CDenseMatrix(aocount, grid_batch_size));
            std::vector<CDenseMatrix> mat_atomvec_chi_z(natoms, CDenseMatrix(aocount, grid_batch_size));

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

                    auto iatom = ao_to_atom_ids[pre_ao_inds[nu]];

                    for (int vecind = 0; vecind < natoms; vecind++)
                    {
                        if (atomIdxVec[vecind] == iatom)
                        {
                            std::memcpy(mat_atomvec_chi_x[vecind].row(idx), submat_x_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                            std::memcpy(mat_atomvec_chi_y[vecind].row(idx), submat_y_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                            std::memcpy(mat_atomvec_chi_z[vecind].row(idx), submat_z_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                        }
                    }

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

            auto v2rho2     = omp_v2rho2_data[thread_id].data();
            auto v2rhosigma = omp_v2rhosigma_data[thread_id].data();
            auto v2sigma2   = omp_v2sigma2_data[thread_id].data();

            sdengridgen::serialGenerateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat_a, gs_sub_dens_mat_b);

            omptimers[thread_id].stop("Generate density grid");

            omptimers[thread_id].start("Density grad. grid matmul");

            auto mat_F_a = sdenblas::serialMultAB(gs_sub_dens_mat_a, mat_chi);
            auto mat_F_b = sdenblas::serialMultAB(gs_sub_dens_mat_b, mat_chi);

            auto mat_F_a_x = sdenblas::serialMultAB(gs_sub_dens_mat_a, mat_chi_x);
            auto mat_F_a_y = sdenblas::serialMultAB(gs_sub_dens_mat_a, mat_chi_y);
            auto mat_F_a_z = sdenblas::serialMultAB(gs_sub_dens_mat_a, mat_chi_z);

            auto mat_F_b_x = sdenblas::serialMultAB(gs_sub_dens_mat_b, mat_chi_x);
            auto mat_F_b_y = sdenblas::serialMultAB(gs_sub_dens_mat_b, mat_chi_y);
            auto mat_F_b_z = sdenblas::serialMultAB(gs_sub_dens_mat_b, mat_chi_z);

            omptimers[thread_id].stop("Density grad. grid matmul");

            auto F_a_val = mat_F_a.values();
            auto F_b_val = mat_F_b.values();

            auto F_a_x_val = mat_F_a_x.values();
            auto F_a_y_val = mat_F_a_y.values();
            auto F_a_z_val = mat_F_a_z.values();

            auto F_b_x_val = mat_F_b_x.values();
            auto F_b_y_val = mat_F_b_y.values();
            auto F_b_z_val = mat_F_b_z.values();

            auto chi_val = mat_chi.values();

            auto chi_x_val = mat_chi_x.values();
            auto chi_y_val = mat_chi_y.values();
            auto chi_z_val = mat_chi_z.values();

            auto chi_xx_val = mat_chi_xx.values();
            auto chi_xy_val = mat_chi_xy.values();
            auto chi_xz_val = mat_chi_xz.values();
            auto chi_yy_val = mat_chi_yy.values();
            auto chi_yz_val = mat_chi_yz.values();
            auto chi_zz_val = mat_chi_zz.values();

            // compute exchange-correlation functional derivative

            omptimers[thread_id].start("XC functional eval.");

            omp_xcfuncs[thread_id].compute_vxc_for_gga(grid_batch_size, rho, sigma, vrho, vsigma);

            omp_xcfuncs[thread_id].compute_fxc_for_gga(grid_batch_size, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

            omptimers[thread_id].stop("XC functional eval.");

            omptimers[thread_id].start("Copy grid weights");

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");

            omptimers[thread_id].start("VxcFockGrad on atoms");

            CDenseMatrix dengrad_a_x(3, grid_batch_size);
            CDenseMatrix dengrad_b_x(3, grid_batch_size);

            CDenseMatrix dengrad_a_xx(9, grid_batch_size);
            CDenseMatrix dengrad_b_xx(9, grid_batch_size);

            CDenseMatrix vxc_a_w(aocount, grid_batch_size);
            CDenseMatrix vxc_a_wx(aocount, grid_batch_size);
            CDenseMatrix vxc_a_wy(aocount, grid_batch_size);
            CDenseMatrix vxc_a_wz(aocount, grid_batch_size);

            CDenseMatrix vxc_b_w(aocount, grid_batch_size);
            CDenseMatrix vxc_b_wx(aocount, grid_batch_size);
            CDenseMatrix vxc_b_wy(aocount, grid_batch_size);
            CDenseMatrix vxc_b_wz(aocount, grid_batch_size);

            for (int vecind = 0; vecind < natoms; vecind++)
            {
                auto atomIdx = atomIdxVec[vecind];

                dengrad_a_x.zero();
                dengrad_b_x.zero();

                auto gden_a_x = dengrad_a_x.row(0);
                auto gden_a_y = dengrad_a_x.row(1);
                auto gden_a_z = dengrad_a_x.row(2);

                auto gden_b_x = dengrad_b_x.row(0);
                auto gden_b_y = dengrad_b_x.row(1);
                auto gden_b_z = dengrad_b_x.row(2);

                dengrad_a_xx.zero();
                dengrad_b_xx.zero();

                auto gden_a_xx = dengrad_a_xx.row(0);
                auto gden_a_xy = dengrad_a_xx.row(1);
                auto gden_a_xz = dengrad_a_xx.row(2);

                auto gden_a_yx = dengrad_a_xx.row(3);
                auto gden_a_yy = dengrad_a_xx.row(4);
                auto gden_a_yz = dengrad_a_xx.row(5);

                auto gden_a_zx = dengrad_a_xx.row(6);
                auto gden_a_zy = dengrad_a_xx.row(7);
                auto gden_a_zz = dengrad_a_xx.row(8);

                auto gden_b_xx = dengrad_b_xx.row(0);
                auto gden_b_xy = dengrad_b_xx.row(1);
                auto gden_b_xz = dengrad_b_xx.row(2);

                auto gden_b_yx = dengrad_b_xx.row(3);
                auto gden_b_yy = dengrad_b_xx.row(4);
                auto gden_b_yz = dengrad_b_xx.row(5);

                auto gden_b_zx = dengrad_b_xx.row(6);
                auto gden_b_zy = dengrad_b_xx.row(7);
                auto gden_b_zz = dengrad_b_xx.row(8);

                vxc_a_w.zero();
                vxc_a_wx.zero();
                vxc_a_wy.zero();
                vxc_a_wz.zero();

                vxc_b_w.zero();
                vxc_b_wx.zero();
                vxc_b_wy.zero();
                vxc_b_wz.zero();

                auto vxc_a_w_val  = vxc_a_w.values();
                auto vxc_a_wx_val = vxc_a_wx.values();
                auto vxc_a_wy_val = vxc_a_wy.values();
                auto vxc_a_wz_val = vxc_a_wz.values();

                auto vxc_b_w_val  = vxc_b_w.values();
                auto vxc_b_wx_val = vxc_b_wx.values();
                auto vxc_b_wy_val = vxc_b_wy.values();
                auto vxc_b_wz_val = vxc_b_wz.values();

                // prepare gradient density

                for (int nu = 0; nu < aocount; nu++)
                {
                    auto iatom = ao_to_atom_ids[aoinds[nu]];

                    if (iatom != atomIdx) continue;

                    auto nu_offset = nu * grid_batch_size;

                    #pragma omp simd 
                    for (int g = 0; g < grid_batch_size; g++)
                    {
                        auto nu_g = nu_offset + g;

                        gden_a_x[g] -= 2.0 * F_a_val[nu_g] * chi_x_val[nu_g];
                        gden_a_y[g] -= 2.0 * F_a_val[nu_g] * chi_y_val[nu_g];
                        gden_a_z[g] -= 2.0 * F_a_val[nu_g] * chi_z_val[nu_g];

                        gden_b_x[g] -= 2.0 * F_b_val[nu_g] * chi_x_val[nu_g];
                        gden_b_y[g] -= 2.0 * F_b_val[nu_g] * chi_y_val[nu_g];
                        gden_b_z[g] -= 2.0 * F_b_val[nu_g] * chi_z_val[nu_g];

                        gden_a_xx[g] -= 2.0 * (F_a_x_val[nu_g] * chi_x_val[nu_g] + F_a_val[nu_g] * chi_xx_val[nu_g]);
                        gden_a_xy[g] -= 2.0 * (F_a_x_val[nu_g] * chi_y_val[nu_g] + F_a_val[nu_g] * chi_xy_val[nu_g]);
                        gden_a_xz[g] -= 2.0 * (F_a_x_val[nu_g] * chi_z_val[nu_g] + F_a_val[nu_g] * chi_xz_val[nu_g]);

                        gden_a_yx[g] -= 2.0 * (F_a_y_val[nu_g] * chi_x_val[nu_g] + F_a_val[nu_g] * chi_xy_val[nu_g]);
                        gden_a_yy[g] -= 2.0 * (F_a_y_val[nu_g] * chi_y_val[nu_g] + F_a_val[nu_g] * chi_yy_val[nu_g]);
                        gden_a_yz[g] -= 2.0 * (F_a_y_val[nu_g] * chi_z_val[nu_g] + F_a_val[nu_g] * chi_yz_val[nu_g]);

                        gden_a_zx[g] -= 2.0 * (F_a_z_val[nu_g] * chi_x_val[nu_g] + F_a_val[nu_g] * chi_xz_val[nu_g]);
                        gden_a_zy[g] -= 2.0 * (F_a_z_val[nu_g] * chi_y_val[nu_g] + F_a_val[nu_g] * chi_yz_val[nu_g]);
                        gden_a_zz[g] -= 2.0 * (F_a_z_val[nu_g] * chi_z_val[nu_g] + F_a_val[nu_g] * chi_zz_val[nu_g]);

                        gden_b_xx[g] -= 2.0 * (F_b_x_val[nu_g] * chi_x_val[nu_g] + F_b_val[nu_g] * chi_xx_val[nu_g]);
                        gden_b_xy[g] -= 2.0 * (F_b_x_val[nu_g] * chi_y_val[nu_g] + F_b_val[nu_g] * chi_xy_val[nu_g]);
                        gden_b_xz[g] -= 2.0 * (F_b_x_val[nu_g] * chi_z_val[nu_g] + F_b_val[nu_g] * chi_xz_val[nu_g]);

                        gden_b_yx[g] -= 2.0 * (F_b_y_val[nu_g] * chi_x_val[nu_g] + F_b_val[nu_g] * chi_xy_val[nu_g]);
                        gden_b_yy[g] -= 2.0 * (F_b_y_val[nu_g] * chi_y_val[nu_g] + F_b_val[nu_g] * chi_yy_val[nu_g]);
                        gden_b_yz[g] -= 2.0 * (F_b_y_val[nu_g] * chi_z_val[nu_g] + F_b_val[nu_g] * chi_yz_val[nu_g]);

                        gden_b_zx[g] -= 2.0 * (F_b_z_val[nu_g] * chi_x_val[nu_g] + F_b_val[nu_g] * chi_xz_val[nu_g]);
                        gden_b_zy[g] -= 2.0 * (F_b_z_val[nu_g] * chi_y_val[nu_g] + F_b_val[nu_g] * chi_yz_val[nu_g]);
                        gden_b_zz[g] -= 2.0 * (F_b_z_val[nu_g] * chi_z_val[nu_g] + F_b_val[nu_g] * chi_zz_val[nu_g]);
                    }
                }

                // Vxc matrix element gradient

                for (int nu = 0; nu < aocount; nu++)
                {
                    auto iatom = ao_to_atom_ids[aoinds[nu]];

                    bool nu_on_atom = (iatom == atomIdx);

                    auto nu_offset = nu * grid_batch_size;

                    #pragma omp simd 
                    for (int g = 0; g < grid_batch_size; g++)
                    {
                        auto nu_g = nu_offset + g;

                        auto w = local_weights[g];

                        // alpha
                        // (exc_rho[a]) * phi[mu]_xi * phi[nu]
                        // (exc_rho[a]) * phi[mu] * phi[nu]_xi
                        // beta
                        // (exc_rho[b]) * phi[mu]_xi * phi[nu]
                        // (exc_rho[b]) * phi[mu] * phi[nu]_xi

                        // note: phi[mu]_xi will be added later (from mat_atomvec_chi_{xyz})

                        auto prefac_a = w * vrho[2 * g + 0];
                        auto prefac_b = w * vrho[2 * g + 1];

                        vxc_a_w_val[nu_g] += -2.0 * prefac_a * chi_val[nu_g];
                        vxc_b_w_val[nu_g] += -2.0 * prefac_b * chi_val[nu_g];

                        // alpha
                        // (exc_rho[a]_rho[a]) * phi[mu] * phi[nu] * rho[a]_xi
                        // (exc_rho[a]_rho[b]) * phi[mu] * phi[nu] * rho[b]_xi
                        // beta
                        // (exc_rho[a]_rho[b]) * phi[mu] * phi[nu] * rho[a]_xi
                        // (exc_rho[b]_rho[b]) * phi[mu] * phi[nu] * rho[b]_xi

                        // note: phi[mu] will be added later (from mat_chi)

                        auto prefac_aa = w * v2rho2[3 * g + 0] * chi_val[nu_g];
                        auto prefac_ab = w * v2rho2[3 * g + 1] * chi_val[nu_g];
                        auto prefac_bb = w * v2rho2[3 * g + 2] * chi_val[nu_g];

                        vxc_a_wx_val[nu_g] += prefac_aa * gden_a_x[g] + prefac_ab * gden_b_x[g];
                        vxc_a_wy_val[nu_g] += prefac_aa * gden_a_y[g] + prefac_ab * gden_b_y[g];
                        vxc_a_wz_val[nu_g] += prefac_aa * gden_a_z[g] + prefac_ab * gden_b_z[g];

                        vxc_b_wx_val[nu_g] += prefac_ab * gden_a_x[g] + prefac_bb * gden_b_x[g];
                        vxc_b_wy_val[nu_g] += prefac_ab * gden_a_y[g] + prefac_bb * gden_b_y[g];
                        vxc_b_wz_val[nu_g] += prefac_ab * gden_a_z[g] + prefac_bb * gden_b_z[g];

                        // alpha
                        // (2 exc_rho[a]_sigma[aa]) * nabla[1]_rho[a]    * nabla[1]_rho[a]_xi * phi[mu] * phi[nu]
                        // (  exc_rho[a]_sigma[ab]) * nabla[1]_rho[a]_xi * nabla[1]_rho[b]    * phi[mu] * phi[nu]
                        // (  exc_rho[a]_sigma[ab]) * nabla[1]_rho[a]    * nabla[1]_rho[b]_xi * phi[mu] * phi[nu]
                        // (2 exc_rho[a]_sigma[bb]) * nabla[1]_rho[b]    * nabla[1]_rho[b]_xi * phi[mu] * phi[nu]
                        // beta
                        // (2 exc_rho[b]_sigma[aa]) * nabla[1]_rho[a]    * nabla[1]_rho[a]_xi * phi[mu] * phi[nu]
                        // (  exc_rho[b]_sigma[ab]) * nabla[1]_rho[a]_xi * nabla[1]_rho[b]    * phi[mu] * phi[nu]
                        // (  exc_rho[b]_sigma[ab]) * nabla[1]_rho[a]    * nabla[1]_rho[b]_xi * phi[mu] * phi[nu]
                        // (2 exc_rho[b]_sigma[bb]) * nabla[1]_rho[b]    * nabla[1]_rho[b]_xi * phi[mu] * phi[nu]

                        // note: phi[mu] will be added later (from mat_chi)

                        auto g_a_x = rhograd[6 * g + 0];
                        auto g_a_y = rhograd[6 * g + 1];
                        auto g_a_z = rhograd[6 * g + 2];

                        auto g_b_x = rhograd[6 * g + 3];
                        auto g_b_y = rhograd[6 * g + 4];
                        auto g_b_z = rhograd[6 * g + 5];

                        auto xcomp_aa = 2.0 * (g_a_x * gden_a_xx[g] + g_a_y * gden_a_yx[g] + g_a_z * gden_a_zx[g]);
                        auto ycomp_aa = 2.0 * (g_a_x * gden_a_xy[g] + g_a_y * gden_a_yy[g] + g_a_z * gden_a_zy[g]);
                        auto zcomp_aa = 2.0 * (g_a_x * gden_a_xz[g] + g_a_y * gden_a_yz[g] + g_a_z * gden_a_zz[g]);

                        auto xcomp_ab = (g_a_x * gden_b_xx[g] + g_a_y * gden_b_yx[g] + g_a_z * gden_b_zx[g]);
                        auto ycomp_ab = (g_a_x * gden_b_xy[g] + g_a_y * gden_b_yy[g] + g_a_z * gden_b_zy[g]);
                        auto zcomp_ab = (g_a_x * gden_b_xz[g] + g_a_y * gden_b_yz[g] + g_a_z * gden_b_zz[g]);

                        xcomp_ab += (g_b_x * gden_a_xx[g] + g_b_y * gden_a_yx[g] + g_b_z * gden_a_zx[g]);
                        ycomp_ab += (g_b_x * gden_a_xy[g] + g_b_y * gden_a_yy[g] + g_b_z * gden_a_zy[g]);
                        zcomp_ab += (g_b_x * gden_a_xz[g] + g_b_y * gden_a_yz[g] + g_b_z * gden_a_zz[g]);

                        auto xcomp_bb = 2.0 * (g_b_x * gden_b_xx[g] + g_b_y * gden_b_yx[g] + g_b_z * gden_b_zx[g]);
                        auto ycomp_bb = 2.0 * (g_b_x * gden_b_xy[g] + g_b_y * gden_b_yy[g] + g_b_z * gden_b_zy[g]);
                        auto zcomp_bb = 2.0 * (g_b_x * gden_b_xz[g] + g_b_y * gden_b_yz[g] + g_b_z * gden_b_zz[g]);

                        auto f_a_aa = v2rhosigma[6 * g + 0];
                        auto f_a_ab = v2rhosigma[6 * g + 1];
                        auto f_a_bb = v2rhosigma[6 * g + 2];
                        auto f_b_aa = v2rhosigma[6 * g + 3];
                        auto f_b_ab = v2rhosigma[6 * g + 4];
                        auto f_b_bb = v2rhosigma[6 * g + 5];

                        prefac_a = w * chi_val[nu_g];
                        prefac_b = w * chi_val[nu_g];

                        vxc_a_wx_val[nu_g] += prefac_a * (f_a_aa * xcomp_aa + f_a_ab * xcomp_ab + f_a_bb * xcomp_bb);
                        vxc_a_wy_val[nu_g] += prefac_a * (f_a_aa * ycomp_aa + f_a_ab * ycomp_ab + f_a_bb * ycomp_bb);
                        vxc_a_wz_val[nu_g] += prefac_a * (f_a_aa * zcomp_aa + f_a_ab * zcomp_ab + f_a_bb * zcomp_bb);

                        vxc_b_wx_val[nu_g] += prefac_b * (f_b_aa * xcomp_aa + f_b_ab * xcomp_ab + f_b_bb * xcomp_bb);
                        vxc_b_wy_val[nu_g] += prefac_b * (f_b_aa * ycomp_aa + f_b_ab * ycomp_ab + f_b_bb * ycomp_bb);
                        vxc_b_wz_val[nu_g] += prefac_b * (f_b_aa * zcomp_aa + f_b_ab * zcomp_ab + f_b_bb * zcomp_bb);

                        // alpha
                        // (2 exc_sigma[aa]) * nabla[0]_phi[nu] * nabla[0]_rho[a]_xi * phi[mu]
                        // (  exc_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[b]_xi * phi[mu]
                        // beta
                        // (2 exc_sigma[bb]) * nabla[0]_phi[nu] * nabla[0]_rho[b]_xi * phi[mu]
                        // (  exc_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[a]_xi * phi[mu]

                        // note: need an extra factor of 2

                        // note: \phi_{\mu} will be added later (from mat_chi)

                        auto f_aa = vsigma[3 * g + 0];
                        auto f_ab = vsigma[3 * g + 1];
                        auto f_bb = vsigma[3 * g + 2];

                        auto xcomp_a = chi_x_val[nu_g] * gden_a_xx[g] + chi_y_val[nu_g] * gden_a_yx[g] + chi_z_val[nu_g] * gden_a_zx[g];
                        auto ycomp_a = chi_x_val[nu_g] * gden_a_xy[g] + chi_y_val[nu_g] * gden_a_yy[g] + chi_z_val[nu_g] * gden_a_zy[g];
                        auto zcomp_a = chi_x_val[nu_g] * gden_a_xz[g] + chi_y_val[nu_g] * gden_a_yz[g] + chi_z_val[nu_g] * gden_a_zz[g];

                        auto xcomp_b = chi_x_val[nu_g] * gden_b_xx[g] + chi_y_val[nu_g] * gden_b_yx[g] + chi_z_val[nu_g] * gden_b_zx[g];
                        auto ycomp_b = chi_x_val[nu_g] * gden_b_xy[g] + chi_y_val[nu_g] * gden_b_yy[g] + chi_z_val[nu_g] * gden_b_zy[g];
                        auto zcomp_b = chi_x_val[nu_g] * gden_b_xz[g] + chi_y_val[nu_g] * gden_b_yz[g] + chi_z_val[nu_g] * gden_b_zz[g];

                        vxc_a_wx_val[nu_g] += 2.0 * w * (2.0 * f_aa * xcomp_a + f_ab * xcomp_b);
                        vxc_a_wy_val[nu_g] += 2.0 * w * (2.0 * f_aa * ycomp_a + f_ab * ycomp_b);
                        vxc_a_wz_val[nu_g] += 2.0 * w * (2.0 * f_aa * zcomp_a + f_ab * zcomp_b);

                        vxc_b_wx_val[nu_g] += 2.0 * w * (2.0 * f_bb * xcomp_b + f_ab * xcomp_a);
                        vxc_b_wy_val[nu_g] += 2.0 * w * (2.0 * f_bb * ycomp_b + f_ab * ycomp_a);
                        vxc_b_wz_val[nu_g] += 2.0 * w * (2.0 * f_bb * zcomp_b + f_ab * zcomp_a);

                        // alpha
                        // 2 (2 exc_sigma[aa] + exc_sigma[ab]) * nabla[0]_phi[nu]_xi * nabla[0]_rho[a] * phi[mu]
                        // beta
                        // 2 (2 exc_sigma[bb] + exc_sigma[ab]) * nabla[0]_phi[nu]_xi * nabla[0]_rho[b] * phi[mu]

                        // alpha
                        // (2 exc_sigma[aa]) * nabla[0]_phi[nu]_xi * nabla[0]_rho[a] * phi[mu]
                        // (  exc_sigma[ab]) * nabla[0]_phi[nu]_xi * nabla[0]_rho[b] * phi[mu]
                        // beta
                        // (2 exc_sigma[bb]) * nabla[0]_phi[nu]_xi * nabla[0]_rho[b] * phi[mu]
                        // (  exc_sigma[ab]) * nabla[0]_phi[nu]_xi * nabla[0]_rho[a] * phi[mu]

                        // note: need an extra factor of 2

                        // note: \phi_{\mu} will be added later (from mat_chi)

                        if (nu_on_atom)
                        {
                            prefac_a = w * (2.0 * f_aa + f_ab);
                            prefac_b = w * (2.0 * f_bb + f_ab);

                            xcomp_a = g_a_x * chi_xx_val[nu_g] + g_a_y * chi_xy_val[nu_g] + g_a_z * chi_xz_val[nu_g];
                            ycomp_a = g_a_x * chi_xy_val[nu_g] + g_a_y * chi_yy_val[nu_g] + g_a_z * chi_yz_val[nu_g];
                            zcomp_a = g_a_x * chi_xz_val[nu_g] + g_a_y * chi_yz_val[nu_g] + g_a_z * chi_zz_val[nu_g];

                            xcomp_b = g_b_x * chi_xx_val[nu_g] + g_b_y * chi_xy_val[nu_g] + g_b_z * chi_xz_val[nu_g];
                            ycomp_b = g_b_x * chi_xy_val[nu_g] + g_b_y * chi_yy_val[nu_g] + g_b_z * chi_yz_val[nu_g];
                            zcomp_b = g_b_x * chi_xz_val[nu_g] + g_b_y * chi_yz_val[nu_g] + g_b_z * chi_zz_val[nu_g];

                            vxc_a_wx_val[nu_g] += -2.0 * w * (2.0 * f_aa * xcomp_a + f_ab * xcomp_b);
                            vxc_a_wy_val[nu_g] += -2.0 * w * (2.0 * f_aa * ycomp_a + f_ab * ycomp_b);
                            vxc_a_wz_val[nu_g] += -2.0 * w * (2.0 * f_aa * zcomp_a + f_ab * zcomp_b);

                            vxc_b_wx_val[nu_g] += -2.0 * w * (2.0 * f_bb * xcomp_b + f_ab * xcomp_a);
                            vxc_b_wy_val[nu_g] += -2.0 * w * (2.0 * f_bb * ycomp_b + f_ab * ycomp_a);
                            vxc_b_wz_val[nu_g] += -2.0 * w * (2.0 * f_bb * zcomp_b + f_ab * zcomp_a);
                        }

                        // alpha
                        // 2 (2 exc_sigma[aa] + exc_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[a] * phi[mu]_xi
                        // beta
                        // 2 (2 exc_sigma[bb] + exc_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[b] * phi[mu]_xi

                        // alpha
                        // (2 exc_sigma[aa]) * nabla[0]_phi[nu] * nabla[0]_rho[a] * phi[mu]_xi
                        // (  exc_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[b] * phi[mu]_xi
                        // beta
                        // (2 exc_sigma[bb]) * nabla[0]_phi[nu] * nabla[0]_rho[b] * phi[mu]_xi
                        // (  exc_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[a] * phi[mu]_xi

                        // note: need an extra factor of 2

                        // note: \phi_{\mu}^{(\xi)} will be added later (from mat_atomvec_chi_{xyz})

                        auto dot_val_a = (g_a_x * chi_x_val[nu_g] + g_a_y * chi_y_val[nu_g] + g_a_z * chi_z_val[nu_g]);
                        auto dot_val_b = (g_b_x * chi_x_val[nu_g] + g_b_y * chi_y_val[nu_g] + g_b_z * chi_z_val[nu_g]);

                        vxc_a_w_val[nu_g] += -2.0 * w * (2.0 * f_aa * dot_val_a + f_ab * dot_val_b);
                        vxc_b_w_val[nu_g] += -2.0 * w * (2.0 * f_bb * dot_val_b + f_ab * dot_val_a);

                        // alpha
                        // (2 exc_rho[a]_sigma[aa]) * nabla[0]_phi[nu] * nabla[0]_rho[a] * phi[mu] * rho[a]_xi
                        // (2 exc_rho[b]_sigma[aa]) * nabla[0]_phi[nu] * nabla[0]_rho[a] * phi[mu] * rho[b]_xi
                        // (  exc_rho[a]_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[b] * phi[mu] * rho[a]_xi
                        // (  exc_rho[b]_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[b] * phi[mu] * rho[b]_xi
                        // beta
                        // (2 exc_rho[a]_sigma[bb]) * nabla[0]_phi[nu] * nabla[0]_rho[b] * phi[mu] * rho[a]_xi
                        // (2 exc_rho[b]_sigma[bb]) * nabla[0]_phi[nu] * nabla[0]_rho[b] * phi[mu] * rho[b]_xi
                        // (  exc_rho[a]_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[a] * phi[mu] * rho[a]_xi
                        // (  exc_rho[b]_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[a] * phi[mu] * rho[b]_xi

                        // note: need an extra factor of 2

                        // note: \phi_{\mu} will be added later (from mat_chi)

                        f_a_aa = v2rhosigma[6 * g + 0];
                        f_a_ab = v2rhosigma[6 * g + 1];
                        f_a_bb = v2rhosigma[6 * g + 2];
                        f_b_aa = v2rhosigma[6 * g + 3];
                        f_b_ab = v2rhosigma[6 * g + 4];
                        f_b_bb = v2rhosigma[6 * g + 5];

                        vxc_a_wx_val[nu_g] += 2.0 * w * (2.0 * f_a_aa * dot_val_a * gden_a_x[g] +
                                                         2.0 * f_b_aa * dot_val_a * gden_b_x[g] +
                                                               f_a_ab * dot_val_b * gden_a_x[g] +
                                                               f_b_ab * dot_val_b * gden_b_x[g]);
                        vxc_a_wy_val[nu_g] += 2.0 * w * (2.0 * f_a_aa * dot_val_a * gden_a_y[g] +
                                                         2.0 * f_b_aa * dot_val_a * gden_b_y[g] +
                                                               f_a_ab * dot_val_b * gden_a_y[g] +
                                                               f_b_ab * dot_val_b * gden_b_y[g]);
                        vxc_a_wz_val[nu_g] += 2.0 * w * (2.0 * f_a_aa * dot_val_a * gden_a_z[g] +
                                                         2.0 * f_b_aa * dot_val_a * gden_b_z[g] +
                                                               f_a_ab * dot_val_b * gden_a_z[g] +
                                                               f_b_ab * dot_val_b * gden_b_z[g]);

                        vxc_b_wx_val[nu_g] += 2.0 * w * (2.0 * f_a_bb * dot_val_b * gden_a_x[g] +
                                                         2.0 * f_b_bb * dot_val_b * gden_b_x[g] +
                                                               f_a_ab * dot_val_a * gden_a_x[g] +
                                                               f_b_ab * dot_val_a * gden_b_x[g]);
                        vxc_b_wy_val[nu_g] += 2.0 * w * (2.0 * f_a_bb * dot_val_b * gden_a_y[g] +
                                                         2.0 * f_b_bb * dot_val_b * gden_b_y[g] +
                                                               f_a_ab * dot_val_a * gden_a_y[g] +
                                                               f_b_ab * dot_val_a * gden_b_y[g]);
                        vxc_b_wz_val[nu_g] += 2.0 * w * (2.0 * f_a_bb * dot_val_b * gden_a_z[g] +
                                                         2.0 * f_b_bb * dot_val_b * gden_b_z[g] +
                                                               f_a_ab * dot_val_a * gden_a_z[g] +
                                                               f_b_ab * dot_val_a * gden_b_z[g]);

                        // alpha
                        // (4 exc_sigma[aa]_sigma[aa]) * nabla[0]_phi[nu] * nabla[0]_rho[a] * nabla[1]_rho[a] * nabla[1]_rho[a]_xi * phi[mu]
                        // (2 exc_sigma[aa]_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[a] * nabla[1]_rho[a]_xi * nabla[1]_rho[b] * phi[mu]
                        // (2 exc_sigma[aa]_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[a] * nabla[1]_rho[a] * nabla[1]_rho[b]_xi * phi[mu]
                        // (4 exc_sigma[aa]_sigma[bb]) * nabla[0]_phi[nu] * nabla[0]_rho[a] * nabla[1]_rho[b] * nabla[1]_rho[b]_xi * phi[mu]
                        // (2 exc_sigma[aa]_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[b] * nabla[1]_rho[a] * nabla[1]_rho[a]_xi * phi[mu]
                        // (exc_sigma[ab]_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[b] * nabla[1]_rho[a]_xi * nabla[1]_rho[b] * phi[mu]
                        // (exc_sigma[ab]_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[b] * nabla[1]_rho[a] * nabla[1]_rho[b]_xi * phi[mu]
                        // (2 exc_sigma[ab]_sigma[bb]) * nabla[0]_phi[nu] * nabla[0]_rho[b] * nabla[1]_rho[b] * nabla[1]_rho[b]_xi * phi[mu]
                        // beta
                        // (4 exc_sigma[aa]_sigma[bb]) * nabla[0]_phi[nu] * nabla[0]_rho[b] * nabla[1]_rho[a] * nabla[1]_rho[a]_xi * phi[mu]
                        // (2 exc_sigma[ab]_sigma[bb]) * nabla[0]_phi[nu] * nabla[0]_rho[b] * nabla[1]_rho[a]_xi * nabla[1]_rho[b] * phi[mu]
                        // (2 exc_sigma[ab]_sigma[bb]) * nabla[0]_phi[nu] * nabla[0]_rho[b] * nabla[1]_rho[a] * nabla[1]_rho[b]_xi * phi[mu]
                        // (4 exc_sigma[bb]_sigma[bb]) * nabla[0]_phi[nu] * nabla[0]_rho[b] * nabla[1]_rho[b] * nabla[1]_rho[b]_xi * phi[mu]
                        // (2 exc_sigma[aa]_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[a] * nabla[1]_rho[a] * nabla[1]_rho[a]_xi * phi[mu]
                        // (  exc_sigma[ab]_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[a] * nabla[1]_rho[a]_xi * nabla[1]_rho[b] * phi[mu]
                        // (  exc_sigma[ab]_sigma[ab]) * nabla[0]_phi[nu] * nabla[0]_rho[a] * nabla[1]_rho[a] * nabla[1]_rho[b]_xi * phi[mu]
                        // (2 exc_sigma[ab]_sigma[bb]) * nabla[0]_phi[nu] * nabla[0]_rho[a] * nabla[1]_rho[b] * nabla[1]_rho[b]_xi * phi[mu]

                        // note: need an extra factor of 2

                        // note: \phi_{\mu} will be added later (from mat_chi)

                        auto f_aa_aa = v2sigma2[6 * g + 0];
                        auto f_aa_ab = v2sigma2[6 * g + 1];
                        auto f_aa_bb = v2sigma2[6 * g + 2];
                        auto f_ab_ab = v2sigma2[6 * g + 3];
                        auto f_ab_bb = v2sigma2[6 * g + 4];
                        auto f_bb_bb = v2sigma2[6 * g + 5];

                        xcomp_aa = 2.0 * (g_a_x * gden_a_xx[g] + g_a_y * gden_a_yx[g] + g_a_z * gden_a_zx[g]);
                        ycomp_aa = 2.0 * (g_a_x * gden_a_xy[g] + g_a_y * gden_a_yy[g] + g_a_z * gden_a_zy[g]);
                        zcomp_aa = 2.0 * (g_a_x * gden_a_xz[g] + g_a_y * gden_a_yz[g] + g_a_z * gden_a_zz[g]);

                        xcomp_ab = g_a_x * gden_b_xx[g] + g_a_y * gden_b_yx[g] + g_a_z * gden_b_zx[g];
                        ycomp_ab = g_a_x * gden_b_xy[g] + g_a_y * gden_b_yy[g] + g_a_z * gden_b_zy[g];
                        zcomp_ab = g_a_x * gden_b_xz[g] + g_a_y * gden_b_yz[g] + g_a_z * gden_b_zz[g];

                        xcomp_ab += g_b_x * gden_a_xx[g] + g_b_y * gden_a_yx[g] + g_b_z * gden_a_zx[g];
                        ycomp_ab += g_b_x * gden_a_xy[g] + g_b_y * gden_a_yy[g] + g_b_z * gden_a_zy[g];
                        zcomp_ab += g_b_x * gden_a_xz[g] + g_b_y * gden_a_yz[g] + g_b_z * gden_a_zz[g];

                        xcomp_bb = 2.0 * (g_b_x * gden_b_xx[g] + g_b_y * gden_b_yx[g] + g_b_z * gden_b_zx[g]);
                        ycomp_bb = 2.0 * (g_b_x * gden_b_xy[g] + g_b_y * gden_b_yy[g] + g_b_z * gden_b_zy[g]);
                        zcomp_bb = 2.0 * (g_b_x * gden_b_xz[g] + g_b_y * gden_b_yz[g] + g_b_z * gden_b_zz[g]);

                        vxc_a_wx_val[nu_g] += 2.0 * w * (2.0 * f_aa_aa * dot_val_a * xcomp_aa +
                                                         2.0 * f_aa_ab * dot_val_a * xcomp_ab +
                                                         2.0 * f_aa_bb * dot_val_a * xcomp_bb +
                                                               f_aa_ab * dot_val_b * xcomp_aa +
                                                               f_ab_ab * dot_val_b * xcomp_ab +
                                                               f_ab_bb * dot_val_b * xcomp_bb);
                        vxc_a_wy_val[nu_g] += 2.0 * w * (2.0 * f_aa_aa * dot_val_a * ycomp_aa +
                                                         2.0 * f_aa_ab * dot_val_a * ycomp_ab +
                                                         2.0 * f_aa_bb * dot_val_a * ycomp_bb +
                                                               f_aa_ab * dot_val_b * ycomp_aa +
                                                               f_ab_ab * dot_val_b * ycomp_ab +
                                                               f_ab_bb * dot_val_b * ycomp_bb);
                        vxc_a_wz_val[nu_g] += 2.0 * w * (2.0 * f_aa_aa * dot_val_a * zcomp_aa +
                                                         2.0 * f_aa_ab * dot_val_a * zcomp_ab +
                                                         2.0 * f_aa_bb * dot_val_a * zcomp_bb +
                                                               f_aa_ab * dot_val_b * zcomp_aa +
                                                               f_ab_ab * dot_val_b * zcomp_ab +
                                                               f_ab_bb * dot_val_b * zcomp_bb);

                        vxc_b_wx_val[nu_g] += 2.0 * w * (2.0 * f_aa_bb * dot_val_b * xcomp_aa +
                                                         2.0 * f_ab_bb * dot_val_b * xcomp_ab +
                                                         2.0 * f_bb_bb * dot_val_b * xcomp_bb +
                                                               f_aa_ab * dot_val_a * xcomp_aa +
                                                               f_ab_ab * dot_val_a * xcomp_ab +
                                                               f_ab_bb * dot_val_a * xcomp_bb);
                        vxc_b_wy_val[nu_g] += 2.0 * w * (2.0 * f_aa_bb * dot_val_b * ycomp_aa +
                                                         2.0 * f_ab_bb * dot_val_b * ycomp_ab +
                                                         2.0 * f_bb_bb * dot_val_b * ycomp_bb +
                                                               f_aa_ab * dot_val_a * ycomp_aa +
                                                               f_ab_ab * dot_val_a * ycomp_ab +
                                                               f_ab_bb * dot_val_a * ycomp_bb);
                        vxc_b_wz_val[nu_g] += 2.0 * w * (2.0 * f_aa_bb * dot_val_b * zcomp_aa +
                                                         2.0 * f_ab_bb * dot_val_b * zcomp_ab +
                                                         2.0 * f_bb_bb * dot_val_b * zcomp_bb +
                                                               f_aa_ab * dot_val_a * zcomp_aa +
                                                               f_ab_ab * dot_val_a * zcomp_ab +
                                                               f_ab_bb * dot_val_a * zcomp_bb);
                    }
                }

                auto vxc_a_gx = sdenblas::serialMultABt(mat_atomvec_chi_x[vecind], vxc_a_w);
                auto vxc_a_gy = sdenblas::serialMultABt(mat_atomvec_chi_y[vecind], vxc_a_w);
                auto vxc_a_gz = sdenblas::serialMultABt(mat_atomvec_chi_z[vecind], vxc_a_w);

                auto vxc_b_gx = sdenblas::serialMultABt(mat_atomvec_chi_x[vecind], vxc_b_w);
                auto vxc_b_gy = sdenblas::serialMultABt(mat_atomvec_chi_y[vecind], vxc_b_w);
                auto vxc_b_gz = sdenblas::serialMultABt(mat_atomvec_chi_z[vecind], vxc_b_w);

                auto vxc_a_gx_2 = sdenblas::serialMultABt(mat_chi, vxc_a_wx);
                auto vxc_a_gy_2 = sdenblas::serialMultABt(mat_chi, vxc_a_wy);
                auto vxc_a_gz_2 = sdenblas::serialMultABt(mat_chi, vxc_a_wz);

                auto vxc_b_gx_2 = sdenblas::serialMultABt(mat_chi, vxc_b_wx);
                auto vxc_b_gy_2 = sdenblas::serialMultABt(mat_chi, vxc_b_wy);
                auto vxc_b_gz_2 = sdenblas::serialMultABt(mat_chi, vxc_b_wz);

                sdenblas::serialInPlaceAddAB(vxc_a_gx, vxc_a_gx_2);
                sdenblas::serialInPlaceAddAB(vxc_a_gy, vxc_a_gy_2);
                sdenblas::serialInPlaceAddAB(vxc_a_gz, vxc_a_gz_2);

                sdenblas::serialInPlaceAddAB(vxc_b_gx, vxc_b_gx_2);
                sdenblas::serialInPlaceAddAB(vxc_b_gy, vxc_b_gy_2);
                sdenblas::serialInPlaceAddAB(vxc_b_gz, vxc_b_gz_2);

                vxc_a_gx.symmetrizeAndScale(0.5);
                vxc_a_gy.symmetrizeAndScale(0.5);
                vxc_a_gz.symmetrizeAndScale(0.5);

                vxc_b_gx.symmetrizeAndScale(0.5);
                vxc_b_gy.symmetrizeAndScale(0.5);
                vxc_b_gz.symmetrizeAndScale(0.5);

                #pragma omp critical
                {
                    sdenblas::serialInPlaceAddAB(sum_partial_vxc_a_gx[vecind], vxc_a_gx);
                    sdenblas::serialInPlaceAddAB(sum_partial_vxc_a_gy[vecind], vxc_a_gy);
                    sdenblas::serialInPlaceAddAB(sum_partial_vxc_a_gz[vecind], vxc_a_gz);

                    sdenblas::serialInPlaceAddAB(sum_partial_vxc_b_gx[vecind], vxc_b_gx);
                    sdenblas::serialInPlaceAddAB(sum_partial_vxc_b_gy[vecind], vxc_b_gy);
                    sdenblas::serialInPlaceAddAB(sum_partial_vxc_b_gz[vecind], vxc_b_gz);
                }
            }

            omptimers[thread_id].stop("VxcFockGrad on atoms");
        }

        timer.start("Vxc grad. dist.");

        for (int vecind = 0; vecind < natoms; vecind++)
        {
            dftsubmat::distributeSubMatrixToDenseMatrix(vxcgrads[vecind * 6 + 0], sum_partial_vxc_a_gx[vecind], aoinds, naos);
            dftsubmat::distributeSubMatrixToDenseMatrix(vxcgrads[vecind * 6 + 1], sum_partial_vxc_a_gy[vecind], aoinds, naos);
            dftsubmat::distributeSubMatrixToDenseMatrix(vxcgrads[vecind * 6 + 2], sum_partial_vxc_a_gz[vecind], aoinds, naos);

            dftsubmat::distributeSubMatrixToDenseMatrix(vxcgrads[vecind * 6 + 3], sum_partial_vxc_b_gx[vecind], aoinds, naos);
            dftsubmat::distributeSubMatrixToDenseMatrix(vxcgrads[vecind * 6 + 4], sum_partial_vxc_b_gy[vecind], aoinds, naos);
            dftsubmat::distributeSubMatrixToDenseMatrix(vxcgrads[vecind * 6 + 5], sum_partial_vxc_b_gz[vecind], aoinds, naos);
        }

        timer.stop("Vxc grad. dist.");
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

    return vxcgrads;
}

}  // namespace xchessgga
