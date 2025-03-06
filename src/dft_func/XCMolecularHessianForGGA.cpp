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

namespace xchessgga {  // xchessgga namespace

auto
integrateExcHessianForGGA(const CMolecule&        molecule,
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

    // density and functional derivatives

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhograd_data(dim->rho * 3 * max_npoints_per_box);
    std::vector<double> sigma_data(dim->sigma * max_npoints_per_box);

    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> vsigma_data(dim->vsigma * max_npoints_per_box);

    std::vector<double> v2rho2_data(dim->v2rho2 * max_npoints_per_box);
    std::vector<double> v2rhosigma_data(dim->v2rhosigma * max_npoints_per_box);
    std::vector<double> v2sigma2_data(dim->v2sigma2 * max_npoints_per_box);

    std::vector<double> weighted_vrho(max_npoints_per_box);
    std::vector<double> weighted_vnabla_x(max_npoints_per_box);
    std::vector<double> weighted_vnabla_y(max_npoints_per_box);
    std::vector<double> weighted_vnabla_z(max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();

    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();

    auto v2rho2     = v2rho2_data.data();
    auto v2rhosigma = v2rhosigma_data.data();
    auto v2sigma2   = v2sigma2_data.data();

    auto w0 = weighted_vrho.data();
    auto wx = weighted_vnabla_x.data();
    auto wy = weighted_vnabla_y.data();
    auto wz = weighted_vnabla_z.data();

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

        CDenseMatrix mat_chi_xxx(aocount, npoints);
        CDenseMatrix mat_chi_xxy(aocount, npoints);
        CDenseMatrix mat_chi_xxz(aocount, npoints);
        CDenseMatrix mat_chi_xyy(aocount, npoints);
        CDenseMatrix mat_chi_xyz(aocount, npoints);
        CDenseMatrix mat_chi_xzz(aocount, npoints);
        CDenseMatrix mat_chi_yyy(aocount, npoints);
        CDenseMatrix mat_chi_yyz(aocount, npoints);
        CDenseMatrix mat_chi_yzz(aocount, npoints);
        CDenseMatrix mat_chi_zzz(aocount, npoints);

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

                    std::memcpy(mat_chi_xxx.row(idx) + grid_batch_offset, submat_xxx_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xxy.row(idx) + grid_batch_offset, submat_xxy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xxz.row(idx) + grid_batch_offset, submat_xxz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xyy.row(idx) + grid_batch_offset, submat_xyy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xyz.row(idx) + grid_batch_offset, submat_xyz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xzz.row(idx) + grid_batch_offset, submat_xzz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yyy.row(idx) + grid_batch_offset, submat_yyy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yyz.row(idx) + grid_batch_offset, submat_yyz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yzz.row(idx) + grid_batch_offset, submat_yzz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_zzz.row(idx) + grid_batch_offset, submat_zzz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");
        }

        timer.stop("OMP GTO evaluation");

        // generate sub density matrix and density grid

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

        timer.stop("Density matrix slicing");

        dengridgen::generateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat, timer);

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

        timer.stop("Density grad. grid prep.");

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(gs_sub_dens_mat, mat_chi);

        auto mat_F_x = denblas::multAB(gs_sub_dens_mat, mat_chi_x);
        auto mat_F_y = denblas::multAB(gs_sub_dens_mat, mat_chi_y);
        auto mat_F_z = denblas::multAB(gs_sub_dens_mat, mat_chi_z);

        timer.stop("Density grad. grid matmul");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

        xcFunctional.compute_fxc_for_gga(npoints, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        timer.start("Accumulate Hessian");

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

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molhess_threads.row(thread_id);

            // prepare w0, wx, wy and wz

            #pragma omp simd 
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

                auto atom_offset = atomidx * npoints;

                auto mu_offset = mu * npoints;

                #pragma omp simd 
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

                auto mu_offset = mu * npoints;

                double gatmxx = 0.0, gatmxy = 0.0, gatmxz = 0.0;
                double gatmyx = 0.0, gatmyy = 0.0, gatmyz = 0.0;
                double gatmzx = 0.0, gatmzy = 0.0, gatmzz = 0.0;

                #pragma omp simd reduction(+ : gatmxx, gatmxy, gatmxz, gatmyx, gatmyy, gatmyz, gatmzx, gatmzy, gatmzz) 
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

                auto mu_offset = mu * npoints;

                for (int nu = 0; nu < aocount; nu++)
                {
                    auto jatom = ao_to_atom_ids[aoinds[nu]];

                    // only consider the upper triangular part, e.g. iatom <= jatom

                    if (iatom > jatom) continue;

                    auto jx = jatom * 3 + 0;
                    auto jy = jatom * 3 + 1;
                    auto jz = jatom * 3 + 2;

                    auto nu_offset = nu * npoints;

                    double gatmxx = 0.0, gatmxy = 0.0, gatmxz = 0.0;
                    double gatmyx = 0.0, gatmyy = 0.0, gatmyz = 0.0;
                    double gatmzx = 0.0, gatmzy = 0.0, gatmzz = 0.0;

                    #pragma omp simd reduction(+ : gatmxx, gatmxy, gatmxz, gatmyx, gatmyy, gatmyz, gatmzx, gatmzy, gatmzz) 
                    for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

                auto i_offset = iatom * npoints;

                for (int jatom = iatom; jatom < natoms; jatom++)
                {
                    auto jx = jatom * 3 + 0;
                    auto jy = jatom * 3 + 1;
                    auto jz = jatom * 3 + 2;

                    auto j_offset = jatom * npoints;

                    double gatmxx = 0.0, gatmyx = 0.0, gatmzx = 0.0;
                    double gatmxy = 0.0, gatmyy = 0.0, gatmzy = 0.0;
                    double gatmxz = 0.0, gatmyz = 0.0, gatmzz = 0.0;

                    #pragma omp simd reduction(+ : gatmxx, gatmxy, gatmxz, gatmyx, gatmyy, gatmyz, gatmzx, gatmzy, gatmzz) 
                    for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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
        }

        timer.stop("Accumulate Hessian");
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
integrateVxcFockGradientForGGA(const CMolecule&        molecule,
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

    // density and functional derivatives

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhograd_data(dim->rho * 3 * max_npoints_per_box);
    std::vector<double> sigma_data(dim->sigma * max_npoints_per_box);

    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> vsigma_data(dim->vsigma * max_npoints_per_box);

    std::vector<double> v2rho2_data(dim->v2rho2 * max_npoints_per_box);
    std::vector<double> v2rhosigma_data(dim->v2rhosigma * max_npoints_per_box);
    std::vector<double> v2sigma2_data(dim->v2sigma2 * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();

    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();

    auto v2rho2     = v2rho2_data.data();
    auto v2rhosigma = v2rhosigma_data.data();
    auto v2sigma2   = v2sigma2_data.data();

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

        std::vector<CDenseMatrix> mat_atomvec_chi_x(natoms, CDenseMatrix(aocount, npoints));
        std::vector<CDenseMatrix> mat_atomvec_chi_y(natoms, CDenseMatrix(aocount, npoints));
        std::vector<CDenseMatrix> mat_atomvec_chi_z(natoms, CDenseMatrix(aocount, npoints));

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

                    auto iatom = ao_to_atom_ids[pre_ao_inds[nu]];

                    for (int vecind = 0; vecind < natoms; vecind++)
                    {
                        if (atomIdxVec[vecind] == iatom)
                        {
                            std::memcpy(mat_atomvec_chi_x[vecind].row(idx) + grid_batch_offset, submat_x_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                            std::memcpy(mat_atomvec_chi_y[vecind].row(idx) + grid_batch_offset, submat_y_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                            std::memcpy(mat_atomvec_chi_z[vecind].row(idx) + grid_batch_offset, submat_z_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                        }
                    }

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

        auto gs_sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

        timer.stop("Density matrix slicing");

        dengridgen::generateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat, timer);

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(gs_sub_dens_mat, mat_chi);

        auto mat_F_x = denblas::multAB(gs_sub_dens_mat, mat_chi_x);
        auto mat_F_y = denblas::multAB(gs_sub_dens_mat, mat_chi_y);
        auto mat_F_z = denblas::multAB(gs_sub_dens_mat, mat_chi_z);

        timer.stop("Density grad. grid matmul");

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

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

        xcFunctional.compute_fxc_for_gga(npoints, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        timer.start("Parallel Vxc grad.");

#pragma omp parallel for schedule(static)
        for (int vecind = 0; vecind < natoms; vecind++)
        {
            auto atomIdx = atomIdxVec[vecind];

            CDenseMatrix dengradx(1, npoints);
            CDenseMatrix dengrady(1, npoints);
            CDenseMatrix dengradz(1, npoints);

            CDenseMatrix dengradxx(1, npoints);
            CDenseMatrix dengradxy(1, npoints);
            CDenseMatrix dengradxz(1, npoints);

            CDenseMatrix dengradyx(1, npoints);
            CDenseMatrix dengradyy(1, npoints);
            CDenseMatrix dengradyz(1, npoints);

            CDenseMatrix dengradzx(1, npoints);
            CDenseMatrix dengradzy(1, npoints);
            CDenseMatrix dengradzz(1, npoints);

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

            CDenseMatrix vxc_w(aocount, npoints);
            CDenseMatrix vxc_wx(aocount, npoints);
            CDenseMatrix vxc_wy(aocount, npoints);
            CDenseMatrix vxc_wz(aocount, npoints);

            auto vxc_w_val  = vxc_w.values();
            auto vxc_wx_val = vxc_wx.values();
            auto vxc_wy_val = vxc_wy.values();
            auto vxc_wz_val = vxc_wz.values();

            // prepare gradient density

            for (int nu = 0; nu < aocount; nu++)
            {
                auto iatom = ao_to_atom_ids[aoinds[nu]];

                if (iatom != atomIdx) continue;

                auto nu_offset = nu * npoints;

                #pragma omp simd 
                for (int g = 0; g < npoints; g++)
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

                auto nu_offset = nu * npoints;

                #pragma omp simd 
                for (int g = 0; g < npoints; g++)
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

            auto vxc_gx_first_contrib = denblas::serialMultABt(mat_atomvec_chi_x[vecind], vxc_w);
            auto vxc_gy_first_contrib = denblas::serialMultABt(mat_atomvec_chi_y[vecind], vxc_w);
            auto vxc_gz_first_contrib = denblas::serialMultABt(mat_atomvec_chi_z[vecind], vxc_w);

            auto vxc_gx_second_contrib = denblas::serialMultABt(mat_chi, vxc_wx);
            auto vxc_gy_second_contrib = denblas::serialMultABt(mat_chi, vxc_wy);
            auto vxc_gz_second_contrib = denblas::serialMultABt(mat_chi, vxc_wz);

            auto vxc_gx = denblas::addAB(vxc_gx_first_contrib, vxc_gx_second_contrib, 1.0);
            auto vxc_gy = denblas::addAB(vxc_gy_first_contrib, vxc_gy_second_contrib, 1.0);
            auto vxc_gz = denblas::addAB(vxc_gz_first_contrib, vxc_gz_second_contrib, 1.0);

            vxc_gx.symmetrizeAndScale(0.5);
            vxc_gy.symmetrizeAndScale(0.5);
            vxc_gz.symmetrizeAndScale(0.5);

            dftsubmat::distributeSubMatrixToDenseMatrix(vxcgrads[vecind * 3 + 0], vxc_gx, aoinds, naos);
            dftsubmat::distributeSubMatrixToDenseMatrix(vxcgrads[vecind * 3 + 1], vxc_gy, aoinds, naos);
            dftsubmat::distributeSubMatrixToDenseMatrix(vxcgrads[vecind * 3 + 2], vxc_gz, aoinds, naos);

        }

        timer.stop("Parallel Vxc grad.");
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
integrateVxcFockGradientForGGA(const CMolecule&        molecule,
                               const CMolecularBasis&  basis,
                               const std::vector<const double*>& gsDensityPointers,
                               const CMolecularGrid&   molecularGrid,
                               const double            screeningThresholdForGTOValues,
                               const CXCFunctional&    xcFunctional,
                               const int               atomIdx) -> std::vector<CDenseMatrix>
{
    std::vector<int> atomIdxVec({atomIdx});

    return integrateVxcFockGradientForGGA(molecule,
                                          basis,
                                          gsDensityPointers,
                                          molecularGrid,
                                          screeningThresholdForGTOValues,
                                          xcFunctional,
                                          atomIdxVec);
}

}  // namespace xchessgga
