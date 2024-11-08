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

#include "XCMolecularHessianForLDA.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "AODensityMatrix.hpp"
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

namespace xchesslda {  // xchesslda namespace

auto
integrateExcHessianForLDA(const CMolecule&        molecule,
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

    xchesslda::computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.number_of_atoms();

    CDenseMatrix molhess_threads(nthreads, (natoms * 3) * (natoms * 3));

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    // density and functional derivatives

    auto       ldafunc = xcFunctional.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> v2rho2_data(dim->v2rho2 * max_npoints_per_box);

    std::vector<double> weighted_vrho(max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho    = rho_data.data();
    auto vrho   = vrho_data.data();
    auto v2rho2 = v2rho2_data.data();

    auto w0 = weighted_vrho.data();

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    timer.stop("Preparation");

    for (int32_t box_id = 0; box_id < counts.size(); box_id++)
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

        auto gs_sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

        timer.stop("Density matrix slicing");

        dengridgen::generateDensityForLDA(rho, mat_chi, gs_sub_dens_mat, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, npoints);
        CDenseMatrix dengrady(natoms, npoints);
        CDenseMatrix dengradz(natoms, npoints);

        auto gdenx = dengradx.values();
        auto gdeny = dengrady.values();
        auto gdenz = dengradz.values();

        timer.stop("Density grad. grid prep.");

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(gs_sub_dens_mat, mat_chi);

        timer.stop("Density grad. grid matmul");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_lda(npoints, rho, vrho);

        xcFunctional.compute_fxc_for_lda(npoints, rho, v2rho2);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        timer.start("Accumulate Hessian");

        auto D_val = gs_sub_dens_mat.values();

        auto F_val = mat_F.values();

        auto chi_x_val = mat_chi_x.values();
        auto chi_y_val = mat_chi_y.values();
        auto chi_z_val = mat_chi_z.values();

        auto chi_xx_val = mat_chi_xx.values();
        auto chi_xy_val = mat_chi_xy.values();
        auto chi_xz_val = mat_chi_xz.values();
        auto chi_yy_val = mat_chi_yy.values();
        auto chi_yz_val = mat_chi_yz.values();
        auto chi_zz_val = mat_chi_zz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molhess_threads.row(thread_id);

            // prepare w0

            #pragma omp simd
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                w0[g] = local_weights[g] * vrho[2 * g + 0];
            }

            // prepare gradient grid

            for (int32_t mu = 0; mu < aocount; mu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[mu]];

                auto atom_offset = atomidx * npoints;

                auto mu_offset = mu * npoints;

                #pragma omp simd
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto mu_g = mu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[mu_g] * chi_x_val[mu_g];
                    gdeny[atom_g] -= 2.0 * F_val[mu_g] * chi_y_val[mu_g];
                    gdenz[atom_g] -= 2.0 * F_val[mu_g] * chi_z_val[mu_g];
                }
            }

            // first contribution to
            // f_{\rho_{\alpha}} \rho_{\alpha}^{(\xi,\zeta)}
            // on the same atom

            for (int32_t mu = 0; mu < aocount; mu++)
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
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto mu_g = mu_offset + g;

                    // f_{\rho_{\alpha}} \rho_{\alpha}^{(\xi,\zeta)}

                    // \rho_{\alpha}^{(\xi,\zeta)} (first contrib.)
                    // = 2 \sum_{\mu\nu} D^{\alpha}_{\mu\nu,\rm{sym}} \phi_{\mu}^{(\xi,\zeta)} \phi_{\nu}
                    // = 2 \sum_{\mu} F_{\mu} \phi_{\mu}^{(\xi,\zeta)}

                    // factor of 2 is added outside of the for loop

                    gatmxx += w0[g] * F_val[mu_g] * chi_xx_val[mu_g];
                    gatmxy += w0[g] * F_val[mu_g] * chi_xy_val[mu_g];
                    gatmxz += w0[g] * F_val[mu_g] * chi_xz_val[mu_g];

                    gatmyx += w0[g] * F_val[mu_g] * chi_xy_val[mu_g];
                    gatmyy += w0[g] * F_val[mu_g] * chi_yy_val[mu_g];
                    gatmyz += w0[g] * F_val[mu_g] * chi_yz_val[mu_g];

                    gatmzx += w0[g] * F_val[mu_g] * chi_xz_val[mu_g];
                    gatmzy += w0[g] * F_val[mu_g] * chi_yz_val[mu_g];
                    gatmzz += w0[g] * F_val[mu_g] * chi_zz_val[mu_g];
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
            // on the same atom and on different atoms

            for (int32_t mu = 0; mu < aocount; mu++)
            {
                auto iatom = ao_to_atom_ids[aoinds[mu]];

                auto ix = iatom * 3 + 0;
                auto iy = iatom * 3 + 1;
                auto iz = iatom * 3 + 2;

                auto mu_offset = mu * npoints;

                for (int32_t nu = 0; nu < aocount; nu++)
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
                    for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                    {
                        auto mu_g = mu_offset + g;
                        auto nu_g = nu_offset + g;

                        // f_{\rho_{\alpha}} \rho_{\alpha}^{(\xi,\zeta)}

                        // \rho_{\alpha}^{(\xi,\zeta)} (second contrib.)
                        // = 2 \sum_{\mu\nu} D^{\alpha}_{\mu\nu,\rm{sym}} \phi_{\mu}^{(\xi)} \phi_{\nu}^{(\zeta)}

                        // factor of 2 and D_{\mu\nu,sym}^{\alpha} are added outside of the for loop

                        gatmxx += w0[g] * chi_x_val[mu_g] * chi_x_val[nu_g];
                        gatmxy += w0[g] * chi_x_val[mu_g] * chi_y_val[nu_g];
                        gatmxz += w0[g] * chi_x_val[mu_g] * chi_z_val[nu_g];

                        gatmyx += w0[g] * chi_y_val[mu_g] * chi_x_val[nu_g];
                        gatmyy += w0[g] * chi_y_val[mu_g] * chi_y_val[nu_g];
                        gatmyz += w0[g] * chi_y_val[mu_g] * chi_z_val[nu_g];

                        gatmzx += w0[g] * chi_z_val[mu_g] * chi_x_val[nu_g];
                        gatmzy += w0[g] * chi_z_val[mu_g] * chi_y_val[nu_g];
                        gatmzz += w0[g] * chi_z_val[mu_g] * chi_z_val[nu_g];
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

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto ix = iatom * 3 + 0;
                auto iy = iatom * 3 + 1;
                auto iz = iatom * 3 + 2;

                auto i_offset = iatom * npoints;

                for (int32_t jatom = iatom; jatom < natoms; jatom++)
                {
                    auto jx = jatom * 3 + 0;
                    auto jy = jatom * 3 + 1;
                    auto jz = jatom * 3 + 2;

                    auto j_offset = jatom * npoints;

                    double gatmxx = 0.0, gatmyx = 0.0, gatmzx = 0.0;
                    double gatmxy = 0.0, gatmyy = 0.0, gatmzy = 0.0;
                    double gatmxz = 0.0, gatmyz = 0.0, gatmzz = 0.0;

                    #pragma omp simd reduction(+ : gatmxx, gatmxy, gatmxz, gatmyx, gatmyy, gatmyz, gatmzx, gatmzy, gatmzz)
                    for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                    {
                        auto ig = i_offset + g;
                        auto jg = j_offset + g;

                        // (f_{\rho_{\alpha} \rho_{\alpha}} + f_{\rho_{\alpha} \rho_{\beta}})
                        // \rho_{\alpha}^{(\xi)} \rho_{\alpha}^{(\zeta)}

                        double prefac = local_weights[g] * (v2rho2[3 * g + 0] + v2rho2[3 * g + 1]);

                        gatmxx += prefac * gdenx[ig] * gdenx[jg];
                        gatmxy += prefac * gdenx[ig] * gdeny[jg];
                        gatmxz += prefac * gdenx[ig] * gdenz[jg];

                        gatmyx += prefac * gdeny[ig] * gdenx[jg];
                        gatmyy += prefac * gdeny[ig] * gdeny[jg];
                        gatmyz += prefac * gdeny[ig] * gdenz[jg];

                        gatmzx += prefac * gdenz[ig] * gdenx[jg];
                        gatmzy += prefac * gdenz[ig] * gdeny[jg];
                        gatmzz += prefac * gdenz[ig] * gdenz[jg];
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

    for (int32_t iatom = 0; iatom < natoms; iatom++)
    {
        auto ix = iatom * 3 + 0;
        auto iy = iatom * 3 + 1;
        auto iz = iatom * 3 + 2;

        for (int32_t jatom = iatom; jatom < natoms; jatom++)
        {
            auto jx = jatom * 3 + 0;
            auto jy = jatom * 3 + 1;
            auto jz = jatom * 3 + 2;

            for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
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
integrateVxcFockGradientForLDA(const CMolecule&        molecule,
                               const CMolecularBasis&  basis,
                               const std::vector<const double*>& gsDensityPointers,
                               const CMolecularGrid&   molecularGrid,
                               const double            screeningThresholdForGTOValues,
                               const CXCFunctional&    xcFunctional,
                               const int32_t           atomIdx) -> std::vector<CDenseMatrix>
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

    xchesslda::computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // Vxc Fock gradeints (in x,y,z directions)

    std::vector<CDenseMatrix> vxcgrads;

    vxcgrads.push_back(CDenseMatrix(naos, naos));
    vxcgrads.push_back(CDenseMatrix(naos, naos));
    vxcgrads.push_back(CDenseMatrix(naos, naos));

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    // density and functional derivatives

    auto       ldafunc = xcFunctional.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> v2rho2_data(dim->v2rho2 * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho    = rho_data.data();
    auto vrho   = vrho_data.data();
    auto v2rho2 = v2rho2_data.data();

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    timer.stop("Preparation");

    for (int32_t box_id = 0; box_id < counts.size(); box_id++)
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

        // GTO values on grid points

        timer.start("OMP GTO evaluation");

        CDenseMatrix mat_chi(aocount, npoints);

        CDenseMatrix mat_chi_x(aocount, npoints);
        CDenseMatrix mat_chi_y(aocount, npoints);
        CDenseMatrix mat_chi_z(aocount, npoints);

        CDenseMatrix mat_atom_chi_x(aocount, npoints);
        CDenseMatrix mat_atom_chi_y(aocount, npoints);
        CDenseMatrix mat_atom_chi_z(aocount, npoints);

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

                    auto iatom = ao_to_atom_ids[pre_ao_inds[nu]];

                    if (iatom == atomIdx)
                    {
                        std::memcpy(mat_atom_chi_x.row(idx) + grid_batch_offset, submat_x_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                        std::memcpy(mat_atom_chi_y.row(idx) + grid_batch_offset, submat_y_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                        std::memcpy(mat_atom_chi_z.row(idx) + grid_batch_offset, submat_z_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    }
                }
            }

            omptimers[thread_id].stop("gtoeval");
        }

        timer.stop("OMP GTO evaluation");

        // generate sub density matrix and density grid

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

        timer.stop("Density matrix slicing");

        dengridgen::generateDensityForLDA(rho, mat_chi, gs_sub_dens_mat, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(1, npoints);
        CDenseMatrix dengrady(1, npoints);
        CDenseMatrix dengradz(1, npoints);

        auto gdenx = dengradx.values();
        auto gdeny = dengrady.values();
        auto gdenz = dengradz.values();

        timer.stop("Density grad. grid prep.");

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(gs_sub_dens_mat, mat_chi);

        timer.stop("Density grad. grid matmul");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_lda(npoints, rho, vrho);

        xcFunctional.compute_fxc_for_lda(npoints, rho, v2rho2);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // compute partial contribution to Vxc Fock gradient and distribute partial
        // Vxc Fock gradient to full matrix

        timer.start("Vxc matrix G");

        auto F_val = mat_F.values();

        auto chi_val = mat_chi.values();

        auto chi_x_val = mat_chi_x.values();
        auto chi_y_val = mat_chi_y.values();
        auto chi_z_val = mat_chi_z.values();

        CDenseMatrix vxc_w(aocount, npoints);
        CDenseMatrix vxc_wx(aocount, npoints);
        CDenseMatrix vxc_wy(aocount, npoints);
        CDenseMatrix vxc_wz(aocount, npoints);

        auto vxc_w_val  = vxc_w.values();
        auto vxc_wx_val = vxc_wx.values();
        auto vxc_wy_val = vxc_wy.values();
        auto vxc_wz_val = vxc_wz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto nthreads = omp_get_max_threads();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            // prepare gradient density

            for (int32_t nu = 0; nu < aocount; nu++)
            {
                auto iatom = ao_to_atom_ids[aoinds[nu]];

                if (iatom != atomIdx) continue;

                auto nu_offset = nu * npoints;

                #pragma omp simd
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto nu_g = nu_offset + g;

                    gdenx[g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];
                    gdeny[g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];
                    gdenz[g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];
                }
            }

            // Vxc matrix element gradient

            for (int32_t nu = 0; nu < aocount; nu++)
            {
                auto nu_offset = nu * npoints;

                #pragma omp simd
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto nu_g = nu_offset + g;

                    auto prefac = local_weights[g] * vrho[2 * g + 0];

                    // 2 f_{\rho_{\alpha}} \phi_{\mu}^{(\xi)} \phi_{\nu}

                    // note: \phi_{\mu}^{(\xi)} will be added later (from mat_atom_chi_{xyz})

                    vxc_w_val[nu_g] = -2.0 * prefac * chi_val[nu_g];

                    // (f_{\rho_{\alpha} \rho_{\alpha}} + f_{\rho_{\alpha} \rho_{\beta}})
                    // \rho_{\alpha}^{(\xi)} \phi_{\mu} \phi_{\nu}

                    // note: \phi_{\mu} will be added later (from mat_chi)

                    prefac = local_weights[g] * (v2rho2[3 * g + 0] + v2rho2[3 * g + 1]) * chi_val[nu_g];

                    vxc_wx_val[nu_g] = prefac * gdenx[g];
                    vxc_wy_val[nu_g] = prefac * gdeny[g];
                    vxc_wz_val[nu_g] = prefac * gdenz[g];
                }
            }
        }

        timer.stop("Vxc matrix G");

        timer.start("Vxc matrix matmul");

        auto vxc_gx_first_contrib = denblas::multABt(mat_atom_chi_x, vxc_w);
        auto vxc_gy_first_contrib = denblas::multABt(mat_atom_chi_y, vxc_w);
        auto vxc_gz_first_contrib = denblas::multABt(mat_atom_chi_z, vxc_w);

        auto vxc_gx_second_contrib = denblas::multABt(mat_chi, vxc_wx);
        auto vxc_gy_second_contrib = denblas::multABt(mat_chi, vxc_wy);
        auto vxc_gz_second_contrib = denblas::multABt(mat_chi, vxc_wz);

        auto vxc_gx = denblas::addAB(vxc_gx_first_contrib, vxc_gx_second_contrib, 1.0);
        auto vxc_gy = denblas::addAB(vxc_gy_first_contrib, vxc_gy_second_contrib, 1.0);
        auto vxc_gz = denblas::addAB(vxc_gz_first_contrib, vxc_gz_second_contrib, 1.0);

        vxc_gx.symmetrizeAndScale(0.5);
        vxc_gy.symmetrizeAndScale(0.5);
        vxc_gz.symmetrizeAndScale(0.5);

        timer.stop("Vxc matrix matmul");

        timer.start("Vxc matrix dist.");

        dftsubmat::distributeSubMatrixToDenseMatrix(vxcgrads[0], vxc_gx, aoinds, naos);
        dftsubmat::distributeSubMatrixToDenseMatrix(vxcgrads[1], vxc_gy, aoinds, naos);
        dftsubmat::distributeSubMatrixToDenseMatrix(vxcgrads[2], vxc_gz, aoinds, naos);

        timer.stop("Vxc matrix dist.");
    }

    timer.stop("Total timing");

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;
    // std::cout << "OpenMP timing" << std::endl;
    // for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
    //{
    //     std::cout << "Thread " << thread_id << std::endl;
    //     std::cout << omptimers[thread_id].getSummary() << std::endl;
    // }

    return vxcgrads;
}

auto
computeAOtoAtomMapping(std::vector<int>& ao_to_atom_ids, const CMolecule& molecule, const CMolecularBasis& basis) -> void
{
    auto natoms = molecule.number_of_atoms();

    auto max_angl = basis.max_angular_momentum();

    // azimuthal quantum number: s,p,d,f,...

    for (int angl = 0, aoidx = 0; angl <= max_angl; angl++)
    {
        auto nsph = angl * 2 + 1;

        // magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...

        for (int isph = 0; isph < nsph; isph++)
        {
            // atoms

            for (int atomidx = 0; atomidx < natoms; atomidx++)
            {
                auto nao = basis.number_of_basis_functions(std::vector<int>({atomidx}), angl);

                // atomic orbitals

                for (int iao = 0; iao < nao; iao++, aoidx++)
                {
                    ao_to_atom_ids[aoidx] = atomidx;
                }
            }
        }
    }
}

}  // namespace xchesslda
