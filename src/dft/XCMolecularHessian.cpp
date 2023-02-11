//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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

#include "XCMolecularHessian.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "AODensityMatrix.hpp"
#include "AngularMomentum.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DensityGridGenerator.hpp"
#include "DensityGridType.hpp"
#include "GridScreener.hpp"
#include "GtoEvaluator.hpp"
#include "NewFunctionalParser.hpp"
#include "SubMatrix.hpp"
#include "XCFuncType.hpp"

CXCMolecularHessian::CXCMolecularHessian(MPI_Comm comm)

    : _screeningThresholdForGTOValues(1.0e-12)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;
}

CDenseMatrix
CXCMolecularHessian::integrateVxcHessian(const CMolecule&        molecule,
                                         const CMolecularBasis&  basis,
                                         const CAODensityMatrix& gsDensityMatrix,
                                         const CMolecularGrid&   molecularGrid,
                                         const std::string&      xcFuncLabel) const
{
    auto newfvxc = newvxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = newfvxc.getFunctionalType();

    if (gsDensityMatrix.isClosedShell())
    {
        if (xcfuntype == xcfun::lda)
        {
            return _integrateVxcHessianForLDA(molecule, basis, gsDensityMatrix, molecularGrid, newfvxc);
        }
        else
        {
            std::string errxcfuntype("XCMolecularHessian.integrateVxcHessian: Only implemented for LDA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCMolecularHessian.integrateVxcHessian: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }

    return CDenseMatrix();
}

CDenseMatrix
CXCMolecularHessian::integrateFxcHessian(const CMolecule&        molecule,
                                         const CMolecularBasis&  basis,
                                         const CAODensityMatrix& gsDensityMatrix,
                                         const CMolecularGrid&   molecularGrid,
                                         const std::string&      xcFuncLabel) const
{
    auto newfvxc = newvxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = newfvxc.getFunctionalType();

    if (gsDensityMatrix.isClosedShell())
    {
        if (xcfuntype == xcfun::lda)
        {
            return _integrateFxcHessianForLDA(molecule, basis, gsDensityMatrix, molecularGrid, newfvxc);
        }
        else
        {
            std::string errxcfuntype("XCMolecularHessian.integrateFxcHessian: Only implemented for LDA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCMolecularHessian.integrateFxcHessian: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }

    return CDenseMatrix();
}

CDenseMatrix
CXCMolecularHessian::_integrateVxcHessianForLDA(const CMolecule&        molecule,
                                                const CMolecularBasis&  basis,
                                                const CAODensityMatrix& gsDensityMatrix,
                                                const CMolecularGrid&   molecularGrid,
                                                const CXCNewFunctional& xcFunctional) const
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // AO-to-atom mapping

    std::vector<int32_t> ao_to_atom_ids(naos);

    _computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.getNumberOfAtoms();

    CMemBlock2D<double> molhess_threads((natoms * 3) * (natoms * 3), nthreads);

    // memory blocks for GTOs on grid points

    CMemBlock2D<double> gaos(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaox(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);
    CMemBlock2D<double> gaoy(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);
    CMemBlock2D<double> gaoz(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaoxx(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);
    CMemBlock2D<double> gaoxy(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);
    CMemBlock2D<double> gaoxz(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);
    CMemBlock2D<double> gaoyy(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);
    CMemBlock2D<double> gaoyz(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);
    CMemBlock2D<double> gaozz(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos);  // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho  = rho_data.data();
    auto vrho = vrho_data.data();

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

        auto boxdim = gtoeval::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // pre-screening of GTOs

        timer.start("GTO pre-screening");

        gtoeval::preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 2, _screeningThresholdForGTOValues, boxdim);  // 2nd order GTO derivative

        timer.stop("GTO pre-screening");

        // GTO values on grid points

        timer.start("OMP GTO evaluation");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForMetaGGA(
                gaos, gaox, gaoy, gaoz, gaoxx, gaoxy, gaoxz, gaoyy, gaoyz, gaozz, gtovec, xcoords, ycoords, zcoords,
                gridblockpos, grid_batch_offset, grid_batch_size, skip_cgto_ids);
        }

        timer.stop("OMP GTO evaluation");

        timer.start("GTO screening");

        int32_t aocount = 0;

        for (int32_t nu = 0; nu < naos; nu++)
        {
            if (skip_ao_ids.data()[nu]) continue;

            bool skip = true;

            auto gaos_nu = gaos.data(nu);

            auto gaox_nu = gaox.data(nu);
            auto gaoy_nu = gaoy.data(nu);
            auto gaoz_nu = gaoz.data(nu);

            auto gaoxx_nu = gaoxx.data(nu);
            auto gaoxy_nu = gaoxy.data(nu);
            auto gaoxz_nu = gaoxz.data(nu);
            auto gaoyy_nu = gaoyy.data(nu);
            auto gaoyz_nu = gaoyz.data(nu);
            auto gaozz_nu = gaozz.data(nu);

            for (int32_t g = 0; g < npoints; g++)
            {
                if ((std::fabs(gaos_nu[g]) > _screeningThresholdForGTOValues) ||
                    (std::fabs(gaox_nu[g]) > _screeningThresholdForGTOValues) ||
                    (std::fabs(gaoy_nu[g]) > _screeningThresholdForGTOValues) ||
                    (std::fabs(gaoz_nu[g]) > _screeningThresholdForGTOValues) ||
                    (std::fabs(gaoxx_nu[g]) > _screeningThresholdForGTOValues) ||
                    (std::fabs(gaoxy_nu[g]) > _screeningThresholdForGTOValues) ||
                    (std::fabs(gaoxz_nu[g]) > _screeningThresholdForGTOValues) ||
                    (std::fabs(gaoyy_nu[g]) > _screeningThresholdForGTOValues) ||
                    (std::fabs(gaoyz_nu[g]) > _screeningThresholdForGTOValues) ||
                    (std::fabs(gaozz_nu[g]) > _screeningThresholdForGTOValues))
                {
                    skip = false;

                    break;
                }
            }

            if (!skip)
            {
                aoinds[aocount] = nu;

                ++aocount;
            }
        }

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

        for (int32_t i = 0; i < aocount; i++)
        {
            std::memcpy(mat_chi.row(i), gaos.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_x.row(i), gaox.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_y.row(i), gaoy.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_z.row(i), gaoz.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_xx.row(i), gaoxx.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_xy.row(i), gaoxy.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_xz.row(i), gaoxz.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_yy.row(i), gaoyy.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_yz.row(i), gaoyz.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_zz.row(i), gaozz.data(aoinds[i]), npoints * sizeof(double));
        }

        timer.stop("GTO screening");

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForLDA(rho, npoints, mat_chi, gs_sub_dens_mat, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        // pairwise hessian contribution: only store the upper triangular part (iatom <= jatom)

        auto n2_upper_triang = natoms * (natoms + 1) / 2;

        CDenseMatrix dengradxx(n2_upper_triang, npoints);
        CDenseMatrix dengradxy(n2_upper_triang, npoints);
        CDenseMatrix dengradxz(n2_upper_triang, npoints);

        CDenseMatrix dengradyx(n2_upper_triang, npoints);
        CDenseMatrix dengradyy(n2_upper_triang, npoints);
        CDenseMatrix dengradyz(n2_upper_triang, npoints);

        CDenseMatrix dengradzx(n2_upper_triang, npoints);
        CDenseMatrix dengradzy(n2_upper_triang, npoints);
        CDenseMatrix dengradzz(n2_upper_triang, npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto D_val = gs_sub_dens_mat.values();

        auto mat_F = denblas::multAB(gs_sub_dens_mat, mat_chi);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto naos = mat_chi.getNumberOfRows();

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

        auto gdenxx = dengradxx.values();
        auto gdenxy = dengradxy.values();
        auto gdenxz = dengradxz.values();

        auto gdenyx = dengradyx.values();
        auto gdenyy = dengradyy.values();
        auto gdenyz = dengradyz.values();

        auto gdenzx = dengradzx.values();
        auto gdenzy = dengradzy.values();
        auto gdenzz = dengradzz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            for (int32_t mu = 0; mu < naos; mu++)
            {
                auto iatom = ao_to_atom_ids[aoinds[mu]];

                // upper triangular index for diagonal i: i * N + i - i * (i + 1) / 2

                auto ii_upper_triang_index = iatom * natoms + iatom - iatom * (iatom + 1) / 2;

                auto ii_offset = ii_upper_triang_index * npoints;

                auto mu_offset = mu * npoints;

                #pragma omp simd aligned(F_val, chi_xx_val, chi_xy_val, chi_xz_val, chi_yy_val, chi_yz_val, chi_zz_val, \
                        gdenxx, gdenxy, gdenxz, gdenyx, gdenyy, gdenyz, gdenzx, gdenzy, gdenzz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto ii_g = ii_offset + g;

                    auto mu_g = mu_offset + g;

                    gdenxx[ii_g] += 2.0 * (F_val[mu_g] * chi_xx_val[mu_g]);
                    gdenxy[ii_g] += 2.0 * (F_val[mu_g] * chi_xy_val[mu_g]);
                    gdenxz[ii_g] += 2.0 * (F_val[mu_g] * chi_xz_val[mu_g]);

                    gdenyx[ii_g] += 2.0 * (F_val[mu_g] * chi_xy_val[mu_g]);
                    gdenyy[ii_g] += 2.0 * (F_val[mu_g] * chi_yy_val[mu_g]);
                    gdenyz[ii_g] += 2.0 * (F_val[mu_g] * chi_yz_val[mu_g]);

                    gdenzx[ii_g] += 2.0 * (F_val[mu_g] * chi_xz_val[mu_g]);
                    gdenzy[ii_g] += 2.0 * (F_val[mu_g] * chi_yz_val[mu_g]);
                    gdenzz[ii_g] += 2.0 * (F_val[mu_g] * chi_zz_val[mu_g]);
                }

                for (int32_t nu = 0; nu < naos; nu++)
                {
                    auto jatom = ao_to_atom_ids[aoinds[nu]];

                    // only consider the upper triangular part, e.g. iatom <= jatom

                    if (iatom > jatom) continue;

                    // upper triangular index for i<=j: i * N + j - i * (i + 1) / 2

                    auto ij_upper_triang_index = iatom * natoms + jatom - iatom * (iatom + 1) / 2;

                    auto ij_offset = ij_upper_triang_index * npoints;

                    auto nu_offset = nu * npoints;

                    auto D_mn = D_val[mu * naos + nu];

                    #pragma omp simd aligned(chi_x_val, chi_y_val, chi_z_val, \
                            gdenxx, gdenxy, gdenxz, gdenyx, gdenyy, gdenyz, gdenzx, gdenzy, gdenzz : VLX_ALIGN)
                    for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                    {
                        auto ij_g = ij_offset + g;

                        auto mu_g = mu_offset + g;
                        auto nu_g = nu_offset + g;

                        gdenxx[ij_g] += 2.0 * (chi_x_val[mu_g] * chi_x_val[nu_g] * D_mn);
                        gdenxy[ij_g] += 2.0 * (chi_x_val[mu_g] * chi_y_val[nu_g] * D_mn);
                        gdenxz[ij_g] += 2.0 * (chi_x_val[mu_g] * chi_z_val[nu_g] * D_mn);

                        gdenyx[ij_g] += 2.0 * (chi_y_val[mu_g] * chi_x_val[nu_g] * D_mn);
                        gdenyy[ij_g] += 2.0 * (chi_y_val[mu_g] * chi_y_val[nu_g] * D_mn);
                        gdenyz[ij_g] += 2.0 * (chi_y_val[mu_g] * chi_z_val[nu_g] * D_mn);

                        gdenzx[ij_g] += 2.0 * (chi_z_val[mu_g] * chi_x_val[nu_g] * D_mn);
                        gdenzy[ij_g] += 2.0 * (chi_z_val[mu_g] * chi_y_val[nu_g] * D_mn);
                        gdenzz[ij_g] += 2.0 * (chi_z_val[mu_g] * chi_z_val[nu_g] * D_mn);
                    }
                }
            }
        }

        timer.stop("Density grad. grid rho");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_lda(npoints, rho, vrho);

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        timer.stop("XC functional eval.");

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molhess_threads.data(thread_id);

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

                    // upper triangular index for i<=j: i * N + j - i * (i + 1) / 2

                    auto ij_upper_triang_index = iatom * natoms + jatom - iatom * (iatom + 1) / 2;

                    auto ij_offset = ij_upper_triang_index * npoints;

                    double gatmxx = 0.0, gatmxy = 0.0, gatmxz = 0.0;
                    double gatmyx = 0.0, gatmyy = 0.0, gatmyz = 0.0;
                    double gatmzx = 0.0, gatmzy = 0.0, gatmzz = 0.0;

                    #pragma omp simd reduction(+ : gatmxx, gatmxy, gatmxz, gatmyx, gatmyy, gatmyz, gatmzx, gatmzy, gatmzz) aligned(local_weights, \
                            vrho, gdenxx, gdenxy, gdenxz, gdenyx, gdenyy, gdenyz, gdenzx, gdenzy, gdenzz : VLX_ALIGN)
                    for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                    {
                        auto ij_g = ij_offset + g;

                        double prefac = local_weights[g] * vrho[2 * g + 0];

                        gatmxx += prefac * gdenxx[ij_g];
                        gatmxy += prefac * gdenxy[ij_g];
                        gatmxz += prefac * gdenxz[ij_g];

                        gatmyx += prefac * gdenyx[ij_g];
                        gatmyy += prefac * gdenyy[ij_g];
                        gatmyz += prefac * gdenyz[ij_g];

                        gatmzx += prefac * gdenzx[ij_g];
                        gatmzy += prefac * gdenzy[ij_g];
                        gatmzz += prefac * gdenzz[ij_g];
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

        timer.stop("Accumulate gradient");
    }

    // destroy GTOs container

    delete gtovec;

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
                molhess.row(ix)[jx] += molhess_threads.data(thread_id)[ix * (natoms * 3) + jx];
                molhess.row(ix)[jy] += molhess_threads.data(thread_id)[ix * (natoms * 3) + jy];
                molhess.row(ix)[jz] += molhess_threads.data(thread_id)[ix * (natoms * 3) + jz];

                molhess.row(iy)[jx] += molhess_threads.data(thread_id)[iy * (natoms * 3) + jx];
                molhess.row(iy)[jy] += molhess_threads.data(thread_id)[iy * (natoms * 3) + jy];
                molhess.row(iy)[jz] += molhess_threads.data(thread_id)[iy * (natoms * 3) + jz];

                molhess.row(iz)[jx] += molhess_threads.data(thread_id)[iz * (natoms * 3) + jx];
                molhess.row(iz)[jy] += molhess_threads.data(thread_id)[iz * (natoms * 3) + jy];
                molhess.row(iz)[jz] += molhess_threads.data(thread_id)[iz * (natoms * 3) + jz];
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

CDenseMatrix
CXCMolecularHessian::_integrateFxcHessianForLDA(const CMolecule&        molecule,
                                                const CMolecularBasis&  basis,
                                                const CAODensityMatrix& gsDensityMatrix,
                                                const CMolecularGrid&   molecularGrid,
                                                const CXCNewFunctional& xcFunctional) const
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // AO-to-atom mapping

    std::vector<int32_t> ao_to_atom_ids(naos);

    _computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.getNumberOfAtoms();

    CMemBlock2D<double> molhess_threads(natoms * 3 * natoms * 3, nthreads);

    // memory blocks for GTOs on grid points

    CMemBlock2D<double> gaos(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaox(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);
    CMemBlock2D<double> gaoy(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);
    CMemBlock2D<double> gaoz(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos);  // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2rho2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();
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

        auto boxdim = gtoeval::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // pre-screening of GTOs

        timer.start("GTO pre-screening");

        gtoeval::preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 1, _screeningThresholdForGTOValues, boxdim);  // 1st order GTO derivative

        timer.stop("GTO pre-screening");

        // GTO values on grid points

        timer.start("OMP GTO evaluation");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForGGA(gaos,
                                             gaox,
                                             gaoy,
                                             gaoz,
                                             gtovec,
                                             xcoords,
                                             ycoords,
                                             zcoords,
                                             gridblockpos,
                                             grid_batch_offset,
                                             grid_batch_size,
                                             skip_cgto_ids);
        }

        timer.stop("OMP GTO evaluation");

        timer.start("GTO screening");

        int32_t aocount = 0;

        for (int32_t nu = 0; nu < naos; nu++)
        {
            if (skip_ao_ids.data()[nu]) continue;

            bool skip = true;

            auto gaos_nu = gaos.data(nu);

            auto gaox_nu = gaox.data(nu);
            auto gaoy_nu = gaoy.data(nu);
            auto gaoz_nu = gaoz.data(nu);

            for (int32_t g = 0; g < npoints; g++)
            {
                if ((std::fabs(gaos_nu[g]) > _screeningThresholdForGTOValues) || (std::fabs(gaox_nu[g]) > _screeningThresholdForGTOValues) ||
                    (std::fabs(gaoy_nu[g]) > _screeningThresholdForGTOValues) || (std::fabs(gaoz_nu[g]) > _screeningThresholdForGTOValues))
                {
                    skip = false;

                    break;
                }
            }

            if (!skip)
            {
                aoinds[aocount] = nu;

                ++aocount;
            }
        }

        CDenseMatrix mat_chi(aocount, npoints);

        CDenseMatrix mat_chi_x(aocount, npoints);
        CDenseMatrix mat_chi_y(aocount, npoints);
        CDenseMatrix mat_chi_z(aocount, npoints);

        for (int32_t i = 0; i < aocount; i++)
        {
            std::memcpy(mat_chi.row(i), gaos.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_x.row(i), gaox.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_y.row(i), gaoy.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_z.row(i), gaoz.data(aoinds[i]), npoints * sizeof(double));
        }

        timer.stop("GTO screening");

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForLDA(rho, npoints, mat_chi, gs_sub_dens_mat, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, npoints);
        CDenseMatrix dengrady(natoms, npoints);
        CDenseMatrix dengradz(natoms, npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(gs_sub_dens_mat, mat_chi);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto naos = mat_chi.getNumberOfRows();

        auto F_val = mat_F.values();

        auto chi_x_val = mat_chi_x.values();
        auto chi_y_val = mat_chi_y.values();
        auto chi_z_val = mat_chi_z.values();

        auto gdenx = dengradx.values();
        auto gdeny = dengrady.values();
        auto gdenz = dengradz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * npoints;

                auto nu_offset = nu * npoints;

                #pragma omp simd aligned(gdenx, gdeny, gdenz, F_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];
                    gdeny[atom_g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];
                    gdenz[atom_g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];
                }
            }
        }

        timer.stop("Density grad. grid rho");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_fxc_for_lda(npoints, rho, v2rho2);

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        timer.stop("XC functional eval.");

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molhess_threads.data(thread_id);

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

                    #pragma omp simd reduction(+ : gatmxx, gatmxy, gatmxz, gatmyx, gatmyy, gatmyz, gatmzx, gatmzy, gatmzz) aligned(local_weights, \
                            v2rho2, gdenx, gdeny, gdenz : VLX_ALIGN)
                    for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                    {
                        auto ig = i_offset + g;
                        auto jg = j_offset + g;

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

        timer.stop("Accumulate gradient");
    }

    // destroy GTOs container

    delete gtovec;

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
                molhess.row(ix)[jx] += molhess_threads.data(thread_id)[ix * (natoms * 3) + jx];
                molhess.row(ix)[jy] += molhess_threads.data(thread_id)[ix * (natoms * 3) + jy];
                molhess.row(ix)[jz] += molhess_threads.data(thread_id)[ix * (natoms * 3) + jz];

                molhess.row(iy)[jx] += molhess_threads.data(thread_id)[iy * (natoms * 3) + jx];
                molhess.row(iy)[jy] += molhess_threads.data(thread_id)[iy * (natoms * 3) + jy];
                molhess.row(iy)[jz] += molhess_threads.data(thread_id)[iy * (natoms * 3) + jz];

                molhess.row(iz)[jx] += molhess_threads.data(thread_id)[iz * (natoms * 3) + jx];
                molhess.row(iz)[jy] += molhess_threads.data(thread_id)[iz * (natoms * 3) + jy];
                molhess.row(iz)[jz] += molhess_threads.data(thread_id)[iz * (natoms * 3) + jz];
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

void
CXCMolecularHessian::_computeAOtoAtomMapping(std::vector<int32_t>& ao_to_atom_ids, const CMolecule& molecule, const CMolecularBasis& basis) const
{
    for (int32_t iatom = 0; iatom < molecule.getNumberOfAtoms(); iatom++)
    {
        CGtoContainer* atomgtovec = new CGtoContainer(molecule, basis, iatom, 1);

        for (int32_t i = 0; i < atomgtovec->getNumberOfGtoBlocks(); i++)
        {
            auto bgtos = atomgtovec->getGtoBlock(i);

            auto bang = bgtos.getAngularMomentum();

            auto bnspher = angmom::to_SphericalComponents(bang);

            for (int32_t j = 0; j < bgtos.getNumberOfContrGtos(); j++)
            {
                for (int32_t k = 0; k < bnspher; k++)
                {
                    auto idx = bgtos.getIdentifiers(k)[j];

                    ao_to_atom_ids[idx] = iatom;
                }
            }
        }

        delete atomgtovec;
    }
}
