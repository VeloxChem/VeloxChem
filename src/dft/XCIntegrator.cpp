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

#include "XCIntegrator.hpp"

#include <omp.h>
#include <xc.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "DenseLinearAlgebra.hpp"
#include "DenseMatrix.hpp"
#include "DensityGridGenerator.hpp"
#include "DftSubMatrix.hpp"
#include "ErrorHandler.hpp"
#include "FunctionalParser.hpp"
#include "GtoFunc.hpp"
#include "GtoValues.hpp"
#include "MathFunc.hpp"
#include "MultiTimer.hpp"
#include "Prescreener.hpp"
#include "StringFormat.hpp"

CXCIntegrator::CXCIntegrator(MPI_Comm comm)

    : _screeningThresholdForGTOValues(1.0e-12)
{
    _locComm = comm;
}

auto
CXCIntegrator::integrateVxcFock(const CMolecule&        molecule,
                                const CMolecularBasis&  basis,
                                const CAODensityMatrix& densityMatrix,
                                const CMolecularGrid&   molecularGrid,
                                const std::string&      xcFuncLabel) const -> CAOKohnShamMatrix
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    auto flag = densityMatrix.isClosedShell() ? std::string("CLOSEDSHELL") : std::string("OPENSHELL");

    if (xcfuntype == xcfun::lda) return _integrateVxcFockForLDA(molecule, basis, densityMatrix, molecularGrid, fvxc, flag);

    if (xcfuntype == xcfun::gga) return _integrateVxcFockForGGA(molecule, basis, densityMatrix, molecularGrid, fvxc, flag);

    std::string errxcfuntype("XCIntegrator.integrateVxcFock: Only implemented for LDA");

    errors::assertMsgCritical(false, errxcfuntype);

    return CAOKohnShamMatrix();
}

auto
CXCIntegrator::_integrateVxcFockForLDA(const CMolecule&        molecule,
                                       const CMolecularBasis&  basis,
                                       const CAODensityMatrix& densityMatrix,
                                       const CMolecularGrid&   molecularGrid,
                                       const CXCFunctional&    xcFunctional,
                                       const std::string&      flag) const -> CAOKohnShamMatrix
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    std::string errnaos("XCIntegrator._integrateVxcFockForLDA: Inconsistent number of AOs");

    errors::assertMsgCritical((naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)), errnaos);

    // Kohn-Sham matrix

    bool closedshell = (fstr::upcase(flag) == std::string("CLOSEDSHELL"));

    CAOKohnShamMatrix mat_Vxc(naos, naos, closedshell);

    mat_Vxc.zero();

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    CDenseMatrix gaos(naos, max_npoints_per_box);

    // density and functional derivatives

    auto       ldafunc = xcFunctional.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);

    std::vector<double> exc_data(dim->zk * max_npoints_per_box);
    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();

    auto exc  = exc_data.data();
    auto vrho = vrho_data.data();

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

    for (int64_t box_id = 0; box_id < counts.size(); box_id++)
    {
        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // prescreening

        timer.start("GTO pre-screening");

        std::vector<std::vector<int64_t>> cgto_masks, ao_masks, pre_ao_inds_blocks;

        for (const auto& gto_block : gto_blocks)
        {
            // 0th order GTO derivative
            auto [cgto_mask, ao_mask] = prescr::preScreenGtoBlock(gto_block, 0, _screeningThresholdForGTOValues, boxdim);

            cgto_masks.push_back(cgto_mask);

            ao_masks.push_back(ao_mask);

            std::vector<int64_t> pre_ao_inds;

            auto gto_ao_inds = gto_block.getAtomicOrbitalsIndexes();

            for (int64_t i = 0; i < static_cast<int64_t>(gto_ao_inds.size()); i++)
            {
                if (ao_mask[i] == 1) pre_ao_inds.push_back(gto_ao_inds[i]);
            }

            pre_ao_inds_blocks.push_back(pre_ao_inds);
        }

        timer.stop("GTO pre-screening");

        // GTO values on grid points

        timer.start("OMP GTO evaluation");

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

            for (size_t i_block = 0; i_block < gto_blocks.size(); i_block++)
            {
                const auto& gto_block = gto_blocks[i_block];

                const auto& cgto_mask = cgto_masks[i_block];

                const auto& ao_mask = ao_masks[i_block];

                const auto& pre_ao_inds = pre_ao_inds_blocks[i_block];

                auto cmat = gtoval::getGtoValuesForLda(gto_block, grid_x, grid_y, grid_z, cgto_mask);

                auto submat_ptr = cmat.getSubMatrix({0, 0});

                auto submat_data = submat_ptr->getData();

                for (int64_t nu = 0; nu < static_cast<int64_t>(pre_ao_inds.size()); nu++)
                {
                    std::memcpy(gaos.row(pre_ao_inds[nu]) + grid_batch_offset, submat_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");
        }

        timer.stop("OMP GTO evaluation");

        // post-screening

        timer.start("GTO screening");

        std::vector<int64_t> aoinds;

        for (const auto& pre_ao_inds : pre_ao_inds_blocks)
        {
            for (const auto nu : pre_ao_inds)
            {
                bool skip = true;

                auto gaos_nu = gaos.row(nu);

                for (int64_t g = 0; g < npoints; g++)
                {
                    if (std::fabs(gaos_nu[g]) > _screeningThresholdForGTOValues)
                    {
                        skip = false;
                        break;
                    }
                }

                if (!skip) aoinds.push_back(nu);
            }
        }

        std::sort(aoinds.begin(), aoinds.end());

        const auto aocount = static_cast<int64_t>(aoinds.size());

        CDenseMatrix mat_chi(aocount, npoints);

        for (int64_t i = 0; i < aocount; i++)
        {
            std::memcpy(mat_chi.row(i), gaos.row(aoinds[i]), npoints * sizeof(double));
        }

        timer.stop("GTO screening");

        if (aocount == 0) continue;

        // generate sub density matrix and density grid

        if (closedshell)
        {
            timer.start("Density matrix slicing");

            auto sub_dens_mat = dftsubmat::getSubDensityMatrix(densityMatrix, 0, std::string("ALPHA"), aoinds);

            timer.stop("Density matrix slicing");

            dengridgen::generateDensityForLDA(rho, mat_chi, sub_dens_mat, timer);
        }
        else
        {
            // TODO: openshell
        }

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_lda(npoints, rho, exc, vrho);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        if (closedshell)
        {
            auto partial_mat_Vxc = _integratePartialVxcFockForLDA(local_weights, mat_chi, vrho, timer);

            timer.start("Vxc matrix dist.");

            dftsubmat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc, aoinds);

            timer.stop("Vxc matrix dist.");
        }
        else
        {
            // TODO: openshell
        }

        // compute partial contribution to XC energy

        timer.start("XC energy");
        for (int64_t g = 0; g < npoints; g++)
        {
            auto rho_total = rho[2 * g + 0] + rho[2 * g + 1];

            nele += local_weights[g] * rho_total;

            xcene += local_weights[g] * exc[g] * rho_total;
        }

        timer.stop("XC energy");
    }

    mat_Vxc.setNumberOfElectrons(nele);

    mat_Vxc.setExchangeCorrelationEnergy(xcene);

    timer.stop("Total timing");

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;
    // std::cout << "OpenMP timing" << std::endl;
    // for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
    // {
    //     std::cout << "Thread " << thread_id << std::endl;
    //     std::cout << omptimers[thread_id].getSummary() << std::endl;
    // }

    return mat_Vxc;
}

auto
CXCIntegrator::_integrateVxcFockForGGA(const CMolecule&        molecule,
                                       const CMolecularBasis&  basis,
                                       const CAODensityMatrix& densityMatrix,
                                       const CMolecularGrid&   molecularGrid,
                                       const CXCFunctional&    xcFunctional,
                                       const std::string&      flag) const -> CAOKohnShamMatrix
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    std::string errnaos("XCIntegrator._integrateVxcFockForGGA: Inconsistent number of AOs");

    errors::assertMsgCritical((naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)), errnaos);

    // Kohn-Sham matrix

    bool closedshell = (fstr::upcase(flag) == std::string("CLOSEDSHELL"));

    CAOKohnShamMatrix mat_Vxc(naos, naos, closedshell);

    mat_Vxc.zero();

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    CDenseMatrix gaos(naos, max_npoints_per_box);

    CDenseMatrix gaox(naos, max_npoints_per_box);
    CDenseMatrix gaoy(naos, max_npoints_per_box);
    CDenseMatrix gaoz(naos, max_npoints_per_box);

    // density and functional derivatives

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhograd_data(dim->rho * 3 * max_npoints_per_box);
    std::vector<double> sigma_data(dim->sigma * max_npoints_per_box);

    std::vector<double> exc_data(dim->zk * max_npoints_per_box);
    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> vsigma_data(dim->vsigma * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();

    auto exc    = exc_data.data();
    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();

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

    for (int32_t box_id = 0; box_id < counts.size(); box_id++)
    {
        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // prescreening

        timer.start("GTO pre-screening");

        std::vector<std::vector<int64_t>> cgto_masks, ao_masks, pre_ao_inds_blocks;

        for (const auto& gto_block : gto_blocks)
        {
            // 1st order GTO derivative
            auto [cgto_mask, ao_mask] = prescr::preScreenGtoBlock(gto_block, 1, _screeningThresholdForGTOValues, boxdim);

            cgto_masks.push_back(cgto_mask);

            ao_masks.push_back(ao_mask);

            std::vector<int64_t> pre_ao_inds;

            auto gto_ao_inds = gto_block.getAtomicOrbitalsIndexes();

            for (int64_t i = 0; i < static_cast<int64_t>(gto_ao_inds.size()); i++)
            {
                if (ao_mask[i] == 1) pre_ao_inds.push_back(gto_ao_inds[i]);
            }

            pre_ao_inds_blocks.push_back(pre_ao_inds);
        }

        timer.stop("GTO pre-screening");

        // GTO values on grid points

        timer.start("OMP GTO evaluation");

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

            for (size_t i_block = 0; i_block < gto_blocks.size(); i_block++)
            {
                const auto& gto_block = gto_blocks[i_block];

                const auto& cgto_mask = cgto_masks[i_block];

                const auto& ao_mask = ao_masks[i_block];

                const auto& pre_ao_inds = pre_ao_inds_blocks[i_block];

                auto cmat = gtoval::getGtoValuesForGga(gto_block, grid_x, grid_y, grid_z, cgto_mask);

                auto submat_0_ptr = cmat.getSubMatrix({0, 0});
                auto submat_x_ptr = cmat.getSubMatrix({1, 0});
                auto submat_y_ptr = cmat.getSubMatrix({1, 1});
                auto submat_z_ptr = cmat.getSubMatrix({1, 2});

                auto submat_0_data = submat_0_ptr->getData();
                auto submat_x_data = submat_x_ptr->getData();
                auto submat_y_data = submat_y_ptr->getData();
                auto submat_z_data = submat_z_ptr->getData();

                for (int64_t nu = 0; nu < static_cast<int64_t>(pre_ao_inds.size()); nu++)
                {
                    std::memcpy(
                        gaos.row(pre_ao_inds[nu]) + grid_batch_offset, submat_0_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(
                        gaox.row(pre_ao_inds[nu]) + grid_batch_offset, submat_x_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(
                        gaoy.row(pre_ao_inds[nu]) + grid_batch_offset, submat_y_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(
                        gaoz.row(pre_ao_inds[nu]) + grid_batch_offset, submat_z_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");
        }

        timer.stop("OMP GTO evaluation");

        // post-screening

        timer.start("GTO screening");

        std::vector<int64_t> aoinds;

        for (const auto& pre_ao_inds : pre_ao_inds_blocks)
        {
            for (const auto nu : pre_ao_inds)
            {
                bool skip = true;

                auto gaos_nu = gaos.row(nu);
                auto gaox_nu = gaox.row(nu);
                auto gaoy_nu = gaoy.row(nu);
                auto gaoz_nu = gaoz.row(nu);

                for (int64_t g = 0; g < npoints; g++)
                {
                    if ((std::fabs(gaos_nu[g]) > _screeningThresholdForGTOValues) || (std::fabs(gaox_nu[g]) > _screeningThresholdForGTOValues) ||
                        (std::fabs(gaoy_nu[g]) > _screeningThresholdForGTOValues) || (std::fabs(gaoz_nu[g]) > _screeningThresholdForGTOValues))
                    {
                        skip = false;
                        break;
                    }
                }

                if (!skip) aoinds.push_back(nu);
            }
        }

        std::sort(aoinds.begin(), aoinds.end());

        const auto aocount = static_cast<int64_t>(aoinds.size());

        CDenseMatrix mat_chi(aocount, npoints);
        CDenseMatrix mat_chi_x(aocount, npoints);
        CDenseMatrix mat_chi_y(aocount, npoints);
        CDenseMatrix mat_chi_z(aocount, npoints);

        for (int64_t i = 0; i < aocount; i++)
        {
            std::memcpy(mat_chi.row(i), gaos.row(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_x.row(i), gaox.row(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_y.row(i), gaoy.row(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_z.row(i), gaoz.row(aoinds[i]), npoints * sizeof(double));
        }

        timer.stop("GTO screening");

        if (aocount == 0) continue;

        // generate sub density matrix and density grid

        if (closedshell)
        {
            timer.start("Density matrix slicing");

            auto sub_dens_mat = dftsubmat::getSubDensityMatrix(densityMatrix, 0, std::string("ALPHA"), aoinds);

            timer.stop("Density matrix slicing");

            dengridgen::generateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat, timer);
        }
        else
        {
            // TODO: openshell
        }

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_gga(npoints, rho, sigma, exc, vrho, vsigma);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        if (closedshell)
        {
            auto partial_mat_Vxc =
                _integratePartialVxcFockForGGA(local_weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rhograd, vrho, vsigma, timer);

            timer.start("Vxc matrix dist.");

            dftsubmat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc, aoinds);

            timer.stop("Vxc matrix dist.");
        }
        else
        {
            // TODO: openshell
        }

        // compute partial contribution to XC energy

        timer.start("XC energy");

        for (int32_t g = 0; g < npoints; g++)
        {
            auto rho_total = rho[2 * g + 0] + rho[2 * g + 1];

            nele += local_weights[g] * rho_total;

            xcene += local_weights[g] * exc[g] * rho_total;
        }

        timer.stop("XC energy");
    }

    mat_Vxc.setNumberOfElectrons(nele);

    mat_Vxc.setExchangeCorrelationEnergy(xcene);

    timer.stop("Total timing");

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;
    // std::cout << "OpenMP timing" << std::endl;
    // for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
    // {
    //     std::cout << "Thread " << thread_id << std::endl;
    //     std::cout << omptimers[thread_id].getSummary() << std::endl;
    // }

    return mat_Vxc;
}

auto
CXCIntegrator::_integratePartialVxcFockForLDA(const double* weights, const CDenseMatrix& gtoValues, const double* vrho, CMultiTimer timer) const
    -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    // GTO values on grid points

    timer.start("Vxc matrix G");

    auto chi_val = gtoValues.values();

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);

    auto G_val = mat_G.values();

#pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

        for (int64_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int64_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                G_val[nu_offset + g] = weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];
            }
        }
    }

    timer.stop("Vxc matrix G");

    timer.start("Vxc matrix matmul");

    auto mat_Vxc = denblas::multABt(gtoValues, mat_G);

    timer.stop("Vxc matrix matmul");

    return mat_Vxc;
}

auto
CXCIntegrator::_integratePartialVxcFockForGGA(const double*       weights,
                                              const CDenseMatrix& gtoValues,
                                              const CDenseMatrix& gtoValuesX,
                                              const CDenseMatrix& gtoValuesY,
                                              const CDenseMatrix& gtoValuesZ,
                                              const double*       rhograd,
                                              const double*       vrho,
                                              const double*       vsigma,
                                              CMultiTimer&        timer) const -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    timer.start("Vxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);
    CDenseMatrix mat_G_gga(naos, npoints);

    auto G_val     = mat_G.values();
    auto G_gga_val = mat_G_gga.values();

    auto chi_val   = gtoValues.values();
    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

#pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                auto vx = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                auto vy = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                auto vz = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                G_val[nu_offset + g] = weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];

                G_gga_val[nu_offset + g] =
                    weights[g] * (vx * chi_x_val[nu_offset + g] + vy * chi_y_val[nu_offset + g] + vz * chi_z_val[nu_offset + g]);
            }
        }
    }

    timer.stop("Vxc matrix G");

    // Note that we use matrix-matrix multiplication only once, and symmetrize
    // the result. This is because the density matrix is symmetric, and the
    // Kohn-Sham matrix from mat_G is also symmetric. Formally only the
    // mat_G_gga contribution should be symmetrized.

    timer.start("Vxc matrix matmul");

    auto mat_Vxc = denblas::multABt(gtoValues, denblas::addAB(mat_G, mat_G_gga, 2.0));

    mat_Vxc.symmetrizeAndScale(0.5);

    timer.stop("Vxc matrix matmul");

    return mat_Vxc;
}

auto
CXCIntegrator::computeGtoValuesOnGridPoints(const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid) const
    -> CDenseMatrix
{
    // number of OpenMP threads

    auto nthreads = omp_get_max_threads();

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // GTO values on grid points

    CDenseMatrix allgtovalues(naos, molecularGrid.getNumberOfGridPoints());

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    for (int64_t box_id = 0; box_id < counts.size(); box_id++)
    {
        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // compute GTO values on grid points

#pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            const auto grid_x_ptr = xcoords + gridblockpos + grid_batch_offset;
            const auto grid_y_ptr = ycoords + gridblockpos + grid_batch_offset;
            const auto grid_z_ptr = zcoords + gridblockpos + grid_batch_offset;

            std::vector<double> grid_x(grid_x_ptr, grid_x_ptr + grid_batch_size);
            std::vector<double> grid_y(grid_y_ptr, grid_y_ptr + grid_batch_size);
            std::vector<double> grid_z(grid_z_ptr, grid_z_ptr + grid_batch_size);

            // go through GTO blocks

            for (const auto& gto_block : gto_blocks)
            {
                auto gto_orb_inds = gto_block.getOrbitalIndexes();

                auto gto_ang = gto_block.getAngularMomentum();

                // prescreen GTO block

                auto [cgto_mask, ao_mask] =
                    prescr::preScreenGtoBlock(gto_block, 0, _screeningThresholdForGTOValues, boxdim);  // 0th order GTO derivative

                auto pre_aocount = mathfunc::countSignificantElements(ao_mask);

                std::vector<int64_t> pre_ao_inds;

                auto gto_ao_inds = gto_block.getAtomicOrbitalsIndexes();

                for (int64_t i = 0; i < static_cast<int64_t>(gto_ao_inds.size()); i++)
                {
                    if (ao_mask[i] == 1) pre_ao_inds.push_back(gto_ao_inds[i]);
                }

                // GTO values on grid points

                auto cmat = gtoval::getGtoValuesForLda(gto_block, grid_x, grid_y, grid_z, cgto_mask);

                auto submat_ptr = cmat.getSubMatrix({0, 0});

                auto subgaos_ptr = submat_ptr->getData();

                for (int64_t nu = 0; nu < pre_aocount; nu++)
                {
                    std::memcpy(allgtovalues.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));
                }
            }
        }
    }

    return allgtovalues;
}
