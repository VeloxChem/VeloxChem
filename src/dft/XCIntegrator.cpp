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
#include "GtoFunc.hpp"
#include "GtoValues.hpp"
#include "MathFunc.hpp"
#include "Prescreener.hpp"
#include "StringFormat.hpp"

CXCIntegrator::CXCIntegrator(MPI_Comm comm)

    : _screeningThresholdForGTOValues(1.0e-12)
{
    _locComm = comm;
}

auto
CXCIntegrator::integrateVxcFock(const CMolecule&       molecule,
                                const CMolecularBasis& basis,
                                const CDenseMatrix&    densityMatrix,
                                const CMolecularGrid&  molecularGrid,
                                const std::string&     flag) const -> CAOKohnShamMatrix
{
    return _integrateVxcFockForLDA(molecule, basis, densityMatrix, molecularGrid, flag);
}

auto
CXCIntegrator::_integrateVxcFockForLDA(const CMolecule&       molecule,
                                       const CMolecularBasis& basis,
                                       const CDenseMatrix&    densityMatrix,
                                       const CMolecularGrid&  molecularGrid,
                                       const std::string&     flag) const -> CAOKohnShamMatrix
{
    // number of OpenMP threads

    auto nthreads = omp_get_max_threads();

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    std::string errnaos("XCIntegrator._integrateVxcFockForLDA: Inconsistent number of AOs");

    errors::assertMsgCritical((naos == densityMatrix.getNumberOfRows()) && (naos == densityMatrix.getNumberOfColumns()), errnaos);

    // Kohn-Sham matrix

    bool closedshell = (fstr::upcase(flag) == std::string("CLOSEDSHELL"));

    CAOKohnShamMatrix mat_Vxc(naos, naos, closedshell);

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    CDenseMatrix gaos(naos, max_npoints_per_box);

    /*
    // density and functional derivatives

    auto       ldafunc = xcFunctional.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    CMemBlock<double> local_weights_data(max_npoints_per_box);

    CMemBlock<double> rho_data(dim->rho * max_npoints_per_box);

    CMemBlock<double> exc_data(dim->zk * max_npoints_per_box);
    CMemBlock<double> vrho_data(dim->vrho * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();

    auto exc  = exc_data.data();
    auto vrho = vrho_data.data();
    */

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(2 * max_npoints_per_box);

    std::vector<double> exc_data(1 * max_npoints_per_box);
    std::vector<double> vrho_data(2 * max_npoints_per_box);

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

    xc_func_type _func;

    auto funcID = xc_functional_get_number(std::string("LDA_X").c_str());

    auto xc_err = xc_func_init(&_func, funcID, XC_POLARIZED);

    for (int64_t box_id = 0; box_id < counts.size(); box_id++)
    {
        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // prescreening

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
        }

        // post-screening

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

        if (aocount == 0) continue;

        // generate sub density matrix and density grid

        auto sub_dens_mat = dftsubmat::getSubDensityMatrix(densityMatrix, aoinds);

        dengridgen::generateDensityForLDA(rho, mat_chi, sub_dens_mat);

        // compute exchange-correlation functional derivative

        xc_lda_exc_vxc(&_func, npoints, rho, exc, vrho);

        for (int64_t g = 0; g < npoints; g++)
        {
            local_weights[g] = weights[gridblockpos + g];
        }

        auto partial_mat_Vxc = _integratePartialVxcFockForLDA(local_weights, mat_chi, vrho);

        dftsubmat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc, aoinds);

        for (int64_t g = 0; g < npoints; g++)
        {
            auto rho_total = rho[2 * g + 0] + rho[2 * g + 1];

            nele += local_weights[g] * rho_total;

            xcene += local_weights[g] * exc[g] * rho_total;
        }

        /*
        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_lda(npoints, rho, exc, vrho);

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        timer.stop("XC functional eval.");
        */

        /*
        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        if (closedshell)
        {
            auto partial_mat_Vxc = _integratePartialVxcFockForLDA(npoints, local_weights, mat_chi, vrho, timer);

            timer.start("Vxc matrix dist.");

            submat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc, aoinds, aocount, naos);

            timer.stop("Vxc matrix dist.");
        }
        else
        {
            auto partial_mat_Vxc_ab = _integratePartialVxcFockForLDAOpenShell(npoints, local_weights, mat_chi, vrho, timer);

            timer.start("Vxc matrix dist.");

            submat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc_ab, aoinds, aocount, naos);

            timer.stop("Vxc matrix dist.");
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
        */
    }

    xc_func_end(&_func);

    mat_Vxc.setNumberOfElectrons(nele);

    mat_Vxc.setExchangeCorrelationEnergy(xcene);

    return mat_Vxc;
}

auto
CXCIntegrator::_integratePartialVxcFockForLDA(const double* weights, const CDenseMatrix& gtoValues, const double* vrho) const -> CDenseMatrix
{
    const auto npoints = gtoValues.getNumberOfColumns();

    // GTO values on grid points

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

    auto mat_Vxc = denblas::multABt(gtoValues, mat_G);

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
