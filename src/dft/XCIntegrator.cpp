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

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "DenseLinearAlgebra.hpp"
#include "DenseMatrix.hpp"
#include "GtoFunc.hpp"
#include "GtoValues.hpp"
#include "MathFunc.hpp"
#include "Prescreener.hpp"

CXCIntegrator::CXCIntegrator()

    : _screeningThresholdForGTOValues(1.0e-12)
{
}

CDenseMatrix
CXCIntegrator::integrateVxcFock(const CMolecule&       molecule,
                                const CMolecularBasis& basis,
                                const CDenseMatrix&    densityMatrix,
                                const CMolecularGrid&  molecularGrid) const
{
    return _integrateVxcFockForLDA(molecule, basis, densityMatrix, molecularGrid);
}

CDenseMatrix
CXCIntegrator::_integrateVxcFockForLDA(const CMolecule&       molecule,
                                       const CMolecularBasis& basis,
                                       const CDenseMatrix&    densityMatrix,
                                       const CMolecularGrid&  molecularGrid) const
{
    // number of OpenMP threads

    auto nthreads = omp_get_max_threads();

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    int64_t naos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();

        const auto ang = gto_block.getAngularMomentum();

        naos += ncgtos * (ang * 2 + 1);
    }

    // Kohn-Sham matrix

    CDenseMatrix mat_Vxc(densityMatrix.getNumberOfRows(), densityMatrix.getNumberOfColumns());

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    CDenseMatrix gaos(naos, max_npoints_per_box);

    /*
    // mapping between AO indices before and after screening

    std::vector<int64_t> aoinds(naos);

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

    for (int64_t box_id = 0; box_id < counts.size(); box_id++)
    {
        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        std::vector<int64_t> post_ao_inds;

        for (const auto& gto_block : gto_blocks)
        {
            auto gto_orb_inds = gto_block.getOrbitalIndexes();

            auto gto_ang = gto_block.getAngularMomentum();

            // prescreen GTO block

            auto [cgto_mask, ao_mask] = prescr::preScreenGtoBlock(gto_block, 0, _screeningThresholdForGTOValues, boxdim);  // 0th order GTO derivative

            auto pre_aocount = mathfunc::countSignificantElements(ao_mask);

            std::vector<int64_t> pre_ao_inds;

            for (int64_t comp = 0, aocount = 0; comp < gto_ang * 2 + 1; comp++)
            {
                for (int64_t ind = 1; ind < static_cast<int64_t>(gto_orb_inds.size()); ind++, aocount++)
                {
                    if (ao_mask[aocount] == 1) pre_ao_inds.push_back(comp * gto_orb_inds[0] + gto_orb_inds[ind]);
                }
            }

            // GTO values on grid points

#pragma omp parallel
            {
                auto thread_id = omp_get_thread_num();

                auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

                auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

                auto cmat = gtoval::getGtoValuesForLda(gto_block,
                                                       grid_batch_size,
                                                       xcoords + gridblockpos + grid_batch_offset,
                                                       ycoords + gridblockpos + grid_batch_offset,
                                                       zcoords + gridblockpos + grid_batch_offset,
                                                       cgto_mask);

                auto submat_ptr = cmat.getSubMatrix({0, 0});

                auto gaos_ptr = submat_ptr->getData();

                for (int64_t nu = 0; nu < pre_aocount; nu++)
                {
                    std::memcpy(gaos.row(pre_ao_inds[nu]) + grid_batch_offset, gaos_ptr + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            // post-screening

            for (int64_t nu = 0; nu < pre_aocount; nu++)
            {
                bool skip = true;

                auto gaos_nu = gaos.row(pre_ao_inds[nu]);

                for (int64_t g = 0; g < npoints; g++)
                {
                    if (std::fabs(gaos_nu[g]) > _screeningThresholdForGTOValues)
                    {
                        skip = false;
                        break;
                    }
                }

                if (!skip) post_ao_inds.push_back(pre_ao_inds[nu]);
            }

            auto aocount = static_cast<int64_t>(post_ao_inds.size());

            if (box_id == 1)
            {
                std::cout << "GTO block:" << std::endl;

                std::cout << "  angular momentum: " << gto_ang << std::endl;

                std::cout << "  n_basis_funcs: " << gto_block.getNumberOfBasisFunctions() << std::endl;

                std::cout << "  orb_inds: ";
                for (const auto& x : gto_orb_inds)
                {
                    std::cout << x << ", ";
                }
                std::cout << std::endl;

                std::cout << "  cgto_mask: ";
                for (const auto x : cgto_mask)
                {
                    std::cout << x << ", ";
                }
                std::cout << std::endl;

                std::cout << "  ao_mask: ";
                for (const auto x : ao_mask)
                {
                    std::cout << x << ", ";
                }
                std::cout << std::endl;

                std::cout << "pre_ao_inds: ";
                for (const auto x : pre_ao_inds)
                {
                    std::cout << x << ", ";
                }
                std::cout << std::endl;

                std::cout << "====" << std::endl;
            }
        }

        std::sort(post_ao_inds.begin(), post_ao_inds.end());

        auto aocount = static_cast<int64_t>(post_ao_inds.size());

        if (box_id == 1)
        {
            std::cout << "npoints: " << npoints << std::endl;
            std::cout << "gridblockpos: " << gridblockpos << std::endl;

            std::cout << "aocount: " << aocount << std::endl;

            std::cout << "post_ao_inds: ";
            for (const auto x : post_ao_inds)
            {
                std::cout << x << ", ";
            }
            std::cout << std::endl;

            std::cout << "gaos: " << std::endl;
            for (const auto i : post_ao_inds)
            {
                auto gaos_i = gaos.row(i);

                std::cout << "  " << std::endl;
                for (int64_t g = 0; g < npoints; g++)
                {
                    std::cout << gaos_i[g] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        /*
        CDenseMatrix mat_chi(aocount, npoints);

        for (int64_t i = 0; i < aocount; i++)
        {
            std::memcpy(mat_chi.row(i), gaos.data(aoinds[i]), npoints * sizeof(double));
        }

        timer.stop("GTO screening");

        if (aocount == 0) continue;

        // generate sub density matrix and density grid

        if (closedshell)
        {
            timer.start("Density matrix slicing");

            auto sub_dens_mat = submat::getSubDensityMatrix(densityMatrix, 0, std::string("ALPHA"), aoinds, aocount, naos);

            timer.stop("Density matrix slicing");

            dengridgen::generateDensityForLDA(rho, npoints, mat_chi, sub_dens_mat, timer);
        }
        else
        {
            timer.start("Density matrix slicing");

            auto sub_dens_mat_a = submat::getSubDensityMatrix(densityMatrix, 0, std::string("ALPHA"), aoinds, aocount, naos);
            auto sub_dens_mat_b = submat::getSubDensityMatrix(densityMatrix, 0, std::string("BETA"), aoinds, aocount, naos);

            timer.stop("Density matrix slicing");

            dengridgen::generateDensityForLDA(rho, npoints, mat_chi, sub_dens_mat_a, sub_dens_mat_b, timer);
        }

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_lda(npoints, rho, exc, vrho);

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        timer.stop("XC functional eval.");

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

    // mat_Vxc.setNumberOfElectrons(nele);

    // mat_Vxc.setExchangeCorrelationEnergy(xcene);

    return mat_Vxc;
}

auto
CXCIntegrator::computeGtoValuesOnGridPoints(const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid) const
    -> CDenseMatrix
{
    auto nthreads = omp_get_max_threads();

    // GTO blocks and number of AOs

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

    for (int32_t box_id = 0; box_id < counts.size(); box_id++)
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

                auto cmat = gtoval::getGtoValuesForLda(gto_block,
                                                       grid_batch_size,
                                                       xcoords + gridblockpos + grid_batch_offset,
                                                       ycoords + gridblockpos + grid_batch_offset,
                                                       zcoords + gridblockpos + grid_batch_offset,
                                                       cgto_mask);

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
