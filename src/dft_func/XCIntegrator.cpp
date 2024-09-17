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

#include "XCIntegrator.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
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
#include "MpiFunc.hpp"
#include "MultiTimer.hpp"
#include "Prescreener.hpp"
#include "StringFormat.hpp"
#include "XCIntegratorForGGA.hpp"
#include "XCIntegratorForLDA.hpp"
#include "XCIntegratorForMGGA.hpp"

CXCIntegrator::CXCIntegrator(MPI_Comm comm)

    : _screeningThresholdForGTOValues(1.0e-12)
{
    _locComm = comm;
}

auto
CXCIntegrator::integrateVxcFock(const CMolecule&                  molecule,
                                const CMolecularBasis&            basis,
                                const std::vector<const double*>& gsDensityPointers,
                                const CMolecularGrid&             molecularGrid,
                                const std::string&                xcFuncLabel) const -> CAOKohnShamMatrix
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    return integrateVxcFock(molecule, basis, gsDensityPointers, molecularGrid, fvxc);
}

auto
CXCIntegrator::integrateVxcFock(const CMolecule&                  molecule,
                                const CMolecularBasis&            basis,
                                const std::vector<const double*>& gsDensityPointers,
                                const CMolecularGrid&             molecularGrid,
                                const CXCFunctional&              fvxc) const -> CAOKohnShamMatrix
{
    auto xcfuntype = fvxc.getFunctionalType();

    auto flag = (gsDensityPointers.size() == 1) ? std::string("CLOSEDSHELL") : std::string("OPENSHELL");

    std::string erropenshell("XCIntegrator.integrateVxcFock: Only implemented for closed-shell");

    if (xcfuntype == xcfun::lda)
    {
        return xcintlda::integrateVxcFockForLDA(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, flag);
    }

    if (xcfuntype == xcfun::gga)
    {
        return xcintgga::integrateVxcFockForGGA(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, flag);
    }

    if (xcfuntype == xcfun::mgga)
    {
        return xcintmgga::integrateVxcFockForMGGA(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, flag);
    }

    std::string errxcfuntype("XCIntegrator.integrateVxcFock: Only implemented for LDA/GGA/meta-GGA");

    errors::assertMsgCritical(false, errxcfuntype);

    return CAOKohnShamMatrix();
}

auto
CXCIntegrator::integrateFxcFock(const std::vector<double*>&       aoFockPointers,
                                const CMolecule&                  molecule,
                                const CMolecularBasis&            basis,
                                const std::vector<const double*>& rwDensityPointers,
                                const std::vector<const double*>& gsDensityPointers,
                                const CMolecularGrid&             molecularGrid,
                                const std::string&                xcFuncLabel) const -> void
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    return integrateFxcFock(aoFockPointers, molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, fvxc);
}

auto
CXCIntegrator::integrateFxcFock(const std::vector<double*>&       aoFockPointers,
                                const CMolecule&                  molecule,
                                const CMolecularBasis&            basis,
                                const std::vector<const double*>& rwDensityPointers,
                                const std::vector<const double*>& gsDensityPointers,
                                const CMolecularGrid&             molecularGrid,
                                const CXCFunctional&              fvxc) const -> void
{
    auto xcfuntype = fvxc.getFunctionalType();

    if (gsDensityPointers.size() == 1)
    {
        if (xcfuntype == xcfun::lda)
        {
            xcintlda::integrateFxcFockForLDA(
                aoFockPointers, molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
            xcintgga::integrateFxcFockForGGA(
                aoFockPointers, molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else if (xcfuntype == xcfun::mgga)
        {
            xcintmgga::integrateFxcFockForMGGA(
                aoFockPointers, molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else
        {
            std::string errxcfuntype("XCIntegrator.integrateFxcFock: Only implemented for LDA/GGA/meta-GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCIntegrator.integrateFxcFock: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }
}

auto
CXCIntegrator::computeGtoValuesOnGridPoints(const CMolecule&       molecule,
                                            const CMolecularBasis& basis,
                                            const CMolecularGrid&  molecularGrid) const -> CDenseMatrix
{
    // number of OpenMP threads

    auto nthreads = omp_get_max_threads();

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

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

    for (size_t box_id = 0; box_id < counts.size(); box_id++)
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
                // prescreen GTO block

                // 0th order GTO derivative
                auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(gto_block, 0, _screeningThresholdForGTOValues, boxdim);

                // GTO values on grid points

                auto cmat = gtoval::get_gto_values_for_lda(gto_block, grid_x, grid_y, grid_z, cgto_mask);

                auto submat_ptr = cmat.sub_matrix({0, 0});

                auto subgaos_ptr = submat_ptr->data();

                for (int nu = 0; nu < static_cast<int>(pre_ao_inds.size()); nu++)
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
