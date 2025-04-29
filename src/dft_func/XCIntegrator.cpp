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

#include "XCIntegrator.hpp"

#include <omp.h>

#include <cmath>
#include <cstring>
#include <vector>

#include "DenseMatrix.hpp"
#include "ErrorHandler.hpp"
#include "FunctionalParser.hpp"
#include "GtoFunc.hpp"
#include "GtoValues.hpp"
#include "MathFunc.hpp"
#include "Prescreener.hpp"
#include "XCIntegratorForGGA.hpp"
#include "XCIntegratorForLDA.hpp"
#include "XCIntegratorForMGGA.hpp"
#include "XCIntegratorForPDFT.hpp"

CXCIntegrator::CXCIntegrator()

    : _screeningThresholdForGTOValues(1.0e-12)
{
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
        if (flag == std::string("CLOSEDSHELL"))
        {
            return xcintlda::integrateVxcFockForLdaClosedShell(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else
        {
            return xcintlda::integrateVxcFockForLdaOpenShell(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
    }

    if (xcfuntype == xcfun::gga)
    {
        if (flag == std::string("CLOSEDSHELL"))
        {
            return xcintgga::integrateVxcFockForGgaClosedShell(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else
        {
            return xcintgga::integrateVxcFockForGgaOpenShell(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
    }

    if (xcfuntype == xcfun::mgga)
    {
        if (flag == std::string("CLOSEDSHELL"))
        {
            return xcintmgga::integrateVxcFockForMetaGgaClosedShell(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else
        {
            return xcintmgga::integrateVxcFockForMetaGgaOpenShell(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
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
            xcintlda::integrateFxcFockForLdaClosedShell(
                aoFockPointers, molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
            xcintgga::integrateFxcFockForGgaClosedShell(
                aoFockPointers, molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else if (xcfuntype == xcfun::mgga)
        {
            xcintmgga::integrateFxcFockForMetaGgaClosedShell(
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
CXCIntegrator::integrateKxcFock(const std::vector<double*>& aoFockPointers,
                                const CMolecule&        molecule,
                                const CMolecularBasis&  basis,
                                const std::vector<const double*>& rwDensityPointers,
                                const std::vector<const double*>& rw2DensityPointers,
                                const std::vector<const double*>& gsDensityPointers,
                                const CMolecularGrid&   molecularGrid,
                                const std::string&      xcFuncLabel,
                                const std::string&      quadMode) const -> void
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    return integrateKxcFock(aoFockPointers, molecule, basis, rwDensityPointers, rw2DensityPointers, gsDensityPointers, molecularGrid, fvxc, quadMode);
}

auto
CXCIntegrator::integrateKxcFock(const std::vector<double*>& aoFockPointers,
                                const CMolecule&        molecule,
                                const CMolecularBasis&  basis,
                                const std::vector<const double*>& rwDensityPointers,
                                const std::vector<const double*>& rw2DensityPointers,
                                const std::vector<const double*>& gsDensityPointers,
                                const CMolecularGrid&   molecularGrid,
                                const CXCFunctional&    fvxc,
                                const std::string&      quadMode) const -> void
{
    auto xcfuntype = fvxc.getFunctionalType();

    if (gsDensityPointers.size() == 1)
    {
        if (xcfuntype == xcfun::lda)
        {
            xcintlda::integrateKxcFockForLdaClosedShell(aoFockPointers, molecule, basis, rwDensityPointers, rw2DensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, quadMode);
        }
        else if (xcfuntype == xcfun::gga)
        {
            xcintgga::integrateKxcFockForGgaClosedShell(aoFockPointers, molecule, basis, rwDensityPointers, rw2DensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, quadMode);
        }
        else if (xcfuntype == xcfun::mgga)
        {
            xcintmgga::integrateKxcFockForMetaGgaClosedShell(aoFockPointers, molecule, basis, rwDensityPointers, rw2DensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, quadMode);
        }
        else
        {
            std::string errxcfuntype("XCIntegrator.integrateKxcFock: Only implemented for LDA/GGA/meta-GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCIntegrator.integrateKxcFock: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }
}

auto
CXCIntegrator::integrateKxcLxcFock(const std::vector<double*>& aoFockPointers,
                                   const CMolecule&        molecule,
                                   const CMolecularBasis&  basis,
                                   const std::vector<const double*>& rwDensityPointers,
                                   const std::vector<const double*>& rw2DensityPointers,
                                   const std::vector<const double*>& rw3DensityPointers,
                                   const std::vector<const double*>& gsDensityPointers,
                                   const CMolecularGrid&   molecularGrid,
                                   const std::string&      xcFuncLabel,
                                   const std::string&      cubeMode) const -> void
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    return integrateKxcLxcFock(aoFockPointers, molecule, basis, rwDensityPointers, rw2DensityPointers, rw3DensityPointers, gsDensityPointers, molecularGrid, fvxc, cubeMode);
}

auto
CXCIntegrator::integrateKxcLxcFock(const std::vector<double*>& aoFockPointers,
                                   const CMolecule&        molecule,
                                   const CMolecularBasis&  basis,
                                   const std::vector<const double*>& rwDensityPointers,
                                   const std::vector<const double*>& rw2DensityPointers,
                                   const std::vector<const double*>& rw3DensityPointers,
                                   const std::vector<const double*>& gsDensityPointers,
                                   const CMolecularGrid&   molecularGrid,
                                   const CXCFunctional&    fvxc,
                                   const std::string&      cubeMode) const -> void
{
    auto xcfuntype = fvxc.getFunctionalType();

    if (gsDensityPointers.size() == 1)
    {
        if (xcfuntype == xcfun::lda)
        {
            xcintlda::integrateKxcLxcFockForLdaClosedShell(
                aoFockPointers, molecule, basis, rwDensityPointers, rw2DensityPointers, rw3DensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, cubeMode);
        }
        else if (xcfuntype == xcfun::gga)
        {
            xcintgga::integrateKxcLxcFockForGgaClosedShell(
                aoFockPointers, molecule, basis, rwDensityPointers, rw2DensityPointers, rw3DensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, cubeMode);
        }
        else if (xcfuntype == xcfun::mgga)
        {
            xcintmgga::integrateKxcLxcFockForMetaGgaClosedShell(
                aoFockPointers, molecule, basis, rwDensityPointers, rw2DensityPointers, rw3DensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, cubeMode);
        }
        else
        {
            std::string errxcfuntype("XCIntegrator.integrateKxcLxcFock: Only implemented for LDA/GGA/meta-GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCIntegrator.integrateKxcLxcFock: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }
}

auto
CXCIntegrator::integrateVxcPDFT(CAOKohnShamMatrix&                  aoFockMatrix,
                                CDense4DTensor&                     tensorWxc,
                                const CMolecule&                    molecule,
                                const CMolecularBasis&              basis,
                                const double*                       densityMatrixPointer,
                                const CDenseMatrix&                 twoBodyDensityMatrix,
                                const CDenseMatrix&                 activeMOs,
                                const CMolecularGrid&               molecularGrid,
                                const CXCPairDensityFunctional&     fvxc,
                                const double                        rs_omega) const -> void
{
    auto xcfuntype = fvxc.getFunctionalType();

    if (xcfuntype == "PLDA")
    {
        xcintpdft::integrateVxcPDFTForLDA(aoFockMatrix,
                                          tensorWxc,
                                          molecule,
                                          basis,
                                          densityMatrixPointer,
                                          twoBodyDensityMatrix,
                                          activeMOs,
                                          molecularGrid,
                                          _screeningThresholdForGTOValues,
                                          fvxc, rs_omega);
    }
    else if (xcfuntype == "PGGA")
    {
        xcintpdft::integrateVxcPDFTForGGA(aoFockMatrix,
                                          tensorWxc,
                                          molecule,
                                          basis,
                                          densityMatrixPointer,
                                          twoBodyDensityMatrix,
                                          activeMOs,
                                          molecularGrid,
                                          _screeningThresholdForGTOValues,
                                          fvxc, rs_omega);
    }
    else
    {
        std::string errxcfuntype("XCIntegrator.integrateVxcPDFT: Only implemented for PLDA/PGGA");

        errors::assertMsgCritical(false, errxcfuntype);
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

                if (cmat.is_empty()) continue;

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

auto
CXCIntegrator::computeGtoValuesAndDerivativesOnGridPoints(const CMolecule&       molecule,
                                                          const CMolecularBasis& basis,
                                                          const CMolecularGrid&  molecularGrid) const -> std::vector<CDenseMatrix>
{
    // number of OpenMP threads

    auto nthreads = omp_get_max_threads();

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // GTO values on grid points

    CDenseMatrix allgtovalues(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_x(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_y(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_z(naos, molecularGrid.getNumberOfGridPoints());

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

                // 1st order GTO derivative
                auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(gto_block, 1, _screeningThresholdForGTOValues, boxdim);

                // GTO values on grid points

                auto cmat = gtoval::get_gto_values_for_gga(gto_block, grid_x, grid_y, grid_z, cgto_mask);

                if (cmat.is_empty()) continue;

                auto submat_0_ptr = cmat.sub_matrix({0, 0});
                auto submat_x_ptr = cmat.sub_matrix({1, 0});
                auto submat_y_ptr = cmat.sub_matrix({1, 1});
                auto submat_z_ptr = cmat.sub_matrix({1, 2});

                auto subgaos_0_ptr = submat_0_ptr->data();
                auto subgaos_x_ptr = submat_x_ptr->data();
                auto subgaos_y_ptr = submat_y_ptr->data();
                auto subgaos_z_ptr = submat_z_ptr->data();

                for (int nu = 0; nu < static_cast<int>(pre_ao_inds.size()); nu++)
                {
                    std::memcpy(allgtovalues.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_0_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));
                    std::memcpy(allgtovalues_x.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_x_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));
                    std::memcpy(allgtovalues_y.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_y_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));
                    std::memcpy(allgtovalues_z.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_z_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));
                }
            }
        }
    }

    return std::vector<CDenseMatrix>({allgtovalues, allgtovalues_x, allgtovalues_y, allgtovalues_z});
}

auto
CXCIntegrator::computeGtoValuesAndSecondOrderDerivativesOnGridPoints(const CMolecule&       molecule,
                                                                     const CMolecularBasis& basis,
                                                                     const CMolecularGrid&  molecularGrid) const -> std::vector<CDenseMatrix>
{
    // number of OpenMP threads

    auto nthreads = omp_get_max_threads();

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // GTO values on grid points

    CDenseMatrix allgtovalues(naos, molecularGrid.getNumberOfGridPoints());

    CDenseMatrix allgtovalues_x(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_y(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_z(naos, molecularGrid.getNumberOfGridPoints());

    CDenseMatrix allgtovalues_xx(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_xy(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_xz(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_yy(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_yz(naos, molecularGrid.getNumberOfGridPoints());
    CDenseMatrix allgtovalues_zz(naos, molecularGrid.getNumberOfGridPoints());

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

                // 1st order GTO derivative
                auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(gto_block, 1, _screeningThresholdForGTOValues, boxdim);

                // GTO values on grid points

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

                auto subgaos_0_ptr = submat_0_ptr->data();

                auto subgaos_x_ptr = submat_x_ptr->data();
                auto subgaos_y_ptr = submat_y_ptr->data();
                auto subgaos_z_ptr = submat_z_ptr->data();

                auto subgaos_xx_ptr = submat_xx_ptr->data();
                auto subgaos_xy_ptr = submat_xy_ptr->data();
                auto subgaos_xz_ptr = submat_xz_ptr->data();
                auto subgaos_yy_ptr = submat_yy_ptr->data();
                auto subgaos_yz_ptr = submat_yz_ptr->data();
                auto subgaos_zz_ptr = submat_zz_ptr->data();

                for (int nu = 0; nu < static_cast<int>(pre_ao_inds.size()); nu++)
                {
                    std::memcpy(allgtovalues.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_0_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));

                    std::memcpy(allgtovalues_x.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_x_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));
                    std::memcpy(allgtovalues_y.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_y_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));
                    std::memcpy(allgtovalues_z.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_z_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));

                    std::memcpy(allgtovalues_xx.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_xx_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));
                    std::memcpy(allgtovalues_xy.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_xy_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));
                    std::memcpy(allgtovalues_xz.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_xz_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));
                    std::memcpy(allgtovalues_yy.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_yy_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));
                    std::memcpy(allgtovalues_yz.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_yz_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));
                    std::memcpy(allgtovalues_zz.row(pre_ao_inds[nu]) + gridblockpos + grid_batch_offset,
                                subgaos_zz_ptr + nu * grid_batch_size,
                                grid_batch_size * sizeof(double));
                }
            }
        }
    }

    return std::vector<CDenseMatrix>({allgtovalues, allgtovalues_x, allgtovalues_y, allgtovalues_z,
                                      allgtovalues_xx, allgtovalues_xy, allgtovalues_xz,
                                      allgtovalues_yy, allgtovalues_yz, allgtovalues_zz});
}
