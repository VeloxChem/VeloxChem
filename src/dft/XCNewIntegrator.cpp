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

#include "XCNewIntegrator.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <omp.h>

#include "AngularMomentum.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DensityGridGenerator.hpp"
#include "DensityGridType.hpp"
#include "FunctionalParser.hpp"
#include "GtoEvaluator.hpp"
#include "GtoFunc.hpp"
#include "XCFuncType.hpp"
#include "XCVarsType.hpp"

CXCNewIntegrator::CXCNewIntegrator(MPI_Comm comm)

    : _screeningThresholdForGTOValues(1.0e-12)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;
}

CXCNewIntegrator::~CXCNewIntegrator()
{
}

CAOKohnShamMatrix
CXCNewIntegrator::integrateVxcFock(const CMolecule&        molecule,
                                   const CMolecularBasis&  basis,
                                   const CAODensityMatrix& densityMatrix,
                                   const CMolecularGrid&   molecularGrid,
                                   const std::string&      xcFuncLabel) const
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (densityMatrix.isClosedShell())
    {
        if (xcfuntype == xcfun::lda)
        {
           return _integrateVxcFockForLDA(molecule, basis, densityMatrix, molecularGrid, fvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
           return _integrateVxcFockForGGA(molecule, basis, densityMatrix, molecularGrid, fvxc);
        }
        else
        {
            std::string errxcfuntype("XCNewIntegrator.integrateVxcFock: Only implemented for LDA/GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCNewIntegrator.integrateVxcFock: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }

    return CAOKohnShamMatrix();
}

CAOKohnShamMatrix
CXCNewIntegrator::_integrateVxcFockForLDA(const CMolecule&        molecule,
                                          const CMolecularBasis&  basis,
                                          const CAODensityMatrix& densityMatrix,
                                          const CMolecularGrid&   molecularGrid,
                                          const CXCFunctional&    xcFunctional) const
{
    CMultiTimer timer;

    // create GTOs container

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    // number of AOs

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    CAOKohnShamMatrix mat_Vxc(densityMatrix.getNumberOfRows(0), densityMatrix.getNumberOfColumns(0), true);

    mat_Vxc.zero();

    double nele = 0.0, xcene = 0.0;

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    timer.start("Total timing");

    for (int32_t box_id = 0; box_id < counts.size(); box_id++)
    {
        auto count = counts.data()[box_id];

        auto displ = displacements.data()[box_id];

        timer.start("GTO evaluation");

        // grid points in box

        auto npoints = count;

        auto xcoords = molecularGrid.getCoordinatesX() + displ;

        auto ycoords = molecularGrid.getCoordinatesY() + displ;

        auto zcoords = molecularGrid.getCoordinatesZ() + displ;

        auto weights = molecularGrid.getWeights() + displ;

        // determine spatial extent of grid points

        double xmin = xcoords[0], ymin = ycoords[0], zmin = zcoords[0];

        double xmax = xcoords[0], ymax = ycoords[0], zmax = zcoords[0];

        for (int32_t g = 0; g < npoints; g++)
        {
            if (xmin > xcoords[g]) xmin = xcoords[g];

            if (ymin > ycoords[g]) ymin = ycoords[g];

            if (zmin > zcoords[g]) zmin = zcoords[g];

            if (xmax < xcoords[g]) xmax = xcoords[g];

            if (ymax < ycoords[g]) ymax = ycoords[g];

            if (zmax < zcoords[g]) zmax = zcoords[g];
        }

        std::array<double, 6> boxdim({xmin, ymin, zmin, xmax, ymax, zmax});

        // GTO values on grid points

        CMemBlock2D<double> gaos(npoints, naos);

        #pragma omp parallel
        {
            auto nthreads = omp_get_max_threads();

            auto thread_id = omp_get_thread_num();

            auto batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            CMemBlock2D<double> local_gaos(batch_size, naos);

            local_gaos.zero();

            gtorec::computeGtosValuesForLDA(local_gaos, gtovec, xcoords, ycoords, zcoords, batch_offset, batch_size,
                                            boxdim, _screeningThresholdForGTOValues);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                std::memcpy(gaos.data(nu) + batch_offset, local_gaos.data(nu), batch_size * sizeof(double));
            }
        }

        timer.stop("GTO evaluation");

        timer.start("GTO screening");

        std::vector<int32_t> aoinds(naos);

        int32_t aocount = 0;

        for (int32_t nu = 0; nu < naos; nu++)
        {
            bool skip = true;

            auto gaos_nu = gaos.data(nu);

            for (int32_t g = 0; g < npoints; g++)
            {
                if (std::fabs(gaos_nu[g]) > _screeningThresholdForGTOValues)
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

        aoinds.resize(aocount);

        CDenseMatrix mat_chi(aocount, npoints);

        for (int32_t i = 0; i < aocount; i++)
        {
            std::memcpy(mat_chi.row(i), gaos.data(aoinds[i]), npoints * sizeof(double));
        }

        timer.stop("GTO screening");

        timer.start("Density matrix");

        // generate sub density matrix

        CAODensityMatrix sub_dens_mat;

        if (aocount < naos)
        {
            CDenseMatrix sub_dens(aocount, aocount);

            const CDenseMatrix& dens = densityMatrix.getReferenceToDensity(0);

            for (int32_t i = 0; i < aocount; i++)
            {
                auto sub_dens_row = sub_dens.row(i);

                auto dens_row = dens.row(aoinds[i]);

                for (int32_t j = 0; j < aocount; j++)
                {
                    sub_dens_row[j] = dens_row[aoinds[j]];
                }
            }

            sub_dens_mat = CAODensityMatrix({sub_dens}, denmat::rest);
        }
        else
        {
            sub_dens_mat = densityMatrix;
        }

        timer.stop("Density matrix");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto dengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, sub_dens_mat, xcfuntype, timer);

        timer.start("XC functional");

        // compute exchange-correlation functional derivative

        CXCGradientGrid vxcgrid(npoints, dengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxcgrid, dengrid);

        timer.stop("XC functional");

        // compute partial contribution to Vxc matrix

        auto partial_mat_Vxc = _integratePartialVxcFockForLDA(npoints, xcoords, ycoords, zcoords, weights, mat_chi, vxcgrid, timer);

        timer.start("Vxc matrix dist.");

        if (aocount < naos)
        {
            for (int32_t row = 0; row < partial_mat_Vxc.getNumberOfRows(); row++)
            {
                auto row_orig = aoinds[row];

                auto mat_Vxc_row_orig = mat_Vxc.getMatrix(0) + row_orig * naos;

                auto partial_mat_Vxc_row = partial_mat_Vxc.row(row);

                for (int32_t col = 0; col < partial_mat_Vxc.getNumberOfColumns(); col++)
                {
                    auto col_orig = aoinds[col];

                    mat_Vxc_row_orig[col_orig] += partial_mat_Vxc_row[col];
                }
            }
        }
        else
        {
            mat_Vxc.addMatrixContribution(partial_mat_Vxc);
        }

        timer.stop("Vxc matrix dist.");

        timer.start("XC energy");

        // compute partial contribution to XC energy

        auto rhoa = dengrid.alphaDensity(0);

        auto efunc = vxcgrid.xcFunctionalValues();

        for (int32_t i = 0; i < npoints; i++)
        {
            nele += weights[i] * rhoa[i];

            xcene += weights[i] * efunc[i];
        }

        timer.stop("XC energy");
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    std::cout << "Timing of new integrator" << std::endl;

    std::cout << "------------------------" << std::endl;

    std::cout << timer.getSummary() << std::endl;

    mat_Vxc.setNumberOfElectrons(nele);

    mat_Vxc.setExchangeCorrelationEnergy(xcene);

    return mat_Vxc;
}

CAOKohnShamMatrix
CXCNewIntegrator::_integrateVxcFockForGGA(const CMolecule&        molecule,
                                          const CMolecularBasis&  basis,
                                          const CAODensityMatrix& densityMatrix,
                                          const CMolecularGrid&   molecularGrid,
                                          const CXCFunctional&    xcFunctional) const
{
    CMultiTimer timer;

    // create GTOs container

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis, "ATOMGTOS");

    // number of AOs

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    CAOKohnShamMatrix mat_Vxc(densityMatrix.getNumberOfRows(0), densityMatrix.getNumberOfColumns(0), true);

    mat_Vxc.zero();

    double nele = 0.0, xcene = 0.0;

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    timer.start("Total timing");

    int32_t skip_count_total = 0, num_blocks_total = 0;

    for (int32_t box_id = 0; box_id < counts.size(); box_id++)
    {
        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        auto xcoords = molecularGrid.getCoordinatesX();

        auto ycoords = molecularGrid.getCoordinatesY();

        auto zcoords = molecularGrid.getCoordinatesZ();

        auto weights = molecularGrid.getWeights();

        timer.start("AO pre-selection");

        // determine spatial extent of grid points

        double xmin = xcoords[gridblockpos], ymin = ycoords[gridblockpos], zmin = zcoords[gridblockpos];

        double xmax = xcoords[gridblockpos], ymax = ycoords[gridblockpos], zmax = zcoords[gridblockpos];

        for (int32_t g = 0; g < npoints; g++)
        {
            xmin = std::min(xmin, xcoords[gridblockpos + g]);

            ymin = std::min(ymin, ycoords[gridblockpos + g]);

            zmin = std::min(zmin, zcoords[gridblockpos + g]);

            xmax = std::max(xmax, xcoords[gridblockpos + g]);

            ymax = std::max(ymax, ycoords[gridblockpos + g]);

            zmax = std::max(zmax, zcoords[gridblockpos + g]);
        }

        std::array<double, 6> boxdim({xmin, ymin, zmin, xmax, ymax, zmax});

        auto numgtoblocks = gtovec->getNumberOfGtoBlocks();

        CMemBlock<int32_t> skip_block_ids(numgtoblocks);

        for (int32_t i = 0; i < numgtoblocks; i++)
        {
            skip_block_ids.data()[i] = 0;

            auto bgtos = gtovec->getGtoBlock(i);

            auto bang = bgtos.getAngularMomentum();

            auto bfnorms = bgtos.getNormFactors();

            auto bfexps = bgtos.getExponents();

            auto bfx = bgtos.getCoordinatesX();

            auto bfy = bgtos.getCoordinatesY();

            auto bfz = bgtos.getCoordinatesZ();

            auto spos = bgtos.getStartPositions();

            auto epos = bgtos.getEndPositions();

            auto firstprim = spos[0];

            double rx = std::max({xmin - bfx[firstprim], bfx[firstprim] - xmax, 0.0});

            double ry = std::max({ymin - bfy[firstprim], bfy[firstprim] - ymax, 0.0});

            double rz = std::max({zmin - bfz[firstprim], bfz[firstprim] - zmax, 0.0});

            auto r2 = rx * rx + ry * ry + rz * rz;

            if (r2 > 1.0)
            {
                auto minexp = bfexps[firstprim];

                auto maxexp = bfexps[firstprim];

                auto maxcoef = std::fabs(bfnorms[firstprim]);

                for (int32_t j = 0; j < bgtos.getNumberOfContrGtos(); j++)
                {
                    for (int32_t iprim = spos[j]; iprim < epos[j]; iprim++)
                    {
                        auto bexp = bfexps[iprim];

                        auto bnorm = std::fabs(bfnorms[iprim]);

                        if (minexp > bexp) minexp = bexp;

                        if (maxexp < bexp) maxexp = bexp;

                        if (maxcoef < bnorm) maxcoef = bnorm;
                    }
                }

                // gto  :           r^{ang}   |C| exp(-alpha r^2)
                // gto_m:           r^{ang-1} |C| exp(-alpha r^2)
                // gto_p: (2 alpha) r^{ang+1} |C| exp(-alpha r^2)

                // Note that gto_m < gto (r > 1)

                auto r = std::sqrt(r2);

                auto gtolimit = maxcoef * std::exp(-minexp * r2);

                for (int32_t ipow = 0; ipow < bang; ipow++) gtolimit *= r;

                auto gtolimit_p = 2.0 * maxexp * r * gtolimit;

                if (std::max(gtolimit, gtolimit_p) < _screeningThresholdForGTOValues)
                {
                    skip_block_ids.data()[i] = 1;

                    ++skip_count_total;
                }
            }

            ++num_blocks_total;
        }

        timer.stop("AO pre-selection");

        timer.start("GTO evaluation");

        // GTO values on grid points

        CMemBlock2D<double> gaos(npoints, naos);

        CMemBlock2D<double> gaox(npoints, naos);

        CMemBlock2D<double> gaoy(npoints, naos);

        CMemBlock2D<double> gaoz(npoints, naos);

        #pragma omp parallel shared(npoints, naos, gtovec, xcoords, ycoords, zcoords, gridblockpos, skip_block_ids, gaos, gaox, gaoy, gaoz)
        {
            auto nthreads = omp_get_max_threads();

            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            CMemBlock2D<double> local_gaos(grid_batch_size, naos);

            CMemBlock2D<double> local_gaox(grid_batch_size, naos);

            CMemBlock2D<double> local_gaoy(grid_batch_size, naos);

            CMemBlock2D<double> local_gaoz(grid_batch_size, naos);

            local_gaos.zero();

            local_gaox.zero();

            local_gaoy.zero();

            local_gaoz.zero();

            gtoeval::computeGtosValuesForGGA(local_gaos, local_gaox, local_gaoy, local_gaoz, gtovec, xcoords, ycoords, zcoords,

                                             gridblockpos, grid_batch_offset, grid_batch_size, skip_block_ids);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                std::memcpy(gaos.data(nu) + grid_batch_offset, local_gaos.data(nu), grid_batch_size * sizeof(double));

                std::memcpy(gaox.data(nu) + grid_batch_offset, local_gaox.data(nu), grid_batch_size * sizeof(double));

                std::memcpy(gaoy.data(nu) + grid_batch_offset, local_gaoy.data(nu), grid_batch_size * sizeof(double));

                std::memcpy(gaoz.data(nu) + grid_batch_offset, local_gaoz.data(nu), grid_batch_size * sizeof(double));
            }
        }

        timer.stop("GTO evaluation");

        timer.start("GTO screening");

        std::vector<int32_t> aoinds(naos);

        int32_t aocount = 0;

        for (int32_t nu = 0; nu < naos; nu++)
        {
            bool skip = true;

            auto gaos_nu = gaos.data(nu);

            auto gaox_nu = gaox.data(nu);

            auto gaoy_nu = gaoy.data(nu);

            auto gaoz_nu = gaoz.data(nu);

            for (int32_t g = 0; g < npoints; g++)
            {
                if ((std::fabs(gaos_nu[g]) > _screeningThresholdForGTOValues) ||
                    (std::fabs(gaox_nu[g]) > _screeningThresholdForGTOValues) ||
                    (std::fabs(gaoy_nu[g]) > _screeningThresholdForGTOValues) ||
                    (std::fabs(gaoz_nu[g]) > _screeningThresholdForGTOValues))
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

        aoinds.resize(aocount);

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

        timer.start("Density matrix");

        // generate sub density matrix

        CAODensityMatrix sub_dens_mat;

        if (aocount < naos)
        {
            CDenseMatrix sub_dens(aocount, aocount);

            const CDenseMatrix& dens = densityMatrix.getReferenceToDensity(0);

            for (int32_t i = 0; i < aocount; i++)
            {
                auto sub_dens_row = sub_dens.row(i);

                auto dens_row = dens.row(aoinds[i]);

                for (int32_t j = 0; j < aocount; j++)
                {
                    sub_dens_row[j] = dens_row[aoinds[j]];
                }
            }

            sub_dens_mat = CAODensityMatrix({sub_dens}, denmat::rest);
        }
        else
        {
            sub_dens_mat = densityMatrix;
        }

        timer.stop("Density matrix");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto dengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat, xcfuntype, timer);

        timer.start("XC functional");

        // compute exchange-correlation functional derivative

        CXCGradientGrid vxcgrid(npoints, dengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxcgrid, dengrid);

        timer.stop("XC functional");

        // compute partial contribution to Vxc matrix

        auto partial_mat_Vxc = _integratePartialVxcFockForGGA(gridblockpos, npoints, weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                              vxcgrid, dengrid, timer);

        timer.start("Vxc matrix dist.");

        if (aocount < naos)
        {
            for (int32_t row = 0; row < partial_mat_Vxc.getNumberOfRows(); row++)
            {
                auto row_orig = aoinds[row];

                auto mat_Vxc_row_orig = mat_Vxc.getMatrix(0) + row_orig * naos;

                auto partial_mat_Vxc_row = partial_mat_Vxc.row(row);

                for (int32_t col = 0; col < partial_mat_Vxc.getNumberOfColumns(); col++)
                {
                    auto col_orig = aoinds[col];

                    mat_Vxc_row_orig[col_orig] += partial_mat_Vxc_row[col];
                }
            }
        }
        else
        {
            mat_Vxc.addMatrixContribution(partial_mat_Vxc);
        }

        timer.stop("Vxc matrix dist.");

        timer.start("XC energy");

        // compute partial contribution to XC energy

        auto rhoa = dengrid.alphaDensity(0);

        auto efunc = vxcgrid.xcFunctionalValues();

        #pragma omp parallel for simd reduction(+ : nele, xcene)
        for (int32_t g = 0; g < npoints; g++)
        {
            nele += weights[gridblockpos + g] * rhoa[g];

            xcene += weights[gridblockpos + g] * efunc[g];
        }

        timer.stop("XC energy");
    }

    std::cout << "Total skip count: " << skip_count_total << " / " << num_blocks_total << std::endl;

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    std::cout << "Timing of new integrator" << std::endl;

    std::cout << "------------------------" << std::endl;

    std::cout << timer.getSummary() << std::endl;

    mat_Vxc.setNumberOfElectrons(nele);

    mat_Vxc.setExchangeCorrelationEnergy(xcene);

    return mat_Vxc;
}

CDenseMatrix
CXCNewIntegrator::_integratePartialVxcFockForLDA(const int32_t          npoints,
                                                 const double*          xcoords,
                                                 const double*          ycoords,
                                                 const double*          zcoords,
                                                 const double*          weights,
                                                 const CDenseMatrix&    gtoValues,
                                                 const CXCGradientGrid& xcGradientGrid,
                                                 CMultiTimer&           timer) const
{
    // GTO values on grid points

    const CDenseMatrix& mat_chi = gtoValues;

    // exchange-correlation functional derivative

    auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix G");

    auto naos = mat_chi.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);

    mat_G.zero();

    for (int32_t nu = 0; nu < naos; nu++)
    {
        auto G_nu = mat_G.row(nu);

        auto chi_nu = mat_chi.row(nu);

        for (int32_t g = 0; g < npoints; g++)
        {
            G_nu[g] = weights[g] * grhoa[g] * chi_nu[g];
        }
    }

    timer.stop("Vxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix matmul");

    auto mat_Vxc = denblas::multABt(mat_chi, mat_G);

    timer.stop("Vxc matrix matmul");

    return mat_Vxc;
}

CDenseMatrix
CXCNewIntegrator::_integratePartialVxcFockForGGA(const int32_t          gridblockpos,
                                                 const int32_t          npoints,
                                                 const double*          weights,
                                                 const CDenseMatrix&    gtoValues,
                                                 const CDenseMatrix&    gtoValuesX,
                                                 const CDenseMatrix&    gtoValuesY,
                                                 const CDenseMatrix&    gtoValuesZ,
                                                 const CXCGradientGrid& xcGradientGrid,
                                                 const CDensityGrid&    densityGrid,
                                                 CMultiTimer&           timer) const
{
    // GTO values on grid points

    const CDenseMatrix& mat_chi = gtoValues;

    const CDenseMatrix& mat_chi_x = gtoValuesX;

    const CDenseMatrix& mat_chi_y = gtoValuesY;

    const CDenseMatrix& mat_chi_z = gtoValuesZ;

    // exchange-correlation functional derivative

    auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);

    auto ggrada = xcGradientGrid.xcGradientValues(xcvars::grada);

    auto ggradab = xcGradientGrid.xcGradientValues(xcvars::gradab);

    // set up pointers to density gradient norms

    auto ngrada = densityGrid.alphaDensityGradient(0);

    auto gradax = densityGrid.alphaDensityGradientX(0);

    auto graday = densityGrid.alphaDensityGradientY(0);

    auto gradaz = densityGrid.alphaDensityGradientZ(0);

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix G");

    auto naos = mat_chi.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);

    CDenseMatrix mat_G_gga(naos, npoints);

    mat_G.zero();

    mat_G_gga.zero();

    for (int32_t nu = 0; nu < naos; nu++)
    {
        auto G_nu = mat_G.row(nu);

        auto G_gga_nu = mat_G_gga.row(nu);

        auto chi_nu = mat_chi.row(nu);

        auto chi_x_nu = mat_chi_x.row(nu);

        auto chi_y_nu = mat_chi_y.row(nu);

        auto chi_z_nu = mat_chi_z.row(nu);

        for (int32_t g = 0; g < npoints; g++)
        {
            G_nu[g] = weights[gridblockpos + g] * grhoa[g] * chi_nu[g];

            G_gga_nu[g] = weights[gridblockpos + g] * (ggrada[g] / ngrada[g] + ggradab[g]) *

                          (gradax[g] * chi_x_nu[g] + graday[g] * chi_y_nu[g] + gradaz[g] * chi_z_nu[g]);
        }
    }

    timer.stop("Vxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix matmul");

    auto mat_Vxc = denblas::multABt(mat_chi, denblas::addAB(mat_G, mat_G_gga, 2.0));

    mat_Vxc.symmetrizeAndScale(0.5);

    timer.stop("Vxc matrix matmul");

    return mat_Vxc;
}
