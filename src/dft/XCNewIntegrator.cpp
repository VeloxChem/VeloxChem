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
#include "AODensityMatrix.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DensityGridGenerator.hpp"
#include "DensityGridType.hpp"
#include "FunctionalParser.hpp"
#include "GridScreener.hpp"
#include "GtoEvaluator.hpp"
#include "SubMatrix.hpp"
#include "XCFuncType.hpp"
#include "XCVarsType.hpp"

CXCNewIntegrator::CXCNewIntegrator(MPI_Comm comm)

    : _screeningThresholdForGTOValues(1.0e-12)

    , _screeningThresholdForDensityValues(1.0e-13)
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

void
CXCNewIntegrator::integrateFxcFock(CAOFockMatrix&          aoFockMatrix,
                                   const CMolecule&        molecule,
                                   const CMolecularBasis&  basis,
                                   const CAODensityMatrix& rwDensityMatrix,
                                   const CAODensityMatrix& gsDensityMatrix,
                                   const CMolecularGrid&   molecularGrid,
                                   const std::string&      xcFuncLabel) const
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (rwDensityMatrix.isClosedShell() && gsDensityMatrix.isClosedShell())
    {
        if (xcfuntype == xcfun::lda)
        {
           _integrateFxcFockForLDA(aoFockMatrix, molecule, basis, rwDensityMatrix, gsDensityMatrix, molecularGrid, fvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
           _integrateFxcFockForGGA(aoFockMatrix, molecule, basis, rwDensityMatrix, gsDensityMatrix, molecularGrid, fvxc);
        }
        else
        {
            std::string errxcfuntype("XCNewIntegrator.integrateFxcFock: Only implemented for LDA/GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCNewIntegrator.integrateFxcFock: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }
}

void
CXCNewIntegrator::integrateKxcFock(CAOFockMatrix&          aoFockMatrix,
                                   const CMolecule&        molecule,
                                   const CMolecularBasis&  basis,
                                   const CAODensityMatrix& rwDensityMatrix,
                                   const CAODensityMatrix& rw2DensityMatrix,
                                   const CAODensityMatrix& gsDensityMatrix,
                                   const CMolecularGrid&   molecularGrid,
                                   const std::string&      xcFuncLabel,
                                   const std::string&      quadMode) const
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (rwDensityMatrix.isClosedShell() && rw2DensityMatrix.isClosedShell() && gsDensityMatrix.isClosedShell())
    {
        if (xcfuntype == xcfun::lda)
        {
           _integrateKxcFockForLDA(aoFockMatrix, molecule, basis, rwDensityMatrix, rw2DensityMatrix, gsDensityMatrix,

                                   molecularGrid, fvxc, quadMode);
        }
        else if (xcfuntype == xcfun::gga)
        {
           _integrateKxcFockForGGA(aoFockMatrix, molecule, basis, rwDensityMatrix, rw2DensityMatrix, gsDensityMatrix,

                                   molecularGrid, fvxc, quadMode);
        }
        else
        {
            std::string errxcfuntype("XCNewIntegrator.integrateKxcFock: Only implemented for LDA/GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCNewIntegrator.integrateKxcFock: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }
}

CAOKohnShamMatrix
CXCNewIntegrator::_integrateVxcFockForLDA(const CMolecule&        molecule,
                                          const CMolecularBasis&  basis,
                                          const CAODensityMatrix& densityMatrix,
                                          const CMolecularGrid&   molecularGrid,
                                          const CXCFunctional&    xcFunctional) const
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // Kohn-Sham matrix

    CAOKohnShamMatrix mat_Vxc(densityMatrix.getNumberOfRows(0), densityMatrix.getNumberOfColumns(0), true);

    mat_Vxc.zero();

    // memory blocks for GTOs on grid points

    CMemBlock2D<double> gaos(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos); // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

    // indices for keeping track of valid grid points

    std::vector<int32_t> screened_point_inds(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> screened_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto screened_weights = screened_weights_data.data();

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

        auto boxdim = gtoeval::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // pre-screening of GTOs

        timer.start("GTO pre-screening");

        gtoeval::preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 0, _screeningThresholdForGTOValues, boxdim);  // 0th order GTO derivative

        timer.stop("GTO pre-screening");

        // GTO values on grid points

        timer.start("OMP GTO evaluation");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            omptimers[thread_id].start("gtoeval");

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForLDA(gaos, gtovec, xcoords, ycoords, zcoords, gridblockpos,

                                             grid_batch_offset, grid_batch_size, skip_cgto_ids);

            omptimers[thread_id].stop("gtoeval");
        }

        timer.stop("OMP GTO evaluation");

        timer.start("GTO screening");

        int32_t aocount = 0;

        for (int32_t nu = 0; nu < naos; nu++)
        {
            if (skip_ao_ids.data()[nu]) continue;

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

        CDenseMatrix mat_chi(aocount, npoints);

        for (int32_t i = 0; i < aocount; i++)
        {
            std::memcpy(mat_chi.row(i), gaos.data(aoinds[i]), npoints * sizeof(double));
        }

        timer.stop("GTO screening");

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto sub_dens_mat = submat::getSubDensityMatrix(densityMatrix, 0, aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto dengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, sub_dens_mat, xcfuntype, timer);

        // screen density grid, weights and GTO matrix

        timer.start("Density screening");

        CDensityGrid screened_dengrid(dengrid);

        gridscreen::screenDensityGridForLDA(screened_point_inds, screened_dengrid, dengrid,

                                            _screeningThresholdForDensityValues);

        auto screened_npoints = screened_dengrid.getNumberOfGridPoints();

        gridscreen::screenWeights(screened_weights, gridblockpos, weights, screened_point_inds, screened_npoints);

        CDenseMatrix screened_mat_chi(mat_chi.getNumberOfRows(), screened_npoints);

        gridscreen::screenGtoMatrixForLDA(screened_mat_chi, mat_chi, screened_point_inds, screened_npoints);

        timer.stop("Density screening");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        CXCGradientGrid vxcgrid(screened_npoints, screened_dengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxcgrid, screened_dengrid);

        timer.stop("XC functional eval.");

        // compute partial contribution to Vxc matrix

        auto partial_mat_Vxc = _integratePartialVxcFockForLDA(screened_npoints, screened_weights, screened_mat_chi,

                                                              vxcgrid, timer);

        // distribute partial Vxc to full Kohn-Sham matrix

        timer.start("Vxc matrix dist.");

        submat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc, aoinds, aocount, naos);

        timer.stop("Vxc matrix dist.");

        // compute partial contribution to XC energy

        timer.start("XC energy");

        auto rhoa = screened_dengrid.alphaDensity(0);

        auto rhob = screened_dengrid.betaDensity(0);

        auto efunc = vxcgrid.xcFunctionalValues();

        for (int32_t g = 0; g < screened_npoints; g++)
        {
            nele += screened_weights[g] * (rhoa[g] + rhob[g]);

            xcene += screened_weights[g] * efunc[g];
        }

        timer.stop("XC energy");
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    //std::cout << "Timing of new integrator" << std::endl;
    //std::cout << "------------------------" << std::endl;
    //std::cout << timer.getSummary() << std::endl;
    //std::cout << "OpenMP timing" << std::endl;
    //for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
    //{
    //    std::cout << "Thread " << thread_id << std::endl;
    //    std::cout << omptimers[thread_id].getSummary() << std::endl;
    //}

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

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // Kohn-Sham matrix

    CAOKohnShamMatrix mat_Vxc(densityMatrix.getNumberOfRows(0), densityMatrix.getNumberOfColumns(0), true);

    mat_Vxc.zero();

    // memory blocks for GTOs on grid points

    CMemBlock2D<double> gaos(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaox(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaoy(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaoz(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos); // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

    // indices for keeping track of valid grid points

    std::vector<int32_t> screened_point_inds(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> screened_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto screened_weights = screened_weights_data.data();

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

            omptimers[thread_id].start("gtoeval");

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtovec, xcoords, ycoords, zcoords, gridblockpos,

                                             grid_batch_offset, grid_batch_size, skip_cgto_ids);

            omptimers[thread_id].stop("gtoeval");
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

        auto sub_dens_mat = submat::getSubDensityMatrix(densityMatrix, 0, aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto dengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                                             sub_dens_mat, xcfuntype, timer);

        // screen density grid, weights and GTO matrix

        timer.start("Density screening");

        CDensityGrid screened_dengrid(dengrid);

        gridscreen::screenDensityGridForGGA(screened_point_inds, screened_dengrid, dengrid,

                                            _screeningThresholdForDensityValues);

        auto screened_npoints = screened_dengrid.getNumberOfGridPoints();

        gridscreen::screenWeights(screened_weights, gridblockpos, weights, screened_point_inds, screened_npoints);

        CDenseMatrix screened_mat_chi(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_x(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_y(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_z(mat_chi.getNumberOfRows(), screened_npoints);

        gridscreen::screenGtoMatrixForGGA(screened_mat_chi, screened_mat_chi_x, screened_mat_chi_y, screened_mat_chi_z,

                                          mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, screened_point_inds, screened_npoints);

        timer.stop("Density screening");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        CXCGradientGrid vxcgrid(screened_npoints, screened_dengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxcgrid, screened_dengrid);

        timer.stop("XC functional eval.");

        // compute partial contribution to Vxc matrix

        auto partial_mat_Vxc = _integratePartialVxcFockForGGA(screened_npoints, screened_weights, screened_mat_chi,

                                                              screened_mat_chi_x, screened_mat_chi_y, screened_mat_chi_z,

                                                              vxcgrid, screened_dengrid, timer);

        // distribute partial Vxc to full Kohn-Sham matrix

        timer.start("Vxc matrix dist.");

        submat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc, aoinds, aocount, naos);

        timer.stop("Vxc matrix dist.");

        // compute partial contribution to XC energy

        timer.start("XC energy");

        auto rhoa = screened_dengrid.alphaDensity(0);

        auto rhob = screened_dengrid.betaDensity(0);

        auto efunc = vxcgrid.xcFunctionalValues();

        for (int32_t g = 0; g < screened_npoints; g++)
        {
            nele += screened_weights[g] * (rhoa[g] + rhob[g]);

            xcene += screened_weights[g] * efunc[g];
        }

        timer.stop("XC energy");
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    //std::cout << "Timing of new integrator" << std::endl;
    //std::cout << "------------------------" << std::endl;
    //std::cout << timer.getSummary() << std::endl;
    //std::cout << "OpenMP timing" << std::endl;
    //for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
    //{
    //    std::cout << "Thread " << thread_id << std::endl;
    //    std::cout << omptimers[thread_id].getSummary() << std::endl;
    //}

    mat_Vxc.setNumberOfElectrons(nele);

    mat_Vxc.setExchangeCorrelationEnergy(xcene);

    return mat_Vxc;
}

void
CXCNewIntegrator::_integrateFxcFockForLDA(CAOFockMatrix&          aoFockMatrix,
                                          const CMolecule&        molecule,
                                          const CMolecularBasis&  basis,
                                          const CAODensityMatrix& rwDensityMatrix,
                                          const CAODensityMatrix& gsDensityMatrix,
                                          const CMolecularGrid&   molecularGrid,
                                          const CXCFunctional&    xcFunctional) const
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // memory blocks for GTOs on grid points

    CMemBlock2D<double> gaos(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos); // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

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

        gtoeval::preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 0, _screeningThresholdForGTOValues, boxdim);  // 0th order GTO derivative

        timer.stop("GTO pre-screening");

        // GTO values on grid points

        timer.start("OMP GTO evaluation");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            omptimers[thread_id].start("gtoeval");

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForLDA(gaos, gtovec, xcoords, ycoords, zcoords, gridblockpos,

                                             grid_batch_offset, grid_batch_size, skip_cgto_ids);

            omptimers[thread_id].stop("gtoeval");
        }

        timer.stop("OMP GTO evaluation");

        timer.start("GTO screening");

        int32_t aocount = 0;

        for (int32_t nu = 0; nu < naos; nu++)
        {
            if (skip_ao_ids.data()[nu]) continue;

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

        CDenseMatrix mat_chi(aocount, npoints);

        for (int32_t i = 0; i < aocount; i++)
        {
            std::memcpy(mat_chi.row(i), gaos.data(aoinds[i]), npoints * sizeof(double));
        }

        timer.stop("GTO screening");

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto gsdengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, sub_dens_mat, xcfuntype, timer);

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        CXCHessianGrid vxc2grid(npoints, gsdengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxc2grid, gsdengrid);

        timer.stop("XC functional eval.");

        // go through rhow density matrices

        for (int32_t idensity = 0; idensity < rwDensityMatrix.getNumberOfDensityMatrices(); idensity++)
        {
            // generate sub density matrix

            timer.start("Density matrix slicing");

            auto sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, idensity, aoinds, aocount, naos);

            timer.stop("Density matrix slicing");

            // generate density grid

            auto xcfuntype = xcFunctional.getFunctionalType();

            auto rwdengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, sub_dens_mat, xcfuntype, timer);

            // compute partial contribution to Fxc matrix

            auto partial_mat_Fxc = _integratePartialFxcFockForLDA(gridblockpos, npoints, weights, mat_chi,

                                                                  vxc2grid, rwdengrid, timer);

            // distribute partial Fxc to full Fock matrix

            timer.start("Fxc matrix dist.");

            submat::distributeSubMatrixToFock(aoFockMatrix, idensity, partial_mat_Fxc, aoinds, aocount, naos);

            timer.stop("Fxc matrix dist.");
        }
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    //std::cout << "Timing of new integrator" << std::endl;
    //std::cout << "------------------------" << std::endl;
    //std::cout << timer.getSummary() << std::endl;
    //std::cout << "OpenMP timing" << std::endl;
    //for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
    //{
    //    std::cout << "Thread " << thread_id << std::endl;
    //    std::cout << omptimers[thread_id].getSummary() << std::endl;
    //}
}

void
CXCNewIntegrator::_integrateFxcFockForGGA(CAOFockMatrix&          aoFockMatrix,
                                          const CMolecule&        molecule,
                                          const CMolecularBasis&  basis,
                                          const CAODensityMatrix& rwDensityMatrix,
                                          const CAODensityMatrix& gsDensityMatrix,
                                          const CMolecularGrid&   molecularGrid,
                                          const CXCFunctional&    xcFunctional) const
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // memory blocks for GTOs on grid points

    CMemBlock2D<double> gaos(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaox(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaoy(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaoz(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos); // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

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

            omptimers[thread_id].start("gtoeval");

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtovec, xcoords, ycoords, zcoords, gridblockpos,

                                             grid_batch_offset, grid_batch_size, skip_cgto_ids);

            omptimers[thread_id].stop("gtoeval");
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

        auto sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto gsdengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                                               sub_dens_mat, xcfuntype, timer);

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        CXCGradientGrid vxcgrid(npoints, gsdengrid.getDensityGridType(), xcfuntype);

        CXCHessianGrid vxc2grid(npoints, gsdengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxcgrid, gsdengrid);

        xcFunctional.compute(vxc2grid, gsdengrid);

        timer.stop("XC functional eval.");

        // go through rhow density matrices

        for (int32_t idensity = 0; idensity < rwDensityMatrix.getNumberOfDensityMatrices(); idensity++)
        {
            // generate sub density matrix

            timer.start("Density matrix slicing");

            auto sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, idensity, aoinds, aocount, naos);

            timer.stop("Density matrix slicing");

            // generate density grid

            auto xcfuntype = xcFunctional.getFunctionalType();

            auto rwdengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                                                   sub_dens_mat, xcfuntype, timer);

            // compute partial contribution to Fxc matrix

            auto partial_mat_Fxc = _integratePartialFxcFockForGGA(gridblockpos, npoints, weights, mat_chi,

                                                                  mat_chi_x, mat_chi_y, mat_chi_z,

                                                                  vxcgrid, vxc2grid, rwdengrid, gsdengrid, timer);

            // distribute partial Fxc to full Fock matrix

            timer.start("Fxc matrix dist.");

            submat::distributeSubMatrixToFock(aoFockMatrix, idensity, partial_mat_Fxc, aoinds, aocount, naos);

            timer.stop("Fxc matrix dist.");
        }
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    //std::cout << "Timing of new integrator" << std::endl;
    //std::cout << "------------------------" << std::endl;
    //std::cout << timer.getSummary() << std::endl;
    //std::cout << "OpenMP timing" << std::endl;
    //for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
    //{
    //    std::cout << "Thread " << thread_id << std::endl;
    //    std::cout << omptimers[thread_id].getSummary() << std::endl;
    //}
}

void
CXCNewIntegrator::_integrateKxcFockForLDA(CAOFockMatrix&          aoFockMatrix,
                                          const CMolecule&        molecule,
                                          const CMolecularBasis&  basis,
                                          const CAODensityMatrix& rwDensityMatrix,
                                          const CAODensityMatrix& rw2DensityMatrix,
                                          const CAODensityMatrix& gsDensityMatrix,
                                          const CMolecularGrid&   molecularGrid,
                                          const CXCFunctional&    xcFunctional,
                                          const std::string&      quadMode) const
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // memory blocks for GTOs on grid points

    CMemBlock2D<double> gaos(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos); // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

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

        gtoeval::preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 0, _screeningThresholdForGTOValues, boxdim);  // 0th order GTO derivative

        timer.stop("GTO pre-screening");

        // GTO values on grid points

        timer.start("OMP GTO evaluation");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            omptimers[thread_id].start("gtoeval");

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForLDA(gaos, gtovec, xcoords, ycoords, zcoords, gridblockpos,

                                             grid_batch_offset, grid_batch_size, skip_cgto_ids);

            omptimers[thread_id].stop("gtoeval");
        }

        timer.stop("OMP GTO evaluation");

        timer.start("GTO screening");

        int32_t aocount = 0;

        for (int32_t nu = 0; nu < naos; nu++)
        {
            if (skip_ao_ids.data()[nu]) continue;

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

        CDenseMatrix mat_chi(aocount, npoints);

        for (int32_t i = 0; i < aocount; i++)
        {
            std::memcpy(mat_chi.row(i), gaos.data(aoinds[i]), npoints * sizeof(double));
        }

        timer.stop("GTO screening");

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, aoinds, aocount, naos);

        auto rw_sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, aoinds, aocount);

        auto rw2_sub_dens_mat = submat::getSubDensityMatrix(rw2DensityMatrix, aoinds, aocount);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto gsdengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, gs_sub_dens_mat, xcfuntype, timer);

        auto rwdengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, rw_sub_dens_mat, xcfuntype, timer);

        auto rw2dengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, rw2_sub_dens_mat, xcfuntype, timer);

        // compute perturbed density

        timer.start("Density grid quad");

        auto numdens_rw2 = rw2DensityMatrix.getNumberOfDensityMatrices();

        CDensityGridQuad rwdengridquad(npoints, numdens_rw2, xcfuntype, dengrid::ab);

        rwdengridquad.DensityProd(rwdengrid, xcfuntype, numdens_rw2, quadMode);

        timer.stop("Density grid quad");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        CXCHessianGrid vxc2grid(npoints, gsdengrid.getDensityGridType(), xcfuntype);

        CXCCubicHessianGrid vxc3grid(npoints, gsdengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxc2grid, gsdengrid);

        xcFunctional.compute(vxc3grid, gsdengrid);

        timer.stop("XC functional eval.");

        // go through density matrices

        for (int32_t idensity = 0; idensity < numdens_rw2; idensity++)
        {
            // compute partial contribution to Kxc matrix

            auto partial_mat_Kxc = _integratePartialKxcFockForLDA(gridblockpos, npoints, weights, mat_chi,

                                                                  vxc2grid, vxc3grid, rwdengridquad, rw2dengrid,

                                                                  idensity, timer);

            // distribute partial Kxc to full Fock matrix

            timer.start("Kxc matrix dist.");

            submat::distributeSubMatrixToFock(aoFockMatrix, idensity, partial_mat_Kxc, aoinds, aocount, naos);

            timer.stop("Kxc matrix dist.");
        }
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    //std::cout << "Timing of new integrator" << std::endl;
    //std::cout << "------------------------" << std::endl;
    //std::cout << timer.getSummary() << std::endl;
    //std::cout << "OpenMP timing" << std::endl;
    //for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
    //{
    //    std::cout << "Thread " << thread_id << std::endl;
    //    std::cout << omptimers[thread_id].getSummary() << std::endl;
    //}
}

void
CXCNewIntegrator::_integrateKxcFockForGGA(CAOFockMatrix&          aoFockMatrix,
                                          const CMolecule&        molecule,
                                          const CMolecularBasis&  basis,
                                          const CAODensityMatrix& rwDensityMatrix,
                                          const CAODensityMatrix& rw2DensityMatrix,
                                          const CAODensityMatrix& gsDensityMatrix,
                                          const CMolecularGrid&   molecularGrid,
                                          const CXCFunctional&    xcFunctional,
                                          const std::string&      quadMode) const
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // memory blocks for GTOs on grid points

    CMemBlock2D<double> gaos(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaox(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaoy(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaoz(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos); // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

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

            omptimers[thread_id].start("gtoeval");

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtovec, xcoords, ycoords, zcoords, gridblockpos,

                                             grid_batch_offset, grid_batch_size, skip_cgto_ids);

            omptimers[thread_id].stop("gtoeval");
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

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, aoinds, aocount, naos);

        auto rw_sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, aoinds, aocount);

        auto rw2_sub_dens_mat = submat::getSubDensityMatrix(rw2DensityMatrix, aoinds, aocount);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto gsdengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                                               gs_sub_dens_mat, xcfuntype, timer);

        auto rwdengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                                               rw_sub_dens_mat, xcfuntype, timer);

        auto rw2dengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                                                rw2_sub_dens_mat, xcfuntype, timer);

        // compute perturbed density

        timer.start("Density grid quad");

        auto numdens_rw2 = rw2DensityMatrix.getNumberOfDensityMatrices();

        CDensityGridQuad rwdengridquad(npoints, numdens_rw2, xcfuntype, dengrid::ab);

        rwdengridquad.DensityProd(rwdengrid, xcfuntype, numdens_rw2, quadMode);

        timer.stop("Density grid quad");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        CXCGradientGrid vxcgrid(npoints, gsdengrid.getDensityGridType(), xcfuntype);

        CXCHessianGrid vxc2grid(npoints, gsdengrid.getDensityGridType(), xcfuntype);

        CXCCubicHessianGrid vxc3grid(npoints, gsdengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxcgrid, gsdengrid);

        xcFunctional.compute(vxc2grid, gsdengrid);

        xcFunctional.compute(vxc3grid, gsdengrid);

        timer.stop("XC functional eval.");

        // go through density matrices

        for (int32_t idensity = 0; idensity < numdens_rw2; idensity++)
        {
            // compute partial contribution to Kxc matrix

            auto partial_mat_Kxc = _integratePartialKxcFockForGGA(gridblockpos, npoints, weights, mat_chi,

                                                                  mat_chi_x, mat_chi_y, mat_chi_z,

                                                                  vxcgrid, vxc2grid, vxc3grid,

                                                                  rwdengridquad, rw2dengrid, gsdengrid,

                                                                  idensity, timer);

            // distribute partial Kxc to full Fock matrix

            timer.start("Kxc matrix dist.");

            submat::distributeSubMatrixToFock(aoFockMatrix, idensity, partial_mat_Kxc, aoinds, aocount, naos);

            timer.stop("Kxc matrix dist.");
        }
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    //std::cout << "Timing of new integrator" << std::endl;
    //std::cout << "------------------------" << std::endl;
    //std::cout << timer.getSummary() << std::endl;
    //std::cout << "OpenMP timing" << std::endl;
    //for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
    //{
    //    std::cout << "Thread " << thread_id << std::endl;
    //    std::cout << omptimers[thread_id].getSummary() << std::endl;
    //}
}

CDenseMatrix
CXCNewIntegrator::_integratePartialVxcFockForLDA(const int32_t          npoints,
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
CXCNewIntegrator::_integratePartialVxcFockForGGA(const int32_t          npoints,
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

    // pointers to density gradient norms

    auto ngrada = densityGrid.alphaDensityGradient(0);

    auto gradax = densityGrid.alphaDensityGradientX(0);

    auto graday = densityGrid.alphaDensityGradientY(0);

    auto gradaz = densityGrid.alphaDensityGradientZ(0);

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix G");

    auto naos = mat_chi.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);

    CDenseMatrix mat_G_gga(naos, npoints);

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
            G_nu[g] = weights[g] * grhoa[g] * chi_nu[g];

            G_gga_nu[g] = weights[g] * (ggrada[g] / ngrada[g] + ggradab[g]) *

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

CDenseMatrix
CXCNewIntegrator::_integratePartialFxcFockForLDA(const int32_t         gridBlockPosition,
                                                 const int32_t         npoints,
                                                 const double*         weights,
                                                 const CDenseMatrix&   gtoValues,
                                                 const CXCHessianGrid& xcHessianGrid,
                                                 const CDensityGrid&   rwDensityGrid,
                                                 CMultiTimer&          timer) const
{
    // GTO values on grid points

    auto chi_val = gtoValues.values();

    // pointers to exchange-correlation functional derivative

    auto grho_aa = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhoa);

    auto grho_ab = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhob);

    // pointers to perturbed density

    auto rhowa = rwDensityGrid.alphaDensity(0);

    auto rhowb = rwDensityGrid.betaDensity(0);

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Fxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);

    auto G_val = mat_G.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(weights, \
                    grho_aa, grho_ab, rhowa, rhowb, \
                    G_val, chi_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                G_val[nu_offset + g] = weights[gridBlockPosition + g] *

                          (grho_aa[g] * rhowa[g] + grho_ab[g] * rhowb[g]) * chi_val[nu_offset + g];
            }
        }
    }

    timer.stop("Fxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Fxc matrix matmul");

    auto mat_Fxc = denblas::multABt(gtoValues, mat_G);

    timer.stop("Fxc matrix matmul");

    return mat_Fxc;
}

CDenseMatrix
CXCNewIntegrator::_integratePartialFxcFockForGGA(const int32_t          gridblockpos,
                                                 const int32_t          npoints,
                                                 const double*          weights,
                                                 const CDenseMatrix&    gtoValues,
                                                 const CDenseMatrix&    gtoValuesX,
                                                 const CDenseMatrix&    gtoValuesY,
                                                 const CDenseMatrix&    gtoValuesZ,
                                                 const CXCGradientGrid& xcGradientGrid,
                                                 const CXCHessianGrid&  xcHessianGrid,
                                                 const CDensityGrid&    rwDensityGrid,
                                                 const CDensityGrid&    gsDensityGrid,
                                                 CMultiTimer&           timer) const
{
    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();

    auto chi_y_val = gtoValuesY.values();

    auto chi_z_val = gtoValuesZ.values();

    // pointers to exchange-correlation functional derrivatives

    auto grho_aa = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhoa);

    auto grho_ab = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhob);

    auto gmix_aa = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::grada);

    auto gmix_ab = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::gradb);

    auto gmix_ac = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::gradab);

    auto gmix_bc = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::gradab);

    auto ggrad_aa = xcHessianGrid.xcHessianValues(xcvars::grada, xcvars::grada);

    auto ggrad_ab = xcHessianGrid.xcHessianValues(xcvars::grada, xcvars::gradb);

    auto ggrad_ac = xcHessianGrid.xcHessianValues(xcvars::grada, xcvars::gradab);

    auto ggrad_bc = xcHessianGrid.xcHessianValues(xcvars::gradb, xcvars::gradab);

    auto ggrad_cc = xcHessianGrid.xcHessianValues(xcvars::gradab, xcvars::gradab);

    auto ggrad_a = xcGradientGrid.xcGradientValues(xcvars::grada);

    auto ggrad_c = xcGradientGrid.xcGradientValues(xcvars::gradab);

    // pointers to ground state density gradient norms

    auto ngrada = gsDensityGrid.alphaDensityGradient(0);

    auto ngradb = gsDensityGrid.betaDensityGradient(0);

    auto grada_x = gsDensityGrid.alphaDensityGradientX(0);

    auto grada_y = gsDensityGrid.alphaDensityGradientY(0);

    auto grada_z = gsDensityGrid.alphaDensityGradientZ(0);

    auto gradb_x = gsDensityGrid.betaDensityGradientX(0);

    auto gradb_y = gsDensityGrid.betaDensityGradientY(0);

    auto gradb_z = gsDensityGrid.betaDensityGradientZ(0);

    // pointers to perturbed density gradient norms

    auto rhowa = rwDensityGrid.alphaDensity(0);

    auto rhowb = rwDensityGrid.betaDensity(0);

    auto gradwa_x = rwDensityGrid.alphaDensityGradientX(0);

    auto gradwa_y = rwDensityGrid.alphaDensityGradientY(0);

    auto gradwa_z = rwDensityGrid.alphaDensityGradientZ(0);

    auto gradwb_x = rwDensityGrid.betaDensityGradientX(0);

    auto gradwb_y = rwDensityGrid.betaDensityGradientY(0);

    auto gradwb_z = rwDensityGrid.betaDensityGradientZ(0);

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Fxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);

    CDenseMatrix mat_G_gga(naos, npoints);

    auto G_val = mat_G.values();

    auto G_gga_val = mat_G_gga.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(weights, ggrad_a, ggrad_c, grho_aa, grho_ab, \
                    gmix_aa, gmix_ab, gmix_ac, gmix_bc, ggrad_aa, ggrad_ab, ggrad_ac, ggrad_bc, ggrad_cc, \
                    ngrada, ngradb, grada_x, grada_y, grada_z, gradb_x, gradb_y, gradb_z, rhowa, rhowb, \
                    gradwa_x, gradwa_y, gradwa_z, gradwb_x, gradwb_y, gradwb_z, \
                    G_val, G_gga_val, chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                double w = weights[gridblockpos + g];

                double znva = 1.0 / ngrada[g];

                double znvb = 1.0 / ngradb[g];

                double rxa = znva * grada_x[g];

                double rya = znva * grada_y[g];

                double rza = znva * grada_z[g];

                double rxb = znvb * gradb_x[g];

                double ryb = znvb * gradb_y[g];

                double rzb = znvb * gradb_z[g];

                double rxwa = gradwa_x[g];

                double rywa = gradwa_y[g];

                double rzwa = gradwa_z[g];

                double rxwb = gradwb_x[g];

                double rywb = gradwb_y[g];

                double rzwb = gradwb_z[g];

                //  variations of functionals variables

                double zetaa = rxwa * rxa + rywa * rya + rzwa * rza;

                double zetab = rxwb * rxb + rywb * ryb + rzwb * rzb;

                double zetac = grada_x[g] * rxwb + grada_y[g] * rywb + grada_z[g] * rzwb +

                               gradb_x[g] * rxwa + gradb_y[g] * rywa + gradb_z[g] * rzwa;

                // first contribution

                double fac0 = gmix_aa[g] * zetaa + gmix_ab[g] * zetab + gmix_ac[g] * zetac +

                              grho_aa[g] * rhowa[g] + grho_ab[g] * rhowb[g];

                G_val[nu_offset + g] = w * fac0 * chi_val[nu_offset + g];

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                // second contribution

                double facr = gmix_aa[g] * rhowa[g] + gmix_ab[g] * rhowb[g] +

                              ggrad_aa[g] * zetaa + ggrad_ab[g] * zetab + ggrad_ac[g] * zetac;

                xcomp += facr * rxa;

                ycomp += facr * rya;

                zcomp += facr * rza;

                // third contribution

                double facz = gmix_ac[g] * rhowa[g] + gmix_bc[g] * rhowb[g] +

                              ggrad_ac[g] * zetaa + ggrad_bc[g] * zetab + ggrad_cc[g] * zetac;

                xcomp += facz * grada_x[g];

                ycomp += facz * grada_y[g];

                zcomp += facz * grada_z[g];

                // fourth contribution

                xcomp += znva * ggrad_a[g] * (rxwa - rxa * zetaa);

                ycomp += znva * ggrad_a[g] * (rywa - rya * zetaa);

                zcomp += znva * ggrad_a[g] * (rzwa - rza * zetaa);

                // fifth contribution

                xcomp += ggrad_c[g] * rxwa;

                ycomp += ggrad_c[g] * rywa;

                zcomp += ggrad_c[g] * rzwa;

                G_gga_val[nu_offset + g] = w * (xcomp * chi_x_val[nu_offset + g] +

                                                ycomp * chi_y_val[nu_offset + g] +

                                                zcomp * chi_z_val[nu_offset + g]);
            }
        }
    }

    timer.stop("Fxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Fxc matrix matmul");

    auto mat_Fxc = denblas::multABt(gtoValues, mat_G);

    auto mat_Fxc_gga = denblas::multABt(gtoValues, mat_G_gga);

    mat_Fxc_gga.symmetrize();  // matrix + matrix.T

    mat_Fxc = denblas::addAB(mat_Fxc, mat_Fxc_gga, 1.0);

    timer.stop("Fxc matrix matmul");

    return mat_Fxc;
}

CDenseMatrix
CXCNewIntegrator::_integratePartialKxcFockForLDA(const int32_t              gridBlockPosition,
                                                 const int32_t              npoints,
                                                 const double*              weights,
                                                 const CDenseMatrix&        gtoValues,
                                                 const CXCHessianGrid&      xcHessianGrid,
                                                 const CXCCubicHessianGrid& xcCubicHessianGrid,
                                                 const CDensityGridQuad&    rwDensityGridQuad,
                                                 const CDensityGrid&        rw2DensityGrid,
                                                 const int32_t              iFock,
                                                 CMultiTimer&               timer) const
{
    timer.start("Kxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    // pointers to exchange-correlation functional derrivatives

    auto grho_aa = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhoa);

    auto grho_ab = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhob);

    auto grho_aaa = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhoa);

    auto grho_aab = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhob);

    auto grho_abb = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::rhob);

    // pointers to perturbed density

    auto rhow1a = rwDensityGridQuad.rhow1rhow2(iFock);

    auto rhow12a = rw2DensityGrid.alphaDensity(iFock);

    auto rhow12b = rw2DensityGrid.betaDensity(iFock);

    timer.stop("Kxc matrix prep.");

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Kxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);

    auto G_val = mat_G.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(weights, \
                    grho_aa, grho_ab, grho_aaa, grho_aab, grho_abb, \
                    rhow1a, rhow12a, rhow12b, \
                    G_val, chi_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                G_val[nu_offset + g] = weights[gridBlockPosition + g] *

                          ((grho_aaa[g] + grho_aab[g] + grho_aab[g] + grho_abb[g]) * rhow1a[g] +

                           grho_aa[g] * rhow12a[g] + grho_ab[g] * rhow12b[g]) *

                          chi_val[nu_offset + g];
            }
        }
    }

    timer.stop("Kxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Kxc matrix matmul");

    auto mat_Kxc = denblas::multABt(gtoValues, mat_G);

    timer.stop("Kxc matrix matmul");

    return mat_Kxc;
}

CDenseMatrix
CXCNewIntegrator::_integratePartialKxcFockForGGA(const int32_t              gridblockpos,
                                                 const int32_t              npoints,
                                                 const double*              weights,
                                                 const CDenseMatrix&        gtoValues,
                                                 const CDenseMatrix&        gtoValuesX,
                                                 const CDenseMatrix&        gtoValuesY,
                                                 const CDenseMatrix&        gtoValuesZ,
                                                 const CXCGradientGrid&     xcGradientGrid,
                                                 const CXCHessianGrid&      xcHessianGrid,
                                                 const CXCCubicHessianGrid& xcCubicHessianGrid,
                                                 const CDensityGridQuad&    rwDensityGridQuad,
                                                 const CDensityGrid&        rw2DensityGrid,
                                                 const CDensityGrid&        gsDensityGrid,
                                                 const int32_t              iFock,
                                                 CMultiTimer&               timer) const
{
    timer.start("Kxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();

    auto chi_y_val = gtoValuesY.values();

    auto chi_z_val = gtoValuesZ.values();

    // pointers to exchange-correlation functional derrivatives

    auto df2000 = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhoa);

    auto df1100 = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhob);

    auto df1010 = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::grada);

    auto df1001 = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::gradb);

    auto df10001 = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::gradab);

    auto df0020 = xcHessianGrid.xcHessianValues(xcvars::grada, xcvars::grada);

    auto df0011 = xcHessianGrid.xcHessianValues(xcvars::grada, xcvars::gradb);

    auto df00101 = xcHessianGrid.xcHessianValues(xcvars::grada, xcvars::gradab);

    auto df00002 = xcHessianGrid.xcHessianValues(xcvars::gradab, xcvars::gradab);

    auto df0010 = xcGradientGrid.xcGradientValues(xcvars::grada);

    auto df00001 = xcGradientGrid.xcGradientValues(xcvars::gradab);

    auto df00011 = xcHessianGrid.xcHessianValues(xcvars::gradb, xcvars::gradab);

    auto df01001 = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::gradab);

    auto df0110 = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::grada);

    auto df3000 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhoa);

    auto df2100 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhob);

    auto df1200 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::rhob);

    auto df2010 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::grada);

    auto df0030 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::grada, xcvars::grada, xcvars::grada);

    auto df0021 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::grada, xcvars::grada, xcvars::gradb);

    auto df0012 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::grada, xcvars::gradb, xcvars::gradb);

    auto df00201 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::grada, xcvars::grada, xcvars::gradab);

    auto df00111 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::grada, xcvars::gradb, xcvars::gradab);

    auto df00102 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::grada, xcvars::gradab, xcvars::gradab);

    auto df00003 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::gradab, xcvars::gradab, xcvars::gradab);

    auto df2001 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::gradb);

    auto df1110 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::grada);

    auto df1101 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::gradb);

    auto df20001 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::gradab);

    auto df11001 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::gradab);

    auto df1020 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::grada, xcvars::grada);

    auto df1011 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::grada, xcvars::gradb);

    auto df1002 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::gradb, xcvars::gradb);

    auto df10101 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::grada, xcvars::gradab);

    auto df10002 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::gradab, xcvars::gradab);

    auto df01002 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhob, xcvars::gradab, xcvars::gradab);

    auto df0120 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhob, xcvars::grada, xcvars::grada);

    auto df0111 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhob, xcvars::grada, xcvars::gradb);

    auto df01101 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhob, xcvars::grada, xcvars::gradab);

    auto df10011 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::gradb, xcvars::gradab);

    auto df01011 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhob, xcvars::gradb, xcvars::gradab);

    auto df0210 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhob, xcvars::rhob, xcvars::grada);

    auto df02001 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhob, xcvars::rhob, xcvars::gradab);

    auto df00021 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::gradb, xcvars::gradb, xcvars::gradab);

    // pointers to ground state density gradient norms

    auto ngrada = gsDensityGrid.alphaDensityGradient(0);

    auto grada_x = gsDensityGrid.alphaDensityGradientX(0);

    auto grada_y = gsDensityGrid.alphaDensityGradientY(0);

    auto grada_z = gsDensityGrid.alphaDensityGradientZ(0);

    // pointers to perturbed density gradient norms

    auto rhow1rhow2 = rwDensityGridQuad.rhow1rhow2(iFock);

    auto rxw1rhow2 = rwDensityGridQuad.rxw1rhow2(iFock);

    auto ryw1rhow2 = rwDensityGridQuad.ryw1rhow2(iFock);

    auto rzw1rhow2 = rwDensityGridQuad.rzw1rhow2(iFock);

    auto rxw1rxw2 = rwDensityGridQuad.rxw1rxw2(iFock);

    auto rxw1ryw2 = rwDensityGridQuad.rxw1ryw2(iFock);

    auto rxw1rzw2 = rwDensityGridQuad.rxw1rzw2(iFock);

    auto ryw1rxw2 = rwDensityGridQuad.ryw1rxw2(iFock);

    auto ryw1ryw2 = rwDensityGridQuad.ryw1ryw2(iFock);

    auto ryw1rzw2 = rwDensityGridQuad.ryw1rzw2(iFock);

    auto rzw1rxw2 = rwDensityGridQuad.rzw1rxw2(iFock);

    auto rzw1ryw2 = rwDensityGridQuad.rzw1ryw2(iFock);

    auto rzw1rzw2 = rwDensityGridQuad.rzw1rzw2(iFock);

    auto rhow12a = rw2DensityGrid.alphaDensity(iFock);

    auto gradw12a_x = rw2DensityGrid.alphaDensityGradientX(iFock);

    auto gradw12a_y = rw2DensityGrid.alphaDensityGradientY(iFock);

    auto gradw12a_z = rw2DensityGrid.alphaDensityGradientZ(iFock);

    timer.stop("Kxc matrix prep.");

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Kxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);

    CDenseMatrix mat_G_gga(naos, npoints);

    auto G_val = mat_G.values();

    auto G_gga_val = mat_G_gga.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(weights, \
                    df2000, df1100, df1010, df1001, df10001, df0020, df0011, df00101, df00002, \
                    df0010, df00001, df00011, df01001, df0110, df3000, df2100, df1200, df2010, \
                    df0030, df0021, df0012, df00201, df00111, df00102, df00003, df2001, df1110, \
                    df1101, df20001, df11001, df1020, df1011, df1002, df10101, df10002, df01002, \
                    df0120, df0111, df01101, df10011, df01011, df0210, df02001, df00021, \
                    ngrada, grada_x, grada_y, grada_z, rhow1rhow2, rxw1rhow2, ryw1rhow2, rzw1rhow2, \
                    rxw1rxw2, rxw1ryw2, rxw1rzw2, ryw1rxw2, ryw1ryw2, ryw1rzw2, rzw1rxw2, rzw1ryw2, rzw1rzw2, \
                    rhow12a, gradw12a_x, gradw12a_y, gradw12a_z, \
                    G_val, G_gga_val, chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                double w = weights[gridblockpos + g];

                double znva = 1.0 / ngrada[g];

                double znva3 = 1.0 / std::pow(ngrada[g], 3.0);

                double znva5 = 1.0 / std::pow(ngrada[g], 5.0);

                double rxw12a = gradw12a_x[g];

                double ryw12a = gradw12a_y[g];

                double rzw12a = gradw12a_z[g];

                double xigrad_x = znva * grada_x[g];

                double xigrad_y = znva * grada_y[g];

                double xigrad_z = znva * grada_z[g];

                double xigrad_xx = (znva - grada_x[g] * grada_x[g] * znva3);

                double xigrad_yy = (znva - grada_y[g] * grada_y[g] * znva3);

                double xigrad_zz = (znva - grada_z[g] * grada_z[g] * znva3);

                double xigrad_xy = -grada_x[g] * grada_y[g] * znva3;

                double xigrad_xz = -grada_x[g] * grada_z[g] * znva3;

                double xigrad_yz = -grada_y[g] * grada_z[g] * znva3;

                double xigrad_xxy = 3.0 * grada_x[g] * grada_x[g] * grada_y[g] * znva5 - grada_y[g] * znva3;

                double xigrad_xxz = 3.0 * grada_x[g] * grada_x[g] * grada_z[g] * znva5 - grada_z[g] * znva3;

                double xigrad_xyy = 3.0 * grada_x[g] * grada_y[g] * grada_y[g] * znva5 - grada_x[g] * znva3;

                double xigrad_xzz = 3.0 * grada_x[g] * grada_z[g] * grada_z[g] * znva5 - grada_x[g] * znva3;

                double xigrad_yzz = 3.0 * grada_y[g] * grada_z[g] * grada_z[g] * znva5 - grada_y[g] * znva3;

                double xigrad_yyz = 3.0 * grada_y[g] * grada_y[g] * grada_z[g] * znva5 - grada_z[g] * znva3;

                double xigrad_xyz = 3.0 * grada_x[g] * grada_y[g] * grada_z[g] * znva5;

                double xigrad_xxx = 3.0 * grada_x[g] * grada_x[g] * grada_x[g] * znva5 - 3.0 * grada_x[g] * znva3;

                double xigrad_yyy = 3.0 * grada_y[g] * grada_y[g] * grada_y[g] * znva5 - 3.0 * grada_y[g] * znva3;

                double xigrad_zzz = 3.0 * grada_z[g] * grada_z[g] * grada_z[g] * znva5 - 3.0 * grada_z[g] * znva3;

                // Various required quantities

                double xigrad_dot_rw12a = (xigrad_x * rxw12a + xigrad_y * ryw12a + xigrad_z * rzw12a);

                double xigrad_dot_rw1rw2 = xigrad_x * rxw1rhow2[g] + xigrad_y * ryw1rhow2[g] + xigrad_z * rzw1rhow2[g];

                double rw1_dot_rw2 = rxw1rxw2[g] + ryw1ryw2[g] + rzw1rzw2[g];

                double xigrad_dot_rw1rhow2 = xigrad_x * rxw1rhow2[g] + xigrad_y * ryw1rhow2[g] + xigrad_z * rzw1rhow2[g];

                double grad_dot_rw12 = grada_x[g] * rxw12a + grada_y[g] * ryw12a + grada_z[g] * rzw12a;

                double grad_dot_rw1rw2 = grada_x[g] * rxw1rhow2[g] + grada_y[g] * ryw1rhow2[g] +
                                         grada_z[g] * rzw1rhow2[g];

                double grad_dot_rw1rhow2 = grada_x[g] * rxw1rhow2[g] + grada_y[g] * ryw1rhow2[g] +
                                           grada_z[g] * rzw1rhow2[g];

                double xigrad_dot_rw1_xigrad_dot_rw2 = xigrad_x * xigrad_x * rxw1rxw2[g] + xigrad_x * xigrad_y * rxw1ryw2[g] +
                                                       xigrad_x * xigrad_z * rxw1rzw2[g] + xigrad_y * xigrad_x * ryw1rxw2[g] +
                                                       xigrad_y * xigrad_y * ryw1ryw2[g] + xigrad_y * xigrad_z * ryw1rzw2[g] +
                                                       xigrad_z * xigrad_x * rzw1rxw2[g] + xigrad_z * xigrad_y * rzw1ryw2[g] +
                                                       xigrad_z * xigrad_z * rzw1rzw2[g];

                // scalar contribution

                double prefac = 0.0;

                // first, second

                prefac += (df2000[g] * rhow12a[g] + df1100[g] * rhow12a[g]) +

                          (df1010[g] + df1001[g]) * xigrad_dot_rw12a +

                          df10001[g] * 2.0 * grad_dot_rw12;

                // fifth

                prefac += (df3000[g] + 2.0 * df2100[g] + df1200[g]) * rhow1rhow2[g];

                // seventh

                prefac += (df2010[g] + df2001[g]) * xigrad_dot_rw1rw2 +

                          (df1110[g] + df1101[g]) * xigrad_dot_rw1rw2 +

                          2.0 * (df20001[g] + df11001[g]) * grad_dot_rw1rw2;

                // eighth

                prefac += (df1020[g] + 2.0 * df1011[g] + df1002[g]) * xigrad_dot_rw1_xigrad_dot_rw2;

                prefac += (df1010[g] + df1001[g]) *

                          (xigrad_xx * rxw1rxw2[g] + xigrad_xy * rxw1ryw2[g] + xigrad_xz * rxw1rzw2[g] +

                           xigrad_xy * ryw1rxw2[g] + xigrad_yy * ryw1ryw2[g] + xigrad_yz * ryw1rzw2[g] +

                           xigrad_xz * rzw1rxw2[g] + xigrad_yz * rzw1ryw2[g] + xigrad_zz * rzw1rzw2[g]);

                prefac += 2.0 * (df10101[g] + df10101[g]) * ngrada[g] * xigrad_dot_rw1_xigrad_dot_rw2;

                prefac += 4.0 * df10002[g] * ngrada[g] * ngrada[g] * xigrad_dot_rw1_xigrad_dot_rw2;

                prefac += 2.0 * df10001[g] * rw1_dot_rw2;

                G_val[nu_offset + g] = w * prefac * chi_val[nu_offset + g];

                // vector contribution

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                // third

                prefac = (df1010[g] + df0110[g]) * rhow12a[g];

                xcomp += prefac * xigrad_x;

                ycomp += prefac * xigrad_y;

                zcomp += prefac * xigrad_z;

                prefac = (df10001[g] + df01001[g]) * rhow12a[g];

                xcomp += prefac * grada_x[g];

                ycomp += prefac * grada_y[g];

                zcomp += prefac * grada_z[g];

                // fourth

                xcomp += df0010[g] * (xigrad_xx * rxw12a + xigrad_xy * ryw12a + xigrad_xz * rzw12a);

                ycomp += df0010[g] * (xigrad_xy * rxw12a + xigrad_yy * ryw12a + xigrad_yz * rzw12a);

                zcomp += df0010[g] * (xigrad_xz * rxw12a + xigrad_yz * ryw12a + xigrad_zz * rzw12a);

                xcomp += df00001[g] * rxw12a;

                ycomp += df00001[g] * ryw12a;

                zcomp += df00001[g] * rzw12a;

                // tenth

                prefac = (df1010[g] + df0110[g]);

                xcomp += prefac * (xigrad_xx * rxw1rhow2[g] + xigrad_xy * ryw1rhow2[g] + xigrad_xz * rzw1rhow2[g]);

                ycomp += prefac * (xigrad_xy * rxw1rhow2[g] + xigrad_yy * ryw1rhow2[g] + xigrad_yz * rzw1rhow2[g]);

                zcomp += prefac * (xigrad_xz * rxw1rhow2[g] + xigrad_yz * ryw1rhow2[g] + xigrad_zz * rzw1rhow2[g]);

                // twelfth first

                prefac = df0010[g];

                xcomp += prefac * (xigrad_xxx * rxw1rxw2[g] + xigrad_xxy * rxw1ryw2[g] + xigrad_xxz * rxw1rzw2[g] +

                                   xigrad_xxy * ryw1rxw2[g] + xigrad_xyy * ryw1ryw2[g] + xigrad_xyz * ryw1rzw2[g] +

                                   xigrad_xxz * rzw1rxw2[g] + xigrad_xyz * rzw1ryw2[g] + xigrad_xzz * rzw1rzw2[g]);

                ycomp += prefac * (xigrad_xxy * rxw1rxw2[g] + xigrad_xyy * rxw1ryw2[g] + xigrad_xyz * rxw1rzw2[g] +

                                   xigrad_xyy * ryw1rxw2[g] + xigrad_yyy * ryw1ryw2[g] + xigrad_yyz * ryw1rzw2[g] +

                                   xigrad_xyz * rzw1rxw2[g] + xigrad_yyz * rzw1ryw2[g] + xigrad_yzz * rzw1rzw2[g]);

                zcomp += prefac * (xigrad_xxz * rxw1rxw2[g] + xigrad_xyz * rxw1ryw2[g] + xigrad_xzz * rxw1rzw2[g] +

                                   xigrad_xyz * ryw1rxw2[g] + xigrad_yyz * ryw1ryw2[g] + xigrad_yzz * ryw1rzw2[g] +

                                   xigrad_xzz * rzw1rxw2[g] + xigrad_yzz * rzw1ryw2[g] + xigrad_zzz * rzw1rzw2[g]);

                // twelfth second

                prefac = (df00101[g] + df00011[g]) * ngrada[g];

                prefac += (df0020[g] + df0011[g]);

                xcomp += prefac * xigrad_x * (xigrad_xx * rxw1rxw2[g] + xigrad_xy * rxw1ryw2[g] + xigrad_xz * rxw1rzw2[g] +

                                              xigrad_xy * ryw1rxw2[g] + xigrad_yy * ryw1ryw2[g] + xigrad_yz * ryw1rzw2[g] +

                                              xigrad_xz * rzw1rxw2[g] + xigrad_yz * rzw1ryw2[g] + xigrad_zz * rzw1rzw2[g]);

                ycomp += prefac * xigrad_y * (xigrad_xx * rxw1rxw2[g] + xigrad_xy * rxw1ryw2[g] + xigrad_xz * rxw1rzw2[g] +

                                              xigrad_xy * ryw1rxw2[g] + xigrad_yy * ryw1ryw2[g] + xigrad_yz * ryw1rzw2[g] +

                                              xigrad_xz * rzw1rxw2[g] + xigrad_yz * rzw1ryw2[g] + xigrad_zz * rzw1rzw2[g]);

                zcomp += prefac * xigrad_z * (xigrad_xx * rxw1rxw2[g] + xigrad_xy * rxw1ryw2[g] + xigrad_xz * rxw1rzw2[g] +

                                              xigrad_xy * ryw1rxw2[g] + xigrad_yy * ryw1ryw2[g] + xigrad_yz * ryw1rzw2[g] +

                                              xigrad_xz * rzw1rxw2[g] + xigrad_yz * rzw1ryw2[g] + xigrad_zz * rzw1rzw2[g]);

                // twelfth third

                prefac = (df0020[g] + df0011[g]);

                prefac += df00101[g] * ngrada[g];

                xcomp += prefac * (xigrad_xx * xigrad_x * (rxw1rxw2[g] + rxw1rxw2[g]) +

                                   xigrad_xx * xigrad_y * (ryw1rxw2[g] + rxw1ryw2[g]) +

                                   xigrad_xx * xigrad_z * (rzw1rxw2[g] + rxw1rzw2[g]) +

                                   xigrad_xy * xigrad_x * (rxw1ryw2[g] + ryw1rxw2[g]) +

                                   xigrad_xy * xigrad_y * (ryw1ryw2[g] + ryw1ryw2[g]) +

                                   xigrad_xy * xigrad_z * (rzw1ryw2[g] + ryw1rzw2[g]) +

                                   xigrad_xz * xigrad_x * (rxw1rzw2[g] + rzw1rxw2[g]) +

                                   xigrad_xz * xigrad_y * (ryw1rzw2[g] + rzw1ryw2[g]) +

                                   xigrad_xz * xigrad_z * (rzw1rzw2[g] + rzw1rzw2[g]));

                ycomp += prefac * (xigrad_xy * xigrad_x * (rxw1rxw2[g] + rxw1rxw2[g]) +

                                   xigrad_xy * xigrad_y * (ryw1rxw2[g] + rxw1ryw2[g]) +

                                   xigrad_xy * xigrad_z * (rzw1rxw2[g] + rxw1rzw2[g]) +

                                   xigrad_yy * xigrad_x * (rxw1ryw2[g] + ryw1rxw2[g]) +

                                   xigrad_yy * xigrad_y * (ryw1ryw2[g] + ryw1ryw2[g]) +

                                   xigrad_yy * xigrad_z * (rzw1ryw2[g] + ryw1rzw2[g]) +

                                   xigrad_yz * xigrad_x * (rxw1rzw2[g] + rzw1rxw2[g]) +

                                   xigrad_yz * xigrad_y * (ryw1rzw2[g] + rzw1ryw2[g]) +

                                   xigrad_yz * xigrad_z * (rzw1rzw2[g] + rzw1rzw2[g]));

                zcomp += prefac * (xigrad_xz * xigrad_x * (rxw1rxw2[g] + rxw1rxw2[g]) +

                                   xigrad_xz * xigrad_y * (ryw1rxw2[g] + rxw1ryw2[g]) +

                                   xigrad_xz * xigrad_z * (rzw1rxw2[g] + rxw1rzw2[g]) +

                                   xigrad_yz * xigrad_x * (rxw1ryw2[g] + ryw1rxw2[g]) +

                                   xigrad_yz * xigrad_y * (ryw1ryw2[g] + ryw1ryw2[g]) +

                                   xigrad_yz * xigrad_z * (rzw1ryw2[g] + ryw1rzw2[g]) +

                                   xigrad_zz * xigrad_x * (rxw1rzw2[g] + rzw1rxw2[g]) +

                                   xigrad_zz * xigrad_y * (ryw1rzw2[g] + rzw1ryw2[g]) +

                                   xigrad_zz * xigrad_z * (rzw1rzw2[g] + rzw1rzw2[g]));

                // twelfth fourth gam

                prefac = (df00101[g] + df00011[g]);

                prefac += df00002[g] * ngrada[g];

                xcomp += prefac * (xigrad_x * rxw1rxw2[g] + xigrad_y * ryw1rxw2[g] + xigrad_z * rzw1rxw2[g]);

                ycomp += prefac * (xigrad_x * rxw1ryw2[g] + xigrad_y * ryw1ryw2[g] + xigrad_z * rzw1ryw2[g]);

                zcomp += prefac * (xigrad_x * rxw1rzw2[g] + xigrad_y * ryw1rzw2[g] + xigrad_z * rzw1rzw2[g]);

                // twelfth fifth gam

                double twelfthfifth_gam = (xigrad_x * grada_x[g] + grada_x[g] * xigrad_x) * rxw1rxw2[g] +

                                          (xigrad_x * grada_y[g] + grada_x[g] * xigrad_y) * rxw1ryw2[g] +

                                          (xigrad_x * grada_z[g] + grada_x[g] * xigrad_z) * rxw1rzw2[g] +

                                          (xigrad_y * grada_x[g] + grada_y[g] * xigrad_x) * rxw1rxw2[g] +

                                          (xigrad_y * grada_y[g] + grada_y[g] * xigrad_y) * rxw1ryw2[g] +

                                          (xigrad_y * grada_z[g] + grada_y[g] * xigrad_z) * rxw1rzw2[g] +

                                          (xigrad_z * grada_x[g] + grada_z[g] * xigrad_x) * rxw1rxw2[g] +

                                          (xigrad_z * grada_y[g] + grada_z[g] * xigrad_y) * rxw1ryw2[g] +

                                          (xigrad_z * grada_z[g] + grada_z[g] * xigrad_z) * rxw1rzw2[g];

                // fourth, ninth, tenth

                // xigrad_dot_omega == (xigrad_x * xomega + xigrad_y * yomega + xigrad_z * zomega);

                prefac = df00101[g] * (2.0 * grad_dot_rw12);

                prefac += (df0020[g] + df0011[g]) * xigrad_dot_rw12a;

                prefac += (df2010[g] + 2.0 * df1110[g] + df0210[g]) * rhow1rhow2[g];

                prefac += (df1020[g] + df1011[g] + df0120[g] + df0111[g]) * xigrad_dot_rw1rhow2;

                prefac += (df10101[g] + df10011[g] + df01101[g] + df0111[g]) * grad_dot_rw1rhow2;

                prefac += (df0030[g] + 2.0 * df0021[g] + df0012[g]) * xigrad_dot_rw1_xigrad_dot_rw2;

                prefac += (df00101[g] + df00011[g]) * rw1_dot_rw2;

                prefac += (df00201[g] + df00111[g]) * twelfthfifth_gam;

                prefac += df00102[g] * xigrad_dot_rw1_xigrad_dot_rw2;

                xcomp += prefac * xigrad_x;

                ycomp += prefac * xigrad_y;

                zcomp += prefac * xigrad_z;

                // fourth, ninth, tenth

                // grad_dot_omega == grada_x[g] * xomega + grada_y[g] * yomega + grada_z[g] * zomega;

                prefac = (df00101[g] * xigrad_dot_rw12a + df00011[g] * xigrad_dot_rw12a);

                prefac += df00002[g] * (2.0 * grad_dot_rw12);

                prefac += (df20001[g] + 2.0 * df11001[g] + df02001[g]) * rhow1rhow2[g];

                prefac += (df10101[g] + df10011[g] + df01101[g] + df0111[g] + df01011[g]) * xigrad_dot_rw1rhow2;

                prefac += (df10002[g] + df01002[g]) * grad_dot_rw1rhow2;

                prefac += df00002[g] * rw1_dot_rw2;

                prefac += (df00201[g] + 2 * df00111[g] + df00021[g]) * xigrad_dot_rw1_xigrad_dot_rw2;

                prefac += (df00102[g] + df00011[g]) * ngrada[g] * xigrad_dot_rw1_xigrad_dot_rw2;

                prefac += df00003[g] * ngrada[g] * ngrada[g];

                xcomp += prefac * grada_x[g];

                ycomp += prefac * grada_y[g];

                zcomp += prefac * grada_z[g];

                // tenth

                // omega_dot_rw1rhow2 == xomega * rxw1rhow2[g] + yomega * ryw1rhow2[g] + zomega * rzw1rhow2[g];

                prefac = (df10001[g] + df01001[g]);

                xcomp += prefac * rxw1rhow2[g];

                ycomp += prefac * ryw1rhow2[g];

                zcomp += prefac * rzw1rhow2[g];

                G_gga_val[nu_offset + g] = w * (xcomp * chi_x_val[nu_offset + g] +

                                                ycomp * chi_y_val[nu_offset + g] +

                                                zcomp * chi_z_val[nu_offset + g]);
            }
        }
    }

    timer.stop("Kxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Kxc matrix matmul");

    auto mat_Kxc = denblas::multABt(gtoValues, mat_G);

    auto mat_Kxc_gga = denblas::multABt(gtoValues, mat_G_gga);

    mat_Kxc_gga.symmetrize();  // matrix + matrix.T

    mat_Kxc = denblas::addAB(mat_Kxc, mat_Kxc_gga, 1.0);

    timer.stop("Kxc matrix matmul");

    return mat_Kxc;
}

CDenseMatrix
CXCNewIntegrator::computeGtoValuesOnGridPoints(const CMolecule&       molecule,
                                               const CMolecularBasis& basis,
                                               const CMolecularGrid&  molecularGrid) const
{
    std::string errnotpartitioned("CXCNewIntegrator.computeGtoValuesOnGridPoints: Cannot use unpartitioned molecular grid");

    errors::assertMsgCritical(molecularGrid.isPartitioned(), errnotpartitioned);

    auto nthreads = omp_get_max_threads();

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // memory blocks for GTOs on grid points

    CDenseMatrix allgtovalues(naos, molecularGrid.getNumberOfGridPoints());

    CMemBlock2D<double> gaos(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos); // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

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

        auto boxdim = gtoeval::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // pre-screening of GTOs

        gtoeval::preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 0, _screeningThresholdForGTOValues, boxdim);  // 0th order GTO derivative

        // GTO values on grid points

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForLDA(gaos, gtovec, xcoords, ycoords, zcoords, gridblockpos,

                                             grid_batch_offset, grid_batch_size, skip_cgto_ids);
        }

        int32_t aocount = 0;

        for (int32_t nu = 0; nu < naos; nu++)
        {
            if (skip_ao_ids.data()[nu]) continue;

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

        for (int32_t i = 0; i < aocount; i++)
        {
            auto aoidx = aoinds[i];

            std::memcpy(allgtovalues.row(aoidx) + gridblockpos, gaos.data(aoidx), npoints * sizeof(double));
        }
    }

    // auto npoints = molecularGrid.getNumberOfGridPoints();
    // CMemBlock2D<double> refgaos(npoints, naos);
    // gtorec::computeGtosValuesForLDA(refgaos, gtovec, xcoords, ycoords, zcoords, 0, 0, npoints);
    // double maxdiff = 0.0;
    // for (int32_t i = 0; i < naos; i++)
    // {
    //     for (int32_t g = 0; g < npoints; g++)
    //     {
    //         maxdiff = std::max(maxdiff, std::fabs(refgaos.data(i)[g] - allgtovalues.row(i)[g]));
    //     }
    // }
    // std::cout << "maxdiff: " << maxdiff << std::endl;

    // destroy GTOs container

    delete gtovec;

    return allgtovalues;
}
