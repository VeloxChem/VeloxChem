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
    // TODO: use _numberOfPointsThreshold

    CMemBlock2D<double> gaos(1024, naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos); // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

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

        auto boxdim = _getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // pre-screening of GTOs

        timer.start("GTO pre-screening");

        _preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 0, boxdim);

        timer.stop("GTO pre-screening");

        // GTO values on grid points

        timer.start("OMP GTO evaluation");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            omptimers[thread_id].start("gtoeval");

            auto batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForLDA(gaos, gtovec, xcoords, ycoords, zcoords, gridblockpos,
                                             batch_offset, batch_size, skip_cgto_ids);

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

        auto sub_dens_mat = _getSubDensityMatrix(densityMatrix, aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto dengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, sub_dens_mat, xcfuntype, timer);

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        CXCGradientGrid vxcgrid(npoints, dengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxcgrid, dengrid);

        timer.stop("XC functional eval.");

        // compute partial contribution to Vxc matrix

        auto partial_mat_Vxc = _integratePartialVxcFockForLDA(gridblockpos, npoints, weights, mat_chi, vxcgrid, timer);

        // distribute partial Vxc to full Vxc

        timer.start("Vxc matrix dist.");

        _distributeVxcMatrix(mat_Vxc, partial_mat_Vxc, aoinds, aocount, naos);

        timer.stop("Vxc matrix dist.");

        // compute partial contribution to XC energy

        timer.start("XC energy");

        auto rhoa = dengrid.alphaDensity(0);

        auto rhob = dengrid.betaDensity(0);

        auto efunc = vxcgrid.xcFunctionalValues();

        for (int32_t g = 0; g < npoints; g++)
        {
            nele += weights[gridblockpos + g] * (rhoa[g] + rhob[g]);

            xcene += weights[gridblockpos + g] * efunc[g];
        }

        timer.stop("XC energy");
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    std::cout << "Timing of new integrator" << std::endl;

    std::cout << "------------------------" << std::endl;

    std::cout << timer.getSummary() << std::endl;

    std::cout << "OpenMP timing" << std::endl;

    for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
    {
        std::cout << "Thread " << thread_id << std::endl;

        std::cout << omptimers[thread_id].getSummary() << std::endl;
    }

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
    // TODO: use _numberOfPointsThreshold

    CMemBlock2D<double> gaos(1024, naos);

    CMemBlock2D<double> gaox(1024, naos);

    CMemBlock2D<double> gaoy(1024, naos);

    CMemBlock2D<double> gaoz(1024, naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos); // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

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

        auto boxdim = _getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // pre-screening of GTOs

        timer.start("GTO pre-screening");

        _preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 1, boxdim);

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

        auto sub_dens_mat = _getSubDensityMatrix(densityMatrix, aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto dengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat, xcfuntype, timer);

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        CXCGradientGrid vxcgrid(npoints, dengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxcgrid, dengrid);

        timer.stop("XC functional eval.");

        // compute partial contribution to Vxc matrix

        auto partial_mat_Vxc = _integratePartialVxcFockForGGA(gridblockpos, npoints, weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                              vxcgrid, dengrid, timer);

        // distribute partial Vxc to full Vxc

        timer.start("Vxc matrix dist.");

        _distributeVxcMatrix(mat_Vxc, partial_mat_Vxc, aoinds, aocount, naos);

        timer.stop("Vxc matrix dist.");

        // compute partial contribution to XC energy

        timer.start("XC energy");

        auto rhoa = dengrid.alphaDensity(0);

        auto rhob = dengrid.betaDensity(0);

        auto efunc = vxcgrid.xcFunctionalValues();

        for (int32_t g = 0; g < npoints; g++)
        {
            nele += weights[gridblockpos + g] * (rhoa[g] + rhob[g]);

            xcene += weights[gridblockpos + g] * efunc[g];
        }

        timer.stop("XC energy");
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    std::cout << "Timing of new integrator" << std::endl;

    std::cout << "------------------------" << std::endl;

    std::cout << timer.getSummary() << std::endl;

    std::cout << "OpenMP timing" << std::endl;

    for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
    {
        std::cout << "Thread " << thread_id << std::endl;

        std::cout << omptimers[thread_id].getSummary() << std::endl;
    }

    mat_Vxc.setNumberOfElectrons(nele);

    mat_Vxc.setExchangeCorrelationEnergy(xcene);

    return mat_Vxc;
}

std::array<double, 6>
CXCNewIntegrator::_getGridBoxDimension(const int32_t gridBlockPosition,
                                       const int32_t nGridPoints,
                                       const double* xcoords,
                                       const double* ycoords,
                                       const double* zcoords) const
{
    double xmin = xcoords[gridBlockPosition], ymin = ycoords[gridBlockPosition], zmin = zcoords[gridBlockPosition];

    double xmax = xcoords[gridBlockPosition], ymax = ycoords[gridBlockPosition], zmax = zcoords[gridBlockPosition];

    for (int32_t g = 0; g < nGridPoints; g++)
    {
        xmin = std::min(xmin, xcoords[gridBlockPosition + g]);

        ymin = std::min(ymin, ycoords[gridBlockPosition + g]);

        zmin = std::min(zmin, zcoords[gridBlockPosition + g]);

        xmax = std::max(xmax, xcoords[gridBlockPosition + g]);

        ymax = std::max(ymax, ycoords[gridBlockPosition + g]);

        zmax = std::max(zmax, zcoords[gridBlockPosition + g]);
    }

    return std::array<double, 6>({xmin, ymin, zmin, xmax, ymax, zmax});
}

void
CXCNewIntegrator::_preScreenGtos(CMemBlock<int32_t>&          skipCgtoIds,
                                 CMemBlock<int32_t>&          skipAOIds,
                                 const CGtoContainer*         gtoContainer,
                                 const int32_t                gtoDeriv,
                                 const std::array<double, 6>& boxDimension) const
{
    skipCgtoIds.zero();

    skipAOIds.zero();

    double xmin = boxDimension[0], ymin = boxDimension[1], zmin = boxDimension[2];

    double xmax = boxDimension[3], ymax = boxDimension[4], zmax = boxDimension[5];

    for (int32_t i = 0, cgto_count = 0; i < gtoContainer->getNumberOfGtoBlocks(); i++)
    {
        auto bgtos = gtoContainer->getGtoBlock(i);

        auto bang = bgtos.getAngularMomentum();

        auto bfnorms = bgtos.getNormFactors();

        auto bfexps = bgtos.getExponents();

        auto bfx = bgtos.getCoordinatesX();

        auto bfy = bgtos.getCoordinatesY();

        auto bfz = bgtos.getCoordinatesZ();

        auto spos = bgtos.getStartPositions();

        auto epos = bgtos.getEndPositions();

        // loop over contracted GTOs

        for (int32_t j = 0; j < bgtos.getNumberOfContrGtos(); j++, cgto_count++)
        {
            // contracted GTO screening

            auto firstprim = spos[j];

            double rx = std::max({xmin - bfx[firstprim], bfx[firstprim] - xmax, 0.0});

            double ry = std::max({ymin - bfy[firstprim], bfy[firstprim] - ymax, 0.0});

            double rz = std::max({zmin - bfz[firstprim], bfz[firstprim] - zmax, 0.0});

            auto r2 = rx * rx + ry * ry + rz * rz;

            if (r2 > 1.0)
            {
                auto minexp = bfexps[firstprim];

                auto maxexp = bfexps[firstprim];

                auto maxcoef = std::fabs(bfnorms[firstprim]);

                for (int32_t iprim = spos[j]; iprim < epos[j]; iprim++)
                {
                    auto bexp = bfexps[iprim];

                    auto bnorm = std::fabs(bfnorms[iprim]);

                    minexp = std::min(minexp, bexp);

                    maxexp = std::max(maxexp, bexp);

                    maxcoef = std::max(maxcoef, bnorm);
                }

                // gto  :           r^{ang}   |C| exp(-alpha r^2)
                // gto_m:           r^{ang-1} |C| exp(-alpha r^2)
                // gto_p: (2 alpha) r^{ang+1} |C| exp(-alpha r^2)

                auto gtolimit = maxcoef * std::exp(-minexp * r2);

                auto r = std::sqrt(r2);

                for (int32_t ipow = 0; ipow < bang; ipow++) gtolimit *= r;

                if (gtoDeriv == 1) gtolimit = std::max(gtolimit, 2.0 * maxexp * r * gtolimit);

                if (gtolimit < _screeningThresholdForGTOValues)
                {
                    skipCgtoIds.data()[cgto_count] = 1;

                    auto bnspher = angmom::to_SphericalComponents(bang);

                    for (int32_t k = 0; k < bnspher; k++)
                    {
                        auto ao_idx = (bgtos.getIdentifiers(k))[j];

                        skipAOIds.data()[ao_idx] = 1;
                    }
                }
            }
        }
    }
}

CAODensityMatrix
CXCNewIntegrator::_getSubDensityMatrix(const CAODensityMatrix&     densityMatrix,
                                       const std::vector<int32_t>& aoIndices,
                                       const int32_t               aoCount,
                                       const int32_t               nAOs) const
{
    if (aoCount < nAOs)
    {
        CDenseMatrix sub_dens(aoCount, aoCount);

        const CDenseMatrix& dens = densityMatrix.getReferenceToDensity(0);

        for (int32_t i = 0; i < aoCount; i++)
        {
            auto sub_dens_row = sub_dens.row(i);

            auto dens_row = dens.row(aoIndices[i]);

            for (int32_t j = 0; j < aoCount; j++)
            {
                sub_dens_row[j] = dens_row[aoIndices[j]];
            }
        }

        return CAODensityMatrix({sub_dens}, denmat::rest);
    }
    else
    {
        return densityMatrix;
    }
}

CDenseMatrix
CXCNewIntegrator::_integratePartialVxcFockForLDA(const int32_t          gridBlockPosition,
                                                 const int32_t          npoints,
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
            G_nu[g] = weights[gridBlockPosition + g] * grhoa[g] * chi_nu[g];
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

void
CXCNewIntegrator::_distributeVxcMatrix(CAOKohnShamMatrix&          matVxc,
                                       const CDenseMatrix&         partialMatVxc,
                                       const std::vector<int32_t>& aoIndices,
                                       const int32_t               aoCount,
                                       const int32_t               nAOs) const
{
    if (aoCount < nAOs)
    {
        for (int32_t row = 0; row < partialMatVxc.getNumberOfRows(); row++)
        {
            auto row_orig = aoIndices[row];

            auto matVxc_row_orig = matVxc.getMatrix(0) + row_orig * nAOs;

            auto partialMatVxc_row = partialMatVxc.row(row);

            for (int32_t col = 0; col < partialMatVxc.getNumberOfColumns(); col++)
            {
                auto col_orig = aoIndices[col];

                matVxc_row_orig[col_orig] += partialMatVxc_row[col];
            }
        }
    }
    else
    {
        matVxc.addMatrixContribution(partialMatVxc);
    }
}
