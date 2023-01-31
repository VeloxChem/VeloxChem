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
#include "GridScreener.hpp"
#include "GtoEvaluator.hpp"
#include "NewFunctionalParser.hpp"
#include "SubMatrix.hpp"
#include "XCFuncType.hpp"

CXCNewIntegrator::CXCNewIntegrator(MPI_Comm comm)

    : _screeningThresholdForGTOValues(1.0e-12)

    , _screeningThresholdForDensityValues(1.0e-13)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;
}

CAOKohnShamMatrix
CXCNewIntegrator::integrateVxcFock(const CMolecule&        molecule,
                                   const CMolecularBasis&  basis,
                                   const CAODensityMatrix& densityMatrix,
                                   const CMolecularGrid&   molecularGrid,
                                   const std::string&      xcFuncLabel) const
{
    auto newfvxc = newvxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = newfvxc.getFunctionalType();

    auto flag = densityMatrix.isClosedShell() ? std::string("CLOSEDSHELL") : std::string("OPENSHELL");

    if (xcfuntype == xcfun::lda)
    {
        return _integrateVxcFockForLDA(molecule, basis, densityMatrix, molecularGrid, newfvxc, flag);
    }
    else if (xcfuntype == xcfun::gga)
    {
        return _integrateVxcFockForGGA(molecule, basis, densityMatrix, molecularGrid, newfvxc, flag);
    }
    else if (xcfuntype == xcfun::mgga)
    {
        return _integrateVxcFockForMGGA(molecule, basis, densityMatrix, molecularGrid, newfvxc, flag);
    }

    std::string errxcfuntype("XCNewIntegrator.integrateVxcFock: Only implemented for LDA/GGA/meta-GGA");

    errors::assertMsgCritical(false, errxcfuntype);

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
    auto fvxc = newvxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

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
        else if (xcfuntype == xcfun::mgga)
        {
            _integrateFxcFockForMGGA(aoFockMatrix, molecule, basis, rwDensityMatrix, gsDensityMatrix, molecularGrid, fvxc);
        }
        else
        {
            std::string errxcfuntype("XCNewIntegrator.integrateFxcFock: Only implemented for LDA/GGA/meta-GGA");

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
    auto fvxc = newvxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

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
        else if (xcfuntype == xcfun::mgga)
        {
            _integrateKxcFockForMGGA(aoFockMatrix, molecule, basis, rwDensityMatrix, rw2DensityMatrix, gsDensityMatrix,
                                     molecularGrid, fvxc, quadMode);
        }
        else
        {
            std::string errxcfuntype("XCNewIntegrator.integrateKxcFock: Only implemented for LDA/GGA/meta-GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCNewIntegrator.integrateKxcFock: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }
}

void
CXCNewIntegrator::integrateLxcFock(CAOFockMatrix&          aoFockMatrix,
                                   const CMolecule&        molecule,
                                   const CMolecularBasis&  basis,
                                   const CAODensityMatrix& rwDensityMatrix,
                                   const CAODensityMatrix& rw2DensityMatrix,
                                   const CAODensityMatrix& rw3DensityMatrix,
                                   const CAODensityMatrix& gsDensityMatrix,
                                   const CMolecularGrid&   molecularGrid,
                                   const std::string&      xcFuncLabel,
                                   const std::string&      cubeMode) const
{
    auto fvxc = newvxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (rwDensityMatrix.isClosedShell() && rw2DensityMatrix.isClosedShell() && rw3DensityMatrix.isClosedShell() && gsDensityMatrix.isClosedShell())
    {
        if (xcfuntype == xcfun::lda)
        {
            _integrateLxcFockForLDA(aoFockMatrix, molecule, basis, rwDensityMatrix, rw2DensityMatrix, rw3DensityMatrix,
                                    gsDensityMatrix, molecularGrid, fvxc, cubeMode);
        }
        else if (xcfuntype == xcfun::gga)
        {
            _integrateLxcFockForGGA(aoFockMatrix, molecule, basis, rwDensityMatrix, rw2DensityMatrix, rw3DensityMatrix,
                                    gsDensityMatrix, molecularGrid, fvxc, cubeMode);
        }
        else
        {
            std::string errxcfuntype("XCNewIntegrator.integrateLxcFock: Only implemented for LDA/GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCNewIntegrator.integrateLxcFock: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }
}

void
CXCNewIntegrator::integrateKxcLxcFock(CAOFockMatrix&          aoFockMatrix,
                                      const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const CAODensityMatrix& rwDensityMatrix,
                                      const CAODensityMatrix& rw2DensityMatrix,
                                      const CAODensityMatrix& rw3DensityMatrix,
                                      const CAODensityMatrix& gsDensityMatrix,
                                      const CMolecularGrid&   molecularGrid,
                                      const std::string&      xcFuncLabel,
                                      const std::string&      cubeMode) const
{
    auto fvxc = newvxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (rwDensityMatrix.isClosedShell() && rw2DensityMatrix.isClosedShell() && rw3DensityMatrix.isClosedShell() && gsDensityMatrix.isClosedShell())
    {
        if (xcfuntype == xcfun::lda)
        {
            _integrateKxcLxcFockForLDA(aoFockMatrix, molecule, basis, rwDensityMatrix, rw2DensityMatrix, rw3DensityMatrix,
                                       gsDensityMatrix, molecularGrid, fvxc, cubeMode);
        }
        else if (xcfuntype == xcfun::gga)
        {
            _integrateKxcLxcFockForGGA(aoFockMatrix, molecule, basis, rwDensityMatrix, rw2DensityMatrix, rw3DensityMatrix,
                                       gsDensityMatrix, molecularGrid, fvxc, cubeMode);
        }
        else if (xcfuntype == xcfun::mgga)
        {
            _integrateKxcLxcFockForMGGA(aoFockMatrix, molecule, basis, rwDensityMatrix, rw2DensityMatrix,rw3DensityMatrix,
                                        gsDensityMatrix, molecularGrid, fvxc, cubeMode);
        }
        else
        {
            std::string errxcfuntype("XCNewIntegrator.integrateKxcLxcFock: Only implemented for LDA/GGA/meta-GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCNewIntegrator.integrateKxcLxcFock: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }
}

void
CXCNewIntegrator::integrateVxcPDFT(CAOKohnShamMatrix&      aoFockMatrix,
                                   CDense4DTensor&         moTwoBodyGradient,
                                   const CMolecule&        molecule,
                                   const CMolecularBasis&  basis,
                                   const CAODensityMatrix& DensityMatrix,
                                   const CDense4DTensor&   TwoBodyDensityMatrix,
                                   const CDenseMatrix&     ActiveMOs,
                                   const CMolecularGrid&   molecularGrid,
                                   const std::string&      xcFuncLabel) const
{
    auto fvxc = newvxcfuncs::getPairDensityExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (xcfuntype == "PLDA")
    {
        _integrateVxcPDFTForLDA(aoFockMatrix, moTwoBodyGradient, molecule, basis,
                                DensityMatrix, TwoBodyDensityMatrix, ActiveMOs, molecularGrid, fvxc);
    }
    else if (xcfuntype == "PGGA")
    {
        _integrateVxcPDFTForGGA(aoFockMatrix, moTwoBodyGradient, molecule, basis,
                                DensityMatrix, TwoBodyDensityMatrix, ActiveMOs, molecularGrid, fvxc);
    }
    else
    {
        std::string errxcfuntype("XCNewIntegrator.integrateVxcPDFT: Only implemented for LDA/GGA");

        errors::assertMsgCritical(false, errxcfuntype);
    }
}

CAOKohnShamMatrix
CXCNewIntegrator::_integrateVxcFockForLDA(const CMolecule&        molecule,
                                          const CMolecularBasis&  basis,
                                          const CAODensityMatrix& densityMatrix,
                                          const CMolecularGrid&   molecularGrid,
                                          const CXCNewFunctional& xcFunctional,
                                          const std::string&      flag) const
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

    bool closedshell = (fstr::upcase(flag) == std::string("CLOSEDSHELL"));

    CAOKohnShamMatrix mat_Vxc(densityMatrix.getNumberOfRows(0), densityMatrix.getNumberOfColumns(0), closedshell);

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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> exc_data(1 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();

    auto exc = exc_data.data();

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

        if (aocount == 0) continue;

        // generate sub density matrix and density grid

        if (closedshell)
        {
            timer.start("Density matrix slicing");

            auto sub_dens_mat = submat::getSubDensityMatrix(densityMatrix, 0, "ALPHA", aoinds, aocount, naos);

            timer.stop("Density matrix slicing");

            dengridgen::generateDensityForLDA(rho, npoints, mat_chi, sub_dens_mat, timer);
        }
        else
        {
            timer.start("Density matrix slicing");

            auto sub_dens_mat_a = submat::getSubDensityMatrix(densityMatrix, 0, "ALPHA", aoinds, aocount, naos);
            auto sub_dens_mat_b = submat::getSubDensityMatrix(densityMatrix, 0, "BETA", aoinds, aocount, naos);

            timer.stop("Density matrix slicing");

            dengridgen::generateDensityForLDA(rho, npoints, mat_chi, sub_dens_mat_a, sub_dens_mat_b, timer);
        }

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_lda(npoints, rho, exc, vrho);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenVxcFockForLDA(rho, exc, vrho, npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

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

        for (int32_t g = 0; g < npoints; g++)
        {
            auto rho_total = rho[2 * g + 0] + rho[2 * g + 1];

            nele += local_weights[g] * rho_total;

            xcene += local_weights[g] * exc[g] * rho_total;
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
                                          const CXCNewFunctional& xcFunctional,
                                          const std::string&      flag) const
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

    bool closedshell = (fstr::upcase(flag) == std::string("CLOSEDSHELL"));

    CAOKohnShamMatrix mat_Vxc(densityMatrix.getNumberOfRows(0), densityMatrix.getNumberOfColumns(0), closedshell);

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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rhograd_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> sigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> exc_data(1 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> vsigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();

    auto rhograd = rhograd_data.data();

    auto sigma = sigma_data.data();

    auto exc = exc_data.data();

    auto vrho = vrho_data.data();

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

        if (aocount == 0) continue;

        // generate sub density matrix and density grid

        if (closedshell)
        {
            timer.start("Density matrix slicing");

            auto sub_dens_mat = submat::getSubDensityMatrix(densityMatrix, 0, "ALPHA", aoinds, aocount, naos);

            dengridgen::generateDensityForGGA(rho, rhograd, sigma, npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                              sub_dens_mat, timer);

            timer.stop("Density matrix slicing");
        }
        else
        {
            timer.start("Density matrix slicing");

            auto sub_dens_mat_a = submat::getSubDensityMatrix(densityMatrix, 0, "ALPHA", aoinds, aocount, naos);
            auto sub_dens_mat_b = submat::getSubDensityMatrix(densityMatrix, 0, "BETA", aoinds, aocount, naos);

            dengridgen::generateDensityForGGA(rho, rhograd, sigma, npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                              sub_dens_mat_a, sub_dens_mat_b, timer);

            timer.stop("Density matrix slicing");
        }

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_gga(npoints, rho, sigma, exc, vrho, vsigma);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenVxcFockForGGA(rho, sigma, exc, vrho, vsigma, npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        if (closedshell)
        {
            auto partial_mat_Vxc = _integratePartialVxcFockForGGA(
                npoints, local_weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rhograd, vrho, vsigma, timer);

            timer.start("Vxc matrix dist.");

            submat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc, aoinds, aocount, naos);

            timer.stop("Vxc matrix dist.");
        }
        else
        {
            auto partial_mat_Vxc_ab = _integratePartialVxcFockForGGAOpenShell(
                npoints, local_weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rhograd, vrho, vsigma, timer);

            timer.start("Vxc matrix dist.");

            submat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc_ab, aoinds, aocount, naos);

            timer.stop("Vxc matrix dist.");
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
CXCNewIntegrator::_integrateVxcFockForMGGA(const CMolecule&        molecule,
                                           const CMolecularBasis&  basis,
                                           const CAODensityMatrix& densityMatrix,
                                           const CMolecularGrid&   molecularGrid,
                                           const CXCNewFunctional& xcFunctional,
                                           const std::string&      flag) const
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

    bool closedshell = (fstr::upcase(flag) == std::string("CLOSEDSHELL"));

    CAOKohnShamMatrix mat_Vxc(densityMatrix.getNumberOfRows(0), densityMatrix.getNumberOfColumns(0), closedshell);

    mat_Vxc.zero();

    // memory blocks for GTOs on grid points
    // TODO implement Laplacian dependence

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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> rhograd_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> sigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> lapl_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> tau_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> exc_data(1 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vsigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vlapl_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vtau_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma = sigma_data.data();
    auto lapl = lapl_data.data();
    auto tau = tau_data.data();

    auto exc = exc_data.data();
    auto vrho = vrho_data.data();
    auto vsigma = vsigma_data.data();
    auto vlapl = vlapl_data.data();
    auto vtau = vtau_data.data();

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

        if (aocount == 0) continue;

        // generate sub density matrix and density grid

        if (closedshell)
        {
            timer.start("Density matrix slicing");

            auto sub_dens_mat = submat::getSubDensityMatrix(densityMatrix, 0, "ALPHA", aoinds, aocount, naos);

            dengridgen::generateDensityForMGGA(rho, rhograd, sigma, lapl, tau, npoints,
                                               mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                               sub_dens_mat, timer);

            timer.stop("Density matrix slicing");
        }
        else
        {
            timer.start("Density matrix slicing");

            auto sub_dens_mat_a = submat::getSubDensityMatrix(densityMatrix, 0, "ALPHA", aoinds, aocount, naos);
            auto sub_dens_mat_b = submat::getSubDensityMatrix(densityMatrix, 0, "BETA", aoinds, aocount, naos);

            dengridgen::generateDensityForMGGA(rho, rhograd, sigma, lapl, tau, npoints,
                                               mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                               sub_dens_mat_a, sub_dens_mat_b, timer);

            timer.stop("Density matrix slicing");
        }

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_mgga(npoints, rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenExcVxcFockForMGGA(rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau,
                                            npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        if (closedshell)
        {
            auto partial_mat_Vxc = _integratePartialVxcFockForMGGA(
                npoints, local_weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                rhograd, vrho, vsigma, vlapl, vtau, timer);

            timer.start("Vxc matrix dist.");

            submat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc, aoinds, aocount, naos);

            timer.stop("Vxc matrix dist.");
        }
        else
        {
            auto partial_mat_Vxc_ab = _integratePartialVxcFockForMGGAOpenShell(
                npoints, local_weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                rhograd, vrho, vsigma, vlapl, vtau, timer);

            timer.start("Vxc matrix dist.");

            submat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc_ab, aoinds, aocount, naos);

            timer.stop("Vxc matrix dist.");
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
                                          const CXCNewFunctional& xcFunctional) const
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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rhow_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> v2rho2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();

    auto rhow = rhow_data.data();

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

        if (aocount == 0) continue;

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForLDA(rho, npoints, mat_chi, sub_dens_mat, timer);

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_fxc_for_lda(npoints, rho, v2rho2);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenFxcFockForLDA(rho, v2rho2, npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // go through rhow density matrices

        for (int32_t idensity = 0; idensity < rwDensityMatrix.getNumberOfDensityMatrices(); idensity++)
        {
            // generate sub density matrix

            timer.start("Density matrix slicing");

            auto sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, idensity, "ALPHA", aoinds, aocount, naos);

            timer.stop("Density matrix slicing");

            // generate density grid

            dengridgen::generateDensityForLDA(rhow, npoints, mat_chi, sub_dens_mat, timer);

            // compute partial contribution to Fxc matrix

            auto partial_mat_Fxc = _integratePartialFxcFockForLDA(npoints, local_weights, mat_chi, rhow, v2rho2, timer);

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
                                          const CXCNewFunctional& xcFunctional) const
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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rhow_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rhograd_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rhowgrad_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> sigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> vsigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> v2rho2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> v2rhosigma_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> v2sigma2_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();

    auto rhow = rhow_data.data();

    auto rhograd = rhograd_data.data();

    auto rhowgrad = rhowgrad_data.data();

    auto sigma = sigma_data.data();

    auto vrho = vrho_data.data();

    auto vsigma = vsigma_data.data();

    auto v2rho2 = v2rho2_data.data();

    auto v2rhosigma = v2rhosigma_data.data();

    auto v2sigma2 = v2sigma2_data.data();

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

        if (aocount == 0) continue;

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForGGA(rho, rhograd, sigma, npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                          sub_dens_mat, timer);

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

        xcFunctional.compute_fxc_for_gga(npoints, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenFxcFockForGGA(rho, sigma, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2,
                                        npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // go through rhow density matrices

        for (int32_t idensity = 0; idensity < rwDensityMatrix.getNumberOfDensityMatrices(); idensity++)
        {
            // generate sub density matrix

            timer.start("Density matrix slicing");

            auto sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, idensity, "ALPHA", aoinds, aocount, naos);

            timer.stop("Density matrix slicing");

            // generate density grid

            dengridgen::generateDensityForGGA(rhow, rhowgrad, nullptr, npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                              sub_dens_mat, timer);

            // compute partial contribution to Fxc matrix

            auto partial_mat_Fxc = _integratePartialFxcFockForGGA(npoints, local_weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                                                  rhow, rhograd, rhowgrad, vsigma, v2rho2, v2rhosigma, v2sigma2,

                                                                  timer);

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
CXCNewIntegrator::_integrateFxcFockForMGGA(CAOFockMatrix&         aoFockMatrix,
                                           const CMolecule&        molecule,
                                           const CMolecularBasis&  basis,
                                           const CAODensityMatrix& rwDensityMatrix,
                                           const CAODensityMatrix& gsDensityMatrix,
                                           const CMolecularGrid&   molecularGrid,
                                           const CXCNewFunctional& xcFunctional) const
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
    // TODO implement Laplacian dependence

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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    // ground-state
    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> lapl_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> tau_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> rhograd_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> sigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    // perturbed 
    CMemBlock<double> rhowgrad_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> rhow_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> tauw_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> laplw_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    // First-order
    CMemBlock<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vsigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vlapl_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vtau_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    // Second-order 
    CMemBlock<double> v2rho2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2lapl2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2tau2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2v2rholapl_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2v2rhotau_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2lapltau_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2rhosigma_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2sigmalapl_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2sigmatau_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2sigma2_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();
    
    // Ground-state 
    auto rho = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma = sigma_data.data();
    auto lapl = lapl_data.data();
    auto tau = tau_data.data();

    // Perturbed
    auto rhow = rhow_data.data();
    auto rhowgrad = rhowgrad_data.data();
    auto tauw = tauw_data.data();
    auto laplw = laplw_data.data();

    // First-order 
    auto vrho = vrho_data.data();
    auto vsigma = vsigma_data.data();
    auto vlapl = vlapl_data.data();
    auto vtau = vtau_data.data();

    // Second-order
    auto v2rho2  = v2rho2_data.data();
    auto v2lapl2 = v2lapl2_data.data();
    auto v2tau2 =  v2tau2_data.data();
    auto v2rholapl = v2v2rholapl_data.data();
    auto v2rhotau = v2v2rhotau_data.data();
    auto v2lapltau = v2lapltau_data.data();
    auto v2rhosigma = v2rhosigma_data.data();
    auto v2sigmalapl = v2sigmalapl_data.data();
    auto v2sigmatau = v2sigmatau_data.data();
    auto v2sigma2 = v2sigma2_data.data();

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

        if (aocount == 0) continue;

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForMGGA(rho, rhograd, sigma, lapl, tau, npoints,
                                           mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                           sub_dens_mat, timer);

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_mgga(npoints, rho, sigma, lapl, tau, vrho,vsigma,vlapl,vtau);

        xcFunctional.compute_fxc_for_mgga(npoints, rho, sigma, lapl, tau, v2rho2,
                                          v2rhosigma,v2rholapl,v2rhotau,v2sigma2,
                                          v2sigmalapl,v2sigmatau,v2lapl2,v2lapltau,v2tau2);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenFxcFockForMGGA(rho, sigma, lapl, tau, v2rho2,
                                         v2rhosigma,v2rholapl,v2rhotau,v2sigma2,
                                         v2sigmalapl,v2sigmatau,v2lapl2,v2lapltau,v2tau2,
                                         npoints, _screeningThresholdForDensityValues);
            
        timer.stop("Density screening");

        // go through rhow density matrices

        for (int32_t idensity = 0; idensity < rwDensityMatrix.getNumberOfDensityMatrices(); idensity++)
        {
            // generate sub density matrix

            timer.start("Density matrix slicing");

            auto sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, idensity, "ALPHA", aoinds, aocount, naos);

            timer.stop("Density matrix slicing");

            // generate density grid

            dengridgen::generateDensityForMGGA(rhow, rhowgrad, sigma, laplw, tauw, npoints,
                                                mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                sub_dens_mat, timer);

            // compute partial contribution to Fxc matrix

            auto partial_mat_Fxc = _integratePartialFxcFockForMGGA(npoints, local_weights, 
                                                                   mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                                   rhow, rhograd, rhowgrad, tauw, laplw,
                                                                   vrho, vsigma, vlapl, vtau, v2rho2,
                                                                   v2lapl2, v2tau2, v2rholapl, v2rhotau,
                                                                   v2lapltau, v2rhosigma, v2sigmalapl, v2sigmatau, v2sigma2, 
                                                                   timer);

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
                                          const CXCNewFunctional& xcFunctional,
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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> v2rho2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rho3_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();

    auto v2rho2 = v2rho2_data.data();
    auto v3rho3 = v3rho3_data.data();

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

        if (aocount == 0) continue;

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForLDA(rho, npoints, mat_chi, gs_sub_dens_mat, timer);

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto rw_sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, aoinds, aocount);

        auto rw2_sub_dens_mat = submat::getSubDensityMatrix(rw2DensityMatrix, aoinds, aocount);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

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

        xcFunctional.compute_fxc_for_lda(npoints, rho, v2rho2);

        xcFunctional.compute_kxc_for_lda(npoints, rho, v3rho3);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenKxcFockForLDA(rho, v2rho2, v3rho3, npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // go through density matrices

        for (int32_t idensity = 0; idensity < numdens_rw2; idensity++)
        {
            // compute partial contribution to Kxc matrix

            auto partial_mat_Kxc = _integratePartialKxcFockForLDA(npoints, local_weights, mat_chi, v2rho2, v3rho3,
                                                                  rwdengridquad, rw2dengrid, idensity, timer);

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
                                          const CXCNewFunctional& xcFunctional,
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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> rhograd_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> sigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vsigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> v2rho2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2rhosigma_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2sigma2_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> v3rho3_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rho2sigma_data(9 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rhosigma2_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3sigma3_data(10 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma = sigma_data.data();

    auto vrho = vrho_data.data();
    auto vsigma = vsigma_data.data();

    auto v2rho2 = v2rho2_data.data();
    auto v2rhosigma = v2rhosigma_data.data();
    auto v2sigma2 = v2sigma2_data.data();

    auto v3rho3 = v3rho3_data.data();
    auto v3rho2sigma = v3rho2sigma_data.data();
    auto v3rhosigma2 = v3rhosigma2_data.data();
    auto v3sigma3 = v3sigma3_data.data();

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

        if (aocount == 0) continue;

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForGGA(rho, rhograd, sigma, npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                          gs_sub_dens_mat, timer);

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto rw_sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, aoinds, aocount);

        auto rw2_sub_dens_mat = submat::getSubDensityMatrix(rw2DensityMatrix, aoinds, aocount);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

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

        xcFunctional.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

        xcFunctional.compute_fxc_for_gga(npoints, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

        xcFunctional.compute_kxc_for_gga(npoints, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenKxcFockForGGA(rho, sigma, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2,
                                        v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,
                                        npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // go through density matrices

        for (int32_t idensity = 0; idensity < numdens_rw2; idensity++)
        {
            // compute partial contribution to Kxc matrix

            auto partial_mat_Kxc = _integratePartialKxcFockForGGA(npoints, local_weights, mat_chi,
                                                                  mat_chi_x, mat_chi_y, mat_chi_z,
                                                                  rhograd, vsigma, v2rho2, v2rhosigma, v2sigma2,
                                                                  v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,
                                                                  rwdengridquad, rw2dengrid, idensity, timer);

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
CXCNewIntegrator::_integrateKxcFockForMGGA(CAOFockMatrix&          aoFockMatrix,
                                           const CMolecule&        molecule,
                                           const CMolecularBasis&  basis,
                                           const CAODensityMatrix& rwDensityMatrix,
                                           const CAODensityMatrix& rw2DensityMatrix,
                                           const CAODensityMatrix& gsDensityMatrix,
                                           const CMolecularGrid&   molecularGrid,
                                           const CXCNewFunctional& xcFunctional,
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
    // TODO implement Laplacian dependence

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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    // ground-state
    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> lapl_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> tau_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> rhograd_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> sigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    // First-order
    CMemBlock<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vsigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vlapl_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vtau_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    // Second-order
    CMemBlock<double> v2rho2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2lapl2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2tau2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2v2rholapl_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2v2rhotau_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2lapltau_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2rhosigma_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2sigmalapl_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2sigmatau_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2sigma2_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    // Third-order
    CMemBlock<double> v3rho3_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rho2sigma_data(9 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rho2lapl_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rho2tau_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rhosigma2_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rhosigmalapl_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rhosigmatau_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rholapl2_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rholapltau_data(8 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rhotau2_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3sigma3_data(10 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3sigma2lapl_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3sigma2tau_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3sigmalapl2_data(9 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3sigmalapltau_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3sigmatau2_data(9 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3lapl3_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3lapl2tau_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3lapltau2_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3tau3_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    // Ground-state
    auto rho = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma = sigma_data.data();
    auto lapl = lapl_data.data();
    auto tau = tau_data.data();

    // First-order
    auto vrho = vrho_data.data();
    auto vsigma = vsigma_data.data();
    auto vlapl = vlapl_data.data();
    auto vtau = vtau_data.data();

    // Second-order
    auto  v2rho2  = v2rho2_data.data();
    auto  v2lapl2 = v2lapl2_data.data();
    auto  v2tau2 =  v2tau2_data.data();
    auto  v2rholapl = v2v2rholapl_data.data();
    auto  v2rhotau = v2v2rhotau_data.data();
    auto  v2lapltau = v2lapltau_data.data();
    auto  v2rhosigma = v2rhosigma_data.data();
    auto  v2sigmalapl = v2sigmalapl_data.data();
    auto  v2sigmatau = v2sigmatau_data.data();
    auto  v2sigma2 = v2sigma2_data.data();

    // Third-order
    auto v3rho3 = v3rho3_data.data();
    auto v3rho2sigma = v3rho2sigma_data.data();
    auto v3rho2lapl = v3rho2lapl_data.data();
    auto v3rho2tau = v3rho2tau_data.data();
    auto v3rhosigma2 = v3rhosigma2_data.data();
    auto v3rhosigmalapl = v3rhosigmalapl_data.data();
    auto v3rhosigmatau = v3rhosigmatau_data.data();
    auto v3rholapl2 = v3rholapl2_data.data();
    auto v3rholapltau = v3rholapltau_data.data();
    auto v3rhotau2 = v3rhotau2_data.data();
    auto v3sigma3 = v3sigma3_data.data();
    auto v3sigma2lapl = v3sigma2lapl_data.data();
    auto v3sigma2tau = v3sigma2tau_data.data();
    auto v3sigmalapl2 = v3sigmalapl2_data.data();
    auto v3sigmalapltau = v3sigmalapltau_data.data();
    auto v3sigmatau2 = v3sigmatau2_data.data();
    auto v3lapl3 = v3lapl3_data.data();
    auto v3lapl2tau = v3lapl2tau_data.data();
    auto v3lapltau2 = v3lapltau2_data.data();
    auto v3tau3 = v3tau3_data.data();

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

        if (aocount == 0) continue;

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForMGGA(rho, rhograd, sigma, lapl, tau, npoints,
                                           mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                           sub_dens_mat, timer);

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto rw_sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, aoinds, aocount);

        auto rw2_sub_dens_mat = submat::getSubDensityMatrix(rw2DensityMatrix, aoinds, aocount);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto rwdengrid = dengridgen::generateDensityGridForMGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                                rw_sub_dens_mat, xcfuntype, timer);

        auto rw2dengrid = dengridgen::generateDensityGridForMGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                                 rw2_sub_dens_mat, xcfuntype, timer);

        // compute perturbed density

        timer.start("Density grid quad");

        auto numdens_rw2 = rw2DensityMatrix.getNumberOfDensityMatrices();

        CDensityGridQuad rwdengridquad(npoints, numdens_rw2, xcfuntype, dengrid::ab);

        rwdengridquad.DensityProd(rwdengrid, xcfuntype, numdens_rw2, quadMode);

        timer.stop("Density grid quad");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_mgga(npoints, rho, sigma, lapl, tau, vrho,vsigma,vlapl,vtau);

        xcFunctional.compute_fxc_for_mgga(npoints, rho, sigma, lapl, tau, v2rho2,
                                          v2rhosigma,v2rholapl,v2rhotau,v2sigma2,
                                          v2sigmalapl,v2sigmatau,v2lapl2,v2lapltau,v2tau2);

        xcFunctional.compute_kxc_for_mgga(npoints, rho, sigma, lapl, tau,
                                          v3rho3,v3rho2sigma,v3rho2lapl,
                                          v3rho2tau,v3rhosigma2,v3rhosigmalapl,
                                          v3rhosigmatau,v3rholapl2,v3rholapltau,
                                          v3rhotau2,v3sigma3,v3sigma2lapl,
                                          v3sigma2tau,v3sigmalapl2,v3sigmalapltau,
                                          v3sigmatau2,v3lapl3,v3lapl2tau,
                                          v3lapltau2,v3tau3);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenFxcFockForMGGA(rho, sigma, lapl, tau, v2rho2,
                                         v2rhosigma,v2rholapl,v2rhotau,v2sigma2,
                                         v2sigmalapl,v2sigmatau,v2lapl2,v2lapltau,v2tau2,
                                         npoints, _screeningThresholdForDensityValues);

        gridscreen::screenKxcFockForMGGA(rho, sigma, lapl, tau,v3rho3,v3rho2sigma,v3rho2lapl,v3rho2tau,v3rhosigma2,
                                         v3rhosigmalapl,v3rhosigmatau,v3rholapl2,v3rholapltau,v3rhotau2,v3sigma3,
                                         v3sigma2lapl,v3sigma2tau,v3sigmalapl2,v3sigmalapltau,v3sigmatau2,v3lapl3,
                                         v3lapl2tau,v3lapltau2,v3tau3,npoints,_screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // go through density matrices

        for (int32_t idensity = 0; idensity < numdens_rw2; idensity++)
        {
            // compute partial contribution to Kxc matrix

            auto partial_mat_Kxc = _integratePartialKxcFockForMGGA(npoints,local_weights,mat_chi,mat_chi_x,mat_chi_y,
                                                                   mat_chi_z,rhograd,vsigma,
                                                                   v2rho2,v2lapl2,v2tau2,v2rholapl,v2rhotau,
                                                                   v2lapltau,v2rhosigma,v2sigmalapl,v2sigmatau,v2sigma2,
                                                                   v3rho3,v3rho2sigma,v3rho2lapl,v3rho2tau,v3rhosigma2,
                                                                   v3rhosigmalapl,v3rhosigmatau,v3rholapl2,v3rholapltau,v3rhotau2,
                                                                   v3sigma3,v3sigma2lapl,v3sigma2tau,v3sigmalapl2,v3sigmalapltau,
                                                                   v3sigmatau2,v3lapl3,v3lapl2tau,v3lapltau2,v3tau3,
                                                                   rwdengridquad,rw2dengrid,idensity,timer);

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
CXCNewIntegrator::_integrateLxcFockForLDA(CAOFockMatrix&          aoFockMatrix,
                                          const CMolecule&        molecule,
                                          const CMolecularBasis&  basis,
                                          const CAODensityMatrix& rwDensityMatrix,
                                          const CAODensityMatrix& rw2DensityMatrix,
                                          const CAODensityMatrix& rw3DensityMatrix,
                                          const CAODensityMatrix& gsDensityMatrix,
                                          const CMolecularGrid&   molecularGrid,
                                          const CXCNewFunctional& xcFunctional,
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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> v2rho2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rho3_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rho4_data(5 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();

    auto v2rho2 = v2rho2_data.data();
    auto v3rho3 = v3rho3_data.data();
    auto v4rho4 = v4rho4_data.data();

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

        if (aocount == 0) continue;

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForLDA(rho, npoints, mat_chi, gs_sub_dens_mat, timer);

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto rw_sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, aoinds, aocount);

        auto rw2_sub_dens_mat = submat::getSubDensityMatrix(rw2DensityMatrix, aoinds, aocount);

        auto rw3_sub_dens_mat = submat::getSubDensityMatrix(rw3DensityMatrix, aoinds, aocount);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto rwdengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, rw_sub_dens_mat, xcfuntype, timer);

        auto rw2dengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, rw2_sub_dens_mat, xcfuntype, timer);

        auto rw3dengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, rw3_sub_dens_mat, xcfuntype, timer);

        // compute perturbed density

        timer.start("Density grid cube");

        auto numdens_rw3 = rw3DensityMatrix.getNumberOfDensityMatrices();

        CDensityGridCubic rwdengridcube(npoints, numdens_rw3, xcfuntype, dengrid::ab);

        rwdengridcube.DensityProd(rwdengrid, rw2dengrid, xcfuntype, numdens_rw3, quadMode);

        timer.stop("Density grid cubic");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_fxc_for_lda(npoints, rho, v2rho2);

        xcFunctional.compute_kxc_for_lda(npoints, rho, v3rho3);

        xcFunctional.compute_lxc_for_lda(npoints, rho, v4rho4);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenLxcFockForLDA(rho, v2rho2, v3rho3, v4rho4, npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // go through density matrices

        for (int32_t idensity = 0; idensity < numdens_rw3; idensity++)
        {
            // compute partial contribution to Kxc matrix

            auto partial_mat_Lxc = _integratePartialLxcFockForLDA(npoints, local_weights, mat_chi, v2rho2, v3rho3,v4rho4,
                                                                  rwdengridcube, rw3dengrid, idensity, timer);

            // distribute partial Kxc to full Fock matrix

            timer.start("Lxc matrix dist.");

            submat::distributeSubMatrixToFock(aoFockMatrix, idensity, partial_mat_Lxc, aoinds, aocount, naos);

            timer.stop("Lxc matrix dist.");
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
CXCNewIntegrator::_integrateLxcFockForGGA(CAOFockMatrix&          aoFockMatrix,
                                          const CMolecule&        molecule,
                                          const CMolecularBasis&  basis,
                                          const CAODensityMatrix& rwDensityMatrix,
                                          const CAODensityMatrix& rw2DensityMatrix,
                                          const CAODensityMatrix& rw3DensityMatrix,
                                          const CAODensityMatrix& gsDensityMatrix,
                                          const CMolecularGrid&   molecularGrid,
                                          const CXCNewFunctional& xcFunctional,
                                          const std::string&      cubeMode) const
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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> rhograd_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> sigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vsigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> v2rho2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2rhosigma_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2sigma2_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> v3rho3_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rho2sigma_data(9 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rhosigma2_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3sigma3_data(10 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> v4rho4_data(5 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rho3sigma_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rho2sigma2_data(18 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rhosigma3_data(20 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4sigma4_data(15 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma = sigma_data.data();

    auto vrho = vrho_data.data();
    auto vsigma = vsigma_data.data();

    auto v2rho2 = v2rho2_data.data();
    auto v2rhosigma = v2rhosigma_data.data();
    auto v2sigma2 = v2sigma2_data.data();

    auto v3rho3 = v3rho3_data.data();
    auto v3rho2sigma = v3rho2sigma_data.data();
    auto v3rhosigma2 = v3rhosigma2_data.data();
    auto v3sigma3 = v3sigma3_data.data();

    auto v4rho4 = v4rho4_data.data();
    auto v4rho3sigma = v4rho3sigma_data.data();
    auto v4rho2sigma2 = v4rho2sigma2_data.data();
    auto v4rhosigma3 = v4rhosigma3_data.data();
    auto v4sigma4 = v4sigma4_data.data();

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

        if (aocount == 0) continue;

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForGGA(rho, rhograd, sigma, npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                          gs_sub_dens_mat, timer);

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto rw_sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, aoinds, aocount);

        auto rw2_sub_dens_mat = submat::getSubDensityMatrix(rw2DensityMatrix, aoinds, aocount);

        auto rw3_sub_dens_mat = submat::getSubDensityMatrix(rw3DensityMatrix, aoinds, aocount);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto rwdengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                               rw_sub_dens_mat, xcfuntype, timer);

        auto rw2dengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                                rw2_sub_dens_mat, xcfuntype, timer);

        auto rw3dengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                                rw3_sub_dens_mat, xcfuntype, timer);
        // compute perturbed density

        timer.start("Density grid cubic");

        auto numdens_rw3 = rw3DensityMatrix.getNumberOfDensityMatrices();

        CDensityGridCubic rwdengridcubic(npoints, numdens_rw3, xcfuntype, dengrid::ab);

        rwdengridcubic.DensityProd(rwdengrid, rw2dengrid, xcfuntype, numdens_rw3, cubeMode);

        timer.stop("Density grid cubic");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

        xcFunctional.compute_fxc_for_gga(npoints, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

        xcFunctional.compute_kxc_for_gga(npoints, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);

        xcFunctional.compute_lxc_for_gga(npoints, rho, sigma, v4rho4, v4rho3sigma, v4rho2sigma2,v4rhosigma3, v4sigma4);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenLxcFockForGGA(rho, sigma, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2,
                                        v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,
                                        v4rho4, v4rho3sigma, v4rho2sigma2,v4rhosigma3, v4sigma4,
                                        npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // go through density matrices

        for (int32_t idensity = 0; idensity < numdens_rw3; idensity++)
        {
            // compute partial contribution to Lxc matrix

            auto partial_mat_Kxc = _integratePartialLxcFockForGGA(npoints, local_weights, mat_chi,
                                                                  mat_chi_x, mat_chi_y, mat_chi_z,
                                                                  rhograd, vsigma, v2rho2, v2rhosigma, v2sigma2,
                                                                  v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,
                                                                  v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4,
                                                                  rwdengridcubic, rw3dengrid, idensity, timer);

            // distribute partial Lxc to full Fock matrix

            timer.start("Lxc matrix dist.");

            submat::distributeSubMatrixToFock(aoFockMatrix, idensity, partial_mat_Kxc, aoinds, aocount, naos);

            timer.stop("Lxc matrix dist.");
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
CXCNewIntegrator::_integrateKxcLxcFockForLDA(CAOFockMatrix&          aoFockMatrix,
                                             const CMolecule&        molecule,
                                             const CMolecularBasis&  basis,
                                             const CAODensityMatrix& rwDensityMatrix,
                                             const CAODensityMatrix& rw2DensityMatrix,
                                             const CAODensityMatrix& rw3DensityMatrix,
                                             const CAODensityMatrix& gsDensityMatrix,
                                             const CMolecularGrid&   molecularGrid,
                                             const CXCNewFunctional& xcFunctional,
                                             const std::string&      cubeMode) const
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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> v2rho2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rho3_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rho4_data(5 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();

    auto v2rho2 = v2rho2_data.data();
    auto v3rho3 = v3rho3_data.data();
    auto v4rho4 = v4rho4_data.data();

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

        if (aocount == 0) continue;

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForLDA(rho, npoints, mat_chi, gs_sub_dens_mat, timer);

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto rw_sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, aoinds, aocount);

        auto rw2_sub_dens_mat = submat::getSubDensityMatrix(rw2DensityMatrix, aoinds, aocount);

        auto rw3_sub_dens_mat = submat::getSubDensityMatrix(rw3DensityMatrix, aoinds, aocount);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto rwdengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, rw_sub_dens_mat, xcfuntype, timer);

        auto rw2dengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, rw2_sub_dens_mat, xcfuntype, timer);

        auto rw3dengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, rw3_sub_dens_mat, xcfuntype, timer);

        // compute perturbed density

        timer.start("Density grid cube");

        auto numdens_rw3 = rw3DensityMatrix.getNumberOfDensityMatrices();

        auto numdens_rw2 = rw2DensityMatrix.getNumberOfDensityMatrices();

        CDensityGridCubic rwdengridcube(npoints, numdens_rw2 + numdens_rw3, xcfuntype, dengrid::ab);

        rwdengridcube.DensityProd(rwdengrid, rw2dengrid, xcfuntype, (numdens_rw2 + numdens_rw3), cubeMode);

        timer.stop("Density grid cubic");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_fxc_for_lda(npoints, rho, v2rho2);

        xcFunctional.compute_kxc_for_lda(npoints, rho, v3rho3);

        xcFunctional.compute_lxc_for_lda(npoints, rho, v4rho4);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenLxcFockForLDA(rho, v2rho2, v3rho3, v4rho4, npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // go through density matrices

        for  (int32_t idensity = 0; idensity < numdens_rw2; idensity++)
        {
            // compute partial contribution to Kxc matrix

            auto partial_mat_Kxc = _integratePartialKxcFockForLDA2(npoints, local_weights, mat_chi, v2rho2, v3rho3,
                                                                  rwdengridcube, rw2dengrid, idensity, timer);

            // distribute partial Kxc to full Fock matrix

            timer.start("Kxc matrix dist.");

            submat::distributeSubMatrixToFock(aoFockMatrix, idensity, partial_mat_Kxc, aoinds, aocount, naos);

            timer.stop("Kxc matrix dist.");
        }

        for (int32_t idensity = 0; idensity < numdens_rw3; idensity++)
        {
            // compute partial contribution to Kxc matrix

            auto partial_mat_Lxc = _integratePartialLxcFockForLDA(npoints, local_weights, mat_chi, v2rho2, v3rho3,v4rho4,
                                                                  rwdengridcube, rw3dengrid, idensity, timer);

            // distribute partial Kxc to full Fock matrix

            timer.start("Lxc matrix dist.");

            submat::distributeSubMatrixToFock(aoFockMatrix, (idensity + numdens_rw2), partial_mat_Lxc, aoinds, aocount, naos);

            timer.stop("Lxc matrix dist.");
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
CXCNewIntegrator::_integrateKxcLxcFockForGGA(CAOFockMatrix&          aoFockMatrix,
                                             const CMolecule&        molecule,
                                             const CMolecularBasis&  basis,
                                             const CAODensityMatrix& rwDensityMatrix,
                                             const CAODensityMatrix& rw2DensityMatrix,
                                             const CAODensityMatrix& rw3DensityMatrix,
                                             const CAODensityMatrix& gsDensityMatrix,
                                             const CMolecularGrid&   molecularGrid,
                                             const CXCNewFunctional& xcFunctional,
                                             const std::string&      cubeMode) const
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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> rhograd_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> sigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vsigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> v2rho2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2rhosigma_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2sigma2_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> v3rho3_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rho2sigma_data(9 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rhosigma2_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3sigma3_data(10 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> v4rho4_data(5 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rho3sigma_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rho2sigma2_data(18 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rhosigma3_data(20 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4sigma4_data(15 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma = sigma_data.data();

    auto vrho = vrho_data.data();
    auto vsigma = vsigma_data.data();

    auto v2rho2 = v2rho2_data.data();
    auto v2rhosigma = v2rhosigma_data.data();
    auto v2sigma2 = v2sigma2_data.data();

    auto v3rho3 = v3rho3_data.data();
    auto v3rho2sigma = v3rho2sigma_data.data();
    auto v3rhosigma2 = v3rhosigma2_data.data();
    auto v3sigma3 = v3sigma3_data.data();

    auto v4rho4 = v4rho4_data.data();
    auto v4rho3sigma = v4rho3sigma_data.data();
    auto v4rho2sigma2 = v4rho2sigma2_data.data();
    auto v4rhosigma3 = v4rhosigma3_data.data();
    auto v4sigma4 = v4sigma4_data.data();

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

        if (aocount == 0) continue;

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForGGA(rho, rhograd, sigma, npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                          gs_sub_dens_mat, timer);

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto rw_sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, aoinds, aocount);

        auto rw2_sub_dens_mat = submat::getSubDensityMatrix(rw2DensityMatrix, aoinds, aocount);

        auto rw3_sub_dens_mat = submat::getSubDensityMatrix(rw3DensityMatrix, aoinds, aocount);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto rwdengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                               rw_sub_dens_mat, xcfuntype, timer);

        auto rw2dengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                                rw2_sub_dens_mat, xcfuntype, timer);

        auto rw3dengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                                rw3_sub_dens_mat, xcfuntype, timer);
        // compute perturbed density

        timer.start("Density grid cubic");

        auto numdens_rw3 = rw3DensityMatrix.getNumberOfDensityMatrices();

        auto numdens_rw2 = rw2DensityMatrix.getNumberOfDensityMatrices();

        CDensityGridCubic rwdengridcube(npoints, (numdens_rw2 + numdens_rw3), xcfuntype, dengrid::ab);

        rwdengridcube.DensityProd(rwdengrid, rw2dengrid, xcfuntype, (numdens_rw2 + numdens_rw3), cubeMode);

        timer.stop("Density grid cubic");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

        xcFunctional.compute_fxc_for_gga(npoints, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

        xcFunctional.compute_kxc_for_gga(npoints, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);

        xcFunctional.compute_lxc_for_gga(npoints, rho, sigma, v4rho4, v4rho3sigma, v4rho2sigma2,v4rhosigma3, v4sigma4);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenLxcFockForGGA(rho, sigma, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2,
                                        v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,
                                        v4rho4, v4rho3sigma, v4rho2sigma2,v4rhosigma3, v4sigma4,
                                        npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // go through density matrices

        for (int32_t idensity = 0; idensity < numdens_rw2; idensity++)
        {
            // compute partial contribution to Kxc matrix

            auto partial_mat_Kxc = _integratePartialKxcFockForGGA2(npoints, local_weights, mat_chi,
                                                                   mat_chi_x, mat_chi_y, mat_chi_z,
                                                                   rhograd, vsigma, v2rho2, v2rhosigma, v2sigma2,
                                                                   v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,
                                                                   rwdengridcube, rw2dengrid, idensity, timer);

            // distribute partial Kxc to full Fock matrix

            timer.start("Kxc matrix dist.");

            submat::distributeSubMatrixToFock(aoFockMatrix, idensity, partial_mat_Kxc, aoinds, aocount, naos);

            timer.stop("Kxc matrix dist.");
        }

        for (int32_t idensity = 0; idensity < numdens_rw3; idensity++)
        {
            // compute partial contribution to Lxc matrix

            auto partial_mat_Lxc = _integratePartialLxcFockForGGA(npoints, local_weights, mat_chi,
                                                                  mat_chi_x, mat_chi_y, mat_chi_z,
                                                                  rhograd, vsigma, v2rho2, v2rhosigma, v2sigma2,
                                                                  v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,
                                                                  v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4,
                                                                  rwdengridcube, rw3dengrid, idensity, timer);

            // distribute partial Lxc to full Fock matrix

            timer.start("Lxc matrix dist.");

            submat::distributeSubMatrixToFock(aoFockMatrix, (idensity + numdens_rw2), partial_mat_Lxc, aoinds, aocount, naos);

            timer.stop("Lxc matrix dist.");
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
CXCNewIntegrator::_integrateKxcLxcFockForMGGA(CAOFockMatrix&      aoFockMatrix,
                                          const CMolecule&        molecule,
                                          const CMolecularBasis&  basis,
                                          const CAODensityMatrix& rwDensityMatrix,
                                          const CAODensityMatrix& rw2DensityMatrix,
                                          const CAODensityMatrix& rw3DensityMatrix,
                                          const CAODensityMatrix& gsDensityMatrix,
                                          const CMolecularGrid&   molecularGrid,
                                          const CXCNewFunctional& xcFunctional,
                                          const std::string&      cubeMode) const
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
    // TODO implement Laplacian dependence

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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    // ground-state
    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> lapl_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> tau_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> rhograd_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> sigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    // First-order
    CMemBlock<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vsigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vlapl_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vtau_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    // Second-order
    CMemBlock<double> v2rho2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2lapl2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2tau2_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2v2rholapl_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2v2rhotau_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2lapltau_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2rhosigma_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2sigmalapl_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2sigmatau_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v2sigma2_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    // Third-order
    CMemBlock<double> v3rho3_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rho2sigma_data(9 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rho2lapl_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rho2tau_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rhosigma2_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rhosigmalapl_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rhosigmatau_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rholapl2_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rholapltau_data(8 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3rhotau2_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3sigma3_data(10 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3sigma2lapl_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3sigma2tau_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3sigmalapl2_data(9 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3sigmalapltau_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3sigmatau2_data(9 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3lapl3_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3lapl2tau_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3lapltau2_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v3tau3_data(4 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    // Fourth-order
    CMemBlock<double> v4rho4_data(5 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rho3sigma_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rho3lapl_data(8 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rho3tau_data(8 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rho2sigma2_data(18 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rho2sigmalapl_data(18 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rho2sigmatau_data(18 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rho2lapl2_data(9 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rho2lapltau_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rho2tau2_data(9 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rhosigma3_data(20 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rhosigma2lapl_data(36 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rhosigma2tau_data(36 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rhosigmalapl2_data(18 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rhosigmalapltau_data(24 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rhosigmatau2_data(36 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rholapl3_data(8 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rholapl2tau_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rholapltau2_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4rhotau3_data(8 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4sigma4_data(15 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4sigma3lapl_data(20 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4sigma3tau_data(30 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4sigma2lapl2_data(18 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4sigma2lapltau_data(24 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4sigma2tau2_data(18 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4sigmalapl3_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4sigmalapl2tau_data(18 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4sigmalapltau2_data(18 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4sigmatau3_data(12 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4lapl4_data(5 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4lapl3tau_data(8 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4lapl2tau2_data(9 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4lapltau3_data(8 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> v4tau4_data(5 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    // Ground-state
    auto rho = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma = sigma_data.data();
    auto lapl = lapl_data.data();
    auto tau = tau_data.data();

    // First-order
    auto vrho = vrho_data.data();
    auto vsigma = vsigma_data.data();
    auto vlapl = vlapl_data.data();
    auto vtau = vtau_data.data();

    // Second-order
    auto  v2rho2  = v2rho2_data.data();
    auto  v2lapl2 = v2lapl2_data.data();
    auto  v2tau2 =  v2tau2_data.data();
    auto  v2rholapl = v2v2rholapl_data.data();
    auto  v2rhotau = v2v2rhotau_data.data();
    auto  v2lapltau = v2lapltau_data.data();
    auto  v2rhosigma = v2rhosigma_data.data();
    auto  v2sigmalapl = v2sigmalapl_data.data();
    auto  v2sigmatau = v2sigmatau_data.data();
    auto  v2sigma2 = v2sigma2_data.data();

    // Third-order
    auto v3rho3 = v3rho3_data.data();
    auto v3rho2sigma = v3rho2sigma_data.data();
    auto v3rho2lapl = v3rho2lapl_data.data();
    auto v3rho2tau = v3rho2tau_data.data();
    auto v3rhosigma2 = v3rhosigma2_data.data();
    auto v3rhosigmalapl = v3rhosigmalapl_data.data();
    auto v3rhosigmatau = v3rhosigmatau_data.data();
    auto v3rholapl2 = v3rholapl2_data.data();
    auto v3rholapltau = v3rholapltau_data.data();
    auto v3rhotau2 = v3rhotau2_data.data();
    auto v3sigma3 = v3sigma3_data.data();
    auto v3sigma2lapl = v3sigma2lapl_data.data();
    auto v3sigma2tau = v3sigma2tau_data.data();
    auto v3sigmalapl2 = v3sigmalapl2_data.data();
    auto v3sigmalapltau = v3sigmalapltau_data.data();
    auto v3sigmatau2 = v3sigmatau2_data.data();
    auto v3lapl3 = v3lapl3_data.data();
    auto v3lapl2tau = v3lapl2tau_data.data();
    auto v3lapltau2 = v3lapltau2_data.data();
    auto v3tau3 = v3tau3_data.data();

    // Fourth-order
    auto v4rho4 = v4rho4_data.data();
    auto v4rho3sigma = v4rho3sigma_data.data();
    auto v4rho3lapl = v4rho3lapl_data.data();
    auto v4rho3tau = v4rho3tau_data.data();
    auto v4rho2sigma2 = v4rho2sigma2_data.data();
    auto v4rho2sigmalapl = v4rho2sigmalapl_data.data();
    auto v4rho2sigmatau = v4rho2sigmatau_data.data();
    auto v4rho2lapl2 = v4rho2lapl2_data.data();
    auto v4rho2lapltau = v4rho2lapltau_data.data();
    auto v4rho2tau2 = v4rho2tau2_data.data();
    auto v4rhosigma3 = v4rhosigma3_data.data();
    auto v4rhosigma2lapl = v4rhosigma2lapl_data.data();
    auto v4rhosigma2tau = v4rhosigma2tau_data.data();
    auto v4rhosigmalapl2 = v4rhosigmalapl2_data.data();
    auto v4rhosigmalapltau = v4rhosigmalapltau_data.data();
    auto v4rhosigmatau2 = v4rhosigmatau2_data.data();
    auto v4rholapl3 = v4rholapl3_data.data();
    auto v4rholapl2tau = v4rholapl2tau_data.data();
    auto v4rholapltau2 = v4rholapltau2_data.data();
    auto v4rhotau3 = v4rhotau3_data.data();
    auto v4sigma4 = v4sigma4_data.data();
    auto v4sigma3lapl = v4sigma3lapl_data.data();
    auto v4sigma3tau = v4sigma3tau_data.data();
    auto v4sigma2lapl2 = v4sigma2lapl2_data.data();
    auto v4sigma2lapltau = v4sigma2lapltau_data.data();
    auto v4sigma2tau2 = v4sigma2tau2_data.data();
    auto v4sigmalapl3 = v4sigmalapl3_data.data();
    auto v4sigmalapl2tau = v4sigmalapl2tau_data.data();
    auto v4sigmalapltau2 = v4sigmalapltau2_data.data();
    auto v4sigmatau3 = v4sigmatau3_data.data();
    auto v4lapl4 = v4lapl4_data.data();
    auto v4lapl3tau = v4lapl3tau_data.data();
    auto v4lapl2tau2 = v4lapl2tau2_data.data();
    auto v4lapltau3 = v4lapltau3_data.data();
    auto v4tau4 = v4tau4_data.data();

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

        if (aocount == 0) continue;

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForMGGA(rho, rhograd, sigma, lapl, tau, npoints,
                                           mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                           sub_dens_mat, timer);

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto rw_sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, aoinds, aocount);

        auto rw2_sub_dens_mat = submat::getSubDensityMatrix(rw2DensityMatrix, aoinds, aocount);

        auto rw3_sub_dens_mat = submat::getSubDensityMatrix(rw3DensityMatrix, aoinds, aocount);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto rwdengrid = dengridgen::generateDensityGridForMGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                                rw_sub_dens_mat, xcfuntype, timer);

        auto rw2dengrid = dengridgen::generateDensityGridForMGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                                rw2_sub_dens_mat, xcfuntype, timer);

        auto rw3dengrid = dengridgen::generateDensityGridForMGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                                 rw3_sub_dens_mat, xcfuntype, timer);

        // compute perturbed density

        timer.start("Density grid quad");

        auto numdens_rw3 = rw3DensityMatrix.getNumberOfDensityMatrices();

        auto numdens_rw2 = rw2DensityMatrix.getNumberOfDensityMatrices();

        CDensityGridCubic rwdengridcube(npoints, (numdens_rw2 + numdens_rw3), xcfuntype, dengrid::ab);

        rwdengridcube.DensityProd(rwdengrid, rw2dengrid, xcfuntype, (numdens_rw2 + numdens_rw3), cubeMode);

        timer.stop("Density grid quad");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_mgga(npoints, rho, sigma, lapl, tau, vrho,vsigma,vlapl,vtau);

        xcFunctional.compute_fxc_for_mgga(npoints, rho, sigma, lapl, tau, v2rho2,
                                          v2rhosigma,v2rholapl,v2rhotau,v2sigma2,
                                          v2sigmalapl,v2sigmatau,v2lapl2,v2lapltau,v2tau2);

        xcFunctional.compute_kxc_for_mgga(npoints, rho, sigma, lapl, tau,
                                          v3rho3,v3rho2sigma,v3rho2lapl,
                                          v3rho2tau,v3rhosigma2,v3rhosigmalapl,
                                          v3rhosigmatau,v3rholapl2,v3rholapltau,
                                          v3rhotau2,v3sigma3,v3sigma2lapl,
                                          v3sigma2tau,v3sigmalapl2,v3sigmalapltau,
                                          v3sigmatau2,v3lapl3,v3lapl2tau,
                                          v3lapltau2,v3tau3);

        xcFunctional.compute_lxc_for_mgga(npoints, rho, sigma, lapl,tau,
                                          v4rho4,v4rho3sigma,v4rho3lapl,v4rho3tau,v4rho2sigma2,
                                          v4rho2sigmalapl,v4rho2sigmatau,v4rho2lapl2,v4rho2lapltau,
                                          v4rho2tau2,v4rhosigma3,v4rhosigma2lapl,v4rhosigma2tau,
                                          v4rhosigmalapl2,v4rhosigmalapltau,v4rhosigmatau2,v4rholapl3,
                                          v4rholapl2tau,v4rholapltau2,v4rhotau3,v4sigma4,v4sigma3lapl,
                                          v4sigma3tau,v4sigma2lapl2,v4sigma2lapltau,v4sigma2tau2,
                                          v4sigmalapl3,v4sigmalapl2tau,v4sigmalapltau2,v4sigmatau3,
                                          v4lapl4,v4lapl3tau,v4lapl2tau2,v4lapltau3,v4tau4);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenFxcFockForMGGA(rho, sigma, lapl, tau, v2rho2,
                                         v2rhosigma,v2rholapl,v2rhotau,v2sigma2,
                                         v2sigmalapl,v2sigmatau,v2lapl2,v2lapltau,v2tau2,
                                         npoints, _screeningThresholdForDensityValues);

        gridscreen::screenKxcFockForMGGA(rho, sigma, lapl, tau,v3rho3,v3rho2sigma,v3rho2lapl,v3rho2tau,v3rhosigma2,
                                         v3rhosigmalapl,v3rhosigmatau,v3rholapl2,v3rholapltau,v3rhotau2,v3sigma3,
                                         v3sigma2lapl,v3sigma2tau,v3sigmalapl2,v3sigmalapltau,v3sigmatau2,v3lapl3,
                                         v3lapl2tau,v3lapltau2,v3tau3,npoints,_screeningThresholdForDensityValues);

        gridscreen::screenLxcFockForMGGA(rho, sigma, lapl, tau,v4rho4,v4rho3sigma,
                                         v4rho3lapl,v4rho3tau,v4rho2sigma2,v4rho2sigmalapl,v4rho2sigmatau,
                                         v4rho2lapl2,v4rho2lapltau,v4rho2tau2,v4rhosigma3,v4rhosigma2lapl,
                                         v4rhosigma2tau,v4rhosigmalapl2,v4rhosigmalapltau,v4rhosigmatau2,
                                         v4rholapl3,v4rholapl2tau,v4rholapltau2,v4rhotau3,v4sigma4,v4sigma3lapl,
                                         v4sigma3tau,v4sigma2lapl2,v4sigma2lapltau,v4sigma2tau2,v4sigmalapl3,
                                         v4sigmalapl2tau,v4sigmalapltau2,v4sigmatau3,v4lapl4,v4lapl3tau,
                                         v4lapl2tau2,v4lapltau3,v4tau4,npoints,_screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // go through density matrices

        for (int32_t idensity = 0; idensity < numdens_rw2; idensity++)
        {
            // compute partial contribution to Kxc matrix

           auto partial_mat_Kxc = _integratePartialKxcFockForMGGA2(npoints,local_weights,mat_chi,mat_chi_x,mat_chi_y,mat_chi_z,
                                                                   rhograd,vsigma,v2rho2,v2lapl2,v2tau2,v2rholapl,v2rhotau,
                                                                   v2lapltau,v2rhosigma,v2sigmalapl,v2sigmatau,v2sigma2,
                                                                   v3rho3,v3rho2sigma,v3rho2lapl,v3rho2tau,v3rhosigma2,
                                                                   v3rhosigmalapl,v3rhosigmatau,v3rholapl2,v3rholapltau,v3rhotau2,
                                                                   v3sigma3,v3sigma2lapl,v3sigma2tau,v3sigmalapl2,v3sigmalapltau,
                                                                   v3sigmatau2,v3lapl3,v3lapl2tau,v3lapltau2,v3tau3,
                                                                   rwdengridcube,rw2dengrid,idensity,timer);

            // distribute partial Kxc to full Fock matrix

            timer.start("Kxc matrix dist.");

            submat::distributeSubMatrixToFock(aoFockMatrix, idensity, partial_mat_Kxc, aoinds, aocount, naos);

            timer.stop("Kxc matrix dist.");
        }

        for (int32_t idensity = 0; idensity < numdens_rw3; idensity++)
        {
            // compute partial contribution to Lxc matrix

            auto partial_mat_Lxc = _integratePartialLxcFockForMGGA(npoints, local_weights, mat_chi,mat_chi_x,mat_chi_y,mat_chi_z,
                                                                   rhograd,vsigma, v2rho2,v2lapl2, v2tau2, v2rholapl, v2rhotau,v2lapltau,
                                                                   v2rhosigma, v2sigmalapl, v2sigmatau,v2sigma2,
                                                                   v3rho3,v3rho2sigma,v3rho2lapl,v3rho2tau,v3rhosigma2,
                                                                   v3rhosigmalapl,v3rhosigmatau,v3rholapl2,v3rholapltau,
                                                                   v3rhotau2,v3sigma3,v3sigma2lapl,v3sigma2tau,v3sigmalapl2,
                                                                   v3sigmalapltau,v3sigmatau2,v3lapl3,v3lapl2tau,v3lapltau2,
                                                                   v3tau3, v4rho4,v4rho3sigma,v4rho3lapl,v4rho3tau,v4rho2sigma2,
                                                                   v4rho2sigmalapl,v4rho2sigmatau,v4rho2lapl2,v4rho2lapltau,
                                                                   v4rho2tau2,v4rhosigma3,v4rhosigma2lapl,v4rhosigma2tau,
                                                                   v4rhosigmalapl2,v4rhosigmalapltau,v4rhosigmatau2,v4rholapl3,
                                                                   v4rholapl2tau,v4rholapltau2,v4rhotau3,v4sigma4,v4sigma3lapl,
                                                                   v4sigma3tau,v4sigma2lapl2,v4sigma2lapltau,v4sigma2tau2,
                                                                   v4sigmalapl3,v4sigmalapl2tau,v4sigmalapltau2,v4sigmatau3,
                                                                   v4lapl4,v4lapl3tau,v4lapl2tau2,v4lapltau3,v4tau4,
                                                                   rwdengridcube, rw3dengrid,idensity,timer);

            // distribute partial Lxc to full Fock matrix

            timer.start("Lxc matrix dist.");

            submat::distributeSubMatrixToFock(aoFockMatrix, (idensity + numdens_rw2), partial_mat_Lxc, aoinds, aocount, naos);

            timer.stop("Lxc matrix dist.");
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
CXCNewIntegrator::_integrateVxcPDFTForLDA(CAOKohnShamMatrix&              aoFockMatrix,
                                          CDense4DTensor&                 moTwoBodyGradient,
                                          const CMolecule&                molecule,
                                          const CMolecularBasis&          basis,
                                          const CAODensityMatrix&         DensityMatrix,
                                          const CDense4DTensor&           TwoBodyDensityMatrix,
                                          const CDenseMatrix&             ActiveMOs,
                                          const CMolecularGrid&           molecularGrid,
                                          const CXCPairDensityFunctional& xcFunctional) const
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // Set up Fock matrix

    aoFockMatrix.zero();

    // memory blocks for GTOs on grid points

    CMemBlock2D<double> gaos(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos); // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> exc_data(1 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();

    auto exc = exc_data.data();
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

        if (aocount == 0) continue;

        // generate sub density matrix and MO coefficients

        timer.start("Density matrix slicing");

        auto sub_dens_mat_a = submat::getSubDensityMatrix(DensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        auto sub_ActiveMOs = submat::getSubMatrixByColumnSlicing(ActiveMOs, aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density and on-top pair density on the grid

        dengridgen::generatePairDensityForLDA(rho, npoints, mat_chi, sub_dens_mat_a, sub_ActiveMOs, TwoBodyDensityMatrix, timer);

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_plda(npoints, rho, exc, vrho);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenVxcFockForPLDA(rho, exc, vrho, npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        auto partial_mat_Vxc = _integratePartialVxcFockForLDA(npoints, local_weights, mat_chi, vrho, timer);

        // TODO (MGD) 2-body gradient

        // distribute partial Vxc to full Kohn-Sham matrix

        timer.start("Vxc matrix dist.");

        submat::distributeSubMatrixToKohnSham(aoFockMatrix, partial_mat_Vxc, aoinds, aocount, naos);

        timer.stop("Vxc matrix dist.");

        // compute partial contribution to XC energy

        timer.start("XC energy");

        for (int32_t g = 0; g < npoints; g++)
        {
            auto rho_total = rho[2 * g + 0];

            nele += local_weights[g] * rho_total;

            xcene += local_weights[g] * exc[g] * rho_total;
        }

        timer.stop("XC energy");
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    aoFockMatrix.setNumberOfElectrons(nele);

    aoFockMatrix.setExchangeCorrelationEnergy(xcene);
}

void
CXCNewIntegrator::_integrateVxcPDFTForGGA(CAOKohnShamMatrix&              aoFockMatrix,
                                          CDense4DTensor&                 moTwoBodyGradient,
                                          const CMolecule&                molecule,
                                          const CMolecularBasis&          basis,
                                          const CAODensityMatrix&         DensityMatrix,
                                          const CDense4DTensor&           TwoBodyDensityMatrix,
                                          const CDenseMatrix&             ActiveMOs,
                                          const CMolecularGrid&           molecularGrid,
                                          const CXCPairDensityFunctional& xcFunctional) const
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // Set up Fock matrix

    aoFockMatrix.zero();

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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> rhograd_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> sigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> exc_data(1 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    CMemBlock<double> vsigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma = sigma_data.data();

    auto exc = exc_data.data();
    auto vrho = vrho_data.data();
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

        if (aocount == 0) continue;

        // generate sub density matrix and MO coefficients

        timer.start("Density matrix slicing");

        auto sub_dens_mat_a = submat::getSubDensityMatrix(DensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        auto sub_ActiveMOs = submat::getSubMatrixByColumnSlicing(ActiveMOs, aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density and on-top pair density on the grid

        dengridgen::generatePairDensityForGGA(rho, rhograd, sigma, npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, sub_dens_mat_a, sub_ActiveMOs, TwoBodyDensityMatrix, timer);

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_pgga(npoints, rho, sigma, exc, vrho, vsigma);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenVxcFockForPGGA(rho, sigma, exc, vrho, vsigma, npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        auto partial_mat_Vxc = _integratePartialVxcFockForGGA(
                npoints, local_weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, rhograd, vrho, vsigma, timer);

        // TODO (MGD) 2-body gradient

        // distribute partial Vxc to full Kohn-Sham matrix

        timer.start("Vxc matrix dist.");

        submat::distributeSubMatrixToKohnSham(aoFockMatrix, partial_mat_Vxc, aoinds, aocount, naos);

        timer.stop("Vxc matrix dist.");

        // compute partial contribution to XC energy

        timer.start("XC energy");

        for (int32_t g = 0; g < npoints; g++)
        {
            auto rho_total = rho[2 * g + 0];

            nele += local_weights[g] * rho_total;

            xcene += local_weights[g] * exc[g] * rho_total;
        }

        timer.stop("XC energy");
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    aoFockMatrix.setNumberOfElectrons(nele);

    aoFockMatrix.setExchangeCorrelationEnergy(xcene);
}

CDenseMatrix
CXCNewIntegrator::_integratePartialVxcFockForLDA(const int32_t          npoints,
                                                 const double*          weights,
                                                 const CDenseMatrix&    gtoValues,
                                                 const double*          vrho,
                                                 CMultiTimer&           timer) const
{
    // GTO values on grid points

    auto chi_val = gtoValues.values();

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix G");

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

            #pragma omp simd aligned(weights, vrho, G_val, chi_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                G_val[nu_offset + g] = weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];
            }
        }
    }

    timer.stop("Vxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix matmul");

    auto mat_Vxc = denblas::multABt(gtoValues, mat_G);

    timer.stop("Vxc matrix matmul");

    return mat_Vxc;
}

std::vector<CDenseMatrix>
CXCNewIntegrator::_integratePartialVxcFockForLDAOpenShell(const int32_t          npoints,
                                                          const double*          weights,
                                                          const CDenseMatrix&    gtoValues,
                                                          const double*          vrho,
                                                          CMultiTimer&           timer) const
{
    // GTO values on grid points

    auto chi_val = gtoValues.values();

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G_a(naos, npoints);

    CDenseMatrix mat_G_b(naos, npoints);

    auto G_a_val = mat_G_a.values();

    auto G_b_val = mat_G_b.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(weights, vrho, G_a_val, G_b_val, chi_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                G_a_val[nu_offset + g] = weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];

                G_b_val[nu_offset + g] = weights[g] * vrho[2 * g + 1] * chi_val[nu_offset + g];
            }
        }
    }

    timer.stop("Vxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix matmul");

    auto mat_Vxc_a = denblas::multABt(gtoValues, mat_G_a);

    auto mat_Vxc_b = denblas::multABt(gtoValues, mat_G_b);

    timer.stop("Vxc matrix matmul");

    return std::vector<CDenseMatrix>{mat_Vxc_a, mat_Vxc_b};
}

CDenseMatrix
CXCNewIntegrator::_integratePartialVxcFockForGGA(const int32_t          npoints,
                                                 const double*          weights,
                                                 const CDenseMatrix&    gtoValues,
                                                 const CDenseMatrix&    gtoValuesX,
                                                 const CDenseMatrix&    gtoValuesY,
                                                 const CDenseMatrix&    gtoValuesZ,
                                                 const double*          rhograd,
                                                 const double*          vrho,
                                                 const double*          vsigma,
                                                 CMultiTimer&           timer) const
{
    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G(naos, npoints);

    CDenseMatrix mat_G_gga(naos, npoints);

    auto G_val = mat_G.values();

    auto G_gga_val = mat_G_gga.values();

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();

    auto chi_y_val = gtoValuesY.values();

    auto chi_z_val = gtoValuesZ.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(weights, vrho, vsigma, rhograd, G_val, G_gga_val, \
                                     chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                auto vx = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                auto vy = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                auto vz = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                G_val[nu_offset + g] = weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];

                G_gga_val[nu_offset + g] = weights[g] * (vx * chi_x_val[nu_offset + g] +
                                                         vy * chi_y_val[nu_offset + g] +
                                                         vz * chi_z_val[nu_offset + g]);
            }
        }
    }

    timer.stop("Vxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

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

std::vector<CDenseMatrix>
CXCNewIntegrator::_integratePartialVxcFockForGGAOpenShell(const int32_t          npoints,
                                                          const double*          weights,
                                                          const CDenseMatrix&    gtoValues,
                                                          const CDenseMatrix&    gtoValuesX,
                                                          const CDenseMatrix&    gtoValuesY,
                                                          const CDenseMatrix&    gtoValuesZ,
                                                          const double*          rhograd,
                                                          const double*          vrho,
                                                          const double*          vsigma,
                                                          CMultiTimer&           timer) const
{
    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    CDenseMatrix mat_G_a(naos, npoints);
    CDenseMatrix mat_G_b(naos, npoints);

    CDenseMatrix mat_G_a_gga(naos, npoints);
    CDenseMatrix mat_G_b_gga(naos, npoints);

    auto G_a_val = mat_G_a.values();
    auto G_b_val = mat_G_b.values();

    auto G_a_gga_val = mat_G_a_gga.values();
    auto G_b_gga_val = mat_G_b_gga.values();

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(weights, vrho, vsigma, rhograd, G_a_val, G_b_val, G_a_gga_val, G_b_gga_val, \
                                     chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                auto vxa = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                auto vya = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                auto vza = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                auto vxb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 3] + vsigma[3 * g + 1] * rhograd[6 * g + 0];
                auto vyb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 4] + vsigma[3 * g + 1] * rhograd[6 * g + 1];
                auto vzb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 5] + vsigma[3 * g + 1] * rhograd[6 * g + 2];

                G_a_val[nu_offset + g] = weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];
                G_b_val[nu_offset + g] = weights[g] * vrho[2 * g + 1] * chi_val[nu_offset + g];

                G_a_gga_val[nu_offset + g] = weights[g] * (vxa * chi_x_val[nu_offset + g] +
                                                           vya * chi_y_val[nu_offset + g] +
                                                           vza * chi_z_val[nu_offset + g]);
                G_b_gga_val[nu_offset + g] = weights[g] * (vxb * chi_x_val[nu_offset + g] +
                                                           vyb * chi_y_val[nu_offset + g] +
                                                           vzb * chi_z_val[nu_offset + g]);
            }
        }
    }

    timer.stop("Vxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    // Note that we use matrix-matrix multiplication only once, and symmetrize
    // the result. This is because the density matrix is symmetric, and the
    // Kohn-Sham matrix from mat_G is also symmetric. Formally only the
    // mat_G_gga contribution should be symmetrized.

    timer.start("Vxc matrix matmul");

    auto mat_Vxc_a = denblas::multABt(gtoValues, denblas::addAB(mat_G_a, mat_G_a_gga, 2.0));
    auto mat_Vxc_b = denblas::multABt(gtoValues, denblas::addAB(mat_G_b, mat_G_b_gga, 2.0));

    mat_Vxc_a.symmetrizeAndScale(0.5);
    mat_Vxc_b.symmetrizeAndScale(0.5);

    timer.stop("Vxc matrix matmul");

    return std::vector<CDenseMatrix>{mat_Vxc_a, mat_Vxc_b};
}

CDenseMatrix
CXCNewIntegrator::_integratePartialVxcFockForMGGA(const int32_t          npoints,
                                                  const double*          weights,
                                                  const CDenseMatrix&    gtoValues,
                                                  const CDenseMatrix&    gtoValuesX,
                                                  const CDenseMatrix&    gtoValuesY,
                                                  const CDenseMatrix&    gtoValuesZ,
                                                  const double*          rhograd,
                                                  const double*          vrho,
                                                  const double*          vsigma,
                                                  const double*          vlapl,
                                                  const double*          vtau,
                                                  CMultiTimer&           timer) const
{
    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    // LDA contribution
    CDenseMatrix mat_G(naos, npoints);

    // GGA contribution
    CDenseMatrix mat_G_gga(naos, npoints);

    // tau contribution
    CDenseMatrix mat_G_gga_x(naos, npoints);
    CDenseMatrix mat_G_gga_y(naos, npoints);
    CDenseMatrix mat_G_gga_z(naos, npoints);

    auto G_val = mat_G.values();

    auto G_gga_val = mat_G_gga.values();

    auto G_gga_x_val = mat_G_gga_x.values();
    auto G_gga_y_val = mat_G_gga_y.values();
    auto G_gga_z_val = mat_G_gga_z.values();

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(weights, vrho, vsigma, rhograd, vtau, \
                                     G_val, G_gga_val, G_gga_x_val, G_gga_y_val, G_gga_z_val, \
                                     chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                auto vx = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                auto vy = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                auto vz = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                // LDA contribution
                G_val[nu_offset + g] = weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];

                // GGA contribution (will be scaled by 2 later)
                G_gga_val[nu_offset + g] = weights[g] * (vx * chi_x_val[nu_offset + g] +
                                                         vy * chi_y_val[nu_offset + g] +
                                                         vz * chi_z_val[nu_offset + g]);

                // TODO implement Laplacian dependence

                // tau contribution (will be scaled by 0.5 later)
                G_gga_x_val[nu_offset + g] = weights[g] * vtau[2 * g + 0] * chi_x_val[nu_offset + g];
                G_gga_y_val[nu_offset + g] = weights[g] * vtau[2 * g + 0] * chi_y_val[nu_offset + g];
                G_gga_z_val[nu_offset + g] = weights[g] * vtau[2 * g + 0] * chi_z_val[nu_offset + g];
            }
        }
    }

    timer.stop("Vxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    // Note that we use matrix-matrix multiplication only once, and symmetrize
    // the result. This is because the density matrix is symmetric, and the
    // Kohn-Sham matrix from mat_G is also symmetric. Formally only the
    // mat_G_gga contribution should be symmetrized.

    timer.start("Vxc matrix matmul");

    // LDA and GGA contribution
    auto mat_Vxc = denblas::multABt(gtoValues, denblas::addAB(mat_G, mat_G_gga, 2.0));

    // tau contribution
    auto mat_Vxc_x = denblas::multABt(gtoValuesX, mat_G_gga_x);
    auto mat_Vxc_y = denblas::multABt(gtoValuesY, mat_G_gga_y);
    auto mat_Vxc_z = denblas::multABt(gtoValuesZ, mat_G_gga_z);

    mat_Vxc = denblas::addAB(mat_Vxc, mat_Vxc_x, 0.5);
    mat_Vxc = denblas::addAB(mat_Vxc, mat_Vxc_y, 0.5);
    mat_Vxc = denblas::addAB(mat_Vxc, mat_Vxc_z, 0.5);

    mat_Vxc.symmetrizeAndScale(0.5);

    timer.stop("Vxc matrix matmul");

    return mat_Vxc;
}

std::vector<CDenseMatrix>
CXCNewIntegrator::_integratePartialVxcFockForMGGAOpenShell(const int32_t          npoints,
                                                           const double*          weights,
                                                           const CDenseMatrix&    gtoValues,
                                                           const CDenseMatrix&    gtoValuesX,
                                                           const CDenseMatrix&    gtoValuesY,
                                                           const CDenseMatrix&    gtoValuesZ,
                                                           const double*          rhograd,
                                                           const double*          vrho,
                                                           const double*          vsigma,
                                                           const double*          vlapl,
                                                           const double*          vtau,
                                                           CMultiTimer&           timer) const
{
    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    // LDA contribution
    CDenseMatrix mat_G_a(naos, npoints);
    CDenseMatrix mat_G_b(naos, npoints);

    // GGA contribution
    CDenseMatrix mat_G_a_gga(naos, npoints);
    CDenseMatrix mat_G_b_gga(naos, npoints);

    // tau contribution
    CDenseMatrix mat_G_a_gga_x(naos, npoints);
    CDenseMatrix mat_G_a_gga_y(naos, npoints);
    CDenseMatrix mat_G_a_gga_z(naos, npoints);

    CDenseMatrix mat_G_b_gga_x(naos, npoints);
    CDenseMatrix mat_G_b_gga_y(naos, npoints);
    CDenseMatrix mat_G_b_gga_z(naos, npoints);

    auto G_a_val = mat_G_a.values();
    auto G_b_val = mat_G_b.values();

    auto G_a_gga_val = mat_G_a_gga.values();
    auto G_b_gga_val = mat_G_b_gga.values();

    auto G_a_gga_x_val = mat_G_a_gga_x.values();
    auto G_a_gga_y_val = mat_G_a_gga_y.values();
    auto G_a_gga_z_val = mat_G_a_gga_z.values();

    auto G_b_gga_x_val = mat_G_b_gga_x.values();
    auto G_b_gga_y_val = mat_G_b_gga_y.values();
    auto G_b_gga_z_val = mat_G_b_gga_z.values();

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(weights, vrho, vsigma, rhograd, vtau, \
                                     G_a_val, G_a_gga_val, G_a_gga_x_val, G_a_gga_y_val, G_a_gga_z_val, \
                                     G_b_val, G_b_gga_val, G_b_gga_x_val, G_b_gga_y_val, G_b_gga_z_val, \
                                     chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                auto vxa = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                auto vya = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                auto vza = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                auto vxb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 3] + vsigma[3 * g + 1] * rhograd[6 * g + 0];
                auto vyb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 4] + vsigma[3 * g + 1] * rhograd[6 * g + 1];
                auto vzb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 5] + vsigma[3 * g + 1] * rhograd[6 * g + 2];

                // LDA contribution
                G_a_val[nu_offset + g] = weights[g] * vrho[2 * g + 0] * chi_val[nu_offset + g];
                G_b_val[nu_offset + g] = weights[g] * vrho[2 * g + 1] * chi_val[nu_offset + g];

                // GGA contribution (will be scaled by 2 later)
                G_a_gga_val[nu_offset + g] = weights[g] * (vxa * chi_x_val[nu_offset + g] +
                                                           vya * chi_y_val[nu_offset + g] +
                                                           vza * chi_z_val[nu_offset + g]);
                G_b_gga_val[nu_offset + g] = weights[g] * (vxb * chi_x_val[nu_offset + g] +
                                                           vyb * chi_y_val[nu_offset + g] +
                                                           vzb * chi_z_val[nu_offset + g]);

                // TODO implement Laplacian dependence

                // tau contribution (will be scaled by 0.5 later)
                G_a_gga_x_val[nu_offset + g] = weights[g] * vtau[2 * g + 0] * chi_x_val[nu_offset + g];
                G_a_gga_y_val[nu_offset + g] = weights[g] * vtau[2 * g + 0] * chi_y_val[nu_offset + g];
                G_a_gga_z_val[nu_offset + g] = weights[g] * vtau[2 * g + 0] * chi_z_val[nu_offset + g];

                G_b_gga_x_val[nu_offset + g] = weights[g] * vtau[2 * g + 1] * chi_x_val[nu_offset + g];
                G_b_gga_y_val[nu_offset + g] = weights[g] * vtau[2 * g + 1] * chi_y_val[nu_offset + g];
                G_b_gga_z_val[nu_offset + g] = weights[g] * vtau[2 * g + 1] * chi_z_val[nu_offset + g];
            }
        }
    }

    timer.stop("Vxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    // Note that we use matrix-matrix multiplication only once, and symmetrize
    // the result. This is because the density matrix is symmetric, and the
    // Kohn-Sham matrix from mat_G is also symmetric. Formally only the
    // mat_G_gga contribution should be symmetrized.

    timer.start("Vxc matrix matmul");

    // LDA and GGA contribution
    auto mat_Vxc_a = denblas::multABt(gtoValues, denblas::addAB(mat_G_a, mat_G_a_gga, 2.0));
    auto mat_Vxc_b = denblas::multABt(gtoValues, denblas::addAB(mat_G_b, mat_G_b_gga, 2.0));

    // tau contribution
    auto mat_Vxc_a_x = denblas::multABt(gtoValuesX, mat_G_a_gga_x);
    auto mat_Vxc_a_y = denblas::multABt(gtoValuesY, mat_G_a_gga_y);
    auto mat_Vxc_a_z = denblas::multABt(gtoValuesZ, mat_G_a_gga_z);

    auto mat_Vxc_b_x = denblas::multABt(gtoValuesX, mat_G_b_gga_x);
    auto mat_Vxc_b_y = denblas::multABt(gtoValuesY, mat_G_b_gga_y);
    auto mat_Vxc_b_z = denblas::multABt(gtoValuesZ, mat_G_b_gga_z);

    mat_Vxc_a = denblas::addAB(mat_Vxc_a, mat_Vxc_a_x, 0.5);
    mat_Vxc_a = denblas::addAB(mat_Vxc_a, mat_Vxc_a_y, 0.5);
    mat_Vxc_a = denblas::addAB(mat_Vxc_a, mat_Vxc_a_z, 0.5);

    mat_Vxc_b = denblas::addAB(mat_Vxc_b, mat_Vxc_b_x, 0.5);
    mat_Vxc_b = denblas::addAB(mat_Vxc_b, mat_Vxc_b_y, 0.5);
    mat_Vxc_b = denblas::addAB(mat_Vxc_b, mat_Vxc_b_z, 0.5);

    mat_Vxc_a.symmetrizeAndScale(0.5);
    mat_Vxc_b.symmetrizeAndScale(0.5);

    timer.stop("Vxc matrix matmul");

    return std::vector<CDenseMatrix>{mat_Vxc_a, mat_Vxc_b};
}

CDenseMatrix
CXCNewIntegrator::_integratePartialFxcFockForLDA(const int32_t         npoints,
                                                 const double*         weights,
                                                 const CDenseMatrix&   gtoValues,
                                                 const double*         rhow,
                                                 const double*         v2rho2,
                                                 CMultiTimer&          timer) const
{
    // GTO values on grid points

    auto chi_val = gtoValues.values();

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

            #pragma omp simd aligned(weights, v2rho2, rhow, G_val, chi_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                G_val[nu_offset + g] = weights[g] * (v2rho2[3 * g + 0] * rhow[2 * g + 0] +
                                                     v2rho2[3 * g + 1] * rhow[2 * g + 1]) * chi_val[nu_offset + g];
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
CXCNewIntegrator::_integratePartialFxcFockForGGA(const int32_t          npoints,
                                                 const double*          weights,
                                                 const CDenseMatrix&    gtoValues,
                                                 const CDenseMatrix&    gtoValuesX,
                                                 const CDenseMatrix&    gtoValuesY,
                                                 const CDenseMatrix&    gtoValuesZ,
                                                 const double*          rhow,
                                                 const double*          rhograd,
                                                 const double*          rhowgrad,
                                                 const double*          vsigma,
                                                 const double*          v2rho2,
                                                 const double*          v2rhosigma,
                                                 const double*          v2sigma2,
                                                 CMultiTimer&           timer) const
{
    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();

    auto chi_y_val = gtoValuesY.values();

    auto chi_z_val = gtoValuesZ.values();

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

            #pragma omp simd aligned(weights, rhow, rhograd, rhowgrad, \
                    vsigma, v2rho2, v2rhosigma, v2sigma2, \
                    G_val, G_gga_val, chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                double w = weights[g];

                auto grhow_grho_aa = 2.0 * (rhowgrad[6 * g + 0] * rhograd[6 * g + 0] +
                                            rhowgrad[6 * g + 1] * rhograd[6 * g + 1] +
                                            rhowgrad[6 * g + 2] * rhograd[6 * g + 2]);

                auto grhow_grho_bb = 2.0 * (rhowgrad[6 * g + 3] * rhograd[6 * g + 3] +
                                            rhowgrad[6 * g + 4] * rhograd[6 * g + 4] +
                                            rhowgrad[6 * g + 5] * rhograd[6 * g + 5]);

                auto grhow_grho_ab = (rhowgrad[6 * g + 0] * rhograd[6 * g + 3] +
                                      rhowgrad[6 * g + 1] * rhograd[6 * g + 4] +
                                      rhowgrad[6 * g + 2] * rhograd[6 * g + 5] +

                                      rhowgrad[6 * g + 3] * rhograd[6 * g + 0] +
                                      rhowgrad[6 * g + 4] * rhograd[6 * g + 1] +
                                      rhowgrad[6 * g + 5] * rhograd[6 * g + 2]);

                // scalar contribution

                double f_0 = v2rho2[3 * g + 0] * rhow[2 * g + 0] +
                             v2rho2[3 * g + 1] * rhow[2 * g + 1] +

                             v2rhosigma[6 * g + 0] * grhow_grho_aa +
                             v2rhosigma[6 * g + 1] * grhow_grho_ab +
                             v2rhosigma[6 * g + 2] * grhow_grho_bb;

                G_val[nu_offset + g] = w * f_0 * chi_val[nu_offset + g];

                // vector contribution

                double f_aa = v2rhosigma[6 * g + 0] * rhow[2 * g + 0] +
                              v2rhosigma[6 * g + 3] * rhow[2 * g + 1] +

                              v2sigma2[6 * g + 0] * grhow_grho_aa +
                              v2sigma2[6 * g + 1] * grhow_grho_ab +
                              v2sigma2[6 * g + 2] * grhow_grho_bb;

                double f_ab = v2rhosigma[6 * g + 1] * rhow[2 * g + 0] +
                              v2rhosigma[6 * g + 4] * rhow[2 * g + 1] +

                              v2sigma2[6 * g + 1] * grhow_grho_aa +
                              v2sigma2[6 * g + 3] * grhow_grho_ab +
                              v2sigma2[6 * g + 4] * grhow_grho_bb;

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                xcomp += 2.0 * f_aa * rhograd[6 * g + 0] + f_ab * rhograd[6 * g + 3];
                ycomp += 2.0 * f_aa * rhograd[6 * g + 1] + f_ab * rhograd[6 * g + 4];
                zcomp += 2.0 * f_aa * rhograd[6 * g + 2] + f_ab * rhograd[6 * g + 5];

                xcomp += 2.0 * vsigma[3 * g + 0] * rhowgrad[6 * g + 0] + vsigma[3 * g + 1] * rhowgrad[6 * g + 3];
                ycomp += 2.0 * vsigma[3 * g + 0] * rhowgrad[6 * g + 1] + vsigma[3 * g + 1] * rhowgrad[6 * g + 4];
                zcomp += 2.0 * vsigma[3 * g + 0] * rhowgrad[6 * g + 2] + vsigma[3 * g + 1] * rhowgrad[6 * g + 5];

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
CXCNewIntegrator::_integratePartialFxcFockForMGGA(const int32_t       npoints, 
                                                  const double*       weights, 
                                                  const CDenseMatrix& gtoValues,
                                                  const CDenseMatrix& gtoValuesX, 
                                                  const CDenseMatrix& gtoValuesY, 
                                                  const CDenseMatrix& gtoValuesZ,
                                                  const double*       rhow, 
                                                  const double*       rhograd, 
                                                  const double*       rhowgrad, 
                                                  const double*       tauw, 
                                                  const double*       laplw,
                                                  const double*       vrho, 
                                                  const double*       vsigma, 
                                                  const double*       vlapl, 
                                                  const double*       vtau, 
                                                  const double*       v2rho2,
                                                  const double*       v2lapl2, 
                                                  const double*       v2tau2, 
                                                  const double*       v2rholapl, 
                                                  const double*       v2rhotau,
                                                  const double*       v2lapltau, 
                                                  const double*       v2rhosigma, 
                                                  const double*       v2sigmalapl, 
                                                  const double*       v2sigmatau,
                                                  const double*       v2sigma2, 
                                                  CMultiTimer&        timer) const
{
    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Fxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    // LDA contribution
    CDenseMatrix mat_G(naos, npoints);

    // GGA contribution
    CDenseMatrix mat_G_gga(naos, npoints);

    // tau contribution
    CDenseMatrix mat_G_gga_x(naos, npoints);
    CDenseMatrix mat_G_gga_y(naos, npoints);
    CDenseMatrix mat_G_gga_z(naos, npoints);

    auto G_val = mat_G.values();

    auto G_gga_val = mat_G_gga.values();

    auto G_gga_x_val = mat_G_gga_x.values();
    auto G_gga_y_val = mat_G_gga_y.values();
    auto G_gga_z_val = mat_G_gga_z.values();

    auto chi_val = gtoValues.values();
    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(weights, rhow, rhograd, rhowgrad, tauw, laplw, \
                                     vrho, vsigma, vlapl, vtau, v2rho2, v2lapl2, v2tau2, v2rholapl, v2rhotau, v2lapltau, \
                                     v2rhosigma, v2sigmalapl, v2sigmatau, v2sigma2, \
                                     G_val, G_gga_val, G_gga_x_val, G_gga_y_val, G_gga_z_val, \
                                     chi_val, chi_x_val, chi_y_val, chi_z_val :VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                double w = weights[g];

                // ground-state gardients
                double grada_x_g = rhograd[6 * g + 0];
                double grada_y_g = rhograd[6 * g + 1];
                double grada_z_g = rhograd[6 * g + 2];

                // perturbed density and its gradients
                double rwa = rhow[2 * g + 0];
                double tauwa = tauw[2 * g + 0];
                double laplwa = laplw[2 * g + 0];

                double rwa_x = rhowgrad[6 * g + 0];
                double rwa_y = rhowgrad[6 * g + 1];
                double rwa_z = rhowgrad[6 * g + 2];

                // functional derivatives 
                // first-order
                double vsigma_a = vsigma[3 * g + 0];
                double vsigma_c = vsigma[3 * g + 1];
                // second-order
                double v2rho2_aa = v2rho2[3 * g + 0];
                double v2rho2_ab = v2rho2[3 * g + 1];

                double v2rhosigma_aa = v2rhosigma[6 * g + 0];
                double v2rhosigma_ac = v2rhosigma[6 * g + 1];
                double v2rhosigma_ab = v2rhosigma[6 * g + 2];
                double v2rhosigma_ba = v2rhosigma[6 * g + 3];
                double v2rhosigma_bc = v2rhosigma[6 * g + 4];
                double v2sigma2_aa = v2sigma2[6 * g + 0];
                double v2sigma2_ac = v2sigma2[6 * g + 1];
                double v2sigma2_ab = v2sigma2[6 * g + 2];
                double v2sigma2_cc = v2sigma2[6 * g + 3];
                double v2sigma2_cb = v2sigma2[6 * g + 4];

                // second-order meta-gga
                double v2rholapl_aa = v2rholapl[4 * g + 0];
                double v2rholapl_ab = v2rholapl[4 * g + 1];
                double v2rholapl_ba = v2rholapl[4 * g + 2];

                double v2rhotau_aa = v2rhotau[4 * g + 0]; 
                double v2rhotau_ab = v2rhotau[4 * g + 1];
                double v2rhotau_ba = v2rhotau[4 * g + 2]; 

                double v2lapltau_aa = v2lapltau[4 * g + 0];
                double v2lapltau_ab = v2lapltau[4 * g + 1];
                double v2lapltau_ba = v2lapltau[4 * g + 2];

                double v2lapl2_aa = v2lapl2[3 * g + 0]; 
                double v2lapl2_ab = v2lapl2[3 * g + 1]; 

                double v2tau2_aa = v2tau2[3 * g + 0]; 
                double v2tau2_ab = v2tau2[3 * g + 1]; 

                double v2sigmalapl_aa = v2sigmalapl[6 * g + 0];
                double v2sigmalapl_ab = v2sigmalapl[6 * g + 1];
                double v2sigmalapl_ca = v2sigmalapl[6 * g + 2];
                double v2sigmalapl_cb = v2sigmalapl[6 * g + 3];
                double v2sigmalapl_ba = v2sigmalapl[6 * g + 4];

                double v2sigmatau_aa = v2sigmatau[6 * g + 0]; 
                double v2sigmatau_ab = v2sigmatau[6 * g + 1]; 
                double v2sigmatau_ca = v2sigmatau[6 * g + 2]; 
                double v2sigmatau_cb = v2sigmatau[6 * g + 3]; 
                double v2sigmatau_ba = v2sigmatau[6 * g + 4]; 

                // sums of functional derivatives that can be used in the restricted case

                // first-order
                double x = vsigma_c + 2.0 * vsigma_a;
                // second-order 
                double rr = v2rho2_aa + v2rho2_ab;
                double rx = 2.0*v2rhosigma_ac + 2.0*v2rhosigma_ab + 2.0*v2rhosigma_aa;
                double rt = v2rhotau_aa + v2rhotau_ab;
                double rl = v2rholapl_aa + v2rholapl_ab;

                // sigma and gamma
                double xr = v2rhosigma_bc + 2.0*v2rhosigma_ba + v2rhosigma_ac + 2.0 * v2rhosigma_aa;
                double xt = v2sigmatau_cb + 2.0*v2sigmatau_ab + v2sigmatau_ca + 2.0 * v2sigmatau_aa;
                double xl = v2sigmalapl_cb + 2.0*v2sigmalapl_ab + v2sigmalapl_ca + 2.0 * v2sigmalapl_aa;
                double xx = 2.0*v2sigma2_cc + 2.0*v2sigma2_cb + 6.0*v2sigma2_ac + 4.0*v2sigma2_ab + 4.0 * v2sigma2_aa;

                // tau
                double tt = v2tau2_aa + v2tau2_ab;
                double tx = 2.0 * v2sigmatau_ca + 2.0 * v2sigmatau_ba + 2.0 * v2sigmatau_aa;
                double tr = v2rhotau_aa + v2rhotau_ba;
                double tl = v2lapltau_aa + v2lapltau_ba;

                // lapl
                //double ll = v2lapl2_aa + v2lapl2_ab;
                //double lx = 2.0 * v2sigmalapl_ca + 2.0 * v2sigmalapl_ba + 2.0 * v2sigmalapl_aa;
                //double lr = v2rholapl_aa + v2rholapl_ba;
                //double lt = v2lapltau_aa + v2lapltau_ab;

                // contraction of perturbed density for restricted case
                double contract = grada_x_g * rwa_x + grada_y_g * rwa_y + grada_z_g * rwa_z;

                // rho-operator
                double r_0 =  rr * rwa 
                            + rx * contract 
                            + rt * tauwa 
                            + rl * laplwa; 

                G_val[nu_offset + g] = w * r_0 * chi_val[nu_offset + g];

                // GGA contribution (will be scaled by 2 later)

                // vector contribution

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                // grad-operator 

                xcomp +=  grada_x_g * ( xr * rwa + xt * tauwa + xl * laplwa)
                        + x * rwa_x 
                        + xx * grada_x_g * contract;

                ycomp += grada_y_g * ( xr * rwa + xt * tauwa + xl * laplwa)
                        + x * rwa_y 
                        + xx * grada_y_g * contract;

                zcomp += grada_z_g * ( xr * rwa + xt * tauwa + xl * laplwa)
                        + x * rwa_z 
                        + xx * grada_z_g * contract;

                G_gga_val[nu_offset + g] = w * (xcomp * chi_x_val[nu_offset + g] +
                                                ycomp * chi_y_val[nu_offset + g] +
                                                zcomp * chi_z_val[nu_offset + g]);

                // TODO implement Laplacian dependence

                //double lap_0 =    lr * rwa 
                //                + lx * contract
                //                + lt * tauwa 
                //                + ll * laplwa; 
                //G_gga_val[nu_offset + g] += w * lap_0 * (chi_xx_val[nu_offset + g] +
                //                                         chi_yy_val[nu_offset + g] +
                //                                         chi_zz_val[nu_offset + g]);

                // tau contribution (will be scaled by 0.5 later)
                double tau_0 =    tr * rwa 
                                + tx * contract
                                + tt * tauwa 
                                + tl * laplwa;   

                G_gga_x_val[nu_offset + g] = w * tau_0 * chi_x_val[nu_offset + g];
                G_gga_y_val[nu_offset + g] = w * tau_0 * chi_y_val[nu_offset + g];
                G_gga_z_val[nu_offset + g] = w * tau_0 * chi_z_val[nu_offset + g];
            }
        }
    }

    timer.stop("Fxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    // Note that we use matrix-matrix multiplication only once, and symmetrize
    // the result. This is because the density matrix is symmetric, and the
    // Kohn-Sham matrix from mat_G is also symmetric. Formally only the
    // mat_G_gga contribution should be symmetrized.

    timer.start("Fxc matrix matmul");

    // LDA and GGA contribution
    auto mat_Fxc = denblas::multABt(gtoValues, denblas::addAB(mat_G, mat_G_gga, 2.0));

    // tau contribution
    auto mat_Fxc_x = denblas::multABt(gtoValuesX, mat_G_gga_x);
    auto mat_Fxc_y = denblas::multABt(gtoValuesY, mat_G_gga_y);
    auto mat_Fxc_z = denblas::multABt(gtoValuesZ, mat_G_gga_z);

    mat_Fxc = denblas::addAB(mat_Fxc, mat_Fxc_x, 0.5);
    mat_Fxc = denblas::addAB(mat_Fxc, mat_Fxc_y, 0.5);
    mat_Fxc = denblas::addAB(mat_Fxc, mat_Fxc_z, 0.5);

    mat_Fxc.symmetrizeAndScale(0.5);

    timer.stop("Fxc matrix matmul");

    return mat_Fxc;
}

CDenseMatrix
CXCNewIntegrator::_integratePartialKxcFockForLDA(const int32_t              npoints,
                                                 const double*              weights,
                                                 const CDenseMatrix&        gtoValues,
                                                 const double*              v2rho2,
                                                 const double*              v3rho3,
                                                 const CDensityGridQuad&    rwDensityGridQuad,
                                                 const CDensityGrid&        rw2DensityGrid,
                                                 const int32_t              iFock,
                                                 CMultiTimer&               timer) const
{
    timer.start("Kxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    // pointers to perturbed density

    auto rhow1a = rwDensityGridQuad.gam(iFock);

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

            #pragma omp simd aligned(weights, v2rho2, v3rho3, rhow1a, rhow12a, rhow12b, G_val, chi_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                G_val[nu_offset + g] = weights[g] *

                          ((v3rho3[4 * g + 0] + 2.0 * v3rho3[4 * g + 1] + v3rho3[4 * g + 2]) * rhow1a[g] +

                           v2rho2[3 * g + 0] * rhow12a[g] + v2rho2[3 * g + 1] * rhow12b[g]) *

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
CXCNewIntegrator::_integratePartialKxcFockForGGA(const int32_t           npoints,
                                                 const double*           weights,
                                                 const CDenseMatrix&     gtoValues,
                                                 const CDenseMatrix&     gtoValuesX,
                                                 const CDenseMatrix&     gtoValuesY,
                                                 const CDenseMatrix&     gtoValuesZ,
                                                 const double*           rhograd,
                                                 const double*           vsigma,
                                                 const double*           v2rho2,
                                                 const double*           v2rhosigma,
                                                 const double*           v2sigma2,
                                                 const double*           v3rho3,
                                                 const double*           v3rho2sigma,
                                                 const double*           v3rhosigma2,
                                                 const double*           v3sigma3,
                                                 const CDensityGridQuad& rwDensityGridQuad,
                                                 const CDensityGrid&     rw2DensityGrid,
                                                 const int32_t           iFock,
                                                 CMultiTimer&            timer) const
{
    timer.start("Kxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    // pointers to perturbed density gradient norms

    auto gam = rwDensityGridQuad.gam(iFock);

    auto gamx = rwDensityGridQuad.gamX(iFock);
    auto gamy = rwDensityGridQuad.gamY(iFock);
    auto gamz = rwDensityGridQuad.gamZ(iFock);

    auto gamxx = rwDensityGridQuad.gamXX(iFock);
    auto gamxy = rwDensityGridQuad.gamXY(iFock);
    auto gamxz = rwDensityGridQuad.gamXZ(iFock);

    auto gamyx = rwDensityGridQuad.gamYX(iFock);
    auto gamyy = rwDensityGridQuad.gamYY(iFock);
    auto gamyz = rwDensityGridQuad.gamYZ(iFock);

    auto gamzx = rwDensityGridQuad.gamZX(iFock);
    auto gamzy = rwDensityGridQuad.gamZY(iFock);
    auto gamzz = rwDensityGridQuad.gamZZ(iFock);

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
                    rhograd, vsigma, v2rho2, v2rhosigma, v2sigma2, \
                    v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3, \
                    gam, gamx, gamy, gamz, \
                    gamxx, gamxy, gamxz, gamyx, gamyy, gamyz, gamzx, gamzy, gamzz, \
                    rhow12a, gradw12a_x, gradw12a_y, gradw12a_z, \
                    G_val, G_gga_val, chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                double w = weights[g];

                double rxw12a = gradw12a_x[g];
                double ryw12a = gradw12a_y[g];
                double rzw12a = gradw12a_z[g];

                double grada_x_g = rhograd[6 * g + 0];
                double grada_y_g = rhograd[6 * g + 1];
                double grada_z_g = rhograd[6 * g + 2];

                double l2contract = grada_x_g * rxw12a + grada_y_g * ryw12a + grada_z_g * rzw12a;
                double l5contract_x = grada_x_g * l2contract;
                double l5contract_y = grada_y_g * l2contract;
                double l5contract_z = grada_z_g * l2contract;
                double q2contract = grada_x_g * gamx[g] + grada_y_g * gamy[g] + grada_z_g * gamz[g];
                double q3contract =   grada_x_g * grada_x_g * gamxx[g]
                                    + grada_x_g * grada_y_g * gamxy[g]
                                    + grada_x_g * grada_z_g * gamxz[g]
                                    + grada_y_g * grada_x_g * gamyx[g]
                                    + grada_y_g * grada_y_g * gamyy[g]
                                    + grada_y_g * grada_z_g * gamyz[g]
                                    + grada_z_g * grada_x_g * gamzx[g]
                                    + grada_z_g * grada_y_g * gamzy[g]
                                    + grada_z_g * grada_z_g * gamzz[g];

                double q4contract = gamxx[g] + gamyy[g] + gamzz[g];
                double q7contract_x =  grada_x_g * (grada_x_g * gamx[g] + grada_y_g * gamy[g] + grada_z_g * gamz[g]);
                double q7contract_y =  grada_y_g * (grada_x_g * gamx[g] + grada_y_g * gamy[g] + grada_z_g * gamz[g]);
                double q7contract_z =  grada_z_g * (grada_x_g * gamx[g] + grada_y_g * gamy[g] + grada_z_g * gamz[g]);

                double q8contract_x =  grada_x_g * gamxx[g] + grada_y_g * gamxy[g] + grada_z_g * gamxz[g];
                double q8contract_y =  grada_x_g * gamyx[g] + grada_y_g * gamyy[g] + grada_z_g * gamyz[g];
                double q8contract_z =  grada_x_g * gamzx[g] + grada_y_g * gamzy[g] + grada_z_g * gamzz[g];

                double q9contract_x =  grada_x_g * q3contract;
                double q9contract_y =  grada_y_g * q3contract;
                double q9contract_z =  grada_z_g * q3contract;

                double q10contract_x =  grada_x_g * gamxx[g] + grada_y_g * gamyx[g] + grada_z_g * gamzx[g];
                double q10contract_y =  grada_x_g * gamxy[g] + grada_y_g * gamyy[g] + grada_z_g * gamzy[g];
                double q10contract_z =  grada_x_g * gamxz[g] + grada_y_g * gamyz[g] + grada_z_g * gamzz[g];

                double q11contract_x =  grada_x_g * gamxx[g] + grada_x_g * gamyy[g] + grada_x_g * gamzz[g];
                double q11contract_y =  grada_y_g * gamxx[g] + grada_y_g * gamyy[g] + grada_y_g * gamzz[g];
                double q11contract_z =  grada_z_g * gamxx[g] + grada_z_g * gamyy[g] + grada_z_g * gamzz[g];

                // functional derivatives in libxc form

                auto vsigma_a = vsigma[3 * g + 0];
                auto vsigma_c = vsigma[3 * g + 1];

                auto v2rho2_aa = v2rho2[3 * g + 0];
                auto v2rho2_ab = v2rho2[3 * g + 1];

                auto v2rhosigma_aa = v2rhosigma[6 * g + 0];
                auto v2rhosigma_ac = v2rhosigma[6 * g + 1];
                auto v2rhosigma_ab = v2rhosigma[6 * g + 2];
                auto v2rhosigma_ba = v2rhosigma[6 * g + 3];
                auto v2rhosigma_bc = v2rhosigma[6 * g + 4];

                auto v2sigma2_aa = v2sigma2[6 * g + 0];
                auto v2sigma2_ac = v2sigma2[6 * g + 1];
                auto v2sigma2_ab = v2sigma2[6 * g + 2];
                auto v2sigma2_cc = v2sigma2[6 * g + 3];
                auto v2sigma2_cb = v2sigma2[6 * g + 4];

                auto v3rho3_aaa = v3rho3[4 * g + 0];
                auto v3rho3_aab = v3rho3[4 * g + 1];
                auto v3rho3_abb = v3rho3[4 * g + 2];

                auto v3rho2sigma_aaa = v3rho2sigma[9 * g + 0];
                auto v3rho2sigma_aac = v3rho2sigma[9 * g + 1];
                auto v3rho2sigma_aab = v3rho2sigma[9 * g + 2];
                auto v3rho2sigma_aba = v3rho2sigma[9 * g + 3];
                auto v3rho2sigma_abc = v3rho2sigma[9 * g + 4];
                auto v3rho2sigma_abb = v3rho2sigma[9 * g + 5];
                auto v3rho2sigma_bba = v3rho2sigma[9 * g + 6];
                auto v3rho2sigma_bbc = v3rho2sigma[9 * g + 7];

                auto v3rhosigma2_aaa = v3rhosigma2[12 * g + 0];
                auto v3rhosigma2_aac = v3rhosigma2[12 * g + 1];
                auto v3rhosigma2_aab = v3rhosigma2[12 * g + 2];
                auto v3rhosigma2_acc = v3rhosigma2[12 * g + 3];
                auto v3rhosigma2_acb = v3rhosigma2[12 * g + 4];
                auto v3rhosigma2_abb = v3rhosigma2[12 * g + 5];
                auto v3rhosigma2_baa = v3rhosigma2[12 * g + 6];
                auto v3rhosigma2_bac = v3rhosigma2[12 * g + 7];
                auto v3rhosigma2_bab = v3rhosigma2[12 * g + 8];
                auto v3rhosigma2_bcc = v3rhosigma2[12 * g + 9];
                auto v3rhosigma2_bcb = v3rhosigma2[12 * g + 10];

                auto v3sigma3_aaa = v3sigma3[10 * g + 0];
                auto v3sigma3_aac = v3sigma3[10 * g + 1];
                auto v3sigma3_aab = v3sigma3[10 * g + 2];
                auto v3sigma3_acc = v3sigma3[10 * g + 3];
                auto v3sigma3_acb = v3sigma3[10 * g + 4];
                auto v3sigma3_abb = v3sigma3[10 * g + 5];
                auto v3sigma3_ccc = v3sigma3[10 * g + 6];
                auto v3sigma3_bcc = v3sigma3[10 * g + 7];
                auto v3sigma3_cbb = v3sigma3[10 * g + 8];

                // functional derivatives
                double rr = (v2rho2_aa + v2rho2_ab);
                double rrr = (v3rho3_aaa + 2.0 * v3rho3_aab + v3rho3_abb);
                double rx = (2.0*v2rhosigma_ac + 2.0*v2rhosigma_ab + 2.0*v2rhosigma_aa);
                double rxr = (2.0*v3rho2sigma_abc + 2.0*v3rho2sigma_abb + 2.0*v3rho2sigma_aba
                             + 2.0*v3rho2sigma_aac + 2.0*v3rho2sigma_aab + 2.0*v3rho2sigma_aaa);
                double rxx = (4.0*v3rhosigma2_acc + 8.0*v3rhosigma2_acb + 4.0*v3rhosigma2_abb
                            + 8.0*v3rhosigma2_aac + 8.0*v3rhosigma2_aab + 4.0*v3rhosigma2_aaa);
                double x = vsigma_c + 2.0*vsigma_a;
                double xr = v2rhosigma_bc + 2.0*v2rhosigma_ba + v2rhosigma_ac + 2.0*v2rhosigma_aa;
                double xx = 2.0*v2sigma2_cc + 2.0*v2sigma2_cb + 6.0*v2sigma2_ac + 4.0*v2sigma2_ab + 4.0*v2sigma2_aa;
                double xrr = v3rho2sigma_bbc + 2.0*v3rho2sigma_bba + 2.0*v3rho2sigma_abc + 4.0*v3rho2sigma_aba
                            + v3rho2sigma_aac + 2.0*v3rho2sigma_aaa;
                double xxr = 2.0*v3rhosigma2_bcc + 2.0*v3rhosigma2_bcb + 6.0*v3rhosigma2_bac
                            + 4.0*v3rhosigma2_bab + 4.0*v3rhosigma2_baa + 2.0*v3rhosigma2_acc
                            + 2.0*v3rhosigma2_acb + 6.0*v3rhosigma2_aac + 4.0*v3rhosigma2_aab + 4.0*v3rhosigma2_aaa;
                double xxx = 4.0*v3sigma3_ccc + 8.0*v3sigma3_bcc + 4.0*v3sigma3_cbb + 16.0*v3sigma3_acc
                            + 24.0*v3sigma3_acb + 8.0*v3sigma3_abb + 20.0*v3sigma3_aac
                            + 16.0*v3sigma3_aab + 8.0*v3sigma3_aaa;

                // Scalar contribution

                double prefac = 0.0;

                // vxc 1 contributions

                prefac += rr * rhow12a[g] // l1
                        + rx * l2contract;

                // vxc 2 contributions

                prefac += rrr * gam[g] // q1
                        + rxr * q2contract
                        + rxx * q3contract
                        + rx * q4contract;

                G_val[nu_offset + g] = w * prefac * chi_val[nu_offset + g];

                // vector contribution

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                // vxc 1 contributions

                xcomp += xr * grada_x_g * rhow12a[g] // l3
                        + x * rxw12a // l4
                        + xx * l5contract_x;

                ycomp += xr * grada_y_g * rhow12a[g] // l3
                        + x * ryw12a // l4
                        + xx * l5contract_y;

                zcomp += xr * grada_z_g * rhow12a[g] // l3
                        + x * rzw12a // l4
                        + xx * l5contract_z;

                // vxc 2 contributions

                xcomp += xrr * grada_x_g * gam[g] // q5
                        + xr * gamx[g] // q6
                        + xxr * q7contract_x
                        + xx * (q8contract_x + q10contract_x + q11contract_x)
                        + xxx * q9contract_x;

                ycomp += xrr * grada_y_g * gam[g] // q5
                        + xr * gamy[g] // q6
                        + xxr * q7contract_y
                        + xx * (q8contract_y + q10contract_y + q11contract_y)
                        + xxx * q9contract_y;

                zcomp += xrr * grada_z_g * gam[g] // q5
                        + xr * gamz[g] // q6
                        + xxr * q7contract_z
                        + xx * (q8contract_z + q10contract_z + q11contract_z)
                        + xxx * q9contract_z;

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
CXCNewIntegrator::_integratePartialKxcFockForMGGA(const int32_t           npoints,
                                                  const double*           weights,
                                                  const CDenseMatrix&     gtoValues,
                                                  const CDenseMatrix&     gtoValuesX,
                                                  const CDenseMatrix&     gtoValuesY,
                                                  const CDenseMatrix&     gtoValuesZ,
                                                  const double*           rhograd,
                                                  const double*           vsigma,
                                                  const double*           v2rho2,
                                                  const double*           v2lapl2,
                                                  const double*           v2tau2,
                                                  const double*           v2rholapl,
                                                  const double*           v2rhotau,
                                                  const double*           v2lapltau,
                                                  const double*           v2rhosigma,
                                                  const double*           v2sigmalapl,
                                                  const double*           v2sigmatau,
                                                  const double*           v2sigma2,
                                                  const double*           v3rho3,
                                                  const double*           v3rho2sigma,
                                                  const double*           v3rho2lapl,
                                                  const double*           v3rho2tau,
                                                  const double*           v3rhosigma2,
                                                  const double*           v3rhosigmalapl,
                                                  const double*           v3rhosigmatau,
                                                  const double*           v3rholapl2,
                                                  const double*           v3rholapltau,
                                                  const double*           v3rhotau2,
                                                  const double*           v3sigma3,
                                                  const double*           v3sigma2lapl,
                                                  const double*           v3sigma2tau,
                                                  const double*           v3sigmalapl2,
                                                  const double*           v3sigmalapltau,
                                                  const double*           v3sigmatau2,
                                                  const double*           v3lapl3,
                                                  const double*           v3lapl2tau,
                                                  const double*           v3lapltau2,
                                                  const double*           v3tau3,
                                                  const CDensityGridQuad& rwDensityGridQuad,
                                                  const CDensityGrid&     rw2DensityGrid,
                                                  const int32_t           iFock,
                                                  CMultiTimer&            timer) const
{
    timer.start("Kxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    // pointers to perturbed densities

    auto gam = rwDensityGridQuad.gam(iFock);
    auto rt_gam = rwDensityGridQuad.rt_gam(iFock);
    auto rl_gam = rwDensityGridQuad.rl_gam(iFock);
    auto tt_gam = rwDensityGridQuad.tt_gam(iFock);
    auto tl_gam = rwDensityGridQuad.tl_gam(iFock);
    auto ll_gam = rwDensityGridQuad.ll_gam(iFock);

    auto gamx = rwDensityGridQuad.gamX(iFock);
    auto gamy = rwDensityGridQuad.gamY(iFock);
    auto gamz = rwDensityGridQuad.gamZ(iFock);
    auto st_gamx = rwDensityGridQuad.st_gamX(iFock);
    auto st_gamy = rwDensityGridQuad.st_gamY(iFock);
    auto st_gamz = rwDensityGridQuad.st_gamZ(iFock);
    auto sl_gamx = rwDensityGridQuad.sl_gamX(iFock);
    auto sl_gamy = rwDensityGridQuad.sl_gamY(iFock);
    auto sl_gamz = rwDensityGridQuad.sl_gamZ(iFock);

    auto gamxx = rwDensityGridQuad.gamXX(iFock);
    auto gamxy = rwDensityGridQuad.gamXY(iFock);
    auto gamxz = rwDensityGridQuad.gamXZ(iFock);
    auto gamyx = rwDensityGridQuad.gamYX(iFock);
    auto gamyy = rwDensityGridQuad.gamYY(iFock);
    auto gamyz = rwDensityGridQuad.gamYZ(iFock);
    auto gamzx = rwDensityGridQuad.gamZX(iFock);
    auto gamzy = rwDensityGridQuad.gamZY(iFock);
    auto gamzz = rwDensityGridQuad.gamZZ(iFock);

    auto rhow12a = rw2DensityGrid.alphaDensity(iFock);
    auto tauw12a = rw2DensityGrid.alphaDensitytau(iFock);
    auto laplw12a = rw2DensityGrid.alphaDensitylapl(iFock);
    auto gradw12a_x = rw2DensityGrid.alphaDensityGradientX(iFock);
    auto gradw12a_y = rw2DensityGrid.alphaDensityGradientY(iFock);
    auto gradw12a_z = rw2DensityGrid.alphaDensityGradientZ(iFock);

    timer.stop("Kxc matrix prep.");

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Kxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    // LDA contribution
    CDenseMatrix mat_G(naos, npoints);

    // GGA contribution
    CDenseMatrix mat_G_gga(naos, npoints);

    // tau contribution
    CDenseMatrix mat_G_gga_x(naos, npoints);
    CDenseMatrix mat_G_gga_y(naos, npoints);
    CDenseMatrix mat_G_gga_z(naos, npoints);

    auto G_val = mat_G.values();

    auto G_gga_val = mat_G_gga.values();

    auto G_gga_x_val = mat_G_gga_x.values();
    auto G_gga_y_val = mat_G_gga_y.values();
    auto G_gga_z_val = mat_G_gga_z.values();

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
                         rhograd, vsigma,\
                         v2rho2,v2lapl2, v2tau2, v2rholapl, v2rhotau,v2lapltau, v2rhosigma,\
                         v2sigmalapl, v2sigmatau,v2sigma2,\
                         v3rho3,v3rho2sigma,v3rho2lapl,v3rho2tau,v3rhosigma2,v3rhosigmalapl,\
                         v3rhosigmatau,v3rholapl2,v3rholapltau,\
                         v3rhotau2,v3sigma3,v3sigma2lapl,v3sigma2tau,v3sigmalapl2,v3sigmalapltau,\
                         v3sigmatau2,v3lapl3,v3lapl2tau,v3lapltau2,v3tau3,\
                         gam,rt_gam ,rl_gam ,tt_gam ,tl_gam ,ll_gam ,gamx ,gamy ,gamz, \
                         st_gamx,st_gamy,st_gamz,sl_gamx,sl_gamy,sl_gamz, \
                         gamxx, gamxy, gamxz, gamyx, gamyy, gamyz, gamzx, gamzy, gamzz, \
                         rhow12a, gradw12a_x, gradw12a_y, gradw12a_z, tauw12a, laplw12a, \
                         G_val, G_gga_val, G_gga_x_val, G_gga_y_val, G_gga_z_val,\
                         chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                // first-order
                double vsigma_a = vsigma[3 * g + 0];
                double vsigma_c = vsigma[3 * g + 1];

                // second-order
                double v2rho2_aa = v2rho2[3 * g + 0];
                double v2rho2_ab = v2rho2[3 * g + 1];

                double v2rhosigma_aa = v2rhosigma[6 * g + 0];
                double v2rhosigma_ac = v2rhosigma[6 * g + 1];
                double v2rhosigma_ab = v2rhosigma[6 * g + 2];
                double v2rhosigma_ba = v2rhosigma[6 * g + 3];
                double v2rhosigma_bc = v2rhosigma[6 * g + 4];

                double v2sigma2_aa = v2sigma2[6 * g + 0];
                double v2sigma2_ac = v2sigma2[6 * g + 1];
                double v2sigma2_ab = v2sigma2[6 * g + 2];
                double v2sigma2_cc = v2sigma2[6 * g + 3];
                double v2sigma2_cb = v2sigma2[6 * g + 4];

                double v2rholapl_aa = v2rholapl[4 * g + 0];
                double v2rholapl_ab = v2rholapl[4 * g + 1];
                double v2rholapl_ba = v2rholapl[4 * g + 2];

                double v2rhotau_aa = v2rhotau[4 * g + 0];
                double v2rhotau_ab = v2rhotau[4 * g + 1];
                double v2rhotau_ba = v2rhotau[4 * g + 2];

                double v2lapltau_aa = v2lapltau[4 * g + 0];
                double v2lapltau_ab = v2lapltau[4 * g + 1];
                double v2lapltau_ba = v2lapltau[4 * g + 2];

                double v2lapl2_aa = v2lapl2[3 * g + 0];
                double v2lapl2_ab = v2lapl2[3 * g + 1];

                double v2tau2_aa = v2tau2[3 * g + 0];
                double v2tau2_ab = v2tau2[3 * g + 1];

                double v2sigmalapl_aa = v2sigmalapl[6 * g + 0];
                double v2sigmalapl_ab = v2sigmalapl[6 * g + 1];
                double v2sigmalapl_ca = v2sigmalapl[6 * g + 2];
                double v2sigmalapl_cb = v2sigmalapl[6 * g + 3];
                double v2sigmalapl_ba = v2sigmalapl[6 * g + 4];

                double v2sigmatau_aa = v2sigmatau[6 * g + 0];
                double v2sigmatau_ab = v2sigmatau[6 * g + 1];
                double v2sigmatau_ca = v2sigmatau[6 * g + 2];
                double v2sigmatau_cb = v2sigmatau[6 * g + 3];
                double v2sigmatau_ba = v2sigmatau[6 * g + 4];

                // Third-order terms

                auto v3rho3_aaa = v3rho3[4 * g + 0];
                auto v3rho3_aab = v3rho3[4 * g + 1];
                auto v3rho3_abb = v3rho3[4 * g + 2];

                auto v3rho2lapl_aaa = v3rho2lapl[6 * g + 0];
                auto v3rho2lapl_aab = v3rho2lapl[6 * g + 1];
                auto v3rho2lapl_aba = v3rho2lapl[6 * g + 2];
                auto v3rho2lapl_abb = v3rho2lapl[6 * g + 3];
                auto v3rho2lapl_bba = v3rho2lapl[6 * g + 4];

                auto v3rho2tau_aaa = v3rho2tau[6 * g + 0];
                auto v3rho2tau_aab = v3rho2tau[6 * g + 1];
                auto v3rho2tau_aba = v3rho2tau[6 * g + 2];
                auto v3rho2tau_abb = v3rho2tau[6 * g + 3];
                auto v3rho2tau_bba = v3rho2tau[6 * g + 4];

                auto v3rholapl2_aaa = v3rholapl2[6 * g + 0];
                auto v3rholapl2_aab = v3rholapl2[6 * g + 1];
                auto v3rholapl2_abb = v3rholapl2[6 * g + 2];
                auto v3rholapl2_baa = v3rholapl2[6 * g + 3];
                auto v3rholapl2_bab = v3rholapl2[6 * g + 4];

                auto v3rholapltau_aaa = v3rholapltau[8 * g + 0];
                auto v3rholapltau_aab = v3rholapltau[8 * g + 1];
                auto v3rholapltau_aba = v3rholapltau[8 * g + 2];
                auto v3rholapltau_abb = v3rholapltau[8 * g + 3];
                auto v3rholapltau_baa = v3rholapltau[8 * g + 4];
                auto v3rholapltau_bab = v3rholapltau[8 * g + 5];
                auto v3rholapltau_bba = v3rholapltau[8 * g + 6];

                auto v3rhotau2_aaa = v3rhotau2[6 * g + 0];
                auto v3rhotau2_aab = v3rhotau2[6 * g + 1];
                auto v3rhotau2_abb = v3rhotau2[6 * g + 2];
                auto v3rhotau2_baa = v3rhotau2[6 * g + 3];
                auto v3rhotau2_bab = v3rhotau2[6 * g + 4];

                auto v3sigma2lapl_aaa = v3sigma2lapl[12 * g + 0];
                auto v3sigma2lapl_aab = v3sigma2lapl[12 * g + 1];
                auto v3sigma2lapl_aca = v3sigma2lapl[12 * g + 2];
                auto v3sigma2lapl_acb = v3sigma2lapl[12 * g + 3];
                auto v3sigma2lapl_aba = v3sigma2lapl[12 * g + 4];
                auto v3sigma2lapl_abb = v3sigma2lapl[12 * g + 5];
                auto v3sigma2lapl_cca = v3sigma2lapl[12 * g + 6];
                auto v3sigma2lapl_ccb = v3sigma2lapl[12 * g + 7];
                auto v3sigma2lapl_cba = v3sigma2lapl[12 * g + 8];
                auto v3sigma2lapl_cbb = v3sigma2lapl[12 * g + 9];
                auto v3sigma2lapl_bba = v3sigma2lapl[12 * g + 10];

                auto v3rhosigmalapl_aaa = v3rhosigmalapl[12 * g + 0];
                auto v3rhosigmalapl_aab = v3rhosigmalapl[12 * g + 1];
                auto v3rhosigmalapl_aca = v3rhosigmalapl[12 * g + 2];
                auto v3rhosigmalapl_acb = v3rhosigmalapl[12 * g + 3];
                auto v3rhosigmalapl_aba = v3rhosigmalapl[12 * g + 4];
                auto v3rhosigmalapl_abb = v3rhosigmalapl[12 * g + 5];
                auto v3rhosigmalapl_baa = v3rhosigmalapl[12 * g + 6];
                auto v3rhosigmalapl_bab = v3rhosigmalapl[12 * g + 7];
                auto v3rhosigmalapl_bca = v3rhosigmalapl[12 * g + 8];
                auto v3rhosigmalapl_bcb = v3rhosigmalapl[12 * g + 9];
                auto v3rhosigmalapl_bba = v3rhosigmalapl[12 * g + 10];

                auto v3rhosigmatau_aaa = v3rhosigmatau[12 * g + 0];
                auto v3rhosigmatau_aab = v3rhosigmatau[12 * g + 1];
                auto v3rhosigmatau_aca = v3rhosigmatau[12 * g + 2];
                auto v3rhosigmatau_acb = v3rhosigmatau[12 * g + 3];
                auto v3rhosigmatau_aba = v3rhosigmatau[12 * g + 4];
                auto v3rhosigmatau_abb = v3rhosigmatau[12 * g + 5];
                auto v3rhosigmatau_baa = v3rhosigmatau[12 * g + 6];
                auto v3rhosigmatau_bab = v3rhosigmatau[12 * g + 7];
                auto v3rhosigmatau_bca = v3rhosigmatau[12 * g + 8];
                auto v3rhosigmatau_bcb = v3rhosigmatau[12 * g + 9];
                auto v3rhosigmatau_bba = v3rhosigmatau[12 * g + 10];

                auto v3sigma2tau_aaa = v3sigma2tau[12 * g + 0];
                auto v3sigma2tau_aab = v3sigma2tau[12 * g + 1];
                auto v3sigma2tau_aca = v3sigma2tau[12 * g + 2];
                auto v3sigma2tau_acb = v3sigma2tau[12 * g + 3];
                auto v3sigma2tau_aba = v3sigma2tau[12 * g + 4];
                auto v3sigma2tau_abb = v3sigma2tau[12 * g + 5];
                auto v3sigma2tau_cca = v3sigma2tau[12 * g + 6];
                auto v3sigma2tau_ccb = v3sigma2tau[12 * g + 7];
                auto v3sigma2tau_cba = v3sigma2tau[12 * g + 8];
                auto v3sigma2tau_cbb = v3sigma2tau[12 * g + 9];
                auto v3sigma2tau_bba = v3sigma2tau[12 * g + 10];

                auto v3sigmalapl2_aaa = v3sigmalapl2[9 * g + 0];
                auto v3sigmalapl2_aab = v3sigmalapl2[9 * g + 1];
                auto v3sigmalapl2_abb = v3sigmalapl2[9 * g + 2];
                auto v3sigmalapl2_caa = v3sigmalapl2[9 * g + 3];
                auto v3sigmalapl2_cab = v3sigmalapl2[9 * g + 4];
                auto v3sigmalapl2_cbb = v3sigmalapl2[9 * g + 5];
                auto v3sigmalapl2_baa = v3sigmalapl2[9 * g + 6];
                auto v3sigmalapl2_bab = v3sigmalapl2[9 * g + 7];

                auto v3sigmalapltau_aaa = v3sigmalapltau[12 * g + 0];
                auto v3sigmalapltau_aab = v3sigmalapltau[12 * g + 1];
                auto v3sigmalapltau_aba = v3sigmalapltau[12 * g + 2];
                auto v3sigmalapltau_abb = v3sigmalapltau[12 * g + 3];
                auto v3sigmalapltau_caa = v3sigmalapltau[12 * g + 4];
                auto v3sigmalapltau_cab = v3sigmalapltau[12 * g + 5];
                auto v3sigmalapltau_cba = v3sigmalapltau[12 * g + 6];
                auto v3sigmalapltau_cbb = v3sigmalapltau[12 * g + 7];
                auto v3sigmalapltau_baa = v3sigmalapltau[12 * g + 8];
                auto v3sigmalapltau_bab = v3sigmalapltau[12 * g + 9];
                auto v3sigmalapltau_bba = v3sigmalapltau[12 * g + 10];

                auto v3sigmatau2_aaa = v3sigmatau2[9 * g + 0];
                auto v3sigmatau2_aab = v3sigmatau2[9 * g + 1];
                auto v3sigmatau2_abb = v3sigmatau2[9 * g + 2];
                auto v3sigmatau2_caa = v3sigmatau2[9 * g + 3];
                auto v3sigmatau2_cab = v3sigmatau2[9 * g + 4];
                auto v3sigmatau2_cbb = v3sigmatau2[9 * g + 5];
                auto v3sigmatau2_baa = v3sigmatau2[9 * g + 6];
                auto v3sigmatau2_bab = v3sigmatau2[9 * g + 7];

                auto v3lapl3_aaa = v3lapl3[4 * g + 0];
                auto v3lapl3_aab = v3lapl3[4 * g + 1];
                auto v3lapl3_abb = v3lapl3[4 * g + 2];

                auto v3lapl2tau_aaa = v3lapl2tau[6 * g + 0];
                auto v3lapl2tau_aab = v3lapl2tau[6 * g + 1];
                auto v3lapl2tau_aba = v3lapl2tau[6 * g + 2];
                auto v3lapl2tau_abb = v3lapl2tau[6 * g + 3];
                auto v3lapl2tau_bba = v3lapl2tau[6 * g + 4];

                auto v3lapltau2_aaa = v3lapltau2[6 * g + 0];
                auto v3lapltau2_aab = v3lapltau2[6 * g + 1];
                auto v3lapltau2_abb = v3lapltau2[6 * g + 2];
                auto v3lapltau2_baa = v3lapltau2[6 * g + 3];
                auto v3lapltau2_bab = v3lapltau2[6 * g + 4];

                auto v3tau3_aaa = v3tau3[4 * g + 0];
                auto v3tau3_aab = v3tau3[4 * g + 1];
                auto v3tau3_abb = v3tau3[4 * g + 2];
                auto v3tau3_bbb = v3tau3[4 * g + 3];

                auto v3rho2sigma_aaa = v3rho2sigma[9 * g + 0];
                auto v3rho2sigma_aac = v3rho2sigma[9 * g + 1];
                auto v3rho2sigma_aab = v3rho2sigma[9 * g + 2];
                auto v3rho2sigma_aba = v3rho2sigma[9 * g + 3];
                auto v3rho2sigma_abc = v3rho2sigma[9 * g + 4];
                auto v3rho2sigma_abb = v3rho2sigma[9 * g + 5];
                auto v3rho2sigma_bba = v3rho2sigma[9 * g + 6];
                auto v3rho2sigma_bbc = v3rho2sigma[9 * g + 7];

                auto v3rhosigma2_aaa = v3rhosigma2[12 * g + 0];
                auto v3rhosigma2_aac = v3rhosigma2[12 * g + 1];
                auto v3rhosigma2_aab = v3rhosigma2[12 * g + 2];
                auto v3rhosigma2_acc = v3rhosigma2[12 * g + 3];
                auto v3rhosigma2_acb = v3rhosigma2[12 * g + 4];
                auto v3rhosigma2_abb = v3rhosigma2[12 * g + 5];
                auto v3rhosigma2_baa = v3rhosigma2[12 * g + 6];
                auto v3rhosigma2_bac = v3rhosigma2[12 * g + 7];
                auto v3rhosigma2_bab = v3rhosigma2[12 * g + 8];
                auto v3rhosigma2_bcc = v3rhosigma2[12 * g + 9];
                auto v3rhosigma2_bcb = v3rhosigma2[12 * g + 10];

                auto v3sigma3_aaa = v3sigma3[10 * g + 0];
                auto v3sigma3_aac = v3sigma3[10 * g + 1];
                auto v3sigma3_aab = v3sigma3[10 * g + 2];
                auto v3sigma3_acc = v3sigma3[10 * g + 3];
                auto v3sigma3_acb = v3sigma3[10 * g + 4];
                auto v3sigma3_abb = v3sigma3[10 * g + 5];
                auto v3sigma3_ccc = v3sigma3[10 * g + 6];
                auto v3sigma3_ccb = v3sigma3[10 * g + 7];
                auto v3sigma3_cbb = v3sigma3[10 * g + 8];

                // functional derivatives

                // first-order
                double x = vsigma_c + 2.0 * vsigma_a;

                // second-order
                // rho
                double rr = v2rho2_aa + v2rho2_ab;
                double rx = 2.0*v2rhosigma_ac + 2.0*v2rhosigma_ab + 2.0*v2rhosigma_aa;
                double rt = v2rhotau_aa + v2rhotau_ab;
                double rl = v2rholapl_aa + v2rholapl_ab;

                // sigma and gamma
                double xr = v2rhosigma_bc + 2.0*v2rhosigma_ba + v2rhosigma_ac + 2.0 * v2rhosigma_aa;
                double xt = v2sigmatau_cb + 2.0*v2sigmatau_ab + v2sigmatau_ca + 2.0 * v2sigmatau_aa;
                double xl = v2sigmalapl_cb + 2.0*v2sigmalapl_ab + v2sigmalapl_ca + 2.0 * v2sigmalapl_aa;
                double xx = 2.0*v2sigma2_cc + 2.0*v2sigma2_cb + 6.0*v2sigma2_ac + 4.0*v2sigma2_ab + 4.0 * v2sigma2_aa;

                // tau
                double tt = v2tau2_aa + v2tau2_ab;
                double tx = 2.0 * v2sigmatau_ca + 2.0 * v2sigmatau_ba + 2.0 * v2sigmatau_aa;
                double tr = v2rhotau_aa + v2rhotau_ba;
                double tl = v2lapltau_aa + v2lapltau_ba;

                // lapl
                //double ll = v2lapl2_aa + v2lapl2_ab;
                //double lx = 2.0 * v2sigmalapl_ca + 2.0 * v2sigmalapl_ba + 2.0 * v2sigmalapl_aa;
                //double lr = v2rholapl_aa + v2rholapl_ba;
                //double lt = v2lapltau_aa + v2lapltau_ab;

                // Third-oder

                // // sigma and gamma
                double xxx  = 4.0 * v3sigma3_ccc + 8.0 * v3sigma3_ccb + 4.0 * v3sigma3_cbb + 16.0 * v3sigma3_acc
                            + 24.0 * v3sigma3_acb + 8.0 * v3sigma3_abb + 20.0 * v3sigma3_aac + 16.0 * v3sigma3_aab + 8.0 * v3sigma3_aaa;

                double xxr  = 2.0 * v3rhosigma2_bcc + 2.0 * v3rhosigma2_bcb + 6.0 * v3rhosigma2_bac
                            + 4.0 * v3rhosigma2_bab + 4.0 * v3rhosigma2_baa + 2.0 * v3rhosigma2_acc
                            + 2.0 * v3rhosigma2_acb + 6.0 * v3rhosigma2_aac + 4.0 * v3rhosigma2_aab + 4.0 * v3rhosigma2_aaa;

                double xrt  = v3rhosigmatau_bcb + v3rhosigmatau_bca
                            + 2.0 * v3rhosigmatau_bab + 2.0 * v3rhosigmatau_baa
                            + 2.0 * v3rhosigmatau_aab + 2.0 * v3rhosigmatau_aaa
                            + v3rhosigmatau_acb + v3rhosigmatau_aca;

                double xxl  = 2.0 * v3sigma2lapl_ccb + 2.0 * v3sigma2lapl_cca + 2.0 * v3sigma2lapl_cbb
                            + 2.0 * v3sigma2lapl_cba + 6.0 * v3sigma2lapl_acb + 6.0 * v3sigma2lapl_aca
                            + 4.0 * v3sigma2lapl_abb + 4.0 * v3sigma2lapl_aba + 4.0 * v3sigma2lapl_aab + 4.0 * v3sigma2lapl_aaa;
                double xrr  = v3rho2sigma_bbc + 2.0 * v3rho2sigma_bba
                            + 2.0 * v3rho2sigma_abc + 4.0 * v3rho2sigma_aba
                            + v3rho2sigma_aac + 2.0 * v3rho2sigma_aaa;

                double xrl  = v3rhosigmalapl_bcb + v3rhosigmalapl_bca + 2.0 * v3rhosigmalapl_bab
                            + 2.0 * v3rhosigmalapl_baa + v3rhosigmalapl_acb + v3rhosigmalapl_aca
                            + 2.0 * v3rhosigmalapl_aab + 2.0 * v3rhosigmalapl_aaa;

                double xtl  = v3sigmalapltau_cbb + v3sigmalapltau_cba + v3sigmalapltau_cab
                            + v3sigmalapltau_caa + 2.0 * v3sigmalapltau_abb + 2.0 * v3sigmalapltau_aba
                            + 2.0 * v3sigmalapltau_aab + 2.0 * v3sigmalapltau_aaa;
                double xll  = v3sigmalapl2_cbb + 2.0 * v3sigmalapl2_cab + v3sigmalapl2_caa
                            + 2.0 * v3sigmalapl2_abb + 4.0 * v3sigmalapl2_aab + 2.0 * v3sigmalapl2_aaa;


                double xxt  = 2.0 * v3sigma2tau_ccb + 2.0 * v3sigma2tau_cca + 2.0 * v3sigma2tau_cbb
                            + 2.0 * v3sigma2tau_cba + 6.0 * v3sigma2tau_acb + 6.0 * v3sigma2tau_aca
                            + 4.0 * v3sigma2tau_abb + 4.0 * v3sigma2tau_aba + 4.0 * v3sigma2tau_aab + 4.0 * v3sigma2tau_aaa;

                double xtt  = v3sigmatau2_cbb + 2.0 * v3sigmatau2_cab + v3sigmatau2_caa
                            + 2.0 * v3sigmatau2_abb + 4.0 * v3sigmatau2_aab + 2.0 * v3sigmatau2_aaa;

                // rho
                double rrr  = v3rho3_abb + 2.0 * v3rho3_aab + v3rho3_aaa;
                double rrt  = v3rho2tau_abb + v3rho2tau_aba + v3rho2tau_aab + v3rho2tau_aaa;
                double rtx  = 2.0 * v3rhosigmatau_acb + 2.0 * v3rhosigmatau_aca
                            + 2.0 * v3rhosigmatau_abb + 2.0 * v3rhosigmatau_aba
                            + 2.0 * v3rhosigmatau_aab + 2.0 * v3rhosigmatau_aaa;
                double rrl  = v3rho2lapl_abb + v3rho2lapl_aba + v3rho2lapl_aab + v3rho2lapl_aaa;
                double rrx  = 2.0 * v3rho2sigma_abc + 2.0 * v3rho2sigma_abb
                            + 2.0 * v3rho2sigma_aba + 2.0 * v3rho2sigma_aac
                            + 2.0 * v3rho2sigma_aab + 2.0 * v3rho2sigma_aaa;
                double rtl  = v3rholapltau_abb + v3rholapltau_aba + v3rholapltau_aab + v3rholapltau_aaa;

                double rtt  = v3rhotau2_abb + 2.0 * v3rhotau2_aab + v3rhotau2_aaa;

                double rll  = v3rholapl2_abb + 2.0 * v3rholapl2_aab + v3rholapl2_aaa;
                double rlx  = 2.0 * v3rhosigmalapl_acb + 2.0 * v3rhosigmalapl_aca + 2.0 * v3rhosigmalapl_abb
                            + 2.0 * v3rhosigmalapl_aba + 2.0 * v3rhosigmalapl_aab + 2.0 * v3rhosigmalapl_aaa;
                double rxx  = 4.0 * v3rhosigma2_acc + 8.0 * v3rhosigma2_acb + 4.0 * v3rhosigma2_abb
                            + 8.0 * v3rhosigma2_aac + 8.0 * v3rhosigma2_aab + 4.0 * v3rhosigma2_aaa;

                // laplacian
                //double lll  = v3lapl3_abb + 2.0 * v3lapl3_aab + v3lapl3_aaa;
                //double llr  = v3rholapl2_bab + v3rholapl2_baa + v3rholapl2_aab + v3rholapl2_aaa;
                //double llt  = v3lapl2tau_abb + v3lapl2tau_aba + v3lapl2tau_aab + v3lapl2tau_aaa;
                //double llx  = 2.0 * v3sigmalapl2_cab + 2.0 * v3sigmalapl2_caa + 2.0 * v3sigmalapl2_bab
                //            + 2.0 * v3sigmalapl2_baa + 2.0 * v3sigmalapl2_aab + 2.0 * v3sigmalapl2_aaa;
                //double lrr  = v3rho2lapl_bba + 2.0 * v3rho2lapl_aba + v3rho2lapl_aaa;
                //double lrt  = v3rholapltau_bab + v3rholapltau_baa + v3rholapltau_aab + v3rholapltau_aaa;
                //double lrx  = 2.0 * v3rhosigmalapl_bca + 2.0 * v3rhosigmalapl_bba + 2.0 * v3rhosigmalapl_baa
                //            + 2.0 * v3rhosigmalapl_aca + 2.0 * v3rhosigmalapl_aba + 2.0 * v3rhosigmalapl_aaa;
                //double ltt  = v3lapltau2_abb + 2.0 * v3lapltau2_aab + v3lapltau2_aaa;
                //double ltx  = 2.0 * v3sigmalapltau_cab + 2.0 * v3sigmalapltau_caa + 2.0 * v3sigmalapltau_bab
                //            + 2.0 * v3sigmalapltau_baa + 2.0 * v3sigmalapltau_aab + 2.0 * v3sigmalapltau_aaa;
                //double lxx  = 4.0 * v3sigma2lapl_cca + 8.0 * v3sigma2lapl_cba + 4.0 * v3sigma2lapl_bba
                //            + 8.0 * v3sigma2lapl_aca + 8.0 * v3sigma2lapl_aba + 4.0 * v3sigma2lapl_aaa;

                // tau
                double trr  = v3rho2tau_bba + 2.0 * v3rho2tau_aba + v3rho2tau_aaa;
                double ttl  = v3lapltau2_bab + v3lapltau2_baa + v3lapltau2_aab + v3lapltau2_aaa;

                double trl  = v3rholapltau_bba + v3rholapltau_baa + v3rholapltau_aba + v3rholapltau_aaa;

                double tll  = v3lapl2tau_bba + 2.0 * v3lapl2tau_aba + v3lapl2tau_aaa;
                double tlx  = 2.0 * v3sigmalapltau_cba + 2.0 * v3sigmalapltau_caa + 2.0 * v3sigmalapltau_bba
                            + 2.0 * v3sigmalapltau_baa + 2.0 * v3sigmalapltau_aba + 2.0 * v3sigmalapltau_aaa;

                double ttt  = v3tau3_abb + 2.0 * v3tau3_aab + v3tau3_aaa;

                double ttr  = v3rhotau2_bab + v3rhotau2_baa + v3rhotau2_aab + v3rhotau2_aaa;

                double trx  = 2.0 * v3rhosigmatau_bca + 2.0 * v3rhosigmatau_bba
                            + 2.0 * v3rhosigmatau_baa + 2.0 * v3rhosigmatau_aca
                            + 2.0 * v3rhosigmatau_aba + 2.0 * v3rhosigmatau_aaa;

                double txx  = 4.0 * v3sigma2tau_cca + 8.0 * v3sigma2tau_cba + 4.0 * v3sigma2tau_bba
                            + 8.0 * v3sigma2tau_aca + 8.0 * v3sigma2tau_aba + 4.0 * v3sigma2tau_aaa;

                double ttx  = 2.0 * v3sigmatau2_cab + 2.0 * v3sigmatau2_caa + 2.0 * v3sigmatau2_bab
                            + 2.0 * v3sigmatau2_baa + 2.0 * v3sigmatau2_aab + 2.0 * v3sigmatau2_aaa;

                double w = weights[g];

                double rxw12a = gradw12a_x[g];
                double ryw12a = gradw12a_y[g];
                double rzw12a = gradw12a_z[g];

                double grada_x_g = rhograd[6 * g + 0];
                double grada_y_g = rhograd[6 * g + 1];
                double grada_z_g = rhograd[6 * g + 2];

                double l2contract =   grada_x_g * rxw12a + grada_y_g * ryw12a + grada_z_g * rzw12a;

                double l5contract_x = grada_x_g * l2contract;
                double l5contract_y = grada_y_g * l2contract;
                double l5contract_z = grada_z_g * l2contract;

                double q2contract = grada_x_g * gamx[g] + grada_y_g * gamy[g] + grada_z_g * gamz[g];
                double sl_q2contract = grada_x_g * sl_gamx[g] + grada_y_g * sl_gamy[g] + grada_z_g * sl_gamz[g];
                double st_q2contract = grada_x_g * st_gamx[g] + grada_y_g * st_gamy[g] + grada_z_g * st_gamz[g];

                double q3contract =   grada_x_g * grada_x_g * gamxx[g]
                                    + grada_x_g * grada_y_g * gamxy[g]
                                    + grada_x_g * grada_z_g * gamxz[g]
                                    + grada_y_g * grada_x_g * gamyx[g]
                                    + grada_y_g * grada_y_g * gamyy[g]
                                    + grada_y_g * grada_z_g * gamyz[g]
                                    + grada_z_g * grada_x_g * gamzx[g]
                                    + grada_z_g * grada_y_g * gamzy[g]
                                    + grada_z_g * grada_z_g * gamzz[g];

                double q4contract = gamxx[g] + gamyy[g] + gamzz[g];

                double q7contract_x =  grada_x_g * grada_x_g *  gamx[g] + grada_x_g * grada_y_g *  gamy[g] + grada_x_g * grada_z_g *  gamz[g];
                double q7contract_y =  grada_y_g * grada_x_g *  gamx[g] + grada_y_g * grada_y_g *  gamy[g] + grada_y_g * grada_z_g *  gamz[g];
                double q7contract_z =  grada_z_g * grada_x_g *  gamx[g] + grada_z_g * grada_y_g *  gamy[g] + grada_z_g * grada_z_g *  gamz[g];

                double sl_q7contract_x =  grada_x_g * grada_x_g * sl_gamx[g] + grada_x_g * grada_y_g * sl_gamy[g] + grada_x_g * grada_z_g * sl_gamz[g];
                double sl_q7contract_y =  grada_y_g * grada_x_g * sl_gamx[g] + grada_y_g * grada_y_g * sl_gamy[g] + grada_y_g * grada_z_g * sl_gamz[g];
                double sl_q7contract_z =  grada_z_g * grada_x_g * sl_gamx[g] + grada_z_g * grada_y_g * sl_gamy[g] + grada_z_g * grada_z_g * sl_gamz[g];

                double st_q7contract_x =  grada_x_g * grada_x_g * st_gamx[g] + grada_x_g * grada_y_g * st_gamy[g] + grada_x_g * grada_z_g * st_gamz[g];
                double st_q7contract_y =  grada_y_g * grada_x_g * st_gamx[g] + grada_y_g * grada_y_g * st_gamy[g] + grada_y_g * grada_z_g * st_gamz[g];
                double st_q7contract_z =  grada_z_g * grada_x_g * st_gamx[g] + grada_z_g * grada_y_g * st_gamy[g] + grada_z_g * grada_z_g * st_gamz[g];

                double q8contract_x =  grada_x_g *  gamxx[g] + grada_y_g *  gamxy[g] + grada_z_g *  gamxz[g];
                double q8contract_y =  grada_x_g *  gamyx[g] + grada_y_g *  gamyy[g] + grada_z_g *  gamyz[g];
                double q8contract_z =  grada_x_g *  gamzx[g] + grada_y_g *  gamzy[g] + grada_z_g *  gamzz[g];

                double q9contract_x =  grada_x_g *  q3contract;
                double q9contract_y =  grada_y_g *  q3contract;
                double q9contract_z =  grada_z_g *  q3contract;

                double q10contract_x =  grada_x_g *  gamxx[g] + grada_y_g *  gamyx[g] + grada_z_g *  gamzx[g];
                double q10contract_y =  grada_x_g *  gamxy[g] + grada_y_g *  gamyy[g] + grada_z_g *  gamzy[g];
                double q10contract_z =  grada_x_g *  gamxz[g] + grada_y_g *  gamyz[g] + grada_z_g *  gamzz[g];

                double q11contract_x =  grada_x_g *  gamxx[g] + grada_x_g *  gamyy[g] + grada_x_g *  gamzz[g];
                double q11contract_y =  grada_y_g *  gamxx[g] + grada_y_g *  gamyy[g] + grada_y_g *  gamzz[g];
                double q11contract_z =  grada_z_g *  gamxx[g] + grada_z_g *  gamyy[g] + grada_z_g *  gamzz[g];

                // Rho operator contributions

                // vxc 1 contributions

                double rho_0 =  rr * rhow12a[g]
                              + rx * l2contract
                              + rt * tauw12a[g]
                              + rl * laplw12a[g];

                //double lap_0 =    lr * rhow12a[g]
                //                + lx * l2contract
                //                + lt * tauw12a[g]
                //                + ll * laplw12a[g];

                double tau_0 =     tr * rhow12a[g]
                                 + tx * l2contract
                                 + tt * tauw12a[g]
                                 + tl * laplw12a[g];

                // vxc 2 contributions

                rho_0 += rrr * gam[g]
                       + rrt * rt_gam[g]
                       + rrl * rl_gam[g]
                       + rll * ll_gam[g]
                       + rtt * tt_gam[g]
                       + rtl * tl_gam[g]
                       + rrx * q2contract
                       + rlx * sl_q2contract
                       + rtx * st_q2contract
                       + rxx * q3contract
                       + rx  * q4contract;

                //lap_0 += lrr * gam[g]
                //       + lrt * rt_gam[g]
                //       + llr * rl_gam[g]
                //       + lll * ll_gam[g]
                //       + ltt * tt_gam[g]
                //       + llt * tl_gam[g]
                //       + lrx * q2contract
                //       + llx * sl_q2contract
                //       + ltx * st_q2contract
                //       + lxx * q3contract
                //       + lx  * q4contract;

                tau_0 += trr * gam[g]
                       + ttr * rt_gam[g]
                       + trl * rl_gam[g]
                       + tll * ll_gam[g]
                       + ttt * tt_gam[g]
                       + ttl * tl_gam[g]
                       + trx * q2contract
                       + tlx * sl_q2contract
                       + ttx * st_q2contract
                       + txx * q3contract
                       + tx  * q4contract;

                // Grad operator contributions

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                // vxc 1 contributions

                xcomp +=  grada_x_g * ( xr * rhow12a[g] + xt * tauw12a[g] + xl * laplw12a[g])
                        + x * rxw12a
                        + xx * l5contract_x;

                ycomp += grada_y_g * ( xr * rhow12a[g] + xt * tauw12a[g] + xl * laplw12a[g])
                        + x * ryw12a
                        + xx * l5contract_y;

                zcomp += grada_z_g * ( xr * rhow12a[g] + xt * tauw12a[g] + xl * laplw12a[g])
                        + x * rzw12a
                        + xx * l5contract_z;

                // vxc 2 contributions

                xcomp +=  xrr * grada_x_g * gam[g]
                        + xrt * grada_x_g * rt_gam[g]
                        + xrl * grada_x_g * rl_gam[g]
                        + xll * grada_x_g * ll_gam[g]
                        + xtt * grada_x_g * tt_gam[g]
                        + xtl * grada_x_g * tl_gam[g]
                        + xr * gamx[g] // q6
                        + xl * sl_gamx[g]
                        + xt * st_gamx[g]
                        + xxr * q7contract_x
                        + xxl * sl_q7contract_x
                        + xxt * st_q7contract_x
                        + xx * (q8contract_x + q10contract_x + q11contract_x)
                        + xxx * q9contract_x;

                ycomp +=  xrr * grada_y_g * gam[g] // q5
                        + xrt * grada_y_g * rt_gam[g]
                        + xrl * grada_y_g * rl_gam[g]
                        + xll * grada_y_g * ll_gam[g]
                        + xtt * grada_y_g * tt_gam[g]
                        + xtl * grada_y_g * tl_gam[g]
                        + xr * gamy[g] // q6
                        + xl * sl_gamy[g]
                        + xt * st_gamy[g]
                        + xxr * q7contract_y
                        + xxl * sl_q7contract_y
                        + xxt * st_q7contract_y
                        + xx * (q8contract_y + q10contract_y + q11contract_y)
                        + xxx * q9contract_y;

                zcomp +=  xrr * grada_z_g * gam[g] // q5
                        + xrt * grada_z_g * rt_gam[g]
                        + xrl * grada_z_g * rl_gam[g]
                        + xll * grada_z_g * ll_gam[g]
                        + xtt * grada_z_g * tt_gam[g]
                        + xtl * grada_z_g * tl_gam[g]
                        + xr * gamz[g] // q6
                        + xl * sl_gamz[g]
                        + xt * st_gamz[g]
                        + xxr * q7contract_z
                        + xxl * sl_q7contract_z
                        + xxt * st_q7contract_z
                        + xx * (q8contract_z + q10contract_z + q11contract_z)
                        + xxx * q9contract_z;

                G_val[nu_offset + g] = w * rho_0 * chi_val[nu_offset + g];

                G_gga_val[nu_offset + g] = w * (xcomp * chi_x_val[nu_offset + g] +
                                                ycomp * chi_y_val[nu_offset + g] +
                                                zcomp * chi_z_val[nu_offset + g]);

                // TODO implement Laplacian dependence

                //G_gga_val[nu_offset + g] += w * lap_0 * (chi_xx_val[nu_offset + g] +
                //                                         chi_yy_val[nu_offset + g] +
                //                                         chi_zz_val[nu_offset + g]);

                G_gga_x_val[nu_offset + g] = w * tau_0 * chi_x_val[nu_offset + g];
                G_gga_y_val[nu_offset + g] = w * tau_0 * chi_y_val[nu_offset + g];
                G_gga_z_val[nu_offset + g] = w * tau_0 * chi_z_val[nu_offset + g];
            }
        }
    }

    timer.stop("Kxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Kxc matrix matmul");

    // LDA and GGA contribution
    auto mat_Kxc = denblas::multABt(gtoValues, denblas::addAB(mat_G, mat_G_gga, 2.0));

    // tau contribution
    auto mat_Kxc_x = denblas::multABt(gtoValuesX, mat_G_gga_x);
    auto mat_Kxc_y = denblas::multABt(gtoValuesY, mat_G_gga_y);
    auto mat_Kxc_z = denblas::multABt(gtoValuesZ, mat_G_gga_z);

    mat_Kxc = denblas::addAB(mat_Kxc, mat_Kxc_x, 0.5);
    mat_Kxc = denblas::addAB(mat_Kxc, mat_Kxc_y, 0.5);
    mat_Kxc = denblas::addAB(mat_Kxc, mat_Kxc_z, 0.5);

    mat_Kxc.symmetrizeAndScale(0.5);

    timer.stop("Kxc matrix matmul");

    return mat_Kxc;
}

CDenseMatrix
CXCNewIntegrator::_integratePartialKxcFockForLDA2(const int32_t            npoints,
                                                  const double*            weights,
                                                  const CDenseMatrix&      gtoValues,
                                                  const double*            v2rho2,
                                                  const double*            v3rho3,
                                                  const CDensityGridCubic& rwDensityGridCubic,
                                                  const CDensityGrid&      rw2DensityGrid,
                                                  const int32_t            iFock,
                                                  CMultiTimer&             timer) const
{
    timer.start("Kxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    // pointers to perturbed density

    auto rhow1a = rwDensityGridCubic.gam2(iFock);

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

            #pragma omp simd aligned(weights, v2rho2, v3rho3, rhow1a, rhow12a, rhow12b, G_val, chi_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                G_val[nu_offset + g] = weights[g] *

                          ((v3rho3[4 * g + 0] + 2.0 * v3rho3[4 * g + 1] + v3rho3[4 * g + 2]) * rhow1a[g] +

                           v2rho2[3 * g + 0] * rhow12a[g] + v2rho2[3 * g + 1] * rhow12b[g]) *

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
CXCNewIntegrator::_integratePartialKxcFockForGGA2(const int32_t            npoints,
                                                  const double*            weights,
                                                  const CDenseMatrix&      gtoValues,
                                                  const CDenseMatrix&      gtoValuesX,
                                                  const CDenseMatrix&      gtoValuesY,
                                                  const CDenseMatrix&      gtoValuesZ,
                                                  const double*            rhograd,
                                                  const double*            vsigma,
                                                  const double*            v2rho2,
                                                  const double*            v2rhosigma,
                                                  const double*            v2sigma2,
                                                  const double*            v3rho3,
                                                  const double*            v3rho2sigma,
                                                  const double*            v3rhosigma2,
                                                  const double*            v3sigma3,
                                                  const CDensityGridCubic& rwDensityGridCubic,
                                                  const CDensityGrid&      rw2DensityGrid,
                                                  const int32_t            iFock,
                                                  CMultiTimer&             timer) const
{
    timer.start("Kxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    // pointers to perturbed density gradient norms

    auto rhow1rhow2 = rwDensityGridCubic.gam2(iFock);

    auto rxw1rhow2 = rwDensityGridCubic.gam2X(iFock);

    auto ryw1rhow2 = rwDensityGridCubic.gam2Y(iFock);

    auto rzw1rhow2 = rwDensityGridCubic.gam2Z(iFock);

    auto rxw1rxw2 = rwDensityGridCubic.gam2XX(iFock);

    auto rxw1ryw2 = rwDensityGridCubic.gam2XY(iFock);

    auto rxw1rzw2 = rwDensityGridCubic.gam2XZ(iFock);

    auto ryw1rxw2 = rwDensityGridCubic.gam2YX(iFock);

    auto ryw1ryw2 = rwDensityGridCubic.gam2YY(iFock);

    auto ryw1rzw2 = rwDensityGridCubic.gam2YZ(iFock);

    auto rzw1rxw2 = rwDensityGridCubic.gam2ZX(iFock);

    auto rzw1ryw2 = rwDensityGridCubic.gam2ZY(iFock);

    auto rzw1rzw2 = rwDensityGridCubic.gam2ZZ(iFock);

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
                    rhograd, vsigma, v2rho2, v2rhosigma, v2sigma2, \
                    v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3, \
                    rhow1rhow2, rxw1rhow2, ryw1rhow2, rzw1rhow2, \
                    rxw1rxw2, rxw1ryw2, rxw1rzw2, ryw1rxw2, ryw1ryw2, ryw1rzw2, rzw1rxw2, rzw1ryw2, rzw1rzw2, \
                    rhow12a, gradw12a_x, gradw12a_y, gradw12a_z, \
                    G_val, G_gga_val, chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {

                double w = weights[g];


                double rxw12a = gradw12a_x[g];
                double ryw12a = gradw12a_y[g];
                double rzw12a = gradw12a_z[g];

                double grada_x_g = rhograd[6 * g + 0];
                double grada_y_g = rhograd[6 * g + 1];
                double grada_z_g = rhograd[6 * g + 2];


                double l2contract = grada_x_g * rxw12a + grada_y_g * ryw12a + grada_z_g * rzw12a;
                double l5contract_x = grada_x_g * l2contract;
                double l5contract_y = grada_y_g * l2contract;
                double l5contract_z = grada_z_g * l2contract;
                double q2contract = grada_x_g * rxw1rhow2[g] + grada_y_g * ryw1rhow2[g] + grada_z_g * rzw1rhow2[g];
                double q3contract =   grada_x_g * grada_x_g * rxw1rxw2[g]
                                    + grada_x_g * grada_y_g * rxw1ryw2[g]
                                    + grada_x_g * grada_z_g * rxw1rzw2[g]
                                    + grada_y_g * grada_x_g * ryw1rxw2[g]
                                    + grada_y_g * grada_y_g * ryw1ryw2[g]
                                    + grada_y_g * grada_z_g * ryw1rzw2[g]
                                    + grada_z_g * grada_x_g * rzw1rxw2[g]
                                    + grada_z_g * grada_y_g * rzw1ryw2[g]
                                    + grada_z_g * grada_z_g * rzw1rzw2[g];

                double q4contract = rxw1rxw2[g] + ryw1ryw2[g] + rzw1rzw2[g];
                double q7contract_x =  grada_x_g * grada_x_g *  rxw1rhow2[g] + grada_x_g * grada_y_g *  ryw1rhow2[g] + grada_x_g * grada_z_g *  rzw1rhow2[g];
                double q7contract_y =  grada_y_g * grada_x_g *  rxw1rhow2[g] + grada_y_g * grada_y_g *  ryw1rhow2[g] + grada_y_g * grada_z_g *  rzw1rhow2[g];
                double q7contract_z =  grada_z_g * grada_x_g *  rxw1rhow2[g] + grada_z_g * grada_y_g *  ryw1rhow2[g] + grada_z_g * grada_z_g *  rzw1rhow2[g];


                double q8contract_x =  grada_x_g *  rxw1rxw2[g] + grada_y_g *  rxw1ryw2[g] + grada_z_g *  rxw1rzw2[g];
                double q8contract_y =  grada_x_g *  ryw1rxw2[g] + grada_y_g *  ryw1ryw2[g] + grada_z_g *  ryw1rzw2[g];
                double q8contract_z =  grada_x_g *  rzw1rxw2[g] + grada_y_g *  rzw1ryw2[g] + grada_z_g *  rzw1rzw2[g];

                double q9contract_x =  grada_x_g *  q3contract;
                double q9contract_y =  grada_y_g *  q3contract;
                double q9contract_z =  grada_z_g *  q3contract;

                double q10contract_x =  grada_x_g *  rxw1rxw2[g] + grada_y_g *  ryw1rxw2[g] + grada_z_g *  rzw1rxw2[g];
                double q10contract_y =  grada_x_g *  rxw1ryw2[g] + grada_y_g *  ryw1ryw2[g] + grada_z_g *  rzw1ryw2[g];
                double q10contract_z =  grada_x_g *  rxw1rzw2[g] + grada_y_g *  ryw1rzw2[g] + grada_z_g *  rzw1rzw2[g];

                double q11contract_x =  grada_x_g *  rxw1rxw2[g] + grada_x_g *  ryw1ryw2[g] + grada_x_g *  rzw1rzw2[g];
                double q11contract_y =  grada_y_g *  rxw1rxw2[g] + grada_y_g *  ryw1ryw2[g] + grada_y_g *  rzw1rzw2[g];
                double q11contract_z =  grada_z_g *  rxw1rxw2[g] + grada_z_g *  ryw1ryw2[g] + grada_z_g *  rzw1rzw2[g];

                // functional derivatives in libxc form

                auto vsigma_a = vsigma[3 * g + 0];
                auto vsigma_c = vsigma[3 * g + 1];

                auto v2rho2_aa = v2rho2[3 * g + 0];
                auto v2rho2_ab = v2rho2[3 * g + 1];

                auto v2rhosigma_aa = v2rhosigma[6 * g + 0];
                auto v2rhosigma_ac = v2rhosigma[6 * g + 1];
                auto v2rhosigma_ab = v2rhosigma[6 * g + 2];
                auto v2rhosigma_ba = v2rhosigma[6 * g + 3];
                auto v2rhosigma_bc = v2rhosigma[6 * g + 4];

                auto v2sigma2_aa = v2sigma2[6 * g + 0];
                auto v2sigma2_ac = v2sigma2[6 * g + 1];
                auto v2sigma2_ab = v2sigma2[6 * g + 2];
                auto v2sigma2_cc = v2sigma2[6 * g + 3];
                auto v2sigma2_cb = v2sigma2[6 * g + 4];

                auto v3rho3_aaa = v3rho3[4 * g + 0];
                auto v3rho3_aab = v3rho3[4 * g + 1];
                auto v3rho3_abb = v3rho3[4 * g + 2];

                auto v3rho2sigma_aaa = v3rho2sigma[9 * g + 0];
                auto v3rho2sigma_aac = v3rho2sigma[9 * g + 1];
                auto v3rho2sigma_aab = v3rho2sigma[9 * g + 2];
                auto v3rho2sigma_aba = v3rho2sigma[9 * g + 3];
                auto v3rho2sigma_abc = v3rho2sigma[9 * g + 4];
                auto v3rho2sigma_abb = v3rho2sigma[9 * g + 5];
                auto v3rho2sigma_bba = v3rho2sigma[9 * g + 6];
                auto v3rho2sigma_bbc = v3rho2sigma[9 * g + 7];

                auto v3rhosigma2_aaa = v3rhosigma2[12 * g + 0];
                auto v3rhosigma2_aac = v3rhosigma2[12 * g + 1];
                auto v3rhosigma2_aab = v3rhosigma2[12 * g + 2];
                auto v3rhosigma2_acc = v3rhosigma2[12 * g + 3];
                auto v3rhosigma2_acb = v3rhosigma2[12 * g + 4];
                auto v3rhosigma2_abb = v3rhosigma2[12 * g + 5];
                auto v3rhosigma2_baa = v3rhosigma2[12 * g + 6];
                auto v3rhosigma2_bac = v3rhosigma2[12 * g + 7];
                auto v3rhosigma2_bab = v3rhosigma2[12 * g + 8];
                auto v3rhosigma2_bcc = v3rhosigma2[12 * g + 9];
                auto v3rhosigma2_bcb = v3rhosigma2[12 * g + 10];

                auto v3sigma3_aaa = v3sigma3[10 * g + 0];
                auto v3sigma3_aac = v3sigma3[10 * g + 1];
                auto v3sigma3_aab = v3sigma3[10 * g + 2];
                auto v3sigma3_acc = v3sigma3[10 * g + 3];
                auto v3sigma3_acb = v3sigma3[10 * g + 4];
                auto v3sigma3_abb = v3sigma3[10 * g + 5];
                auto v3sigma3_ccc = v3sigma3[10 * g + 6];
                auto v3sigma3_bcc = v3sigma3[10 * g + 7];
                auto v3sigma3_cbb = v3sigma3[10 * g + 8];


                // functional derivatives
                double rr = (v2rho2_aa + v2rho2_ab);
                double rrr = (v3rho3_aaa + 2.0 * v3rho3_aab + v3rho3_abb);
                double rx = (2.0*v2rhosigma_ac + 2.0*v2rhosigma_ab + 2.0*v2rhosigma_aa);
                double rxr = (2.0*v3rho2sigma_abc + 2.0*v3rho2sigma_abb + 2.0*v3rho2sigma_aba
                             + 2.0*v3rho2sigma_aac + 2.0*v3rho2sigma_aab + 2.0*v3rho2sigma_aaa);
                double rxx = (4.0*v3rhosigma2_acc + 8.0*v3rhosigma2_acb + 4.0*v3rhosigma2_abb
                            + 8.0*v3rhosigma2_aac + 8.0*v3rhosigma2_aab + 4.0*v3rhosigma2_aaa);
                double x = vsigma_c + 2.0*vsigma_a;
                double xr = v2rhosigma_bc + 2.0*v2rhosigma_ba + v2rhosigma_ac + 2.0*v2rhosigma_aa;
                double xx = 2.0*v2sigma2_cc + 2.0*v2sigma2_cb + 6.0*v2sigma2_ac + 4.0*v2sigma2_ab + 4.0*v2sigma2_aa;
                double xrr = v3rho2sigma_bbc + 2.0*v3rho2sigma_bba + 2.0*v3rho2sigma_abc + 4.0*v3rho2sigma_aba
                            + v3rho2sigma_aac + 2.0*v3rho2sigma_aaa;
                double xxr = 2.0*v3rhosigma2_bcc + 2.0*v3rhosigma2_bcb + 6.0*v3rhosigma2_bac
                            + 4.0*v3rhosigma2_bab + 4.0*v3rhosigma2_baa + 2.0*v3rhosigma2_acc
                            + 2.0*v3rhosigma2_acb + 6.0*v3rhosigma2_aac + 4.0*v3rhosigma2_aab + 4.0*v3rhosigma2_aaa;
                double xxx = 4.0*v3sigma3_ccc + 8.0*v3sigma3_bcc + 4.0*v3sigma3_cbb + 16.0*v3sigma3_acc
                            + 24.0*v3sigma3_acb + 8.0*v3sigma3_abb + 20.0*v3sigma3_aac
                            + 16.0*v3sigma3_aab + 8.0*v3sigma3_aaa;

                // Scalar contribution

                double prefac = 0.0;

                // vxc 1 contributions

                prefac += rr * rhow12a[g] // l1
                        + rx * l2contract;

                // vxc 2 contributions

                prefac += rrr * rhow1rhow2[g] // q1
                        + rxr * q2contract
                        + rxx * q3contract
                        + rx * q4contract;

                G_val[nu_offset + g] = w * prefac * chi_val[nu_offset + g];

                // vector contribution

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                // vxc 1 contributions

                xcomp += xr * grada_x_g * rhow12a[g] // l3
                        + x * rxw12a // l4
                        + xx * l5contract_x;

                ycomp += xr * grada_y_g * rhow12a[g] // l3
                        + x * ryw12a // l4
                        + xx * l5contract_y;

                zcomp += xr * grada_z_g * rhow12a[g] // l3
                        + x * rzw12a // l4
                        + xx * l5contract_z;

                // vxc 2 contributions

                xcomp += xrr * grada_x_g * rhow1rhow2[g] // q5
                        + xr * rxw1rhow2[g] // q6
                        + xxr * q7contract_x
                        + xx * (q8contract_x + q10contract_x + q11contract_x)
                        + xxx * q9contract_x;

                ycomp += xrr * grada_y_g * rhow1rhow2[g] // q5
                        + xr * ryw1rhow2[g] // q6
                        + xxr * q7contract_y
                        + xx * (q8contract_y + q10contract_y + q11contract_y)
                        + xxx * q9contract_y;

                zcomp += xrr * grada_z_g * rhow1rhow2[g] // q5
                        + xr * rzw1rhow2[g] // q6
                        + xxr * q7contract_z
                        + xx * (q8contract_z + q10contract_z + q11contract_z)
                        + xxx * q9contract_z;

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
CXCNewIntegrator::_integratePartialKxcFockForMGGA2(const int32_t            npoints,
                                                   const double*            weights,
                                                   const CDenseMatrix&      gtoValues,
                                                   const CDenseMatrix&      gtoValuesX,
                                                   const CDenseMatrix&      gtoValuesY,
                                                   const CDenseMatrix&      gtoValuesZ,
                                                   const double*            rhograd,
                                                   const double*            vsigma,
                                                   const double*            v2rho2,
                                                   const double*            v2lapl2,
                                                   const double*            v2tau2,
                                                   const double*            v2rholapl,
                                                   const double*            v2rhotau,
                                                   const double*            v2lapltau,
                                                   const double*            v2rhosigma,
                                                   const double*            v2sigmalapl,
                                                   const double*            v2sigmatau,
                                                   const double*            v2sigma2,
                                                   const double*            v3rho3,
                                                   const double*            v3rho2sigma,
                                                   const double*            v3rho2lapl,
                                                   const double*            v3rho2tau,
                                                   const double*            v3rhosigma2,
                                                   const double*            v3rhosigmalapl,
                                                   const double*            v3rhosigmatau,
                                                   const double*            v3rholapl2,
                                                   const double*            v3rholapltau,
                                                   const double*            v3rhotau2,
                                                   const double*            v3sigma3,
                                                   const double*            v3sigma2lapl,
                                                   const double*            v3sigma2tau,
                                                   const double*            v3sigmalapl2,
                                                   const double*            v3sigmalapltau,
                                                   const double*            v3sigmatau2,
                                                   const double*            v3lapl3,
                                                   const double*            v3lapl2tau,
                                                   const double*            v3lapltau2,
                                                   const double*            v3tau3,
                                                   const CDensityGridCubic& rwDensityGridCubic,
                                                   const CDensityGrid&      rw2DensityGrid,
                                                   const int32_t            iFock,
                                                   CMultiTimer&             timer) const
{
    timer.start("Kxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    // pointers to perturbed densities

    auto gam = rwDensityGridCubic.gam2(iFock);
    auto rt_gam = rwDensityGridCubic.rt_gam2(iFock);
    auto rl_gam = rwDensityGridCubic.rl_gam2(iFock);
    auto tt_gam = rwDensityGridCubic.tt_gam2(iFock);
    auto tl_gam = rwDensityGridCubic.tl_gam2(iFock);
    auto ll_gam = rwDensityGridCubic.ll_gam2(iFock);

    auto gamx = rwDensityGridCubic.gam2X(iFock);
    auto gamy = rwDensityGridCubic.gam2Y(iFock);
    auto gamz = rwDensityGridCubic.gam2Z(iFock);
    auto st_gamx = rwDensityGridCubic.st_gam2X(iFock);
    auto st_gamy = rwDensityGridCubic.st_gam2Y(iFock);
    auto st_gamz = rwDensityGridCubic.st_gam2Z(iFock);
    auto sl_gamx = rwDensityGridCubic.sl_gam2X(iFock);
    auto sl_gamy = rwDensityGridCubic.sl_gam2Y(iFock);
    auto sl_gamz = rwDensityGridCubic.sl_gam2Z(iFock);

    auto gamxx = rwDensityGridCubic.gam2XX(iFock);
    auto gamxy = rwDensityGridCubic.gam2XY(iFock);
    auto gamxz = rwDensityGridCubic.gam2XZ(iFock);
    auto gamyx = rwDensityGridCubic.gam2YX(iFock);
    auto gamyy = rwDensityGridCubic.gam2YY(iFock);
    auto gamyz = rwDensityGridCubic.gam2YZ(iFock);
    auto gamzx = rwDensityGridCubic.gam2ZX(iFock);
    auto gamzy = rwDensityGridCubic.gam2ZY(iFock);
    auto gamzz = rwDensityGridCubic.gam2ZZ(iFock);

    auto rhow = rw2DensityGrid.alphaDensity(iFock);
    auto tauw = rw2DensityGrid.alphaDensitytau(iFock);
    auto laplw = rw2DensityGrid.alphaDensitylapl(iFock);
    auto gradw_x = rw2DensityGrid.alphaDensityGradientX(iFock);
    auto gradw_y = rw2DensityGrid.alphaDensityGradientY(iFock);
    auto gradw_z = rw2DensityGrid.alphaDensityGradientZ(iFock);

    timer.stop("Kxc matrix prep.");

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Kxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    // LDA contribution
    CDenseMatrix mat_G(naos, npoints);

    // GGA contribution
    CDenseMatrix mat_G_gga(naos, npoints);

    // tau contribution
    CDenseMatrix mat_G_gga_x(naos, npoints);
    CDenseMatrix mat_G_gga_y(naos, npoints);
    CDenseMatrix mat_G_gga_z(naos, npoints);

    auto G_val = mat_G.values();

    auto G_gga_val = mat_G_gga.values();

    auto G_gga_x_val = mat_G_gga_x.values();
    auto G_gga_y_val = mat_G_gga_y.values();
    auto G_gga_z_val = mat_G_gga_z.values();

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
                         rhograd, vsigma,\
                         v2rho2,v2lapl2, v2tau2, v2rholapl, v2rhotau,v2lapltau, v2rhosigma,\
                         v2sigmalapl, v2sigmatau,v2sigma2,\
                         v3rho3,v3rho2sigma,v3rho2lapl,v3rho2tau,v3rhosigma2,v3rhosigmalapl,\
                         v3rhosigmatau,v3rholapl2,v3rholapltau,\
                         v3rhotau2,v3sigma3,v3sigma2lapl,v3sigma2tau,v3sigmalapl2,v3sigmalapltau,\
                         v3sigmatau2,v3lapl3,v3lapl2tau,v3lapltau2,v3tau3,\
                         gam,rt_gam ,rl_gam ,tt_gam ,tl_gam ,ll_gam ,gamx ,gamy ,gamz, \
                         st_gamx,st_gamy,st_gamz,sl_gamx,sl_gamy,sl_gamz, \
                         gamxx, gamxy, gamxz, gamyx, gamyy, gamyz, gamzx, gamzy, gamzz, \
                         rhow, gradw_x, gradw_y, gradw_z, tauw, laplw, \
                         G_val, G_gga_val, G_gga_x_val, G_gga_y_val, G_gga_z_val,\
                         chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                // first-order
                double vsigma_a = vsigma[3 * g + 0];
                double vsigma_c = vsigma[3 * g + 1];

                // second-order
                double v2rho2_aa = v2rho2[3 * g + 0];
                double v2rho2_ab = v2rho2[3 * g + 1];

                double v2rhosigma_aa = v2rhosigma[6 * g + 0];
                double v2rhosigma_ac = v2rhosigma[6 * g + 1];
                double v2rhosigma_ab = v2rhosigma[6 * g + 2];
                double v2rhosigma_ba = v2rhosigma[6 * g + 3];
                double v2rhosigma_bc = v2rhosigma[6 * g + 4];

                double v2sigma2_aa = v2sigma2[6 * g + 0];
                double v2sigma2_ac = v2sigma2[6 * g + 1];
                double v2sigma2_ab = v2sigma2[6 * g + 2];
                double v2sigma2_cc = v2sigma2[6 * g + 3];
                double v2sigma2_cb = v2sigma2[6 * g + 4];

                double v2rholapl_aa = v2rholapl[4 * g + 0];
                double v2rholapl_ab = v2rholapl[4 * g + 1];
                double v2rholapl_ba = v2rholapl[4 * g + 2];

                double v2rhotau_aa = v2rhotau[4 * g + 0];
                double v2rhotau_ab = v2rhotau[4 * g + 1];
                double v2rhotau_ba = v2rhotau[4 * g + 2];

                double v2lapltau_aa = v2lapltau[4 * g + 0];
                double v2lapltau_ab = v2lapltau[4 * g + 1];
                double v2lapltau_ba = v2lapltau[4 * g + 2];

                double v2lapl2_aa = v2lapl2[3 * g + 0];
                double v2lapl2_ab = v2lapl2[3 * g + 1];

                double v2tau2_aa = v2tau2[3 * g + 0];
                double v2tau2_ab = v2tau2[3 * g + 1];

                double v2sigmalapl_aa = v2sigmalapl[6 * g + 0];
                double v2sigmalapl_ab = v2sigmalapl[6 * g + 1];
                double v2sigmalapl_ca = v2sigmalapl[6 * g + 2];
                double v2sigmalapl_cb = v2sigmalapl[6 * g + 3];
                double v2sigmalapl_ba = v2sigmalapl[6 * g + 4];

                double v2sigmatau_aa = v2sigmatau[6 * g + 0];
                double v2sigmatau_ab = v2sigmatau[6 * g + 1];
                double v2sigmatau_ca = v2sigmatau[6 * g + 2];
                double v2sigmatau_cb = v2sigmatau[6 * g + 3];
                double v2sigmatau_ba = v2sigmatau[6 * g + 4];

                // Third-order terms

                auto v3rho3_aaa = v3rho3[4 * g + 0];
                auto v3rho3_aab = v3rho3[4 * g + 1];
                auto v3rho3_abb = v3rho3[4 * g + 2];

                auto v3rho2lapl_aaa = v3rho2lapl[6 * g + 0];
                auto v3rho2lapl_aab = v3rho2lapl[6 * g + 1];
                auto v3rho2lapl_aba = v3rho2lapl[6 * g + 2];
                auto v3rho2lapl_abb = v3rho2lapl[6 * g + 3];
                auto v3rho2lapl_bba = v3rho2lapl[6 * g + 4];

                auto v3rho2tau_aaa = v3rho2tau[6 * g + 0];
                auto v3rho2tau_aab = v3rho2tau[6 * g + 1];
                auto v3rho2tau_aba = v3rho2tau[6 * g + 2];
                auto v3rho2tau_abb = v3rho2tau[6 * g + 3];
                auto v3rho2tau_bba = v3rho2tau[6 * g + 4];

                auto v3rholapl2_aaa = v3rholapl2[6 * g + 0];
                auto v3rholapl2_aab = v3rholapl2[6 * g + 1];
                auto v3rholapl2_abb = v3rholapl2[6 * g + 2];
                auto v3rholapl2_baa = v3rholapl2[6 * g + 3];
                auto v3rholapl2_bab = v3rholapl2[6 * g + 4];

                auto v3rholapltau_aaa = v3rholapltau[8 * g + 0];
                auto v3rholapltau_aab = v3rholapltau[8 * g + 1];
                auto v3rholapltau_aba = v3rholapltau[8 * g + 2];
                auto v3rholapltau_abb = v3rholapltau[8 * g + 3];
                auto v3rholapltau_baa = v3rholapltau[8 * g + 4];
                auto v3rholapltau_bab = v3rholapltau[8 * g + 5];
                auto v3rholapltau_bba = v3rholapltau[8 * g + 6];

                auto v3rhotau2_aaa = v3rhotau2[6 * g + 0];
                auto v3rhotau2_aab = v3rhotau2[6 * g + 1];
                auto v3rhotau2_abb = v3rhotau2[6 * g + 2];
                auto v3rhotau2_baa = v3rhotau2[6 * g + 3];
                auto v3rhotau2_bab = v3rhotau2[6 * g + 4];

                auto v3sigma2lapl_aaa = v3sigma2lapl[12 * g + 0];
                auto v3sigma2lapl_aab = v3sigma2lapl[12 * g + 1];
                auto v3sigma2lapl_aca = v3sigma2lapl[12 * g + 2];
                auto v3sigma2lapl_acb = v3sigma2lapl[12 * g + 3];
                auto v3sigma2lapl_aba = v3sigma2lapl[12 * g + 4];
                auto v3sigma2lapl_abb = v3sigma2lapl[12 * g + 5];
                auto v3sigma2lapl_cca = v3sigma2lapl[12 * g + 6];
                auto v3sigma2lapl_ccb = v3sigma2lapl[12 * g + 7];
                auto v3sigma2lapl_cba = v3sigma2lapl[12 * g + 8];
                auto v3sigma2lapl_cbb = v3sigma2lapl[12 * g + 9];
                auto v3sigma2lapl_bba = v3sigma2lapl[12 * g + 10];

                auto v3rhosigmalapl_aaa = v3rhosigmalapl[12 * g + 0];
                auto v3rhosigmalapl_aab = v3rhosigmalapl[12 * g + 1];
                auto v3rhosigmalapl_aca = v3rhosigmalapl[12 * g + 2];
                auto v3rhosigmalapl_acb = v3rhosigmalapl[12 * g + 3];
                auto v3rhosigmalapl_aba = v3rhosigmalapl[12 * g + 4];
                auto v3rhosigmalapl_abb = v3rhosigmalapl[12 * g + 5];
                auto v3rhosigmalapl_baa = v3rhosigmalapl[12 * g + 6];
                auto v3rhosigmalapl_bab = v3rhosigmalapl[12 * g + 7];
                auto v3rhosigmalapl_bca = v3rhosigmalapl[12 * g + 8];
                auto v3rhosigmalapl_bcb = v3rhosigmalapl[12 * g + 9];
                auto v3rhosigmalapl_bba = v3rhosigmalapl[12 * g + 10];

                auto v3rhosigmatau_aaa = v3rhosigmatau[12 * g + 0];
                auto v3rhosigmatau_aab = v3rhosigmatau[12 * g + 1];
                auto v3rhosigmatau_aca = v3rhosigmatau[12 * g + 2];
                auto v3rhosigmatau_acb = v3rhosigmatau[12 * g + 3];
                auto v3rhosigmatau_aba = v3rhosigmatau[12 * g + 4];
                auto v3rhosigmatau_abb = v3rhosigmatau[12 * g + 5];
                auto v3rhosigmatau_baa = v3rhosigmatau[12 * g + 6];
                auto v3rhosigmatau_bab = v3rhosigmatau[12 * g + 7];
                auto v3rhosigmatau_bca = v3rhosigmatau[12 * g + 8];
                auto v3rhosigmatau_bcb = v3rhosigmatau[12 * g + 9];
                auto v3rhosigmatau_bba = v3rhosigmatau[12 * g + 10];

                auto v3sigma2tau_aaa = v3sigma2tau[12 * g + 0];
                auto v3sigma2tau_aab = v3sigma2tau[12 * g + 1];
                auto v3sigma2tau_aca = v3sigma2tau[12 * g + 2];
                auto v3sigma2tau_acb = v3sigma2tau[12 * g + 3];
                auto v3sigma2tau_aba = v3sigma2tau[12 * g + 4];
                auto v3sigma2tau_abb = v3sigma2tau[12 * g + 5];
                auto v3sigma2tau_cca = v3sigma2tau[12 * g + 6];
                auto v3sigma2tau_ccb = v3sigma2tau[12 * g + 7];
                auto v3sigma2tau_cba = v3sigma2tau[12 * g + 8];
                auto v3sigma2tau_cbb = v3sigma2tau[12 * g + 9];
                auto v3sigma2tau_bba = v3sigma2tau[12 * g + 10];

                auto v3sigmalapl2_aaa = v3sigmalapl2[9 * g + 0];
                auto v3sigmalapl2_aab = v3sigmalapl2[9 * g + 1];
                auto v3sigmalapl2_abb = v3sigmalapl2[9 * g + 2];
                auto v3sigmalapl2_caa = v3sigmalapl2[9 * g + 3];
                auto v3sigmalapl2_cab = v3sigmalapl2[9 * g + 4];
                auto v3sigmalapl2_cbb = v3sigmalapl2[9 * g + 5];
                auto v3sigmalapl2_baa = v3sigmalapl2[9 * g + 6];
                auto v3sigmalapl2_bab = v3sigmalapl2[9 * g + 7];

                auto v3sigmalapltau_aaa = v3sigmalapltau[12 * g + 0];
                auto v3sigmalapltau_aab = v3sigmalapltau[12 * g + 1];
                auto v3sigmalapltau_aba = v3sigmalapltau[12 * g + 2];
                auto v3sigmalapltau_abb = v3sigmalapltau[12 * g + 3];
                auto v3sigmalapltau_caa = v3sigmalapltau[12 * g + 4];
                auto v3sigmalapltau_cab = v3sigmalapltau[12 * g + 5];
                auto v3sigmalapltau_cba = v3sigmalapltau[12 * g + 6];
                auto v3sigmalapltau_cbb = v3sigmalapltau[12 * g + 7];
                auto v3sigmalapltau_baa = v3sigmalapltau[12 * g + 8];
                auto v3sigmalapltau_bab = v3sigmalapltau[12 * g + 9];
                auto v3sigmalapltau_bba = v3sigmalapltau[12 * g + 10];

                auto v3sigmatau2_aaa = v3sigmatau2[9 * g + 0];
                auto v3sigmatau2_aab = v3sigmatau2[9 * g + 1];
                auto v3sigmatau2_abb = v3sigmatau2[9 * g + 2];
                auto v3sigmatau2_caa = v3sigmatau2[9 * g + 3];
                auto v3sigmatau2_cab = v3sigmatau2[9 * g + 4];
                auto v3sigmatau2_cbb = v3sigmatau2[9 * g + 5];
                auto v3sigmatau2_baa = v3sigmatau2[9 * g + 6];
                auto v3sigmatau2_bab = v3sigmatau2[9 * g + 7];

                auto v3lapl3_aaa = v3lapl3[4 * g + 0];
                auto v3lapl3_aab = v3lapl3[4 * g + 1];
                auto v3lapl3_abb = v3lapl3[4 * g + 2];

                auto v3lapl2tau_aaa = v3lapl2tau[6 * g + 0];
                auto v3lapl2tau_aab = v3lapl2tau[6 * g + 1];
                auto v3lapl2tau_aba = v3lapl2tau[6 * g + 2];
                auto v3lapl2tau_abb = v3lapl2tau[6 * g + 3];
                auto v3lapl2tau_bba = v3lapl2tau[6 * g + 4];

                auto v3lapltau2_aaa = v3lapltau2[6 * g + 0];
                auto v3lapltau2_aab = v3lapltau2[6 * g + 1];
                auto v3lapltau2_abb = v3lapltau2[6 * g + 2];
                auto v3lapltau2_baa = v3lapltau2[6 * g + 3];
                auto v3lapltau2_bab = v3lapltau2[6 * g + 4];

                auto v3tau3_aaa = v3tau3[4 * g + 0];
                auto v3tau3_aab = v3tau3[4 * g + 1];
                auto v3tau3_abb = v3tau3[4 * g + 2];
                auto v3tau3_bbb = v3tau3[4 * g + 3];

                auto v3rho2sigma_aaa = v3rho2sigma[9 * g + 0];
                auto v3rho2sigma_aac = v3rho2sigma[9 * g + 1];
                auto v3rho2sigma_aab = v3rho2sigma[9 * g + 2];
                auto v3rho2sigma_aba = v3rho2sigma[9 * g + 3];
                auto v3rho2sigma_abc = v3rho2sigma[9 * g + 4];
                auto v3rho2sigma_abb = v3rho2sigma[9 * g + 5];
                auto v3rho2sigma_bba = v3rho2sigma[9 * g + 6];
                auto v3rho2sigma_bbc = v3rho2sigma[9 * g + 7];

                auto v3rhosigma2_aaa = v3rhosigma2[12 * g + 0];
                auto v3rhosigma2_aac = v3rhosigma2[12 * g + 1];
                auto v3rhosigma2_aab = v3rhosigma2[12 * g + 2];
                auto v3rhosigma2_acc = v3rhosigma2[12 * g + 3];
                auto v3rhosigma2_acb = v3rhosigma2[12 * g + 4];
                auto v3rhosigma2_abb = v3rhosigma2[12 * g + 5];
                auto v3rhosigma2_baa = v3rhosigma2[12 * g + 6];
                auto v3rhosigma2_bac = v3rhosigma2[12 * g + 7];
                auto v3rhosigma2_bab = v3rhosigma2[12 * g + 8];
                auto v3rhosigma2_bcc = v3rhosigma2[12 * g + 9];
                auto v3rhosigma2_bcb = v3rhosigma2[12 * g + 10];

                auto v3sigma3_aaa = v3sigma3[10 * g + 0];
                auto v3sigma3_aac = v3sigma3[10 * g + 1];
                auto v3sigma3_aab = v3sigma3[10 * g + 2];
                auto v3sigma3_acc = v3sigma3[10 * g + 3];
                auto v3sigma3_acb = v3sigma3[10 * g + 4];
                auto v3sigma3_abb = v3sigma3[10 * g + 5];
                auto v3sigma3_ccc = v3sigma3[10 * g + 6];
                auto v3sigma3_ccb = v3sigma3[10 * g + 7];
                auto v3sigma3_cbb = v3sigma3[10 * g + 8];

                // functional derivatives

                // first-order
                double x = vsigma_c + 2.0 * vsigma_a;

                // second-order
                // rho
                double rr = v2rho2_aa + v2rho2_ab;
                double rx = 2.0*v2rhosigma_ac + 2.0*v2rhosigma_ab + 2.0*v2rhosigma_aa;
                double rt = v2rhotau_aa + v2rhotau_ab;
                double rl = v2rholapl_aa + v2rholapl_ab;

                // sigma and gamma
                double xr = v2rhosigma_bc + 2.0*v2rhosigma_ba + v2rhosigma_ac + 2.0 * v2rhosigma_aa;
                double xt = v2sigmatau_cb + 2.0*v2sigmatau_ab + v2sigmatau_ca + 2.0 * v2sigmatau_aa;
                double xl = v2sigmalapl_cb + 2.0*v2sigmalapl_ab + v2sigmalapl_ca + 2.0 * v2sigmalapl_aa;
                double xx = 2.0*v2sigma2_cc + 2.0*v2sigma2_cb + 6.0*v2sigma2_ac + 4.0*v2sigma2_ab + 4.0 * v2sigma2_aa;

                // tau
                double tt = v2tau2_aa + v2tau2_ab;
                double tx = 2.0 * v2sigmatau_ca + 2.0 * v2sigmatau_ba + 2.0 * v2sigmatau_aa;
                double tr = v2rhotau_aa + v2rhotau_ba;
                double tl = v2lapltau_aa + v2lapltau_ba;

                // lapl
                //double ll = v2lapl2_aa + v2lapl2_ab;
                //double lx = 2.0 * v2sigmalapl_ca + 2.0 * v2sigmalapl_ba + 2.0 * v2sigmalapl_aa;
                //double lr = v2rholapl_aa + v2rholapl_ba;
                //double lt = v2lapltau_aa + v2lapltau_ab;

                // Third-oder

                // sigma and gamma
                double xxx  = 4.0 * v3sigma3_ccc + 8.0 * v3sigma3_ccb + 4.0 * v3sigma3_cbb + 16.0 * v3sigma3_acc
                            + 24.0 * v3sigma3_acb + 8.0 * v3sigma3_abb + 20.0 * v3sigma3_aac + 16.0 * v3sigma3_aab + 8.0 * v3sigma3_aaa;

                double xxr  = 2.0 * v3rhosigma2_bcc + 2.0 * v3rhosigma2_bcb + 6.0 * v3rhosigma2_bac
                            + 4.0 * v3rhosigma2_bab + 4.0 * v3rhosigma2_baa + 2.0 * v3rhosigma2_acc
                            + 2.0 * v3rhosigma2_acb + 6.0 * v3rhosigma2_aac + 4.0 * v3rhosigma2_aab + 4.0 * v3rhosigma2_aaa;

                double xrt  = v3rhosigmatau_bcb + v3rhosigmatau_bca
                            + 2.0 * v3rhosigmatau_bab + 2.0 * v3rhosigmatau_baa
                            + 2.0 * v3rhosigmatau_aab + 2.0 * v3rhosigmatau_aaa
                            + v3rhosigmatau_acb + v3rhosigmatau_aca;

                double xxl  = 2.0 * v3sigma2lapl_ccb + 2.0 * v3sigma2lapl_cca + 2.0 * v3sigma2lapl_cbb
                            + 2.0 * v3sigma2lapl_cba + 6.0 * v3sigma2lapl_acb + 6.0 * v3sigma2lapl_aca
                            + 4.0 * v3sigma2lapl_abb + 4.0 * v3sigma2lapl_aba + 4.0 * v3sigma2lapl_aab + 4.0 * v3sigma2lapl_aaa;
                double xrr  = v3rho2sigma_bbc + 2.0 * v3rho2sigma_bba
                            + 2.0 * v3rho2sigma_abc + 4.0 * v3rho2sigma_aba
                            + v3rho2sigma_aac + 2.0 * v3rho2sigma_aaa;

                double xrl  = v3rhosigmalapl_bcb + v3rhosigmalapl_bca + 2.0 * v3rhosigmalapl_bab
                            + 2.0 * v3rhosigmalapl_baa + v3rhosigmalapl_acb + v3rhosigmalapl_aca
                            + 2.0 * v3rhosigmalapl_aab + 2.0 * v3rhosigmalapl_aaa;

                double xtl  = v3sigmalapltau_cbb + v3sigmalapltau_cba + v3sigmalapltau_cab
                            + v3sigmalapltau_caa + 2.0 * v3sigmalapltau_abb + 2.0 * v3sigmalapltau_aba
                            + 2.0 * v3sigmalapltau_aab + 2.0 * v3sigmalapltau_aaa;
                double xll  = v3sigmalapl2_cbb + 2.0 * v3sigmalapl2_cab + v3sigmalapl2_caa
                            + 2.0 * v3sigmalapl2_abb + 4.0 * v3sigmalapl2_aab + 2.0 * v3sigmalapl2_aaa;


                double xxt  = 2.0 * v3sigma2tau_ccb + 2.0 * v3sigma2tau_cca + 2.0 * v3sigma2tau_cbb
                            + 2.0 * v3sigma2tau_cba + 6.0 * v3sigma2tau_acb + 6.0 * v3sigma2tau_aca
                            + 4.0 * v3sigma2tau_abb + 4.0 * v3sigma2tau_aba + 4.0 * v3sigma2tau_aab + 4.0 * v3sigma2tau_aaa;


                double xtt  = v3sigmatau2_cbb + 2.0 * v3sigmatau2_cab + v3sigmatau2_caa
                            + 2.0 * v3sigmatau2_abb + 4.0 * v3sigmatau2_aab + 2.0 * v3sigmatau2_aaa;

                // rho
                double rrr  = v3rho3_abb + 2.0 * v3rho3_aab + v3rho3_aaa;
                double rrt  = v3rho2tau_abb + v3rho2tau_aba + v3rho2tau_aab + v3rho2tau_aaa;
                double rtx  = 2.0 * v3rhosigmatau_acb + 2.0 * v3rhosigmatau_aca
                            + 2.0 * v3rhosigmatau_abb + 2.0 * v3rhosigmatau_aba
                            + 2.0 * v3rhosigmatau_aab + 2.0 * v3rhosigmatau_aaa;
                double rrl  = v3rho2lapl_abb + v3rho2lapl_aba + v3rho2lapl_aab + v3rho2lapl_aaa;
                double rrx  = 2.0 * v3rho2sigma_abc + 2.0 * v3rho2sigma_abb
                            + 2.0 * v3rho2sigma_aba + 2.0 * v3rho2sigma_aac
                            + 2.0 * v3rho2sigma_aab + 2.0 * v3rho2sigma_aaa;
                double rtl  = v3rholapltau_abb + v3rholapltau_aba + v3rholapltau_aab + v3rholapltau_aaa;

                double rtt  = v3rhotau2_abb + 2.0 * v3rhotau2_aab + v3rhotau2_aaa;

                double rll  = v3rholapl2_abb + 2.0 * v3rholapl2_aab + v3rholapl2_aaa;
                double rlx  = 2.0 * v3rhosigmalapl_acb + 2.0 * v3rhosigmalapl_aca + 2.0 * v3rhosigmalapl_abb
                            + 2.0 * v3rhosigmalapl_aba + 2.0 * v3rhosigmalapl_aab + 2.0 * v3rhosigmalapl_aaa;
                double rxx  = 4.0 * v3rhosigma2_acc + 8.0 * v3rhosigma2_acb + 4.0 * v3rhosigma2_abb
                            + 8.0 * v3rhosigma2_aac + 8.0 * v3rhosigma2_aab + 4.0 * v3rhosigma2_aaa;

                // laplacian
                //double lll  = v3lapl3_abb + 2.0 * v3lapl3_aab + v3lapl3_aaa;
                //double llr  = v3rholapl2_bab + v3rholapl2_baa + v3rholapl2_aab + v3rholapl2_aaa;
                //double llt  = v3lapl2tau_abb + v3lapl2tau_aba + v3lapl2tau_aab + v3lapl2tau_aaa;
                //double llx  = 2.0 * v3sigmalapl2_cab + 2.0 * v3sigmalapl2_caa + 2.0 * v3sigmalapl2_bab
                //            + 2.0 * v3sigmalapl2_baa + 2.0 * v3sigmalapl2_aab + 2.0 * v3sigmalapl2_aaa;
                //double lrr  = v3rho2lapl_bba + 2.0 * v3rho2lapl_aba + v3rho2lapl_aaa;
                //double lrt  = v3rholapltau_bab + v3rholapltau_baa + v3rholapltau_aab + v3rholapltau_aaa;
                //double lrx  = 2.0 * v3rhosigmalapl_bca + 2.0 * v3rhosigmalapl_bba + 2.0 * v3rhosigmalapl_baa
                //            + 2.0 * v3rhosigmalapl_aca + 2.0 * v3rhosigmalapl_aba + 2.0 * v3rhosigmalapl_aaa;
                //double ltt  = v3lapltau2_abb + 2.0 * v3lapltau2_aab + v3lapltau2_aaa;
                //double ltx  = 2.0 * v3sigmalapltau_cab + 2.0 * v3sigmalapltau_caa + 2.0 * v3sigmalapltau_bab
                //            + 2.0 * v3sigmalapltau_baa + 2.0 * v3sigmalapltau_aab + 2.0 * v3sigmalapltau_aaa;
                //double lxx  = 4.0 * v3sigma2lapl_cca + 8.0 * v3sigma2lapl_cba + 4.0 * v3sigma2lapl_bba
                //            + 8.0 * v3sigma2lapl_aca + 8.0 * v3sigma2lapl_aba + 4.0 * v3sigma2lapl_aaa;

                // tau
                double trr  = v3rho2tau_bba + 2.0 * v3rho2tau_aba + v3rho2tau_aaa;
                double ttl  = v3lapltau2_bab + v3lapltau2_baa + v3lapltau2_aab + v3lapltau2_aaa;

                double trl  = v3rholapltau_bba + v3rholapltau_baa + v3rholapltau_aba + v3rholapltau_aaa;

                double tll  = v3lapl2tau_bba + 2.0 * v3lapl2tau_aba + v3lapl2tau_aaa;
                double tlx  = 2.0 * v3sigmalapltau_cba + 2.0 * v3sigmalapltau_caa + 2.0 * v3sigmalapltau_bba
                            + 2.0 * v3sigmalapltau_baa + 2.0 * v3sigmalapltau_aba + 2.0 * v3sigmalapltau_aaa;

                double ttt  = v3tau3_abb + 2.0 * v3tau3_aab + v3tau3_aaa;

                double ttr  = v3rhotau2_bab + v3rhotau2_baa + v3rhotau2_aab + v3rhotau2_aaa;

                double trx  = 2.0 * v3rhosigmatau_bca + 2.0 * v3rhosigmatau_bba
                            + 2.0 * v3rhosigmatau_baa + 2.0 * v3rhosigmatau_aca
                            + 2.0 * v3rhosigmatau_aba + 2.0 * v3rhosigmatau_aaa;

                double txx  = 4.0 * v3sigma2tau_cca + 8.0 * v3sigma2tau_cba + 4.0 * v3sigma2tau_bba
                            + 8.0 * v3sigma2tau_aca + 8.0 * v3sigma2tau_aba + 4.0 * v3sigma2tau_aaa;

                double ttx  = 2.0 * v3sigmatau2_cab + 2.0 * v3sigmatau2_caa + 2.0 * v3sigmatau2_bab
                            + 2.0 * v3sigmatau2_baa + 2.0 * v3sigmatau2_aab + 2.0 * v3sigmatau2_aaa;


                double w = weights[g];

                double grada_x_g = rhograd[6 * g + 0];
                double grada_y_g = rhograd[6 * g + 1];
                double grada_z_g = rhograd[6 * g + 2];

                double l2contract =   grada_x_g * gradw_x[g] + grada_y_g * gradw_y[g] + grada_z_g * gradw_z[g];

                double l5contract_x = grada_x_g * l2contract;
                double l5contract_y = grada_y_g * l2contract;
                double l5contract_z = grada_z_g * l2contract;

                double q2contract = grada_x_g * gamx[g] + grada_y_g * gamy[g] + grada_z_g * gamz[g];
                double sl_q2contract = grada_x_g * sl_gamx[g] + grada_y_g * sl_gamy[g] + grada_z_g * sl_gamz[g];
                double st_q2contract = grada_x_g * st_gamx[g] + grada_y_g * st_gamy[g] + grada_z_g * st_gamz[g];


                double q3contract =   grada_x_g * grada_x_g * gamxx[g]
                                    + grada_x_g * grada_y_g * gamxy[g]
                                    + grada_x_g * grada_z_g * gamxz[g]
                                    + grada_y_g * grada_x_g * gamyx[g]
                                    + grada_y_g * grada_y_g * gamyy[g]
                                    + grada_y_g * grada_z_g * gamyz[g]
                                    + grada_z_g * grada_x_g * gamzx[g]
                                    + grada_z_g * grada_y_g * gamzy[g]
                                    + grada_z_g * grada_z_g * gamzz[g];

                double q4contract = gamxx[g] + gamyy[g] + gamzz[g];

                double q7contract_x =  grada_x_g * grada_x_g *  gamx[g] + grada_x_g * grada_y_g *  gamy[g] + grada_x_g * grada_z_g *  gamz[g];
                double q7contract_y =  grada_y_g * grada_x_g *  gamx[g] + grada_y_g * grada_y_g *  gamy[g] + grada_y_g * grada_z_g *  gamz[g];
                double q7contract_z =  grada_z_g * grada_x_g *  gamx[g] + grada_z_g * grada_y_g *  gamy[g] + grada_z_g * grada_z_g *  gamz[g];

                double sl_q7contract_x =  grada_x_g * grada_x_g * sl_gamx[g] + grada_x_g * grada_y_g * sl_gamy[g] + grada_x_g * grada_z_g * sl_gamz[g];
                double sl_q7contract_y =  grada_y_g * grada_x_g * sl_gamx[g] + grada_y_g * grada_y_g * sl_gamy[g] + grada_y_g * grada_z_g * sl_gamz[g];
                double sl_q7contract_z =  grada_z_g * grada_x_g * sl_gamx[g] + grada_z_g * grada_y_g * sl_gamy[g] + grada_z_g * grada_z_g * sl_gamz[g];

                double st_q7contract_x =  grada_x_g * grada_x_g * st_gamx[g] + grada_x_g * grada_y_g * st_gamy[g] + grada_x_g * grada_z_g * st_gamz[g];
                double st_q7contract_y =  grada_y_g * grada_x_g * st_gamx[g] + grada_y_g * grada_y_g * st_gamy[g] + grada_y_g * grada_z_g * st_gamz[g];
                double st_q7contract_z =  grada_z_g * grada_x_g * st_gamx[g] + grada_z_g * grada_y_g * st_gamy[g] + grada_z_g * grada_z_g * st_gamz[g];

                double q8contract_x =  grada_x_g *  gamxx[g] + grada_y_g *  gamxy[g] + grada_z_g *  gamxz[g];
                double q8contract_y =  grada_x_g *  gamyx[g] + grada_y_g *  gamyy[g] + grada_z_g *  gamyz[g];
                double q8contract_z =  grada_x_g *  gamzx[g] + grada_y_g *  gamzy[g] + grada_z_g *  gamzz[g];

                double q9contract_x =  grada_x_g *  q3contract;
                double q9contract_y =  grada_y_g *  q3contract;
                double q9contract_z =  grada_z_g *  q3contract;

                double q10contract_x =  grada_x_g *  gamxx[g] + grada_y_g *  gamyx[g] + grada_z_g *  gamzx[g];
                double q10contract_y =  grada_x_g *  gamxy[g] + grada_y_g *  gamyy[g] + grada_z_g *  gamzy[g];
                double q10contract_z =  grada_x_g *  gamxz[g] + grada_y_g *  gamyz[g] + grada_z_g *  gamzz[g];

                double q11contract_x =  grada_x_g *  gamxx[g] + grada_x_g *  gamyy[g] + grada_x_g *  gamzz[g];
                double q11contract_y =  grada_y_g *  gamxx[g] + grada_y_g *  gamyy[g] + grada_y_g *  gamzz[g];
                double q11contract_z =  grada_z_g *  gamxx[g] + grada_z_g *  gamyy[g] + grada_z_g *  gamzz[g];

                // Rho operator contributions

                // vxc 1 contributions

                double rho_0 =  rr * rhow[g]
                              + rx * l2contract
                              + rt * tauw[g]
                              + rl * laplw[g];

                //double lap_0 =    lr * rhow[g]
                //                + lx * l2contract
                //                + lt * tauw[g]
                //                + ll * laplw[g];

                double tau_0 =     tr * rhow[g]
                                 + tx * l2contract
                                 + tt * tauw[g]
                                 + tl * laplw[g];

                // vxc 2 contributions

                rho_0 += rrr * gam[g]
                       + rrt * rt_gam[g]
                       + rrl * rl_gam[g]
                       + rll * ll_gam[g]
                       + rtt * tt_gam[g]
                       + rtl * tl_gam[g]
                       + rrx * q2contract
                       + rlx * sl_q2contract
                       + rtx * st_q2contract
                       + rxx * q3contract
                       + rx  * q4contract;

                //lap_0 += lrr * gam[g]
                //       + lrt * rt_gam[g]
                //       + llr * rl_gam[g]
                //       + lll * ll_gam[g]
                //       + ltt * tt_gam[g]
                //       + llt * tl_gam[g]
                //       + lrx * q2contract
                //       + llx * sl_q2contract
                //       + ltx * st_q2contract
                //       + lxx * q3contract
                //       + lx  * q4contract;

                tau_0 += trr * gam[g]
                       + ttr * rt_gam[g]
                       + trl * rl_gam[g]
                       + tll * ll_gam[g]
                       + ttt * tt_gam[g]
                       + ttl * tl_gam[g]
                       + trx * q2contract
                       + tlx * sl_q2contract
                       + ttx * st_q2contract
                       + txx * q3contract
                       + tx  * q4contract;


                // Grad operator contributions

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                // vxc 1 contributions

                xcomp +=  grada_x_g * ( xr * rhow[g] + xt * tauw[g] + xl * laplw[g])
                        + x * gradw_x[g]
                        + xx * l5contract_x;

                ycomp += grada_y_g * ( xr * rhow[g] + xt * tauw[g] + xl * laplw[g])
                        + x * gradw_y[g]
                        + xx * l5contract_y;

                zcomp += grada_z_g * ( xr * rhow[g] + xt * tauw[g] + xl * laplw[g])
                        + x * gradw_z[g]
                        + xx * l5contract_z;

                // vxc 2 contributions

                xcomp +=  xrr * grada_x_g * gam[g]
                        + xrt * grada_x_g * rt_gam[g]
                        + xrl * grada_x_g * rl_gam[g]
                        + xll * grada_x_g * ll_gam[g]
                        + xtt * grada_x_g * tt_gam[g]
                        + xtl * grada_x_g * tl_gam[g]
                        + xr * gamx[g] // q6
                        + xl * sl_gamx[g]
                        + xt * st_gamx[g]
                        + xxr * q7contract_x
                        + xxl * sl_q7contract_x
                        + xxt * st_q7contract_x
                        + xx * (q8contract_x + q10contract_x + q11contract_x)
                        + xxx * q9contract_x;

                ycomp +=  xrr * grada_y_g * gam[g] // q5
                        + xrt * grada_y_g * rt_gam[g]
                        + xrl * grada_y_g * rl_gam[g]
                        + xll * grada_y_g * ll_gam[g]
                        + xtt * grada_y_g * tt_gam[g]
                        + xtl * grada_y_g * tl_gam[g]
                        + xr * gamy[g] // q6
                        + xl * sl_gamy[g]
                        + xt * st_gamy[g]
                        + xxr * q7contract_y
                        + xxl * sl_q7contract_y
                        + xxt * st_q7contract_y
                        + xx * (q8contract_y + q10contract_y + q11contract_y)
                        + xxx * q9contract_y;

                zcomp +=  xrr * grada_z_g * gam[g] // q5
                        + xrt * grada_z_g * rt_gam[g]
                        + xrl * grada_z_g * rl_gam[g]
                        + xll * grada_z_g * ll_gam[g]
                        + xtt * grada_z_g * tt_gam[g]
                        + xtl * grada_z_g * tl_gam[g]
                        + xr * gamz[g] // q6
                        + xl * sl_gamz[g]
                        + xt * st_gamz[g]
                        + xxr * q7contract_z
                        + xxl * sl_q7contract_z
                        + xxt * st_q7contract_z
                        + xx * (q8contract_z + q10contract_z + q11contract_z)
                        + xxx * q9contract_z;

                G_val[nu_offset + g] = w * rho_0 * chi_val[nu_offset + g];

                G_gga_val[nu_offset + g] = w * (xcomp * chi_x_val[nu_offset + g] +
                                                ycomp * chi_y_val[nu_offset + g] +
                                                zcomp * chi_z_val[nu_offset + g]);

                // TODO implement Laplacian dependence

                //G_gga_val[nu_offset + g] += w * lap_0 * (chi_xx_val[nu_offset + g] +
                //                                         chi_yy_val[nu_offset + g] +
                //                                         chi_zz_val[nu_offset + g]);

                G_gga_x_val[nu_offset + g] = w * tau_0 * chi_x_val[nu_offset + g];
                G_gga_y_val[nu_offset + g] = w * tau_0 * chi_y_val[nu_offset + g];
                G_gga_z_val[nu_offset + g] = w * tau_0 * chi_z_val[nu_offset + g];
            }
        }
    }

    timer.stop("Kxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Kxc matrix matmul");

    // LDA and GGA contribution
    auto mat_Kxc = denblas::multABt(gtoValues, denblas::addAB(mat_G, mat_G_gga, 2.0));

    // tau contribution
    auto mat_Kxc_x = denblas::multABt(gtoValuesX, mat_G_gga_x);
    auto mat_Kxc_y = denblas::multABt(gtoValuesY, mat_G_gga_y);
    auto mat_Kxc_z = denblas::multABt(gtoValuesZ, mat_G_gga_z);

    mat_Kxc = denblas::addAB(mat_Kxc, mat_Kxc_x, 0.5);
    mat_Kxc = denblas::addAB(mat_Kxc, mat_Kxc_y, 0.5);
    mat_Kxc = denblas::addAB(mat_Kxc, mat_Kxc_z, 0.5);

    mat_Kxc.symmetrizeAndScale(0.5);

    timer.stop("Kxc matrix matmul");

    return mat_Kxc;
}

CDenseMatrix
CXCNewIntegrator::_integratePartialLxcFockForLDA(const int32_t              npoints,
                                                 const double*              weights,
                                                 const CDenseMatrix&        gtoValues,
                                                 const double*              v2rho2,
                                                 const double*              v3rho3,
                                                 const double*              v4rho4,
                                                 const CDensityGridCubic&   rwDensityGridCubic,
                                                 const CDensityGrid&        rw3DensityGrid,
                                                 const int32_t              iFock,
                                                 CMultiTimer&               timer) const
{
    timer.start("Lxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    // pointers to perturbed density

    auto rx_ry_rz = rwDensityGridCubic.pi(iFock);
    auto rxy_rz = rwDensityGridCubic.gam(iFock);
    auto r_xyz = rw3DensityGrid.alphaDensity(iFock);

    timer.stop("Lxc matrix prep.");

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Lxc matrix G");

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

            #pragma omp simd aligned(weights, v2rho2, v3rho3, v4rho4, rx_ry_rz, rxy_rz, r_xyz, G_val, chi_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                auto v2rho2_aa = v2rho2[3 * g + 0];
                auto v2rho2_ab = v2rho2[3 * g + 1];
                auto v3rho3_aaa = v3rho3[4 * g + 0];
                auto v3rho3_aab = v3rho3[4 * g + 1];
                auto v3rho3_abb = v3rho3[4 * g + 2];
                auto v4rho4_aaaa = v4rho4[5 * g + 0];
                auto v4rho4_aaab = v4rho4[5 * g + 1];
                auto v4rho4_aabb = v4rho4[5 * g + 2];
                auto v4rho4_abbb = v4rho4[5 * g + 3];

                double rr = (v2rho2_aa + v2rho2_ab);
                double rrr = (v3rho3_aaa + 2.0 * v3rho3_aab + v3rho3_abb);
                double rrrr = (v4rho4_aaaa +  3.0 * v4rho4_aaab + 3.0 * v4rho4_aabb + v4rho4_abbb);

                G_val[nu_offset + g] = weights[g] *
                            (    rr * r_xyz[g]
                            +   rrr * rxy_rz[g]
                            +   rrrr * rx_ry_rz[g]) * chi_val[nu_offset + g];
            }
        }
    }

    timer.stop("Lxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Lxc matrix matmul");

    auto mat_Kxc = denblas::multABt(gtoValues, mat_G);

    timer.stop("Lxc matrix matmul");

    return mat_Kxc;
}

CDenseMatrix
CXCNewIntegrator::_integratePartialLxcFockForGGA(const int32_t              npoints,
                                                 const double*              weights,
                                                 const CDenseMatrix&        gtoValues,
                                                 const CDenseMatrix&        gtoValuesX,
                                                 const CDenseMatrix&        gtoValuesY,
                                                 const CDenseMatrix&        gtoValuesZ,
                                                 const double*              rhograd,
                                                 const double*              vsigma,
                                                 const double*              v2rho2,
                                                 const double*              v2rhosigma,
                                                 const double*              v2sigma2,
                                                 const double*              v3rho3,
                                                 const double*              v3rho2sigma,
                                                 const double*              v3rhosigma2,
                                                 const double*              v3sigma3,
                                                 const double*              v4rho4,
                                                 const double*              v4rho3sigma,
                                                 const double*              v4rho2sigma2,
                                                 const double*              v4rhosigma3,
                                                 const double*              v4sigma4,                   
                                                 const CDensityGridCubic&   rwDensityGridCubic,
                                                 const CDensityGrid&        rw3DensityGrid,
                                                 const int32_t              iFock,
                                                 CMultiTimer&               timer) const
{
    timer.start("Lxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    // pointers to perturbed density gradient norms

    auto gam = rwDensityGridCubic.gam(iFock);
    
    auto gamx = rwDensityGridCubic.gamX(iFock);
    
    auto gamy = rwDensityGridCubic.gamY(iFock);
    
    auto gamz = rwDensityGridCubic.gamZ(iFock);

    auto gamxx = rwDensityGridCubic.gamXX(iFock);

    auto gamxy = rwDensityGridCubic.gamXY(iFock); 

    auto gamxz = rwDensityGridCubic.gamXZ(iFock);
    
    auto gamyx = rwDensityGridCubic.gamYX(iFock);
    
    auto gamyy = rwDensityGridCubic.gamYY(iFock);
    
    auto gamyz = rwDensityGridCubic.gamYZ(iFock);
    
    auto gamzx = rwDensityGridCubic.gamZX(iFock);
    
    auto gamzy = rwDensityGridCubic.gamZY(iFock);

    auto gamzz = rwDensityGridCubic.gamZZ(iFock);

    auto pi = rwDensityGridCubic.pi(iFock);

    auto pix = rwDensityGridCubic.piX(iFock);

    auto piy = rwDensityGridCubic.piY(iFock);

    auto piz = rwDensityGridCubic.piZ(iFock);

    auto pixx = rwDensityGridCubic.piXX(iFock);

    auto pixy = rwDensityGridCubic.piXY(iFock);

    auto pixz = rwDensityGridCubic.piXZ(iFock);

    auto piyx = rwDensityGridCubic.piYX(iFock);

    auto piyy = rwDensityGridCubic.piYY(iFock);

    auto piyz = rwDensityGridCubic.piYZ(iFock);

    auto pizx = rwDensityGridCubic.piZX(iFock);

    auto pizy = rwDensityGridCubic.piZY(iFock);

    auto pizz = rwDensityGridCubic.piZZ(iFock);

    auto pixxx = rwDensityGridCubic.piXXX(iFock);

    auto pixxy = rwDensityGridCubic.piXXY(iFock);

    auto pixxz = rwDensityGridCubic.piXXZ(iFock);

    auto pixyx = rwDensityGridCubic.piXYX(iFock);

    auto pixyy = rwDensityGridCubic.piXYY(iFock);

    auto pixyz = rwDensityGridCubic.piXYZ(iFock);

    auto pixzx = rwDensityGridCubic.piXZX(iFock);

    auto pixzy = rwDensityGridCubic.piXZY(iFock);

    auto pixzz = rwDensityGridCubic.piXZZ(iFock);

    auto piyxx = rwDensityGridCubic.piYXX(iFock);

    auto piyxy = rwDensityGridCubic.piYXY(iFock);

    auto piyxz = rwDensityGridCubic.piYXZ(iFock);

    auto piyyx = rwDensityGridCubic.piYYX(iFock);

    auto piyyy = rwDensityGridCubic.piYYY(iFock);

    auto piyyz = rwDensityGridCubic.piYYZ(iFock);

    auto piyzx = rwDensityGridCubic.piYZX(iFock);

    auto piyzy = rwDensityGridCubic.piYZY(iFock);

    auto piyzz = rwDensityGridCubic.piYZZ(iFock);

    auto pizxx = rwDensityGridCubic.piZXX(iFock);

    auto pizxy = rwDensityGridCubic.piZXY(iFock);

    auto pizxz = rwDensityGridCubic.piZXZ(iFock);

    auto pizyx = rwDensityGridCubic.piZYX(iFock);

    auto pizyy = rwDensityGridCubic.piZYY(iFock);

    auto pizyz = rwDensityGridCubic.piZYZ(iFock);

    auto pizzx = rwDensityGridCubic.piZZX(iFock);

    auto pizzy = rwDensityGridCubic.piZZY(iFock);

    auto pizzz = rwDensityGridCubic.piZZZ(iFock);

    auto rho3 = rw3DensityGrid.alphaDensity(iFock);
    
    auto grad3_x = rw3DensityGrid.alphaDensityGradientX(iFock);
    
    auto grad3_y = rw3DensityGrid.alphaDensityGradientY(iFock);
    
    auto grad3_z = rw3DensityGrid.alphaDensityGradientZ(iFock);


    timer.stop("Lxc matrix prep.");

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Lxc matrix G");

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
                    gam, gamx, gamy, gamz  ,gamxx ,gamxy ,gamxz ,gamyx , \
                    gamyy ,gamyz ,gamzx ,gamzy ,gamzz ,pi  ,pix  ,piy  ,piz  ,pixx  , pixy  ,\
                    pixz  ,piyx  ,piyy  ,piyz  ,pizx  ,pizy  ,pizz  ,pixxx ,pixxy ,pixxz ,pixyx , \
                    pixyy ,pixyz ,pixzx ,pixzy ,pixzz ,piyxx ,piyxy ,piyxz ,piyyx ,piyyy ,piyyz ,piyzx , \
                    piyzy ,piyzz ,pizxx ,pizxy ,pizxz ,pizyx ,pizyy ,pizyz ,pizzx ,pizzy ,pizzz ,\
                    rho3  ,grad3_x,grad3_y,grad3_z, vsigma, v2rho2, v2rhosigma, v2sigma2, rhograd, \
                    v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3, \
                    v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4,\
                    G_val, G_gga_val, chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                double w = weights[g];

                double rxw123a = grad3_x[g];
                double ryw123a = grad3_y[g];
                double rzw123a = grad3_z[g];

                double grada_x_g = rhograd[6 * g + 0];
                double grada_y_g = rhograd[6 * g + 1];
                double grada_z_g = rhograd[6 * g + 2];

                double l2contract = grada_x_g * rxw123a + grada_y_g * ryw123a + grada_z_g * rzw123a;
                double l5contract_x = grada_x_g * l2contract;
                double l5contract_y = grada_y_g * l2contract;
                double l5contract_z = grada_z_g * l2contract; 

                // vx2 terms
                double q2contract = grada_x_g * gamx[g] + grada_y_g * gamy[g] + grada_z_g * gamz[g];

                double q3contract =   grada_x_g * grada_x_g * gamxx[g]
                                    + grada_x_g * grada_y_g * gamxy[g]
                                    + grada_x_g * grada_z_g * gamxz[g]
                                    + grada_y_g * grada_x_g * gamyx[g]
                                    + grada_y_g * grada_y_g * gamyy[g]
                                    + grada_y_g * grada_z_g * gamyz[g]
                                    + grada_z_g * grada_x_g * gamzx[g]
                                    + grada_z_g * grada_y_g * gamzy[g]
                                    + grada_z_g * grada_z_g * gamzz[g];

                double q4contract = gamxx[g] + gamyy[g] + gamzz[g];

                double q7contract_x =  grada_x_g * q2contract;
                double q7contract_y =  grada_y_g * q2contract;
                double q7contract_z =  grada_z_g * q2contract;

                double q8contract_x =  grada_x_g *  gamxx[g] + grada_y_g *  gamxy[g] + grada_z_g *  gamxz[g];
                double q8contract_y =  grada_x_g *  gamyx[g] + grada_y_g *  gamyy[g] + grada_z_g *  gamyz[g];
                double q8contract_z =  grada_x_g *  gamzx[g] + grada_y_g *  gamzy[g] + grada_z_g *  gamzz[g];

                double q10contract_x =  grada_x_g *  gamxx[g] + grada_y_g *  gamyx[g] + grada_z_g *  gamzx[g];
                double q10contract_y =  grada_x_g *  gamxy[g] + grada_y_g *  gamyy[g] + grada_z_g *  gamzy[g];
                double q10contract_z =  grada_x_g *  gamxz[g] + grada_y_g *  gamyz[g] + grada_z_g *  gamzz[g];

                double q11contract_x =  grada_x_g *  gamxx[g] + grada_x_g *  gamyy[g] + grada_x_g *  gamzz[g];
                double q11contract_y =  grada_y_g *  gamxx[g] + grada_y_g *  gamyy[g] + grada_y_g *  gamzz[g];
                double q11contract_z =  grada_z_g *  gamxx[g] + grada_z_g *  gamyy[g] + grada_z_g *  gamzz[g];

                double q9contract_x =  grada_x_g *  q3contract;
                double q9contract_y =  grada_y_g *  q3contract;
                double q9contract_z =  grada_z_g *  q3contract;

                // vx3 terms
                double c1 = pi[g];

                double c2 = grada_x_g * pix[g] + grada_y_g * piy[g] + grada_z_g * piz[g];

                double c3 =  grada_x_g * grada_x_g * pixx[g]
                            +grada_x_g * grada_y_g * pixy[g]
                            +grada_x_g * grada_z_g * pixz[g]
                            +grada_y_g * grada_x_g * piyx[g]
                            +grada_y_g * grada_y_g * piyy[g]
                            +grada_y_g * grada_z_g * piyz[g]
                            +grada_z_g * grada_x_g * pizx[g]
                            +grada_z_g * grada_y_g * pizy[g]
                            +grada_z_g * grada_z_g * pizz[g];

                double c4 = pixx[g] + piyy[g] + pizz[g];

                double c5_6 = grada_x_g *(pixxx[g] + pixxx[g])
                            +grada_x_g * (piyxy[g] + pixyy[g])
                            +grada_x_g * (pizxz[g] + pixzz[g])
                            +grada_y_g * (pixyx[g] + piyxx[g])
                            +grada_y_g * (piyyy[g] + piyyy[g])
                            +grada_y_g * (pizyz[g] + piyzz[g])
                            +grada_z_g * (pixzx[g] + pizxx[g])
                            +grada_z_g * (piyzy[g] + pizyy[g])
                            +grada_z_g * (pizzz[g] + pizzz[g]);
                
                double c7 =   grada_x_g * grada_x_g * grada_x_g * pixxx[g]
                            + grada_x_g * grada_x_g * grada_y_g * pixxy[g]
                            + grada_x_g * grada_x_g * grada_z_g * pixxz[g]
                            + grada_x_g * grada_y_g * grada_x_g * pixyx[g]
                            + grada_x_g * grada_y_g * grada_y_g * pixyy[g]
                            + grada_x_g * grada_y_g * grada_z_g * pixyz[g]
                            + grada_x_g * grada_z_g * grada_x_g * pixzx[g]
                            + grada_x_g * grada_z_g * grada_y_g * pixzy[g]
                            + grada_x_g * grada_z_g * grada_z_g * pixzz[g]
                            + grada_y_g * grada_x_g * grada_x_g * piyxx[g]
                            + grada_y_g * grada_x_g * grada_y_g * piyxy[g]
                            + grada_y_g * grada_x_g * grada_z_g * piyxz[g]
                            + grada_y_g * grada_y_g * grada_x_g * piyyx[g]
                            + grada_y_g * grada_y_g * grada_y_g * piyyy[g]
                            + grada_y_g * grada_y_g * grada_z_g * piyyz[g]
                            + grada_y_g * grada_z_g * grada_x_g * piyzx[g]
                            + grada_y_g * grada_z_g * grada_y_g * piyzy[g]
                            + grada_y_g * grada_z_g * grada_z_g * piyzz[g]
                            + grada_z_g * grada_x_g * grada_x_g * pizxx[g]
                            + grada_z_g * grada_x_g * grada_y_g * pizxy[g]
                            + grada_z_g * grada_x_g * grada_z_g * pizxz[g]
                            + grada_z_g * grada_y_g * grada_x_g * pizyx[g]
                            + grada_z_g * grada_y_g * grada_y_g * pizyy[g]
                            + grada_z_g * grada_y_g * grada_z_g * pizyz[g]
                            + grada_z_g * grada_z_g * grada_x_g * pizzx[g]
                            + grada_z_g * grada_z_g * grada_y_g * pizzy[g]
                            + grada_z_g * grada_z_g * grada_z_g * pizzz[g];

                double c8 =  grada_x_g * pixxx[g]
                            +grada_y_g * pixxy[g]
                            +grada_z_g * pixxz[g]
                            +grada_x_g * piyyx[g]
                            +grada_y_g * piyyy[g]
                            +grada_z_g * piyyz[g]
                            +grada_x_g * pizzx[g]
                            +grada_y_g * pizzy[g]
                            +grada_z_g * pizzz[g];
                
                double c9_x = grada_x_g * pi[g];
                double c9_y = grada_y_g * pi[g];
                double c9_z = grada_z_g * pi[g];

                double c10_x = pix[g];
                double c10_y = piy[g];
                double c10_z = piz[g];

                double c11_x = c2 * grada_x_g;
                double c11_y = c2 * grada_y_g;
                double c11_z = c2 * grada_z_g;

                double c12_c14_x =   grada_x_g * (pixx[g] + pixx[g])
                                    +grada_y_g * (pixy[g] + piyx[g])
                                    +grada_z_g * (pixz[g] + pizx[g]);

                double c12_c14_y= grada_x_g *(piyx[g] + pixy[g])
                                 + grada_y_g * (piyy[g] + piyy[g])
                                 + grada_z_g * (piyz[g] + pizy[g]);

                double c12_c14_z= grada_x_g * (pizx[g] + pixz[g])
                                 + grada_y_g * (pizy[g] + piyz[g])
                                 + grada_z_g * (pizz[g] + pizz[g]);

                double c13 = grada_x_g * grada_x_g * pixx[g]
                            +grada_x_g * grada_y_g * pixy[g]
                            +grada_x_g * grada_z_g * pixz[g]
                            +grada_y_g * grada_x_g * piyx[g]
                            +grada_y_g * grada_y_g * piyy[g]
                            +grada_y_g * grada_z_g * piyz[g]
                            +grada_z_g * grada_x_g * pizx[g]
                            +grada_z_g * grada_y_g * pizy[g]
                            +grada_z_g * grada_z_g * pizz[g];
                            
                double c13_x = c13 * grada_x_g;
                double c13_y = c13 * grada_y_g;
                double c13_z = c13 * grada_z_g;

                double c15_x = grada_x_g * c4;
                double c15_y = grada_y_g * c4;
                double c15_z = grada_z_g * c4;

                double c16_19_22_x = grada_x_g * grada_x_g * pixxx[g] + grada_x_g * grada_x_g * pixxx[g] + grada_x_g * grada_x_g * pixxx[g]
                                    +grada_x_g * grada_y_g * pixxy[g] + grada_x_g * grada_y_g * pixyx[g] + grada_x_g * grada_y_g * pixxy[g]
                                    +grada_x_g * grada_z_g * pixxz[g] + grada_x_g * grada_z_g * pixzx[g] + grada_x_g * grada_z_g * pixxz[g]
                                    +grada_y_g * grada_x_g * pixyx[g] + grada_y_g * grada_x_g * piyxx[g] + grada_y_g * grada_x_g * piyxx[g]
                                    +grada_y_g * grada_y_g * pixyy[g] + grada_y_g * grada_y_g * piyyx[g] + grada_y_g * grada_y_g * piyxy[g]
                                    +grada_y_g * grada_z_g * pixyz[g] + grada_y_g * grada_z_g * piyzx[g] + grada_y_g * grada_z_g * piyxz[g]
                                    +grada_z_g * grada_x_g * pixzx[g] + grada_z_g * grada_x_g * pizxx[g] + grada_z_g * grada_x_g * pizxx[g]
                                    +grada_z_g * grada_y_g * pixzy[g] + grada_z_g * grada_y_g * pizyx[g] + grada_z_g * grada_y_g * pizxy[g]
                                    +grada_z_g * grada_z_g * pixzz[g] + grada_z_g * grada_z_g * pizzx[g] + grada_z_g * grada_z_g * pizxz[g];

                double c16_19_22_y = grada_x_g * grada_x_g * piyxx[g] + grada_x_g * grada_x_g * pixxy[g] + grada_x_g * grada_x_g * pixyx[g]
                                    +grada_x_g * grada_y_g * piyxy[g] + grada_x_g * grada_y_g * pixyy[g] + grada_x_g * grada_y_g * pixyy[g]
                                    +grada_x_g * grada_z_g * piyxz[g] + grada_x_g * grada_z_g * pixzy[g] + grada_x_g * grada_z_g * pixyz[g]
                                    +grada_y_g * grada_x_g * piyyx[g] + grada_y_g * grada_x_g * piyxy[g] + grada_y_g * grada_x_g * piyyx[g]
                                    +grada_y_g * grada_y_g * piyyy[g] + grada_y_g * grada_y_g * piyyy[g] + grada_y_g * grada_y_g * piyyy[g]
                                    +grada_y_g * grada_z_g * piyyz[g] + grada_y_g * grada_z_g * piyzy[g] + grada_y_g * grada_z_g * piyyz[g]
                                    +grada_z_g * grada_x_g * piyzx[g] + grada_z_g * grada_x_g * pizxy[g] + grada_z_g * grada_x_g * pizyx[g]
                                    +grada_z_g * grada_y_g * piyzy[g] + grada_z_g * grada_y_g * pizyy[g] + grada_z_g * grada_y_g * pizyy[g]
                                    +grada_z_g * grada_z_g * piyzz[g] + grada_z_g * grada_z_g * pizzy[g] + grada_z_g * grada_z_g * pizyz[g];

                double c16_19_22_z = grada_x_g * grada_x_g * pizxx[g] + grada_x_g * grada_x_g * pixxz[g] + grada_x_g * grada_x_g * pixzx[g]
                                    +grada_x_g * grada_y_g * pizxy[g] + grada_x_g * grada_y_g * pixyz[g] + grada_x_g * grada_y_g * pixzy[g]
                                    +grada_x_g * grada_z_g * pizxz[g] + grada_x_g * grada_z_g * pixzz[g] + grada_x_g * grada_z_g * pixzz[g]
                                    +grada_y_g * grada_x_g * pizyx[g] + grada_y_g * grada_x_g * piyxz[g] + grada_y_g * grada_x_g * piyzx[g]
                                    +grada_y_g * grada_y_g * pizyy[g] + grada_y_g * grada_y_g * piyyz[g] + grada_y_g * grada_y_g * piyzy[g]
                                    +grada_y_g * grada_z_g * pizyz[g] + grada_y_g * grada_z_g * piyzz[g] + grada_y_g * grada_z_g * piyzz[g]
                                    +grada_z_g * grada_x_g * pizzx[g] + grada_z_g * grada_x_g * pizxz[g] + grada_z_g * grada_x_g * pizzx[g]
                                    +grada_z_g * grada_y_g * pizzy[g] + grada_z_g * grada_y_g * pizyz[g] + grada_z_g * grada_y_g * pizzy[g]
                                    +grada_z_g * grada_z_g * pizzz[g] + grada_z_g * grada_z_g * pizzz[g] + grada_z_g * grada_z_g * pizzz[g];

                double c17_24_25_x =  pixxx[g] + pixxx[g] + pixxx[g]
                                    + pixyy[g] + piyxy[g] + piyyx[g]
                                    + pixzz[g] + pizxz[g] + pizzx[g];

                double c17_24_25_y =  piyxx[g] + pixyx[g] + pixxy[g]
                                    + piyyy[g] + piyyy[g] + piyyy[g]
                                    + piyzz[g] + pizyz[g] + pizzy[g];

                double c17_24_25_z =  pizxx[g] + pixzx[g] + pixxz[g]
                                    + pizyy[g] + piyzy[g] + piyyz[g]
                                    + pizzz[g] + pizzz[g] + pizzz[g];

                double c18_x = c7 * grada_x_g;
                double c18_y = c7 * grada_y_g;
                double c18_z = c7 * grada_z_g;


                double c20_21_23 =   grada_x_g * (pixxx[g] + pixxx[g] + pixxx[g]) 
                                    +grada_x_g * (piyxy[g] + pixyy[g] + piyyx[g]) 
                                    +grada_x_g * (pizxz[g] + pixzz[g] + pizzx[g]) 
                                    +grada_y_g * (pixyx[g] + piyxx[g] + pixxy[g]) 
                                    +grada_y_g * (piyyy[g] + piyyy[g] + piyyy[g]) 
                                    +grada_y_g * (pizyz[g] + piyzz[g] + pizzy[g]) 
                                    +grada_z_g * (pixzx[g] + pizxx[g] + pixxz[g]) 
                                    +grada_z_g * (piyzy[g] + pizyy[g] + piyyz[g]) 
                                    +grada_z_g * (pizzz[g] + pizzz[g] + pizzz[g]);

                double c20_21_23_x = grada_x_g * c20_21_23; 
                double c20_21_23_y = grada_y_g * c20_21_23; 
                double c20_21_23_z = grada_z_g * c20_21_23;

                // functional derivatives in libxc form

                auto vsigma_a = vsigma[3 * g + 0];
                auto vsigma_c = vsigma[3 * g + 1];

                auto v2rho2_aa = v2rho2[3 * g + 0];
                auto v2rho2_ab = v2rho2[3 * g + 1];

                auto v2rhosigma_aa = v2rhosigma[6 * g + 0];
                auto v2rhosigma_ac = v2rhosigma[6 * g + 1];
                auto v2rhosigma_ab = v2rhosigma[6 * g + 2];
                auto v2rhosigma_ba = v2rhosigma[6 * g + 3];
                auto v2rhosigma_bc = v2rhosigma[6 * g + 4];

                auto v2sigma2_aa = v2sigma2[6 * g + 0];
                auto v2sigma2_ac = v2sigma2[6 * g + 1];
                auto v2sigma2_ab = v2sigma2[6 * g + 2];
                auto v2sigma2_cc = v2sigma2[6 * g + 3];
                auto v2sigma2_cb = v2sigma2[6 * g + 4];

                auto v3rho3_aaa = v3rho3[4 * g + 0];
                auto v3rho3_aab = v3rho3[4 * g + 1];
                auto v3rho3_abb = v3rho3[4 * g + 2];

                auto v3rho2sigma_aaa = v3rho2sigma[9 * g + 0];
                auto v3rho2sigma_aac = v3rho2sigma[9 * g + 1];
                auto v3rho2sigma_aab = v3rho2sigma[9 * g + 2];
                auto v3rho2sigma_aba = v3rho2sigma[9 * g + 3];
                auto v3rho2sigma_abc = v3rho2sigma[9 * g + 4];
                auto v3rho2sigma_abb = v3rho2sigma[9 * g + 5];
                auto v3rho2sigma_bba = v3rho2sigma[9 * g + 6];
                auto v3rho2sigma_bbc = v3rho2sigma[9 * g + 7];

                auto v3rhosigma2_aaa = v3rhosigma2[12 * g + 0];
                auto v3rhosigma2_aac = v3rhosigma2[12 * g + 1];
                auto v3rhosigma2_aab = v3rhosigma2[12 * g + 2];
                auto v3rhosigma2_acc = v3rhosigma2[12 * g + 3];
                auto v3rhosigma2_acb = v3rhosigma2[12 * g + 4];
                auto v3rhosigma2_abb = v3rhosigma2[12 * g + 5];
                auto v3rhosigma2_baa = v3rhosigma2[12 * g + 6];
                auto v3rhosigma2_bac = v3rhosigma2[12 * g + 7];
                auto v3rhosigma2_bab = v3rhosigma2[12 * g + 8];
                auto v3rhosigma2_bcc = v3rhosigma2[12 * g + 9];
                auto v3rhosigma2_bcb = v3rhosigma2[12 * g + 10];

                auto v3sigma3_aaa = v3sigma3[10 * g + 0];
                auto v3sigma3_aac = v3sigma3[10 * g + 1];
                auto v3sigma3_aab = v3sigma3[10 * g + 2];
                auto v3sigma3_acc = v3sigma3[10 * g + 3];
                auto v3sigma3_acb = v3sigma3[10 * g + 4];
                auto v3sigma3_abb = v3sigma3[10 * g + 5];
                auto v3sigma3_ccc = v3sigma3[10 * g + 6];
                auto v3sigma3_bcc = v3sigma3[10 * g + 7];
                auto v3sigma3_cbb = v3sigma3[10 * g + 8];

                auto v4rho4_aaaa = v4rho4[5 * g + 0];
                auto v4rho4_aaab = v4rho4[5 * g + 1];
                auto v4rho4_aabb = v4rho4[5 * g + 2];
                auto v4rho4_abbb = v4rho4[5 * g + 3];
                auto v4rho4_bbbb = v4rho4[5 * g + 4];

                auto v4rho3sigma_aaaa = v4rho3sigma[12 * g + 0];
                auto v4rho3sigma_aaac = v4rho3sigma[12 * g + 1];
                auto v4rho3sigma_aaab = v4rho3sigma[12 * g + 2];
                auto v4rho3sigma_aaba = v4rho3sigma[12 * g + 3];
                auto v4rho3sigma_aabc = v4rho3sigma[12 * g + 4];
                auto v4rho3sigma_aabb = v4rho3sigma[12 * g + 5];
                auto v4rho3sigma_abba = v4rho3sigma[12 * g + 6];
                auto v4rho3sigma_abbc = v4rho3sigma[12 * g + 7];
                auto v4rho3sigma_abbb = v4rho3sigma[12 * g + 8];
                auto v4rho3sigma_bbba = v4rho3sigma[12 * g + 9];
                auto v4rho3sigma_bbbc = v4rho3sigma[12 * g + 10];

                auto v4rho2sigma2_aaaa = v4rho2sigma2[18 * g + 0];
                auto v4rho2sigma2_aaac = v4rho2sigma2[18 * g + 1];
                auto v4rho2sigma2_aaab = v4rho2sigma2[18 * g + 2];
                auto v4rho2sigma2_aacc = v4rho2sigma2[18 * g + 3];
                auto v4rho2sigma2_aacb = v4rho2sigma2[18 * g + 4];
                auto v4rho2sigma2_aabb = v4rho2sigma2[18 * g + 5];
                auto v4rho2sigma2_abaa = v4rho2sigma2[18 * g + 6];
                auto v4rho2sigma2_abac = v4rho2sigma2[18 * g + 7];
                auto v4rho2sigma2_abab = v4rho2sigma2[18 * g + 8];
                auto v4rho2sigma2_abcc = v4rho2sigma2[18 * g + 9];
                auto v4rho2sigma2_abcb = v4rho2sigma2[18 * g + 10];
                auto v4rho2sigma2_abbb = v4rho2sigma2[18 * g + 11];
                auto v4rho2sigma2_bbaa = v4rho2sigma2[18 * g + 12];
                auto v4rho2sigma2_bbac = v4rho2sigma2[18 * g + 13];
                auto v4rho2sigma2_bbab = v4rho2sigma2[18 * g + 14];
                auto v4rho2sigma2_bbcc = v4rho2sigma2[18 * g + 15];
                auto v4rho2sigma2_bbcb = v4rho2sigma2[18 * g + 16];
                auto v4rho2sigma2_bbbb = v4rho2sigma2[18 * g + 17];

                auto v4rhosigma3_aaaa = v4rhosigma3[20 * g + 0];
                auto v4rhosigma3_aaac = v4rhosigma3[20 * g + 1];
                auto v4rhosigma3_aaab = v4rhosigma3[20 * g + 2];
                auto v4rhosigma3_aacc = v4rhosigma3[20 * g + 3];
                auto v4rhosigma3_aacb = v4rhosigma3[20 * g + 4];
                auto v4rhosigma3_aabb = v4rhosigma3[20 * g + 5];
                auto v4rhosigma3_accc = v4rhosigma3[20 * g + 6];
                auto v4rhosigma3_accb = v4rhosigma3[20 * g + 7];
                auto v4rhosigma3_acbb = v4rhosigma3[20 * g + 8];
                auto v4rhosigma3_abbb = v4rhosigma3[20 * g + 9];
                auto v4rhosigma3_baaa = v4rhosigma3[20 * g + 10];
                auto v4rhosigma3_baac = v4rhosigma3[20 * g + 11];
                auto v4rhosigma3_baab = v4rhosigma3[20 * g + 12];
                auto v4rhosigma3_bacc = v4rhosigma3[20 * g + 13];
                auto v4rhosigma3_bacb = v4rhosigma3[20 * g + 14];
                auto v4rhosigma3_babb = v4rhosigma3[20 * g + 15];
                auto v4rhosigma3_bccc = v4rhosigma3[20 * g + 16];
                auto v4rhosigma3_bccb = v4rhosigma3[20 * g + 17];
                auto v4rhosigma3_bcbb = v4rhosigma3[20 * g + 18];
                auto v4rhosigma3_bbbb = v4rhosigma3[20 * g + 19];

                auto v4sigma4_aaaa = v4sigma4[15 * g + 0];
                auto v4sigma4_aaac = v4sigma4[15 * g + 1];
                auto v4sigma4_aaab = v4sigma4[15 * g + 2];
                auto v4sigma4_aacc = v4sigma4[15 * g + 3];
                auto v4sigma4_aacb = v4sigma4[15 * g + 4];
                auto v4sigma4_aabb = v4sigma4[15 * g + 5];
                auto v4sigma4_accc = v4sigma4[15 * g + 6];
                auto v4sigma4_accb = v4sigma4[15 * g + 7];
                auto v4sigma4_acbb = v4sigma4[15 * g + 8];
                auto v4sigma4_abbb = v4sigma4[15 * g + 9];
                auto v4sigma4_cccc = v4sigma4[15 * g + 10];
                auto v4sigma4_cccb = v4sigma4[15 * g + 11];
                auto v4sigma4_ccbb = v4sigma4[15 * g + 12];
                auto v4sigma4_cbbb = v4sigma4[15 * g + 13];
                auto v4sigma4_bbbb = v4sigma4[15 * g + 14];

                // Transformation of derivatives 

                // first-order 
                double x = vsigma_c + 2.0*vsigma_a;
                
                // second-order
                double rr = (v2rho2_aa + v2rho2_ab);
                double rx = (2.0*v2rhosigma_ac + 2.0*v2rhosigma_ab + 2.0*v2rhosigma_aa);
                double xr = v2rhosigma_bc + 2.0*v2rhosigma_ba + v2rhosigma_ac + 2.0*v2rhosigma_aa;
                double xx = 2.0*v2sigma2_cc + 2.0*v2sigma2_cb + 6.0*v2sigma2_ac + 4.0*v2sigma2_ab + 4.0*v2sigma2_aa;

                // third-order
                double rrr = (v3rho3_aaa + 2.0 * v3rho3_aab + v3rho3_abb);
                double rxr = (2.0*v3rho2sigma_abc + 2.0*v3rho2sigma_abb + 2.0*v3rho2sigma_aba 
                             + 2.0*v3rho2sigma_aac + 2.0*v3rho2sigma_aab + 2.0*v3rho2sigma_aaa);
                double rxx = (4.0*v3rhosigma2_acc + 8.0*v3rhosigma2_acb + 4.0*v3rhosigma2_abb 
                            + 8.0*v3rhosigma2_aac + 8.0*v3rhosigma2_aab + 4.0*v3rhosigma2_aaa);
                double xrr = v3rho2sigma_bbc + 2.0*v3rho2sigma_bba + 2.0*v3rho2sigma_abc + 4.0*v3rho2sigma_aba 
                            + v3rho2sigma_aac + 2.0*v3rho2sigma_aaa;
                double xxr = 2.0*v3rhosigma2_bcc + 2.0*v3rhosigma2_bcb + 6.0*v3rhosigma2_bac 
                            + 4.0*v3rhosigma2_bab + 4.0*v3rhosigma2_baa + 2.0*v3rhosigma2_acc 
                            + 2.0*v3rhosigma2_acb + 6.0*v3rhosigma2_aac + 4.0*v3rhosigma2_aab + 4.0*v3rhosigma2_aaa;
                double xxx = 4.0*v3sigma3_ccc + 8.0*v3sigma3_bcc + 4.0*v3sigma3_cbb + 16.0*v3sigma3_acc 
                            + 24.0*v3sigma3_acb + 8.0*v3sigma3_abb + 20.0*v3sigma3_aac 
                            + 16.0*v3sigma3_aab + 8.0*v3sigma3_aaa;

                // fourth-order 

                double xrrr =   v4rho3sigma_bbbc + 2.0 * v4rho3sigma_bbba + 3.0 * v4rho3sigma_abbc 
                                + 6.0 * v4rho3sigma_abba + 3.0 * v4rho3sigma_aabc + 6.0 * v4rho3sigma_aaba 
                                + v4rho3sigma_aaac + 2.0 * v4rho3sigma_aaaa;

                double xxrr = 2.0 * v4rho2sigma2_bbcc + 2.0 * v4rho2sigma2_bbcb + 6.0 * v4rho2sigma2_bbac 
                            + 4.0 * v4rho2sigma2_bbab + 4.0 * v4rho2sigma2_bbaa + 4.0 * v4rho2sigma2_abcc 
                            + 4.0 * v4rho2sigma2_abcb + 12.0 * v4rho2sigma2_abac + 8.0 * v4rho2sigma2_abab 
                            + 8.0 * v4rho2sigma2_abaa + 2.0 * v4rho2sigma2_aacc + 2.0 * v4rho2sigma2_aacb 
                            + 6.0 * v4rho2sigma2_aaac + 4.0 * v4rho2sigma2_aaab + 4.0 * v4rho2sigma2_aaaa; 

                double xxxr = 4.0 * v4rhosigma3_bccc + 8.0 * v4rhosigma3_bccb + 4.0 * v4rhosigma3_bcbb 
                            + 16.0 * v4rhosigma3_bacc + 24.0 * v4rhosigma3_bacb + 8.0 * v4rhosigma3_babb 
                            + 20.0 * v4rhosigma3_baac + 16.0 * v4rhosigma3_baab + 8.0 * v4rhosigma3_baaa 
                            + 4.0 * v4rhosigma3_accc + 8.0 * v4rhosigma3_accb + 4.0 * v4rhosigma3_acbb 
                            + 16.0 * v4rhosigma3_aacc + 24.0 * v4rhosigma3_aacb + 8.0 * v4rhosigma3_aabb 
                            + 20.0 * v4rhosigma3_aaac + 16.0 * v4rhosigma3_aaab + 8.0 * v4rhosigma3_aaaa;

                double xxxx = 8.0 * v4sigma4_cccc + 24.0 * v4sigma4_cccb + 24.0 * v4sigma4_ccbb + 8.0 * v4sigma4_cbbb 
                            + 40.0 * v4sigma4_accc + 96.0 * v4sigma4_accb + 72.0 * v4sigma4_acbb + 16.0 * v4sigma4_abbb 
                            + 72.0 * v4sigma4_aacc + 120.0 * v4sigma4_aacb + 48.0 * v4sigma4_aabb + 56.0 * v4sigma4_aaac 
                            + 48.0 * v4sigma4_aaab + 16.0 * v4sigma4_aaaa;

                double rrrr =  v4rho4_aaaa + 3.0  *  v4rho4_aaab + 3.0  *  v4rho4_aabb + v4rho4_abbb;

                double rxrr =  2.0 * v4rho3sigma_abbc + 2.0 * v4rho3sigma_abbb + 2.0 * v4rho3sigma_abba + 4.0 * v4rho3sigma_aabc 
                                + 4.0 * v4rho3sigma_aabb + 4.0 * v4rho3sigma_aaba + 2.0 * v4rho3sigma_aaac + 2.0 * v4rho3sigma_aaab + 2.0 * v4rho3sigma_aaaa;

                double rxxr = 4.0 * v4rho2sigma2_abcc + 8.0 * v4rho2sigma2_abcb + 4.0 * v4rho2sigma2_abbb + 8.0 * v4rho2sigma2_abac 
                            + 8.0 * v4rho2sigma2_abab + 4.0 * v4rho2sigma2_abaa + 4.0 * v4rho2sigma2_aacc + 8.0 * v4rho2sigma2_aacb 
                            + 4.0 * v4rho2sigma2_aabb + 8.0 * v4rho2sigma2_aaac + 8.0 * v4rho2sigma2_aaab + 4.0 * v4rho2sigma2_aaaa;

                double rxxx = 8.0 * v4rhosigma3_accc + 24.0 * v4rhosigma3_accb + 24.0 * v4rhosigma3_acbb + 8.0 * v4rhosigma3_abbb 
                            + 24.0 * v4rhosigma3_aacc + 48.0 * v4rhosigma3_aacb + 24.0 * v4rhosigma3_aabb + 24.0 * v4rhosigma3_aaac 
                            + 24.0 * v4rhosigma3_aaab + 8.0 * v4rhosigma3_aaaa; 
                            
                // Scalar contribution

                double prefac = 0.0;

                // vxc 1 contributions

                prefac += rr * rho3[g] // l1 
                        + rx * l2contract; 
                
                // // vxc 2 contributions
                
                prefac += rrr * gam[g] // q1 
                        + rxr * q2contract 
                        + rxx * q3contract 
                        + rx * q4contract; 
                
                // // // vxc 3 contributions
                prefac += rrrr * c1
                        + rxrr * c2
                        + rxxr * c3 
                        + rxr * c4 
                        + rxx * (c5_6 + c8) 
                        + rxxx * c7; 

                G_val[nu_offset + g] = w * prefac * chi_val[nu_offset + g];

                // vector contribution

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                // vxc 1 contributions

                xcomp += xr * grada_x_g * rho3[g] 
                        + x * rxw123a
                        + xx * l5contract_x;

                ycomp += xr * grada_y_g * rho3[g] 
                        + x * ryw123a
                        + xx * l5contract_y;

                zcomp += xr * grada_z_g * rho3[g] 
                        + x * rzw123a
                        + xx * l5contract_z;
                
                // // vxc 2 contributions
                
                xcomp += xrr * grada_x_g * gam[g]
                        + xr * gamx[g] 
                        + xxr * q7contract_x 
                        + xx * (q8contract_x + q10contract_x + q11contract_x) 
                        + xxx * q9contract_x; 

                ycomp += xrr * grada_y_g * gam[g] 
                        + xr * gamy[g] 
                        + xxr * q7contract_y
                        + xx * (q8contract_y + q10contract_y + q11contract_y)
                        + xxx * q9contract_y;
 
                zcomp += xrr * grada_z_g * gam[g] 
                        + xr * gamz[g] 
                        + xxr * q7contract_z
                        + xx * (q8contract_z + q10contract_z + q11contract_z)
                        + xxx * q9contract_z;

                // vxc 3 contributions
                xcomp +=  xrrr * c9_x
                        + xrr * c10_x 
                        + xxrr * c11_x
                        + xxr * (c12_c14_x + c15_x)
                        + xxxr * c13_x 
                        + xx * c17_24_25_x
                        + xxxx * c18_x
                        + xxx * (c16_19_22_x + c20_21_23_x); 

                ycomp +=  xrrr * c9_y
                        + xrr * c10_y
                        + xxrr * c11_y
                        + xxr * (c12_c14_y + c15_y)
                        + xxxr * c13_y
                        + xx * c17_24_25_y
                        + xxxx * c18_y
                        + xxx * (c16_19_22_y + c20_21_23_y);
             
                zcomp +=  xrrr * c9_z
                        + xrr * c10_z
                        + xxrr * c11_z
                        + xxr * (c12_c14_z + c15_z)
                        + xxxr * c13_z
                        + xx * c17_24_25_z
                        + xxxx * c18_z
                        + xxx * (c16_19_22_z + c20_21_23_z);


                G_gga_val[nu_offset + g] = w * (xcomp * chi_x_val[nu_offset + g] +
                                                ycomp * chi_y_val[nu_offset + g] +
                                                zcomp * chi_z_val[nu_offset + g]);
            }
        }
    }

    timer.stop("Lxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Lxc matrix matmul");

    auto mat_Lxc = denblas::multABt(gtoValues, mat_G);

    auto mat_Lxc_gga = denblas::multABt(gtoValues, mat_G_gga);

    mat_Lxc_gga.symmetrize();  // matrix + matrix.T

    mat_Lxc = denblas::addAB(mat_Lxc, mat_Lxc_gga, 1.0);

    timer.stop("Lxc matrix matmul");

    return mat_Lxc;
}

CDenseMatrix
CXCNewIntegrator::_integratePartialLxcFockForMGGA(const int32_t            npoints,
                                                  const double*            weights,
                                                  const CDenseMatrix&      gtoValues,
                                                  const CDenseMatrix&      gtoValuesX,
                                                  const CDenseMatrix&      gtoValuesY,
                                                  const CDenseMatrix&      gtoValuesZ,
                                                  const double*            rhograd,
                                                  const double*            vsigma,
                                                  const double*            v2rho2,
                                                  const double*            v2lapl2,
                                                  const double*            v2tau2,
                                                  const double*            v2rholapl,
                                                  const double*            v2rhotau,
                                                  const double*            v2lapltau,
                                                  const double*            v2rhosigma,
                                                  const double*            v2sigmalapl,
                                                  const double*            v2sigmatau,
                                                  const double*            v2sigma2,
                                                  const double*            v3rho3,
                                                  const double*            v3rho2sigma,
                                                  const double*            v3rho2lapl,
                                                  const double*            v3rho2tau,
                                                  const double*            v3rhosigma2,
                                                  const double*            v3rhosigmalapl,
                                                  const double*            v3rhosigmatau,
                                                  const double*            v3rholapl2,
                                                  const double*            v3rholapltau,
                                                  const double*            v3rhotau2,
                                                  const double*            v3sigma3,
                                                  const double*            v3sigma2lapl,
                                                  const double*            v3sigma2tau,
                                                  const double*            v3sigmalapl2,
                                                  const double*            v3sigmalapltau,
                                                  const double*            v3sigmatau2,
                                                  const double*            v3lapl3,
                                                  const double*            v3lapl2tau,
                                                  const double*            v3lapltau2,
                                                  const double*            v3tau3,
                                                  const double*            v4rho4,
                                                  const double*            v4rho3sigma,
                                                  const double*            v4rho3lapl,
                                                  const double*            v4rho3tau,
                                                  const double*            v4rho2sigma2,
                                                  const double*            v4rho2sigmalapl,
                                                  const double*            v4rho2sigmatau,
                                                  const double*            v4rho2lapl2,
                                                  const double*            v4rho2lapltau,
                                                  const double*            v4rho2tau2,
                                                  const double*            v4rhosigma3,
                                                  const double*            v4rhosigma2lapl,
                                                  const double*            v4rhosigma2tau,
                                                  const double*            v4rhosigmalapl2,
                                                  const double*            v4rhosigmalapltau,
                                                  const double*            v4rhosigmatau2,
                                                  const double*            v4rholapl3,
                                                  const double*            v4rholapl2tau,
                                                  const double*            v4rholapltau2,
                                                  const double*            v4rhotau3,
                                                  const double*            v4sigma4,
                                                  const double*            v4sigma3lapl,
                                                  const double*            v4sigma3tau,
                                                  const double*            v4sigma2lapl2,
                                                  const double*            v4sigma2lapltau,
                                                  const double*            v4sigma2tau2,
                                                  const double*            v4sigmalapl3,
                                                  const double*            v4sigmalapl2tau,
                                                  const double*            v4sigmalapltau2,
                                                  const double*            v4sigmatau3,
                                                  const double*            v4lapl4,
                                                  const double*            v4lapl3tau,
                                                  const double*            v4lapl2tau2,
                                                  const double*            v4lapltau3,
                                                  const double*            v4tau4,
                                                  const CDensityGridCubic& rwDensityGridCubic,
                                                  const CDensityGrid&      rw3DensityGrid,
                                                  const int32_t            iFock,
                                                  CMultiTimer&             timer) const
{
    timer.start("Lxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    // pointers to perturbed density gradient norms

    auto gam = rwDensityGridCubic.gam(iFock);
    auto rt_gam = rwDensityGridCubic.rt_gam(iFock);
    auto rl_gam = rwDensityGridCubic.rl_gam(iFock);
    auto tt_gam = rwDensityGridCubic.tt_gam(iFock);
    auto tl_gam = rwDensityGridCubic.tl_gam(iFock);
    auto ll_gam = rwDensityGridCubic.ll_gam(iFock);

    auto gamx = rwDensityGridCubic.gamX(iFock);
    auto gamy = rwDensityGridCubic.gamY(iFock);
    auto gamz = rwDensityGridCubic.gamZ(iFock);
    auto st_gamx = rwDensityGridCubic.st_gamX(iFock);
    auto st_gamy = rwDensityGridCubic.st_gamY(iFock);
    auto st_gamz = rwDensityGridCubic.st_gamZ(iFock);
    auto sl_gamx = rwDensityGridCubic.sl_gamX(iFock);
    auto sl_gamy = rwDensityGridCubic.sl_gamY(iFock);
    auto sl_gamz = rwDensityGridCubic.sl_gamZ(iFock);

    auto gamxx = rwDensityGridCubic.gamXX(iFock);
    auto gamxy = rwDensityGridCubic.gamXY(iFock);
    auto gamxz = rwDensityGridCubic.gamXZ(iFock);
    auto gamyx = rwDensityGridCubic.gamYX(iFock);
    auto gamyy = rwDensityGridCubic.gamYY(iFock);
    auto gamyz = rwDensityGridCubic.gamYZ(iFock);
    auto gamzx = rwDensityGridCubic.gamZX(iFock);
    auto gamzy = rwDensityGridCubic.gamZY(iFock);
    auto gamzz = rwDensityGridCubic.gamZZ(iFock);

    auto pi = rwDensityGridCubic.pi(iFock);
    auto rrt_pi = rwDensityGridCubic.rrt_pi(iFock);
    auto rrl_pi = rwDensityGridCubic.rrl_pi(iFock);
    auto rtt_pi = rwDensityGridCubic.rtt_pi(iFock);
    auto rtl_pi = rwDensityGridCubic.rtl_pi(iFock);
    auto rll_pi = rwDensityGridCubic.rll_pi(iFock);
    auto ttt_pi = rwDensityGridCubic.ttt_pi(iFock);
    auto ttl_pi = rwDensityGridCubic.ttl_pi(iFock);
    auto tll_pi = rwDensityGridCubic.tll_pi(iFock);
    auto lll_pi = rwDensityGridCubic.lll_pi(iFock);

    auto pix = rwDensityGridCubic.piX(iFock);
    auto rt_pix = rwDensityGridCubic.rt_piX(iFock);
    auto rl_pix = rwDensityGridCubic.rl_piX(iFock);
    auto ll_pix = rwDensityGridCubic.ll_piX(iFock);
    auto tt_pix = rwDensityGridCubic.tt_piX(iFock);
    auto tl_pix = rwDensityGridCubic.tl_piX(iFock);

    auto piy = rwDensityGridCubic.piY(iFock);
    auto rt_piy = rwDensityGridCubic.rt_piY(iFock);
    auto rl_piy = rwDensityGridCubic.rl_piY(iFock);
    auto ll_piy = rwDensityGridCubic.ll_piY(iFock);
    auto tt_piy = rwDensityGridCubic.tt_piY(iFock);
    auto tl_piy = rwDensityGridCubic.tl_piY(iFock);

    auto piz = rwDensityGridCubic.piZ(iFock);
    auto rt_piz = rwDensityGridCubic.rt_piZ(iFock);
    auto rl_piz = rwDensityGridCubic.rl_piZ(iFock);
    auto ll_piz = rwDensityGridCubic.ll_piZ(iFock);
    auto tt_piz = rwDensityGridCubic.tt_piZ(iFock);
    auto tl_piz = rwDensityGridCubic.tl_piZ(iFock);

    auto pixx = rwDensityGridCubic.piXX(iFock);
    auto pixy = rwDensityGridCubic.piXY(iFock);
    auto pixz = rwDensityGridCubic.piXZ(iFock);
    auto piyx = rwDensityGridCubic.piYX(iFock);
    auto piyy = rwDensityGridCubic.piYY(iFock);
    auto piyz = rwDensityGridCubic.piYZ(iFock);
    auto pizx = rwDensityGridCubic.piZX(iFock);
    auto pizy = rwDensityGridCubic.piZY(iFock);
    auto pizz = rwDensityGridCubic.piZZ(iFock);

    auto l_pixx = rwDensityGridCubic.l_piXX(iFock);
    auto l_pixy = rwDensityGridCubic.l_piXY(iFock);
    auto l_pixz = rwDensityGridCubic.l_piXZ(iFock);
    auto l_piyx = rwDensityGridCubic.l_piYX(iFock);
    auto l_piyy = rwDensityGridCubic.l_piYY(iFock);
    auto l_piyz = rwDensityGridCubic.l_piYZ(iFock);
    auto l_pizx = rwDensityGridCubic.l_piZX(iFock);
    auto l_pizy = rwDensityGridCubic.l_piZY(iFock);
    auto l_pizz = rwDensityGridCubic.l_piZZ(iFock);

    auto t_pixx = rwDensityGridCubic.t_piXX(iFock);
    auto t_pixy = rwDensityGridCubic.t_piXY(iFock);
    auto t_pixz = rwDensityGridCubic.t_piXZ(iFock);
    auto t_piyx = rwDensityGridCubic.t_piYX(iFock);
    auto t_piyy = rwDensityGridCubic.t_piYY(iFock);
    auto t_piyz = rwDensityGridCubic.t_piYZ(iFock);
    auto t_pizx = rwDensityGridCubic.t_piZX(iFock);
    auto t_pizy = rwDensityGridCubic.t_piZY(iFock);
    auto t_pizz = rwDensityGridCubic.t_piZZ(iFock);

    auto pixxx = rwDensityGridCubic.piXXX(iFock);
    auto pixxy = rwDensityGridCubic.piXXY(iFock);
    auto pixxz = rwDensityGridCubic.piXXZ(iFock);
    auto pixyx = rwDensityGridCubic.piXYX(iFock);
    auto pixyy = rwDensityGridCubic.piXYY(iFock);
    auto pixyz = rwDensityGridCubic.piXYZ(iFock);
    auto pixzx = rwDensityGridCubic.piXZX(iFock);
    auto pixzy = rwDensityGridCubic.piXZY(iFock);
    auto pixzz = rwDensityGridCubic.piXZZ(iFock);
    auto piyxx = rwDensityGridCubic.piYXX(iFock);
    auto piyxy = rwDensityGridCubic.piYXY(iFock);
    auto piyxz = rwDensityGridCubic.piYXZ(iFock);
    auto piyyx = rwDensityGridCubic.piYYX(iFock);
    auto piyyy = rwDensityGridCubic.piYYY(iFock);
    auto piyyz = rwDensityGridCubic.piYYZ(iFock);
    auto piyzx = rwDensityGridCubic.piYZX(iFock);
    auto piyzy = rwDensityGridCubic.piYZY(iFock);
    auto piyzz = rwDensityGridCubic.piYZZ(iFock);
    auto pizxx = rwDensityGridCubic.piZXX(iFock);
    auto pizxy = rwDensityGridCubic.piZXY(iFock);
    auto pizxz = rwDensityGridCubic.piZXZ(iFock);
    auto pizyx = rwDensityGridCubic.piZYX(iFock);
    auto pizyy = rwDensityGridCubic.piZYY(iFock);
    auto pizyz = rwDensityGridCubic.piZYZ(iFock);
    auto pizzx = rwDensityGridCubic.piZZX(iFock);
    auto pizzy = rwDensityGridCubic.piZZY(iFock);
    auto pizzz = rwDensityGridCubic.piZZZ(iFock);

    auto rhow = rw3DensityGrid.alphaDensity(iFock);
    auto tauw = rw3DensityGrid.alphaDensitytau(iFock);
    auto laplw = rw3DensityGrid.alphaDensitylapl(iFock);
    auto gradw_x = rw3DensityGrid.alphaDensityGradientX(iFock);
    auto gradw_y = rw3DensityGrid.alphaDensityGradientY(iFock);
    auto gradw_z = rw3DensityGrid.alphaDensityGradientZ(iFock);

    timer.stop("Lxc matrix prep.");

    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Lxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    // LDA contribution
    CDenseMatrix mat_G(naos, npoints);

    // GGA contribution
    CDenseMatrix mat_G_gga(naos, npoints);

    // tau contribution
    CDenseMatrix mat_G_gga_x(naos, npoints);
    CDenseMatrix mat_G_gga_y(naos, npoints);
    CDenseMatrix mat_G_gga_z(naos, npoints);

    auto G_val = mat_G.values();
    auto G_gga_val = mat_G_gga.values();

    auto G_gga_x_val = mat_G_gga_x.values();
    auto G_gga_y_val = mat_G_gga_y.values();
    auto G_gga_z_val = mat_G_gga_z.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;
        #pragma omp simd aligned(gam,rt_gam,rl_gam,tt_gam,tl_gam,ll_gam,gamx,gamy,gamz,\
                                 st_gamx,st_gamy,st_gamz,sl_gamx,sl_gamy,sl_gamz,gamxx,\
                                 gamxy,gamxz,gamyx,gamyy,gamyz,gamzx,gamzy,gamzz,pi,rrt_pi,\
                                 rrl_pi,rtt_pi,rtl_pi,rll_pi,ttt_pi,ttl_pi,tll_pi,lll_pi,pix,\
                                 rt_pix,rl_pix,ll_pix,tt_pix,tl_pix,piy,rt_piy,rl_piy,ll_piy,\
                                 tt_piy,tl_piy,piz,rt_piz,rl_piz,ll_piz,tt_piz,tl_piz,pixx,pixy,\
                                 pixz,piyx,piyy,piyz,pizx,pizy,pizz,l_pixx,l_pixy,l_pixz,l_piyx,l_piyy,\
                                 l_piyz,l_pizx,l_pizy,l_pizz,t_pixx,t_pixy,t_pixz,t_piyx,t_piyy,t_piyz,\
                                 t_pizx,t_pizy,t_pizz,pixxx,pixxy,pixxz,pixyx,pixyy,pixyz,pixzx,pixzy,\
                                 pixzz,piyxx,piyxy,piyxz,piyyx,piyyy,piyyz,piyzx,piyzy,piyzz,pizxx,pizxy,\
                                 pizxz,pizyx,pizyy,pizyz,pizzx,pizzy,pizzz,rhow,tauw,laplw, gradw_x,gradw_y,\
                                 gradw_z ,chi_x_val ,chi_y_val ,chi_z_val ,\
                                 rhograd,vsigma ,v2rho2,v2lapl2 ,v2tau2 ,v2rholapl ,v2rhotau,v2lapltau ,v2rhosigma,\
                                 v2sigmalapl ,v2sigmatau,v2sigma2,v3rho3,v3rho2sigma,v3rho2lapl,v3rho2tau,v3rhosigma2,\
                                 v3rhosigmalapl,v3rhosigmatau,v3rholapl2,v3rholapltau,v3rhotau2,v3sigma3,v3sigma2lapl,\
                                 v3sigma2tau,v3sigmalapl2,v3sigmalapltau,v3sigmatau2,v3lapl3,v3lapl2tau,v3lapltau2,v3tau3,\
                                 v4rho4,v4rho3sigma,v4rho3lapl,v4rho3tau,v4rho2sigma2,v4rho2sigmalapl,v4rho2sigmatau,v4rho2lapl2,\
                                 v4rho2lapltau,v4rho2tau2,v4rhosigma3,v4rhosigma2lapl,v4rhosigma2tau,v4rhosigmalapl2,v4rhosigmalapltau,\
                                 v4rhosigmatau2,v4rholapl3,v4rholapl2tau,v4rholapltau2,v4rhotau3,v4sigma4,v4sigma3lapl,v4sigma3tau,\
                                 v4sigma2lapl2,v4sigma2lapltau,v4sigma2tau2,v4sigmalapl3,v4sigmalapl2tau,v4sigmalapltau2,v4sigmatau3,\
                                 v4lapl4,v4lapl3tau,v4lapl2tau2,v4lapltau3,v4tau4:VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                double w = weights[g];

                double grada_x_g = rhograd[6 * g + 0];
                double grada_y_g = rhograd[6 * g + 1];
                double grada_z_g = rhograd[6 * g + 2];

                // vx1 contractions

                double l2contract =   grada_x_g * gradw_x[g] + grada_y_g * gradw_y[g] + grada_z_g * gradw_z[g];

                double l5contract_x = grada_x_g * l2contract;
                double l5contract_y = grada_y_g * l2contract;
                double l5contract_z = grada_z_g * l2contract;

                // vx2 contractions

                double q2contract = grada_x_g * gamx[g] + grada_y_g * gamy[g] + grada_z_g * gamz[g];
                double sl_q2contract = grada_x_g * sl_gamx[g] + grada_y_g * sl_gamy[g] + grada_z_g * sl_gamz[g];
                double st_q2contract = grada_x_g * st_gamx[g] + grada_y_g * st_gamy[g] + grada_z_g * st_gamz[g];

                double q3contract =   grada_x_g * grada_x_g * gamxx[g]
                                    + grada_x_g * grada_y_g * gamxy[g]
                                    + grada_x_g * grada_z_g * gamxz[g]
                                    + grada_y_g * grada_x_g * gamyx[g]
                                    + grada_y_g * grada_y_g * gamyy[g]
                                    + grada_y_g * grada_z_g * gamyz[g]
                                    + grada_z_g * grada_x_g * gamzx[g]
                                    + grada_z_g * grada_y_g * gamzy[g]
                                    + grada_z_g * grada_z_g * gamzz[g];

                double q4contract = gamxx[g] + gamyy[g] + gamzz[g];

                double q7contract_x =  grada_x_g * grada_x_g *  gamx[g] + grada_x_g * grada_y_g *  gamy[g] + grada_x_g * grada_z_g *  gamz[g];
                double q7contract_y =  grada_y_g * grada_x_g *  gamx[g] + grada_y_g * grada_y_g *  gamy[g] + grada_y_g * grada_z_g *  gamz[g];
                double q7contract_z =  grada_z_g * grada_x_g *  gamx[g] + grada_z_g * grada_y_g *  gamy[g] + grada_z_g * grada_z_g *  gamz[g];

                double sl_q7contract_x =  grada_x_g * grada_x_g * sl_gamx[g] + grada_x_g * grada_y_g * sl_gamy[g] + grada_x_g * grada_z_g * sl_gamz[g];
                double sl_q7contract_y =  grada_y_g * grada_x_g * sl_gamx[g] + grada_y_g * grada_y_g * sl_gamy[g] + grada_y_g * grada_z_g * sl_gamz[g];
                double sl_q7contract_z =  grada_z_g * grada_x_g * sl_gamx[g] + grada_z_g * grada_y_g * sl_gamy[g] + grada_z_g * grada_z_g * sl_gamz[g];

                double st_q7contract_x =  grada_x_g * grada_x_g * st_gamx[g] + grada_x_g * grada_y_g * st_gamy[g] + grada_x_g * grada_z_g * st_gamz[g];
                double st_q7contract_y =  grada_y_g * grada_x_g * st_gamx[g] + grada_y_g * grada_y_g * st_gamy[g] + grada_y_g * grada_z_g * st_gamz[g];
                double st_q7contract_z =  grada_z_g * grada_x_g * st_gamx[g] + grada_z_g * grada_y_g * st_gamy[g] + grada_z_g * grada_z_g * st_gamz[g];

                double q8contract_x =  grada_x_g *  gamxx[g] + grada_y_g *  gamxy[g] + grada_z_g *  gamxz[g];
                double q8contract_y =  grada_x_g *  gamyx[g] + grada_y_g *  gamyy[g] + grada_z_g *  gamyz[g];
                double q8contract_z =  grada_x_g *  gamzx[g] + grada_y_g *  gamzy[g] + grada_z_g *  gamzz[g];

                double q9contract_x =  grada_x_g *  q3contract;
                double q9contract_y =  grada_y_g *  q3contract;
                double q9contract_z =  grada_z_g *  q3contract;

                double q10contract_x =  grada_x_g *  gamxx[g] + grada_y_g *  gamyx[g] + grada_z_g *  gamzx[g];
                double q10contract_y =  grada_x_g *  gamxy[g] + grada_y_g *  gamyy[g] + grada_z_g *  gamzy[g];
                double q10contract_z =  grada_x_g *  gamxz[g] + grada_y_g *  gamyz[g] + grada_z_g *  gamzz[g];

                double q11contract_x =  grada_x_g *  gamxx[g] + grada_x_g *  gamyy[g] + grada_x_g *  gamzz[g];
                double q11contract_y =  grada_y_g *  gamxx[g] + grada_y_g *  gamyy[g] + grada_y_g *  gamzz[g];
                double q11contract_z =  grada_z_g *  gamxx[g] + grada_z_g *  gamyy[g] + grada_z_g *  gamzz[g];

                // vx3 contractions

                double    c2 = grada_x_g *    pix[g] + grada_y_g *    piy[g] + grada_z_g *    piz[g];
                double rt_c2 = grada_x_g * rt_pix[g] + grada_y_g * rt_piy[g] + grada_z_g * rt_piz[g];
                double rl_c2 = grada_x_g * rl_pix[g] + grada_y_g * rl_piy[g] + grada_z_g * rl_piz[g];
                double ll_c2 = grada_x_g * ll_pix[g] + grada_y_g * ll_piy[g] + grada_z_g * ll_piz[g];
                double tt_c2 = grada_x_g * tt_pix[g] + grada_y_g * tt_piy[g] + grada_z_g * tt_piz[g];
                double tl_c2 = grada_x_g * tl_pix[g] + grada_y_g * tl_piy[g] + grada_z_g * tl_piz[g];

                double c3 =  grada_x_g * grada_x_g * pixx[g]
                           + grada_x_g * grada_y_g * pixy[g]
                           + grada_x_g * grada_z_g * pixz[g]
                           + grada_y_g * grada_x_g * piyx[g]
                           + grada_y_g * grada_y_g * piyy[g]
                           + grada_y_g * grada_z_g * piyz[g]
                           + grada_z_g * grada_x_g * pizx[g]
                           + grada_z_g * grada_y_g * pizy[g]
                           + grada_z_g * grada_z_g * pizz[g];

                double l_c3 =  grada_x_g * grada_x_g * l_pixx[g]
                             + grada_x_g * grada_y_g * l_pixy[g]
                             + grada_x_g * grada_z_g * l_pixz[g]
                             + grada_y_g * grada_x_g * l_piyx[g]
                             + grada_y_g * grada_y_g * l_piyy[g]
                             + grada_y_g * grada_z_g * l_piyz[g]
                             + grada_z_g * grada_x_g * l_pizx[g]
                             + grada_z_g * grada_y_g * l_pizy[g]
                             + grada_z_g * grada_z_g * l_pizz[g];

                double t_c3 =  grada_x_g * grada_x_g * t_pixx[g]
                             + grada_x_g * grada_y_g * t_pixy[g]
                             + grada_x_g * grada_z_g * t_pixz[g]
                             + grada_y_g * grada_x_g * t_piyx[g]
                             + grada_y_g * grada_y_g * t_piyy[g]
                             + grada_y_g * grada_z_g * t_piyz[g]
                             + grada_z_g * grada_x_g * t_pizx[g]
                             + grada_z_g * grada_y_g * t_pizy[g]
                             + grada_z_g * grada_z_g * t_pizz[g];

                double   c4 =   pixx[g] +   piyy[g] +   pizz[g];
                double l_c4 = l_pixx[g] + l_piyy[g] + l_pizz[g];
                double t_c4 = t_pixx[g] + t_piyy[g] + t_pizz[g];

                double c5_6 = grada_x_g * (pixxx[g] + pixxx[g])
                            + grada_x_g * (piyxy[g] + pixyy[g])
                            + grada_x_g * (pizxz[g] + pixzz[g])
                            + grada_y_g * (pixyx[g] + piyxx[g])
                            + grada_y_g * (piyyy[g] + piyyy[g])
                            + grada_y_g * (pizyz[g] + piyzz[g])
                            + grada_z_g * (pixzx[g] + pizxx[g])
                            + grada_z_g * (piyzy[g] + pizyy[g])
                            + grada_z_g * (pizzz[g] + pizzz[g]);

                double c7 =   grada_x_g * grada_x_g * grada_x_g * pixxx[g]
                            + grada_x_g * grada_x_g * grada_y_g * pixxy[g]
                            + grada_x_g * grada_x_g * grada_z_g * pixxz[g]
                            + grada_x_g * grada_y_g * grada_x_g * pixyx[g]
                            + grada_x_g * grada_y_g * grada_y_g * pixyy[g]
                            + grada_x_g * grada_y_g * grada_z_g * pixyz[g]
                            + grada_x_g * grada_z_g * grada_x_g * pixzx[g]
                            + grada_x_g * grada_z_g * grada_y_g * pixzy[g]
                            + grada_x_g * grada_z_g * grada_z_g * pixzz[g]
                            + grada_y_g * grada_x_g * grada_x_g * piyxx[g]
                            + grada_y_g * grada_x_g * grada_y_g * piyxy[g]
                            + grada_y_g * grada_x_g * grada_z_g * piyxz[g]
                            + grada_y_g * grada_y_g * grada_x_g * piyyx[g]
                            + grada_y_g * grada_y_g * grada_y_g * piyyy[g]
                            + grada_y_g * grada_y_g * grada_z_g * piyyz[g]
                            + grada_y_g * grada_z_g * grada_x_g * piyzx[g]
                            + grada_y_g * grada_z_g * grada_y_g * piyzy[g]
                            + grada_y_g * grada_z_g * grada_z_g * piyzz[g]
                            + grada_z_g * grada_x_g * grada_x_g * pizxx[g]
                            + grada_z_g * grada_x_g * grada_y_g * pizxy[g]
                            + grada_z_g * grada_x_g * grada_z_g * pizxz[g]
                            + grada_z_g * grada_y_g * grada_x_g * pizyx[g]
                            + grada_z_g * grada_y_g * grada_y_g * pizyy[g]
                            + grada_z_g * grada_y_g * grada_z_g * pizyz[g]
                            + grada_z_g * grada_z_g * grada_x_g * pizzx[g]
                            + grada_z_g * grada_z_g * grada_y_g * pizzy[g]
                            + grada_z_g * grada_z_g * grada_z_g * pizzz[g];

                double c8 =  grada_x_g * pixxx[g]
                           + grada_y_g * pixxy[g]
                           + grada_z_g * pixxz[g]
                           + grada_x_g * piyyx[g]
                           + grada_y_g * piyyy[g]
                           + grada_z_g * piyyz[g]
                           + grada_x_g * pizzx[g]
                           + grada_y_g * pizzy[g]
                           + grada_z_g * pizzz[g];

                double c9_x = grada_x_g * pi[g];
                double c9_y = grada_y_g * pi[g];
                double c9_z = grada_z_g * pi[g];

                double c10_x = pix[g];
                double c10_y = piy[g];
                double c10_z = piz[g];

                double c11_x = c2 * grada_x_g;
                double c11_y = c2 * grada_y_g;
                double c11_z = c2 * grada_z_g;

                double c12_c14_x =   grada_x_g * (pixx[g] + pixx[g])
                                   + grada_y_g * (pixy[g] + piyx[g])
                                   + grada_z_g * (pixz[g] + pizx[g]);

                double l_c12_c14_x =   grada_x_g * (l_pixx[g] + l_pixx[g])
                                     + grada_y_g * (l_pixy[g] + l_piyx[g])
                                     + grada_z_g * (l_pixz[g] + l_pizx[g]);

                double t_c12_c14_x =   grada_x_g * (t_pixx[g] + t_pixx[g])
                                     + grada_y_g * (t_pixy[g] + t_piyx[g])
                                     + grada_z_g * (t_pixz[g] + t_pizx[g]);

                double c12_c14_y=  grada_x_g * (piyx[g] + pixy[g])
                                 + grada_y_g * (piyy[g] + piyy[g])
                                 + grada_z_g * (piyz[g] + pizy[g]);

                double l_c12_c14_y=  grada_x_g * (l_piyx[g] + l_pixy[g])
                                   + grada_y_g * (l_piyy[g] + l_piyy[g])
                                   + grada_z_g * (l_piyz[g] + l_pizy[g]);

                double t_c12_c14_y=  grada_x_g * (t_piyx[g] + t_pixy[g])
                                   + grada_y_g * (t_piyy[g] + t_piyy[g])
                                   + grada_z_g * (t_piyz[g] + t_pizy[g]);

                double c12_c14_z=  grada_x_g * (pizx[g] + pixz[g])
                                 + grada_y_g * (pizy[g] + piyz[g])
                                 + grada_z_g * (pizz[g] + pizz[g]);

                double l_c12_c14_z=  grada_x_g * (l_pizx[g] + l_pixz[g])
                                   + grada_y_g * (l_pizy[g] + l_piyz[g])
                                   + grada_z_g * (l_pizz[g] + l_pizz[g]);

                double t_c12_c14_z=  grada_x_g * (t_pizx[g] + t_pixz[g])
                                   + grada_y_g * (t_pizy[g] + t_piyz[g])
                                   + grada_z_g * (t_pizz[g] + t_pizz[g]);

                double c13 = grada_x_g * grada_x_g * pixx[g]
                           + grada_x_g * grada_y_g * pixy[g]
                           + grada_x_g * grada_z_g * pixz[g]
                           + grada_y_g * grada_x_g * piyx[g]
                           + grada_y_g * grada_y_g * piyy[g]
                           + grada_y_g * grada_z_g * piyz[g]
                           + grada_z_g * grada_x_g * pizx[g]
                           + grada_z_g * grada_y_g * pizy[g]
                           + grada_z_g * grada_z_g * pizz[g];

                double l_c13 = grada_x_g * grada_x_g * l_pixx[g]
                             + grada_x_g * grada_y_g * l_pixy[g]
                             + grada_x_g * grada_z_g * l_pixz[g]
                             + grada_y_g * grada_x_g * l_piyx[g]
                             + grada_y_g * grada_y_g * l_piyy[g]
                             + grada_y_g * grada_z_g * l_piyz[g]
                             + grada_z_g * grada_x_g * l_pizx[g]
                             + grada_z_g * grada_y_g * l_pizy[g]
                             + grada_z_g * grada_z_g * l_pizz[g];

                double t_c13 = grada_x_g * grada_x_g * t_pixx[g]
                             + grada_x_g * grada_y_g * t_pixy[g]
                             + grada_x_g * grada_z_g * t_pixz[g]
                             + grada_y_g * grada_x_g * t_piyx[g]
                             + grada_y_g * grada_y_g * t_piyy[g]
                             + grada_y_g * grada_z_g * t_piyz[g]
                             + grada_z_g * grada_x_g * t_pizx[g]
                             + grada_z_g * grada_y_g * t_pizy[g]
                             + grada_z_g * grada_z_g * t_pizz[g];

                double   c13_x =   c13 * grada_x_g;
                double l_c13_x = l_c13 * grada_x_g;
                double t_c13_x = t_c13 * grada_x_g;

                double   c13_y = c13 * grada_y_g;
                double l_c13_y = l_c13 * grada_y_g;
                double t_c13_y = t_c13 * grada_y_g;

                double   c13_z = c13 * grada_z_g;
                double l_c13_z = l_c13 * grada_z_g;
                double t_c13_z = t_c13 * grada_z_g;

                double   c15_x = grada_x_g * c4;
                double l_c15_x = grada_x_g * l_c4;
                double t_c15_x = grada_x_g * t_c4;

                double   c15_y = grada_y_g * c4;
                double l_c15_y = grada_y_g * l_c4;
                double t_c15_y = grada_y_g * t_c4;

                double   c15_z = grada_z_g * c4;
                double l_c15_z = grada_z_g * l_c4;
                double t_c15_z = grada_z_g * t_c4;

                double c16_19_22_x = grada_x_g * grada_x_g * pixxx[g] + grada_x_g * grada_x_g * pixxx[g] + grada_x_g * grada_x_g * pixxx[g]
                                    +grada_x_g * grada_y_g * pixxy[g] + grada_x_g * grada_y_g * pixyx[g] + grada_x_g * grada_y_g * pixxy[g]
                                    +grada_x_g * grada_z_g * pixxz[g] + grada_x_g * grada_z_g * pixzx[g] + grada_x_g * grada_z_g * pixxz[g]
                                    +grada_y_g * grada_x_g * pixyx[g] + grada_y_g * grada_x_g * piyxx[g] + grada_y_g * grada_x_g * piyxx[g]
                                    +grada_y_g * grada_y_g * pixyy[g] + grada_y_g * grada_y_g * piyyx[g] + grada_y_g * grada_y_g * piyxy[g]
                                    +grada_y_g * grada_z_g * pixyz[g] + grada_y_g * grada_z_g * piyzx[g] + grada_y_g * grada_z_g * piyxz[g]
                                    +grada_z_g * grada_x_g * pixzx[g] + grada_z_g * grada_x_g * pizxx[g] + grada_z_g * grada_x_g * pizxx[g]
                                    +grada_z_g * grada_y_g * pixzy[g] + grada_z_g * grada_y_g * pizyx[g] + grada_z_g * grada_y_g * pizxy[g]
                                    +grada_z_g * grada_z_g * pixzz[g] + grada_z_g * grada_z_g * pizzx[g] + grada_z_g * grada_z_g * pizxz[g];

                double c16_19_22_y = grada_x_g * grada_x_g * piyxx[g] + grada_x_g * grada_x_g * pixxy[g] + grada_x_g * grada_x_g * pixyx[g]
                                    +grada_x_g * grada_y_g * piyxy[g] + grada_x_g * grada_y_g * pixyy[g] + grada_x_g * grada_y_g * pixyy[g]
                                    +grada_x_g * grada_z_g * piyxz[g] + grada_x_g * grada_z_g * pixzy[g] + grada_x_g * grada_z_g * pixyz[g]
                                    +grada_y_g * grada_x_g * piyyx[g] + grada_y_g * grada_x_g * piyxy[g] + grada_y_g * grada_x_g * piyyx[g]
                                    +grada_y_g * grada_y_g * piyyy[g] + grada_y_g * grada_y_g * piyyy[g] + grada_y_g * grada_y_g * piyyy[g]
                                    +grada_y_g * grada_z_g * piyyz[g] + grada_y_g * grada_z_g * piyzy[g] + grada_y_g * grada_z_g * piyyz[g]
                                    +grada_z_g * grada_x_g * piyzx[g] + grada_z_g * grada_x_g * pizxy[g] + grada_z_g * grada_x_g * pizyx[g]
                                    +grada_z_g * grada_y_g * piyzy[g] + grada_z_g * grada_y_g * pizyy[g] + grada_z_g * grada_y_g * pizyy[g]
                                    +grada_z_g * grada_z_g * piyzz[g] + grada_z_g * grada_z_g * pizzy[g] + grada_z_g * grada_z_g * pizyz[g];

                double c16_19_22_z = grada_x_g * grada_x_g * pizxx[g] + grada_x_g * grada_x_g * pixxz[g] + grada_x_g * grada_x_g * pixzx[g]
                                    +grada_x_g * grada_y_g * pizxy[g] + grada_x_g * grada_y_g * pixyz[g] + grada_x_g * grada_y_g * pixzy[g]
                                    +grada_x_g * grada_z_g * pizxz[g] + grada_x_g * grada_z_g * pixzz[g] + grada_x_g * grada_z_g * pixzz[g]
                                    +grada_y_g * grada_x_g * pizyx[g] + grada_y_g * grada_x_g * piyxz[g] + grada_y_g * grada_x_g * piyzx[g]
                                    +grada_y_g * grada_y_g * pizyy[g] + grada_y_g * grada_y_g * piyyz[g] + grada_y_g * grada_y_g * piyzy[g]
                                    +grada_y_g * grada_z_g * pizyz[g] + grada_y_g * grada_z_g * piyzz[g] + grada_y_g * grada_z_g * piyzz[g]
                                    +grada_z_g * grada_x_g * pizzx[g] + grada_z_g * grada_x_g * pizxz[g] + grada_z_g * grada_x_g * pizzx[g]
                                    +grada_z_g * grada_y_g * pizzy[g] + grada_z_g * grada_y_g * pizyz[g] + grada_z_g * grada_y_g * pizzy[g]
                                    +grada_z_g * grada_z_g * pizzz[g] + grada_z_g * grada_z_g * pizzz[g] + grada_z_g * grada_z_g * pizzz[g];

                double c17_24_25_x =  pixxx[g] + pixxx[g] + pixxx[g]
                                    + pixyy[g] + piyxy[g] + piyyx[g]
                                    + pixzz[g] + pizxz[g] + pizzx[g];

                double c17_24_25_y =  piyxx[g] + pixyx[g] + pixxy[g]
                                    + piyyy[g] + piyyy[g] + piyyy[g]
                                    + piyzz[g] + pizyz[g] + pizzy[g];

                double c17_24_25_z =  pizxx[g] + pixzx[g] + pixxz[g]
                                    + pizyy[g] + piyzy[g] + piyyz[g]
                                    + pizzz[g] + pizzz[g] + pizzz[g];

                double c18_x = c7 * grada_x_g;
                double c18_y = c7 * grada_y_g;
                double c18_z = c7 * grada_z_g;

                double c20_21_23 =   grada_x_g * (pixxx[g] + pixxx[g] + pixxx[g])
                                   + grada_x_g * (piyxy[g] + pixyy[g] + piyyx[g])
                                   + grada_x_g * (pizxz[g] + pixzz[g] + pizzx[g])
                                   + grada_y_g * (pixyx[g] + piyxx[g] + pixxy[g])
                                   + grada_y_g * (piyyy[g] + piyyy[g] + piyyy[g])
                                   + grada_y_g * (pizyz[g] + piyzz[g] + pizzy[g])
                                   + grada_z_g * (pixzx[g] + pizxx[g] + pixxz[g])
                                   + grada_z_g * (piyzy[g] + pizyy[g] + piyyz[g])
                                   + grada_z_g * (pizzz[g] + pizzz[g] + pizzz[g]);

                double c20_21_23_x = grada_x_g * c20_21_23;
                double c20_21_23_y = grada_y_g * c20_21_23;
                double c20_21_23_z = grada_z_g * c20_21_23;

                // functional derivatives in libxc form

                // first-order
                double vsigma_a = vsigma[3 * g + 0];
                double vsigma_c = vsigma[3 * g + 1];

                // second-order
                double v2rho2_aa = v2rho2[3 * g + 0];
                double v2rho2_ab = v2rho2[3 * g + 1];

                double v2rhosigma_aa = v2rhosigma[6 * g + 0];
                double v2rhosigma_ac = v2rhosigma[6 * g + 1];
                double v2rhosigma_ab = v2rhosigma[6 * g + 2];
                double v2rhosigma_ba = v2rhosigma[6 * g + 3];
                double v2rhosigma_bc = v2rhosigma[6 * g + 4];

                double v2sigma2_aa = v2sigma2[6 * g + 0];
                double v2sigma2_ac = v2sigma2[6 * g + 1];
                double v2sigma2_ab = v2sigma2[6 * g + 2];
                double v2sigma2_cc = v2sigma2[6 * g + 3];
                double v2sigma2_cb = v2sigma2[6 * g + 4];

                double v2rholapl_aa = v2rholapl[4 * g + 0];
                double v2rholapl_ab = v2rholapl[4 * g + 1];
                double v2rholapl_ba = v2rholapl[4 * g + 2];

                double v2rhotau_aa = v2rhotau[4 * g + 0];
                double v2rhotau_ab = v2rhotau[4 * g + 1];
                double v2rhotau_ba = v2rhotau[4 * g + 2];

                double v2lapltau_aa = v2lapltau[4 * g + 0];
                double v2lapltau_ab = v2lapltau[4 * g + 1];
                double v2lapltau_ba = v2lapltau[4 * g + 2];

                double v2lapl2_aa = v2lapl2[3 * g + 0];
                double v2lapl2_ab = v2lapl2[3 * g + 1];

                double v2tau2_aa = v2tau2[3 * g + 0];
                double v2tau2_ab = v2tau2[3 * g + 1];

                double v2sigmalapl_aa = v2sigmalapl[6 * g + 0];
                double v2sigmalapl_ab = v2sigmalapl[6 * g + 1];
                double v2sigmalapl_ca = v2sigmalapl[6 * g + 2];
                double v2sigmalapl_cb = v2sigmalapl[6 * g + 3];
                double v2sigmalapl_ba = v2sigmalapl[6 * g + 4];

                double v2sigmatau_aa = v2sigmatau[6 * g + 0];
                double v2sigmatau_ab = v2sigmatau[6 * g + 1];
                double v2sigmatau_ca = v2sigmatau[6 * g + 2];
                double v2sigmatau_cb = v2sigmatau[6 * g + 3];
                double v2sigmatau_ba = v2sigmatau[6 * g + 4];

                // Third-order terms

                auto v3rho3_aaa = v3rho3[4 * g + 0];
                auto v3rho3_aab = v3rho3[4 * g + 1];
                auto v3rho3_abb = v3rho3[4 * g + 2];

                auto v3rho2lapl_aaa = v3rho2lapl[6 * g + 0];
                auto v3rho2lapl_aab = v3rho2lapl[6 * g + 1];
                auto v3rho2lapl_aba = v3rho2lapl[6 * g + 2];
                auto v3rho2lapl_abb = v3rho2lapl[6 * g + 3];
                auto v3rho2lapl_bba = v3rho2lapl[6 * g + 4];

                auto v3rho2tau_aaa = v3rho2tau[6 * g + 0];
                auto v3rho2tau_aab = v3rho2tau[6 * g + 1];
                auto v3rho2tau_aba = v3rho2tau[6 * g + 2];
                auto v3rho2tau_abb = v3rho2tau[6 * g + 3];
                auto v3rho2tau_bba = v3rho2tau[6 * g + 4];

                auto v3rholapl2_aaa = v3rholapl2[6 * g + 0];
                auto v3rholapl2_aab = v3rholapl2[6 * g + 1];
                auto v3rholapl2_abb = v3rholapl2[6 * g + 2];
                auto v3rholapl2_baa = v3rholapl2[6 * g + 3];
                auto v3rholapl2_bab = v3rholapl2[6 * g + 4];

                auto v3rholapltau_aaa = v3rholapltau[8 * g + 0];
                auto v3rholapltau_aab = v3rholapltau[8 * g + 1];
                auto v3rholapltau_aba = v3rholapltau[8 * g + 2];
                auto v3rholapltau_abb = v3rholapltau[8 * g + 3];
                auto v3rholapltau_baa = v3rholapltau[8 * g + 4];
                auto v3rholapltau_bab = v3rholapltau[8 * g + 5];
                auto v3rholapltau_bba = v3rholapltau[8 * g + 6];

                auto v3rhotau2_aaa = v3rhotau2[6 * g + 0];
                auto v3rhotau2_aab = v3rhotau2[6 * g + 1];
                auto v3rhotau2_abb = v3rhotau2[6 * g + 2];
                auto v3rhotau2_baa = v3rhotau2[6 * g + 3];
                auto v3rhotau2_bab = v3rhotau2[6 * g + 4];

                auto v3sigma2lapl_aaa = v3sigma2lapl[12 * g + 0];
                auto v3sigma2lapl_aab = v3sigma2lapl[12 * g + 1];
                auto v3sigma2lapl_aca = v3sigma2lapl[12 * g + 2];
                auto v3sigma2lapl_acb = v3sigma2lapl[12 * g + 3];
                auto v3sigma2lapl_aba = v3sigma2lapl[12 * g + 4];
                auto v3sigma2lapl_abb = v3sigma2lapl[12 * g + 5];
                auto v3sigma2lapl_cca = v3sigma2lapl[12 * g + 6];
                auto v3sigma2lapl_ccb = v3sigma2lapl[12 * g + 7];
                auto v3sigma2lapl_cba = v3sigma2lapl[12 * g + 8];
                auto v3sigma2lapl_cbb = v3sigma2lapl[12 * g + 9];
                auto v3sigma2lapl_bba = v3sigma2lapl[12 * g + 10];

                auto v3rhosigmalapl_aaa = v3rhosigmalapl[12 * g + 0];
                auto v3rhosigmalapl_aab = v3rhosigmalapl[12 * g + 1];
                auto v3rhosigmalapl_aca = v3rhosigmalapl[12 * g + 2];
                auto v3rhosigmalapl_acb = v3rhosigmalapl[12 * g + 3];
                auto v3rhosigmalapl_aba = v3rhosigmalapl[12 * g + 4];
                auto v3rhosigmalapl_abb = v3rhosigmalapl[12 * g + 5];
                auto v3rhosigmalapl_baa = v3rhosigmalapl[12 * g + 6];
                auto v3rhosigmalapl_bab = v3rhosigmalapl[12 * g + 7];
                auto v3rhosigmalapl_bca = v3rhosigmalapl[12 * g + 8];
                auto v3rhosigmalapl_bcb = v3rhosigmalapl[12 * g + 9];
                auto v3rhosigmalapl_bba = v3rhosigmalapl[12 * g + 10];

                auto v3rhosigmatau_aaa = v3rhosigmatau[12 * g + 0];
                auto v3rhosigmatau_aab = v3rhosigmatau[12 * g + 1];
                auto v3rhosigmatau_aca = v3rhosigmatau[12 * g + 2];
                auto v3rhosigmatau_acb = v3rhosigmatau[12 * g + 3];
                auto v3rhosigmatau_aba = v3rhosigmatau[12 * g + 4];
                auto v3rhosigmatau_abb = v3rhosigmatau[12 * g + 5];
                auto v3rhosigmatau_baa = v3rhosigmatau[12 * g + 6];
                auto v3rhosigmatau_bab = v3rhosigmatau[12 * g + 7];
                auto v3rhosigmatau_bca = v3rhosigmatau[12 * g + 8];
                auto v3rhosigmatau_bcb = v3rhosigmatau[12 * g + 9];
                auto v3rhosigmatau_bba = v3rhosigmatau[12 * g + 10];

                auto v3sigma2tau_aaa = v3sigma2tau[12 * g + 0];
                auto v3sigma2tau_aab = v3sigma2tau[12 * g + 1];
                auto v3sigma2tau_aca = v3sigma2tau[12 * g + 2];
                auto v3sigma2tau_acb = v3sigma2tau[12 * g + 3];
                auto v3sigma2tau_aba = v3sigma2tau[12 * g + 4];
                auto v3sigma2tau_abb = v3sigma2tau[12 * g + 5];
                auto v3sigma2tau_cca = v3sigma2tau[12 * g + 6];
                auto v3sigma2tau_ccb = v3sigma2tau[12 * g + 7];
                auto v3sigma2tau_cba = v3sigma2tau[12 * g + 8];
                auto v3sigma2tau_cbb = v3sigma2tau[12 * g + 9];
                auto v3sigma2tau_bba = v3sigma2tau[12 * g + 10];

                auto v3sigmalapl2_aaa = v3sigmalapl2[9 * g + 0];
                auto v3sigmalapl2_aab = v3sigmalapl2[9 * g + 1];
                auto v3sigmalapl2_abb = v3sigmalapl2[9 * g + 2];
                auto v3sigmalapl2_caa = v3sigmalapl2[9 * g + 3];
                auto v3sigmalapl2_cab = v3sigmalapl2[9 * g + 4];
                auto v3sigmalapl2_cbb = v3sigmalapl2[9 * g + 5];
                auto v3sigmalapl2_baa = v3sigmalapl2[9 * g + 6];
                auto v3sigmalapl2_bab = v3sigmalapl2[9 * g + 7];

                auto v3sigmalapltau_aaa = v3sigmalapltau[12 * g + 0];
                auto v3sigmalapltau_aab = v3sigmalapltau[12 * g + 1];
                auto v3sigmalapltau_aba = v3sigmalapltau[12 * g + 2];
                auto v3sigmalapltau_abb = v3sigmalapltau[12 * g + 3];
                auto v3sigmalapltau_caa = v3sigmalapltau[12 * g + 4];
                auto v3sigmalapltau_cab = v3sigmalapltau[12 * g + 5];
                auto v3sigmalapltau_cba = v3sigmalapltau[12 * g + 6];
                auto v3sigmalapltau_cbb = v3sigmalapltau[12 * g + 7];
                auto v3sigmalapltau_baa = v3sigmalapltau[12 * g + 8];
                auto v3sigmalapltau_bab = v3sigmalapltau[12 * g + 9];
                auto v3sigmalapltau_bba = v3sigmalapltau[12 * g + 10];

                auto v3sigmatau2_aaa = v3sigmatau2[9 * g + 0];
                auto v3sigmatau2_aab = v3sigmatau2[9 * g + 1];
                auto v3sigmatau2_abb = v3sigmatau2[9 * g + 2];
                auto v3sigmatau2_caa = v3sigmatau2[9 * g + 3];
                auto v3sigmatau2_cab = v3sigmatau2[9 * g + 4];
                auto v3sigmatau2_cbb = v3sigmatau2[9 * g + 5];
                auto v3sigmatau2_baa = v3sigmatau2[9 * g + 6];
                auto v3sigmatau2_bab = v3sigmatau2[9 * g + 7];

                auto v3lapl3_aaa = v3lapl3[4 * g + 0];
                auto v3lapl3_aab = v3lapl3[4 * g + 1];
                auto v3lapl3_abb = v3lapl3[4 * g + 2];

                auto v3lapl2tau_aaa = v3lapl2tau[6 * g + 0];
                auto v3lapl2tau_aab = v3lapl2tau[6 * g + 1];
                auto v3lapl2tau_aba = v3lapl2tau[6 * g + 2];
                auto v3lapl2tau_abb = v3lapl2tau[6 * g + 3];
                auto v3lapl2tau_bba = v3lapl2tau[6 * g + 4];

                auto v3lapltau2_aaa = v3lapltau2[6 * g + 0];
                auto v3lapltau2_aab = v3lapltau2[6 * g + 1];
                auto v3lapltau2_abb = v3lapltau2[6 * g + 2];
                auto v3lapltau2_baa = v3lapltau2[6 * g + 3];
                auto v3lapltau2_bab = v3lapltau2[6 * g + 4];

                auto v3tau3_aaa = v3tau3[4 * g + 0];
                auto v3tau3_aab = v3tau3[4 * g + 1];
                auto v3tau3_abb = v3tau3[4 * g + 2];
                auto v3tau3_bbb = v3tau3[4 * g + 3];

                auto v3rho2sigma_aaa = v3rho2sigma[9 * g + 0];
                auto v3rho2sigma_aac = v3rho2sigma[9 * g + 1];
                auto v3rho2sigma_aab = v3rho2sigma[9 * g + 2];
                auto v3rho2sigma_aba = v3rho2sigma[9 * g + 3];
                auto v3rho2sigma_abc = v3rho2sigma[9 * g + 4];
                auto v3rho2sigma_abb = v3rho2sigma[9 * g + 5];
                auto v3rho2sigma_bba = v3rho2sigma[9 * g + 6];
                auto v3rho2sigma_bbc = v3rho2sigma[9 * g + 7];

                auto v3rhosigma2_aaa = v3rhosigma2[12 * g + 0];
                auto v3rhosigma2_aac = v3rhosigma2[12 * g + 1];
                auto v3rhosigma2_aab = v3rhosigma2[12 * g + 2];
                auto v3rhosigma2_acc = v3rhosigma2[12 * g + 3];
                auto v3rhosigma2_acb = v3rhosigma2[12 * g + 4];
                auto v3rhosigma2_abb = v3rhosigma2[12 * g + 5];
                auto v3rhosigma2_baa = v3rhosigma2[12 * g + 6];
                auto v3rhosigma2_bac = v3rhosigma2[12 * g + 7];
                auto v3rhosigma2_bab = v3rhosigma2[12 * g + 8];
                auto v3rhosigma2_bcc = v3rhosigma2[12 * g + 9];
                auto v3rhosigma2_bcb = v3rhosigma2[12 * g + 10];

                auto v3sigma3_aaa = v3sigma3[10 * g + 0];
                auto v3sigma3_aac = v3sigma3[10 * g + 1];
                auto v3sigma3_aab = v3sigma3[10 * g + 2];
                auto v3sigma3_acc = v3sigma3[10 * g + 3];
                auto v3sigma3_acb = v3sigma3[10 * g + 4];
                auto v3sigma3_abb = v3sigma3[10 * g + 5];
                auto v3sigma3_ccc = v3sigma3[10 * g + 6];
                auto v3sigma3_ccb = v3sigma3[10 * g + 7];
                auto v3sigma3_cbb = v3sigma3[10 * g + 8];

                // fourth-order terms
                auto v4rho4_aaaa = v4rho4[5 * g + 0];
                auto v4rho4_aaab = v4rho4[5 * g + 1];
                auto v4rho4_aabb = v4rho4[5 * g + 2];
                auto v4rho4_abbb = v4rho4[5 * g + 3];
                auto v4rho4_bbbb = v4rho4[5 * g + 4];

                auto v4rho3sigma_aaaa = v4rho3sigma[12 * g + 0];
                auto v4rho3sigma_aaac = v4rho3sigma[12 * g + 1];
                auto v4rho3sigma_aaab = v4rho3sigma[12 * g + 2];
                auto v4rho3sigma_aaba = v4rho3sigma[12 * g + 3];
                auto v4rho3sigma_aabc = v4rho3sigma[12 * g + 4];
                auto v4rho3sigma_aabb = v4rho3sigma[12 * g + 5];
                auto v4rho3sigma_abba = v4rho3sigma[12 * g + 6];
                auto v4rho3sigma_abbc = v4rho3sigma[12 * g + 7];
                auto v4rho3sigma_abbb = v4rho3sigma[12 * g + 8];
                auto v4rho3sigma_bbba = v4rho3sigma[12 * g + 9];
                auto v4rho3sigma_bbbc = v4rho3sigma[12 * g + 10];
                auto v4rho3sigma_bbbb = v4rho3sigma[12 * g + 11];

                auto v4rho3lapl_aaaa = v4rho3lapl[8 * g + 0];
                auto v4rho3lapl_aaab = v4rho3lapl[8 * g + 1];
                auto v4rho3lapl_aaba = v4rho3lapl[8 * g + 2];
                auto v4rho3lapl_aabb = v4rho3lapl[8 * g + 3];
                auto v4rho3lapl_abba = v4rho3lapl[8 * g + 4];
                auto v4rho3lapl_abbb = v4rho3lapl[8 * g + 5];
                auto v4rho3lapl_bbba = v4rho3lapl[8 * g + 6];
                auto v4rho3lapl_bbbb = v4rho3lapl[8 * g + 7];

                auto v4rho3tau_aaaa = v4rho3tau[8 * g + 0];
                auto v4rho3tau_aaab = v4rho3tau[8 * g + 1];
                auto v4rho3tau_aaba = v4rho3tau[8 * g + 2];
                auto v4rho3tau_aabb = v4rho3tau[8 * g + 3];
                auto v4rho3tau_abba = v4rho3tau[8 * g + 4];
                auto v4rho3tau_abbb = v4rho3tau[8 * g + 5];
                auto v4rho3tau_bbba = v4rho3tau[8 * g + 6];
                auto v4rho3tau_bbbb = v4rho3tau[8 * g + 7];

                auto v4rho2sigma2_aaaa = v4rho2sigma2[18 * g + 0];
                auto v4rho2sigma2_aaac = v4rho2sigma2[18 * g + 1];
                auto v4rho2sigma2_aaab = v4rho2sigma2[18 * g + 2];
                auto v4rho2sigma2_aacc = v4rho2sigma2[18 * g + 3];
                auto v4rho2sigma2_aacb = v4rho2sigma2[18 * g + 4];
                auto v4rho2sigma2_aabb = v4rho2sigma2[18 * g + 5];
                auto v4rho2sigma2_abaa = v4rho2sigma2[18 * g + 6];
                auto v4rho2sigma2_abac = v4rho2sigma2[18 * g + 7];
                auto v4rho2sigma2_abab = v4rho2sigma2[18 * g + 8];
                auto v4rho2sigma2_abcc = v4rho2sigma2[18 * g + 9];
                auto v4rho2sigma2_abcb = v4rho2sigma2[18 * g + 10];
                auto v4rho2sigma2_abbb = v4rho2sigma2[18 * g + 11];
                auto v4rho2sigma2_bbaa = v4rho2sigma2[18 * g + 12];
                auto v4rho2sigma2_bbac = v4rho2sigma2[18 * g + 13];
                auto v4rho2sigma2_bbab = v4rho2sigma2[18 * g + 14];
                auto v4rho2sigma2_bbcc = v4rho2sigma2[18 * g + 15];
                auto v4rho2sigma2_bbcb = v4rho2sigma2[18 * g + 16];
                auto v4rho2sigma2_bbbb = v4rho2sigma2[18 * g + 17];

                auto v4rho2sigmalapl_aaaa = v4rho2sigmalapl[18 * g + 0];
                auto v4rho2sigmalapl_aaab = v4rho2sigmalapl[18 * g + 1];
                auto v4rho2sigmalapl_aaca = v4rho2sigmalapl[18 * g + 2];
                auto v4rho2sigmalapl_aacb = v4rho2sigmalapl[18 * g + 3];
                auto v4rho2sigmalapl_aaba = v4rho2sigmalapl[18 * g + 4];
                auto v4rho2sigmalapl_aabb = v4rho2sigmalapl[18 * g + 5];
                auto v4rho2sigmalapl_abaa = v4rho2sigmalapl[18 * g + 6];
                auto v4rho2sigmalapl_abab = v4rho2sigmalapl[18 * g + 7];
                auto v4rho2sigmalapl_abca = v4rho2sigmalapl[18 * g + 8];
                auto v4rho2sigmalapl_abcb = v4rho2sigmalapl[18 * g + 9];
                auto v4rho2sigmalapl_abba = v4rho2sigmalapl[18 * g + 10];
                auto v4rho2sigmalapl_abbb = v4rho2sigmalapl[18 * g + 11];
                auto v4rho2sigmalapl_bbaa = v4rho2sigmalapl[18 * g + 12];
                auto v4rho2sigmalapl_bbab = v4rho2sigmalapl[18 * g + 13];
                auto v4rho2sigmalapl_bbca = v4rho2sigmalapl[18 * g + 14];
                auto v4rho2sigmalapl_bbcb = v4rho2sigmalapl[18 * g + 15];
                auto v4rho2sigmalapl_bbba = v4rho2sigmalapl[18 * g + 16];
                auto v4rho2sigmalapl_bbbb = v4rho2sigmalapl[18 * g + 17];

                auto v4rho2sigmatau_aaaa = v4rho2sigmatau[18 * g + 0];
                auto v4rho2sigmatau_aaab = v4rho2sigmatau[18 * g + 1];
                auto v4rho2sigmatau_aaca = v4rho2sigmatau[18 * g + 2];
                auto v4rho2sigmatau_aacb = v4rho2sigmatau[18 * g + 3];
                auto v4rho2sigmatau_aaba = v4rho2sigmatau[18 * g + 4];
                auto v4rho2sigmatau_aabb = v4rho2sigmatau[18 * g + 5];
                auto v4rho2sigmatau_abaa = v4rho2sigmatau[18 * g + 6];
                auto v4rho2sigmatau_abab = v4rho2sigmatau[18 * g + 7];
                auto v4rho2sigmatau_abca = v4rho2sigmatau[18 * g + 8];
                auto v4rho2sigmatau_abcb = v4rho2sigmatau[18 * g + 9];
                auto v4rho2sigmatau_abba = v4rho2sigmatau[18 * g + 10];
                auto v4rho2sigmatau_abbb = v4rho2sigmatau[18 * g + 11];
                auto v4rho2sigmatau_bbaa = v4rho2sigmatau[18 * g + 12];
                auto v4rho2sigmatau_bbab = v4rho2sigmatau[18 * g + 13];
                auto v4rho2sigmatau_bbca = v4rho2sigmatau[18 * g + 14];
                auto v4rho2sigmatau_bbcb = v4rho2sigmatau[18 * g + 15];
                auto v4rho2sigmatau_bbba = v4rho2sigmatau[18 * g + 16];
                auto v4rho2sigmatau_bbbb = v4rho2sigmatau[18 * g + 17];

                auto v4rho2lapl2_aaaa = v4rho2lapl2[9 * g + 0];
                auto v4rho2lapl2_aaab = v4rho2lapl2[9 * g + 1];
                auto v4rho2lapl2_aabb = v4rho2lapl2[9 * g + 2];
                auto v4rho2lapl2_abaa = v4rho2lapl2[9 * g + 3];
                auto v4rho2lapl2_abab = v4rho2lapl2[9 * g + 4];
                auto v4rho2lapl2_abbb = v4rho2lapl2[9 * g + 5];
                auto v4rho2lapl2_bbaa = v4rho2lapl2[9 * g + 6];
                auto v4rho2lapl2_bbab = v4rho2lapl2[9 * g + 7];
                auto v4rho2lapl2_bbbb = v4rho2lapl2[9 * g + 8];

                auto v4rho2lapltau_aaaa = v4rho2lapltau[12 * g + 0];
                auto v4rho2lapltau_aaab = v4rho2lapltau[12 * g + 1];
                auto v4rho2lapltau_aaba = v4rho2lapltau[12 * g + 2];
                auto v4rho2lapltau_aabb = v4rho2lapltau[12 * g + 3];
                auto v4rho2lapltau_abaa = v4rho2lapltau[12 * g + 4];
                auto v4rho2lapltau_abab = v4rho2lapltau[12 * g + 5];
                auto v4rho2lapltau_abba = v4rho2lapltau[12 * g + 6];
                auto v4rho2lapltau_abbb = v4rho2lapltau[12 * g + 7];
                auto v4rho2lapltau_bbaa = v4rho2lapltau[12 * g + 8];
                auto v4rho2lapltau_bbab = v4rho2lapltau[12 * g + 9];
                auto v4rho2lapltau_bbba = v4rho2lapltau[12 * g + 10];
                auto v4rho2lapltau_bbbb = v4rho2lapltau[12 * g + 11];

                auto v4rho2tau2_aaaa = v4rho2tau2[9 * g + 0];
                auto v4rho2tau2_aaab = v4rho2tau2[9 * g + 1];
                auto v4rho2tau2_aabb = v4rho2tau2[9 * g + 2];
                auto v4rho2tau2_abaa = v4rho2tau2[9 * g + 3];
                auto v4rho2tau2_abab = v4rho2tau2[9 * g + 4];
                auto v4rho2tau2_abbb = v4rho2tau2[9 * g + 5];
                auto v4rho2tau2_bbaa = v4rho2tau2[9 * g + 6];
                auto v4rho2tau2_bbab = v4rho2tau2[9 * g + 7];
                auto v4rho2tau2_bbbb = v4rho2tau2[9 * g + 8];

                auto v4rhosigma3_aaaa = v4rhosigma3[20 * g + 0];
                auto v4rhosigma3_aaac = v4rhosigma3[20 * g + 1];
                auto v4rhosigma3_aaab = v4rhosigma3[20 * g + 2];
                auto v4rhosigma3_aacc = v4rhosigma3[20 * g + 3];
                auto v4rhosigma3_aacb = v4rhosigma3[20 * g + 4];
                auto v4rhosigma3_aabb = v4rhosigma3[20 * g + 5];
                auto v4rhosigma3_accc = v4rhosigma3[20 * g + 6];
                auto v4rhosigma3_accb = v4rhosigma3[20 * g + 7];
                auto v4rhosigma3_acbb = v4rhosigma3[20 * g + 8];
                auto v4rhosigma3_abbb = v4rhosigma3[20 * g + 9];
                auto v4rhosigma3_baaa = v4rhosigma3[20 * g + 10];
                auto v4rhosigma3_baac = v4rhosigma3[20 * g + 11];
                auto v4rhosigma3_baab = v4rhosigma3[20 * g + 12];
                auto v4rhosigma3_bacc = v4rhosigma3[20 * g + 13];
                auto v4rhosigma3_bacb = v4rhosigma3[20 * g + 14];
                auto v4rhosigma3_babb = v4rhosigma3[20 * g + 15];
                auto v4rhosigma3_bccc = v4rhosigma3[20 * g + 16];
                auto v4rhosigma3_bccb = v4rhosigma3[20 * g + 17];
                auto v4rhosigma3_bcbb = v4rhosigma3[20 * g + 18];
                auto v4rhosigma3_bbbb = v4rhosigma3[20 * g + 19];

                auto v4rhosigma2lapl_aaaa = v4rhosigma2lapl[36 * g + 0];
                auto v4rhosigma2lapl_aaab = v4rhosigma2lapl[36 * g + 1];
                auto v4rhosigma2lapl_aaca = v4rhosigma2lapl[36 * g + 2];
                auto v4rhosigma2lapl_aacb = v4rhosigma2lapl[36 * g + 3];
                auto v4rhosigma2lapl_aaba = v4rhosigma2lapl[36 * g + 4];
                auto v4rhosigma2lapl_aabb = v4rhosigma2lapl[36 * g + 5];
                auto v4rhosigma2lapl_acca = v4rhosigma2lapl[36 * g + 6];
                auto v4rhosigma2lapl_accb = v4rhosigma2lapl[36 * g + 7];
                auto v4rhosigma2lapl_acba = v4rhosigma2lapl[36 * g + 8];
                auto v4rhosigma2lapl_acbb = v4rhosigma2lapl[36 * g + 9];
                auto v4rhosigma2lapl_abba = v4rhosigma2lapl[36 * g + 10];
                auto v4rhosigma2lapl_abbb = v4rhosigma2lapl[36 * g + 11];
                auto v4rhosigma2lapl_baaa = v4rhosigma2lapl[36 * g + 12];
                auto v4rhosigma2lapl_baab = v4rhosigma2lapl[36 * g + 13];
                auto v4rhosigma2lapl_baca = v4rhosigma2lapl[36 * g + 14];
                auto v4rhosigma2lapl_bacb = v4rhosigma2lapl[36 * g + 15];
                auto v4rhosigma2lapl_baba = v4rhosigma2lapl[36 * g + 16];
                auto v4rhosigma2lapl_babb = v4rhosigma2lapl[36 * g + 17];
                auto v4rhosigma2lapl_bcca = v4rhosigma2lapl[36 * g + 18];
                auto v4rhosigma2lapl_bccb = v4rhosigma2lapl[36 * g + 19];
                auto v4rhosigma2lapl_bcba = v4rhosigma2lapl[36 * g + 20];
                auto v4rhosigma2lapl_bcbb = v4rhosigma2lapl[36 * g + 21];
                auto v4rhosigma2lapl_bbba = v4rhosigma2lapl[36 * g + 22];
                auto v4rhosigma2lapl_bbbb = v4rhosigma2lapl[36 * g + 23];

                auto v4rhosigma2tau_aaaa = v4rhosigma2tau[36 * g + 0];
                auto v4rhosigma2tau_aaab = v4rhosigma2tau[36 * g + 1];
                auto v4rhosigma2tau_aaca = v4rhosigma2tau[36 * g + 2];
                auto v4rhosigma2tau_aacb = v4rhosigma2tau[36 * g + 3];
                auto v4rhosigma2tau_aaba = v4rhosigma2tau[36 * g + 4];
                auto v4rhosigma2tau_aabb = v4rhosigma2tau[36 * g + 5];
                auto v4rhosigma2tau_acca = v4rhosigma2tau[36 * g + 6];
                auto v4rhosigma2tau_accb = v4rhosigma2tau[36 * g + 7];
                auto v4rhosigma2tau_acba = v4rhosigma2tau[36 * g + 8];
                auto v4rhosigma2tau_acbb = v4rhosigma2tau[36 * g + 9];
                auto v4rhosigma2tau_abba = v4rhosigma2tau[36 * g + 10];
                auto v4rhosigma2tau_abbb = v4rhosigma2tau[36 * g + 11];
                auto v4rhosigma2tau_baaa = v4rhosigma2tau[36 * g + 12];
                auto v4rhosigma2tau_baab = v4rhosigma2tau[36 * g + 13];
                auto v4rhosigma2tau_baca = v4rhosigma2tau[36 * g + 14];
                auto v4rhosigma2tau_bacb = v4rhosigma2tau[36 * g + 15];
                auto v4rhosigma2tau_baba = v4rhosigma2tau[36 * g + 16];
                auto v4rhosigma2tau_babb = v4rhosigma2tau[36 * g + 17];
                auto v4rhosigma2tau_bcca = v4rhosigma2tau[36 * g + 18];
                auto v4rhosigma2tau_bccb = v4rhosigma2tau[36 * g + 19];
                auto v4rhosigma2tau_bcba = v4rhosigma2tau[36 * g + 20];
                auto v4rhosigma2tau_bcbb = v4rhosigma2tau[36 * g + 21];
                auto v4rhosigma2tau_bbba = v4rhosigma2tau[36 * g + 22];
                auto v4rhosigma2tau_bbbb = v4rhosigma2tau[36 * g + 23];

                auto v4rhosigmalapl2_aaaa = v4rhosigmalapl2[18 * g + 0];
                auto v4rhosigmalapl2_aaab = v4rhosigmalapl2[18 * g + 1];
                auto v4rhosigmalapl2_aabb = v4rhosigmalapl2[18 * g + 2];
                auto v4rhosigmalapl2_acaa = v4rhosigmalapl2[18 * g + 3];
                auto v4rhosigmalapl2_acab = v4rhosigmalapl2[18 * g + 4];
                auto v4rhosigmalapl2_acbb = v4rhosigmalapl2[18 * g + 5];
                auto v4rhosigmalapl2_abaa = v4rhosigmalapl2[18 * g + 6];
                auto v4rhosigmalapl2_abab = v4rhosigmalapl2[18 * g + 7];
                auto v4rhosigmalapl2_abbb = v4rhosigmalapl2[18 * g + 8];
                auto v4rhosigmalapl2_baaa = v4rhosigmalapl2[18 * g + 9];
                auto v4rhosigmalapl2_baab = v4rhosigmalapl2[18 * g + 10];
                auto v4rhosigmalapl2_babb = v4rhosigmalapl2[18 * g + 11];
                auto v4rhosigmalapl2_bcaa = v4rhosigmalapl2[18 * g + 12];
                auto v4rhosigmalapl2_bcab = v4rhosigmalapl2[18 * g + 13];
                auto v4rhosigmalapl2_bcbb = v4rhosigmalapl2[18 * g + 14];
                auto v4rhosigmalapl2_bbaa = v4rhosigmalapl2[18 * g + 15];
                auto v4rhosigmalapl2_bbab = v4rhosigmalapl2[18 * g + 16];
                auto v4rhosigmalapl2_bbbb = v4rhosigmalapl2[18 * g + 17];

                auto v4rhosigmalapltau_aaaa = v4rhosigmalapltau[24 * g + 0];
                auto v4rhosigmalapltau_aaab = v4rhosigmalapltau[24 * g + 1];
                auto v4rhosigmalapltau_aaba = v4rhosigmalapltau[24 * g + 2];
                auto v4rhosigmalapltau_aabb = v4rhosigmalapltau[24 * g + 3];
                auto v4rhosigmalapltau_acaa = v4rhosigmalapltau[24 * g + 4];
                auto v4rhosigmalapltau_acab = v4rhosigmalapltau[24 * g + 5];
                auto v4rhosigmalapltau_acba = v4rhosigmalapltau[24 * g + 6];
                auto v4rhosigmalapltau_acbb = v4rhosigmalapltau[24 * g + 7];
                auto v4rhosigmalapltau_abaa = v4rhosigmalapltau[24 * g + 8];
                auto v4rhosigmalapltau_abab = v4rhosigmalapltau[24 * g + 9];
                auto v4rhosigmalapltau_abba = v4rhosigmalapltau[24 * g + 10];
                auto v4rhosigmalapltau_abbb = v4rhosigmalapltau[24 * g + 11];
                auto v4rhosigmalapltau_baaa = v4rhosigmalapltau[24 * g + 12];
                auto v4rhosigmalapltau_baab = v4rhosigmalapltau[24 * g + 13];
                auto v4rhosigmalapltau_baba = v4rhosigmalapltau[24 * g + 14];
                auto v4rhosigmalapltau_babb = v4rhosigmalapltau[24 * g + 15];
                auto v4rhosigmalapltau_bcaa = v4rhosigmalapltau[24 * g + 16];
                auto v4rhosigmalapltau_bcab = v4rhosigmalapltau[24 * g + 17];
                auto v4rhosigmalapltau_bcba = v4rhosigmalapltau[24 * g + 18];
                auto v4rhosigmalapltau_bcbb = v4rhosigmalapltau[24 * g + 19];
                auto v4rhosigmalapltau_bbaa = v4rhosigmalapltau[24 * g + 20];
                auto v4rhosigmalapltau_bbab = v4rhosigmalapltau[24 * g + 21];
                auto v4rhosigmalapltau_bbba = v4rhosigmalapltau[24 * g + 22];
                auto v4rhosigmalapltau_bbbb = v4rhosigmalapltau[24 * g + 23];

                auto v4rhosigmatau2_aaaa = v4rhosigmatau2[36 * g + 0];
                auto v4rhosigmatau2_aaab = v4rhosigmatau2[36 * g + 1];
                auto v4rhosigmatau2_aabb = v4rhosigmatau2[36 * g + 2];
                auto v4rhosigmatau2_acaa = v4rhosigmatau2[36 * g + 3];
                auto v4rhosigmatau2_acab = v4rhosigmatau2[36 * g + 4];
                auto v4rhosigmatau2_acbb = v4rhosigmatau2[36 * g + 5];
                auto v4rhosigmatau2_abaa = v4rhosigmatau2[36 * g + 6];
                auto v4rhosigmatau2_abab = v4rhosigmatau2[36 * g + 7];
                auto v4rhosigmatau2_abbb = v4rhosigmatau2[36 * g + 8];
                auto v4rhosigmatau2_baaa = v4rhosigmatau2[36 * g + 9];
                auto v4rhosigmatau2_baab = v4rhosigmatau2[36 * g + 10];
                auto v4rhosigmatau2_babb = v4rhosigmatau2[36 * g + 11];
                auto v4rhosigmatau2_bcaa = v4rhosigmatau2[36 * g + 12];
                auto v4rhosigmatau2_bcab = v4rhosigmatau2[36 * g + 13];
                auto v4rhosigmatau2_bcbb = v4rhosigmatau2[36 * g + 14];
                auto v4rhosigmatau2_bbaa = v4rhosigmatau2[36 * g + 15];
                auto v4rhosigmatau2_bbab = v4rhosigmatau2[36 * g + 16];
                auto v4rhosigmatau2_bbbb = v4rhosigmatau2[36 * g + 17];

                auto v4rholapl3_aaaa = v4rholapl3[8 * g + 0];
                auto v4rholapl3_aaab = v4rholapl3[8 * g + 1];
                auto v4rholapl3_aabb = v4rholapl3[8 * g + 2];
                auto v4rholapl3_abbb = v4rholapl3[8 * g + 3];
                auto v4rholapl3_baaa = v4rholapl3[8 * g + 4];
                auto v4rholapl3_baab = v4rholapl3[8 * g + 5];
                auto v4rholapl3_babb = v4rholapl3[8 * g + 6];
                auto v4rholapl3_bbbb = v4rholapl3[8 * g + 7];

                auto v4rholapl2tau_aaaa = v4rholapl2tau[12 * g + 0];
                auto v4rholapl2tau_aaab = v4rholapl2tau[12 * g + 1];
                auto v4rholapl2tau_aaba = v4rholapl2tau[12 * g + 2];
                auto v4rholapl2tau_aabb = v4rholapl2tau[12 * g + 3];
                auto v4rholapl2tau_abba = v4rholapl2tau[12 * g + 4];
                auto v4rholapl2tau_abbb = v4rholapl2tau[12 * g + 5];
                auto v4rholapl2tau_baaa = v4rholapl2tau[12 * g + 6];
                auto v4rholapl2tau_baab = v4rholapl2tau[12 * g + 7];
                auto v4rholapl2tau_baba = v4rholapl2tau[12 * g + 8];
                auto v4rholapl2tau_babb = v4rholapl2tau[12 * g + 9];
                auto v4rholapl2tau_bbba = v4rholapl2tau[12 * g + 10];
                auto v4rholapl2tau_bbbb = v4rholapl2tau[12 * g + 11];

                auto v4rholapltau2_aaaa = v4rholapltau2[12 * g + 0];
                auto v4rholapltau2_aaab = v4rholapltau2[12 * g + 1];
                auto v4rholapltau2_aabb = v4rholapltau2[12 * g + 2];
                auto v4rholapltau2_abaa = v4rholapltau2[12 * g + 3];
                auto v4rholapltau2_abab = v4rholapltau2[12 * g + 4];
                auto v4rholapltau2_abbb = v4rholapltau2[12 * g + 5];
                auto v4rholapltau2_baaa = v4rholapltau2[12 * g + 6];
                auto v4rholapltau2_baab = v4rholapltau2[12 * g + 7];
                auto v4rholapltau2_babb = v4rholapltau2[12 * g + 8];
                auto v4rholapltau2_bbaa = v4rholapltau2[12 * g + 9];
                auto v4rholapltau2_bbab = v4rholapltau2[12 * g + 10];
                auto v4rholapltau2_bbbb = v4rholapltau2[12 * g + 11];

                auto v4rhotau3_aaaa = v4rhotau3[8 * g + 0];
                auto v4rhotau3_aaab = v4rhotau3[8 * g + 1];
                auto v4rhotau3_aabb = v4rhotau3[8 * g + 2];
                auto v4rhotau3_abbb = v4rhotau3[8 * g + 3];
                auto v4rhotau3_baaa = v4rhotau3[8 * g + 4];
                auto v4rhotau3_baab = v4rhotau3[8 * g + 5];
                auto v4rhotau3_babb = v4rhotau3[8 * g + 6];
                auto v4rhotau3_bbbb = v4rhotau3[8 * g + 7];

                auto v4sigma4_aaaa = v4sigma4[15 * g + 0];
                auto v4sigma4_aaac = v4sigma4[15 * g + 1];
                auto v4sigma4_aaab = v4sigma4[15 * g + 2];
                auto v4sigma4_aacc = v4sigma4[15 * g + 3];
                auto v4sigma4_aacb = v4sigma4[15 * g + 4];
                auto v4sigma4_aabb = v4sigma4[15 * g + 5];
                auto v4sigma4_accc = v4sigma4[15 * g + 6];
                auto v4sigma4_accb = v4sigma4[15 * g + 7];
                auto v4sigma4_acbb = v4sigma4[15 * g + 8];
                auto v4sigma4_abbb = v4sigma4[15 * g + 9];
                auto v4sigma4_cccc = v4sigma4[15 * g + 10];
                auto v4sigma4_cccb = v4sigma4[15 * g + 11];
                auto v4sigma4_ccbb = v4sigma4[15 * g + 12];
                auto v4sigma4_cbbb = v4sigma4[15 * g + 13];
                auto v4sigma4_bbbb = v4sigma4[15 * g + 14];

                auto v4sigma3lapl_aaaa = v4sigma3lapl[20 * g + 0];
                auto v4sigma3lapl_aaab = v4sigma3lapl[20 * g + 1];
                auto v4sigma3lapl_aaca = v4sigma3lapl[20 * g + 2];
                auto v4sigma3lapl_aacb = v4sigma3lapl[20 * g + 3];
                auto v4sigma3lapl_aaba = v4sigma3lapl[20 * g + 4];
                auto v4sigma3lapl_aabb = v4sigma3lapl[20 * g + 5];
                auto v4sigma3lapl_acca = v4sigma3lapl[20 * g + 6];
                auto v4sigma3lapl_accb = v4sigma3lapl[20 * g + 7];
                auto v4sigma3lapl_acba = v4sigma3lapl[20 * g + 8];
                auto v4sigma3lapl_acbb = v4sigma3lapl[20 * g + 9];
                auto v4sigma3lapl_abba = v4sigma3lapl[20 * g + 10];
                auto v4sigma3lapl_abbb = v4sigma3lapl[20 * g + 11];
                auto v4sigma3lapl_ccca = v4sigma3lapl[20 * g + 12];
                auto v4sigma3lapl_cccb = v4sigma3lapl[20 * g + 13];
                auto v4sigma3lapl_ccba = v4sigma3lapl[20 * g + 14];
                auto v4sigma3lapl_ccbb = v4sigma3lapl[20 * g + 15];
                auto v4sigma3lapl_cbba = v4sigma3lapl[20 * g + 16];
                auto v4sigma3lapl_cbbb = v4sigma3lapl[20 * g + 17];
                auto v4sigma3lapl_bbba = v4sigma3lapl[20 * g + 18];
                auto v4sigma3lapl_bbbb = v4sigma3lapl[20 * g + 19];

                auto v4sigma3tau_aaaa = v4sigma3tau[30 * g + 0];
                auto v4sigma3tau_aaab = v4sigma3tau[30 * g + 1];
                auto v4sigma3tau_aaca = v4sigma3tau[30 * g + 2];
                auto v4sigma3tau_aacb = v4sigma3tau[30 * g + 3];
                auto v4sigma3tau_aaba = v4sigma3tau[30 * g + 4];
                auto v4sigma3tau_aabb = v4sigma3tau[30 * g + 5];
                auto v4sigma3tau_acca = v4sigma3tau[30 * g + 6];
                auto v4sigma3tau_accb = v4sigma3tau[30 * g + 7];
                auto v4sigma3tau_acba = v4sigma3tau[30 * g + 8];
                auto v4sigma3tau_acbb = v4sigma3tau[30 * g + 9];
                auto v4sigma3tau_abba = v4sigma3tau[30 * g + 10];
                auto v4sigma3tau_abbb = v4sigma3tau[30 * g + 11];
                auto v4sigma3tau_ccca = v4sigma3tau[30 * g + 12];
                auto v4sigma3tau_cccb = v4sigma3tau[30 * g + 13];
                auto v4sigma3tau_ccba = v4sigma3tau[30 * g + 14];
                auto v4sigma3tau_ccbb = v4sigma3tau[30 * g + 15];
                auto v4sigma3tau_cbba = v4sigma3tau[30 * g + 16];
                auto v4sigma3tau_cbbb = v4sigma3tau[30 * g + 17];
                auto v4sigma3tau_bbba = v4sigma3tau[30 * g + 18];
                auto v4sigma3tau_bbbb = v4sigma3tau[30 * g + 19];

                auto v4sigma2lapl2_aaaa = v4sigma2lapl2[18 * g + 0];
                auto v4sigma2lapl2_aaab = v4sigma2lapl2[18 * g + 1];
                auto v4sigma2lapl2_aabb = v4sigma2lapl2[18 * g + 2];
                auto v4sigma2lapl2_acaa = v4sigma2lapl2[18 * g + 3];
                auto v4sigma2lapl2_acab = v4sigma2lapl2[18 * g + 4];
                auto v4sigma2lapl2_acbb = v4sigma2lapl2[18 * g + 5];
                auto v4sigma2lapl2_abaa = v4sigma2lapl2[18 * g + 6];
                auto v4sigma2lapl2_abab = v4sigma2lapl2[18 * g + 7];
                auto v4sigma2lapl2_abbb = v4sigma2lapl2[18 * g + 8];
                auto v4sigma2lapl2_ccaa = v4sigma2lapl2[18 * g + 9];
                auto v4sigma2lapl2_ccab = v4sigma2lapl2[18 * g + 10];
                auto v4sigma2lapl2_ccbb = v4sigma2lapl2[18 * g + 11];
                auto v4sigma2lapl2_cbaa = v4sigma2lapl2[18 * g + 12];
                auto v4sigma2lapl2_cbab = v4sigma2lapl2[18 * g + 13];
                auto v4sigma2lapl2_cbbb = v4sigma2lapl2[18 * g + 14];
                auto v4sigma2lapl2_bbaa = v4sigma2lapl2[18 * g + 15];
                auto v4sigma2lapl2_bbab = v4sigma2lapl2[18 * g + 16];
                auto v4sigma2lapl2_bbbb = v4sigma2lapl2[18 * g + 17];

                auto v4sigma2lapltau_aaaa = v4sigma2lapltau[24 * g + 0];
                auto v4sigma2lapltau_aaab = v4sigma2lapltau[24 * g + 1];
                auto v4sigma2lapltau_aaba = v4sigma2lapltau[24 * g + 2];
                auto v4sigma2lapltau_aabb = v4sigma2lapltau[24 * g + 3];
                auto v4sigma2lapltau_acaa = v4sigma2lapltau[24 * g + 4];
                auto v4sigma2lapltau_acab = v4sigma2lapltau[24 * g + 5];
                auto v4sigma2lapltau_acba = v4sigma2lapltau[24 * g + 6];
                auto v4sigma2lapltau_acbb = v4sigma2lapltau[24 * g + 7];
                auto v4sigma2lapltau_abaa = v4sigma2lapltau[24 * g + 8];
                auto v4sigma2lapltau_abab = v4sigma2lapltau[24 * g + 9];
                auto v4sigma2lapltau_abba = v4sigma2lapltau[24 * g + 10];
                auto v4sigma2lapltau_abbb = v4sigma2lapltau[24 * g + 11];
                auto v4sigma2lapltau_ccaa = v4sigma2lapltau[24 * g + 12];
                auto v4sigma2lapltau_ccab = v4sigma2lapltau[24 * g + 13];
                auto v4sigma2lapltau_ccba = v4sigma2lapltau[24 * g + 14];
                auto v4sigma2lapltau_ccbb = v4sigma2lapltau[24 * g + 15];
                auto v4sigma2lapltau_cbaa = v4sigma2lapltau[24 * g + 16];
                auto v4sigma2lapltau_cbab = v4sigma2lapltau[24 * g + 17];
                auto v4sigma2lapltau_cbba = v4sigma2lapltau[24 * g + 18];
                auto v4sigma2lapltau_cbbb = v4sigma2lapltau[24 * g + 19];
                auto v4sigma2lapltau_bbaa = v4sigma2lapltau[24 * g + 20];
                auto v4sigma2lapltau_bbab = v4sigma2lapltau[24 * g + 21];
                auto v4sigma2lapltau_bbba = v4sigma2lapltau[24 * g + 22];
                auto v4sigma2lapltau_bbbb = v4sigma2lapltau[24 * g + 23];

                auto v4sigma2tau2_aaaa = v4sigma2tau2[18 * g + 0];
                auto v4sigma2tau2_aaab = v4sigma2tau2[18 * g + 1];
                auto v4sigma2tau2_aabb = v4sigma2tau2[18 * g + 2];
                auto v4sigma2tau2_acaa = v4sigma2tau2[18 * g + 3];
                auto v4sigma2tau2_acab = v4sigma2tau2[18 * g + 4];
                auto v4sigma2tau2_acbb = v4sigma2tau2[18 * g + 5];
                auto v4sigma2tau2_abaa = v4sigma2tau2[18 * g + 6];
                auto v4sigma2tau2_abab = v4sigma2tau2[18 * g + 7];
                auto v4sigma2tau2_abbb = v4sigma2tau2[18 * g + 8];
                auto v4sigma2tau2_ccaa = v4sigma2tau2[18 * g + 9];
                auto v4sigma2tau2_ccab = v4sigma2tau2[18 * g + 10];
                auto v4sigma2tau2_ccbb = v4sigma2tau2[18 * g + 11];
                auto v4sigma2tau2_cbaa = v4sigma2tau2[18 * g + 12];
                auto v4sigma2tau2_cbab = v4sigma2tau2[18 * g + 13];
                auto v4sigma2tau2_cbbb = v4sigma2tau2[18 * g + 14];
                auto v4sigma2tau2_bbaa = v4sigma2tau2[18 * g + 15];
                auto v4sigma2tau2_bbab = v4sigma2tau2[18 * g + 16];
                auto v4sigma2tau2_bbbb = v4sigma2tau2[18 * g + 17];

                auto v4sigmalapl3_aaaa = v4sigmalapl3[12 * g + 0];
                auto v4sigmalapl3_aaab = v4sigmalapl3[12 * g + 1];
                auto v4sigmalapl3_aabb = v4sigmalapl3[12 * g + 2];
                auto v4sigmalapl3_abbb = v4sigmalapl3[12 * g + 3];
                auto v4sigmalapl3_caaa = v4sigmalapl3[12 * g + 4];
                auto v4sigmalapl3_caab = v4sigmalapl3[12 * g + 5];
                auto v4sigmalapl3_cabb = v4sigmalapl3[12 * g + 6];
                auto v4sigmalapl3_cbbb = v4sigmalapl3[12 * g + 7];
                auto v4sigmalapl3_baaa = v4sigmalapl3[12 * g + 8];
                auto v4sigmalapl3_baab = v4sigmalapl3[12 * g + 9];
                auto v4sigmalapl3_babb = v4sigmalapl3[12 * g + 10];
                auto v4sigmalapl3_bbbb = v4sigmalapl3[12 * g + 11];

                auto v4sigmalapl2tau_aaaa = v4sigmalapl2tau[18 * g + 0];
                auto v4sigmalapl2tau_aaab = v4sigmalapl2tau[18 * g + 1];
                auto v4sigmalapl2tau_aaba = v4sigmalapl2tau[18 * g + 2];
                auto v4sigmalapl2tau_aabb = v4sigmalapl2tau[18 * g + 3];
                auto v4sigmalapl2tau_abba = v4sigmalapl2tau[18 * g + 4];
                auto v4sigmalapl2tau_abbb = v4sigmalapl2tau[18 * g + 5];
                auto v4sigmalapl2tau_caaa = v4sigmalapl2tau[18 * g + 6];
                auto v4sigmalapl2tau_caab = v4sigmalapl2tau[18 * g + 7];
                auto v4sigmalapl2tau_caba = v4sigmalapl2tau[18 * g + 8];
                auto v4sigmalapl2tau_cabb = v4sigmalapl2tau[18 * g + 9];
                auto v4sigmalapl2tau_cbba = v4sigmalapl2tau[18 * g + 10];
                auto v4sigmalapl2tau_cbbb = v4sigmalapl2tau[18 * g + 11];
                auto v4sigmalapl2tau_baaa = v4sigmalapl2tau[18 * g + 12];
                auto v4sigmalapl2tau_baab = v4sigmalapl2tau[18 * g + 13];
                auto v4sigmalapl2tau_baba = v4sigmalapl2tau[18 * g + 14];
                auto v4sigmalapl2tau_babb = v4sigmalapl2tau[18 * g + 15];
                auto v4sigmalapl2tau_bbba = v4sigmalapl2tau[18 * g + 16];
                auto v4sigmalapl2tau_bbbb = v4sigmalapl2tau[18 * g + 17];

                auto v4sigmalapltau2_aaaa = v4sigmalapltau2[18 * g + 0];
                auto v4sigmalapltau2_aaab = v4sigmalapltau2[18 * g + 1];
                auto v4sigmalapltau2_aabb = v4sigmalapltau2[18 * g + 2];
                auto v4sigmalapltau2_abaa = v4sigmalapltau2[18 * g + 3];
                auto v4sigmalapltau2_abab = v4sigmalapltau2[18 * g + 4];
                auto v4sigmalapltau2_abbb = v4sigmalapltau2[18 * g + 5];
                auto v4sigmalapltau2_caaa = v4sigmalapltau2[18 * g + 6];
                auto v4sigmalapltau2_caab = v4sigmalapltau2[18 * g + 7];
                auto v4sigmalapltau2_cabb = v4sigmalapltau2[18 * g + 8];
                auto v4sigmalapltau2_cbaa = v4sigmalapltau2[18 * g + 9];
                auto v4sigmalapltau2_cbab = v4sigmalapltau2[18 * g + 10];
                auto v4sigmalapltau2_cbbb = v4sigmalapltau2[18 * g + 11];
                auto v4sigmalapltau2_baaa = v4sigmalapltau2[18 * g + 12];
                auto v4sigmalapltau2_baab = v4sigmalapltau2[18 * g + 13];
                auto v4sigmalapltau2_babb = v4sigmalapltau2[18 * g + 14];
                auto v4sigmalapltau2_bbaa = v4sigmalapltau2[18 * g + 15];
                auto v4sigmalapltau2_bbab = v4sigmalapltau2[18 * g + 16];
                auto v4sigmalapltau2_bbbb = v4sigmalapltau2[18 * g + 17];

                auto v4sigmatau3_aaaa = v4sigmatau3[12 * g + 0];
                auto v4sigmatau3_aaab = v4sigmatau3[12 * g + 1];
                auto v4sigmatau3_aabb = v4sigmatau3[12 * g + 2];
                auto v4sigmatau3_abbb = v4sigmatau3[12 * g + 3];
                auto v4sigmatau3_caaa = v4sigmatau3[12 * g + 4];
                auto v4sigmatau3_caab = v4sigmatau3[12 * g + 5];
                auto v4sigmatau3_cabb = v4sigmatau3[12 * g + 6];
                auto v4sigmatau3_cbbb = v4sigmatau3[12 * g + 7];
                auto v4sigmatau3_baaa = v4sigmatau3[12 * g + 8];
                auto v4sigmatau3_baab = v4sigmatau3[12 * g + 9];
                auto v4sigmatau3_babb = v4sigmatau3[12 * g + 10];
                auto v4sigmatau3_bbbb = v4sigmatau3[12 * g + 11];

                auto v4lapl4_aaaa = v4lapl4[5 * g + 0];
                auto v4lapl4_aaab = v4lapl4[5 * g + 1];
                auto v4lapl4_aabb = v4lapl4[5 * g + 2];
                auto v4lapl4_abbb = v4lapl4[5 * g + 3];
                auto v4lapl4_bbbb = v4lapl4[5 * g + 4];

                auto v4lapl3tau_aaaa = v4lapl3tau[8 * g + 0];
                auto v4lapl3tau_aaab = v4lapl3tau[8 * g + 1];
                auto v4lapl3tau_aaba = v4lapl3tau[8 * g + 2];
                auto v4lapl3tau_aabb = v4lapl3tau[8 * g + 3];
                auto v4lapl3tau_abba = v4lapl3tau[8 * g + 4];
                auto v4lapl3tau_abbb = v4lapl3tau[8 * g + 5];
                auto v4lapl3tau_bbba = v4lapl3tau[8 * g + 6];
                auto v4lapl3tau_bbbb = v4lapl3tau[8 * g + 7];

                auto v4lapl2tau2_aaaa = v4lapl2tau2[9 * g + 0];
                auto v4lapl2tau2_aaab = v4lapl2tau2[9 * g + 1];
                auto v4lapl2tau2_aabb = v4lapl2tau2[9 * g + 2];
                auto v4lapl2tau2_abaa = v4lapl2tau2[9 * g + 3];
                auto v4lapl2tau2_abab = v4lapl2tau2[9 * g + 4];
                auto v4lapl2tau2_abbb = v4lapl2tau2[9 * g + 5];
                auto v4lapl2tau2_bbaa = v4lapl2tau2[9 * g + 6];
                auto v4lapl2tau2_bbab = v4lapl2tau2[9 * g + 7];
                auto v4lapl2tau2_bbbb = v4lapl2tau2[9 * g + 8];

                auto v4lapltau3_aaaa = v4lapltau3[8 * g + 0];
                auto v4lapltau3_aaab = v4lapltau3[8 * g + 1];
                auto v4lapltau3_aabb = v4lapltau3[8 * g + 2];
                auto v4lapltau3_abbb = v4lapltau3[8 * g + 3];
                auto v4lapltau3_baaa = v4lapltau3[8 * g + 4];
                auto v4lapltau3_baab = v4lapltau3[8 * g + 5];
                auto v4lapltau3_babb = v4lapltau3[8 * g + 6];
                auto v4lapltau3_bbbb = v4lapltau3[8 * g + 7];

                auto v4tau4_aaaa = v4tau4[5 * g + 0];
                auto v4tau4_aaab = v4tau4[5 * g + 1];
                auto v4tau4_aabb = v4tau4[5 * g + 2];
                auto v4tau4_abbb = v4tau4[5 * g + 3];
                auto v4tau4_bbbb = v4tau4[5 * g + 4];

                // sums of functional derivatives

                // first-order
                double x = vsigma_c + 2.0 * vsigma_a;

                // second-order
                // rho
                double rr = v2rho2_aa + v2rho2_ab;
                double rx = 2.0*v2rhosigma_ac + 2.0*v2rhosigma_ab + 2.0*v2rhosigma_aa;
                double rt = v2rhotau_aa + v2rhotau_ab;
                double rl = v2rholapl_aa + v2rholapl_ab;

                // sigma and gamma
                double xr = v2rhosigma_bc + 2.0*v2rhosigma_ba + v2rhosigma_ac + 2.0 * v2rhosigma_aa;
                double xt = v2sigmatau_cb + 2.0*v2sigmatau_ab + v2sigmatau_ca + 2.0 * v2sigmatau_aa;
                double xl = v2sigmalapl_cb + 2.0*v2sigmalapl_ab + v2sigmalapl_ca + 2.0 * v2sigmalapl_aa;
                double xx = 2.0*v2sigma2_cc + 2.0*v2sigma2_cb + 6.0*v2sigma2_ac + 4.0*v2sigma2_ab + 4.0 * v2sigma2_aa;

                // tau
                double tt = v2tau2_aa + v2tau2_ab;
                double tx = 2.0 * v2sigmatau_ca + 2.0 * v2sigmatau_ba + 2.0 * v2sigmatau_aa;
                double tr = v2rhotau_aa + v2rhotau_ba;
                double tl = v2lapltau_aa + v2lapltau_ba;

                // lapl
                //double ll = v2lapl2_aa + v2lapl2_ab;
                //double lx = 2.0 * v2sigmalapl_ca + 2.0 * v2sigmalapl_ba + 2.0 * v2sigmalapl_aa;
                //double lr = v2rholapl_aa + v2rholapl_ba;
                //double lt = v2lapltau_aa + v2lapltau_ab;

                // Third-oder

                // sigma and gamma
                double xxx  = 4.0 * v3sigma3_ccc + 8.0 * v3sigma3_ccb + 4.0 * v3sigma3_cbb + 16.0 * v3sigma3_acc
                            + 24.0 * v3sigma3_acb + 8.0 * v3sigma3_abb + 20.0 * v3sigma3_aac + 16.0 * v3sigma3_aab + 8.0 * v3sigma3_aaa;
                double xxr  = 2.0 * v3rhosigma2_bcc + 2.0 * v3rhosigma2_bcb + 6.0 * v3rhosigma2_bac
                            + 4.0 * v3rhosigma2_bab + 4.0 * v3rhosigma2_baa + 2.0 * v3rhosigma2_acc
                            + 2.0 * v3rhosigma2_acb + 6.0 * v3rhosigma2_aac + 4.0 * v3rhosigma2_aab + 4.0 * v3rhosigma2_aaa;
                double xrt  = v3rhosigmatau_bcb + v3rhosigmatau_bca
                            + 2.0 * v3rhosigmatau_bab + 2.0 * v3rhosigmatau_baa
                            + 2.0 * v3rhosigmatau_aab + 2.0 * v3rhosigmatau_aaa
                            + v3rhosigmatau_acb + v3rhosigmatau_aca;
                double xxl  = 2.0 * v3sigma2lapl_ccb + 2.0 * v3sigma2lapl_cca + 2.0 * v3sigma2lapl_cbb
                            + 2.0 * v3sigma2lapl_cba + 6.0 * v3sigma2lapl_acb + 6.0 * v3sigma2lapl_aca
                            + 4.0 * v3sigma2lapl_abb + 4.0 * v3sigma2lapl_aba + 4.0 * v3sigma2lapl_aab + 4.0 * v3sigma2lapl_aaa;
                double xrr  = v3rho2sigma_bbc + 2.0 * v3rho2sigma_bba
                            + 2.0 * v3rho2sigma_abc + 4.0 * v3rho2sigma_aba
                            + v3rho2sigma_aac + 2.0 * v3rho2sigma_aaa;
                double xrl  = v3rhosigmalapl_bcb + v3rhosigmalapl_bca + 2.0 * v3rhosigmalapl_bab
                            + 2.0 * v3rhosigmalapl_baa + v3rhosigmalapl_acb + v3rhosigmalapl_aca
                            + 2.0 * v3rhosigmalapl_aab + 2.0 * v3rhosigmalapl_aaa;
                double xtl  = v3sigmalapltau_cbb + v3sigmalapltau_cba + v3sigmalapltau_cab
                            + v3sigmalapltau_caa + 2.0 * v3sigmalapltau_abb + 2.0 * v3sigmalapltau_aba
                            + 2.0 * v3sigmalapltau_aab + 2.0 * v3sigmalapltau_aaa;
                double xll  = v3sigmalapl2_cbb + 2.0 * v3sigmalapl2_cab + v3sigmalapl2_caa
                            + 2.0 * v3sigmalapl2_abb + 4.0 * v3sigmalapl2_aab + 2.0 * v3sigmalapl2_aaa;
                double xxt  = 2.0 * v3sigma2tau_ccb + 2.0 * v3sigma2tau_cca + 2.0 * v3sigma2tau_cbb
                            + 2.0 * v3sigma2tau_cba + 6.0 * v3sigma2tau_acb + 6.0 * v3sigma2tau_aca
                            + 4.0 * v3sigma2tau_abb + 4.0 * v3sigma2tau_aba + 4.0 * v3sigma2tau_aab + 4.0 * v3sigma2tau_aaa;
                double xtt  = v3sigmatau2_cbb + 2.0 * v3sigmatau2_cab + v3sigmatau2_caa
                            + 2.0 * v3sigmatau2_abb + 4.0 * v3sigmatau2_aab + 2.0 * v3sigmatau2_aaa;

                // rho
                double rrr  = v3rho3_abb + 2.0 * v3rho3_aab + v3rho3_aaa;
                double rrt  = v3rho2tau_abb + v3rho2tau_aba + v3rho2tau_aab + v3rho2tau_aaa;
                double rtx  = 2.0 * v3rhosigmatau_acb + 2.0 * v3rhosigmatau_aca
                            + 2.0 * v3rhosigmatau_abb + 2.0 * v3rhosigmatau_aba
                            + 2.0 * v3rhosigmatau_aab + 2.0 * v3rhosigmatau_aaa;
                double rrl  = v3rho2lapl_abb + v3rho2lapl_aba + v3rho2lapl_aab + v3rho2lapl_aaa;
                double rrx  = 2.0 * v3rho2sigma_abc + 2.0 * v3rho2sigma_abb
                            + 2.0 * v3rho2sigma_aba + 2.0 * v3rho2sigma_aac
                            + 2.0 * v3rho2sigma_aab + 2.0 * v3rho2sigma_aaa;
                double rtl  = v3rholapltau_abb + v3rholapltau_aba + v3rholapltau_aab + v3rholapltau_aaa;
                double rtt  = v3rhotau2_abb + 2.0 * v3rhotau2_aab + v3rhotau2_aaa;
                double rll  = v3rholapl2_abb + 2.0 * v3rholapl2_aab + v3rholapl2_aaa;
                double rlx  = 2.0 * v3rhosigmalapl_acb + 2.0 * v3rhosigmalapl_aca + 2.0 * v3rhosigmalapl_abb
                            + 2.0 * v3rhosigmalapl_aba + 2.0 * v3rhosigmalapl_aab + 2.0 * v3rhosigmalapl_aaa;
                double rxx  = 4.0 * v3rhosigma2_acc + 8.0 * v3rhosigma2_acb + 4.0 * v3rhosigma2_abb
                            + 8.0 * v3rhosigma2_aac + 8.0 * v3rhosigma2_aab + 4.0 * v3rhosigma2_aaa;

                // laplacian
                //double lll  = v3lapl3_abb + 2.0 * v3lapl3_aab + v3lapl3_aaa;
                //double llr  = v3rholapl2_bab + v3rholapl2_baa + v3rholapl2_aab + v3rholapl2_aaa;
                //double llt  = v3lapl2tau_abb + v3lapl2tau_aba + v3lapl2tau_aab + v3lapl2tau_aaa;
                //double llx  = 2.0 * v3sigmalapl2_cab + 2.0 * v3sigmalapl2_caa + 2.0 * v3sigmalapl2_bab
                //            + 2.0 * v3sigmalapl2_baa + 2.0 * v3sigmalapl2_aab + 2.0 * v3sigmalapl2_aaa;
                //double lrr  = v3rho2lapl_bba + 2.0 * v3rho2lapl_aba + v3rho2lapl_aaa;
                //double lrt  = v3rholapltau_bab + v3rholapltau_baa + v3rholapltau_aab + v3rholapltau_aaa;
                //double lrx  = 2.0 * v3rhosigmalapl_bca + 2.0 * v3rhosigmalapl_bba + 2.0 * v3rhosigmalapl_baa
                //            + 2.0 * v3rhosigmalapl_aca + 2.0 * v3rhosigmalapl_aba + 2.0 * v3rhosigmalapl_aaa;
                //double ltt  = v3lapltau2_abb + 2.0 * v3lapltau2_aab + v3lapltau2_aaa;
                //double ltx  = 2.0 * v3sigmalapltau_cab + 2.0 * v3sigmalapltau_caa + 2.0 * v3sigmalapltau_bab
                //            + 2.0 * v3sigmalapltau_baa + 2.0 * v3sigmalapltau_aab + 2.0 * v3sigmalapltau_aaa;
                //double lxx  = 4.0 * v3sigma2lapl_cca + 8.0 * v3sigma2lapl_cba + 4.0 * v3sigma2lapl_bba
                //            + 8.0 * v3sigma2lapl_aca + 8.0 * v3sigma2lapl_aba + 4.0 * v3sigma2lapl_aaa;

                // tau
                double trr  = v3rho2tau_bba + 2.0 * v3rho2tau_aba + v3rho2tau_aaa;
                double ttl  = v3lapltau2_bab + v3lapltau2_baa + v3lapltau2_aab + v3lapltau2_aaa;
                double trl  = v3rholapltau_bba + v3rholapltau_baa + v3rholapltau_aba + v3rholapltau_aaa;
                double tll  = v3lapl2tau_bba + 2.0 * v3lapl2tau_aba + v3lapl2tau_aaa;
                double tlx  = 2.0 * v3sigmalapltau_cba + 2.0 * v3sigmalapltau_caa + 2.0 * v3sigmalapltau_bba
                            + 2.0 * v3sigmalapltau_baa + 2.0 * v3sigmalapltau_aba + 2.0 * v3sigmalapltau_aaa;
                double ttt  = v3tau3_abb + 2.0 * v3tau3_aab + v3tau3_aaa;
                double ttr  = v3rhotau2_bab + v3rhotau2_baa + v3rhotau2_aab + v3rhotau2_aaa;
                double trx  = 2.0 * v3rhosigmatau_bca + 2.0 * v3rhosigmatau_bba
                            + 2.0 * v3rhosigmatau_baa + 2.0 * v3rhosigmatau_aca
                            + 2.0 * v3rhosigmatau_aba + 2.0 * v3rhosigmatau_aaa;
                double txx  = 4.0 * v3sigma2tau_cca + 8.0 * v3sigma2tau_cba + 4.0 * v3sigma2tau_bba
                            + 8.0 * v3sigma2tau_aca + 8.0 * v3sigma2tau_aba + 4.0 * v3sigma2tau_aaa;
                double ttx  = 2.0 * v3sigmatau2_cab + 2.0 * v3sigmatau2_caa + 2.0 * v3sigmatau2_bab
                            + 2.0 * v3sigmatau2_baa + 2.0 * v3sigmatau2_aab + 2.0 * v3sigmatau2_aaa;

                // fourth-order

                double rrrr  = v4rho4_abbb + 3.0 * v4rho4_aabb + 3.0 * v4rho4_aaab + v4rho4_aaaa;
                double rrrx  = 2.0 * v4rho3sigma_abbc + 2.0 * v4rho3sigma_abbb + 2.0 * v4rho3sigma_abba + 4.0 * v4rho3sigma_aabc
                             + 4.0 * v4rho3sigma_aabb + 4.0 * v4rho3sigma_aaba + 2.0 * v4rho3sigma_aaac + 2.0 * v4rho3sigma_aaab + 2.0 * v4rho3sigma_aaaa;
                double rrrl  = v4rho3lapl_abbb + v4rho3lapl_abba + 2.0 * v4rho3lapl_aabb + 2.0 * v4rho3lapl_aaba + v4rho3lapl_aaab + v4rho3lapl_aaaa;
                double rrrt  = v4rho3tau_abbb + v4rho3tau_abba + 2.0 * v4rho3tau_aabb + 2.0 * v4rho3tau_aaba + v4rho3tau_aaab + v4rho3tau_aaaa;
                double rrxx  = 4.0 * v4rho2sigma2_abcc + 8.0 * v4rho2sigma2_abcb + 4.0 * v4rho2sigma2_abbb + 8.0 * v4rho2sigma2_abac
                             + 8.0 * v4rho2sigma2_abab + 4.0 * v4rho2sigma2_abaa + 4.0 * v4rho2sigma2_aacc + 8.0 * v4rho2sigma2_aacb
                             + 4.0 * v4rho2sigma2_aabb + 8.0 * v4rho2sigma2_aaac + 8.0 * v4rho2sigma2_aaab + 4.0 * v4rho2sigma2_aaaa;
                double rrxl  = 2.0 * v4rho2sigmalapl_abcb + 2.0 * v4rho2sigmalapl_abca + 2.0 * v4rho2sigmalapl_abbb + 2.0 * v4rho2sigmalapl_abba + 2.0 * v4rho2sigmalapl_abab + 2.0 * v4rho2sigmalapl_abaa + 2.0 * v4rho2sigmalapl_aacb + 2.0 * v4rho2sigmalapl_aaca + 2.0 * v4rho2sigmalapl_aabb + 2.0 * v4rho2sigmalapl_aaba + 2.0 * v4rho2sigmalapl_aaab + 2.0 * v4rho2sigmalapl_aaaa;
                double rrxt  = 2.0 * v4rho2sigmatau_abcb + 2.0 * v4rho2sigmatau_abca + 2.0 * v4rho2sigmatau_abbb + 2.0 * v4rho2sigmatau_abba + 2.0 * v4rho2sigmatau_abab + 2.0 * v4rho2sigmatau_abaa + 2.0 * v4rho2sigmatau_aacb + 2.0 * v4rho2sigmatau_aaca + 2.0 * v4rho2sigmatau_aabb + 2.0 * v4rho2sigmatau_aaba + 2.0 * v4rho2sigmatau_aaab + 2.0 * v4rho2sigmatau_aaaa;
                double rrll  = v4rho2lapl2_abbb + 2.0 * v4rho2lapl2_abab + v4rho2lapl2_abaa + v4rho2lapl2_aabb + 2.0 * v4rho2lapl2_aaab + v4rho2lapl2_aaaa;
                double rrlt  = v4rho2lapltau_abbb + v4rho2lapltau_abba + v4rho2lapltau_abab + v4rho2lapltau_abaa + v4rho2lapltau_aabb + v4rho2lapltau_aaba + v4rho2lapltau_aaab + v4rho2lapltau_aaaa;
                double rrtt  = v4rho2tau2_abbb + 2.0 * v4rho2tau2_abab + v4rho2tau2_abaa + v4rho2tau2_aabb + 2.0 * v4rho2tau2_aaab + v4rho2tau2_aaaa;
                double rxxx  = 8.0 * v4rhosigma3_accc + 24.0 * v4rhosigma3_accb + 24.0 * v4rhosigma3_acbb + 8.0 * v4rhosigma3_abbb
                             + 24.0 * v4rhosigma3_aacc + 48.0 * v4rhosigma3_aacb + 24.0 * v4rhosigma3_aabb + 24.0 * v4rhosigma3_aaac
                             + 24.0 * v4rhosigma3_aaab + 8.0 * v4rhosigma3_aaaa;
                double rxxl  = 4.0 * v4rhosigma2lapl_accb + 4.0 * v4rhosigma2lapl_acca + 8.0 * v4rhosigma2lapl_acbb
                             + 8.0 * v4rhosigma2lapl_acba + 4.0 * v4rhosigma2lapl_abbb + 4.0 * v4rhosigma2lapl_abba
                             + 8.0 * v4rhosigma2lapl_aacb + 8.0 * v4rhosigma2lapl_aaca + 8.0 * v4rhosigma2lapl_aabb
                             + 8.0 * v4rhosigma2lapl_aaba + 4.0 * v4rhosigma2lapl_aaab + 4.0 * v4rhosigma2lapl_aaaa;
                double rxxt  = 4.0 * v4rhosigma2tau_accb + 4.0 * v4rhosigma2tau_acca + 8.0 * v4rhosigma2tau_acbb
                             + 8.0 * v4rhosigma2tau_acba + 4.0 * v4rhosigma2tau_abbb + 4.0 * v4rhosigma2tau_abba
                             + 8.0 * v4rhosigma2tau_aacb + 8.0 * v4rhosigma2tau_aaca + 8.0 * v4rhosigma2tau_aabb
                             + 8.0 * v4rhosigma2tau_aaba + 4.0 * v4rhosigma2tau_aaab + 4.0 * v4rhosigma2tau_aaaa;
                double rxll  = 2.0 * v4rhosigmalapl2_acbb + 4.0 * v4rhosigmalapl2_acab + 2.0 * v4rhosigmalapl2_acaa
                             + 2.0 * v4rhosigmalapl2_abbb + 4.0 * v4rhosigmalapl2_abab + 2.0 * v4rhosigmalapl2_abaa
                             + 2.0 * v4rhosigmalapl2_aabb + 4.0 * v4rhosigmalapl2_aaab + 2.0 * v4rhosigmalapl2_aaaa;
                double rxlt  = 2.0 * v4rhosigmalapltau_acbb + 2.0 * v4rhosigmalapltau_acba + 2.0 * v4rhosigmalapltau_acab
                             + 2.0 * v4rhosigmalapltau_acaa + 2.0 * v4rhosigmalapltau_abbb + 2.0 * v4rhosigmalapltau_abba
                             + 2.0 * v4rhosigmalapltau_abab + 2.0 * v4rhosigmalapltau_abaa + 2.0 * v4rhosigmalapltau_aabb
                             + 2.0 * v4rhosigmalapltau_aaba + 2.0 * v4rhosigmalapltau_aaab + 2.0 * v4rhosigmalapltau_aaaa;
                double rxtt  = 2.0 * v4rhosigmatau2_acbb + 4.0 * v4rhosigmatau2_acab + 2.0 * v4rhosigmatau2_acaa
                             + 2.0 * v4rhosigmatau2_abbb + 4.0 * v4rhosigmatau2_abab + 2.0 * v4rhosigmatau2_abaa
                             + 2.0 * v4rhosigmatau2_aabb + 4.0 * v4rhosigmatau2_aaab + 2.0 * v4rhosigmatau2_aaaa;
                double rlll  = v4rholapl3_abbb + 3.0 * v4rholapl3_aabb + 3.0 * v4rholapl3_aaab + v4rholapl3_aaaa;
                double rllt  = v4rholapl2tau_abbb + v4rholapl2tau_abba
                             + 2.0 * v4rholapl2tau_aabb + 2.0 * v4rholapl2tau_aaba
                            + v4rholapl2tau_aaab + v4rholapl2tau_aaaa;
                double rltt  = v4rholapltau2_abbb + 2.0 * v4rholapltau2_abab + v4rholapltau2_abaa
                             + v4rholapltau2_aabb + 2.0 * v4rholapltau2_aaab + v4rholapltau2_aaaa;
                double rttt  = v4rhotau3_abbb + 3.0 * v4rhotau3_aabb + 3.0 * v4rhotau3_aaab + v4rhotau3_aaaa;
                double xxxx  = 8.0*v4sigma4_cccc + 24.0*v4sigma4_cccb + 24.0*v4sigma4_ccbb + 8.0*v4sigma4_cbbb + 40.0*v4sigma4_accc
                             + 96.0*v4sigma4_accb + 72.0*v4sigma4_acbb + 16.0*v4sigma4_abbb + 72.0*v4sigma4_aacc + 120.0*v4sigma4_aacb
                             + 48.0*v4sigma4_aabb + 56.0*v4sigma4_aaac + 48.0*v4sigma4_aaab + 16.0*v4sigma4_aaaa;
                double xxxr  = 4.0 * v4rhosigma3_bccc + 8.0 * v4rhosigma3_bccb + 4.0 * v4rhosigma3_bcbb + 16.0 * v4rhosigma3_bacc
                            + 24.0 * v4rhosigma3_bacb + 8.0 * v4rhosigma3_babb + 20.0 * v4rhosigma3_baac + 16.0 * v4rhosigma3_baab
                            + 8.0 * v4rhosigma3_baaa + 4.0 * v4rhosigma3_accc + 8.0 * v4rhosigma3_accb + 4.0 * v4rhosigma3_acbb
                            + 16.0 * v4rhosigma3_aacc + 24.0 * v4rhosigma3_aacb + 8.0 * v4rhosigma3_aabb + 20.0 * v4rhosigma3_aaac
                            + 16.0 * v4rhosigma3_aaab + 8.0 * v4rhosigma3_aaaa;
                double xxxl  = 4.0 * v4sigma3lapl_cccb + 4.0 * v4sigma3lapl_ccca + 8.0 * v4sigma3lapl_ccbb + 8.0 * v4sigma3lapl_ccba
                             + 4.0 * v4sigma3lapl_cbbb + 4.0 * v4sigma3lapl_cbba + 16.0 * v4sigma3lapl_accb + 16.0 * v4sigma3lapl_acca
                             + 24.0 * v4sigma3lapl_acbb + 24.0 * v4sigma3lapl_acba + 8.0 * v4sigma3lapl_abbb + 8.0 * v4sigma3lapl_abba
                             + 20.0 * v4sigma3lapl_aacb + 20.0 * v4sigma3lapl_aaca + 16.0 * v4sigma3lapl_aabb + 16.0 * v4sigma3lapl_aaba
                             + 8.0 * v4sigma3lapl_aaab + 8.0 * v4sigma3lapl_aaaa;
                double xxxt  = 4.0 * v4sigma3tau_cccb + 4.0 * v4sigma3tau_ccca + 8.0 * v4sigma3tau_ccbb
                             + 8.0 * v4sigma3tau_ccba + 4.0 * v4sigma3tau_cbbb + 4.0 * v4sigma3tau_cbba
                            + 16.0 * v4sigma3tau_accb + 16.0 * v4sigma3tau_acca + 24.0 * v4sigma3tau_acbb
                            + 24.0 * v4sigma3tau_acba + 8.0 * v4sigma3tau_abbb + 8.0 * v4sigma3tau_abba
                            + 20.0 * v4sigma3tau_aacb + 20.0 * v4sigma3tau_aaca + 16.0 * v4sigma3tau_aabb
                            + 16.0 * v4sigma3tau_aaba + 8.0 * v4sigma3tau_aaab + 8.0 * v4sigma3tau_aaaa;
                double xxrr  = 2.0 * v4rho2sigma2_bbcc + 2.0 * v4rho2sigma2_bbcb + 6.0 * v4rho2sigma2_bbac
                             + 4.0 * v4rho2sigma2_bbab + 4.0 * v4rho2sigma2_bbaa + 4.0 * v4rho2sigma2_abcc
                             + 4.0 * v4rho2sigma2_abcb + 12.0 * v4rho2sigma2_abac + 8.0 * v4rho2sigma2_abab
                             + 8.0 * v4rho2sigma2_abaa + 2.0 * v4rho2sigma2_aacc + 2.0 * v4rho2sigma2_aacb
                             + 6.0 * v4rho2sigma2_aaac + 4.0 * v4rho2sigma2_aaab + 4.0 * v4rho2sigma2_aaaa;
                double xxrl  = 2.0 * v4rhosigma2lapl_bccb + 2.0 * v4rhosigma2lapl_bcca + 2.0 * v4rhosigma2lapl_bcbb
                             + 2.0 * v4rhosigma2lapl_bcba + 6.0 * v4rhosigma2lapl_bacb + 6.0 * v4rhosigma2lapl_baca
                             + 4.0 * v4rhosigma2lapl_babb + 4.0 * v4rhosigma2lapl_baba + 4.0 * v4rhosigma2lapl_baab
                             + 4.0 * v4rhosigma2lapl_baaa + 2.0 * v4rhosigma2lapl_accb + 2.0 * v4rhosigma2lapl_acca
                             + 2.0 * v4rhosigma2lapl_acbb + 2.0 * v4rhosigma2lapl_acba + 6.0 * v4rhosigma2lapl_aacb
                             + 6.0 * v4rhosigma2lapl_aaca + 4.0 * v4rhosigma2lapl_aabb + 4.0 * v4rhosigma2lapl_aaba
                             + 4.0 * v4rhosigma2lapl_aaab + 4.0 * v4rhosigma2lapl_aaaa;
                double xxrt  = 2.0 * v4rhosigma2tau_bccb + 2.0 * v4rhosigma2tau_bcca + 2.0 * v4rhosigma2tau_bcbb + 2.0 * v4rhosigma2tau_bcba + 6.0 * v4rhosigma2tau_bacb + 6.0 * v4rhosigma2tau_baca + 4.0 * v4rhosigma2tau_babb + 4.0 * v4rhosigma2tau_baba + 4.0 * v4rhosigma2tau_baab + 4.0 * v4rhosigma2tau_baaa + 2.0 * v4rhosigma2tau_accb + 2.0 * v4rhosigma2tau_acca + 2.0 * v4rhosigma2tau_acbb + 2.0 * v4rhosigma2tau_acba + 6.0 * v4rhosigma2tau_aacb + 6.0 * v4rhosigma2tau_aaca + 4.0 * v4rhosigma2tau_aabb + 4.0 * v4rhosigma2tau_aaba + 4.0 * v4rhosigma2tau_aaab + 4.0 * v4rhosigma2tau_aaaa;
                double xxll  = 2.0 * v4sigma2lapl2_ccbb + 4.0 * v4sigma2lapl2_ccab + 2.0 * v4sigma2lapl2_ccaa + 2.0 * v4sigma2lapl2_cbbb + 4.0 * v4sigma2lapl2_cbab + 2.0 * v4sigma2lapl2_cbaa + 6.0 * v4sigma2lapl2_acbb + 12.0 * v4sigma2lapl2_acab + 6.0 * v4sigma2lapl2_acaa + 4.0 * v4sigma2lapl2_abbb + 8.0 * v4sigma2lapl2_abab + 4.0 * v4sigma2lapl2_abaa + 4.0 * v4sigma2lapl2_aabb + 8.0 * v4sigma2lapl2_aaab + 4.0 * v4sigma2lapl2_aaaa;
                double xxlt  = 2.0 * v4sigma2lapltau_ccbb + 2.0 * v4sigma2lapltau_ccba + 2.0 * v4sigma2lapltau_ccab + 2.0 * v4sigma2lapltau_ccaa + 2.0 * v4sigma2lapltau_cbbb + 2.0 * v4sigma2lapltau_cbba + 2.0 * v4sigma2lapltau_cbab + 2.0 * v4sigma2lapltau_cbaa + 6.0 * v4sigma2lapltau_acbb + 6.0 * v4sigma2lapltau_acba + 6.0 * v4sigma2lapltau_acab + 6.0 * v4sigma2lapltau_acaa + 4.0 * v4sigma2lapltau_abbb + 4.0 * v4sigma2lapltau_abba + 4.0 * v4sigma2lapltau_abab + 4.0 * v4sigma2lapltau_abaa + 4.0 * v4sigma2lapltau_aabb + 4.0 * v4sigma2lapltau_aaba + 4.0 * v4sigma2lapltau_aaab + 4.0 * v4sigma2lapltau_aaaa;
                double xxtt  = 2.0 * v4sigma2tau2_ccbb + 4.0 * v4sigma2tau2_ccab + 2.0 * v4sigma2tau2_ccaa + 2.0 * v4sigma2tau2_cbbb + 4.0 * v4sigma2tau2_cbab + 2.0 * v4sigma2tau2_cbaa + 6.0 * v4sigma2tau2_acbb + 12.0 * v4sigma2tau2_acab + 6.0 * v4sigma2tau2_acaa + 4.0 * v4sigma2tau2_abbb + 8.0 * v4sigma2tau2_abab + 4.0 * v4sigma2tau2_abaa + 4.0 * v4sigma2tau2_aabb + 8.0 * v4sigma2tau2_aaab + 4.0 * v4sigma2tau2_aaaa;
                double xrrr  = v4rho3sigma_bbbc + 2.0 * v4rho3sigma_bbba + 3.0 * v4rho3sigma_abbc + 6.0 * v4rho3sigma_abba + 3.0 * v4rho3sigma_aabc
                             + 6.0 * v4rho3sigma_aaba + v4rho3sigma_aaac + 2.0 * v4rho3sigma_aaaa;
                double xrrl  = v4rho2sigmalapl_bbcb + v4rho2sigmalapl_bbca + 2.0 * v4rho2sigmalapl_bbab + 2.0 * v4rho2sigmalapl_bbaa + 2.0 * v4rho2sigmalapl_abcb + 2.0 * v4rho2sigmalapl_abca + 4.0 * v4rho2sigmalapl_abab + 4.0 * v4rho2sigmalapl_abaa + v4rho2sigmalapl_aacb + v4rho2sigmalapl_aaca + 2.0 * v4rho2sigmalapl_aaab + 2.0 * v4rho2sigmalapl_aaaa;
                double xrrt  = v4rho2sigmatau_bbcb + v4rho2sigmatau_bbca + 2.0 * v4rho2sigmatau_bbab + 2.0 * v4rho2sigmatau_bbaa + 2.0 * v4rho2sigmatau_abcb + 2.0 * v4rho2sigmatau_abca + 4.0 * v4rho2sigmatau_abab + 4.0 * v4rho2sigmatau_abaa + v4rho2sigmatau_aacb + v4rho2sigmatau_aaca + 2.0 * v4rho2sigmatau_aaab + 2.0 * v4rho2sigmatau_aaaa;
                double xrll  = v4rhosigmalapl2_bcbb + 2.0 * v4rhosigmalapl2_bcab + v4rhosigmalapl2_bcaa + 2.0 * v4rhosigmalapl2_babb + 4.0 * v4rhosigmalapl2_baab + 2.0 * v4rhosigmalapl2_baaa + v4rhosigmalapl2_acbb + 2.0 * v4rhosigmalapl2_acab + v4rhosigmalapl2_acaa + 2.0 * v4rhosigmalapl2_aabb + 4.0 * v4rhosigmalapl2_aaab + 2.0 * v4rhosigmalapl2_aaaa;
                double xrlt  = v4rhosigmalapltau_bcbb + v4rhosigmalapltau_bcba + v4rhosigmalapltau_bcab + v4rhosigmalapltau_bcaa + 2.0 * v4rhosigmalapltau_babb + 2.0 * v4rhosigmalapltau_baba + 2.0 * v4rhosigmalapltau_baab + 2.0 * v4rhosigmalapltau_baaa + v4rhosigmalapltau_acbb + v4rhosigmalapltau_acba + v4rhosigmalapltau_acab + v4rhosigmalapltau_acaa + 2.0 * v4rhosigmalapltau_aabb + 2.0 * v4rhosigmalapltau_aaba + 2.0 * v4rhosigmalapltau_aaab + 2.0 * v4rhosigmalapltau_aaaa;
                double xrtt  = v4rhosigmatau2_bcbb + 2.0 * v4rhosigmatau2_bcab + v4rhosigmatau2_bcaa + 2.0 * v4rhosigmatau2_babb + 4.0 * v4rhosigmatau2_baab + 2.0 * v4rhosigmatau2_baaa + v4rhosigmatau2_acbb + 2.0 * v4rhosigmatau2_acab + v4rhosigmatau2_acaa + 2.0 * v4rhosigmatau2_aabb + 4.0 * v4rhosigmatau2_aaab + 2.0 * v4rhosigmatau2_aaaa;
                double xlll  = v4sigmalapl3_cbbb + 3.0 * v4sigmalapl3_cabb + 3.0 * v4sigmalapl3_caab + v4sigmalapl3_caaa + 2.0 * v4sigmalapl3_abbb + 6.0 * v4sigmalapl3_aabb + 6.0 * v4sigmalapl3_aaab + 2.0 * v4sigmalapl3_aaaa;
                double xllt  = v4sigmalapl2tau_cbbb + v4sigmalapl2tau_cbba + 2.0 * v4sigmalapl2tau_cabb + 2.0 * v4sigmalapl2tau_caba
                             + v4sigmalapl2tau_caab + v4sigmalapl2tau_caaa + 2.0 * v4sigmalapl2tau_abbb + 2.0 * v4sigmalapl2tau_abba
                             + 4.0 * v4sigmalapl2tau_aabb + 4.0 * v4sigmalapl2tau_aaba + 2.0 * v4sigmalapl2tau_aaab + 2.0 * v4sigmalapl2tau_aaaa;
                double xltt  = v4sigmalapltau2_cbbb + 2.0 * v4sigmalapltau2_cbab + v4sigmalapltau2_cbaa + v4sigmalapltau2_cabb + 2.0 * v4sigmalapltau2_caab + v4sigmalapltau2_caaa + 2.0 * v4sigmalapltau2_abbb + 4.0 * v4sigmalapltau2_abab + 2.0 * v4sigmalapltau2_abaa + 2.0 * v4sigmalapltau2_aabb + 4.0 * v4sigmalapltau2_aaab + 2.0 * v4sigmalapltau2_aaaa;
                double xttt  = v4sigmatau3_cbbb + 3.0 * v4sigmatau3_cabb + 3.0 * v4sigmatau3_caab + v4sigmatau3_caaa + 2.0 * v4sigmatau3_abbb + 6.0 * v4sigmatau3_aabb + 6.0 * v4sigmatau3_aaab + 2.0 * v4sigmatau3_aaaa;
                double llll  = v4lapl4_abbb + 3.0 * v4lapl4_aabb + 3.0 * v4lapl4_aaab + v4lapl4_aaaa;
                double lllr  = v4rholapl3_babb + 2.0 * v4rholapl3_baab + v4rholapl3_baaa + v4rholapl3_aabb + 2.0 * v4rholapl3_aaab + v4rholapl3_aaaa;
                double lllx  = 2.0 * v4sigmalapl3_cabb + 4.0 * v4sigmalapl3_caab + 2.0 * v4sigmalapl3_caaa + 2.0 * v4sigmalapl3_babb + 4.0 * v4sigmalapl3_baab + 2.0 * v4sigmalapl3_baaa + 2.0 * v4sigmalapl3_aabb + 4.0 * v4sigmalapl3_aaab + 2.0 * v4sigmalapl3_aaaa;
                double lllt  = v4lapl3tau_abbb + v4lapl3tau_abba + 2.0 * v4lapl3tau_aabb + 2.0 * v4lapl3tau_aaba + v4lapl3tau_aaab + v4lapl3tau_aaaa;
                double llrr  = v4rho2lapl2_bbab + v4rho2lapl2_bbaa + 2.0 * v4rho2lapl2_abab + 2.0 * v4rho2lapl2_abaa + v4rho2lapl2_aaab + v4rho2lapl2_aaaa;
                double llrx  = 2.0 * v4rhosigmalapl2_bcab + 2.0 * v4rhosigmalapl2_bcaa + 2.0 * v4rhosigmalapl2_bbab + 2.0 * v4rhosigmalapl2_bbaa + 2.0 * v4rhosigmalapl2_baab + 2.0 * v4rhosigmalapl2_baaa + 2.0 * v4rhosigmalapl2_acab + 2.0 * v4rhosigmalapl2_acaa + 2.0 * v4rhosigmalapl2_abab + 2.0 * v4rhosigmalapl2_abaa + 2.0 * v4rhosigmalapl2_aaab + 2.0 * v4rhosigmalapl2_aaaa;
                double llrt  = v4rholapl2tau_babb + v4rholapl2tau_baba + v4rholapl2tau_baab + v4rholapl2tau_baaa + v4rholapl2tau_aabb + v4rholapl2tau_aaba + v4rholapl2tau_aaab + v4rholapl2tau_aaaa;
                double llxx  = 4.0 * v4sigma2lapl2_ccab + 4.0 * v4sigma2lapl2_ccaa + 8.0 * v4sigma2lapl2_cbab + 8.0 * v4sigma2lapl2_cbaa + 4.0 * v4sigma2lapl2_bbab + 4.0 * v4sigma2lapl2_bbaa + 8.0 * v4sigma2lapl2_acab + 8.0 * v4sigma2lapl2_acaa + 8.0 * v4sigma2lapl2_abab + 8.0 * v4sigma2lapl2_abaa + 4.0 * v4sigma2lapl2_aaab + 4.0 * v4sigma2lapl2_aaaa;
                double llxt  = 2.0 * v4sigmalapl2tau_cabb + 2.0 * v4sigmalapl2tau_caba + 2.0 * v4sigmalapl2tau_caab + 2.0 * v4sigmalapl2tau_caaa + 2.0 * v4sigmalapl2tau_babb + 2.0 * v4sigmalapl2tau_baba + 2.0 * v4sigmalapl2tau_baab + 2.0 * v4sigmalapl2tau_baaa + 2.0 * v4sigmalapl2tau_aabb + 2.0 * v4sigmalapl2tau_aaba + 2.0 * v4sigmalapl2tau_aaab + 2.0 * v4sigmalapl2tau_aaaa;
                double lltt  = v4lapl2tau2_abbb + 2.0 * v4lapl2tau2_abab + v4lapl2tau2_abaa + v4lapl2tau2_aabb + 2.0 * v4lapl2tau2_aaab + v4lapl2tau2_aaaa;
                double lrrr  = v4rho3lapl_bbba + 3.0 * v4rho3lapl_abba + 3.0 * v4rho3lapl_aaba + v4rho3lapl_aaaa;
                double lrrx  = 2.0 * v4rho2sigmalapl_bbca + 2.0 * v4rho2sigmalapl_bbba + 2.0 * v4rho2sigmalapl_bbaa + 4.0 * v4rho2sigmalapl_abca + 4.0 * v4rho2sigmalapl_abba + 4.0 * v4rho2sigmalapl_abaa + 2.0 * v4rho2sigmalapl_aaca + 2.0 * v4rho2sigmalapl_aaba + 2.0 * v4rho2sigmalapl_aaaa;
                double lrrt  = v4rho2lapltau_bbab + v4rho2lapltau_bbaa + 2.0 * v4rho2lapltau_abab + 2.0 * v4rho2lapltau_abaa + v4rho2lapltau_aaab + v4rho2lapltau_aaaa;
                double lrxx  = 4.0 * v4rhosigma2lapl_bcca + 8.0 * v4rhosigma2lapl_bcba + 4.0 * v4rhosigma2lapl_bbba + 8.0 * v4rhosigma2lapl_baca
                             + 8.0 * v4rhosigma2lapl_baba + 4.0 * v4rhosigma2lapl_baaa + 4.0 * v4rhosigma2lapl_acca + 8.0 * v4rhosigma2lapl_acba
                             + 4.0 * v4rhosigma2lapl_abba + 8.0 * v4rhosigma2lapl_aaca + 8.0 * v4rhosigma2lapl_aaba + 4.0 * v4rhosigma2lapl_aaaa;
                double lrxt  = 2.0 * v4rhosigmalapltau_bcab + 2.0 * v4rhosigmalapltau_bcaa + 2.0 * v4rhosigmalapltau_bbab + 2.0 * v4rhosigmalapltau_bbaa + 2.0 * v4rhosigmalapltau_baab + 2.0 * v4rhosigmalapltau_baaa + 2.0 * v4rhosigmalapltau_acab + 2.0 * v4rhosigmalapltau_acaa + 2.0 * v4rhosigmalapltau_abab + 2.0 * v4rhosigmalapltau_abaa + 2.0 * v4rhosigmalapltau_aaab + 2.0 * v4rhosigmalapltau_aaaa;
                double lrtt  = v4rholapltau2_babb + 2.0 * v4rholapltau2_baab + v4rholapltau2_baaa + v4rholapltau2_aabb + 2.0 * v4rholapltau2_aaab + v4rholapltau2_aaaa;
                double lxxx  = 8.0 * v4sigma3lapl_ccca + 24.0 * v4sigma3lapl_ccba + 24.0 * v4sigma3lapl_cbba + 8.0 * v4sigma3lapl_bbba + 24.0 * v4sigma3lapl_acca + 48.0 * v4sigma3lapl_acba + 24.0 * v4sigma3lapl_abba + 24.0 * v4sigma3lapl_aaca + 24.0 * v4sigma3lapl_aaba + 8.0 * v4sigma3lapl_aaaa;
                double lxxt  = 4.0 * v4sigma2lapltau_ccab + 4.0 * v4sigma2lapltau_ccaa + 8.0 * v4sigma2lapltau_cbab + 8.0 * v4sigma2lapltau_cbaa + 4.0 * v4sigma2lapltau_bbab + 4.0 * v4sigma2lapltau_bbaa + 8.0 * v4sigma2lapltau_acab + 8.0 * v4sigma2lapltau_acaa + 8.0 * v4sigma2lapltau_abab + 8.0 * v4sigma2lapltau_abaa + 4.0 * v4sigma2lapltau_aaab + 4.0 * v4sigma2lapltau_aaaa;
                double lxtt  = 2.0 * v4sigmalapltau2_cabb + 4.0 * v4sigmalapltau2_caab + 2.0 * v4sigmalapltau2_caaa + 2.0 * v4sigmalapltau2_babb
                             + 4.0 * v4sigmalapltau2_baab + 2.0 * v4sigmalapltau2_baaa + 2.0 * v4sigmalapltau2_aabb + 4.0 * v4sigmalapltau2_aaab
                             + 2.0 * v4sigmalapltau2_aaaa;
                double lttt  = v4lapltau3_abbb + 3.0 * v4lapltau3_aabb + 3.0 * v4lapltau3_aaab + v4lapltau3_aaaa;
                double tttt  = v4tau4_abbb + 3.0 * v4tau4_aabb + 3.0 * v4tau4_aaab + v4tau4_aaaa;
                double tttr  = v4rhotau3_babb + 2.0 * v4rhotau3_baab + v4rhotau3_baaa + v4rhotau3_aabb + 2.0 * v4rhotau3_aaab + v4rhotau3_aaaa;
                double tttx  = 2.0 * v4sigmatau3_cabb + 4.0 * v4sigmatau3_caab + 2.0 * v4sigmatau3_caaa + 2.0 * v4sigmatau3_babb + 4.0 * v4sigmatau3_baab + 2.0 * v4sigmatau3_baaa + 2.0 * v4sigmatau3_aabb + 4.0 * v4sigmatau3_aaab + 2.0 * v4sigmatau3_aaaa;
                double tttl  = v4lapltau3_babb + 2.0 * v4lapltau3_baab + v4lapltau3_baaa + v4lapltau3_aabb + 2.0 * v4lapltau3_aaab + v4lapltau3_aaaa;
                double ttrr  = v4rho2tau2_bbab + v4rho2tau2_bbaa + 2.0 * v4rho2tau2_abab + 2.0 * v4rho2tau2_abaa + v4rho2tau2_aaab + v4rho2tau2_aaaa;
                double ttrx  = 2.0 * v4rhosigmatau2_bcab + 2.0 * v4rhosigmatau2_bcaa + 2.0 * v4rhosigmatau2_bbab + 2.0 * v4rhosigmatau2_bbaa + 2.0 * v4rhosigmatau2_baab + 2.0 * v4rhosigmatau2_baaa + 2.0 * v4rhosigmatau2_acab + 2.0 * v4rhosigmatau2_acaa + 2.0 * v4rhosigmatau2_abab + 2.0 * v4rhosigmatau2_abaa + 2.0 * v4rhosigmatau2_aaab + 2.0 * v4rhosigmatau2_aaaa;
                double ttrl  = v4rholapltau2_bbab + v4rholapltau2_bbaa + v4rholapltau2_baab + v4rholapltau2_baaa + v4rholapltau2_abab + v4rholapltau2_abaa + v4rholapltau2_aaab + v4rholapltau2_aaaa;
                double ttxx  = 4.0 * v4sigma2tau2_ccab + 4.0 * v4sigma2tau2_ccaa + 8.0 * v4sigma2tau2_cbab + 8.0 * v4sigma2tau2_cbaa + 4.0 * v4sigma2tau2_bbab + 4.0 * v4sigma2tau2_bbaa + 8.0 * v4sigma2tau2_acab + 8.0 * v4sigma2tau2_acaa + 8.0 * v4sigma2tau2_abab + 8.0 * v4sigma2tau2_abaa + 4.0 * v4sigma2tau2_aaab + 4.0 * v4sigma2tau2_aaaa;
                double ttxl  = 2.0 * v4sigmalapltau2_cbab + 2.0 * v4sigmalapltau2_cbaa + 2.0 * v4sigmalapltau2_caab + 2.0 * v4sigmalapltau2_caaa + 2.0 * v4sigmalapltau2_bbab + 2.0 * v4sigmalapltau2_bbaa + 2.0 * v4sigmalapltau2_baab + 2.0 * v4sigmalapltau2_baaa + 2.0 * v4sigmalapltau2_abab + 2.0 * v4sigmalapltau2_abaa + 2.0 * v4sigmalapltau2_aaab + 2.0 * v4sigmalapltau2_aaaa;
                double ttll  = v4lapl2tau2_bbab + v4lapl2tau2_bbaa + 2.0 * v4lapl2tau2_abab + 2.0 * v4lapl2tau2_abaa + v4lapl2tau2_aaab + v4lapl2tau2_aaaa;
                double trrr  = v4rho3tau_bbba + 3.0 * v4rho3tau_abba + 3.0 * v4rho3tau_aaba + v4rho3tau_aaaa;
                double trrx  = 2.0 * v4rho2sigmatau_bbca + 2.0 * v4rho2sigmatau_bbba + 2.0 * v4rho2sigmatau_bbaa + 4.0 * v4rho2sigmatau_abca + 4.0 * v4rho2sigmatau_abba + 4.0 * v4rho2sigmatau_abaa + 2.0 * v4rho2sigmatau_aaca + 2.0 * v4rho2sigmatau_aaba + 2.0 * v4rho2sigmatau_aaaa;
                double trrl  = v4rho2lapltau_bbba + v4rho2lapltau_bbaa + 2.0 * v4rho2lapltau_abba + 2.0 * v4rho2lapltau_abaa + v4rho2lapltau_aaba + v4rho2lapltau_aaaa;
                double trxx  = 4.0 * v4rhosigma2tau_bcca + 8.0 * v4rhosigma2tau_bcba + 4.0 * v4rhosigma2tau_bbba + 8.0 * v4rhosigma2tau_baca + 8.0 * v4rhosigma2tau_baba + 4.0 * v4rhosigma2tau_baaa + 4.0 * v4rhosigma2tau_acca + 8.0 * v4rhosigma2tau_acba + 4.0 * v4rhosigma2tau_abba + 8.0 * v4rhosigma2tau_aaca + 8.0 * v4rhosigma2tau_aaba + 4.0 * v4rhosigma2tau_aaaa;
                double trxl  = 2.0 * v4rhosigmalapltau_bcba + 2.0 * v4rhosigmalapltau_bcaa + 2.0 * v4rhosigmalapltau_bbba + 2.0 * v4rhosigmalapltau_bbaa + 2.0 * v4rhosigmalapltau_baba + 2.0 * v4rhosigmalapltau_baaa + 2.0 * v4rhosigmalapltau_acba + 2.0 * v4rhosigmalapltau_acaa + 2.0 * v4rhosigmalapltau_abba + 2.0 * v4rhosigmalapltau_abaa + 2.0 * v4rhosigmalapltau_aaba + 2.0 * v4rhosigmalapltau_aaaa;
                double trll  = v4rholapl2tau_bbba + 2.0 * v4rholapl2tau_baba + v4rholapl2tau_baaa + v4rholapl2tau_abba + 2.0 * v4rholapl2tau_aaba + v4rholapl2tau_aaaa;
                double txxx  = 8.0 * v4sigma3tau_ccca + 24.0 * v4sigma3tau_ccba + 24.0 * v4sigma3tau_cbba + 8.0 * v4sigma3tau_bbba + 24.0 * v4sigma3tau_acca + 48.0 * v4sigma3tau_acba + 24.0 * v4sigma3tau_abba + 24.0 * v4sigma3tau_aaca + 24.0 * v4sigma3tau_aaba + 8.0 * v4sigma3tau_aaaa;
                double txxl  = 4.0 * v4sigma2lapltau_ccba + 4.0 * v4sigma2lapltau_ccaa + 8.0 * v4sigma2lapltau_cbba + 8.0 * v4sigma2lapltau_cbaa + 4.0 * v4sigma2lapltau_bbba + 4.0 * v4sigma2lapltau_bbaa + 8.0 * v4sigma2lapltau_acba + 8.0 * v4sigma2lapltau_acaa + 8.0 * v4sigma2lapltau_abba + 8.0 * v4sigma2lapltau_abaa + 4.0 * v4sigma2lapltau_aaba + 4.0 * v4sigma2lapltau_aaaa;
                double txll  = 2.0 * v4sigmalapl2tau_cbba + 4.0 * v4sigmalapl2tau_caba + 2.0 * v4sigmalapl2tau_caaa + 2.0 * v4sigmalapl2tau_bbba + 4.0 * v4sigmalapl2tau_baba + 2.0 * v4sigmalapl2tau_baaa + 2.0 * v4sigmalapl2tau_abba + 4.0 * v4sigmalapl2tau_aaba + 2.0 * v4sigmalapl2tau_aaaa;
                double tlll  = v4lapl3tau_bbba + 3.0 * v4lapl3tau_abba + 3.0 * v4lapl3tau_aaba + v4lapl3tau_aaaa;

                // Scalar contribution

                double tau_0 = 0.0;
                double rho_0 = 0.0;
                //double lap_0 = 0.0;

                // vxc 1 contributions

                rho_0 +=  rr * rhow[g]
                        + rx * l2contract
                        + rt * tauw[g]
                        + rl * laplw[g];

                //lap_0 +=   lr * rhow[g]
                //         + lx * l2contract
                //         + lt * tauw[g]
                //         + ll * laplw[g];

                tau_0 +=    tr * rhow[g]
                          + tx * l2contract
                          + tt * tauw[g]
                          + tl * laplw[g];

                // vxc 2 contributions

                rho_0 += rrr * gam[g]
                       + rrt * rt_gam[g]
                       + rrl * rl_gam[g]
                       + rll * ll_gam[g]
                       + rtt * tt_gam[g]
                       + rtl * tl_gam[g]
                       + rrx * q2contract
                       + rlx * sl_q2contract
                       + rtx * st_q2contract
                       + rxx * q3contract
                       + rx  * q4contract;

                //lap_0 += lrr * gam[g]
                //       + lrt * rt_gam[g]
                //       + llr * rl_gam[g]
                //       + lll * ll_gam[g]
                //       + ltt * tt_gam[g]
                //       + llt * tl_gam[g]
                //       + lrx * q2contract
                //       + llx * sl_q2contract
                //       + ltx * st_q2contract
                //       + lxx * q3contract
                //       + lx  * q4contract;

                tau_0 += trr * gam[g]
                       + ttr * rt_gam[g]
                       + trl * rl_gam[g]
                       + tll * ll_gam[g]
                       + ttt * tt_gam[g]
                       + ttl * tl_gam[g]
                       + trx * q2contract
                       + tlx * sl_q2contract
                       + ttx * st_q2contract
                       + txx * q3contract
                       + tx  * q4contract;

                // vxc 3 contributions

                rho_0 += + rrrr * pi[g]
                         + rrrt * rrt_pi[g]
                         + rrrl * rrl_pi[g]
                         + rrtt * rtt_pi[g]
                         + rrlt * rtl_pi[g]
                         + rrll * rll_pi[g]
                         + rttt * ttt_pi[g]
                         + rltt * ttl_pi[g]
                         + rllt * tll_pi[g]
                         + rlll * lll_pi[g]

                         + rrrx * c2
                         + rrxt * rt_c2
                         + rrxl * rl_c2
                         + rxll * ll_c2
                         + rxtt * tt_c2
                         + rxlt * tl_c2

                         + rrxx * c3
                         + rxxl * l_c3
                         + rxxt * t_c3

                         + rrx * c4
                         + rlx * l_c4
                         + rtx * t_c4

                         + rxx * (c5_6 + c8)
                         + rxxx * c7;

                //lap_0 += + lrrr * pi[g]
                //         + lrrt * rrt_pi[g]
                //         + llrr * rrl_pi[g]
                //         + lrtt * rtt_pi[g]
                //         + llrt * rtl_pi[g]
                //         + lllr * rll_pi[g]
                //         + lttt * ttt_pi[g]
                //         + lltt * ttl_pi[g]
                //         + lllt * tll_pi[g]
                //         + llll * lll_pi[g]
                //
                //         + lrrx * c2
                //         + lrxt * rt_c2
                //         + llrx * rl_c2
                //         + lllx * ll_c2
                //         + lxtt * tt_c2
                //         + llxt * tl_c2
                //
                //         + lrxx * c3
                //         + llxx * l_c3
                //         + lxxt * t_c3
                //
                //         + lrx * c4
                //         + llx * l_c4
                //         + ltx * t_c4
                //
                //         + lxx * (c5_6 + c8)
                //         + lxxx * c7;

                tau_0 += tttt * ttt_pi[g]
                       + tttr * rtt_pi[g]
                       + tttl * ttl_pi[g]
                       + ttrr * rrt_pi[g]
                       + ttrl * rtl_pi[g]
                       + ttll * tll_pi[g]
                       + trrr * pi[g]
                       + trrl * rrl_pi[g]
                       + trll * rll_pi[g]

                       + trrx *c2
                       + ttrx *rt_c2
                       + trxl *rl_c2
                       + txll *ll_c2
                       + tttx *tt_c2
                       + ttxl *tl_c2

                       + trxx * c3
                       + txxl * l_c3
                       + ttxx * t_c3

                       + trx * c4
                       + tlx * l_c4
                       + ttx * t_c4

                       + txx * (c5_6 + c8)
                       + txxx * c7;


                // vector contribution

                double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                // vxc 1 contribution

                xcomp +=  xr * grada_x_g * rhow[g]
                        + xt * grada_x_g * tauw[g]
                        + xl * grada_x_g * laplw[g]
                        +  x * gradw_x[g]
                        + xx * l5contract_x;

                ycomp +=  xr * grada_y_g * rhow[g]
                        + xt * grada_y_g * tauw[g]
                        + xl * grada_y_g * laplw[g]
                        +  x * gradw_y[g]
                        + xx * l5contract_y;

                zcomp +=  xr * grada_z_g * rhow[g]
                        + xt * grada_z_g * tauw[g]
                        + xl * grada_z_g * laplw[g]
                        +  x * gradw_z[g]
                        + xx * l5contract_z;

                // // vxc 2 contributions

                xcomp +=  xrr * grada_x_g *    gam[g]
                        + xrt * grada_x_g * rt_gam[g]
                        + xrl * grada_x_g * rl_gam[g]
                        + xll * grada_x_g * ll_gam[g]
                        + xtt * grada_x_g * tt_gam[g]
                        + xtl * grada_x_g * tl_gam[g]
                        + xr * gamx[g] // q6
                        + xl * sl_gamx[g]
                        + xt * st_gamx[g]
                        + xxr * q7contract_x
                        + xxl * sl_q7contract_x
                        + xxt * st_q7contract_x
                        + xx * (q8contract_x + q10contract_x + q11contract_x)
                        + xxx * q9contract_x;

                ycomp +=  xrr * grada_y_g *    gam[g] // q5
                        + xrt * grada_y_g * rt_gam[g]
                        + xrl * grada_y_g * rl_gam[g]
                        + xll * grada_y_g * ll_gam[g]
                        + xtt * grada_y_g * tt_gam[g]
                        + xtl * grada_y_g * tl_gam[g]
                        + xr * gamy[g] // q6
                        + xl * sl_gamy[g]
                        + xt * st_gamy[g]
                        + xxr * q7contract_y
                        + xxl * sl_q7contract_y
                        + xxt * st_q7contract_y
                        + xx * (q8contract_y + q10contract_y + q11contract_y)
                        + xxx * q9contract_y;

                zcomp +=  xrr * grada_z_g *    gam[g] // q5
                        + xrt * grada_z_g * rt_gam[g]
                        + xrl * grada_z_g * rl_gam[g]
                        + xll * grada_z_g * ll_gam[g]
                        + xtt * grada_z_g * tt_gam[g]
                        + xtl * grada_z_g * tl_gam[g]
                        + xr * gamz[g] // q6
                        + xl * sl_gamz[g]
                        + xt * st_gamz[g]
                        + xxr * q7contract_z
                        + xxl * sl_q7contract_z
                        + xxt * st_q7contract_z
                        + xx * (q8contract_z + q10contract_z + q11contract_z)
                        + xxx * q9contract_z;

                // vxc 3 contributions

                xcomp +=  xrrr * grada_x_g * pi[g]  // c9 terms
                        + xrrt * grada_x_g * rrt_pi[g]
                        + xrrl * grada_x_g * rrl_pi[g]
                        + xrtt * grada_x_g * rtt_pi[g]
                        + xrlt * grada_x_g * rtl_pi[g]
                        + xrll * grada_x_g * rll_pi[g]
                        + xttt * grada_x_g * ttt_pi[g]
                        + xltt * grada_x_g * ttl_pi[g]
                        + xllt * grada_x_g * tll_pi[g]
                        + xlll * grada_x_g * lll_pi[g]

                        + xrr * pix[g] // c10 terms
                        + xrt * rt_pix[g]
                        + xrl * rl_pix[g]
                        + xll * ll_pix[g]
                        + xtt * tt_pix[g]
                        + xtl * tl_pix[g]

                        + xxrr *    c2 * grada_x_g   // c11 terms
                        + xxrt * rt_c2 * grada_x_g
                        + xxrl * rl_c2 * grada_x_g
                        + xxll * ll_c2 * grada_x_g
                        + xxtt * tt_c2 * grada_x_g
                        + xxlt * tl_c2 * grada_x_g

                        + xxr * (c12_c14_x + c15_x)
                        + xxl * (l_c12_c14_x + l_c15_x)
                        + xxt * (t_c12_c14_x + t_c15_x)

                        + xxxr * c13_x
                        + xxxl * l_c13_x
                        + xxxt * t_c13_x

                        + xx * c17_24_25_x
                        + xxxx * c18_x
                        + xxx * (c16_19_22_x + c20_21_23_x);

                ycomp +=  xrrr * grada_y_g * pi[g]  // c9 terms
                        + xrrt * grada_y_g * rrt_pi[g]
                        + xrrl * grada_y_g * rrl_pi[g]
                        + xrtt * grada_y_g * rtt_pi[g]
                        + xrlt * grada_y_g * rtl_pi[g]
                        + xrll * grada_y_g * rll_pi[g]
                        + xttt * grada_y_g * ttt_pi[g]
                        + xltt * grada_y_g * ttl_pi[g]
                        + xllt * grada_y_g * tll_pi[g]
                        + xlll * grada_y_g * lll_pi[g]

                        + xrr * piy[g] // c10 terms
                        + xrt * rt_piy[g]
                        + xrl * rl_piy[g]
                        + xll * ll_piy[g]
                        + xtt * tt_piy[g]
                        + xtl * tl_piy[g]

                        + xxrr *    c2 * grada_y_g   // c11 terms
                        + xxrt * rt_c2 * grada_y_g
                        + xxrl * rl_c2 * grada_y_g
                        + xxll * ll_c2 * grada_y_g
                        + xxtt * tt_c2 * grada_y_g
                        + xxlt * tl_c2 * grada_y_g

                        + xxr * (c12_c14_y + c15_y)
                        + xxl * (l_c12_c14_y + l_c15_y)
                        + xxt * (t_c12_c14_y + t_c15_y)

                        + xxxr * c13_y
                        + xxxl * l_c13_y
                        + xxxt * t_c13_y


                        + xx * c17_24_25_y
                        + xxxx * c18_y
                        + xxx * (c16_19_22_y + c20_21_23_y);

                zcomp +=  xrrr * grada_z_g * pi[g]  // c9 terms
                        + xrrt * grada_z_g * rrt_pi[g]
                        + xrrl * grada_z_g * rrl_pi[g]
                        + xrtt * grada_z_g * rtt_pi[g]
                        + xrlt * grada_z_g * rtl_pi[g]
                        + xrll * grada_z_g * rll_pi[g]
                        + xttt * grada_z_g * ttt_pi[g]
                        + xltt * grada_z_g * ttl_pi[g]
                        + xllt * grada_z_g * tll_pi[g]
                        + xlll * grada_z_g * lll_pi[g]

                        + xrr * piz[g] // c10 terms
                        + xrt * rt_piz[g]
                        + xrl * rl_piz[g]
                        + xll * ll_piz[g]
                        + xtt * tt_piz[g]
                        + xtl * tl_piz[g]

                        + xxrr *    c2 * grada_z_g   // c11 terms
                        + xxrt * rt_c2 * grada_z_g
                        + xxrl * rl_c2 * grada_z_g
                        + xxll * ll_c2 * grada_z_g
                        + xxtt * tt_c2 * grada_z_g
                        + xxlt * tl_c2 * grada_z_g

                        + xxr * (c12_c14_z + c15_z)
                        + xxl * (l_c12_c14_z + l_c15_z)
                        + xxt * (t_c12_c14_z + t_c15_z)

                        + xxxr * c13_z
                        + xxxl * l_c13_z
                        + xxxt * t_c13_z

                        + xx * c17_24_25_z
                        + xxxx * c18_z
                        + xxx * (c16_19_22_z + c20_21_23_z);

                G_val[nu_offset + g] = w * rho_0 * chi_val[nu_offset + g];

                G_gga_val[nu_offset + g] = w * (xcomp * chi_x_val[nu_offset + g] +
                                                ycomp * chi_y_val[nu_offset + g] +
                                                zcomp * chi_z_val[nu_offset + g]);

                // TODO implement Laplacian dependence

                //G_gga_val[nu_offset + g] += w * lap_0 * (chi_xx_val[nu_offset + g] +
                //                                         chi_yy_val[nu_offset + g] +
                //                                         chi_zz_val[nu_offset + g]);

                G_gga_x_val[nu_offset + g] = w * tau_0 * chi_x_val[nu_offset + g];
                G_gga_y_val[nu_offset + g] = w * tau_0 * chi_y_val[nu_offset + g];
                G_gga_z_val[nu_offset + g] = w * tau_0 * chi_z_val[nu_offset + g];
            }
        }
    }

    timer.stop("Lxc matrix G");

    // eq.(31), JCTC 2021, 17, 1512-1521

    timer.start("Lxc matrix matmul");

    // LDA and GGA contribution
    auto mat_Lxc = denblas::multABt(gtoValues, denblas::addAB(mat_G, mat_G_gga, 2.0));

    // tau contribution
    auto mat_Lxc_x = denblas::multABt(gtoValuesX, mat_G_gga_x);
    auto mat_Lxc_y = denblas::multABt(gtoValuesY, mat_G_gga_y);
    auto mat_Lxc_z = denblas::multABt(gtoValuesZ, mat_G_gga_z);

    mat_Lxc = denblas::addAB(mat_Lxc, mat_Lxc_x, 0.5);
    mat_Lxc = denblas::addAB(mat_Lxc, mat_Lxc_y, 0.5);
    mat_Lxc = denblas::addAB(mat_Lxc, mat_Lxc_z, 0.5);

    mat_Lxc.symmetrizeAndScale(0.5);

    timer.stop("Lxc matrix matmul");

    return mat_Lxc;
}

CDenseMatrix
CXCNewIntegrator::computeGtoValuesOnGridPoints(const CMolecule&       molecule,
                                               const CMolecularBasis& basis,
                                               const CMolecularGrid&  molecularGrid) const
{
    // std::string errnotpartitioned("CXCNewIntegrator.computeGtoValuesOnGridPoints: Cannot use unpartitioned molecular grid");

    // errors::assertMsgCritical(molecularGrid.isPartitioned(), errnotpartitioned);

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

    // destroy GTOs container

    delete gtovec;

    return allgtovalues;
}

std::vector<CDenseMatrix>
CXCNewIntegrator::computeGtoValuesAndDerivativesOnGridPoints(const CMolecule&       molecule,
                                                             const CMolecularBasis& basis,
                                                             const CMolecularGrid&  molecularGrid) const
{
    // std::string errnotpartitioned("CXCNewIntegrator.computeGtoValuesAndDerivativesOnGridPoints: Cannot use unpartitioned molecular grid");

    // errors::assertMsgCritical(molecularGrid.isPartitioned(), errnotpartitioned);

    auto nthreads = omp_get_max_threads();

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // memory blocks for GTOs on grid points

    CDenseMatrix allgtovalues(naos, molecularGrid.getNumberOfGridPoints());

    CDenseMatrix allgtoderivx(naos, molecularGrid.getNumberOfGridPoints());

    CDenseMatrix allgtoderivy(naos, molecularGrid.getNumberOfGridPoints());

    CDenseMatrix allgtoderivz(naos, molecularGrid.getNumberOfGridPoints());

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

        gtoeval::preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 1, _screeningThresholdForGTOValues, boxdim);  // 1st order GTO derivative

        // GTO values on grid points

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtovec, xcoords, ycoords, zcoords, gridblockpos,

                                             grid_batch_offset, grid_batch_size, skip_cgto_ids);
        }

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

        for (int32_t i = 0; i < aocount; i++)
        {
            auto aoidx = aoinds[i];

            std::memcpy(allgtovalues.row(aoidx) + gridblockpos, gaos.data(aoidx), npoints * sizeof(double));

            std::memcpy(allgtoderivx.row(aoidx) + gridblockpos, gaox.data(aoidx), npoints * sizeof(double));

            std::memcpy(allgtoderivy.row(aoidx) + gridblockpos, gaoy.data(aoidx), npoints * sizeof(double));

            std::memcpy(allgtoderivz.row(aoidx) + gridblockpos, gaoz.data(aoidx), npoints * sizeof(double));
        }
    }

    // destroy GTOs container

    delete gtovec;

    return std::vector<CDenseMatrix>({allgtovalues, allgtoderivx, allgtoderivy, allgtoderivz});
}

std::vector<CDenseMatrix>
CXCNewIntegrator::computeGtoValuesAndDerivativesOnGridPoints(const CMolecule&       molecule,
                                                             const CMolecularBasis& basis,
                                                             const int32_t          npoints,
                                                             const double*          xcoords,
                                                             const double*          ycoords,
                                                             const double*          zcoords) const
{
    auto nthreads = omp_get_max_threads();

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // memory blocks for GTOs on grid points

    CDenseMatrix allgtovalues(naos, npoints);

    CDenseMatrix allgtoderivx(naos, npoints);
    CDenseMatrix allgtoderivy(naos, npoints);
    CDenseMatrix allgtoderivz(naos, npoints);

    CMemBlock2D<double> gaos(npoints, naos);

    CMemBlock2D<double> gaox(npoints, naos);
    CMemBlock2D<double> gaoy(npoints, naos);
    CMemBlock2D<double> gaoz(npoints, naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos); // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

    {
        auto gridblockpos = 0;

        // dimension of grid box

        auto boxdim = gtoeval::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // pre-screening of GTOs

        gtoeval::preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 1, _screeningThresholdForGTOValues, boxdim);  // 1st order GTO derivative

        // GTO values on grid points

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtovec, xcoords, ycoords, zcoords, gridblockpos,
                                             grid_batch_offset, grid_batch_size, skip_cgto_ids);
        }

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

        for (int32_t i = 0; i < aocount; i++)
        {
            auto aoidx = aoinds[i];

            std::memcpy(allgtovalues.row(aoidx) + gridblockpos, gaos.data(aoidx), npoints * sizeof(double));

            std::memcpy(allgtoderivx.row(aoidx) + gridblockpos, gaox.data(aoidx), npoints * sizeof(double));
            std::memcpy(allgtoderivy.row(aoidx) + gridblockpos, gaoy.data(aoidx), npoints * sizeof(double));
            std::memcpy(allgtoderivz.row(aoidx) + gridblockpos, gaoz.data(aoidx), npoints * sizeof(double));
        }
    }

    // destroy GTOs container

    delete gtovec;

    return std::vector<CDenseMatrix>({allgtovalues, allgtoderivx, allgtoderivy, allgtoderivz});
}
