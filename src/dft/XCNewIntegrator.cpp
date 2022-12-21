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
#include "NewFunctionalParser.hpp"
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

CAOKohnShamMatrix
CXCNewIntegrator::integrateVxcFock(const CMolecule&        molecule,
                                   const CMolecularBasis&  basis,
                                   const CAODensityMatrix& densityMatrix,
                                   CMolecularGrid&         molecularGrid,
                                   const std::string&      xcFuncLabel) const
{
    molecularGrid.partitionGridPoints();

    molecularGrid.distributeCountsAndDisplacements(_locRank, _locNodes, _locComm);

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
        std::string erropenshell("XCNewIntegrator.integrateVxcFock: Not implemented for open-shell meta-GGA");

        errors::assertMsgCritical(fstr::upcase(flag) == "CLOSEDSHELL", erropenshell);

        return _integrateVxcFockForMGGA(molecule, basis, densityMatrix, molecularGrid, newfvxc, flag);
    }

    std::string errxcfuntype("XCNewIntegrator.integrateVxcFock: Only implemented for LDA/GGA/MGGA");

    errors::assertMsgCritical(false, errxcfuntype);

    return CAOKohnShamMatrix();
}

void
CXCNewIntegrator::integrateFxcFock(CAOFockMatrix&          aoFockMatrix,
                                   const CMolecule&        molecule,
                                   const CMolecularBasis&  basis,
                                   const CAODensityMatrix& rwDensityMatrix,
                                   const CAODensityMatrix& gsDensityMatrix,
                                   CMolecularGrid&         molecularGrid,
                                   const std::string&      xcFuncLabel) const
{
    molecularGrid.partitionGridPoints();

    molecularGrid.distributeCountsAndDisplacements(_locRank, _locNodes, _locComm);

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
            std::string errxcfuntype("XCNewIntegrator.integrateFxcFock: Only implemented for LDA/GGA/MGGA");

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
                                   CMolecularGrid&         molecularGrid,
                                   const std::string&      xcFuncLabel,
                                   const std::string&      quadMode) const
{
    molecularGrid.partitionGridPoints();

    molecularGrid.distributeCountsAndDisplacements(_locRank, _locNodes, _locComm);

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

void
CXCNewIntegrator::integrateLxcFock(CAOFockMatrix&          aoFockMatrix,
                                   const CMolecule&        molecule,
                                   const CMolecularBasis&  basis,
                                   const CAODensityMatrix& rwDensityMatrix,
                                   const CAODensityMatrix& rw2DensityMatrix,
                                   const CAODensityMatrix& rw3DensityMatrix,
                                   const CAODensityMatrix& gsDensityMatrix,
                                   CMolecularGrid&         molecularGrid,
                                   const std::string&      xcFuncLabel,
                                   const std::string&      cubeMode) const
{
    molecularGrid.partitionGridPoints();

    molecularGrid.distributeCountsAndDisplacements(_locRank, _locNodes, _locComm);

    auto fvxc = newvxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (rwDensityMatrix.isClosedShell() && rw2DensityMatrix.isClosedShell() && rw3DensityMatrix.isClosedShell() && gsDensityMatrix.isClosedShell())
    {
        if (xcfuntype == xcfun::lda)
        {
            _integrateLxcFockForLDA(aoFockMatrix, molecule, basis, rwDensityMatrix, rw2DensityMatrix,rw3DensityMatrix, gsDensityMatrix,

                                    molecularGrid, fvxc, cubeMode);
        }
        else if (xcfuntype == xcfun::gga)
        {
            _integrateLxcFockForGGA(aoFockMatrix, molecule, basis, rwDensityMatrix, rw2DensityMatrix,rw3DensityMatrix, gsDensityMatrix,

                                     molecularGrid, fvxc, cubeMode);
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
CXCNewIntegrator::integrateVxcPDFT(CAOKohnShamMatrix&      aoFockMatrix,
                                   CDense4DTensor&         moTwoBodyGradient,
                                   const CMolecule&        molecule,
                                   const CMolecularBasis&  basis,
                                   const CAODensityMatrix& DensityMatrix,
                                   const CDense4DTensor&   TwoBodyDensityMatrix,
                                   const CDenseMatrix&     ActiveMOs,
                                   CMolecularGrid&         molecularGrid,
                                   const std::string&      xcFuncLabel) const
{
    molecularGrid.partitionGridPoints();

    molecularGrid.distributeCountsAndDisplacements(_locRank, _locNodes, _locComm);

    auto fvxc = newvxcfuncs::getPairDensityExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (xcfuntype == "PLDA")
    {
        _integrateVxcPDFTForLDA(aoFockMatrix, moTwoBodyGradient, molecule, basis,
                                DensityMatrix, TwoBodyDensityMatrix, ActiveMOs, molecularGrid, fvxc);
    }
    else if (xcfuntype == "PGGA")
    {
        // _integrateVxcPDFTForGGA(aoFockMatrix, moTwoBodyGradient, molecule, basis,
        //                         DensityMatrix, TwoBodyDensityMatrix, ActiveMOs, molecularGrid, fvxc);
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

        gtoeval::preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 2, _screeningThresholdForGTOValues, boxdim);  // 2nd order GTO derivative

        timer.stop("GTO pre-screening");

        // GTO values on grid points

        timer.start("OMP GTO evaluation");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            omptimers[thread_id].start("gtoeval");

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForMetaGGA(gaos, gaox, gaoy, gaoz, gaoxx, gaoxy, gaoxz, gaoyy, gaoyz, gaozz,
                                                 gtovec, xcoords, ycoords, zcoords, gridblockpos,
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
        CDenseMatrix mat_chi_yy(aocount, npoints);
        CDenseMatrix mat_chi_zz(aocount, npoints);

        for (int32_t i = 0; i < aocount; i++)
        {
            std::memcpy(mat_chi.row(i), gaos.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_x.row(i), gaox.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_y.row(i), gaoy.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_z.row(i), gaoz.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_xx.row(i), gaoxx.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_yy.row(i), gaoyy.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_zz.row(i), gaozz.data(aoinds[i]), npoints * sizeof(double));
        }

        timer.stop("GTO screening");

        if (aocount == 0) continue;

        // generate sub density matrix and density grid

        timer.start("Density matrix slicing");

        auto sub_dens_mat = submat::getSubDensityMatrix(densityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        // TODO open-shell meta-GGA
        dengridgen::generateDensityForMGGA(rho, rhograd, sigma, lapl, tau, npoints,
                                           mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                           mat_chi_xx, mat_chi_yy, mat_chi_zz,
                                           sub_dens_mat, timer);

        timer.stop("Density matrix slicing");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_mgga(npoints, rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        // TODO screening
        // gridscreen::screenVxcFockForMGGA(rho, sigma, exc, vrho, vsigma, vlapl, vtau,
        //                                  npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // compute partial contribution to Vxc matrix and distribute partial
        // Vxc to full Kohn-Sham matrix

        // TODO open-shell meta-GGA
        auto partial_mat_Vxc = _integratePartialVxcFockForMGGA(
            npoints, local_weights, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
            mat_chi_xx, mat_chi_yy, mat_chi_zz, rhograd, vrho, vsigma, vlapl, vtau, timer);

        timer.start("Vxc matrix dist.");

        submat::distributeSubMatrixToKohnSham(mat_Vxc, partial_mat_Vxc, aoinds, aocount, naos);

        timer.stop("Vxc matrix dist.");

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

            gtoeval::computeGtosValuesForMetaGGA(gaos, gaox, gaoy, gaoz, gaoxx, gaoxy, gaoxz, gaoyy, gaoyz, gaozz,
                                                 gtovec, xcoords, ycoords, zcoords, gridblockpos,
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
        CDenseMatrix mat_chi_yy(aocount, npoints);
        CDenseMatrix mat_chi_zz(aocount, npoints);

        for (int32_t i = 0; i < aocount; i++)
        {
            std::memcpy(mat_chi.row(i), gaos.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_x.row(i), gaox.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_y.row(i), gaoy.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_z.row(i), gaoz.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_xx.row(i), gaoxx.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_yy.row(i), gaoyy.data(aoinds[i]), npoints * sizeof(double));
            std::memcpy(mat_chi_zz.row(i), gaozz.data(aoinds[i]), npoints * sizeof(double));
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
                                           mat_chi_xx, mat_chi_yy, mat_chi_zz,
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
                                                mat_chi_xx, mat_chi_yy, mat_chi_zz,
                                                sub_dens_mat, timer);

            // compute partial contribution to Fxc matrix


            auto partial_mat_Fxc = _integratePartialFxcFockForMGGA(npoints, local_weights, 
                                                                  mat_chi, mat_chi_x, 
                                                                  mat_chi_y, mat_chi_z,
                                                                  mat_chi_xx, mat_chi_yy, 
                                                                  mat_chi_zz, 
                                                                  rhow, rhograd, rhowgrad, tauw, laplw,
                                                                  vrho, vsigma, vlapl, vtau, v2rho2,
                                                                  v2lapl2, v2tau2, v2rholapl, v2rhotau,
                                                                  v2lapltau, v2rhosigma, v2sigmalapl, v2sigmatau,v2sigma2, 
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
CXCNewIntegrator::_integrateVxcPDFTForGGA(CAOKohnShamMatrix&      aoFockMatrix,
                                          CDense4DTensor&         moTwoBodyGradient,
                                          const CMolecule&        molecule,
                                          const CMolecularBasis&  basis,
                                          const CAODensityMatrix& DensityMatrix,
                                          const CDense4DTensor&   TwoBodyDensityMatrix,
                                          const CDenseMatrix&     ActiveMOs,
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

        // TODO (MGD) implement and use own functionals
        xcFunctional.compute_exc_vxc_for_gga(npoints, rho, sigma, exc, vrho, vsigma);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenVxcFockForGGA(rho, sigma, exc, vrho, vsigma, npoints, _screeningThresholdForDensityValues);

        // TODO (MGD) screening

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
            auto rho_total = rho[2 * g + 0] + rho[2 * g + 1];

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

    auto G_gga_a_val = mat_G_a_gga.values();
    auto G_gga_b_val = mat_G_b_gga.values();

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

                G_gga_a_val[nu_offset + g] = weights[g] * (vxa * chi_x_val[nu_offset + g] +
                                                           vya * chi_y_val[nu_offset + g] +
                                                           vza * chi_z_val[nu_offset + g]);
                G_gga_b_val[nu_offset + g] = weights[g] * (vxb * chi_x_val[nu_offset + g] +
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
                                                  const CDenseMatrix&    gtoValuesXX,
                                                  const CDenseMatrix&    gtoValuesYY,
                                                  const CDenseMatrix&    gtoValuesZZ,
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

    // GGA and laplacian contribution
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

    auto chi_xx_val = gtoValuesXX.values();
    auto chi_yy_val = gtoValuesYY.values();
    auto chi_zz_val = gtoValuesZZ.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

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

                // laplacian contribution (will be scaled by 2 later)
                G_gga_val[nu_offset + g] += weights[g] * vlapl[2 * g + 0] * (chi_xx_val[nu_offset + g] +
                                                                             chi_yy_val[nu_offset + g] +
                                                                             chi_zz_val[nu_offset + g]);

                // tau contribution (will be scaled by 0.5 later)
                G_gga_x_val[nu_offset + g] += weights[g] * vtau[2 * g + 0] * chi_x_val[nu_offset + g];
                G_gga_y_val[nu_offset + g] += weights[g] * vtau[2 * g + 0] * chi_y_val[nu_offset + g];
                G_gga_z_val[nu_offset + g] += weights[g] * vtau[2 * g + 0] * chi_z_val[nu_offset + g];
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

    // LDA, GGA and laplacian contribution
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
CXCNewIntegrator::_integratePartialFxcFockForMGGA(const int32_t           npoints, 
                                                  const double*           weights, 
                                                  const CDenseMatrix&     gtoValues,
                                                  const CDenseMatrix&     gtoValuesX, 
                                                  const CDenseMatrix&     gtoValuesY, 
                                                  const CDenseMatrix&     gtoValuesZ,
                                                  const CDenseMatrix&     gtoValuesXX, 
                                                  const CDenseMatrix&     gtoValuesYY, 
                                                  const CDenseMatrix&     gtoValuesZZ, 
                                                  const double*           rhow, 
                                                  const double*           rhograd, 
                                                  const double*           rhowgrad, 
                                                  const double*           tauw, 
                                                  const double*           laplw,
                                                  const double*           vrho, 
                                                  const double*           vsigma, 
                                                  const double*           vlapl, 
                                                  const double*           vtau, 
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
                                                  CMultiTimer&            timer) const                                                                                              
{
   
    // eq.(30), JCTC 2021, 17, 1512-1521

    timer.start("Vxc matrix G");

    auto naos = gtoValues.getNumberOfRows();

    // LDA contribution
    CDenseMatrix mat_G(naos, npoints);

    // GGA and laplacian contribution
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

    auto chi_xx_val = gtoValuesXX.values();
    auto chi_yy_val = gtoValuesYY.values();
    auto chi_zz_val = gtoValuesZZ.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto nthreads = omp_get_max_threads();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

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
                // second-order gga
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
                double ll = v2lapl2_aa + v2lapl2_ab;
                double lx = 2.0 * v2sigmalapl_ca + 2.0 * v2sigmalapl_ba + 2.0 * v2sigmalapl_aa;
                double lr = v2rholapl_aa + v2rholapl_ba;
                double lt = v2lapltau_aa + v2lapltau_ab;

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


                // laplacian contribution (will be scaled by 2 later)
                double lap_0 =    lr * rwa 
                                + lx * contract
                                + lt * tauwa 
                                + ll * laplwa; 

                G_gga_val[nu_offset + g] += w * lap_0 * (chi_xx_val[nu_offset + g] +
                                                        chi_yy_val[nu_offset + g] +
                                                        chi_zz_val[nu_offset + g]);

                // tau contribution (will be scaled by 0.5 later)
                double tau_0 =    tr * rwa 
                                + tx * contract
                                + tt * tauwa 
                                + tl * laplwa;   

                G_gga_x_val[nu_offset + g] += w * tau_0 * chi_x_val[nu_offset + g];
                G_gga_y_val[nu_offset + g] += w * tau_0 * chi_y_val[nu_offset + g];
                G_gga_z_val[nu_offset + g] += w * tau_0 * chi_z_val[nu_offset + g];
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

    // LDA, GGA and laplacian contribution
    auto mat_Fxc = denblas::multABt(gtoValues, denblas::addAB(mat_G, mat_G_gga, 2.0));

    // tau contribution
    auto mat_Fxc_x = denblas::multABt(gtoValuesX, mat_G_gga_x);
    auto mat_Fxc_y = denblas::multABt(gtoValuesY, mat_G_gga_y);
    auto mat_Fxc_z = denblas::multABt(gtoValuesZ, mat_G_gga_z);

    mat_Fxc = denblas::addAB(mat_Fxc, mat_Fxc_x, 0.5);
    mat_Fxc = denblas::addAB(mat_Fxc, mat_Fxc_y, 0.5);
    mat_Fxc = denblas::addAB(mat_Fxc, mat_Fxc_z, 0.5);

    mat_Fxc.symmetrizeAndScale(0.5);

    timer.stop("Vxc matrix matmul");

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
CXCNewIntegrator::_integratePartialKxcFockForGGA(const int32_t              npoints,
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
                                                 const CDensityGridQuad&    rwDensityGridQuad,
                                                 const CDensityGrid&        rw2DensityGrid,
                                                 const int32_t              iFock,
                                                 CMultiTimer&               timer) const
{
    timer.start("Kxc matrix prep.");

    // GTO values on grid points

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

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

                double c3 =  grada_x_g * grada_x_g *  pixx[g]
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
CXCNewIntegrator::computeGtoValuesOnGridPoints(const CMolecule&       molecule,
                                               const CMolecularBasis& basis,
                                               CMolecularGrid&        molecularGrid) const
{
    molecularGrid.partitionGridPoints();

    molecularGrid.distributeCountsAndDisplacements(_locRank, _locNodes, _locComm);

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
                                                             CMolecularGrid&        molecularGrid) const
{
    molecularGrid.partitionGridPoints();

    molecularGrid.distributeCountsAndDisplacements(_locRank, _locNodes, _locComm);

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

void
CXCNewIntegrator::computeExcVxcForLDA(const std::string& xcFuncLabel,
                                      const int32_t      npoints,
                                      const double*      rho,
                                      double*            exc,
                                      double*            vrho) const
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    // form density grid

    CDensityGrid dgrid(npoints, 1, xcfuntype, dengrid::ab);

    auto rhoa = dgrid.alphaDensity(0);

    auto rhob = dgrid.betaDensity(0);

    for (int32_t g = 0; g < npoints; g++)
    {
        rhoa[g] = rho[2 * g + 0];

        rhob[g] = rho[2 * g + 1];
    }

    // compute functional derivative

    CXCGradientGrid vxcgrid(npoints, dengrid::ab, xcfuntype);

    fvxc.compute(vxcgrid, dgrid);

    auto efunc = vxcgrid.xcFunctionalValues();

    auto grhoa = vxcgrid.xcGradientValues(xcvars::rhoa);

    auto grhob = vxcgrid.xcGradientValues(xcvars::rhob);

    for (int32_t g = 0; g < npoints; g++)
    {
        exc[g] = efunc[g] / (rho[2 * g + 0] + rho[2 * g + 1]);

        vrho[2 * g + 0] = grhoa[g];

        vrho[2 * g + 1] = grhob[g];
    }
}

void
CXCNewIntegrator::computeExcVxcForGGA(const std::string& xcFuncLabel,
                                      const int32_t      npoints,
                                      const double*      rho,
                                      const double*      sigma,
                                      double*            exc,
                                      double*            vrho,
                                      double*            vsigma) const
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    // form density grid

    CDensityGrid dgrid(npoints, 1, xcfuntype, dengrid::ab);

    auto rhoa = dgrid.alphaDensity(0);

    auto rhob = dgrid.betaDensity(0);

    auto grada = dgrid.alphaDensityGradient(0);

    auto gradb = dgrid.betaDensityGradient(0);

    auto gradab = dgrid.mixedDensityGradient(0);

    for (int32_t g = 0; g < npoints; g++)
    {
        rhoa[g] = rho[2 * g + 0];

        rhob[g] = rho[2 * g + 1];

        grada[g] = std::sqrt(sigma[3 * g + 0]);

        gradab[g] = sigma[3 * g + 1];

        gradb[g] = std::sqrt(sigma[3 * g + 2]);
    }

    // compute functional derivative

    CXCGradientGrid vxcgrid(npoints, dengrid::ab, xcfuntype);

    fvxc.compute(vxcgrid, dgrid);

    auto efunc = vxcgrid.xcFunctionalValues();

    auto grhoa = vxcgrid.xcGradientValues(xcvars::rhoa);

    auto grhob = vxcgrid.xcGradientValues(xcvars::rhob);

    auto ggrada = vxcgrid.xcGradientValues(xcvars::grada);

    auto ggradab = vxcgrid.xcGradientValues(xcvars::gradab);

    auto ggradb = vxcgrid.xcGradientValues(xcvars::gradb);

    auto ngrada = grada;

    auto ngradb = gradb;

    for (int32_t g = 0; g < npoints; g++)
    {
        exc[g] = efunc[g] / (rho[2 * g + 0] + rho[2 * g + 1]);

        vrho[2 * g + 0] = grhoa[g];

        vrho[2 * g + 1] = grhob[g];

        vsigma[3 * g + 0] = 0.5 * ggrada[g] / ngrada[g];

        vsigma[3 * g + 1] = ggradab[g];

        vsigma[3 * g + 2] = 0.5 * ggradb[g] / ngradb[g];
    }
}

void
CXCNewIntegrator::computeFxcForLDA(const std::string& xcFuncLabel,
                                   const int32_t      npoints,
                                   const double*      rho,
                                   double*            v2rho2) const
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    // form density grid

    CDensityGrid dgrid(npoints, 1, xcfuntype, dengrid::ab);

    auto rhoa = dgrid.alphaDensity(0);

    auto rhob = dgrid.betaDensity(0);

    for (int32_t g = 0; g < npoints; g++)
    {
        rhoa[g] = rho[2 * g + 0];

        rhob[g] = rho[2 * g + 1];
    }

    // compute functional derivative

    CXCHessianGrid vxc2grid(npoints, dengrid::ab, xcfuntype);

    fvxc.compute(vxc2grid, dgrid);

    auto grho_aa = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::rhoa);

    auto grho_ab = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::rhob);

    auto grho_bb = vxc2grid.xcHessianValues(xcvars::rhob, xcvars::rhob);

    for (int32_t g = 0; g < npoints; g++)
    {
        v2rho2[3 * g + 0] = grho_aa[g];

        v2rho2[3 * g + 1] = grho_ab[g];

        v2rho2[3 * g + 2] = grho_bb[g];
    }
}

void
CXCNewIntegrator::computeFxcForGGA(const std::string& xcFuncLabel,
                                   const int32_t      npoints,
                                   const double*      rho,
                                   const double*      sigma,
                                   double*            v2rho2,
                                   double*            v2rhosigma,
                                   double*            v2sigma2) const
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    // form density grid

    CDensityGrid dgrid(npoints, 1, xcfuntype, dengrid::ab);

    auto rhoa = dgrid.alphaDensity(0);

    auto rhob = dgrid.betaDensity(0);

    auto grada = dgrid.alphaDensityGradient(0);

    auto gradb = dgrid.betaDensityGradient(0);

    auto gradab = dgrid.mixedDensityGradient(0);

    for (int32_t g = 0; g < npoints; g++)
    {
        rhoa[g] = rho[2 * g + 0];

        rhob[g] = rho[2 * g + 1];

        grada[g] = std::sqrt(sigma[3 * g + 0]);

        gradab[g] = sigma[3 * g + 1];

        gradb[g] = std::sqrt(sigma[3 * g + 2]);
    }

    // compute functional derivative

    CXCGradientGrid vxcgrid(npoints, dengrid::ab, xcfuntype);

    CXCHessianGrid vxc2grid(npoints, dengrid::ab, xcfuntype);

    fvxc.compute(vxcgrid, dgrid);

    fvxc.compute(vxc2grid, dgrid);

    auto ngrad_a = grada;

    auto ngrad_b = gradb;

    auto ggrad_a = vxcgrid.xcGradientValues(xcvars::grada);

    auto ggrad_b = vxcgrid.xcGradientValues(xcvars::gradb);

    auto grho_aa = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::rhoa);

    auto grho_ab = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::rhob);

    auto grho_bb = vxc2grid.xcHessianValues(xcvars::rhob, xcvars::rhob);

    auto gmix_aa = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::grada);

    auto gmix_ab = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::gradb);

    auto gmix_ba = vxc2grid.xcHessianValues(xcvars::rhob, xcvars::grada);

    auto gmix_bb = vxc2grid.xcHessianValues(xcvars::rhob, xcvars::gradb);

    auto gmix_ac = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::gradab);

    auto gmix_bc = vxc2grid.xcHessianValues(xcvars::rhob, xcvars::gradab);

    auto ggrad_aa = vxc2grid.xcHessianValues(xcvars::grada, xcvars::grada);

    auto ggrad_ab = vxc2grid.xcHessianValues(xcvars::grada, xcvars::gradb);

    auto ggrad_bb = vxc2grid.xcHessianValues(xcvars::gradb, xcvars::gradb);

    auto ggrad_ac = vxc2grid.xcHessianValues(xcvars::grada, xcvars::gradab);

    auto ggrad_bc = vxc2grid.xcHessianValues(xcvars::gradb, xcvars::gradab);

    auto ggrad_cc = vxc2grid.xcHessianValues(xcvars::gradab, xcvars::gradab);

    for (int32_t g = 0; g < npoints; g++)
    {
        v2rho2[3 * g + 0] = grho_aa[g];

        v2rho2[3 * g + 1] = grho_ab[g];

        v2rho2[3 * g + 2] = grho_bb[g];

        v2rhosigma[6 * g + 0] = 0.5 * gmix_aa[g] / ngrad_a[g];

        v2rhosigma[6 * g + 1] = gmix_ac[g];

        v2rhosigma[6 * g + 2] = 0.5 * gmix_ab[g] / ngrad_b[g];

        v2rhosigma[6 * g + 3] = 0.5 * gmix_ba[g] / ngrad_a[g];

        v2rhosigma[6 * g + 4] = gmix_bc[g];

        v2rhosigma[6 * g + 5] = 0.5 * gmix_bb[g] / ngrad_b[g];

        v2sigma2[6 * g + 0] = 0.25 * (ggrad_aa[g] - ggrad_a[g] / ngrad_a[g]) / (ngrad_a[g] * ngrad_a[g]);

        v2sigma2[6 * g + 1] = 0.5 * ggrad_ac[g] / ngrad_a[g];

        v2sigma2[6 * g + 2] = 0.25 * ggrad_ab[g] / (ngrad_a[g] * ngrad_b[g]);

        v2sigma2[6 * g + 3] = ggrad_cc[g];

        v2sigma2[6 * g + 4] = 0.5 * ggrad_bc[g] / ngrad_b[g];

        v2sigma2[6 * g + 5] = 0.25 * (ggrad_bb[g] - ggrad_b[g] / ngrad_b[g]) / (ngrad_b[g] * ngrad_b[g]);
    }
}
