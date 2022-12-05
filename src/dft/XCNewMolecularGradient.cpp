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

#include "XCNewMolecularGradient.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "AODensityMatrix.hpp"
#include "AngularMomentum.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DensityGridGenerator.hpp"
#include "DensityGridType.hpp"
#include "FunctionalParser.hpp"
#include "GridScreener.hpp"
#include "GtoEvaluator.hpp"
#include "SubMatrix.hpp"
#include "XCFuncType.hpp"
#include "XCVarsType.hpp"

CXCNewMolecularGradient::CXCNewMolecularGradient(MPI_Comm comm)

    : _screeningThresholdForGTOValues(1.0e-12)

    , _screeningThresholdForDensityValues(1.0e-13)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;
}

CXCNewMolecularGradient::~CXCNewMolecularGradient()
{
}

CDenseMatrix
CXCNewMolecularGradient::integrateVxcGradient(const CMolecule&        molecule,
                                              const CMolecularBasis&  basis,
                                              const CAODensityMatrix& gsDensityMatrix,
                                              CMolecularGrid&         molecularGrid,
                                              const std::string&      xcFuncLabel) const
{
    return integrateVxcGradient(molecule, basis, gsDensityMatrix, gsDensityMatrix, molecularGrid, xcFuncLabel);
}

CDenseMatrix
CXCNewMolecularGradient::integrateVxcGradient(const CMolecule&        molecule,
                                              const CMolecularBasis&  basis,
                                              const CAODensityMatrix& rwDensityMatrix,
                                              const CAODensityMatrix& gsDensityMatrix,
                                              CMolecularGrid&         molecularGrid,
                                              const std::string&      xcFuncLabel) const
{
    molecularGrid.partitionGridPoints();

    molecularGrid.distributeCountsAndDisplacements(_locRank, _locNodes, _locComm);

    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (rwDensityMatrix.isClosedShell() && gsDensityMatrix.isClosedShell())
    {
        if (xcfuntype == xcfun::lda)
        {
            return _integrateVxcGradientForLDA(molecule, basis, rwDensityMatrix, gsDensityMatrix, molecularGrid, fvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
            return _integrateVxcGradientForGGA(molecule, basis, rwDensityMatrix, gsDensityMatrix, molecularGrid, fvxc);
        }
        else
        {
            std::string errxcfuntype("XCNewMolecularGradient.integrateVxcGradient: Only implemented for LDA/GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCNewMolecularGradient.integrateVxcGradient: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }

    return CDenseMatrix();
}

CDenseMatrix
CXCNewMolecularGradient::integrateFxcGradient(const CMolecule&        molecule,
                                              const CMolecularBasis&  basis,
                                              const CAODensityMatrix& rwDensityMatrixOne,
                                              const CAODensityMatrix& rwDensityMatrixTwo,
                                              const CAODensityMatrix& gsDensityMatrix,
                                              CMolecularGrid&         molecularGrid,
                                              const std::string&      xcFuncLabel) const
{
    molecularGrid.partitionGridPoints();

    molecularGrid.distributeCountsAndDisplacements(_locRank, _locNodes, _locComm);

    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (rwDensityMatrixOne.isClosedShell() && rwDensityMatrixTwo.isClosedShell() && gsDensityMatrix.isClosedShell())
    {
        if (xcfuntype == xcfun::lda)
        {
            return _integrateFxcGradientForLDA(molecule, basis, rwDensityMatrixOne, rwDensityMatrixTwo,

                                               gsDensityMatrix, molecularGrid, fvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
            return _integrateFxcGradientForGGA(molecule, basis, rwDensityMatrixOne, rwDensityMatrixTwo,

                                               gsDensityMatrix, molecularGrid, fvxc);
        }
        else
        {
            std::string errxcfuntype("XCNewMolecularGradient.integrateFxcGradient: Only implemented for LDA/GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCNewMolecularGradient.integrateFxcGradient: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }

    return CDenseMatrix();
}

CDenseMatrix
CXCNewMolecularGradient::integrateKxcGradient(const CMolecule&        molecule,
                                              const CMolecularBasis&  basis,
                                              const CAODensityMatrix& rwDensityMatrixOne,
                                              const CAODensityMatrix& rwDensityMatrixTwo,
                                              const CAODensityMatrix& gsDensityMatrix,
                                              CMolecularGrid&         molecularGrid,
                                              const std::string&      xcFuncLabel) const
{
    molecularGrid.partitionGridPoints();

    molecularGrid.distributeCountsAndDisplacements(_locRank, _locNodes, _locComm);

    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (rwDensityMatrixOne.isClosedShell() && rwDensityMatrixTwo.isClosedShell() && gsDensityMatrix.isClosedShell())
    {
        if (xcfuntype == xcfun::lda)
        {
            return _integrateKxcGradientForLDA(molecule, basis, rwDensityMatrixOne, rwDensityMatrixTwo,

                                               gsDensityMatrix, molecularGrid, fvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
            return _integrateKxcGradientForGGA(molecule, basis, rwDensityMatrixOne, rwDensityMatrixTwo,

                                               gsDensityMatrix, molecularGrid, fvxc);
        }
        else
        {
            std::string errxcfuntype("XCNewMolecularGradient.integrateKxcGradient: Only implemented for LDA/GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCNewMolecularGradient.integrateKxcGradient: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }

    return CDenseMatrix();
}

CDenseMatrix
CXCNewMolecularGradient::_integrateVxcGradientForLDA(const CMolecule&        molecule,
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

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // AO-to-atom mapping

    std::vector<int32_t> ao_to_atom_ids(naos);

    _computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.getNumberOfAtoms();

    CMemBlock2D<double> molgrad_threads(natoms * 3,  nthreads);

    // memory blocks for GTOs on grid points

    CMemBlock2D<double> gaos(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaox(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaoy(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaoz(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos);  // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

    // indices for keeping track of valid grid points

    std::vector<int32_t> screened_point_inds(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> screened_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto screened_weights = screened_weights_data.data();

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

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtovec, xcoords, ycoords, zcoords,

                                             gridblockpos, grid_batch_offset, grid_batch_size, skip_cgto_ids);
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

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        auto rw_sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto gsdengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, gs_sub_dens_mat, xcfuntype, timer);

        // screen density grid, weights and GTO matrix

        timer.start("Density screening");

        CDensityGrid screened_gsdengrid(gsdengrid);

        gridscreen::screenDensityGridForLDA(screened_point_inds, screened_gsdengrid, gsdengrid,

                                            _screeningThresholdForDensityValues);

        auto screened_npoints = screened_gsdengrid.getNumberOfGridPoints();

        gridscreen::screenWeights(screened_weights, gridblockpos, weights, screened_point_inds, screened_npoints);

        CDenseMatrix screened_mat_chi(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_x(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_y(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_z(mat_chi.getNumberOfRows(), screened_npoints);

        gridscreen::screenGtoMatrixForGGA(screened_mat_chi, screened_mat_chi_x, screened_mat_chi_y, screened_mat_chi_z,

                                          mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, screened_point_inds, screened_npoints);

        timer.stop("Density screening");

        if (screened_npoints == 0) continue;

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, screened_npoints);

        CDenseMatrix dengrady(natoms, screened_npoints);

        CDenseMatrix dengradz(natoms, screened_npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(rw_sub_dens_mat, screened_mat_chi);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto naos = screened_mat_chi.getNumberOfRows();

        auto F_val = mat_F.values();

        auto chi_x_val = screened_mat_chi_x.values();

        auto chi_y_val = screened_mat_chi_y.values();

        auto chi_z_val = screened_mat_chi_z.values();

        auto gdenx = dengradx.values();

        auto gdeny = dengrady.values();

        auto gdenz = dengradz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(screened_npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(screened_npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * screened_npoints;

                auto nu_offset = nu * screened_npoints;

                #pragma omp simd aligned(gdenx, gdeny, gdenz, \
                        F_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];

                    gdeny[atom_g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];

                    gdenz[atom_g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];
                }
            }
        }

        timer.stop("Density grad. grid rho");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        CXCGradientGrid vxcgrid(screened_npoints, screened_gsdengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxcgrid, screened_gsdengrid);

        timer.stop("XC functional eval.");

        // exchange-correlation functional derivative

        auto grhoa = vxcgrid.xcGradientValues(xcvars::rhoa);

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(screened_npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(screened_npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.data(thread_id);

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * screened_npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(screened_weights, grhoa, \
                        gdenx, gdeny, gdenz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac = screened_weights[g] * grhoa[g];

                    gatmx += prefac * gdenx[atom_g];

                    gatmy += prefac * gdeny[atom_g];

                    gatmz += prefac * gdenz[atom_g];
                }

                // factor of 2 from sum of alpha and beta contributions

                gatm[iatom * 3 + 0] += 2.0 * gatmx;

                gatm[iatom * 3 + 1] += 2.0 * gatmy;

                gatm[iatom * 3 + 2] += 2.0 * gatmz;
            }
        }

        timer.stop("Accumulate gradient");
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    std::cout << "Timing of new integrator" << std::endl;
    std::cout << "------------------------" << std::endl;
    std::cout << timer.getSummary() << std::endl;

    CDenseMatrix molgrad(natoms, 3);

    for (int32_t iatom = 0; iatom < natoms; iatom++)
    {
        for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
        {
            molgrad.row(iatom)[0] += molgrad_threads.data(thread_id)[iatom * 3 + 0];

            molgrad.row(iatom)[1] += molgrad_threads.data(thread_id)[iatom * 3 + 1];

            molgrad.row(iatom)[2] += molgrad_threads.data(thread_id)[iatom * 3 + 2];
        }
    }

    return molgrad;
}

CDenseMatrix
CXCNewMolecularGradient::_integrateVxcGradientForGGA(const CMolecule&        molecule,
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

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // AO-to-atom mapping

    std::vector<int32_t> ao_to_atom_ids(naos);

    _computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.getNumberOfAtoms();

    CMemBlock2D<double> molgrad_threads(natoms * 3,  nthreads);

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

    CMemBlock<int32_t> skip_cgto_ids(naos);  // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

    // indices for keeping track of valid grid points

    std::vector<int32_t> screened_point_inds(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> screened_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto screened_weights = screened_weights_data.data();

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

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForMetaGGA(gaos, gaox, gaoy, gaoz, gaoxx, gaoxy, gaoxz, gaoyy, gaoyz, gaozz,

                                                 gtovec, xcoords, ycoords, zcoords,

                                                 gridblockpos, grid_batch_offset, grid_batch_size, skip_cgto_ids);
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

        CDenseMatrix mat_chi_xy(aocount, npoints);

        CDenseMatrix mat_chi_xz(aocount, npoints);

        CDenseMatrix mat_chi_yy(aocount, npoints);

        CDenseMatrix mat_chi_yz(aocount, npoints);

        CDenseMatrix mat_chi_zz(aocount, npoints);

        for (int32_t i = 0; i < aocount; i++)
        {
            std::memcpy(mat_chi.row(i), gaos.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_x.row(i), gaox.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_y.row(i), gaoy.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_z.row(i), gaoz.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_xx.row(i), gaoxx.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_xy.row(i), gaoxy.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_xz.row(i), gaoxz.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_yy.row(i), gaoyy.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_yz.row(i), gaoyz.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_zz.row(i), gaozz.data(aoinds[i]), npoints * sizeof(double));
        }

        timer.stop("GTO screening");

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        auto rw_sub_dens_mat = submat::getSubDensityMatrix(rwDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto gsdengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                                               gs_sub_dens_mat, xcfuntype, timer);

        // screen density grid, weights and GTO matrix

        timer.start("Density screening");

        CDensityGrid screened_gsdengrid(gsdengrid);

        gridscreen::screenDensityGridForGGA(screened_point_inds, screened_gsdengrid, gsdengrid,

                                            _screeningThresholdForDensityValues);

        auto screened_npoints = screened_gsdengrid.getNumberOfGridPoints();

        gridscreen::screenWeights(screened_weights, gridblockpos, weights, screened_point_inds, screened_npoints);

        CDenseMatrix screened_mat_chi(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_x(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_y(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_z(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_xx(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_xy(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_xz(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_yy(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_yz(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_zz(mat_chi.getNumberOfRows(), screened_npoints);

        gridscreen::screenGtoMatrixForMetaGGA(screened_mat_chi, screened_mat_chi_x, screened_mat_chi_y, screened_mat_chi_z,

                                              screened_mat_chi_xx, screened_mat_chi_xy, screened_mat_chi_xz,

                                              screened_mat_chi_yy, screened_mat_chi_yz, screened_mat_chi_zz,

                                              mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                              mat_chi_xx, mat_chi_xy, mat_chi_xz, mat_chi_yy, mat_chi_yz, mat_chi_zz,

                                              screened_point_inds, screened_npoints);

        timer.stop("Density screening");

        if (screened_npoints == 0) continue;

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, screened_npoints);

        CDenseMatrix dengrady(natoms, screened_npoints);

        CDenseMatrix dengradz(natoms, screened_npoints);

        CDenseMatrix dengradxx(natoms, screened_npoints);

        CDenseMatrix dengradxy(natoms, screened_npoints);

        CDenseMatrix dengradxz(natoms, screened_npoints);

        CDenseMatrix dengradyx(natoms, screened_npoints);

        CDenseMatrix dengradyy(natoms, screened_npoints);

        CDenseMatrix dengradyz(natoms, screened_npoints);

        CDenseMatrix dengradzx(natoms, screened_npoints);

        CDenseMatrix dengradzy(natoms, screened_npoints);

        CDenseMatrix dengradzz(natoms, screened_npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(rw_sub_dens_mat, screened_mat_chi);

        auto mat_F_x = denblas::multAB(rw_sub_dens_mat, screened_mat_chi_x);

        auto mat_F_y = denblas::multAB(rw_sub_dens_mat, screened_mat_chi_y);

        auto mat_F_z = denblas::multAB(rw_sub_dens_mat, screened_mat_chi_z);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto naos = screened_mat_chi.getNumberOfRows();

        auto F_val = mat_F.values();

        auto F_x_val = mat_F_x.values();

        auto F_y_val = mat_F_y.values();

        auto F_z_val = mat_F_z.values();

        auto chi_x_val = screened_mat_chi_x.values();

        auto chi_y_val = screened_mat_chi_y.values();

        auto chi_z_val = screened_mat_chi_z.values();

        auto chi_xx_val = screened_mat_chi_xx.values();

        auto chi_xy_val = screened_mat_chi_xy.values();

        auto chi_xz_val = screened_mat_chi_xz.values();

        auto chi_yy_val = screened_mat_chi_yy.values();

        auto chi_yz_val = screened_mat_chi_yz.values();

        auto chi_zz_val = screened_mat_chi_zz.values();

        auto gdenx = dengradx.values();

        auto gdeny = dengrady.values();

        auto gdenz = dengradz.values();

        auto gdenxx = dengradxx.values();

        auto gdenxy = dengradxy.values();

        auto gdenxz = dengradxz.values();

        auto gdenyx = dengradyx.values();

        auto gdenyy = dengradyy.values();

        auto gdenyz = dengradyz.values();

        auto gdenzx = dengradzx.values();

        auto gdenzy = dengradzy.values();

        auto gdenzz = dengradzz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(screened_npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(screened_npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * screened_npoints;

                auto nu_offset = nu * screened_npoints;

                #pragma omp simd aligned(gdenx, gdeny, gdenz, \
                        gdenxx, gdenxy, gdenxz, gdenyx, gdenyy, gdenyz, gdenzx, gdenzy, gdenzz, \
                        F_val, F_x_val, F_y_val, F_z_val, chi_x_val, chi_y_val, chi_z_val, \
                        chi_xx_val, chi_xy_val, chi_xz_val, chi_yy_val, chi_yz_val, chi_zz_val : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];

                    gdeny[atom_g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];

                    gdenz[atom_g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];

                    gdenxx[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xx_val[nu_g]);

                    gdenxy[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_xy_val[nu_g]);

                    gdenxz[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_xz_val[nu_g]);

                    gdenyx[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xy_val[nu_g]);

                    gdenyy[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_yy_val[nu_g]);

                    gdenyz[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_yz_val[nu_g]);

                    gdenzx[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xz_val[nu_g]);

                    gdenzy[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_yz_val[nu_g]);

                    gdenzz[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_zz_val[nu_g]);
                }
            }
        }

        timer.stop("Density grad. grid rho");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        CXCGradientGrid vxcgrid(screened_npoints, screened_gsdengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxcgrid, screened_gsdengrid);

        timer.stop("XC functional eval.");

        // pointers to exchange-correlation functional derrivatives

        auto grhoa = vxcgrid.xcGradientValues(xcvars::rhoa);

        auto ggrada = vxcgrid.xcGradientValues(xcvars::grada);

        auto ggradab = vxcgrid.xcGradientValues(xcvars::gradab);

        // pointers to density gradient norms

        auto ngrada = screened_gsdengrid.alphaDensityGradient(0);

        auto gradax = screened_gsdengrid.alphaDensityGradientX(0);

        auto graday = screened_gsdengrid.alphaDensityGradientY(0);

        auto gradaz = screened_gsdengrid.alphaDensityGradientZ(0);

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(screened_npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(screened_npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.data(thread_id);

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * screened_npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(screened_weights, \
                        grhoa, ggrada, ggradab, ngrada, gradax, graday, gradaz, gdenx, gdeny, gdenz, \
                        gdenxx, gdenxy, gdenxz, gdenyx, gdenyy, gdenyz, gdenzx, gdenzy, gdenzz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto prefac = screened_weights[g] * grhoa[g];

                    gatmx += prefac * gdenx[atom_g];

                    gatmy += prefac * gdeny[atom_g];

                    gatmz += prefac * gdenz[atom_g];

                    prefac = screened_weights[g] * (ggrada[g] / ngrada[g] + ggradab[g]);

                    gatmx += prefac * (gradax[g] * gdenxx[atom_g] + graday[g] * gdenyx[atom_g] + gradaz[g] * gdenzx[atom_g]);

                    gatmy += prefac * (gradax[g] * gdenxy[atom_g] + graday[g] * gdenyy[atom_g] + gradaz[g] * gdenzy[atom_g]);

                    gatmz += prefac * (gradax[g] * gdenxz[atom_g] + graday[g] * gdenyz[atom_g] + gradaz[g] * gdenzz[atom_g]);
                }

                // factor of 2 from sum of alpha and beta contributions

                gatm[iatom * 3 + 0] += 2.0 * gatmx;

                gatm[iatom * 3 + 1] += 2.0 * gatmy;

                gatm[iatom * 3 + 2] += 2.0 * gatmz;
            }
        }

        timer.stop("Accumulate gradient");
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    std::cout << "Timing of new integrator" << std::endl;
    std::cout << "------------------------" << std::endl;
    std::cout << timer.getSummary() << std::endl;

    CDenseMatrix molgrad(natoms, 3);

    for (int32_t iatom = 0; iatom < natoms; iatom++)
    {
        for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
        {
            molgrad.row(iatom)[0] += molgrad_threads.data(thread_id)[iatom * 3 + 0];

            molgrad.row(iatom)[1] += molgrad_threads.data(thread_id)[iatom * 3 + 1];

            molgrad.row(iatom)[2] += molgrad_threads.data(thread_id)[iatom * 3 + 2];
        }
    }

    return molgrad;
}

CDenseMatrix
CXCNewMolecularGradient::_integrateFxcGradientForLDA(const CMolecule&        molecule,
                                                     const CMolecularBasis&  basis,
                                                     const CAODensityMatrix& rwDensityMatrixOne,
                                                     const CAODensityMatrix& rwDensityMatrixTwo,
                                                     const CAODensityMatrix& gsDensityMatrix,
                                                     const CMolecularGrid&   molecularGrid,
                                                     const CXCFunctional&    xcFunctional) const
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // AO-to-atom mapping

    std::vector<int32_t> ao_to_atom_ids(naos);

    _computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.getNumberOfAtoms();

    CMemBlock2D<double> molgrad_threads(natoms * 3,  nthreads);

    // memory blocks for GTOs on grid points

    CMemBlock2D<double> gaos(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaox(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaoy(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaoz(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos);  // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

    // indices for keeping track of valid grid points

    std::vector<int32_t> screened_point_inds(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> screened_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto screened_weights = screened_weights_data.data();

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

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtovec, xcoords, ycoords, zcoords,

                                             gridblockpos, grid_batch_offset, grid_batch_size, skip_cgto_ids);
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

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        auto rw_sub_dens_mat_one = submat::getSubDensityMatrix(rwDensityMatrixOne, 0, "ALPHA", aoinds, aocount, naos);

        auto rw_sub_dens_mat_two = submat::getSubDensityMatrix(rwDensityMatrixTwo, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto gsdengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, gs_sub_dens_mat, xcfuntype, timer);

        // screen density grid, weights and GTO matrix

        timer.start("Density screening");

        CDensityGrid screened_gsdengrid(gsdengrid);

        gridscreen::screenDensityGridForLDA(screened_point_inds, screened_gsdengrid, gsdengrid,

                                            _screeningThresholdForDensityValues);

        auto screened_npoints = screened_gsdengrid.getNumberOfGridPoints();

        gridscreen::screenWeights(screened_weights, gridblockpos, weights, screened_point_inds, screened_npoints);

        CDenseMatrix screened_mat_chi(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_x(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_y(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_z(mat_chi.getNumberOfRows(), screened_npoints);

        gridscreen::screenGtoMatrixForGGA(screened_mat_chi, screened_mat_chi_x, screened_mat_chi_y, screened_mat_chi_z,

                                          mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, screened_point_inds, screened_npoints);

        timer.stop("Density screening");

        if (screened_npoints == 0) continue;

        // generate density grid

        auto rwdengrid = dengridgen::generateDensityGridForLDA(screened_npoints, screened_mat_chi, rw_sub_dens_mat_one,

                                                               xcfuntype, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, screened_npoints);

        CDenseMatrix dengrady(natoms, screened_npoints);

        CDenseMatrix dengradz(natoms, screened_npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(rw_sub_dens_mat_two, screened_mat_chi);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto naos = screened_mat_chi.getNumberOfRows();

        auto F_val = mat_F.values();

        auto chi_x_val = screened_mat_chi_x.values();

        auto chi_y_val = screened_mat_chi_y.values();

        auto chi_z_val = screened_mat_chi_z.values();

        auto gdenx = dengradx.values();

        auto gdeny = dengrady.values();

        auto gdenz = dengradz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(screened_npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(screened_npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * screened_npoints;

                auto nu_offset = nu * screened_npoints;

                #pragma omp simd aligned(gdenx, gdeny, gdenz, \
                        F_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];

                    gdeny[atom_g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];

                    gdenz[atom_g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];
                }
            }
        }

        timer.stop("Density grad. grid rho");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        CXCHessianGrid vxc2grid(screened_npoints, screened_gsdengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxc2grid, screened_gsdengrid);

        timer.stop("XC functional eval.");

        // pointers to exchange-correlation functional derivative

        auto grho_aa = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::rhoa);

        auto grho_ab = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::rhob);

        // pointers to perturbed density gradient norms

        auto rhowa = rwdengrid.alphaDensity(0);

        auto rhowb = rwdengrid.betaDensity(0);

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(screened_npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(screened_npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.data(thread_id);

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * screened_npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(screened_weights, \
                        grho_aa, grho_ab, rhowa, rhowb, gdenx, gdeny, gdenz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac = screened_weights[g] * (grho_aa[g] * rhowa[g] + grho_ab[g] * rhowb[g]);

                    gatmx += prefac * gdenx[atom_g];

                    gatmy += prefac * gdeny[atom_g];

                    gatmz += prefac * gdenz[atom_g];
                }

                // factor of 2 from sum of alpha and beta contributions

                gatm[iatom * 3 + 0] += 2.0 * gatmx;

                gatm[iatom * 3 + 1] += 2.0 * gatmy;

                gatm[iatom * 3 + 2] += 2.0 * gatmz;
            }
        }

        timer.stop("Accumulate gradient");
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    std::cout << "Timing of new integrator" << std::endl;
    std::cout << "------------------------" << std::endl;
    std::cout << timer.getSummary() << std::endl;

    CDenseMatrix molgrad(natoms, 3);

    for (int32_t iatom = 0; iatom < natoms; iatom++)
    {
        for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
        {
            molgrad.row(iatom)[0] += molgrad_threads.data(thread_id)[iatom * 3 + 0];

            molgrad.row(iatom)[1] += molgrad_threads.data(thread_id)[iatom * 3 + 1];

            molgrad.row(iatom)[2] += molgrad_threads.data(thread_id)[iatom * 3 + 2];
        }
    }

    return molgrad;
}

CDenseMatrix
CXCNewMolecularGradient::_integrateFxcGradientForGGA(const CMolecule&        molecule,
                                                     const CMolecularBasis&  basis,
                                                     const CAODensityMatrix& rwDensityMatrixOne,
                                                     const CAODensityMatrix& rwDensityMatrixTwo,
                                                     const CAODensityMatrix& gsDensityMatrix,
                                                     const CMolecularGrid&   molecularGrid,
                                                     const CXCFunctional&    xcFunctional) const
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // AO-to-atom mapping

    std::vector<int32_t> ao_to_atom_ids(naos);

    _computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.getNumberOfAtoms();

    CMemBlock2D<double> molgrad_threads(natoms * 3,  nthreads);

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

    CMemBlock<int32_t> skip_cgto_ids(naos);  // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

    // indices for keeping track of valid grid points

    std::vector<int32_t> screened_point_inds(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> screened_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto screened_weights = screened_weights_data.data();

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

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForMetaGGA(gaos, gaox, gaoy, gaoz, gaoxx, gaoxy, gaoxz, gaoyy, gaoyz, gaozz,

                                                 gtovec, xcoords, ycoords, zcoords,

                                                 gridblockpos, grid_batch_offset, grid_batch_size, skip_cgto_ids);
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

        CDenseMatrix mat_chi_xy(aocount, npoints);

        CDenseMatrix mat_chi_xz(aocount, npoints);

        CDenseMatrix mat_chi_yy(aocount, npoints);

        CDenseMatrix mat_chi_yz(aocount, npoints);

        CDenseMatrix mat_chi_zz(aocount, npoints);

        for (int32_t i = 0; i < aocount; i++)
        {
            std::memcpy(mat_chi.row(i), gaos.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_x.row(i), gaox.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_y.row(i), gaoy.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_z.row(i), gaoz.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_xx.row(i), gaoxx.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_xy.row(i), gaoxy.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_xz.row(i), gaoxz.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_yy.row(i), gaoyy.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_yz.row(i), gaoyz.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_zz.row(i), gaozz.data(aoinds[i]), npoints * sizeof(double));
        }

        timer.stop("GTO screening");

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        auto rw_sub_dens_mat_one = submat::getSubDensityMatrix(rwDensityMatrixOne, 0, "ALPHA", aoinds, aocount, naos);

        auto rw_sub_dens_mat_two = submat::getSubDensityMatrix(rwDensityMatrixTwo, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto gsdengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                                               gs_sub_dens_mat, xcfuntype, timer);

        // screen density grid, weights and GTO matrix

        timer.start("Density screening");

        CDensityGrid screened_gsdengrid(gsdengrid);

        gridscreen::screenDensityGridForGGA(screened_point_inds, screened_gsdengrid, gsdengrid,

                                            _screeningThresholdForDensityValues);

        auto screened_npoints = screened_gsdengrid.getNumberOfGridPoints();

        gridscreen::screenWeights(screened_weights, gridblockpos, weights, screened_point_inds, screened_npoints);

        CDenseMatrix screened_mat_chi(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_x(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_y(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_z(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_xx(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_xy(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_xz(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_yy(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_yz(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_zz(mat_chi.getNumberOfRows(), screened_npoints);

        gridscreen::screenGtoMatrixForMetaGGA(screened_mat_chi, screened_mat_chi_x, screened_mat_chi_y, screened_mat_chi_z,

                                              screened_mat_chi_xx, screened_mat_chi_xy, screened_mat_chi_xz,

                                              screened_mat_chi_yy, screened_mat_chi_yz, screened_mat_chi_zz,

                                              mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                              mat_chi_xx, mat_chi_xy, mat_chi_xz, mat_chi_yy, mat_chi_yz, mat_chi_zz,

                                              screened_point_inds, screened_npoints);

        timer.stop("Density screening");

        if (screened_npoints == 0) continue;

        // generate density grid

        auto rwdengrid = dengridgen::generateDensityGridForGGA(screened_npoints, screened_mat_chi,

                                                               screened_mat_chi_x, screened_mat_chi_y, screened_mat_chi_z,

                                                               rw_sub_dens_mat_one, xcfuntype, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, screened_npoints);

        CDenseMatrix dengrady(natoms, screened_npoints);

        CDenseMatrix dengradz(natoms, screened_npoints);

        CDenseMatrix dengradxx(natoms, screened_npoints);

        CDenseMatrix dengradxy(natoms, screened_npoints);

        CDenseMatrix dengradxz(natoms, screened_npoints);

        CDenseMatrix dengradyx(natoms, screened_npoints);

        CDenseMatrix dengradyy(natoms, screened_npoints);

        CDenseMatrix dengradyz(natoms, screened_npoints);

        CDenseMatrix dengradzx(natoms, screened_npoints);

        CDenseMatrix dengradzy(natoms, screened_npoints);

        CDenseMatrix dengradzz(natoms, screened_npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(rw_sub_dens_mat_two, screened_mat_chi);

        auto mat_F_x = denblas::multAB(rw_sub_dens_mat_two, screened_mat_chi_x);

        auto mat_F_y = denblas::multAB(rw_sub_dens_mat_two, screened_mat_chi_y);

        auto mat_F_z = denblas::multAB(rw_sub_dens_mat_two, screened_mat_chi_z);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto naos = screened_mat_chi.getNumberOfRows();

        auto F_val = mat_F.values();

        auto F_x_val = mat_F_x.values();

        auto F_y_val = mat_F_y.values();

        auto F_z_val = mat_F_z.values();

        auto chi_x_val = screened_mat_chi_x.values();

        auto chi_y_val = screened_mat_chi_y.values();

        auto chi_z_val = screened_mat_chi_z.values();

        auto chi_xx_val = screened_mat_chi_xx.values();

        auto chi_xy_val = screened_mat_chi_xy.values();

        auto chi_xz_val = screened_mat_chi_xz.values();

        auto chi_yy_val = screened_mat_chi_yy.values();

        auto chi_yz_val = screened_mat_chi_yz.values();

        auto chi_zz_val = screened_mat_chi_zz.values();

        auto gdenx = dengradx.values();

        auto gdeny = dengrady.values();

        auto gdenz = dengradz.values();

        auto gdenxx = dengradxx.values();

        auto gdenxy = dengradxy.values();

        auto gdenxz = dengradxz.values();

        auto gdenyx = dengradyx.values();

        auto gdenyy = dengradyy.values();

        auto gdenyz = dengradyz.values();

        auto gdenzx = dengradzx.values();

        auto gdenzy = dengradzy.values();

        auto gdenzz = dengradzz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(screened_npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(screened_npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * screened_npoints;

                auto nu_offset = nu * screened_npoints;

                #pragma omp simd aligned(gdenx, gdeny, gdenz, \
                        gdenxx, gdenxy, gdenxz, gdenyx, gdenyy, gdenyz, gdenzx, gdenzy, gdenzz, \
                        F_val, F_x_val, F_y_val, F_z_val, chi_x_val, chi_y_val, chi_z_val, \
                        chi_xx_val, chi_xy_val, chi_xz_val, chi_yy_val, chi_yz_val, chi_zz_val : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];

                    gdeny[atom_g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];

                    gdenz[atom_g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];

                    gdenxx[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xx_val[nu_g]);

                    gdenxy[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_xy_val[nu_g]);

                    gdenxz[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_xz_val[nu_g]);

                    gdenyx[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xy_val[nu_g]);

                    gdenyy[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_yy_val[nu_g]);

                    gdenyz[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_yz_val[nu_g]);

                    gdenzx[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xz_val[nu_g]);

                    gdenzy[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_yz_val[nu_g]);

                    gdenzz[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_zz_val[nu_g]);
                }
            }
        }

        timer.stop("Density grad. grid rho");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        CXCGradientGrid vxcgrid(screened_npoints, screened_gsdengrid.getDensityGridType(), xcfuntype);

        CXCHessianGrid vxc2grid(screened_npoints, screened_gsdengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxcgrid, screened_gsdengrid);

        xcFunctional.compute(vxc2grid, screened_gsdengrid);

        timer.stop("XC functional eval.");

        // pointers to exchange-correlation functional derrivatives

        auto ggrad_a = vxcgrid.xcGradientValues(xcvars::grada);

        auto ggrad_c = vxcgrid.xcGradientValues(xcvars::gradab);

        auto grho_aa = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::rhoa);

        auto grho_ab = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::rhob);

        auto gmix_aa = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::grada);

        auto gmix_ab = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::gradb);

        auto gmix_ac = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::gradab);

        auto gmix_bc = vxc2grid.xcHessianValues(xcvars::rhob, xcvars::gradab);

        auto ggrad_aa = vxc2grid.xcHessianValues(xcvars::grada, xcvars::grada);

        auto ggrad_ab = vxc2grid.xcHessianValues(xcvars::grada, xcvars::gradb);

        auto ggrad_ac = vxc2grid.xcHessianValues(xcvars::grada, xcvars::gradab);

        auto ggrad_bc = vxc2grid.xcHessianValues(xcvars::gradb, xcvars::gradab);

        auto ggrad_cc = vxc2grid.xcHessianValues(xcvars::gradab, xcvars::gradab);

        // pointers to density gradient norms

        auto ngrada = screened_gsdengrid.alphaDensityGradient(0);

        auto ngradb = screened_gsdengrid.betaDensityGradient(0);

        auto grada_x = screened_gsdengrid.alphaDensityGradientX(0);

        auto grada_y = screened_gsdengrid.alphaDensityGradientY(0);

        auto grada_z = screened_gsdengrid.alphaDensityGradientZ(0);

        auto gradb_x = screened_gsdengrid.betaDensityGradientX(0);

        auto gradb_y = screened_gsdengrid.betaDensityGradientY(0);

        auto gradb_z = screened_gsdengrid.betaDensityGradientZ(0);

        // pointers to perturbed density gradient norms

        auto rhowa = rwdengrid.alphaDensity(0);

        auto rhowb = rwdengrid.betaDensity(0);

        auto gradwa_x = rwdengrid.alphaDensityGradientX(0);

        auto gradwa_y = rwdengrid.alphaDensityGradientY(0);

        auto gradwa_z = rwdengrid.alphaDensityGradientZ(0);

        auto gradwb_x = rwdengrid.betaDensityGradientX(0);

        auto gradwb_y = rwdengrid.betaDensityGradientY(0);

        auto gradwb_z = rwdengrid.betaDensityGradientZ(0);

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(screened_npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(screened_npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.data(thread_id);

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * screened_npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(screened_weights, \
                        ggrad_a, ggrad_c, grho_aa, grho_ab, gmix_aa, gmix_ab, gmix_ac, gmix_bc, \
                        ggrad_aa, ggrad_ab, ggrad_ac, ggrad_bc, ggrad_cc, ngrada, ngradb, \
                        grada_x, grada_y, grada_z, gradb_x, gradb_y, gradb_z, rhowa, rhowb, \
                        gradwa_x, gradwa_y, gradwa_z, gradwb_x, gradwb_y, gradwb_z, gdenx, gdeny, gdenz, \
                        gdenxx, gdenxy, gdenxz, gdenyx, gdenyy, gdenyz, gdenzx, gdenzy, gdenzz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double w = screened_weights[g];

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

                    double zetac = grada_x[g] * rxwb + grada_y[g] * rywb

                                   + grada_z[g] * rzwb + gradb_x[g] * rxwa

                                   + gradb_y[g] * rywa + gradb_z[g] * rzwa;

                    // contribution from \nabla_A (\phi_mu \phi_nu)

                    // first contribution

                    double fac0 = gmix_aa[g] * zetaa + gmix_ab[g] * zetab

                                  + gmix_ac[g] * zetac + grho_aa[g] * rhowa[g]

                                  + grho_ab[g] * rhowb[g];

                    gatmx += w * fac0 * gdenx[atom_g];

                    gatmy += w * fac0 * gdeny[atom_g];

                    gatmz += w * fac0 * gdenz[atom_g];

                    // contribution from \nabla_A (\nabla (\phi_mu \phi_nu))

                    double xcomp = 0.0;

                    double ycomp = 0.0;

                    double zcomp = 0.0;

                    // second contribution

                    double facr = gmix_aa[g] * rhowa[g] + gmix_ab[g] * rhowb[g]

                                  + ggrad_aa[g] * zetaa + ggrad_ab[g] * zetab + ggrad_ac[g] * zetac;

                    xcomp += facr * rxa;

                    ycomp += facr * rya;

                    zcomp += facr * rza;

                    // third contribution

                    double facz = gmix_ac[g] * rhowa[g] + gmix_bc[g] * rhowb[g]

                                  + ggrad_ac[g] * zetaa + ggrad_bc[g] * zetab + ggrad_cc[g] * zetac;

                    xcomp += facz * grada_x[g];

                    ycomp += facz * grada_y[g];

                    zcomp += facz * grada_z[g];

                    // fourth contribution

                    double prefac = znva * ggrad_a[g];

                    xcomp += prefac * (rxwa - rxa * zetaa);

                    ycomp += prefac * (rywa - rya * zetaa);

                    zcomp += prefac * (rzwa - rza * zetaa);

                    // fifth contribution

                    xcomp += ggrad_c[g] * rxwa;

                    ycomp += ggrad_c[g] * rywa;

                    zcomp += ggrad_c[g] * rzwa;

                    gatmx += w * (xcomp * gdenxx[atom_g] + ycomp * gdenyx[atom_g] + zcomp * gdenzx[atom_g]);

                    gatmy += w * (xcomp * gdenxy[atom_g] + ycomp * gdenyy[atom_g] + zcomp * gdenzy[atom_g]);

                    gatmz += w * (xcomp * gdenxz[atom_g] + ycomp * gdenyz[atom_g] + zcomp * gdenzz[atom_g]);
                }

                // factor of 2 from sum of alpha and beta contributions

                gatm[iatom * 3 + 0] += 2.0 * gatmx;

                gatm[iatom * 3 + 1] += 2.0 * gatmy;

                gatm[iatom * 3 + 2] += 2.0 * gatmz;
            }
        }

        timer.stop("Accumulate gradient");
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    std::cout << "Timing of new integrator" << std::endl;
    std::cout << "------------------------" << std::endl;
    std::cout << timer.getSummary() << std::endl;

    CDenseMatrix molgrad(natoms, 3);

    for (int32_t iatom = 0; iatom < natoms; iatom++)
    {
        for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
        {
            molgrad.row(iatom)[0] += molgrad_threads.data(thread_id)[iatom * 3 + 0];

            molgrad.row(iatom)[1] += molgrad_threads.data(thread_id)[iatom * 3 + 1];

            molgrad.row(iatom)[2] += molgrad_threads.data(thread_id)[iatom * 3 + 2];
        }
    }

    return molgrad;
}

CDenseMatrix
CXCNewMolecularGradient::_integrateKxcGradientForLDA(const CMolecule&        molecule,
                                                     const CMolecularBasis&  basis,
                                                     const CAODensityMatrix& rwDensityMatrixOne,
                                                     const CAODensityMatrix& rwDensityMatrixTwo,
                                                     const CAODensityMatrix& gsDensityMatrix,
                                                     const CMolecularGrid&   molecularGrid,
                                                     const CXCFunctional&    xcFunctional) const
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // AO-to-atom mapping

    std::vector<int32_t> ao_to_atom_ids(naos);

    _computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.getNumberOfAtoms();

    CMemBlock2D<double> molgrad_threads(natoms * 3,  nthreads);

    // memory blocks for GTOs on grid points

    CMemBlock2D<double> gaos(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaox(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaoy(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    CMemBlock2D<double> gaoz(molecularGrid.getMaxNumberOfGridPointsPerBox(), naos);

    // indices for keeping track of GTOs

    // skip_cgto_ids: whether a CGTO should be skipped
    // skip_ao_ids: whether an AO should be skipped
    // aoinds: mapping between AO indices before and after screening

    CMemBlock<int32_t> skip_cgto_ids(naos);  // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

    // indices for keeping track of valid grid points

    std::vector<int32_t> screened_point_inds(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> screened_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto screened_weights = screened_weights_data.data();

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

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtovec, xcoords, ycoords, zcoords,

                                             gridblockpos, grid_batch_offset, grid_batch_size, skip_cgto_ids);
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

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        auto rw_sub_dens_mat_one = submat::getSubDensityMatrix(rwDensityMatrixOne, 0, "ALPHA", aoinds, aocount, naos);

        auto rw_sub_dens_mat_two = submat::getSubDensityMatrix(rwDensityMatrixTwo, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto gsdengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, gs_sub_dens_mat, xcfuntype, timer);

        // screen density grid, weights and GTO matrix

        timer.start("Density screening");

        CDensityGrid screened_gsdengrid(gsdengrid);

        gridscreen::screenDensityGridForLDA(screened_point_inds, screened_gsdengrid, gsdengrid,

                                            _screeningThresholdForDensityValues);

        auto screened_npoints = screened_gsdengrid.getNumberOfGridPoints();

        gridscreen::screenWeights(screened_weights, gridblockpos, weights, screened_point_inds, screened_npoints);

        CDenseMatrix screened_mat_chi(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_x(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_y(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_z(mat_chi.getNumberOfRows(), screened_npoints);

        gridscreen::screenGtoMatrixForGGA(screened_mat_chi, screened_mat_chi_x, screened_mat_chi_y, screened_mat_chi_z,

                                          mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, screened_point_inds, screened_npoints);

        timer.stop("Density screening");

        if (screened_npoints == 0) continue;

        // compute perturbed density

        // prepare rwdenmat for quadratic response

        timer.start("Density grid quad");

        CDenseMatrix zero_sub_den_mat_one(rw_sub_dens_mat_one);

        CDenseMatrix zero_sub_den_mat_two(rw_sub_dens_mat_two);

        zero_sub_den_mat_one.zero();

        zero_sub_den_mat_two.zero();

        CAODensityMatrix rwdenmat(std::vector<CDenseMatrix>({rw_sub_dens_mat_one, zero_sub_den_mat_one,

                                                             rw_sub_dens_mat_two, zero_sub_den_mat_two}), denmat::rest);

        // Note: We use quadratic response (quadMode == "QRF") to calculate
        // third-order functional derivative contribution. The rw2DensityMatrix
        // contains zero matrices and is therefore removed from the following code.
        // Same for rw2dengrid.

        // For "QRF" we have rwDensityMatrix.getNumberOfDensityMatrices() ==
        // 2 * rw2DensityMatrix.getNumberOfDensityMatrices()

        std::string quadMode("QRF");

        auto numdens_rw2 = rwdenmat.getNumberOfDensityMatrices() / 2;

        auto rwdengrid = dengridgen::generateDensityGridForLDA(screened_npoints, screened_mat_chi, rwdenmat, xcfuntype, timer);

        CDensityGridQuad rwdengridquad(screened_npoints, numdens_rw2, xcfuntype, dengrid::ab);

        rwdengridquad.DensityProd(rwdengrid, xcfuntype, numdens_rw2, quadMode);

        timer.stop("Density grid quad");

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, screened_npoints);

        CDenseMatrix dengrady(natoms, screened_npoints);

        CDenseMatrix dengradz(natoms, screened_npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(gs_sub_dens_mat, screened_mat_chi);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto naos = screened_mat_chi.getNumberOfRows();

        auto F_val = mat_F.values();

        auto chi_x_val = screened_mat_chi_x.values();

        auto chi_y_val = screened_mat_chi_y.values();

        auto chi_z_val = screened_mat_chi_z.values();

        auto gdenx = dengradx.values();

        auto gdeny = dengrady.values();

        auto gdenz = dengradz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(screened_npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(screened_npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * screened_npoints;

                auto nu_offset = nu * screened_npoints;

                #pragma omp simd aligned(gdenx, gdeny, gdenz, \
                        F_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];

                    gdeny[atom_g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];

                    gdenz[atom_g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];
                }
            }
        }

        timer.stop("Density grad. grid rho");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        CXCCubicHessianGrid vxc3grid(screened_npoints, screened_gsdengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxc3grid, screened_gsdengrid);

        timer.stop("XC functional eval.");

        // pointers to exchange-correlation functional derivative

        auto df3000 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhoa);

        auto df2100 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhob);

        auto df1200 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::rhob);

        // pointers to perturbed density gradient norms

        auto rhow1a = rwdengridquad.rhow1rhow2(0);

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(screened_npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(screened_npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.data(thread_id);

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * screened_npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(screened_weights, \
                        df3000, df2100, df1200, rhow1a, gdenx, gdeny, gdenz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac = screened_weights[g] * (df3000[g] + df2100[g] + df2100[g] + df1200[g]) * rhow1a[g];

                    gatmx += prefac * gdenx[atom_g];

                    gatmy += prefac * gdeny[atom_g];

                    gatmz += prefac * gdenz[atom_g];
                }

                // factor of 2 from sum of alpha and beta contributions
                // factor of 0.25 from quadratic response

                gatm[iatom * 3 + 0] += 0.25 * (2.0 * gatmx);

                gatm[iatom * 3 + 1] += 0.25 * (2.0 * gatmy);

                gatm[iatom * 3 + 2] += 0.25 * (2.0 * gatmz);
            }
        }

        timer.stop("Accumulate gradient");
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    std::cout << "Timing of new integrator" << std::endl;
    std::cout << "------------------------" << std::endl;
    std::cout << timer.getSummary() << std::endl;

    CDenseMatrix molgrad(natoms, 3);

    for (int32_t iatom = 0; iatom < natoms; iatom++)
    {
        for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
        {
            molgrad.row(iatom)[0] += molgrad_threads.data(thread_id)[iatom * 3 + 0];

            molgrad.row(iatom)[1] += molgrad_threads.data(thread_id)[iatom * 3 + 1];

            molgrad.row(iatom)[2] += molgrad_threads.data(thread_id)[iatom * 3 + 2];
        }
    }

    return molgrad;
}

CDenseMatrix
CXCNewMolecularGradient::_integrateKxcGradientForGGA(const CMolecule&        molecule,
                                                     const CMolecularBasis&  basis,
                                                     const CAODensityMatrix& rwDensityMatrixOne,
                                                     const CAODensityMatrix& rwDensityMatrixTwo,
                                                     const CAODensityMatrix& gsDensityMatrix,
                                                     const CMolecularGrid&   molecularGrid,
                                                     const CXCFunctional&    xcFunctional) const
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    // GTOs container and number of AOs

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    // AO-to-atom mapping

    std::vector<int32_t> ao_to_atom_ids(naos);

    _computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.getNumberOfAtoms();

    CMemBlock2D<double> molgrad_threads(natoms * 3,  nthreads);

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

    CMemBlock<int32_t> skip_cgto_ids(naos);  // note: naos >= ncgtos

    CMemBlock<int32_t> skip_ao_ids(naos);

    std::vector<int32_t> aoinds(naos);

    // indices for keeping track of valid grid points

    std::vector<int32_t> screened_point_inds(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> screened_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto screened_weights = screened_weights_data.data();

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

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForMetaGGA(gaos, gaox, gaoy, gaoz, gaoxx, gaoxy, gaoxz, gaoyy, gaoyz, gaozz,

                                                 gtovec, xcoords, ycoords, zcoords,

                                                 gridblockpos, grid_batch_offset, grid_batch_size, skip_cgto_ids);
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

        CDenseMatrix mat_chi_xy(aocount, npoints);

        CDenseMatrix mat_chi_xz(aocount, npoints);

        CDenseMatrix mat_chi_yy(aocount, npoints);

        CDenseMatrix mat_chi_yz(aocount, npoints);

        CDenseMatrix mat_chi_zz(aocount, npoints);

        for (int32_t i = 0; i < aocount; i++)
        {
            std::memcpy(mat_chi.row(i), gaos.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_x.row(i), gaox.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_y.row(i), gaoy.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_z.row(i), gaoz.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_xx.row(i), gaoxx.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_xy.row(i), gaoxy.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_xz.row(i), gaoxz.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_yy.row(i), gaoyy.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_yz.row(i), gaoyz.data(aoinds[i]), npoints * sizeof(double));

            std::memcpy(mat_chi_zz.row(i), gaozz.data(aoinds[i]), npoints * sizeof(double));
        }

        timer.stop("GTO screening");

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = submat::getSubDensityMatrix(gsDensityMatrix, 0, "ALPHA", aoinds, aocount, naos);

        auto rw_sub_dens_mat_one = submat::getSubDensityMatrix(rwDensityMatrixOne, 0, "ALPHA", aoinds, aocount, naos);

        auto rw_sub_dens_mat_two = submat::getSubDensityMatrix(rwDensityMatrixTwo, 0, "ALPHA", aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto gsdengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                                               gs_sub_dens_mat, xcfuntype, timer);

        // screen density grid, weights and GTO matrix

        timer.start("Density screening");

        CDensityGrid screened_gsdengrid(gsdengrid);

        gridscreen::screenDensityGridForGGA(screened_point_inds, screened_gsdengrid, gsdengrid,

                                            _screeningThresholdForDensityValues);

        auto screened_npoints = screened_gsdengrid.getNumberOfGridPoints();

        gridscreen::screenWeights(screened_weights, gridblockpos, weights, screened_point_inds, screened_npoints);

        CDenseMatrix screened_mat_chi(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_x(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_y(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_z(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_xx(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_xy(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_xz(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_yy(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_yz(mat_chi.getNumberOfRows(), screened_npoints);

        CDenseMatrix screened_mat_chi_zz(mat_chi.getNumberOfRows(), screened_npoints);

        gridscreen::screenGtoMatrixForMetaGGA(screened_mat_chi, screened_mat_chi_x, screened_mat_chi_y, screened_mat_chi_z,

                                              screened_mat_chi_xx, screened_mat_chi_xy, screened_mat_chi_xz,

                                              screened_mat_chi_yy, screened_mat_chi_yz, screened_mat_chi_zz,

                                              mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                              mat_chi_xx, mat_chi_xy, mat_chi_xz, mat_chi_yy, mat_chi_yz, mat_chi_zz,

                                              screened_point_inds, screened_npoints);

        timer.stop("Density screening");

        if (screened_npoints == 0) continue;

        // compute perturbed density

        // prepare rwdenmat for quadratic response

        timer.start("Density grid quad");

        CDenseMatrix zero_sub_den_mat_one(rw_sub_dens_mat_one);

        CDenseMatrix zero_sub_den_mat_two(rw_sub_dens_mat_two);

        zero_sub_den_mat_one.zero();

        zero_sub_den_mat_two.zero();

        CAODensityMatrix rwdenmat(std::vector<CDenseMatrix>({rw_sub_dens_mat_one, zero_sub_den_mat_one,

                                                             rw_sub_dens_mat_two, zero_sub_den_mat_two}), denmat::rest);

        // Note: We use quadratic response (quadMode == "QRF") to calculate
        // third-order functional derivative contribution. The rw2DensityMatrix
        // contains zero matrices and is therefore removed from the following code.
        // Same for rw2dengrid.

        // For "QRF" we have rwDensityMatrix.getNumberOfDensityMatrices() ==
        // 2 * rw2DensityMatrix.getNumberOfDensityMatrices()

        std::string quadMode("QRF");

        auto numdens_rw2 = rwdenmat.getNumberOfDensityMatrices() / 2;

        auto rwdengrid = dengridgen::generateDensityGridForGGA(screened_npoints, screened_mat_chi,

                                                               screened_mat_chi_x, screened_mat_chi_y, screened_mat_chi_z,

                                                               rwdenmat, xcfuntype, timer);

        CDensityGridQuad rwdengridquad(screened_npoints, numdens_rw2, xcfuntype, dengrid::ab);

        rwdengridquad.DensityProd(rwdengrid, xcfuntype, numdens_rw2, quadMode);

        timer.stop("Density grid quad");

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, screened_npoints);

        CDenseMatrix dengrady(natoms, screened_npoints);

        CDenseMatrix dengradz(natoms, screened_npoints);

        CDenseMatrix dengradxx(natoms, screened_npoints);

        CDenseMatrix dengradxy(natoms, screened_npoints);

        CDenseMatrix dengradxz(natoms, screened_npoints);

        CDenseMatrix dengradyx(natoms, screened_npoints);

        CDenseMatrix dengradyy(natoms, screened_npoints);

        CDenseMatrix dengradyz(natoms, screened_npoints);

        CDenseMatrix dengradzx(natoms, screened_npoints);

        CDenseMatrix dengradzy(natoms, screened_npoints);

        CDenseMatrix dengradzz(natoms, screened_npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(gs_sub_dens_mat, screened_mat_chi);

        auto mat_F_x = denblas::multAB(gs_sub_dens_mat, screened_mat_chi_x);

        auto mat_F_y = denblas::multAB(gs_sub_dens_mat, screened_mat_chi_y);

        auto mat_F_z = denblas::multAB(gs_sub_dens_mat, screened_mat_chi_z);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto naos = screened_mat_chi.getNumberOfRows();

        auto F_val = mat_F.values();

        auto F_x_val = mat_F_x.values();

        auto F_y_val = mat_F_y.values();

        auto F_z_val = mat_F_z.values();

        auto chi_x_val = screened_mat_chi_x.values();

        auto chi_y_val = screened_mat_chi_y.values();

        auto chi_z_val = screened_mat_chi_z.values();

        auto chi_xx_val = screened_mat_chi_xx.values();

        auto chi_xy_val = screened_mat_chi_xy.values();

        auto chi_xz_val = screened_mat_chi_xz.values();

        auto chi_yy_val = screened_mat_chi_yy.values();

        auto chi_yz_val = screened_mat_chi_yz.values();

        auto chi_zz_val = screened_mat_chi_zz.values();

        auto gdenx = dengradx.values();

        auto gdeny = dengrady.values();

        auto gdenz = dengradz.values();

        auto gdenxx = dengradxx.values();

        auto gdenxy = dengradxy.values();

        auto gdenxz = dengradxz.values();

        auto gdenyx = dengradyx.values();

        auto gdenyy = dengradyy.values();

        auto gdenyz = dengradyz.values();

        auto gdenzx = dengradzx.values();

        auto gdenzy = dengradzy.values();

        auto gdenzz = dengradzz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(screened_npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(screened_npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * screened_npoints;

                auto nu_offset = nu * screened_npoints;

                #pragma omp simd aligned(gdenx, gdeny, gdenz, \
                        gdenxx, gdenxy, gdenxz, gdenyx, gdenyy, gdenyz, gdenzx, gdenzy, gdenzz, \
                        F_val, F_x_val, F_y_val, F_z_val, chi_x_val, chi_y_val, chi_z_val, \
                        chi_xx_val, chi_xy_val, chi_xz_val, chi_yy_val, chi_yz_val, chi_zz_val : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];

                    gdeny[atom_g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];

                    gdenz[atom_g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];

                    gdenxx[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xx_val[nu_g]);

                    gdenxy[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_xy_val[nu_g]);

                    gdenxz[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_xz_val[nu_g]);

                    gdenyx[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xy_val[nu_g]);

                    gdenyy[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_yy_val[nu_g]);

                    gdenyz[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_yz_val[nu_g]);

                    gdenzx[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xz_val[nu_g]);

                    gdenzy[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_yz_val[nu_g]);

                    gdenzz[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_zz_val[nu_g]);
                }
            }
        }

        timer.stop("Density grad. grid rho");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        CXCGradientGrid vxcgrid(screened_npoints, screened_gsdengrid.getDensityGridType(), xcfuntype);

        CXCHessianGrid vxc2grid(screened_npoints, screened_gsdengrid.getDensityGridType(), xcfuntype);

        CXCCubicHessianGrid vxc3grid(screened_npoints, screened_gsdengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxcgrid, screened_gsdengrid);

        xcFunctional.compute(vxc2grid, screened_gsdengrid);

        xcFunctional.compute(vxc3grid, screened_gsdengrid);

        timer.stop("XC functional eval.");

        // pointers to exchange-correlation functional derrivatives

        auto df0010 = vxcgrid.xcGradientValues(xcvars::grada);

        auto df1010 = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::grada);

        auto df1001 = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::gradb);

        auto df10001 = vxc2grid.xcHessianValues(xcvars::rhoa, xcvars::gradab);

        auto df0020 = vxc2grid.xcHessianValues(xcvars::grada, xcvars::grada);

        auto df0011 = vxc2grid.xcHessianValues(xcvars::grada, xcvars::gradb);

        auto df00101 = vxc2grid.xcHessianValues(xcvars::grada, xcvars::gradab);

        auto df00002 = vxc2grid.xcHessianValues(xcvars::gradab, xcvars::gradab);

        auto df00011 = vxc2grid.xcHessianValues(xcvars::gradb, xcvars::gradab);

        auto df01001 = vxc2grid.xcHessianValues(xcvars::rhob, xcvars::gradab);

        auto df0110 = vxc2grid.xcHessianValues(xcvars::rhob, xcvars::grada);

        auto df3000 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhoa);

        auto df2100 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhob);

        auto df1200 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::rhob);

        auto df2010 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::grada);

        auto df0030 = vxc3grid.xcCubicHessianValues(xcvars::grada, xcvars::grada, xcvars::grada);

        auto df0021 = vxc3grid.xcCubicHessianValues(xcvars::grada, xcvars::grada, xcvars::gradb);

        auto df0012 = vxc3grid.xcCubicHessianValues(xcvars::grada, xcvars::gradb, xcvars::gradb);

        auto df00201 = vxc3grid.xcCubicHessianValues(xcvars::grada, xcvars::grada, xcvars::gradab);

        auto df00111 = vxc3grid.xcCubicHessianValues(xcvars::grada, xcvars::gradb, xcvars::gradab);

        auto df00102 = vxc3grid.xcCubicHessianValues(xcvars::grada, xcvars::gradab, xcvars::gradab);

        auto df00003 = vxc3grid.xcCubicHessianValues(xcvars::gradab, xcvars::gradab, xcvars::gradab);

        auto df2001 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::gradb);

        auto df1110 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::grada);

        auto df1101 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::gradb);

        auto df20001 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::gradab);

        auto df11001 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::gradab);

        auto df1020 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::grada, xcvars::grada);

        auto df1011 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::grada, xcvars::gradb);

        auto df1002 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::gradb, xcvars::gradb);

        auto df10101 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::grada, xcvars::gradab);

        auto df10002 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::gradab, xcvars::gradab);

        auto df01002 = vxc3grid.xcCubicHessianValues(xcvars::rhob, xcvars::gradab, xcvars::gradab);

        auto df0120 = vxc3grid.xcCubicHessianValues(xcvars::rhob, xcvars::grada, xcvars::grada);

        auto df0111 = vxc3grid.xcCubicHessianValues(xcvars::rhob, xcvars::grada, xcvars::gradb);

        auto df01101 = vxc3grid.xcCubicHessianValues(xcvars::rhob, xcvars::grada, xcvars::gradab);

        auto df10011 = vxc3grid.xcCubicHessianValues(xcvars::rhoa, xcvars::gradb, xcvars::gradab);

        auto df01011 = vxc3grid.xcCubicHessianValues(xcvars::rhob, xcvars::gradb, xcvars::gradab);

        auto df0210 = vxc3grid.xcCubicHessianValues(xcvars::rhob, xcvars::rhob, xcvars::grada);

        auto df02001 = vxc3grid.xcCubicHessianValues(xcvars::rhob, xcvars::rhob, xcvars::gradab);

        auto df00021 = vxc3grid.xcCubicHessianValues(xcvars::gradb, xcvars::gradb, xcvars::gradab);

        // pointers to density gradient norms

        auto ngrada = screened_gsdengrid.alphaDensityGradient(0);

        auto grada_x = screened_gsdengrid.alphaDensityGradientX(0);

        auto grada_y = screened_gsdengrid.alphaDensityGradientY(0);

        auto grada_z = screened_gsdengrid.alphaDensityGradientZ(0);

        // pointers to perturbed density gradient norms

        auto rhow1rhow2 = rwdengridquad.rhow1rhow2(0);

        auto rxw1rhow2 = rwdengridquad.rxw1rhow2(0);

        auto ryw1rhow2 = rwdengridquad.ryw1rhow2(0);

        auto rzw1rhow2 = rwdengridquad.rzw1rhow2(0);

        auto rxw1rxw2 = rwdengridquad.rxw1rxw2(0);

        auto rxw1ryw2 = rwdengridquad.rxw1ryw2(0);

        auto rxw1rzw2 = rwdengridquad.rxw1rzw2(0);

        auto ryw1rxw2 = rwdengridquad.ryw1rxw2(0);

        auto ryw1ryw2 = rwdengridquad.ryw1ryw2(0);

        auto ryw1rzw2 = rwdengridquad.ryw1rzw2(0);

        auto rzw1rxw2 = rwdengridquad.rzw1rxw2(0);

        auto rzw1ryw2 = rwdengridquad.rzw1ryw2(0);

        auto rzw1rzw2 = rwdengridquad.rzw1rzw2(0);

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(screened_npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(screened_npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.data(thread_id);

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * screened_npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(screened_weights, df0010, \
                        df1010, df1001, df10001, df0020, df0011, df00101, df00002, df00011, \
                        df01001, df0110, df3000, df2100, df1200, df2010, df0030, df0021, \
                        df0012, df00201, df00111, df00102, df00003, df2001, df1110, df1101, \
                        df20001, df11001, df1020, df1011, df1002, df10101, df10002, \
                        df01002, df0120, df0111, df01101, df10011, df01011, df0210, \
                        df02001, df00021, ngrada, grada_x, grada_y, grada_z, rhow1rhow2, \
                        rxw1rhow2, ryw1rhow2, rzw1rhow2, rxw1rxw2, rxw1ryw2, rxw1rzw2, ryw1rxw2, \
                        ryw1ryw2, ryw1rzw2, rzw1rxw2, rzw1ryw2, rzw1rzw2, gdenx, gdeny, gdenz, \
                        gdenxx, gdenxy, gdenxz, gdenyx, gdenyy, gdenyz, gdenzx, gdenzy, gdenzz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double w = screened_weights[g];

                    double znva = 1.0 / ngrada[g];

                    double znva3 = 1.0 / std::pow(ngrada[g], 3.0);

                    double znva5 = 1.0 / std::pow(ngrada[g], 5.0);

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

                    double xigrad_dot_rw1rw2 = xigrad_x * rxw1rhow2[g] + xigrad_y * ryw1rhow2[g] + xigrad_z * rzw1rhow2[g];

                    double rw1_dot_rw2 = rxw1rxw2[g] + ryw1ryw2[g] + rzw1rzw2[g];

                    double xigrad_dot_rw1rhow2 = xigrad_x * rxw1rhow2[g] + xigrad_y * ryw1rhow2[g] + xigrad_z * rzw1rhow2[g];

                    double grad_dot_rw1rw2 = grada_x[g] * rxw1rhow2[g] + grada_y[g] * ryw1rhow2[g] + grada_z[g] * rzw1rhow2[g];

                    double grad_dot_rw1rhow2 = grada_x[g] * rxw1rhow2[g] + grada_y[g] * ryw1rhow2[g] + grada_z[g] * rzw1rhow2[g];

                    double xigrad_dot_rw1_xigrad_dot_rw2 =
                        xigrad_x * xigrad_x * rxw1rxw2[g] + xigrad_x * xigrad_y * rxw1ryw2[g] + xigrad_x * xigrad_z * rxw1rzw2[g] +
                        xigrad_y * xigrad_x * ryw1rxw2[g] + xigrad_y * xigrad_y * ryw1ryw2[g] + xigrad_y * xigrad_z * ryw1rzw2[g] +
                        xigrad_z * xigrad_x * rzw1rxw2[g] + xigrad_z * xigrad_y * rzw1ryw2[g] + xigrad_z * xigrad_z * rzw1rzw2[g];

                    double twelthfifth_gam =
                        (xigrad_x * grada_x[g] + grada_x[g] * xigrad_x) * rxw1rxw2[g] +
                        (xigrad_x * grada_y[g] + grada_x[g] * xigrad_y) * rxw1ryw2[g] +
                        (xigrad_x * grada_z[g] + grada_x[g] * xigrad_z) * rxw1rzw2[g] +
                        (xigrad_y * grada_x[g] + grada_y[g] * xigrad_x) * rxw1rxw2[g] +
                        (xigrad_y * grada_y[g] + grada_y[g] * xigrad_y) * rxw1ryw2[g] +
                        (xigrad_y * grada_z[g] + grada_y[g] * xigrad_z) * rxw1rzw2[g] +
                        (xigrad_z * grada_x[g] + grada_z[g] * xigrad_x) * rxw1rxw2[g] +
                        (xigrad_z * grada_y[g] + grada_z[g] * xigrad_y) * rxw1ryw2[g] +
                        (xigrad_z * grada_z[g] + grada_z[g] * xigrad_z) * rxw1rzw2[g];

                    // contribution from \nabla_A (\phi_mu \phi_nu)

                    double prefac;

                    // fifth
                    // seventh
                    // eighth

                    prefac = (df3000[g] + 2.0 * df2100[g] + df1200[g]) * rhow1rhow2[g]

                             + (df2010[g] + df2001[g]) * xigrad_dot_rw1rw2

                             + (df1110[g] + df1101[g]) * xigrad_dot_rw1rw2

                             + 2.0 * (df20001[g] + df11001[g]) * grad_dot_rw1rw2

                             + (df1020[g] + 2.0 * df1011[g] + df1002[g]) * xigrad_dot_rw1_xigrad_dot_rw2

                             + (df1010[g] + df1001[g]) *

                               (xigrad_xx * rxw1rxw2[g] + xigrad_xy * rxw1ryw2[g] + xigrad_xz * rxw1rzw2[g]

                                + xigrad_xy * ryw1rxw2[g] + xigrad_yy * ryw1ryw2[g] + xigrad_yz * ryw1rzw2[g]

                                + xigrad_xz * rzw1rxw2[g] + xigrad_yz * rzw1ryw2[g] + xigrad_zz * rzw1rzw2[g])

                             + 2.0 * (df10101[g] + df10101[g]) * ngrada[g] * xigrad_dot_rw1_xigrad_dot_rw2

                             + 4.0 * df10002[g] * ngrada[g] * ngrada[g] * xigrad_dot_rw1_xigrad_dot_rw2

                             + 2.0 * df10001[g] * rw1_dot_rw2;

                    gatmx += w * prefac * gdenx[atom_g];

                    gatmy += w * prefac * gdeny[atom_g];

                    gatmz += w * prefac * gdenz[atom_g];

                    // contribution from \nabla_A (\nabla (\phi_mu \phi_nu))

                    double xcomp = 0.0;

                    double ycomp = 0.0;

                    double zcomp = 0.0;

                    // ninth
                    // tenth
                    // twelfth

                    prefac = (df2010[g] + 2.0 * df1110[g] + df0210[g]) * rhow1rhow2[g]

                             + (df1020[g] + df1011[g] + df0120[g] + df0111[g]) * xigrad_dot_rw1rhow2

                             + (df10101[g] + df10011[g] + df01101[g] + df0111[g]) * grad_dot_rw1rhow2

                             + (df0030[g] + 2.0 * df0021[g] + df0012[g]) * xigrad_dot_rw1_xigrad_dot_rw2

                             + (df00101[g] + df00011[g]) * rw1_dot_rw2

                             + (df00201[g] + df00111[g]) * twelthfifth_gam

                             + df00102[g] * xigrad_dot_rw1_xigrad_dot_rw2;

                    xcomp += prefac * xigrad_x;

                    ycomp += prefac * xigrad_y;

                    zcomp += prefac * xigrad_z;

                    // ninth
                    // tenth
                    // twelfth

                    prefac = (df20001[g] + 2.0 * df11001[g] + df02001[g]) * rhow1rhow2[g]

                             + (df10101[g] + df10011[g] + df01101[g] + df0111[g] + df01011[g]) * xigrad_dot_rw1rhow2

                             + (df10002[g] + df01002[g]) * grad_dot_rw1rhow2

                             + df00002[g] * rw1_dot_rw2

                             + (df00201[g] + 2 * df00111[g] + df00021[g]) * xigrad_dot_rw1_xigrad_dot_rw2

                             + (df00102[g] + df00011[g]) * ngrada[g] * xigrad_dot_rw1_xigrad_dot_rw2

                             + df00003[g] * ngrada[g] * ngrada[g];

                    xcomp += prefac * grada_x[g];

                    ycomp += prefac * grada_y[g];

                    zcomp += prefac * grada_z[g];

                    // tenth

                    prefac = (df10001[g] + df01001[g]);

                    xcomp += prefac * rxw1rhow2[g];

                    ycomp += prefac * ryw1rhow2[g];

                    zcomp += prefac * rzw1rhow2[g];

                    // tenth

                    prefac = (df1010[g] + df0110[g]);

                    xcomp += prefac * (xigrad_xx * rxw1rhow2[g] + xigrad_xy * ryw1rhow2[g] + xigrad_xz * rzw1rhow2[g]);

                    ycomp += prefac * (xigrad_xy * rxw1rhow2[g] + xigrad_yy * ryw1rhow2[g] + xigrad_yz * rzw1rhow2[g]);

                    zcomp += prefac * (xigrad_xz * rxw1rhow2[g] + xigrad_yz * ryw1rhow2[g] + xigrad_zz * rzw1rhow2[g]);

                    // twelfth
                    // twelthfirst

                    prefac = df0010[g];

                    xcomp += prefac * (xigrad_xxx * rxw1rxw2[g] + xigrad_xxy * rxw1ryw2[g] + xigrad_xxz * rxw1rzw2[g]

                                       + xigrad_xxy * ryw1rxw2[g] + xigrad_xyy * ryw1ryw2[g] + xigrad_xyz * ryw1rzw2[g]

                                       + xigrad_xxz * rzw1rxw2[g] + xigrad_xyz * rzw1ryw2[g] + xigrad_xzz * rzw1rzw2[g]);

                    ycomp += prefac * (xigrad_xxy * rxw1rxw2[g] + xigrad_xyy * rxw1ryw2[g] + xigrad_xyz * rxw1rzw2[g]

                                       + xigrad_xyy * ryw1rxw2[g] + xigrad_yyy * ryw1ryw2[g] + xigrad_yyz * ryw1rzw2[g]

                                       + xigrad_xyz * rzw1rxw2[g] + xigrad_yyz * rzw1ryw2[g] + xigrad_yzz * rzw1rzw2[g]);

                    zcomp += prefac * (xigrad_xxz * rxw1rxw2[g] + xigrad_xyz * rxw1ryw2[g] + xigrad_xzz * rxw1rzw2[g]

                                       + xigrad_xyz * ryw1rxw2[g] + xigrad_yyz * ryw1ryw2[g] + xigrad_yzz * ryw1rzw2[g]

                                       + xigrad_xzz * rzw1rxw2[g] + xigrad_yzz * rzw1ryw2[g] + xigrad_zzz * rzw1rzw2[g]);

                    // twelfth
                    // twelthsecond

                    prefac = (df0020[g] + df0011[g]) + (df00101[g] + df00011[g]) * ngrada[g];

                    xcomp += prefac * (xigrad_xx * xigrad_x * rxw1rxw2[g]

                                       + xigrad_xy * xigrad_x * rxw1ryw2[g]

                                       + xigrad_xz * xigrad_x * rxw1rzw2[g]

                                       + xigrad_xy * xigrad_x * ryw1rxw2[g]

                                       + xigrad_yy * xigrad_x * ryw1ryw2[g]

                                       + xigrad_yz * xigrad_x * ryw1rzw2[g]

                                       + xigrad_xz * xigrad_x * rzw1rxw2[g]

                                       + xigrad_yz * xigrad_x * rzw1ryw2[g]

                                       + xigrad_zz * xigrad_x * rzw1rzw2[g]);

                    ycomp += prefac * (xigrad_xx * xigrad_y * rxw1rxw2[g]

                                       + xigrad_xy * xigrad_y * rxw1ryw2[g]

                                       + xigrad_xz * xigrad_y * rxw1rzw2[g]

                                       + xigrad_xy * xigrad_y * ryw1rxw2[g]

                                       + xigrad_yy * xigrad_y * ryw1ryw2[g]

                                       + xigrad_yz * xigrad_y * ryw1rzw2[g]

                                       + xigrad_xz * xigrad_y * rzw1rxw2[g]

                                       + xigrad_yz * xigrad_y * rzw1ryw2[g]

                                       + xigrad_zz * xigrad_y * rzw1rzw2[g]);

                    zcomp += prefac * (xigrad_xx * xigrad_z * rxw1rxw2[g]

                                       + xigrad_xy * xigrad_z * rxw1ryw2[g]

                                       + xigrad_xz * xigrad_z * rxw1rzw2[g]

                                       + xigrad_xy * xigrad_z * ryw1rxw2[g]

                                       + xigrad_yy * xigrad_z * ryw1ryw2[g]

                                       + xigrad_yz * xigrad_z * ryw1rzw2[g]

                                       + xigrad_xz * xigrad_z * rzw1rxw2[g]

                                       + xigrad_yz * xigrad_z * rzw1ryw2[g]

                                       + xigrad_zz * xigrad_z * rzw1rzw2[g]);

                    // twelfth
                    // twelththird

                    prefac = (df0020[g] + df0011[g]) + df00101[g] * ngrada[g];

                    xcomp += prefac * (xigrad_xx * xigrad_x * (rxw1rxw2[g] + rxw1rxw2[g])

                                       + xigrad_xx * xigrad_y * (ryw1rxw2[g] + rxw1ryw2[g])

                                       + xigrad_xx * xigrad_z * (rzw1rxw2[g] + rxw1rzw2[g])

                                       + xigrad_xy * xigrad_x * (rxw1ryw2[g] + ryw1rxw2[g])

                                       + xigrad_xy * xigrad_y * (ryw1ryw2[g] + ryw1ryw2[g])

                                       + xigrad_xy * xigrad_z * (rzw1ryw2[g] + ryw1rzw2[g])

                                       + xigrad_xz * xigrad_x * (rxw1rzw2[g] + rzw1rxw2[g])

                                       + xigrad_xz * xigrad_y * (ryw1rzw2[g] + rzw1ryw2[g])

                                       + xigrad_xz * xigrad_z * (rzw1rzw2[g] + rzw1rzw2[g]));

                    ycomp += prefac * (xigrad_xy * xigrad_x * (rxw1rxw2[g] + rxw1rxw2[g])

                                       + xigrad_xy * xigrad_y * (ryw1rxw2[g] + rxw1ryw2[g])

                                       + xigrad_xy * xigrad_z * (rzw1rxw2[g] + rxw1rzw2[g])

                                       + xigrad_yy * xigrad_x * (rxw1ryw2[g] + ryw1rxw2[g])

                                       + xigrad_yy * xigrad_y * (ryw1ryw2[g] + ryw1ryw2[g])

                                       + xigrad_yy * xigrad_z * (rzw1ryw2[g] + ryw1rzw2[g])

                                       + xigrad_yz * xigrad_x * (rxw1rzw2[g] + rzw1rxw2[g])

                                       + xigrad_yz * xigrad_y * (ryw1rzw2[g] + rzw1ryw2[g])

                                       + xigrad_yz * xigrad_z * (rzw1rzw2[g] + rzw1rzw2[g]));

                    zcomp += prefac * (xigrad_xz * xigrad_x * (rxw1rxw2[g] + rxw1rxw2[g])

                                       + xigrad_xz * xigrad_y * (ryw1rxw2[g] + rxw1ryw2[g])

                                       + xigrad_xz * xigrad_z * (rzw1rxw2[g] + rxw1rzw2[g])

                                       + xigrad_yz * xigrad_x * (rxw1ryw2[g] + ryw1rxw2[g])

                                       + xigrad_yz * xigrad_y * (ryw1ryw2[g] + ryw1ryw2[g])

                                       + xigrad_yz * xigrad_z * (rzw1ryw2[g] + ryw1rzw2[g])

                                       + xigrad_zz * xigrad_x * (rxw1rzw2[g] + rzw1rxw2[g])

                                       + xigrad_zz * xigrad_y * (ryw1rzw2[g] + rzw1ryw2[g])

                                       + xigrad_zz * xigrad_z * (rzw1rzw2[g] + rzw1rzw2[g]));

                    // twelfth
                    // twelthfourth_gam =

                    prefac = (df00101[g] + df00011[g]) + df00002[g] * ngrada[g];

                    xcomp += prefac * (xigrad_x * rxw1rxw2[g] + xigrad_y * ryw1rxw2[g] + xigrad_z * rzw1rxw2[g]);

                    ycomp += prefac * (xigrad_x * rxw1ryw2[g] + xigrad_y * ryw1ryw2[g] + xigrad_z * rzw1ryw2[g]);

                    zcomp += prefac * (xigrad_x * rxw1rzw2[g] + xigrad_y * ryw1rzw2[g] + xigrad_z * rzw1rzw2[g]);

                    gatmx += w * (xcomp * gdenxx[atom_g] + ycomp * gdenyx[atom_g] + zcomp * gdenzx[atom_g]);

                    gatmy += w * (xcomp * gdenxy[atom_g] + ycomp * gdenyy[atom_g] + zcomp * gdenzy[atom_g]);

                    gatmz += w * (xcomp * gdenxz[atom_g] + ycomp * gdenyz[atom_g] + zcomp * gdenzz[atom_g]);
                }

                // factor of 2 from sum of alpha and beta contributions
                // factor of 0.25 from quadratic response

                gatm[iatom * 3 + 0] += 0.25 * (2.0 * gatmx);

                gatm[iatom * 3 + 1] += 0.25 * (2.0 * gatmy);

                gatm[iatom * 3 + 2] += 0.25 * (2.0 * gatmz);
            }
        }

        timer.stop("Accumulate gradient");
    }

    // destroy GTOs container

    delete gtovec;

    timer.stop("Total timing");

    std::cout << "Timing of new integrator" << std::endl;
    std::cout << "------------------------" << std::endl;
    std::cout << timer.getSummary() << std::endl;

    CDenseMatrix molgrad(natoms, 3);

    for (int32_t iatom = 0; iatom < natoms; iatom++)
    {
        for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
        {
            molgrad.row(iatom)[0] += molgrad_threads.data(thread_id)[iatom * 3 + 0];

            molgrad.row(iatom)[1] += molgrad_threads.data(thread_id)[iatom * 3 + 1];

            molgrad.row(iatom)[2] += molgrad_threads.data(thread_id)[iatom * 3 + 2];
        }
    }

    return molgrad;
}

void
CXCNewMolecularGradient::_computeAOtoAtomMapping(std::vector<int32_t>&  ao_to_atom_ids,
                                                 const CMolecule&       molecule,
                                                 const CMolecularBasis& basis) const
{
    for (int32_t iatom = 0; iatom < molecule.getNumberOfAtoms(); iatom++)
    {
        CGtoContainer* atomgtovec = new CGtoContainer(molecule, basis, iatom, 1);

        for (int32_t i = 0; i < atomgtovec->getNumberOfGtoBlocks(); i++)
        {
            auto bgtos = atomgtovec->getGtoBlock(i);

            auto bang = bgtos.getAngularMomentum();

            auto bnspher = angmom::to_SphericalComponents(bang);

            for (int32_t j = 0; j < bgtos.getNumberOfContrGtos(); j++)
            {
                for (int32_t k = 0; k < bnspher; k++)
                {
                    auto idx = bgtos.getIdentifiers(k)[j];

                    ao_to_atom_ids[idx] = iatom;
                }
            }
        }

        delete atomgtovec;
    }
}
