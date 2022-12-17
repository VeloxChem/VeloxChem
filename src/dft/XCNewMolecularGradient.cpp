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
#include "NewFunctionalParser.hpp"
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

    auto newfvxc = newvxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = newfvxc.getFunctionalType();

    if (rwDensityMatrix.isClosedShell() && gsDensityMatrix.isClosedShell())
    {
        if (xcfuntype == xcfun::lda)
        {
            return _integrateVxcGradientForLDA(molecule, basis, rwDensityMatrix, gsDensityMatrix, molecularGrid, newfvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
            return _integrateVxcGradientForGGA(molecule, basis, rwDensityMatrix, gsDensityMatrix, molecularGrid, newfvxc);
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

    auto newfvxc = newvxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = newfvxc.getFunctionalType();

    if (rwDensityMatrixOne.isClosedShell() && rwDensityMatrixTwo.isClosedShell() && gsDensityMatrix.isClosedShell())
    {
        if (xcfuntype == xcfun::lda)
        {
            return _integrateFxcGradientForLDA(molecule, basis, rwDensityMatrixOne, rwDensityMatrixTwo,

                                               gsDensityMatrix, molecularGrid, newfvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
            return _integrateFxcGradientForGGA(molecule, basis, rwDensityMatrixOne, rwDensityMatrixTwo,

                                               gsDensityMatrix, molecularGrid, newfvxc);
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

    auto fvxc = newvxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

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
                                                     const CXCNewFunctional& xcFunctional) const
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

    // density and functional derivatives

    CMemBlock<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> exc_data(1 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    CMemBlock<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();

    auto exc = exc_data.data();

    auto vrho = vrho_data.data();

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

        dengridgen::generateDensityForLDA(rho, npoints, mat_chi, gs_sub_dens_mat, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, npoints);

        CDenseMatrix dengrady(natoms, npoints);

        CDenseMatrix dengradz(natoms, npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(rw_sub_dens_mat, mat_chi);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto naos = mat_chi.getNumberOfRows();

        auto F_val = mat_F.values();

        auto chi_x_val = mat_chi_x.values();

        auto chi_y_val = mat_chi_y.values();

        auto chi_z_val = mat_chi_z.values();

        auto gdenx = dengradx.values();

        auto gdeny = dengrady.values();

        auto gdenz = dengradz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * npoints;

                auto nu_offset = nu * npoints;

                #pragma omp simd aligned(gdenx, gdeny, gdenz, F_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
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

        xcFunctional.compute_exc_vxc_for_lda(npoints, rho, exc, vrho);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenVxcFockForLDA(rho, exc, vrho, npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.data(thread_id);

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(local_weights, vrho, gdenx, gdeny, gdenz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac = local_weights[g] * vrho[2 * g + 0];

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

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;

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
                                                     const CXCNewFunctional& xcFunctional) const
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

        dengridgen::generateDensityForGGA(rho, rhograd, sigma, npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                          gs_sub_dens_mat, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, npoints);
        CDenseMatrix dengrady(natoms, npoints);
        CDenseMatrix dengradz(natoms, npoints);

        CDenseMatrix dengradxx(natoms, npoints);
        CDenseMatrix dengradxy(natoms, npoints);
        CDenseMatrix dengradxz(natoms, npoints);

        CDenseMatrix dengradyx(natoms, npoints);
        CDenseMatrix dengradyy(natoms, npoints);
        CDenseMatrix dengradyz(natoms, npoints);

        CDenseMatrix dengradzx(natoms, npoints);
        CDenseMatrix dengradzy(natoms, npoints);
        CDenseMatrix dengradzz(natoms, npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(rw_sub_dens_mat, mat_chi);

        auto mat_F_x = denblas::multAB(rw_sub_dens_mat, mat_chi_x);
        auto mat_F_y = denblas::multAB(rw_sub_dens_mat, mat_chi_y);
        auto mat_F_z = denblas::multAB(rw_sub_dens_mat, mat_chi_z);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto naos = mat_chi.getNumberOfRows();

        auto F_val = mat_F.values();

        auto F_x_val = mat_F_x.values();
        auto F_y_val = mat_F_y.values();
        auto F_z_val = mat_F_z.values();

        auto chi_x_val = mat_chi_x.values();
        auto chi_y_val = mat_chi_y.values();
        auto chi_z_val = mat_chi_z.values();

        auto chi_xx_val = mat_chi_xx.values();
        auto chi_xy_val = mat_chi_xy.values();
        auto chi_xz_val = mat_chi_xz.values();

        auto chi_yy_val = mat_chi_yy.values();
        auto chi_yz_val = mat_chi_yz.values();
        auto chi_zz_val = mat_chi_zz.values();

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

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * npoints;

                auto nu_offset = nu * npoints;

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

        xcFunctional.compute_exc_vxc_for_gga(npoints, rho, sigma, exc, vrho, vsigma);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenVxcFockForGGA(rho, sigma, exc, vrho, vsigma, npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.data(thread_id);

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(local_weights, \
                        rhograd, vrho, vsigma, gdenx, gdeny, gdenz, \
                        gdenxx, gdenxy, gdenxz, gdenyx, gdenyy, gdenyz, gdenzx, gdenzy, gdenzz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto vx = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                    auto vy = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                    auto vz = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                    gatmx += local_weights[g] * vrho[2 * g + 0] * gdenx[atom_g];
                    gatmy += local_weights[g] * vrho[2 * g + 0] * gdeny[atom_g];
                    gatmz += local_weights[g] * vrho[2 * g + 0] * gdenz[atom_g];

                    gatmx += local_weights[g] * (vx * gdenxx[atom_g] + vy * gdenyx[atom_g] + vz * gdenzx[atom_g]);
                    gatmy += local_weights[g] * (vx * gdenxy[atom_g] + vy * gdenyy[atom_g] + vz * gdenzy[atom_g]);
                    gatmz += local_weights[g] * (vx * gdenxz[atom_g] + vy * gdenyz[atom_g] + vz * gdenzz[atom_g]);
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

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;

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
                                                     const CXCNewFunctional& xcFunctional) const
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

        dengridgen::generateDensityForLDA(rho, npoints, mat_chi, gs_sub_dens_mat, timer);

        dengridgen::generateDensityForLDA(rhow, npoints, mat_chi, rw_sub_dens_mat_one, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, npoints);
        CDenseMatrix dengrady(natoms, npoints);
        CDenseMatrix dengradz(natoms, npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(rw_sub_dens_mat_two, mat_chi);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto naos = mat_chi.getNumberOfRows();

        auto F_val = mat_F.values();

        auto chi_x_val = mat_chi_x.values();
        auto chi_y_val = mat_chi_y.values();
        auto chi_z_val = mat_chi_z.values();

        auto gdenx = dengradx.values();
        auto gdeny = dengrady.values();
        auto gdenz = dengradz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * npoints;

                auto nu_offset = nu * npoints;

                #pragma omp simd aligned(gdenx, gdeny, gdenz, F_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
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

        xcFunctional.compute_fxc_for_lda(npoints, rho, v2rho2);

        timer.stop("XC functional eval.");

        // screen density and functional derivatives

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenFxcFockForLDA(rho, v2rho2, npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.data(thread_id);

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(local_weights, \
                        rhow, v2rho2, gdenx, gdeny, gdenz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac = local_weights[g] * (v2rho2[3 * g + 0] * rhow[2 * g + 0] +
                                                        v2rho2[3 * g + 1] * rhow[2 * g + 1]);

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

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;

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
                                                     const CXCNewFunctional& xcFunctional) const
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

        dengridgen::generateDensityForGGA(rho, rhograd, sigma, npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                          gs_sub_dens_mat, timer);

        dengridgen::generateDensityForGGA(rhow, rhowgrad, nullptr, npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                          rw_sub_dens_mat_one, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, npoints);
        CDenseMatrix dengrady(natoms, npoints);
        CDenseMatrix dengradz(natoms, npoints);

        CDenseMatrix dengradxx(natoms, npoints);
        CDenseMatrix dengradxy(natoms, npoints);
        CDenseMatrix dengradxz(natoms, npoints);

        CDenseMatrix dengradyx(natoms, npoints);
        CDenseMatrix dengradyy(natoms, npoints);
        CDenseMatrix dengradyz(natoms, npoints);

        CDenseMatrix dengradzx(natoms, npoints);
        CDenseMatrix dengradzy(natoms, npoints);
        CDenseMatrix dengradzz(natoms, npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(rw_sub_dens_mat_two, mat_chi);

        auto mat_F_x = denblas::multAB(rw_sub_dens_mat_two, mat_chi_x);
        auto mat_F_y = denblas::multAB(rw_sub_dens_mat_two, mat_chi_y);
        auto mat_F_z = denblas::multAB(rw_sub_dens_mat_two, mat_chi_z);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto naos = mat_chi.getNumberOfRows();

        auto F_val = mat_F.values();

        auto F_x_val = mat_F_x.values();
        auto F_y_val = mat_F_y.values();
        auto F_z_val = mat_F_z.values();

        auto chi_x_val = mat_chi_x.values();
        auto chi_y_val = mat_chi_y.values();
        auto chi_z_val = mat_chi_z.values();

        auto chi_xx_val = mat_chi_xx.values();
        auto chi_xy_val = mat_chi_xy.values();
        auto chi_xz_val = mat_chi_xz.values();
        auto chi_yy_val = mat_chi_yy.values();
        auto chi_yz_val = mat_chi_yz.values();
        auto chi_zz_val = mat_chi_zz.values();

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

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * npoints;

                auto nu_offset = nu * npoints;

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

        xcFunctional.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

        xcFunctional.compute_fxc_for_gga(npoints, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

        timer.stop("XC functional eval.");

        // screen density grid, weights and GTO matrix

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenFxcFockForGGA(rho, sigma, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2,

                                        npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.data(thread_id);

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(local_weights, \
                        rhow, rhograd, rhowgrad, vsigma, v2rho2, v2rhosigma, v2sigma2, gdenx, gdeny, gdenz, \
                        gdenxx, gdenxy, gdenxz, gdenyx, gdenyy, gdenyz, gdenzx, gdenzy, gdenzz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double w = local_weights[g];

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

                    // scalar contribution, \nabla_A (\phi_mu \phi_nu)

                    double f_0 = v2rho2[3 * g + 0] * rhow[2 * g + 0] +
                                 v2rho2[3 * g + 1] * rhow[2 * g + 1] +

                                 v2rhosigma[6 * g + 0] * grhow_grho_aa +
                                 v2rhosigma[6 * g + 1] * grhow_grho_ab +
                                 v2rhosigma[6 * g + 2] * grhow_grho_bb;

                    gatmx += w * f_0 * gdenx[atom_g];
                    gatmy += w * f_0 * gdeny[atom_g];
                    gatmz += w * f_0 * gdenz[atom_g];

                    // vector contribution, \nabla_A (\nabla (\phi_mu \phi_nu))

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

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;

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
                                                     const CXCNewFunctional& xcFunctional) const
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

        gtoeval::preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 1, _screeningThresholdForGTOValues, boxdim);  // 1st order GTO derivative

        timer.stop("GTO pre-screening");

        // GTO values on grid points

        timer.start("OMP GTO evaluation");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            gtoeval::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtovec, xcoords, ycoords, zcoords, gridblockpos,
                                             grid_batch_offset, grid_batch_size, skip_cgto_ids);
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

        dengridgen::generateDensityForLDA(rho, npoints, mat_chi, gs_sub_dens_mat, timer);

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

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto rwdengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, rwdenmat, xcfuntype, timer);

        CDensityGridQuad rwdengridquad(npoints, numdens_rw2, xcfuntype, dengrid::ab);

        rwdengridquad.DensityProd(rwdengrid, xcfuntype, numdens_rw2, quadMode);

        timer.stop("Density grid quad");

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, npoints);
        CDenseMatrix dengrady(natoms, npoints);
        CDenseMatrix dengradz(natoms, npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(gs_sub_dens_mat, mat_chi);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto naos = mat_chi.getNumberOfRows();

        auto F_val = mat_F.values();

        auto chi_x_val = mat_chi_x.values();
        auto chi_y_val = mat_chi_y.values();
        auto chi_z_val = mat_chi_z.values();

        auto gdenx = dengradx.values();
        auto gdeny = dengrady.values();
        auto gdenz = dengradz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * npoints;

                auto nu_offset = nu * npoints;

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

        xcFunctional.compute_fxc_for_lda(npoints, rho, v2rho2);

        xcFunctional.compute_kxc_for_lda(npoints, rho, v3rho3);

        timer.stop("XC functional eval.");

        // screen density grid, weights and GTO matrix

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenKxcFockForLDA(rho, v2rho2, v3rho3, npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // pointers to perturbed density gradient norms

        auto rhow1a = rwdengridquad.rhow1rhow2(0);

        // Note: rw2DensityMatrix is zero in KxcGradientForLDA
        // auto rhow12a = rw2DensityGrid.alphaDensity(iFock);
        // auto rhow12b = rw2DensityGrid.betaDensity(iFock);

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.data(thread_id);

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(local_weights, \
                        v3rho3, rhow1a, gdenx, gdeny, gdenz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac = local_weights[g] * (v3rho3[4 * g + 0] + 2.0 * v3rho3[4 * g + 1] + v3rho3[4 * g + 2]) * rhow1a[g];

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

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;

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
                                                     const CXCNewFunctional& xcFunctional) const
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
                                                 gtovec, xcoords, ycoords, zcoords, gridblockpos,
                                                 grid_batch_offset, grid_batch_size, skip_cgto_ids);
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

        dengridgen::generateDensityForGGA(rho, rhograd, sigma, npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                          gs_sub_dens_mat, timer);

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

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto rwdengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,
                                                               rwdenmat, xcfuntype, timer);

        CDensityGridQuad rwdengridquad(npoints, numdens_rw2, xcfuntype, dengrid::ab);

        rwdengridquad.DensityProd(rwdengrid, xcfuntype, numdens_rw2, quadMode);

        timer.stop("Density grid quad");

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, npoints);
        CDenseMatrix dengrady(natoms, npoints);
        CDenseMatrix dengradz(natoms, npoints);

        CDenseMatrix dengradxx(natoms, npoints);
        CDenseMatrix dengradxy(natoms, npoints);
        CDenseMatrix dengradxz(natoms, npoints);

        CDenseMatrix dengradyx(natoms, npoints);
        CDenseMatrix dengradyy(natoms, npoints);
        CDenseMatrix dengradyz(natoms, npoints);

        CDenseMatrix dengradzx(natoms, npoints);
        CDenseMatrix dengradzy(natoms, npoints);
        CDenseMatrix dengradzz(natoms, npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(gs_sub_dens_mat, mat_chi);

        auto mat_F_x = denblas::multAB(gs_sub_dens_mat, mat_chi_x);
        auto mat_F_y = denblas::multAB(gs_sub_dens_mat, mat_chi_y);
        auto mat_F_z = denblas::multAB(gs_sub_dens_mat, mat_chi_z);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto naos = mat_chi.getNumberOfRows();

        auto F_val = mat_F.values();

        auto F_x_val = mat_F_x.values();
        auto F_y_val = mat_F_y.values();
        auto F_z_val = mat_F_z.values();

        auto chi_x_val = mat_chi_x.values();
        auto chi_y_val = mat_chi_y.values();
        auto chi_z_val = mat_chi_z.values();

        auto chi_xx_val = mat_chi_xx.values();
        auto chi_xy_val = mat_chi_xy.values();
        auto chi_xz_val = mat_chi_xz.values();
        auto chi_yy_val = mat_chi_yy.values();
        auto chi_yz_val = mat_chi_yz.values();
        auto chi_zz_val = mat_chi_zz.values();

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

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * npoints;

                auto nu_offset = nu * npoints;

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

        xcFunctional.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

        xcFunctional.compute_fxc_for_gga(npoints, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

        xcFunctional.compute_kxc_for_gga(npoints, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);

        timer.stop("XC functional eval.");

        // screen density grid, weights and GTO matrix

        timer.start("Density screening");

        gridscreen::copyWeights(local_weights, gridblockpos, weights, npoints);

        gridscreen::screenKxcFockForGGA(rho, sigma, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2,
                                        v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,
                                        npoints, _screeningThresholdForDensityValues);

        timer.stop("Density screening");

        // pointers to perturbed densities

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

        // Note: rw2DensityMatrix is zero in KxcGradientForGGA
        // auto rhow12a = rw2DensityGrid.alphaDensity(iFock);
        // auto gradw12a_x = rw2DensityGrid.alphaDensityGradientX(iFock);
        // auto gradw12a_y = rw2DensityGrid.alphaDensityGradientY(iFock);
        // auto gradw12a_z = rw2DensityGrid.alphaDensityGradientZ(iFock);

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.data(thread_id);

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(local_weights, \
                        rhograd, vsigma, v2rho2, v2rhosigma, v2sigma2, \
                        v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3, \
                        rhow1rhow2, rxw1rhow2, ryw1rhow2, rzw1rhow2, \
                        rxw1rxw2, rxw1ryw2, rxw1rzw2, ryw1rxw2, ryw1ryw2, ryw1rzw2, rzw1rxw2, rzw1ryw2, rzw1rzw2, \
                        gdenx, gdeny, gdenz, gdenxx, gdenxy, gdenxz, gdenyx, gdenyy, gdenyz, gdenzx, gdenzy, gdenzz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double w = local_weights[g];

                    // double rxw12a = gradw12a_x[g];
                    // double ryw12a = gradw12a_y[g];
                    // double rzw12a = gradw12a_z[g];

                    double grada_x_g = rhograd[6 * g + 0];
                    double grada_y_g = rhograd[6 * g + 1];
                    double grada_z_g = rhograd[6 * g + 2];

                    // double l2contract = grada_x_g * rxw12a + grada_y_g * ryw12a + grada_z_g * rzw12a;
                    // double l5contract_x = grada_x_g * l2contract;
                    // double l5contract_y = grada_y_g * l2contract;
                    // double l5contract_z = grada_z_g * l2contract;

                    double q2contract = grada_x_g * rxw1rhow2[g] + grada_y_g * ryw1rhow2[g] + grada_z_g * rzw1rhow2[g];

                    double q3contract = grada_x_g * grada_x_g * rxw1rxw2[g] +
                                        grada_x_g * grada_y_g * rxw1ryw2[g] +
                                        grada_x_g * grada_z_g * rxw1rzw2[g] +
                                        grada_y_g * grada_x_g * ryw1rxw2[g] +
                                        grada_y_g * grada_y_g * ryw1ryw2[g] +
                                        grada_y_g * grada_z_g * ryw1rzw2[g] +
                                        grada_z_g * grada_x_g * rzw1rxw2[g] +
                                        grada_z_g * grada_y_g * rzw1ryw2[g] +
                                        grada_z_g * grada_z_g * rzw1rzw2[g];

                    double q4contract = rxw1rxw2[g] + ryw1ryw2[g] + rzw1rzw2[g];

                    double q7contract_x = grada_x_g * (grada_x_g * rxw1rhow2[g] + grada_y_g * ryw1rhow2[g] + grada_z_g * rzw1rhow2[g]);
                    double q7contract_y = grada_y_g * (grada_x_g * rxw1rhow2[g] + grada_y_g * ryw1rhow2[g] + grada_z_g * rzw1rhow2[g]);
                    double q7contract_z = grada_z_g * (grada_x_g * rxw1rhow2[g] + grada_y_g * ryw1rhow2[g] + grada_z_g * rzw1rhow2[g]);

                    double q8contract_x = grada_x_g * rxw1rxw2[g] + grada_y_g * rxw1ryw2[g] + grada_z_g * rxw1rzw2[g];
                    double q8contract_y = grada_x_g * ryw1rxw2[g] + grada_y_g * ryw1ryw2[g] + grada_z_g * ryw1rzw2[g];
                    double q8contract_z = grada_x_g * rzw1rxw2[g] + grada_y_g * rzw1ryw2[g] + grada_z_g * rzw1rzw2[g];

                    double q9contract_x = grada_x_g * q3contract;
                    double q9contract_y = grada_y_g * q3contract;
                    double q9contract_z = grada_z_g * q3contract;

                    double q10contract_x = grada_x_g * rxw1rxw2[g] + grada_y_g * ryw1rxw2[g] + grada_z_g * rzw1rxw2[g];
                    double q10contract_y = grada_x_g * rxw1ryw2[g] + grada_y_g * ryw1ryw2[g] + grada_z_g * rzw1ryw2[g];
                    double q10contract_z = grada_x_g * rxw1rzw2[g] + grada_y_g * ryw1rzw2[g] + grada_z_g * rzw1rzw2[g];

                    double q11contract_x = grada_x_g * rxw1rxw2[g] + grada_x_g * ryw1ryw2[g] + grada_x_g * rzw1rzw2[g];
                    double q11contract_y = grada_y_g * rxw1rxw2[g] + grada_y_g * ryw1ryw2[g] + grada_y_g * rzw1rzw2[g];
                    double q11contract_z = grada_z_g * rxw1rxw2[g] + grada_z_g * ryw1ryw2[g] + grada_z_g * rzw1rzw2[g];

                    // functional derivatives

                    // auto vsigma_a = vsigma[3 * g + 0];
                    // auto vsigma_c = vsigma[3 * g + 1];

                    // auto v2rho2_aa = v2rho2[3 * g + 0];
                    // auto v2rho2_ab = v2rho2[3 * g + 1];

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
                    auto v3sigma3_ccb = v3sigma3[10 * g + 7];
                    auto v3sigma3_cbb = v3sigma3[10 * g + 8];

                    // Scalar contribution

                    double prefac = 0.0;

                    // vxc 1 contributions

                    // L1
                    // prefac += (v2rho2_aa + v2rho2_ab) * rhow12a[g];

                    // L2
                    // prefac += (2.0*v2rhosigma_ac + 2.0*v2rhosigma_ab + 2.0*v2rhosigma_aa) * l2contract;

                    // vxc 2 contributions

                    // Q1
                    prefac += (v3rho3_aaa + 2.0 * v3rho3_aab + v3rho3_abb) * rhow1rhow2[g];

                    // Q2
                    prefac += (2.0*v3rho2sigma_abc + 2.0*v3rho2sigma_abb + 2.0*v3rho2sigma_aba 
                             + 2.0*v3rho2sigma_aac + 2.0*v3rho2sigma_aab + 2.0*v3rho2sigma_aaa) * q2contract;

                    // Q3
                    prefac += (4.0*v3rhosigma2_acc + 8.0*v3rhosigma2_acb + 4.0*v3rhosigma2_abb 
                             + 8.0*v3rhosigma2_aac + 8.0*v3rhosigma2_aab + 4.0*v3rhosigma2_aaa) * q3contract;

                    // Q4
                    prefac += (2.0*v2rhosigma_ac + 2.0*v2rhosigma_ab + 2.0*v2rhosigma_aa) * q4contract;

                    gatmx += w * prefac * gdenx[atom_g];
                    gatmy += w * prefac * gdeny[atom_g];
                    gatmz += w * prefac * gdenz[atom_g];

                    // vector contribution

                    double xcomp = 0.0, ycomp = 0.0, zcomp = 0.0;

                    // vxc 1 contributions

                    // L3
                    // double l3 = v2rhosigma_bc + 2.0*v2rhosigma_ba + v2rhosigma_ac + 2.0*v2rhosigma_aa;
                    // xcomp += l3 * grada_x_g * rhow12a[g];
                    // ycomp += l3 * grada_y_g * rhow12a[g];
                    // zcomp += l3 * grada_z_g * rhow12a[g];

                    // L4
                    // double l4 = vsigma_c + 2.0*vsigma_a;
                    // xcomp += l4 * rxw12a;
                    // ycomp += l4 * ryw12a;
                    // zcomp += l4 * rzw12a;

                    // L5
                    // double l5 = 2.0*v2sigma2_cc + 2.0*v2sigma2_cb + 6.0*v2sigma2_ac + 4.0*v2sigma2_ab + 4.0*v2sigma2_aa;
                    // xcomp += l5 * l5contract_x;
                    // ycomp += l5 * l5contract_y;
                    // zcomp += l5 * l5contract_z;

                    // vxc 2 contributions

                    // Q5
                    double q5 = v3rho2sigma_bbc + 2.0*v3rho2sigma_bba + 2.0*v3rho2sigma_abc + 4.0*v3rho2sigma_aba
                              + v3rho2sigma_aac + 2.0*v3rho2sigma_aaa;

                    xcomp += q5 * grada_x_g * rhow1rhow2[g];
                    ycomp += q5 * grada_y_g * rhow1rhow2[g];
                    zcomp += q5 * grada_z_g * rhow1rhow2[g];

                    // Q6
                    double q6 = v2rhosigma_bc + 2.0*v2rhosigma_ba + v2rhosigma_ac + 2.0*v2rhosigma_aa;

                    xcomp += q6 * rxw1rhow2[g];
                    ycomp += q6 * ryw1rhow2[g];
                    zcomp += q6 * rzw1rhow2[g];

                    // Q7
                    double q7 = 2.0*v3rhosigma2_bcc + 2.0*v3rhosigma2_bcb + 6.0*v3rhosigma2_bac + 4.0*v3rhosigma2_bab + 4.0*v3rhosigma2_baa
                              + 2.0*v3rhosigma2_acc + 2.0*v3rhosigma2_acb + 6.0*v3rhosigma2_aac + 4.0*v3rhosigma2_aab + 4.0*v3rhosigma2_aaa;

                    xcomp += q7 * q7contract_x;
                    ycomp += q7 * q7contract_y;
                    zcomp += q7 * q7contract_z;

                    // Q8
                    double q8 = 2.0*v2sigma2_cc + 2.0*v2sigma2_cb + 6.0*v2sigma2_ac + 4.0*v2sigma2_ab + 4.0*v2sigma2_aa;

                    xcomp += q8 * (q8contract_x + q10contract_x + q11contract_x);
                    ycomp += q8 * (q8contract_y + q10contract_y + q11contract_y);
                    zcomp += q8 * (q8contract_z + q10contract_z + q11contract_z);

                    // Q9
                    double q9 =  4.0*v3sigma3_ccc + 8.0*v3sigma3_ccb +  4.0*v3sigma3_cbb + 16.0*v3sigma3_acc 
                              + 24.0*v3sigma3_acb + 8.0*v3sigma3_abb + 20.0*v3sigma3_aac 
                              + 16.0*v3sigma3_aab + 8.0*v3sigma3_aaa;

                    xcomp += q9 * q9contract_x;
                    ycomp += q9 * q9contract_y;
                    zcomp += q9 * q9contract_z;

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

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;

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
