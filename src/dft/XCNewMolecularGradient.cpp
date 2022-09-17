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
#include "GtoEvaluator.hpp"
#include "XCFuncType.hpp"
#include "XCVarsType.hpp"

CXCNewMolecularGradient::CXCNewMolecularGradient(MPI_Comm comm)

    : _screeningThresholdForGTOValues(1.0e-12)
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
                                              const CMolecularGrid&   molecularGrid,
                                              const std::string&      xcFuncLabel) const
{
    return integrateVxcGradient(molecule, basis, gsDensityMatrix, gsDensityMatrix, molecularGrid, xcFuncLabel);
}

CDenseMatrix
CXCNewMolecularGradient::integrateVxcGradient(const CMolecule&        molecule,
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
                                              const CMolecularGrid&   molecularGrid,
                                              const std::string&      xcFuncLabel) const
{
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
                                              const CMolecularGrid&   molecularGrid,
                                              const std::string&      xcFuncLabel) const
{
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
            //return _integrateKxcGradientForGGA(molecule, basis, rwDensityMatrixOne, rwDensityMatrixTwo,

            //                                   gsDensityMatrix, molecularGrid, fvxc);
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

    molgrad_threads.zero();

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

        _preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 1, boxdim);  // 1st order GTO derivative

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

        auto gs_sub_dens_mat = _getSubDensityMatrix(gsDensityMatrix, 0, aoinds, aocount, naos);

        auto rw_sub_dens_mat = _getSubDensityMatrix(rwDensityMatrix, 0, aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto gsdengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, gs_sub_dens_mat, xcfuntype, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, npoints);

        CDenseMatrix dengrady(natoms, npoints);

        CDenseMatrix dengradz(natoms, npoints);

        dengradx.zero();

        dengrady.zero();

        dengradz.zero();

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

        CXCGradientGrid vxcgrid(npoints, gsdengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxcgrid, gsdengrid);

        timer.stop("XC functional eval.");

        // exchange-correlation functional derivative

        auto grhoa = vxcgrid.xcGradientValues(xcvars::rhoa);

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

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(weights, grhoa, \
                        gdenx, gdeny, gdenz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac = weights[gridblockpos + g] * grhoa[g];

                    gatmx += 2.0 * prefac * gdenx[atom_g];

                    gatmy += 2.0 * prefac * gdeny[atom_g];

                    gatmz += 2.0 * prefac * gdenz[atom_g];
                }

                gatm[iatom * 3 + 0] += gatmx;

                gatm[iatom * 3 + 1] += gatmy;

                gatm[iatom * 3 + 2] += gatmz;
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

    molgrad.zero();

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

    molgrad_threads.zero();

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

        _preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 2, boxdim);  // 2nd order GTO derivative

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

        auto gs_sub_dens_mat = _getSubDensityMatrix(gsDensityMatrix, 0, aoinds, aocount, naos);

        auto rw_sub_dens_mat = _getSubDensityMatrix(rwDensityMatrix, 0, aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto gsdengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                                               gs_sub_dens_mat, xcfuntype, timer);

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

        dengradx.zero();

        dengrady.zero();

        dengradz.zero();

        dengradxx.zero();

        dengradxy.zero();

        dengradxz.zero();

        dengradyx.zero();

        dengradyy.zero();

        dengradyz.zero();

        dengradzx.zero();

        dengradzy.zero();

        dengradzz.zero();

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

        CXCGradientGrid vxcgrid(npoints, gsdengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxcgrid, gsdengrid);

        timer.stop("XC functional eval.");

        // pointers to exchange-correlation functional derrivatives

        auto grhoa = vxcgrid.xcGradientValues(xcvars::rhoa);

        auto ggrada = vxcgrid.xcGradientValues(xcvars::grada);

        auto ggradab = vxcgrid.xcGradientValues(xcvars::gradab);

        // pointers to density gradient norms

        auto ngrada = gsdengrid.alphaDensityGradient(0);

        auto gradax = gsdengrid.alphaDensityGradientX(0);

        auto graday = gsdengrid.alphaDensityGradientY(0);

        auto gradaz = gsdengrid.alphaDensityGradientZ(0);

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

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(weights, \
                        grhoa, ggrada, ggradab, ngrada, gradax, graday, gradaz, gdenx, gdeny, gdenz, \
                        gdenxx, gdenxy, gdenxz, gdenyx, gdenyy, gdenyz, gdenzx, gdenzy, gdenzz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto prefac = weights[gridblockpos + g] * grhoa[g];

                    gatmx += 2.0 * prefac * gdenx[atom_g];

                    gatmy += 2.0 * prefac * gdeny[atom_g];

                    gatmz += 2.0 * prefac * gdenz[atom_g];

                    prefac =  weights[gridblockpos + g] * (ggrada[g] / ngrada[g] + ggradab[g]);

                    gatmx += 2.0 * prefac * (gradax[g] * gdenxx[atom_g] + graday[g] * gdenyx[atom_g] + gradaz[g] * gdenzx[atom_g]);

                    gatmy += 2.0 * prefac * (gradax[g] * gdenxy[atom_g] + graday[g] * gdenyy[atom_g] + gradaz[g] * gdenzy[atom_g]);

                    gatmz += 2.0 * prefac * (gradax[g] * gdenxz[atom_g] + graday[g] * gdenyz[atom_g] + gradaz[g] * gdenzz[atom_g]);
                }

                gatm[iatom * 3 + 0] += gatmx;

                gatm[iatom * 3 + 1] += gatmy;

                gatm[iatom * 3 + 2] += gatmz;
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

    molgrad.zero();

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

    molgrad_threads.zero();

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

        _preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 1, boxdim);  // 1st order GTO derivative

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

        auto gs_sub_dens_mat = _getSubDensityMatrix(gsDensityMatrix, 0, aoinds, aocount, naos);

        auto rw_sub_dens_mat_one = _getSubDensityMatrix(rwDensityMatrixOne, 0, aoinds, aocount, naos);

        auto rw_sub_dens_mat_two = _getSubDensityMatrix(rwDensityMatrixTwo, 0, aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto gsdengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, gs_sub_dens_mat, xcfuntype, timer);

        auto rwdengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, rw_sub_dens_mat_one, xcfuntype, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, npoints);

        CDenseMatrix dengrady(natoms, npoints);

        CDenseMatrix dengradz(natoms, npoints);

        dengradx.zero();

        dengrady.zero();

        dengradz.zero();

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

        CXCHessianGrid vxc2grid(npoints, gsdengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxc2grid, gsdengrid);

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

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.data(thread_id);

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(weights, \
                        grho_aa, grho_ab, rhowa, rhowb, gdenx, gdeny, gdenz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac = weights[gridblockpos + g] * (grho_aa[g] * rhowa[g] + grho_ab[g] * rhowb[g]);

                    gatmx += 2.0 * prefac * gdenx[atom_g];

                    gatmy += 2.0 * prefac * gdeny[atom_g];

                    gatmz += 2.0 * prefac * gdenz[atom_g];
                }

                gatm[iatom * 3 + 0] += gatmx;

                gatm[iatom * 3 + 1] += gatmy;

                gatm[iatom * 3 + 2] += gatmz;
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

    molgrad.zero();

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

    molgrad_threads.zero();

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

        _preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 2, boxdim);  // 2nd order GTO derivative

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

        auto gs_sub_dens_mat = _getSubDensityMatrix(gsDensityMatrix, 0, aoinds, aocount, naos);

        auto rw_sub_dens_mat_one = _getSubDensityMatrix(rwDensityMatrixOne, 0, aoinds, aocount, naos);

        auto rw_sub_dens_mat_two = _getSubDensityMatrix(rwDensityMatrixTwo, 0, aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto gsdengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                                               gs_sub_dens_mat, xcfuntype, timer);

        auto rwdengrid = dengridgen::generateDensityGridForGGA(npoints, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z,

                                                               rw_sub_dens_mat_one, xcfuntype, timer);

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

        dengradx.zero();

        dengrady.zero();

        dengradz.zero();

        dengradxx.zero();

        dengradxy.zero();

        dengradxz.zero();

        dengradyx.zero();

        dengradyy.zero();

        dengradyz.zero();

        dengradzx.zero();

        dengradzy.zero();

        dengradzz.zero();

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

        CXCGradientGrid vxcgrid(npoints, gsdengrid.getDensityGridType(), xcfuntype);

        CXCHessianGrid vxc2grid(npoints, gsdengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxcgrid, gsdengrid);

        xcFunctional.compute(vxc2grid, gsdengrid);

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

        auto ngrada = gsdengrid.alphaDensityGradient(0);

        auto ngradb = gsdengrid.betaDensityGradient(0);

        auto grada_x = gsdengrid.alphaDensityGradientX(0);

        auto grada_y = gsdengrid.alphaDensityGradientY(0);

        auto grada_z = gsdengrid.alphaDensityGradientZ(0);

        auto gradb_x = gsdengrid.betaDensityGradientX(0);

        auto gradb_y = gsdengrid.betaDensityGradientY(0);

        auto gradb_z = gsdengrid.betaDensityGradientZ(0);

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

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.data(thread_id);

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(weights, \
                        ggrad_a, ggrad_c, grho_aa, grho_ab, gmix_aa, gmix_ab, gmix_ac, gmix_bc, \
                        ggrad_aa, ggrad_ab, ggrad_ac, ggrad_bc, ggrad_cc, ngrada, ngradb, \
                        grada_x, grada_y, grada_z, gradb_x, gradb_y, gradb_z, rhowa, rhowb, \
                        gradwa_x, gradwa_y, gradwa_z, gradwb_x, gradwb_y, gradwb_z, gdenx, gdeny, gdenz, \
                        gdenxx, gdenxy, gdenxz, gdenyx, gdenyy, gdenyz, gdenzx, gdenzy, gdenzz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

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

                    double zetac = grada_x[g] * rxwb + grada_y[g] * rywb

                                   + grada_z[g] * rzwb + gradb_x[g] * rxwa

                                   + gradb_y[g] * rywa + gradb_z[g] * rzwa;

                    // contribution from \nabla_A (\phi_mu \phi_nu)

                    // first contribution

                    double fac0 = gmix_aa[g] * zetaa + gmix_ab[g] * zetab

                                  + gmix_ac[g] * zetac + grho_aa[g] * rhowa[g]

                                  + grho_ab[g] * rhowb[g];

                    // w * a0 * fac0;   a0 = bgaos[m] * kgaos[m];

                    double prefac = w * fac0;

                    gatmx += prefac * gdenx[atom_g];

                    gatmy += prefac * gdeny[atom_g];

                    gatmz += prefac * gdenz[atom_g];

                    // contribution from \nabla_A (\nabla (\phi_mu \phi_nu))

                    double xcomp = 0.0;

                    double ycomp = 0.0;

                    double zcomp = 0.0;

                    // second contribution

                    double facr = gmix_aa[g] * rhowa[g] + gmix_ab[g] * rhowb[g]

                                  + ggrad_aa[g] * zetaa + ggrad_ab[g] * zetab + ggrad_ac[g] * zetac;

                    // w * facr * ar;   ar = ax * rxa + ay * rya + az * rza;

                    prefac = w * facr;

                    xcomp += prefac * rxa;

                    ycomp += prefac * rya;

                    zcomp += prefac * rza;

                    // third contribution

                    double facz = gmix_ac[g] * rhowa[g] + gmix_bc[g] * rhowb[g]

                                  + ggrad_ac[g] * zetaa + ggrad_bc[g] * zetab + ggrad_cc[g] * zetac;

                    // w * facz * arb;   arb = ax * grada_x[g] + ay * grada_y[g] + az * grada_z[g];

                    prefac = w * facz;

                    xcomp += prefac * grada_x[g];

                    ycomp += prefac * grada_y[g];

                    zcomp += prefac * grada_z[g];

                    // fourth contribution

                    // w * znva * ggrad_a[g] * ab;
                    // ab = ax * rxwa + ay * rywa + az * rzwa - ar * zetaa;
                    // ar = ax * rxa + ay * rya + az * rza;

                    prefac = w * znva * ggrad_a[g];

                    xcomp += prefac * (rxwa - rxa * zetaa);

                    ycomp += prefac * (rywa - rya * zetaa);

                    zcomp += prefac * (rzwa - rza * zetaa);

                    // fifth contribution

                    // w * ggrad_c[g] * abw;
                    // abw = ax * rxwa + ay * rywa + az * rzwa;

                    prefac = w * ggrad_c[g];

                    xcomp += prefac * rxwa;

                    ycomp += prefac * rywa;

                    zcomp += prefac * rzwa;

                    gatmx += (xcomp * gdenxx[atom_g] + ycomp * gdenyx[atom_g] + zcomp * gdenzx[atom_g]);

                    gatmy += (xcomp * gdenxy[atom_g] + ycomp * gdenyy[atom_g] + zcomp * gdenzy[atom_g]);

                    gatmz += (xcomp * gdenxz[atom_g] + ycomp * gdenyz[atom_g] + zcomp * gdenzz[atom_g]);
                }

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

    molgrad.zero();

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

    molgrad_threads.zero();

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

        _preScreenGtos(skip_cgto_ids, skip_ao_ids, gtovec, 1, boxdim);  // 1st order GTO derivative

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

        auto gs_sub_dens_mat = _getSubDensityMatrix(gsDensityMatrix, 0, aoinds, aocount, naos);

        auto rw_sub_dens_mat_one = _getSubDensityMatrix(rwDensityMatrixOne, 0, aoinds, aocount, naos);

        auto rw_sub_dens_mat_two = _getSubDensityMatrix(rwDensityMatrixTwo, 0, aoinds, aocount, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        auto xcfuntype = xcFunctional.getFunctionalType();

        auto gsdengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, gs_sub_dens_mat, xcfuntype, timer);

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

        auto rwdengrid = dengridgen::generateDensityGridForLDA(npoints, mat_chi, rwdenmat, xcfuntype, timer);

        CDensityGridQuad rwdengridquad(npoints, numdens_rw2, xcfuntype, dengrid::ab);

        rwdengridquad.DensityProd(rwdengrid, xcfuntype, numdens_rw2, quadMode);

        timer.stop("Density grid quad");

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, npoints);

        CDenseMatrix dengrady(natoms, npoints);

        CDenseMatrix dengradz(natoms, npoints);

        dengradx.zero();

        dengrady.zero();

        dengradz.zero();

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

        CXCCubicHessianGrid vxc3grid(npoints, gsdengrid.getDensityGridType(), xcfuntype);

        xcFunctional.compute(vxc3grid, gsdengrid);

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

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.data(thread_id);

            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(weights, \
                        df3000, df2100, df1200, rhow1a, gdenx, gdeny, gdenz : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac = weights[gridblockpos + g] * (df3000[g] + df2100[g] + df2100[g] + df1200[g]) * rhow1a[g];

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

    molgrad.zero();

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

std::array<double, 6>
CXCNewMolecularGradient::_getGridBoxDimension(const int32_t gridBlockPosition,
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
CXCNewMolecularGradient::_preScreenGtos(CMemBlock<int32_t>&          skipCgtoIds,
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

                // 0th-order derivative
                // gto:                    r^{ang}   |C| exp(-alpha r^2)

                // 1st-order derivative
                // gto_m:              ang r^{ang-1} |C| exp(-alpha r^2)
                // gto_p:           2alpha r^{ang+1} |C| exp(-alpha r^2)

                // 2nd-order derivative
                // gto_m2:     ang (ang-1) r^{ang-2} |C| exp(-alpha r^2)
                // gto   : 2alpha (2ang+1) r^{ang}   |C| exp(-alpha r^2)
                // gto_p2:        4alpha^2 r^{ang+2} |C| exp(-alpha r^2)

                auto gtolimit = maxcoef * std::exp(-minexp * r2);

                auto r = std::sqrt(r2);

                for (int32_t ipow = 0; ipow < bang; ipow++) gtolimit *= r;

                if (gtoDeriv > 0)
                {
                    gtolimit = std::max(gtolimit, 2.0 * maxexp * r * gtolimit);

                    if (bang > 0) gtolimit = std::max(gtolimit, gtolimit / r * bang);
                }

                if (gtoDeriv > 1)
                {
                    gtolimit = std::max({gtolimit, 4.0 * maxexp * maxexp * r2 * gtolimit,

                                         2.0 * maxexp * (2 * bang + 1) * gtolimit});

                    if (bang > 1) gtolimit = std::max(gtolimit, gtolimit / r2 * bang * (bang - 1));
                }

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

CDenseMatrix
CXCNewMolecularGradient::_getSubDensityMatrix(const CAODensityMatrix&     densityMatrix,
                                              const int32_t               densityIndex,
                                              const std::vector<int32_t>& aoIndices,
                                              const int32_t               aoCount,
                                              const int32_t               nAOs) const
{
    if (aoCount <= nAOs)
    {
        CDenseMatrix sub_dens(aoCount, aoCount);

        const CDenseMatrix& dens = densityMatrix.getReferenceToDensity(densityIndex);

        for (int32_t i = 0; i < aoCount; i++)
        {
            auto sub_dens_row = sub_dens.row(i);

            auto dens_row = dens.row(aoIndices[i]);

            for (int32_t j = 0; j < aoCount; j++)
            {
                sub_dens_row[j] = dens_row[aoIndices[j]];
            }
        }

        return sub_dens;
    }

    return CDenseMatrix();
}
