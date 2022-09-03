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
#include <iomanip>
#include <sstream>
#include <iostream>

#include <omp.h>

#include "DenseLinearAlgebra.hpp"
#include "DensityGridType.hpp"
#include "FunctionalParser.hpp"
#include "GtoFunc.hpp"
#include "XCFuncType.hpp"
#include "XCVarsType.hpp"

CXCNewIntegrator::CXCNewIntegrator(MPI_Comm comm)
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
        auto count = counts.data(box_id)[0];

        auto displ = displacements.data(box_id)[0];

        timer.start("GTO evaluation");

        // grid points in box

        auto npoints = count;

        auto xcoords = molecularGrid.getCoordinatesX() + displ;

        auto ycoords = molecularGrid.getCoordinatesY() + displ;

        auto zcoords = molecularGrid.getCoordinatesZ() + displ;

        auto weights = molecularGrid.getWeights() + displ;

        // GTO values on grid points

        CMemBlock2D<double> gaos(npoints, naos);

        #pragma omp parallel
        {
            auto nthreads = omp_get_max_threads();

            auto thread_id = omp_get_thread_num();

            auto b_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto b_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            CMemBlock2D<double> local_gaos(b_size, naos);

            local_gaos.zero();

            gtorec::computeGtosValuesForLDA(local_gaos, gtovec, xcoords, ycoords, zcoords, b_offset, 0, b_size);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                std::memcpy(gaos.data(nu) + b_offset, local_gaos.data(nu), b_size * sizeof(double));
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
                if (std::fabs(gaos_nu[g]) > 1.0e-12)
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

        CAODensityMatrix sub_dens_mat({sub_dens}, denmat::rest);

        timer.stop("Density matrix");

        // generate density grid

        auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

        auto xcfuntype = fvxc.getFunctionalType();

        auto dengrid = _generateDensityGrid(npoints, mat_chi, sub_dens_mat, xcfuntype, timer);

        timer.start("XC functional");

        // compute exchange-correlation functional derivative

        CXCGradientGrid vxcgrid(npoints, dengrid.getDensityGridType(), xcfuntype);

        fvxc.compute(vxcgrid, dengrid);

        timer.stop("XC functional");

        // compute partial contribution to Vxc matrix

        CDenseMatrix partial_mat_Vxc;

        if (xcfuntype == xcfun::lda)
        {
            partial_mat_Vxc = _integratePartialVxcFockForLDA(npoints, xcoords, ycoords, zcoords, weights, mat_chi, vxcgrid, timer);
        }
        else if (xcfuntype == xcfun::gga)
        {
            //partial_mat_Vxc = _integratePartialVxcFockForGGA(npoints, xcoords, ycoords, zcoords, weights, mat_chi, vxcgrid);
        }

        timer.start("Vxc matrix dist.");

        for (int32_t row = 0; row < partial_mat_Vxc.getNumberOfRows(); row++)
        {
            auto row_orig = aoinds[row];

            for (int32_t col = 0; col < partial_mat_Vxc.getNumberOfColumns(); col++)
            {
                auto col_orig = aoinds[col];

                mat_Vxc.getMatrix(0)[row_orig * naos + col_orig] += partial_mat_Vxc.row(row)[col];
            }
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

CDensityGrid
CXCNewIntegrator::_generateDensityGrid(const int32_t           npoints,
                                       const CDenseMatrix&     gtoValuesOnGridPoints,
                                       const CAODensityMatrix& densityMatrix,
                                       const xcfun             xcFunType,
                                       CMultiTimer&            timer) const
{
    CDensityGrid dengrid(npoints, densityMatrix.getNumberOfDensityMatrices(), xcFunType, dengrid::ab);

    dengrid.zero();

    auto rhoa = dengrid.alphaDensity(0);

    auto rhob = dengrid.betaDensity(0);

    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    const CDenseMatrix& mat_chi = gtoValuesOnGridPoints;

    auto mat_F = denblas::multAB(densityMatrix.getReferenceToDensity(0), mat_chi);

    timer.stop("Density grid matmul");

    // eq.(27), JCTC 2021, 17, 1512-1521

    timer.start("Density grid rho");

    auto naos = mat_chi.getNumberOfRows();

    for (int32_t nu = 0; nu < naos; nu++)
    {
        auto F_nu = mat_F.row(nu);

        auto chi_nu = mat_chi.row(nu);

        for (int32_t g = 0; g < npoints; g++)
        {
            rhoa[g] += F_nu[g] * chi_nu[g];
        }
    }

    std::memcpy(rhob, rhoa, npoints * sizeof(double));

    timer.stop("Density grid rho");

    return dengrid;
}

CDenseMatrix
CXCNewIntegrator::_integratePartialVxcFockForLDA(const int32_t          npoints,
                                                 const double*          xcoords,
                                                 const double*          ycoords,
                                                 const double*          zcoords,
                                                 const double*          weights,
                                                 const CDenseMatrix&    gtoValuesOnGridPoints,
                                                 const CXCGradientGrid& xcGradientGrid,
                                                 CMultiTimer&           timer) const
{
    // GTO values on grid points

    const CDenseMatrix& mat_chi = gtoValuesOnGridPoints;

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
