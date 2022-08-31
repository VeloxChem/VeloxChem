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
#include <sstream>

#include "DenseLinearAlgebra.hpp"
#include "FunctionalParser.hpp"
#include "GtoFunc.hpp"
#include "XCFuncType.hpp"
#include "XCVarsType.hpp"

CXCNewIntegrator::CXCNewIntegrator(MPI_Comm comm)

    : _numberOfPointsThreshold(800)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;
}

CXCNewIntegrator::~CXCNewIntegrator()
{
}

void
CXCNewIntegrator::setNumberOfPointsThreshold(const int32_t thresh)
{
    _numberOfPointsThreshold = thresh;
}

int32_t
CXCNewIntegrator::getNumberOfBoxes() const
{
    return static_cast<int32_t>(_boxes.size());
}

void
CXCNewIntegrator::initializeGrid(const CMolecularGrid& molgrid)
{
    auto boxdim = molgrid.getSpatialExtent();

    boxdim[0] -= 0.0001;

    boxdim[1] -= 0.0001;

    boxdim[2] -= 0.0001;

    boxdim[3] += 0.0001;

    boxdim[4] += 0.0001;

    boxdim[5] += 0.0001;

    _boxes.clear();

    _boxes.push_back(CGridBox(boxdim, molgrid.getGridPoints()));
}

void
CXCNewIntegrator::partitionGrid()
{
   auto box = _boxes.begin();

   while (box != _boxes.end())
   {
        auto npoints = box->getNumberOfGridPoints();

        if (npoints <= _numberOfPointsThreshold)
        {
            ++box;
        }
        else
        {
            // divide the box

            // box info

            auto boxdim = box->getBoxDimension();

            auto xmin = boxdim[0];

            auto ymin = boxdim[1];

            auto zmin = boxdim[2];

            auto xmax = boxdim[3];

            auto ymax = boxdim[4];

            auto zmax = boxdim[5];

            // grid points info

            auto xcoords = box->getCoordinatesX();

            auto ycoords = box->getCoordinatesY();

            auto zcoords = box->getCoordinatesZ();

            auto weights = box->getWeights();

            // find the center for dividing the box

            double xhalf = 0.0;

            double yhalf = 0.0;

            double zhalf = 0.0;

            for (int32_t g = 0; g < npoints; g++)
            {
                xhalf += xcoords[g];

                yhalf += ycoords[g];

                zhalf += zcoords[g];
            }

            xhalf /= (double)npoints;

            yhalf /= (double)npoints;

            zhalf /= (double)npoints;

            // sub boxes and grid points

            std::vector<std::array<double, 6>> subboxdims;

            std::vector<CMemBlock2D<double>> subgridpoints;

            std::vector<int32_t> subnumpoints;

            for (int32_t xval = 0; xval < 2; xval++)
            {
                std::array<double, 2> xrange = (!xval) ? std::array<double, 2>({xmin, xhalf}) :
                                                         std::array<double, 2>({xhalf, xmax});

                for (int32_t yval = 0; yval < 2; yval++)
                {
                    std::array<double, 2> yrange = (!yval) ? std::array<double, 2>({ymin, yhalf}) :
                                                             std::array<double, 2>({yhalf, ymax});

                    for (int32_t zval = 0; zval < 2; zval++)
                    {
                        std::array<double, 2> zrange = (!zval) ? std::array<double, 2>({zmin, zhalf}) :
                                                                 std::array<double, 2>({zhalf, zmax});

                        subboxdims.push_back(std::array<double, 6>(
                            {xrange[0], yrange[0], zrange[0], xrange[1], yrange[1], zrange[1]}));

                        subgridpoints.push_back(CMemBlock2D<double>(npoints, 4));

                        subnumpoints.push_back(0);
                    }
                }
            }

            for (int32_t g = 0; g < npoints; g++)
            {
                int32_t xval = (xcoords[g] < xhalf) ? 0 : 1;

                int32_t yval = (ycoords[g] < yhalf) ? 0 : 1;

                int32_t zval = (zcoords[g] < zhalf) ? 0 : 1;

                // x y z  id
                // ---------
                // 0 0 0   0
                // 0 0 1   1
                // 0 1 0   2
                // 0 1 1   3
                // 1 0 0   4
                // 1 0 1   5
                // 1 1 0   6
                // 1 1 1   7

                int32_t box_id = (xval << 2) | (yval << 1) | (zval << 0);

                auto count = subnumpoints[box_id];

                subgridpoints[box_id].data(0)[count] = xcoords[g];

                subgridpoints[box_id].data(1)[count] = ycoords[g];

                subgridpoints[box_id].data(2)[count] = zcoords[g];

                subgridpoints[box_id].data(3)[count] = weights[g];

                subnumpoints[box_id]++;
            }

            for (int32_t box_id = 0; box_id < static_cast<int32_t>(subboxdims.size()); box_id++)
            {
                auto count = subnumpoints[box_id];

                if (count > 0)
                {
                    _boxes.push_back(CGridBox(subboxdims[box_id], subgridpoints[box_id].slice(0, count)));
                }
            }

            box = _boxes.erase(box);
        }
    }
}

std::string
CXCNewIntegrator::getGridInformation() const
{
    std::stringstream ss;

    ss << "Number of grid boxes: " << getNumberOfBoxes() << "\n";

    int32_t boxind = 0;

    for (auto box = _boxes.cbegin(); box != _boxes.cend(); ++box)
    {
        auto boxdim = box->getBoxDimension();

        auto npoints = box->getNumberOfGridPoints();

        ss << "  Grid boxdim " << boxind << ", number of points: " << npoints << ", xyz: ";

        ss << "(" << boxdim[0] << ", " << boxdim[1] << ", " << boxdim[2] << "), ";

        ss << "(" << boxdim[3] << ", " << boxdim[4] << ", " << boxdim[5] << ")\n";

        ++boxind;
    }

    return ss.str();
}

CDenseMatrix
CXCNewIntegrator::integrateVxcFock(const CMolecule&        molecule,
                                   const CMolecularBasis&  basis,
                                   const CAODensityMatrix& densityMatrix,
                                   const std::string&      xcFuncLabel) const
{
    // create GTOs container

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    // number of AOs

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    CDenseMatrix mat_Vxc(naos, naos);

    mat_Vxc.zero();

    for (auto box = _boxes.cbegin(); box != _boxes.cend(); ++box)
    {
        // grid points in box

        const CMemBlock2D<double>& points = box->getGridPoints();

        // generate reference density grid and compute exchange-correlation functional derivative

        auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

        auto xcfuntype = fvxc.getFunctionalType();

        CDensityGridDriver dgdrv(_locComm);

        CMolecularGrid mgrid(points);

        auto dengrid = dgdrv.generate(densityMatrix, molecule, basis, mgrid, xcfuntype);

        CXCGradientGrid vxcgrid(mgrid.getNumberOfGridPoints(), dengrid.getDensityGridType(), xcfuntype);

        fvxc.compute(vxcgrid, dengrid);

        // compute partial contribution to Vxc matrix

        CDenseMatrix partial_mat_Vxc;

        if (xcfuntype == xcfun::lda)
        {
            partial_mat_Vxc = _integratePartialVxcFockForLDA(points, gtovec, vxcgrid);
        }
        else if (xcfuntype == xcfun::gga)
        {
            //partial_mat_Vxc = _integratePartialVxcFockForGGA(points, gtovec, vxcgrid);
        }

        mat_Vxc = denblas::addAB(mat_Vxc, partial_mat_Vxc, 1.0);
    }

    // destroy GTOs container

    delete gtovec;

    return mat_Vxc;
}

CDenseMatrix
CXCNewIntegrator::_integratePartialVxcFockForLDA(const CMemBlock2D<double>& points,
                                                 const CGtoContainer*       gtoContainer,
                                                 const CXCGradientGrid&     xcGradientGrid) const
{
    // grid points info

    auto npoints = points.size(0);

    auto xcoords = points.data(0);

    auto ycoords = points.data(1);

    auto zcoords = points.data(2);

    auto weights = points.data(3);

    // GTO values on grid points

    auto naos = gtoContainer->getNumberOfAtomicOrbitals();

    CMemBlock2D<double> gaos(npoints, naos);

    gaos.zero();

    gtorec::computeGtosValuesForLDA(gaos, gtoContainer, xcoords, ycoords, zcoords, 0, 0, npoints);

    CDenseMatrix mat_chi(naos, npoints);

    for (int32_t nu = 0; nu < naos; nu++)
    {
        std::memcpy(mat_chi.row(nu), gaos.data(nu), npoints * sizeof(double));
    }

    // exchange-correlation functional derivative

    auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);

    // eq.(30), JCTC 2021, 17, 1512-1521

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

    // eq.(31), JCTC 2021, 17, 1512-1521

    return denblas::multABt(mat_chi, mat_G);
}
