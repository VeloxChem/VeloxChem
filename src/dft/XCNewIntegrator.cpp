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

    : _boxes(std::vector<std::array<double, 6>>())

    , _pointsInBoxes(std::vector<CMemBlock2D<double>>())

    , _numberOfPointsThreshold(800)

    , _maximumLevel(100)
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

void
CXCNewIntegrator::setMaximumLevel(const int32_t maxlevel)
{
    _maximumLevel = maxlevel;
}

int32_t
CXCNewIntegrator::getNumberOfBoxes() const
{
    return static_cast<int32_t>(_boxes.size());
}

void
CXCNewIntegrator::initializeGrid(const CMolecularGrid& molgrid)
{
    auto box = molgrid.getSpatialExtent();

    double xmin = box[0] - 0.0001;

    double ymin = box[1] - 0.0001;

    double zmin = box[2] - 0.0001;

    double xmax = box[3] + 0.0001;

    double ymax = box[4] + 0.0001;

    double zmax = box[5] + 0.0001;

    double xhalf = 0.5 * (xmax + xmin);

    double yhalf = 0.5 * (ymax + ymin);

    double zhalf = 0.5 * (zmax + zmin);

    double maxlen = std::max({xmax - xmin, ymax - ymin, zmax - zmin});

    _boxes.clear();

    _boxes.push_back(std::array<double, 6>({xhalf - 0.5 * maxlen, yhalf - 0.5 * maxlen, zhalf - 0.5 * maxlen,
                                            xhalf + 0.5 * maxlen, yhalf + 0.5 * maxlen, zhalf + 0.5 * maxlen}));

    _pointsInBoxes.clear();

    _pointsInBoxes.push_back(molgrid.getGridPoints());
}

void
CXCNewIntegrator::partitionGrid()
{
    // go through the levels

    for (int32_t level = 0; level < _maximumLevel; level++)
    {
        // check number of grid points in the boxes

        bool done = true;

        for (int32_t i = 0; i < getNumberOfBoxes(); i++)
        {
            if (_pointsInBoxes[i].size(0) > _numberOfPointsThreshold)
            {
                done = false;

                break;
            }
        }

        if (done) break;

        // need to further divide the boxes

        std::vector<std::array<double ,6>> new_boxes;

        std::vector<CMemBlock2D<double>> new_points_in_boxes;

        for (int32_t i = 0; i < getNumberOfBoxes(); i++)
        {
            const std::array<double, 6>& box = _boxes[i];

            const CMemBlock2D<double>& points = _pointsInBoxes[i];

            auto npoints = points.size(0);

            if (npoints > _numberOfPointsThreshold)
            {
                // divide the box

                // box info

                auto xmin = box[0];

                auto ymin = box[1];

                auto zmin = box[2];

                auto xmax = box[3];

                auto ymax = box[4];

                auto zmax = box[5];

                double xhalf = 0.5 * (xmax + xmin);

                double yhalf = 0.5 * (ymax + ymin);

                double zhalf = 0.5 * (zmax + zmin);

                // grid points info

                auto xcoords = points.data(0);

                auto ycoords = points.data(1);

                auto zcoords = points.data(2);

                auto weights = points.data(3);

                // sub boxes and grid points

                std::vector<std::array<double, 6>> subboxes;

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

                            subboxes.push_back(std::array<double, 6>({xrange[0], yrange[0], zrange[0],
                                                                      xrange[1], yrange[1], zrange[1]}));

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

                for (int32_t box_id = 0; box_id < static_cast<int32_t>(subboxes.size()); box_id++)
                {
                    auto count = subnumpoints[box_id];

                    if (count > 0)
                    {
                        new_boxes.push_back(subboxes[box_id]);

                        new_points_in_boxes.push_back(subgridpoints[box_id].slice(0, count));
                    }
                }
            }
            else
            {
                // keep the box

                new_boxes.push_back(box);

                new_points_in_boxes.push_back(points);
            }
        }

        _boxes = new_boxes;

        _pointsInBoxes = new_points_in_boxes;
    }
}

std::string
CXCNewIntegrator::getGridInformation() const
{
    std::stringstream ss;

    ss << "Number of grid boxes: " << getNumberOfBoxes() << "\n";

    for (int32_t i = 0; i < getNumberOfBoxes(); i++)
    {
        const std::array<double, 6>& box = _boxes[i];

        auto npoints = _pointsInBoxes[i].size(0);

        ss << "  Grid box " << i << ", number of points: " << npoints << ", xyz: ";

        ss << "(" << box[0] << ", " << box[1] << ", " << box[2] << "), ";

        ss << "(" << box[3] << ", " << box[4] << ", " << box[5] << ")\n";
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

    for (int32_t i = 0; i < getNumberOfBoxes(); i++)
    {
        // grid points in box

        const CMemBlock2D<double>& points = _pointsInBoxes[i];

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
