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

#include "DenseLinearAlgebra.hpp"
#include "DensityGridType.hpp"
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

        if (npoints > _numberOfPointsThreshold)
        {
            auto newboxes = divideBoxIntoTwo(*box);

            _boxes.insert(_boxes.end(), newboxes.begin(), newboxes.end());

            box = _boxes.erase(box);
        }
        else
        {
            ++box;
        }
    }
}

std::list<CGridBox>
CXCNewIntegrator::divideBoxIntoEight(const CGridBox& box) const
{
    // box info

    auto boxdim = box.getBoxDimension();

    auto xmin = boxdim[0];

    auto ymin = boxdim[1];

    auto zmin = boxdim[2];

    auto xmax = boxdim[3];

    auto ymax = boxdim[4];

    auto zmax = boxdim[5];

    // grid points info

    auto npoints = box.getNumberOfGridPoints();

    auto xcoords = box.getCoordinatesX();

    auto ycoords = box.getCoordinatesY();

    auto zcoords = box.getCoordinatesZ();

    auto weights = box.getWeights();

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

    std::vector<std::array<double, 2>> xrange({std::array<double, 2>({xmin, xhalf}), std::array<double, 2>({xhalf, xmax})});

    std::vector<std::array<double, 2>> yrange({std::array<double, 2>({ymin, yhalf}), std::array<double, 2>({yhalf, ymax})});

    std::vector<std::array<double, 2>> zrange({std::array<double, 2>({zmin, zhalf}), std::array<double, 2>({zhalf, zmax})});

    for (const auto& x : xrange)
    {
        for (const auto& y : yrange)
        {
            for (const auto& z : zrange)
            {
                subboxdims.push_back(std::array<double, 6>({x[0], y[0], z[0], x[1], y[1], z[1]}));

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

    std::list<CGridBox> newboxes;

    for (int32_t box_id = 0; box_id < static_cast<int32_t>(subboxdims.size()); box_id++)
    {
        auto count = subnumpoints[box_id];

        if (count > 0) newboxes.push_back(CGridBox(subboxdims[box_id], subgridpoints[box_id].slice(0, count)));
    }

    return newboxes;
}

std::list<CGridBox>
CXCNewIntegrator::divideBoxIntoTwo(const CGridBox& box) const
{
    // box info

    auto boxdim = box.getBoxDimension();

    auto xmin = boxdim[0];

    auto ymin = boxdim[1];

    auto zmin = boxdim[2];

    auto xmax = boxdim[3];

    auto ymax = boxdim[4];

    auto zmax = boxdim[5];

    auto xlen = xmax - xmin;

    auto ylen = ymax - ymin;

    auto zlen = zmax - zmin;

    // the dimension that will be divided

    int32_t icart = 0;

    if ((ylen >= xlen) && (ylen >= zlen)) icart = 1;

    else if ((zlen >= xlen) && (zlen >= ylen)) icart = 2;

    // grid points info

    auto npoints = box.getNumberOfGridPoints();

    auto xcoords = box.getCoordinatesX();

    auto ycoords = box.getCoordinatesY();

    auto zcoords = box.getCoordinatesZ();

    auto weights = box.getWeights();

    std::vector<const double*> coords({xcoords, ycoords, zcoords});

    // find the center for dividing the box

    double center = 0.0;

    for (int32_t g = 0; g < npoints; g++)
    {
        center += coords[icart][g];
    }

    center /= (double)npoints;

    // sub boxes and grid points

    std::vector<std::array<double, 6>> subboxdims;

    std::vector<CMemBlock2D<double>> subgridpoints;

    std::vector<int32_t> subnumpoints;

    for (int32_t box_id = 0; box_id < 2; box_id++)
    {
        subboxdims.push_back(std::array<double, 6>(boxdim));

        subgridpoints.push_back(CMemBlock2D<double>(npoints, 4));

        subnumpoints.push_back(0);
    }

    // update subboxdims

    // icart == 0
    // std::array<double, 6>({xmin, ymin, zmin, center, ymax, zmax})
    // std::array<double, 6>({center, ymin, zmin, xmax, ymax, zmax})

    // icart == 1
    // std::array<double, 6>({xmin, ymin, zmin, xmax, center, zmax})
    // std::array<double, 6>({xmin, center, zmin, xmax, ymax, zmax})

    // icart == 2
    // std::array<double, 6>({xmin, ymin, zmin, xmax, ymax, center})
    // std::array<double, 6>({xmin, ymin, center, xmax, ymax, zmax})

    subboxdims[0][3 + icart] = center;

    subboxdims[1][icart] = center;

    // update subgridpoints

    for (int32_t g = 0; g < npoints; g++)
    {
        int32_t box_id = (coords[icart][g] < center) ? 0 : 1;

        auto count = subnumpoints[box_id];

        subgridpoints[box_id].data(0)[count] = xcoords[g];

        subgridpoints[box_id].data(1)[count] = ycoords[g];

        subgridpoints[box_id].data(2)[count] = zcoords[g];

        subgridpoints[box_id].data(3)[count] = weights[g];

        ++subnumpoints[box_id];
    }

    std::list<CGridBox> newboxes;

    for (int32_t box_id = 0; box_id < static_cast<int32_t>(subboxdims.size()); box_id++)
    {
        auto count = subnumpoints[box_id];

        if (count > 0) newboxes.push_back(CGridBox(subboxdims[box_id], subgridpoints[box_id].slice(0, count)));
    }

    return newboxes;
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

        ss << "  Grid box " << boxind << ", number of points: " << npoints << ", xyz: ";

        ss << "(" << boxdim[0] << ", " << boxdim[1] << ", " << boxdim[2] << "), ";

        ss << "(" << boxdim[3] << ", " << boxdim[4] << ", " << boxdim[5] << ")\n";

        ++boxind;
    }

    return ss.str();
}

std::string
CXCNewIntegrator::getGridStatistics() const
{
    std::stringstream ss;

    ss << "Threshold for number of points per box: " << _numberOfPointsThreshold << "\n";

    ss << "Total number of boxes: " << getNumberOfBoxes() << "\n";

    int32_t npoints_max = -1, npoints_min = -1, npoints_sum = 0;

    int32_t nboxes_dense = 0, nboxes_sparse = 0;

    int32_t npoints_bin = 100;

    std::vector<int32_t> nboxes_count(_numberOfPointsThreshold / npoints_bin + 1, 0);

    for (auto box = _boxes.cbegin(); box != _boxes.cend(); ++box)
    {
        auto npoints = box->getNumberOfGridPoints();

        if ((npoints_max == -1) || (npoints > npoints_max)) npoints_max = npoints;

        if ((npoints_min == -1) || (npoints < npoints_min)) npoints_min = npoints;

        npoints_sum += npoints;

        if (npoints > 0.8 * _numberOfPointsThreshold) ++nboxes_dense;

        if (npoints < 0.1 * _numberOfPointsThreshold) ++nboxes_sparse;

        ++nboxes_count[npoints / npoints_bin];
    }

    ss << "Maximum number of points per box: " << npoints_max << "\n";

    ss << "Minimum number of points per box: " << npoints_min << "\n";

    ss << "Total number of points in all boxes: " << npoints_sum << "\n";

    ss << "Average number of points per box: " << npoints_sum / getNumberOfBoxes() << "\n";

    ss << nboxes_dense << " boxes with npoints > 80\% of threshold " << "\n";

    ss << nboxes_sparse << " boxes with npoints < 10\% of threshold " << "\n";

    ss << "-----------------\n";

    ss << " NPoints  NBoxes\n";

    ss << "-----------------\n";

    for (auto i = 0; i < static_cast<int32_t>(nboxes_count.size()); i++)
    {
        ss << std::right << std::setfill(' ') << std::setw(4) << i * npoints_bin << "-";

        ss << std::left << std::setfill(' ') << std::setw(4) << (i + 1) * npoints_bin - 1 << ":";

        ss << std::right << std::setfill(' ') << std::setw(6) << nboxes_count[i] << "\n";
    }

    ss << "-----------------\n";

    return ss.str();
}

CAOKohnShamMatrix
CXCNewIntegrator::integrateVxcFock(const CMolecule&        molecule,
                                   const CMolecularBasis&  basis,
                                   const CAODensityMatrix& densityMatrix,
                                   const std::string&      xcFuncLabel) const
{
    // create GTOs container

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    // number of AOs

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    CAOKohnShamMatrix mat_Vxc(densityMatrix.getNumberOfRows(0), densityMatrix.getNumberOfColumns(0), true);

    mat_Vxc.zero();

    double nele = 0.0, xcene = 0.0;

    for (auto box = _boxes.cbegin(); box != _boxes.cend(); ++box)
    {
        // grid points in box

        auto npoints = box->getNumberOfGridPoints();

        auto xcoords = box->getCoordinatesX();

        auto ycoords = box->getCoordinatesY();

        auto zcoords = box->getCoordinatesZ();

        auto weights = box->getWeights();

        // GTO values on grid points

        CMemBlock2D<double> gaos(npoints, naos);

        gaos.zero();

        gtorec::computeGtosValuesForLDA(gaos, gtovec, xcoords, ycoords, zcoords, 0, 0, npoints);

        CDenseMatrix mat_chi(naos, npoints);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            std::memcpy(mat_chi.row(nu), gaos.data(nu), npoints * sizeof(double));
        }

        // generate density grid

        auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

        auto xcfuntype = fvxc.getFunctionalType();

        auto dengrid = _generateDensityGrid(npoints, mat_chi, densityMatrix, xcfuntype);

        // compute exchange-correlation functional derivative

        CXCGradientGrid vxcgrid(npoints, dengrid.getDensityGridType(), xcfuntype);

        fvxc.compute(vxcgrid, dengrid);

        // compute partial contribution to Vxc matrix

        CDenseMatrix partial_mat_Vxc;

        if (xcfuntype == xcfun::lda)
        {
            partial_mat_Vxc = _integratePartialVxcFockForLDA(npoints, xcoords, ycoords, zcoords, weights, mat_chi, vxcgrid);
        }
        else if (xcfuntype == xcfun::gga)
        {
            //partial_mat_Vxc = _integratePartialVxcFockForGGA(npoints, xcoords, ycoords, zcoords, weights, mat_chi, vxcgrid);
        }

        for (int32_t i = 0; i < partial_mat_Vxc.getNumberOfElements(); i++)
        {
            mat_Vxc.getMatrix(0)[i] += partial_mat_Vxc.values()[i];
        }

        // compute partial contribution to XC energy

        auto rhoa = dengrid.alphaDensity(0);

        auto efunc = vxcgrid.xcFunctionalValues();

        for (int32_t i = 0; i < npoints; i++)
        {
            nele += weights[i] * rhoa[i];

            xcene += weights[i] * efunc[i];
        }
    }

    // destroy GTOs container

    delete gtovec;

    mat_Vxc.setNumberOfElectrons(nele);

    mat_Vxc.setExchangeCorrelationEnergy(xcene);

    return mat_Vxc;
}

CDensityGrid
CXCNewIntegrator::_generateDensityGrid(const int32_t           npoints,
                                       const CDenseMatrix&     gtoValuesOnGridPoints,
                                       const CAODensityMatrix& densityMatrix,
                                       const xcfun             xcFunType) const
{
    CDensityGrid dengrid(npoints, densityMatrix.getNumberOfDensityMatrices(), xcFunType, dengrid::ab);

    dengrid.zero();

    auto rhoa = dengrid.alphaDensity(0);

    auto rhob = dengrid.betaDensity(0);

    // eq.(26), JCTC 2021, 17, 1512-1521

    const CDenseMatrix& mat_chi = gtoValuesOnGridPoints;

    auto mat_F = denblas::multAB(densityMatrix.getReferenceToDensity(0), mat_chi);

    // eq.(27), JCTC 2021, 17, 1512-1521

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

    return dengrid;
}

CDenseMatrix
CXCNewIntegrator::_integratePartialVxcFockForLDA(const int32_t          npoints,
                                                 const double*          xcoords,
                                                 const double*          ycoords,
                                                 const double*          zcoords,
                                                 const double*          weights,
                                                 const CDenseMatrix&    gtoValuesOnGridPoints,
                                                 const CXCGradientGrid& xcGradientGrid) const
{
    // GTO values on grid points

    const CDenseMatrix& mat_chi = gtoValuesOnGridPoints;

    // exchange-correlation functional derivative

    auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);

    // eq.(30), JCTC 2021, 17, 1512-1521

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

    // eq.(31), JCTC 2021, 17, 1512-1521

    return denblas::multABt(mat_chi, mat_G);
}
