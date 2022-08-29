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

#include "MolecularGrid.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

#include <mpi.h>

#include "DenseLinearAlgebra.hpp"
#include "FunctionalParser.hpp"
#include "GtoContainer.hpp"
#include "GtoFunc.hpp"
#include "XCVarsType.hpp"

CMolecularGrid::CMolecularGrid()

    : _isDistributed(false)
{
}

CMolecularGrid::CMolecularGrid(const CMemBlock2D<double>& gridPoints)

    : _isDistributed(false)

    , _gridPoints(gridPoints)
{
}

CMolecularGrid::CMolecularGrid(const CMolecularGrid& source)

    : _isDistributed(source._isDistributed)

    , _gridPoints(source._gridPoints)
{
}

CMolecularGrid::CMolecularGrid(CMolecularGrid&& source) noexcept

    : _isDistributed(std::move(source._isDistributed))

    , _gridPoints(std::move(source._gridPoints))
{
}

CMolecularGrid::~CMolecularGrid()
{
}

CMolecularGrid&
CMolecularGrid::operator=(const CMolecularGrid& source)
{
    if (this == &source) return *this;

    _isDistributed = source._isDistributed;

    _gridPoints = source._gridPoints;

    return *this;
}

CMolecularGrid&
CMolecularGrid::operator=(CMolecularGrid&& source) noexcept
{
    if (this == &source) return *this;

    _isDistributed = std::move(source._isDistributed);

    _gridPoints = std::move(source._gridPoints);

    return *this;
}

bool
CMolecularGrid::operator==(const CMolecularGrid& other) const
{
    if (_isDistributed != other._isDistributed) return false;

    if (_gridPoints != other._gridPoints) return false;

    return true;
}

bool
CMolecularGrid::operator!=(const CMolecularGrid& other) const
{
    return !(*this == other);
}

void
CMolecularGrid::slice(const int32_t nGridPoints)
{
    if (nGridPoints < getNumberOfGridPoints())
    {
        _gridPoints = _gridPoints.slice(0, nGridPoints); 
    }
}

CMemBlock2D<double>
CMolecularGrid::getGridPoints() const
{
    return _gridPoints;
}

int32_t
CMolecularGrid::getNumberOfGridPoints() const
{
    return _gridPoints.size(0);
}

const double*
CMolecularGrid::getCoordinatesX() const
{
    return _gridPoints.data(0);
}

double*
CMolecularGrid::getCoordinatesX()
{
    return _gridPoints.data(0);
}

const double*
CMolecularGrid::getCoordinatesY() const
{
    return _gridPoints.data(1);
}

double*
CMolecularGrid::getCoordinatesY()
{
    return _gridPoints.data(1);
}

const double*
CMolecularGrid::getCoordinatesZ() const
{
    return _gridPoints.data(2);
}

double*
CMolecularGrid::getCoordinatesZ()
{
    return _gridPoints.data(2);
}

const double*
CMolecularGrid::getWeights() const
{
    return _gridPoints.data(3);
}

double*
CMolecularGrid::getWeights()
{
    return _gridPoints.data(3);
}

void
CMolecularGrid::distribute(int32_t rank, int32_t nodes, MPI_Comm comm)
{
    if (!_isDistributed)
    {
        _isDistributed = true;

        _gridPoints.scatter(rank, nodes, comm);
    }
}

void
CMolecularGrid::broadcast(int32_t rank, MPI_Comm comm)
{
    _gridPoints.broadcast(rank, comm);
}

void
CMolecularGrid::read_blocked_grid(const std::string& fileName)
{
    const int32_t mpoints = 20000000;
    
    CMemBlock2D<double> rgrid(mpoints, 4);
    
    auto f = std::fopen(fileName.c_str(), "rb");
    
    if (f != nullptr)
    {
        int32_t npnt = 0;
        
        int32_t cpnt = 0;
        
        while(true)
        {
            // read number of grid points in block
            
            auto nblock = std::fread(&npnt, sizeof(int32_t), 1u, f);
            
            if (nblock < 1u) break;
            
            // read shell block data
            
            int32_t nshells = 0;
            
            nblock = std::fread(&nshells, sizeof(int32_t), 1u, f);
            
            if (nblock < 1u) break;
            
            for(int32_t i = 0; i < 2 * nshells; i++)
            {
                int32_t bf = 0;
                
                nblock = std::fread(&bf, sizeof(int32_t), 1u, f);
                
                if (nblock < 1u) break;
            }
            
            // read coordinates and weights

            auto rx = rgrid.data(0, cpnt);
            
            auto ry = rgrid.data(1, cpnt);
            
            auto rz = rgrid.data(2, cpnt);
            
            auto rw = rgrid.data(3, cpnt);
            
            nblock = std::fread(rx, sizeof(double), npnt, f);
            
            if (nblock < 1u) break;
            
            nblock = std::fread(ry, sizeof(double), npnt, f);
            
            if (nblock < 1u) break;
            
            nblock = std::fread(rz, sizeof(double), npnt, f);
            
            if (nblock < 1u) break;
            
            nblock = std::fread(rw, sizeof(double), npnt, f);
            
            if (nblock < 1u) break;
            
            cpnt += npnt;
        }
        
        _isDistributed = false;
        
        _gridPoints = rgrid.slice(0, cpnt);
    }
    
    std::fclose(f);
}

void
CMolecularGrid::read_raw_grid(const std::string& fileName)
{
    // allocate raw grid
    
    const int32_t mpoints = 20000000;
    
    CMemBlock2D<double> rgrid(mpoints, 4);
    
    // set up pointers to grid data
    
    auto rx = rgrid.data(0);
    
    auto ry = rgrid.data(1);
    
    auto rz = rgrid.data(2);
    
    auto rw = rgrid.data(3);
    
    // read grid data from text file
    
    std::ifstream fst(fileName.c_str());
    
    int32_t npoints = 0;
    
    while(true)
    {
        std::string str;
        
        std::getline(fst, str);

        if (fst.eof()) break;
        
        std::istringstream iss(str);
        
        int32_t idx = 0;
        
        iss >> idx >> rx[npoints] >> ry[npoints] >> rz[npoints] >> rw[npoints];
        
        npoints++;
    }
    
    _isDistributed = false;
    
    _gridPoints = rgrid.slice(0, npoints);
}

void
CMolecularGrid::write_raw_grid(const std::string fileName) const
{
    auto f = std::fopen(fileName.c_str(), "w");
    
    if (f != nullptr)
    {
        auto rx = getCoordinatesX();
        
        auto ry = getCoordinatesY();
        
        auto rz = getCoordinatesZ();
        
        auto rw = getWeights();
        
        for (int32_t i = 0; i < getNumberOfGridPoints(); i++)
        {
            std::fprintf(f, "%8i %20.16lf %20.16lf %20.16lf %20.16lf\n", i, rx[i], ry[i], rz[i], rw[i]);
        }
    }
    
    std::fclose(f);
}

std::array<double, 6>
CMolecularGrid::getSpatialExtent() const
{
    // initialize min, max coordinates
    
    double xmin = 0.0, ymin = 0.0, zmin = 0.0;
    
    double xmax = 0.0, ymax = 0.0, zmax = 0.0;
    
    auto ngpoints = getNumberOfGridPoints();
    
    if (ngpoints > 0)
    {
        // set up pointers to grid data
        
        auto rx = getCoordinatesX();
        
        auto ry = getCoordinatesY();
        
        auto rz = getCoordinatesZ();
        
        // update min, max coordinates
        
        xmin = rx[0]; xmax = rx[0];
        
        ymin = ry[0]; ymax = ry[0];
        
        zmin = rz[0]; zmax = rz[0];
        
        for (int32_t i = 0; i < ngpoints; i++)
        {
            if (xmin > rx[i]) xmin = rx[i];
            
            if (ymin > ry[i]) ymin = ry[i];
            
            if (zmin > rz[i]) zmin = rz[i];
            
            if (xmax < rx[i]) xmax = rx[i];
            
            if (ymax < ry[i]) ymax = ry[i];
            
            if (zmax < rz[i]) zmax = rz[i];
        }
    }
    
    std::array<double, 6> rext = {xmin, ymin, zmin, xmax, ymax, zmax};
    
    return rext;
}

CDenseMatrix
CMolecularGrid::testPartition(const CMolecule&        molecule,
                              const CMolecularBasis&  basis,
                              const CAODensityMatrix& density,
                              const std::string&      xcFuncLabel,
                              CDensityGridDriver      dgdrv) const
{
    // determine the boundary of grid points

    auto spatial_extent = getSpatialExtent();

    // xmin, ymin, zmin

    spatial_extent[0] -= 0.0001;
    spatial_extent[1] -= 0.0001;
    spatial_extent[2] -= 0.0001;

    // xmax, ymax, zmax

    spatial_extent[3] += 0.0001;
    spatial_extent[4] += 0.0001;
    spatial_extent[5] += 0.0001;

    // grid point info

    auto totalgridpoints = getNumberOfGridPoints();
    auto allgridpoints = getGridPoints();

    // start with one box

    std::vector<std::array<double ,6>> boxes;
    std::vector<CMemBlock2D<double>> points_in_boxes;

    boxes.push_back(spatial_extent);
    points_in_boxes.push_back(allgridpoints);

    // divide the boxes

    const int32_t count_thresh = 800;

    for(int32_t level = 0; level < 100; level++)
    {
        auto numboxes = static_cast<int32_t>(boxes.size());

        // check number of grid points in the boxes

        bool done = true;

        for (int32_t k = 0; k < numboxes; k++)
        {
            if (points_in_boxes[k].size(0) > count_thresh)
            {
                done = false;
                break;
            }
        }

        if (done) break;

        // start dividing the boxes

        std::cout << "Level: " << level << std::endl;

        std::vector<std::array<double ,6>> new_boxes;
        std::vector<CMemBlock2D<double>> new_points_in_boxes;

        for (int32_t k = 0; k < numboxes; k++)
        {
            if (points_in_boxes[k].size(0) <= count_thresh)
            {
                new_boxes.push_back(boxes[k]);
                new_points_in_boxes.push_back(points_in_boxes[k]);
            }
            else
            {
                auto xmin = boxes[k][0];
                auto ymin = boxes[k][1];
                auto zmin = boxes[k][2];

                auto xmax = boxes[k][3];
                auto ymax = boxes[k][4];
                auto zmax = boxes[k][5];

                double xhalf = 0.5 * (xmax + xmin);
                double yhalf = 0.5 * (ymax + ymin);
                double zhalf = 0.5 * (zmax + zmin);

                std::vector<std::array<double, 6>> newboxdims = {
                    std::array<double, 6>({xmin, ymin, zmin, xhalf, yhalf, zhalf}),
                    std::array<double, 6>({xhalf, ymin, zmin, xmax, yhalf, zhalf}),
                    std::array<double, 6>({xmin, yhalf, zmin, xhalf, ymax, zhalf}),
                    std::array<double, 6>({xmin, ymin, zhalf, xhalf, yhalf, zmax}),
                    std::array<double, 6>({xhalf, yhalf, zmin, xmax, ymax, zhalf}),
                    std::array<double, 6>({xhalf, ymin, zhalf, xmax, yhalf, zmax}),
                    std::array<double, 6>({xmin, yhalf, zhalf, xhalf, ymax, zmax}),
                    std::array<double, 6>({xhalf, yhalf, zhalf, xmax, ymax, zmax}),
                };

                auto points = points_in_boxes[k];

                auto numdims = static_cast<int32_t>(newboxdims.size());

                for (int32_t d = 0; d < numdims; d++)
                {
                    auto newpoints = getGridPointsInBox(newboxdims[d], points);
                    new_boxes.push_back(newboxdims[d]);
                    new_points_in_boxes.push_back(newpoints);
                }
            }
        }

        boxes = new_boxes;
        points_in_boxes = new_points_in_boxes;
    }

    auto numboxes = static_cast<int32_t>(boxes.size());
    int32_t count_sum = 0;

    for (int32_t k = 0; k < numboxes; k++)
    {
        std::cout << "box " << k << ": ";
        std::cout << "(" << boxes[k][0] << ", " << boxes[k][1] << ", " << boxes[k][2] << "), ";
        std::cout << "(" << boxes[k][3] << ", " << boxes[k][4] << ", " << boxes[k][5] << ")  ";
        std::cout << "n=" << points_in_boxes[k].size(0) << std::endl;
        count_sum += points_in_boxes[k].size(0);
    }

    std::cout << "----" << std::endl;
    std::cout << "Sum of grid points in boxes: " << count_sum << std::endl;
    std::cout << "Total number of grid points: " << totalgridpoints << std::endl;

    // create GTOs container

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    // set up number of AOs

    auto naos = gtovec->getNumberOfAtomicOrbitals();

    CDenseMatrix mat_Vxc(naos, naos);

    mat_Vxc.zero();

    for (int32_t k = 0; k < numboxes; k++)
    {
        auto points = points_in_boxes[k];

        auto npoints = points.size(0);
        auto xcoords = points.data(0);
        auto ycoords = points.data(1);
        auto zcoords = points.data(2);
        auto weights = points.data(3);

        std::cout << "Computing GTOs on " << npoints << " grid points..." << std::endl;

        CMemBlock2D<double> gaos(npoints, naos);
        CMemBlock2D<double> gaox(npoints, naos);
        CMemBlock2D<double> gaoy(npoints, naos);
        CMemBlock2D<double> gaoz(npoints, naos);

        gaos.zero();
        gaox.zero();
        gaoy.zero();
        gaoz.zero();

        gtorec::computeGtosValuesForGGA(gaos, gaox, gaoy, gaoz, gtovec, xcoords, ycoords, zcoords, 0, 0, npoints);

        // chi: GTO at grid points

        CDenseMatrix mat_chi(naos, npoints);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            for (int32_t g = 0; g < npoints; g++)
            {
                mat_chi.row(nu)[g] = gaos.data(nu)[g];
            }
        }

        // eq.(26), JCTC 2021, 17, 1512-1521

        auto mat_F = denblas::multAB(density.getReferenceToDensity(0), mat_chi);

        // eq.(27), JCTC 2021, 17, 1512-1521

        CMemBlock<double> rho(npoints);

        rho.zero();

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto F_nu = mat_F.row(nu);
            auto chi_nu = mat_chi.row(nu);

            for (int32_t g = 0; g < npoints; g++)
            {
                rho.at(g) += F_nu[g] * chi_nu[g];
            }
        }

        // generate reference density grid and compute exchange-correlation functional derivative

        CMolecularGrid mgrid(points);

        auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

        auto dengrid = dgdrv.generate(density, molecule, basis, mgrid, fvxc.getFunctionalType());

        CXCGradientGrid vxcgrid(mgrid.getNumberOfGridPoints(), dengrid.getDensityGridType(), fvxc.getFunctionalType());

        fvxc.compute(vxcgrid, dengrid);

        auto grhoa = vxcgrid.xcGradientValues(xcvars::rhoa);

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

        auto mat_V = denblas::multABt(mat_chi, mat_G);
        mat_Vxc = denblas::addAB(mat_Vxc, mat_V, 1.0);
    }

    // destroy GTOs container

    delete gtovec;

    return mat_Vxc;
}

CMemBlock2D<double>
CMolecularGrid::getGridPointsInBox(const std::array<double, 6>& boxdim,
                                   const CMemBlock2D<double>&   points) const
{
    auto npoints = points.size(0);

    auto xcoords = points.data(0);
    auto ycoords = points.data(1);
    auto zcoords = points.data(2);
    auto weights = points.data(3);

    CMemBlock2D<double> new_points (npoints, 4);

    auto xmin = boxdim[0];
    auto ymin = boxdim[1];
    auto zmin = boxdim[2];

    auto xmax = boxdim[3];
    auto ymax = boxdim[4];
    auto zmax = boxdim[5];

    int32_t count = 0;

    for(int32_t i = 0; i < npoints; i++)
    {
        if (xmin <= xcoords[i] && xcoords[i] < xmax &&
            ymin <= ycoords[i] && ycoords[i] < ymax &&
            zmin <= zcoords[i] && zcoords[i] < zmax)
        {
            new_points.data(0)[count] = xcoords[i];
            new_points.data(1)[count] = ycoords[i];
            new_points.data(2)[count] = zcoords[i];
            new_points.data(3)[count] = weights[i];

            count++;
        }
    }

    return new_points.slice(0, count);
}

std::ostream&
operator<<(std::ostream& output, const CMolecularGrid& source)
{
    output << std::endl;
    
    output << "[CMolecularGrid (Object):" << &source << "]" << std::endl;

    output << "_isDistributed: " << source._isDistributed << std::endl;

    output << "_gridPoints: " << std::endl;

    output << source._gridPoints << std::endl;

    return output;
}
