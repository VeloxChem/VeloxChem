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

#include "DensityGrid.hpp"

#include <cmath>
#include <iostream>

CDensityGrid::CDensityGrid()

    : _gridType(dengrid::undefined)

    , _nDensityMatrices(0)

    , _densityValues(CMemBlock2D<double>())
{
}

CDensityGrid::CDensityGrid(const CMemBlock2D<double>& densityValues, const dengrid gridType)

    : _gridType(gridType)

    , _nDensityMatrices(1)

    , _densityValues(densityValues)
{
}

CDensityGrid::CDensityGrid(const int32_t nGridPoints, const int32_t nDensityMatrices, const xcfun xcFuncType, const dengrid gridType)
{
    _gridType = gridType;

    _nDensityMatrices = nDensityMatrices;

    int32_t ncomp = 0;

    if (xcFuncType == xcfun::lda) ncomp = (_gridType == dengrid::ab) ? 2 : 1;

    if (xcFuncType == xcfun::gga) ncomp = (_gridType == dengrid::ab) ? 11 : 5;

    // NOTE: this needs to be checked with mgga functionals implementation

    if (xcFuncType == xcfun::mgga) ncomp = (_gridType == dengrid::ab) ? 13 : 6;

    _densityValues = CMemBlock2D<double>(nGridPoints, _nDensityMatrices * ncomp);
}

CDensityGrid::CDensityGrid(const int32_t nGridPoints,
                           const int32_t nDensityMatrices,
                           const int32_t nComponents,
                           const dengrid gridType)
{
    _gridType = gridType;
    
    _nDensityMatrices = nDensityMatrices;
    
    _densityValues = CMemBlock2D<double>(nGridPoints, _nDensityMatrices * nComponents);
}

CDensityGrid::CDensityGrid(const CDensityGrid& source)

    : _gridType(source._gridType)

    , _nDensityMatrices(source._nDensityMatrices)

    , _densityValues(source._densityValues)
{
}

CDensityGrid::CDensityGrid(CDensityGrid&& source) noexcept

    : _gridType(std::move(source._gridType))

    , _nDensityMatrices(std::move(source._nDensityMatrices))

    , _densityValues(std::move(source._densityValues))
{
}

CDensityGrid::~CDensityGrid()
{
}

CDensityGrid&
CDensityGrid::operator=(const CDensityGrid& source)
{
    if (this == &source) return *this;

    _gridType = source._gridType;

    _nDensityMatrices = source._nDensityMatrices;

    _densityValues = source._densityValues;

    return *this;
}

CDensityGrid&
CDensityGrid::operator=(CDensityGrid&& source) noexcept
{
    if (this == &source) return *this;

    _gridType = std::move(source._gridType);

    _nDensityMatrices = std::move(source._nDensityMatrices);

    _densityValues = std::move(source._densityValues);

    return *this;
}

bool
CDensityGrid::operator==(const CDensityGrid& other) const
{
    if (_gridType != other._gridType) return false;

    if (_nDensityMatrices != other._nDensityMatrices) return false;

    if (_densityValues != other._densityValues) return false;

    return true;
}

bool
CDensityGrid::operator!=(const CDensityGrid& other) const
{
    return !(*this == other);
}

void
CDensityGrid::zero()
{
    _densityValues.zero();
}

void
CDensityGrid::slice(const int32_t nGridPoints)
{
    if (nGridPoints < getNumberOfGridPoints())
    {
        _densityValues = _densityValues.slice(0, nGridPoints);
    }
}

void
CDensityGrid::updateBetaDensities()
{
    if (_gridType != dengrid::ab) return;

    auto ngpoints = getNumberOfGridPoints();

    if ((2 * _nDensityMatrices) == _densityValues.blocks())
    {
        for (int32_t i = 0; i < _nDensityMatrices; i++)
        {
            auto rhoa = alphaDensity(i);

            auto rhob = betaDensity(i);

            #pragma omp simd aligned(rhoa, rhob : VLX_ALIGN)
            for (int32_t j = 0; j < ngpoints; j++)
            {
                rhob[j] = rhoa[j];
            }
        }
    }

    if ((11 * _nDensityMatrices) == _densityValues.blocks())
    {
        for (int32_t i = 0; i < _nDensityMatrices; i++)
        {
            auto rhoa = alphaDensity(i);

            auto rhob = betaDensity(i);

            auto grada_x = alphaDensityGradientX(i);

            auto gradb_x = betaDensityGradientX(i);

            auto grada_y = alphaDensityGradientY(i);

            auto gradb_y = betaDensityGradientY(i);

            auto grada_z = alphaDensityGradientZ(i);

            auto gradb_z = betaDensityGradientZ(i);

            #pragma omp simd aligned(rhoa, rhob, grada_x, gradb_x, grada_y, gradb_y, grada_z, gradb_z : VLX_ALIGN)
            for (int32_t j = 0; j < ngpoints; j++)
            {
                rhob[j] = rhoa[j];

                gradb_x[j] = grada_x[j];

                gradb_y[j] = grada_y[j];

                gradb_z[j] = grada_z[j];
            }
        }
    }

    if ((13 * _nDensityMatrices) == _densityValues.blocks())
    {
        for (int32_t i = 0; i < _nDensityMatrices; i++)
        {
            auto rhoa = alphaDensity(i);

            auto rhob = betaDensity(i);

            auto grada_x = alphaDensityGradientX(i);

            auto gradb_x = betaDensityGradientX(i);

            auto grada_y = alphaDensityGradientY(i);

            auto gradb_y = betaDensityGradientY(i);

            auto grada_z = alphaDensityGradientZ(i);

            auto gradb_z = betaDensityGradientZ(i);

            auto taua = alphaDensitytau(i);

            auto taub = betaDensitytau(i);

            #pragma omp simd aligned(rhoa, rhob, grada_x, gradb_x, grada_y, gradb_y, grada_z, gradb_z, taua, taub : VLX_ALIGN)
            for (int32_t j = 0; j < ngpoints; j++)
            {
                rhob[j] = rhoa[j];

                gradb_x[j] = grada_x[j];

                gradb_y[j] = grada_y[j];

                gradb_z[j] = grada_z[j];

                taub[j] = taua[j];
            }
        }
    }
}

void
CDensityGrid::computeDensityNorms()
{
    auto ngpoints = getNumberOfGridPoints();

    if (_gridType == dengrid::ab)
    {
        if ((11 * _nDensityMatrices) == _densityValues.blocks())
        {
            for (int32_t i = 0; i < _nDensityMatrices; i++)
            {
                auto grada_x = alphaDensityGradientX(i);

                auto gradb_x = betaDensityGradientX(i);

                auto grada_y = alphaDensityGradientY(i);

                auto gradb_y = betaDensityGradientY(i);

                auto grada_z = alphaDensityGradientZ(i);

                auto gradb_z = betaDensityGradientZ(i);

                auto grada = alphaDensityGradient(i);

                auto gradb = betaDensityGradient(i);

                auto gradab = mixedDensityGradient(i);

                #pragma omp simd aligned(grada_x, gradb_x, grada_y, gradb_y, grada_z, gradb_z, grada, gradb, gradab : VLX_ALIGN)
                for (int32_t j = 0; j < ngpoints; j++)
                {
                    grada[j] = std::sqrt(grada_x[j] * grada_x[j] + grada_y[j] * grada_y[j] + grada_z[j] * grada_z[j]);

                    gradb[j] = std::sqrt(gradb_x[j] * gradb_x[j] + gradb_y[j] * gradb_y[j] + gradb_z[j] * gradb_z[j]);

                    gradab[j] = grada_x[j] * gradb_x[j] + grada_y[j] * gradb_y[j] + grada_z[j] * gradb_z[j];
                }
            }
        }

        if ((13 * _nDensityMatrices) == _densityValues.blocks())
        {
            for (int32_t i = 0; i < _nDensityMatrices; i++)
            {
                auto grada_x = alphaDensityGradientX(i);

                auto gradb_x = betaDensityGradientX(i);

                auto grada_y = alphaDensityGradientY(i);

                auto gradb_y = betaDensityGradientY(i);

                auto grada_z = alphaDensityGradientZ(i);

                auto gradb_z = betaDensityGradientZ(i);

                auto grada = alphaDensityGradient(i);

                auto gradb = betaDensityGradient(i);

                auto gradab = mixedDensityGradient(i);

                #pragma omp simd aligned(grada_x, gradb_x, grada_y, gradb_y, grada_z, gradb_z, grada, gradb, gradab : VLX_ALIGN)
                for (int32_t j = 0; j < ngpoints; j++)
                {
                    grada[j] = std::sqrt(grada_x[j] * grada_x[j] + grada_y[j] * grada_y[j] + grada_z[j] * grada_z[j]);

                    gradb[j] = std::sqrt(gradb_x[j] * gradb_x[j] + gradb_y[j] * gradb_y[j] + gradb_z[j] * gradb_z[j]);

                    gradab[j] = grada_x[j] * gradb_x[j] + grada_y[j] * gradb_y[j] + grada_z[j] * gradb_z[j];
                }
            }
        }
    }

    if (_gridType == dengrid::lima)
    {
        if ((5 * _nDensityMatrices) == _densityValues.blocks())
        {
            for (int32_t i = 0; i < _nDensityMatrices; i++)
            {
                auto gradb_x = betaDensityGradientX(i);

                auto gradb_y = betaDensityGradientY(i);

                auto gradb_z = betaDensityGradientZ(i);

                auto gradb = betaDensityGradient(i);

                #pragma omp simd aligned(gradb_x, gradb_y, gradb_z, gradb : VLX_ALIGN)
                for (int32_t j = 0; j < ngpoints; j++)
                {
                    gradb[j] = std::sqrt(gradb_x[j] * gradb_x[j] + gradb_y[j] * gradb_y[j] + gradb_z[j] * gradb_z[j]);
                }
            }
        }
    }

    if (_gridType == dengrid::limb)
    {
        if ((5 * _nDensityMatrices) == _densityValues.blocks())
        {
            for (int32_t i = 0; i < _nDensityMatrices; i++)
            {
                auto grada_x = alphaDensityGradientX(i);

                auto grada_y = alphaDensityGradientY(i);

                auto grada_z = alphaDensityGradientZ(i);

                auto grada = alphaDensityGradient(i);

                #pragma omp simd aligned(grada_x, grada_y, grada_z, grada : VLX_ALIGN)
                for (int32_t j = 0; j < ngpoints; j++)
                {
                    grada[j] = std::sqrt(grada_x[j] * grada_x[j] + grada_y[j] * grada_y[j] + grada_z[j] * grada_z[j]);
                }
            }
        }
    }
}

int32_t
CDensityGrid::getNumberOfGridPoints() const
{
    return _densityValues.size(0);
}

int32_t
CDensityGrid::getNumberOfDensityMatrices() const
{
    return _nDensityMatrices;
}

dengrid
CDensityGrid::getDensityGridType() const
{
    return _gridType;
}

const double*
CDensityGrid::alphaDensity(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::lima) return nullptr;

    return _densityValues.data(iDensityMatrix);
}

double*
CDensityGrid::alphaDensity(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::lima) return nullptr;

    return _densityValues.data(iDensityMatrix);
}

const double*
CDensityGrid::betaDensity(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(_nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::lima) return _densityValues.data(iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::betaDensity(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(_nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::lima) return _densityValues.data(iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::alphaDensityGradient(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::limb) return _densityValues.data(_nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::alphaDensityGradient(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::limb) return _densityValues.data(_nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::betaDensityGradient(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::lima) return _densityValues.data(_nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::betaDensityGradient(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::lima) return _densityValues.data(_nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::mixedDensityGradient(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::mixedDensityGradient(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::alphaDensityGradientX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::limb) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::alphaDensityGradientX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::limb) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::alphaDensityGradientY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::limb) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::alphaDensityGradientY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::limb) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::alphaDensityGradientZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::limb) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::alphaDensityGradientZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::limb) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::betaDensityGradientX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::lima) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::betaDensityGradientX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::lima) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::betaDensityGradientY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::lima) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::betaDensityGradientY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::lima) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::betaDensityGradientZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::lima) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::betaDensityGradientZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);

    if (_gridType == dengrid::lima) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::alphaDensitytau(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(11 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::alphaDensitytau(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(11 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::betaDensitytau(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(12 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGrid::betaDensitytau(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(12 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGrid::getComponent(const int32_t iComponent) const
{
    return _densityValues.data(iComponent);
}

double*
CDensityGrid::getComponent(const int32_t iComponent)
{
    return _densityValues.data(iComponent);
}

void
CDensityGrid::getScreenedGridsPair(CDensityGrid&   densityGridAB,
                                   CMolecularGrid& molecularGridab,
                                   const int32_t   iDensityMatrix,
                                   const double    densityThreshold,
                                   const xcfun     xcFuncType) const
{
    if (_gridType != dengrid::ab) return;

    // create density grid

    densityGridAB = CDensityGrid(getNumberOfGridPoints(), 1, xcFuncType, _gridType);

    // generate screened molecular grid

    molecularGridab = getScreenedGrid(molecularGridab, iDensityMatrix, densityThreshold, xcFuncType);

    // set grid points data

    auto npoints = getNumberOfGridPoints();

    int32_t ipoints = 0;

    // set up pointers to source density

    auto srhoa = alphaDensity(iDensityMatrix);

    auto srhob = betaDensity(iDensityMatrix);

    // set up pointers to destination density

    auto drhoa = densityGridAB.alphaDensity(0);

    auto drhob = densityGridAB.betaDensity(0);

    // density screening for LDA

    if (xcFuncType == xcfun::lda)
    {
        for (int32_t i = 0; i < npoints; i++)
        {
            if (_isValidGridPointForLda(srhoa[i], srhob[i], densityThreshold))
            {
                drhoa[ipoints] = srhoa[i];
                drhob[ipoints] = srhob[i];

                ipoints++;
            }
        }
    }

    // set up pointers to source density gradient

    auto sgrada = alphaDensityGradient(iDensityMatrix);

    auto sgradb = betaDensityGradient(iDensityMatrix);

    auto sgradab = mixedDensityGradient(iDensityMatrix);

    auto sgrada_x = alphaDensityGradientX(iDensityMatrix);

    auto sgrada_y = alphaDensityGradientY(iDensityMatrix);

    auto sgrada_z = alphaDensityGradientZ(iDensityMatrix);

    auto sgradb_x = betaDensityGradientX(iDensityMatrix);

    auto sgradb_y = betaDensityGradientY(iDensityMatrix);

    auto sgradb_z = betaDensityGradientZ(iDensityMatrix);

    // set up pointers to destination density gradient

    auto dgrada = densityGridAB.alphaDensityGradient(0);

    auto dgradb = densityGridAB.betaDensityGradient(0);

    auto dgradab = densityGridAB.mixedDensityGradient(0);

    auto dgrada_x = densityGridAB.alphaDensityGradientX(0);

    auto dgrada_y = densityGridAB.alphaDensityGradientY(0);

    auto dgrada_z = densityGridAB.alphaDensityGradientZ(0);

    auto dgradb_x = densityGridAB.betaDensityGradientX(0);

    auto dgradb_y = densityGridAB.betaDensityGradientY(0);

    auto dgradb_z = densityGridAB.betaDensityGradientZ(0);

    // density screening for GGA

    if (xcFuncType == xcfun::gga)
    {
        for (int32_t i = 0; i < npoints; i++)
        {
            if (_isValidGridPointForGga(srhoa[i], srhob[i], sgrada[i], sgradb[i], densityThreshold))
            {
                drhoa[ipoints] = srhoa[i];
                drhob[ipoints] = srhob[i];

                dgrada[ipoints]  = sgrada[i];
                dgradb[ipoints]  = sgradb[i];
                dgradab[ipoints] = sgradab[i];

                dgrada_x[ipoints] = sgrada_x[i];
                dgrada_y[ipoints] = sgrada_y[i];
                dgrada_z[ipoints] = sgrada_z[i];

                dgradb_x[ipoints] = sgradb_x[i];
                dgradb_y[ipoints] = sgradb_y[i];
                dgradb_z[ipoints] = sgradb_z[i];

                ipoints++;
            }
        }
    }

    // density screening for m-GGA

    // set up pointers to source  tau

    auto staua = alphaDensitytau(iDensityMatrix);

    auto staub = betaDensitytau(iDensityMatrix);

    // set up pointers to destination  tau

    auto dtaua = densityGridAB.alphaDensitytau(0);

    auto dtaub = densityGridAB.betaDensitytau(0);

    if (xcFuncType == xcfun::mgga)
    {
        for (int32_t i = 0; i < npoints; i++)
        {
            if (_isValidGridPointForMgga(srhoa[i], srhob[i], sgrada[i], sgradb[i], staua[i], staub[i], densityThreshold))
            {
                drhoa[ipoints] = srhoa[i];
                drhob[ipoints] = srhob[i];

                dgrada[ipoints]  = sgrada[i];
                dgradb[ipoints]  = sgradb[i];
                dgradab[ipoints] = sgradab[i];

                dgrada_x[ipoints] = sgrada_x[i];
                dgrada_y[ipoints] = sgrada_y[i];
                dgrada_z[ipoints] = sgrada_z[i];

                dgradb_x[ipoints] = sgradb_x[i];
                dgradb_y[ipoints] = sgradb_y[i];
                dgradb_z[ipoints] = sgradb_z[i];

                dtaua[ipoints] = staua[i];
                dtaub[ipoints] = staub[i];

                ipoints++;
            }
        }
    }

    // compress screened density grid size

    densityGridAB.slice(ipoints);
}

void
CDensityGrid::getScreenedGridPairUnrestricted(CDensityGrid&   densityGridAB,
                                              CDensityGrid&   densityGridA,
                                              CDensityGrid&   densityGridB,
                                              CMolecularGrid& molecularGridab,
                                              CMolecularGrid& molecularGrida,
                                              CMolecularGrid& molecularGridb,
                                              const int32_t   iDensityMatrix,
                                              const double    densityThreshold,
                                              const xcfun     xcFuncType) const
{
    // create density grid

    densityGridAB = CDensityGrid(getNumberOfGridPoints(), 1, xcFuncType, dengrid::ab);

    densityGridA = CDensityGrid(getNumberOfGridPoints(), 1, xcFuncType, dengrid::lima);

    densityGridB = CDensityGrid(getNumberOfGridPoints(), 1, xcFuncType, dengrid::limb);

    // generate screened molecular grid

    getScreenedGridUnrestricted(molecularGridab, molecularGrida, molecularGridb, iDensityMatrix, densityThreshold, xcFuncType);

    // set grid points data

    auto npoints = getNumberOfGridPoints();

    int32_t ipoints = 0;

    int32_t ipointsa = 0;

    int32_t ipointsb = 0;

    // set up pointers to source density

    auto srhoa = alphaDensity(iDensityMatrix);

    auto srhob = betaDensity(iDensityMatrix);

    // set up pointers to destination density

    auto drhoa = densityGridAB.alphaDensity(0);

    auto drhob = densityGridAB.betaDensity(0);

    auto drhobA = densityGridA.betaDensity(0);

    auto drhoaB = densityGridB.alphaDensity(0);

    // density screening for LDA

    if (xcFuncType == xcfun::lda)
    {
        for (int32_t i = 0; i < npoints; i++)
        {
            double srhoA = srhoa[i];

            double srhoB = srhob[i];

            if (_isValidGridPointForLdaUnrestricted(srhoa[i], srhob[i], densityThreshold))
            {
                drhoa[ipoints] = srhoA;

                drhob[ipoints] = srhoB;

                ipoints++;

                continue;
            }

            if (_isValidGridPointForLdaA(srhoa[i], srhob[i], densityThreshold))
            {
                drhobA[ipointsa] = srhoB;

                ipointsa++;

                continue;
            }

            if (_isValidGridPointForLdaB(srhoa[i], srhob[i], densityThreshold))
            {
                drhoaB[ipointsb] = srhoA;

                ipointsb++;
            }
        }
    }

    // set up pointers to source density gradient

    auto sgrada = alphaDensityGradient(iDensityMatrix);

    auto sgradb = betaDensityGradient(iDensityMatrix);

    auto sgradab = mixedDensityGradient(iDensityMatrix);

    auto sgrada_x = alphaDensityGradientX(iDensityMatrix);

    auto sgrada_y = alphaDensityGradientY(iDensityMatrix);

    auto sgrada_z = alphaDensityGradientZ(iDensityMatrix);

    auto sgradb_x = betaDensityGradientX(iDensityMatrix);

    auto sgradb_y = betaDensityGradientY(iDensityMatrix);

    auto sgradb_z = betaDensityGradientZ(iDensityMatrix);

    // set up pointers to destination density gradient

    auto dgrada = densityGridAB.alphaDensityGradient(0);

    auto dgradb = densityGridAB.betaDensityGradient(0);

    auto dgradab = densityGridAB.mixedDensityGradient(0);

    auto dgrada_x = densityGridAB.alphaDensityGradientX(0);

    auto dgrada_y = densityGridAB.alphaDensityGradientY(0);

    auto dgrada_z = densityGridAB.alphaDensityGradientZ(0);

    auto dgradb_x = densityGridAB.betaDensityGradientX(0);

    auto dgradb_y = densityGridAB.betaDensityGradientY(0);

    auto dgradb_z = densityGridAB.betaDensityGradientZ(0);

    auto dgradbA = densityGridA.betaDensityGradient(0);

    auto dgradbA_x = densityGridA.betaDensityGradientX(0);

    auto dgradbA_y = densityGridA.betaDensityGradientY(0);

    auto dgradbA_z = densityGridA.betaDensityGradientZ(0);

    auto dgradbB = densityGridB.alphaDensityGradient(0);

    auto dgradbB_x = densityGridB.alphaDensityGradientX(0);

    auto dgradbB_y = densityGridB.alphaDensityGradientY(0);

    auto dgradbB_z = densityGridB.alphaDensityGradientZ(0);

    // density screening for GGA

    if (xcFuncType == xcfun::gga)
    {
        for (int32_t i = 0; i < npoints; i++)
        {
            if (_isValidGridPointForGgaUnrestrictedAB(srhoa[i], srhob[i], sgrada[i], sgradb[i], densityThreshold))
            {
                drhoa[ipoints] = srhoa[i];

                drhob[ipoints] = srhob[i];

                dgrada[ipoints] = sgrada[i];

                dgradb[ipoints] = sgradb[i];

                dgradab[ipoints] = sgradab[i];

                dgrada_x[ipoints] = sgrada_x[i];

                dgrada_y[ipoints] = sgrada_y[i];

                dgrada_z[ipoints] = sgrada_z[i];

                dgradb_x[ipoints] = sgradb_x[i];

                dgradb_y[ipoints] = sgradb_y[i];

                dgradb_z[ipoints] = sgradb_z[i];

                ipoints++;
            }

            if (_isValidGridPointForGgaUnrestrictedA(srhoa[i], srhob[i], sgrada[i], sgradb[i], densityThreshold))
            {
                drhobA[ipointsa] = srhob[i];

                dgradbA[ipointsa] = sgradb[i];

                dgradbA_x[ipointsa] = sgradb_x[i];

                dgradbA_y[ipointsa] = sgradb_y[i];

                dgradbA_z[ipointsa] = sgradb_z[i];

                ipointsa++;
            }

            if (_isValidGridPointForGgaUnrestrictedB(srhoa[i], srhob[i], sgrada[i], sgradb[i], densityThreshold))
            {
                drhoaB[ipointsb] = srhoa[i];

                dgradbB[ipointsb] = sgrada[i];

                dgradbB_x[ipointsb] = sgrada_x[i];

                dgradbB_y[ipointsb] = sgrada_y[i];

                dgradbB_z[ipointsb] = sgrada_z[i];

                ipointsb++;
            }
        }
    }

    // compress screened density grid size

    densityGridAB.slice(ipoints);

    densityGridA.slice(ipointsa);

    densityGridB.slice(ipointsb);
}

void
CDensityGrid::getScreenedGridUnrestricted(CMolecularGrid& molecularGridsAB,
                                          CMolecularGrid& molecularGridsA,
                                          CMolecularGrid& molecularGridsB,
                                          const int32_t   iDensityMatrix,
                                          const double    densityThreshold,
                                          const xcfun     xcFuncType) const
{
    // set up pointers to molecular grid data

    auto gx = molecularGridsAB.getCoordinatesX();

    auto gy = molecularGridsAB.getCoordinatesY();

    auto gz = molecularGridsAB.getCoordinatesZ();

    auto gw = molecularGridsAB.getWeights();

    auto gxA = molecularGridsA.getCoordinatesX();

    auto gyA = molecularGridsA.getCoordinatesY();

    auto gzA = molecularGridsA.getCoordinatesZ();

    auto gwA = molecularGridsA.getWeights();

    auto gxB = molecularGridsB.getCoordinatesX();

    auto gyB = molecularGridsB.getCoordinatesY();

    auto gzB = molecularGridsB.getCoordinatesZ();

    auto gwB = molecularGridsB.getWeights();

    // set grid points data

    auto npoints = getNumberOfGridPoints();

    int32_t ipoints = 0;

    int32_t ipointsA = 0;

    int32_t ipointsB = 0;

    // set up pointers to density data

    auto rhoa = alphaDensity(iDensityMatrix);

    auto rhob = betaDensity(iDensityMatrix);

    // screening for LDA

    if (xcFuncType == xcfun::lda)
    {
        for (int32_t i = 0; i < npoints; i++)
        {
            if (_isValidGridPointForLdaUnrestricted(rhoa[i], rhob[i], densityThreshold))
            {
                gx[ipoints] = gx[i];

                gy[ipoints] = gy[i];

                gz[ipoints] = gz[i];

                gw[ipoints] = gw[i];

                ipoints++;

                continue;
            }

            if (_isValidGridPointForLdaA(rhoa[i], rhob[i], densityThreshold))
            {
                gxA[ipointsA] = gxA[i];

                gyA[ipointsA] = gyA[i];

                gzA[ipointsA] = gzA[i];

                gwA[ipointsA] = gwA[i];

                ipointsA++;

                continue;
            }

            if (_isValidGridPointForLdaB(rhoa[i], rhob[i], densityThreshold))
            {
                gxB[ipointsB] = gxB[i];

                gyB[ipointsB] = gyB[i];

                gzB[ipointsB] = gzB[i];

                gwB[ipointsB] = gwB[i];

                ipointsB++;
            }
        }
    }

    // set up pointers to density gradient data

    auto grada = alphaDensityGradient(iDensityMatrix);

    auto gradb = betaDensityGradient(iDensityMatrix);

    // screening for GGA

    if (xcFuncType == xcfun::gga)
    {
        for (int32_t i = 0; i < npoints; i++)
        {
            if (_isValidGridPointForGgaUnrestrictedAB(rhoa[i], rhob[i], grada[i], gradb[i], densityThreshold))
            {
                gx[ipoints] = gx[i];

                gy[ipoints] = gy[i];

                gz[ipoints] = gz[i];

                gw[ipoints] = gw[i];

                ipoints++;

                continue;
            }

            if (_isValidGridPointForGgaUnrestrictedA(rhoa[i], rhob[i], grada[i], gradb[i], densityThreshold))
            {
                gxA[ipointsA] = gxA[i];

                gyA[ipointsA] = gyA[i];

                gzA[ipointsA] = gzA[i];

                gwA[ipointsA] = gwA[i];

                ipointsA++;

                continue;
            }

            if (_isValidGridPointForGgaUnrestrictedB(rhoa[i], rhob[i], grada[i], gradb[i], densityThreshold))
            {
                gxB[ipointsB] = gxB[i];

                gyB[ipointsB] = gyB[i];

                gzB[ipointsB] = gzB[i];

                gwB[ipointsB] = gwB[i];

                ipointsB++;
            }
        }
    }

    // compress molecular grid size

    molecularGridsAB.slice(ipoints);

    molecularGridsA.slice(ipointsA);

    molecularGridsB.slice(ipointsB);
}

CMolecularGrid
CDensityGrid::getScreenedGrid(CMolecularGrid& molecularGridsAB,
                              const int32_t   iDensityMatrix,
                              const double    densityThreshold,
                              const xcfun     xcFuncType) const
{
    auto mgrid = molecularGridsAB;

    // FIX ME: Implement for general case

    if (_gridType != dengrid::ab) return mgrid;

    // set up pointers to molecular grid data

    auto gx = mgrid.getCoordinatesX();

    auto gy = mgrid.getCoordinatesY();

    auto gz = mgrid.getCoordinatesZ();

    auto gw = mgrid.getWeights();

    // set grid points data

    auto npoints = getNumberOfGridPoints();

    int32_t ipoints = 0;

    // set up pointers to density data

    auto rhoa = alphaDensity(iDensityMatrix);

    auto rhob = betaDensity(iDensityMatrix);

    // screening for LDA

    if (xcFuncType == xcfun::lda)
    {
        for (int32_t i = 0; i < npoints; i++)
        {
            if (_isValidGridPointForLda(rhoa[i], rhob[i], densityThreshold))
            {
                gx[ipoints] = gx[i];
                gy[ipoints] = gy[i];
                gz[ipoints] = gz[i];
                gw[ipoints] = gw[i];

                ipoints++;
            }
        }
    }

    // set up pointers to density gradient data

    auto grada = alphaDensityGradient(iDensityMatrix);

    auto gradb = betaDensityGradient(iDensityMatrix);

    // screening for GGA

    if (xcFuncType == xcfun::gga)
    {
        for (int32_t i = 0; i < npoints; i++)
        {
            if (_isValidGridPointForGga(rhoa[i], rhob[i], grada[i], gradb[i], densityThreshold))
            {
                gx[ipoints] = gx[i];
                gy[ipoints] = gy[i];
                gz[ipoints] = gz[i];
                gw[ipoints] = gw[i];

                ipoints++;
            }
        }
    }

    auto taua = alphaDensitytau(iDensityMatrix);

    auto taub = betaDensitytau(iDensityMatrix);

    if (xcFuncType == xcfun::mgga)
    {
        for (int32_t i = 0; i < npoints; i++)
        {
            if (_isValidGridPointForMgga(rhoa[i], rhob[i], grada[i], gradb[i], taua[i], taub[i], densityThreshold))
            {
                gx[ipoints] = gx[i];
                gy[ipoints] = gy[i];
                gz[ipoints] = gz[i];
                gw[ipoints] = gw[i];

                ipoints++;
            }
        }
    }

    // compress molecular grid size

    mgrid.slice(ipoints);

    return mgrid;
}

bool
CDensityGrid::_isValidGridPointForLda(const double alphaDensity, const double betaDensity, const double densityThreshold) const
{
    if (std::fabs(alphaDensity) < densityThreshold) return false;

    if (std::fabs(betaDensity) < densityThreshold) return false;

    return true;
}

bool
CDensityGrid::_isValidGridPointForLdaUnrestricted(const double alphaDensity, const double betaDensity, const double densityThreshold) const
{
    if ((std::fabs(alphaDensity) > densityThreshold) && (std::fabs(betaDensity) > densityThreshold)) return true;

    return false;
}

bool
CDensityGrid::_isValidGridPointForLdaA(const double alphaDensity, const double betaDensity, const double densityThreshold) const
{
    if ((std::fabs(alphaDensity) < densityThreshold) && (std::fabs(betaDensity) > densityThreshold)) return true;

    return false;
}

bool
CDensityGrid::_isValidGridPointForLdaB(const double alphaDensity, const double betaDensity, const double densityThreshold) const
{
    if ((std::fabs(alphaDensity) > densityThreshold) && (std::fabs(betaDensity) < densityThreshold)) return true;

    return false;
}

bool
CDensityGrid::_isValidGridPointForGga(const double alphaDensity,
                                      const double betaDensity,
                                      const double alphaDensityGradient,
                                      const double betaDensityGradient,
                                      const double densityThreshold) const
{
    if (std::fabs(alphaDensity) < densityThreshold) return false;

    if (std::fabs(betaDensity) < densityThreshold) return false;

    if (std::fabs(alphaDensityGradient) < densityThreshold) return false;

    if (std::fabs(betaDensityGradient) < densityThreshold) return false;

    return true;
}

bool
CDensityGrid::_isValidGridPointForMgga(const double alphaDensity,
                                       const double betaDensity,
                                       const double alphaDensityGradient,
                                       const double betaDensityGradient,
                                       const double alphatau,
                                       const double betatau,
                                       const double densityThreshold) const
{
    if (std::fabs(alphaDensity) < densityThreshold) return false;

    if (std::fabs(betaDensity) < densityThreshold) return false;

    if (std::fabs(alphaDensityGradient) < densityThreshold) return false;

    if (std::fabs(betaDensityGradient) < densityThreshold) return false;

    if (std::fabs(alphatau) < densityThreshold) return false;

    if (std::fabs(betatau) < densityThreshold) return false;

    return true;
}

bool
CDensityGrid::_isValidGridPointForGgaUnrestrictedAB(const double alphaDensity,
                                                    const double betaDensity,
                                                    const double alphaDensityGradient,
                                                    const double betaDensityGradient,
                                                    const double densityThreshold) const
{
    if ((std::fabs(alphaDensity) > densityThreshold) && (std::fabs(alphaDensityGradient) > densityThreshold) &&
        (std::fabs(betaDensity) > densityThreshold) && (std::fabs(betaDensityGradient) > densityThreshold))
        return true;

    return false;
}

bool
CDensityGrid::_isValidGridPointForGgaUnrestrictedA(const double alphaDensity,
                                                   const double betaDensity,
                                                   const double alphaDensityGradient,
                                                   const double betaDensityGradient,
                                                   const double densityThreshold) const
{
    if ((std::fabs(alphaDensity) < densityThreshold) && (std::fabs(alphaDensityGradient) < densityThreshold) &&
        (std::fabs(betaDensity) > densityThreshold) && (std::fabs(betaDensityGradient) > densityThreshold))
        return true;
    return false;
}

bool
CDensityGrid::_isValidGridPointForGgaUnrestrictedB(const double alphaDensity,
                                                   const double betaDensity,
                                                   const double alphaDensityGradient,
                                                   const double betaDensityGradient,
                                                   const double densityThreshold) const
{
    if ((std::fabs(alphaDensity) > densityThreshold) && (std::fabs(alphaDensityGradient) > densityThreshold) &&
        (std::fabs(betaDensity) < densityThreshold) && (std::fabs(betaDensityGradient) < densityThreshold))
        return true;

    return false;
}

std::ostream&
operator<<(std::ostream& output, const CDensityGrid& source)
{
    output << std::endl;

    output << "[CDensityGrid (Object):" << &source << "]" << std::endl;

    output << "_gridType: " << to_string(source._gridType) << std::endl;

    output << "_densityValues: " << std::endl;

    output << source._densityValues << std::endl;

    return output;
}
