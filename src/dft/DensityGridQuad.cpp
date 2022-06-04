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

#include "DensityGridQuad.hpp"

#include <cmath>
#include <iostream>

#include "StringFormat.hpp"

CDensityGridQuad::CDensityGridQuad()

    : _gridType(dengrid::undefined)

    , _nDensityMatrices(0)

    , _densityValues(CMemBlock2D<double>())
{
}

CDensityGridQuad::CDensityGridQuad(const CMemBlock2D<double>& densityValues, const dengrid gridType)

    : _gridType(gridType)

    , _nDensityMatrices(1)

    , _densityValues(densityValues)
{
}

CDensityGridQuad::CDensityGridQuad(const int32_t nGridPoints, const int32_t nDensityMatrices, const xcfun xcFuncType, const dengrid gridType)
{
    _gridType = gridType;

    _nDensityMatrices = nDensityMatrices;

    int32_t ncomp = 0;

    if (xcFuncType == xcfun::lda) ncomp = (_gridType == dengrid::ab) ? 1 : 1;

    if (xcFuncType == xcfun::gga) ncomp = (_gridType == dengrid::ab) ? 13 : 5;

    // NOTE: this needs to be checked with mgga functionals implementation

    if (xcFuncType == xcfun::mgga) ncomp = (_gridType == dengrid::ab) ? 13 : 6;

    _densityValues = CMemBlock2D<double>(nGridPoints, _nDensityMatrices * ncomp);
}

CDensityGridQuad::CDensityGridQuad(const CDensityGridQuad& source)

    : _gridType(source._gridType)

    , _nDensityMatrices(source._nDensityMatrices)

    , _densityValues(source._densityValues)
{
}

CDensityGridQuad::CDensityGridQuad(CDensityGridQuad&& source) noexcept

    : _gridType(std::move(source._gridType))

    , _nDensityMatrices(std::move(source._nDensityMatrices))

    , _densityValues(std::move(source._densityValues))
{
}

CDensityGridQuad::~CDensityGridQuad()
{
}

CDensityGridQuad&
CDensityGridQuad::operator=(const CDensityGridQuad& source)
{
    if (this == &source) return *this;

    _gridType = source._gridType;

    _nDensityMatrices = source._nDensityMatrices;

    _densityValues = source._densityValues;

    return *this;
}

CDensityGridQuad&
CDensityGridQuad::operator=(CDensityGridQuad&& source) noexcept
{
    if (this == &source) return *this;

    _gridType = std::move(source._gridType);

    _nDensityMatrices = std::move(source._nDensityMatrices);

    _densityValues = std::move(source._densityValues);

    return *this;
}

bool
CDensityGridQuad::operator==(const CDensityGridQuad& other) const
{
    if (_gridType != other._gridType) return false;

    if (_nDensityMatrices != other._nDensityMatrices) return false;

    if (_densityValues != other._densityValues) return false;

    return true;
}

bool
CDensityGridQuad::operator!=(const CDensityGridQuad& other) const
{
    return !(*this == other);
}

void
CDensityGridQuad::zero()
{
    _densityValues.zero();
}

int32_t
CDensityGridQuad::getNumberOfGridPoints() const
{
    return _densityValues.size(0);
}

int32_t
CDensityGridQuad::getNumberOfDensityMatrices() const
{
    return _nDensityMatrices;
}

dengrid
CDensityGridQuad::getDensityGridType() const
{
    return _gridType;
}

const double*
CDensityGridQuad::rhow1rhow2(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::lima) return nullptr;

    return _densityValues.data(iDensityMatrix);
}

double*
CDensityGridQuad::rhow1rhow2(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::lima) return nullptr;

    return _densityValues.data(iDensityMatrix);
}

const double*
CDensityGridQuad::rxw1rhow2(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(1 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::rxw1rhow2(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(1 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::ryw1rhow2(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::ryw1rhow2(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::rzw1rhow2(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::rzw1rhow2(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::rxw1rxw2(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::rxw1rxw2(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::rxw1ryw2(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::rxw1ryw2(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::rxw1rzw2(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::rxw1rzw2(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::ryw1rxw2(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::ryw1rxw2(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::ryw1ryw2(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::ryw1ryw2(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::ryw1rzw2(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::ryw1rzw2(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::rzw1rxw2(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::rzw1rxw2(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::rzw1ryw2(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(11 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::rzw1ryw2(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(11 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::rzw1rzw2(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(12 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::rzw1rzw2(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(12 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

void
CDensityGridQuad::DensityProd(CDensityGridQuad&   densityGridAB,
                              CMolecularGrid&     molecularGridab,
                              const CDensityGrid& rwDensityGrid,
                              const xcfun         xcFuncType,
                              int32_t             numdens,
                              const std::string&  quadMode) const

{
    if (_gridType != dengrid::ab) return;

    // set grid points data

    auto npoints = getNumberOfGridPoints();

    // set up pointers to source density

    auto rwdenptr = &rwDensityGrid;

    // set up pointers to destination density

    if (xcFuncType == xcfun::lda)
    {
        if (fstr::upcase(quadMode) == "SHG")
        {
            for (int32_t j = 0; j < numdens / 12; j++)
            {
                // Density products to be stored

                auto rho_sig_x_r = densityGridAB.rhow1rhow2(12 * j);

                auto rho_sig_x_i = densityGridAB.rhow1rhow2(12 * j + 1);

                auto rho_sig_y_r = densityGridAB.rhow1rhow2(12 * j + 2);

                auto rho_sig_y_i = densityGridAB.rhow1rhow2(12 * j + 3);

                auto rho_sig_z_r = densityGridAB.rhow1rhow2(12 * j + 4);

                auto rho_sig_z_i = densityGridAB.rhow1rhow2(12 * j + 5);

                auto rho_lam_xy_r = densityGridAB.rhow1rhow2(12 * j + 6);

                auto rho_lam_xy_i = densityGridAB.rhow1rhow2(12 * j + 7);

                auto rho_lam_xz_r = densityGridAB.rhow1rhow2(12 * j + 8);

                auto rho_lam_xz_i = densityGridAB.rhow1rhow2(12 * j + 9);

                auto rho_lam_yz_r = densityGridAB.rhow1rhow2(12 * j + 10);

                auto rho_lam_yz_i = densityGridAB.rhow1rhow2(12 * j + 11);

                auto rhow_kx_r = rwdenptr->alphaDensity(6 * j);

                auto rhow_kx_i = rwdenptr->alphaDensity(6 * j + 1);

                auto rhow_ky_r = rwdenptr->alphaDensity(6 * j + 2);

                auto rhow_ky_i = rwdenptr->alphaDensity(6 * j + 3);

                auto rhow_kz_r = rwdenptr->alphaDensity(6 * j + 4);

                auto rhow_kz_i = rwdenptr->alphaDensity(6 * j + 5);

                for (int32_t i = 0; i < npoints; i++)
                {
                    double jj_r = 2.0 * (rhow_kx_r[i] * rhow_kx_r[i] + rhow_ky_r[i] * rhow_ky_r[i] + rhow_kz_r[i] * rhow_kz_r[i])

                                  - 2.0 * (rhow_kx_i[i] * rhow_kx_i[i] + rhow_ky_i[i] * rhow_ky_i[i] + rhow_kz_i[i] * rhow_kz_i[i]);

                    double jj_i = 2.0 * (rhow_kx_r[i] * rhow_kx_i[i] + rhow_ky_r[i] * rhow_ky_i[i] + rhow_kz_r[i] * rhow_kz_i[i]

                                         + rhow_kx_i[i] * rhow_kx_r[i] + rhow_ky_i[i] * rhow_ky_r[i] + rhow_kz_i[i] * rhow_kz_r[i]);

                    rho_sig_x_r[i] = 4.0 * (rhow_kx_r[i] * rhow_kx_r[i] - rhow_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    rho_sig_x_i[i] = 4.0 * (rhow_kx_r[i] * rhow_kx_i[i] + rhow_kx_i[i] * rhow_kx_r[i]) + jj_i;

                    rho_sig_y_r[i] = 4.0 * (rhow_ky_r[i] * rhow_ky_r[i] - rhow_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    rho_sig_y_i[i] = 4.0 * (rhow_ky_r[i] * rhow_ky_i[i] + rhow_ky_i[i] * rhow_ky_r[i]) + jj_i;

                    rho_sig_z_r[i] = 4.0 * (rhow_kz_r[i] * rhow_kz_r[i] - rhow_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    rho_sig_z_i[i] = 4.0 * (rhow_kz_r[i] * rhow_kz_i[i] + rhow_kz_i[i] * rhow_kz_r[i]) + jj_i;

                    rho_lam_xy_r[i] = (rhow_kx_r[i] * rhow_ky_r[i] - rhow_kx_i[i] * rhow_ky_i[i]

                                       + rhow_ky_r[i] * rhow_kx_r[i] - rhow_ky_i[i] * rhow_kx_i[i]);

                    rho_lam_xy_i[i] = (rhow_kx_r[i] * rhow_ky_i[i] + rhow_kx_i[i] * rhow_ky_r[i]

                                       + rhow_ky_r[i] * rhow_kx_i[i] + rhow_ky_i[i] * rhow_kx_r[i]);

                    rho_lam_xz_r[i] = (rhow_kx_r[i] * rhow_kz_r[i] - rhow_kx_i[i] * rhow_kz_i[i]

                                       + rhow_kz_r[i] * rhow_kx_r[i] - rhow_kz_i[i] * rhow_kx_i[i]);

                    rho_lam_xz_i[i] = (rhow_kx_r[i] * rhow_kz_i[i] + rhow_kx_i[i] * rhow_kz_r[i]

                                       + rhow_kz_r[i] * rhow_kx_i[i] + rhow_kz_i[i] * rhow_kx_r[i]);

                    rho_lam_yz_r[i] = (rhow_ky_r[i] * rhow_kz_r[i] - rhow_ky_i[i] * rhow_kz_i[i]

                                       + rhow_kz_r[i] * rhow_ky_r[i] - rhow_kz_i[i] * rhow_ky_i[i]);

                    rho_lam_yz_i[i] = (rhow_ky_r[i] * rhow_kz_i[i] + rhow_ky_i[i] * rhow_kz_r[i]

                                       + rhow_kz_r[i] * rhow_ky_i[i] + rhow_kz_i[i] * rhow_ky_r[i]);
                }
            }
        }
        if (fstr::upcase(quadMode) == "SHG_RED")
        {
            for (int32_t j = 0; j < numdens / 6; j++)
            {
                // Density products to be stored

                auto rho_sig_x_r = densityGridAB.rhow1rhow2(6 * j);

                auto rho_sig_y_r = densityGridAB.rhow1rhow2(6 * j + 1);

                auto rho_sig_z_r = densityGridAB.rhow1rhow2(6 * j + 2);

                auto rho_lam_xy_r = densityGridAB.rhow1rhow2(6 * j + 3);

                auto rho_lam_xz_r = densityGridAB.rhow1rhow2(6 * j + 4);

                auto rho_lam_yz_r = densityGridAB.rhow1rhow2(6 * j + 5);

                auto rhow_kx_r = rwdenptr->alphaDensity(6 * j);

                auto rhow_kx_i = rwdenptr->alphaDensity(6 * j + 1);

                auto rhow_ky_r = rwdenptr->alphaDensity(6 * j + 2);

                auto rhow_ky_i = rwdenptr->alphaDensity(6 * j + 3);

                auto rhow_kz_r = rwdenptr->alphaDensity(6 * j + 4);

                auto rhow_kz_i = rwdenptr->alphaDensity(6 * j + 5);

                for (int32_t i = 0; i < npoints; i++)
                {
                    double jj_r = 2.0 * (rhow_kx_r[i] * rhow_kx_r[i] + rhow_ky_r[i] * rhow_ky_r[i] + rhow_kz_r[i] * rhow_kz_r[i])

                                  - 2.0 * (rhow_kx_i[i] * rhow_kx_i[i] + rhow_ky_i[i] * rhow_ky_i[i] + rhow_kz_i[i] * rhow_kz_i[i]);

                    rho_sig_x_r[i] = 4.0 * (rhow_kx_r[i] * rhow_kx_r[i] - rhow_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    rho_sig_y_r[i] = 4.0 * (rhow_ky_r[i] * rhow_ky_r[i] - rhow_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    rho_sig_z_r[i] = 4.0 * (rhow_kz_r[i] * rhow_kz_r[i] - rhow_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    rho_lam_xy_r[i] = (rhow_kx_r[i] * rhow_ky_r[i] - rhow_kx_i[i] * rhow_ky_i[i]

                                       + rhow_ky_r[i] * rhow_kx_r[i] - rhow_ky_i[i] * rhow_kx_i[i]);

                    rho_lam_xz_r[i] = (rhow_kx_r[i] * rhow_kz_r[i] - rhow_kx_i[i] * rhow_kz_i[i]

                                       + rhow_kz_r[i] * rhow_kx_r[i] - rhow_kz_i[i] * rhow_kx_i[i]);

                    rho_lam_yz_r[i] = (rhow_ky_r[i] * rhow_kz_r[i] - rhow_ky_i[i] * rhow_kz_i[i]

                                       + rhow_kz_r[i] * rhow_ky_r[i] - rhow_kz_i[i] * rhow_ky_i[i]);
                }
            }
        }
        if (fstr::upcase(quadMode) == "QRF")
        {
            for (int32_t j = 0; j < numdens / 2; j++)
            {
                auto rhorho_r = densityGridAB.rhow1rhow2(2 * j);

                auto rhorho_i = densityGridAB.rhow1rhow2(2 * j + 1);

                auto rhow1a_r = rwDensityGrid.alphaDensity(4 * j);

                auto rhow1a_i = rwDensityGrid.alphaDensity(4 * j + 1);

                auto rhow2a_r = rwDensityGrid.alphaDensity(4 * j + 2);

                auto rhow2a_i = rwDensityGrid.alphaDensity(4 * j + 3);

                for (int32_t i = 0; i < npoints; i++)
                {
                    rhorho_r[i] = rhow1a_r[i] * rhow2a_r[i] - rhow1a_i[i] * rhow2a_i[i]

                                  + rhow2a_r[i] * rhow1a_r[i] - rhow2a_i[i] * rhow1a_i[i];

                    rhorho_i[i] = rhow1a_r[i] * rhow2a_i[i] + rhow1a_i[i] * rhow2a_r[i]

                                  + rhow2a_r[i] * rhow1a_i[i] + rhow2a_i[i] * rhow1a_r[i];
                }
            }
        }
    }
    if (xcFuncType == xcfun::gga)
    {
        if (fstr::upcase(quadMode) == "SHG")
        {
            for (int32_t j = 0; j < numdens / 12; j++)

            {
                // Density products to be stored

                // RhoRho part

                auto rho_sig_x_r = densityGridAB.rhow1rhow2(12 * j);

                auto rho_sig_x_i = densityGridAB.rhow1rhow2(12 * j + 1);

                auto rho_sig_y_r = densityGridAB.rhow1rhow2(12 * j + 2);

                auto rho_sig_y_i = densityGridAB.rhow1rhow2(12 * j + 3);

                auto rho_sig_z_r = densityGridAB.rhow1rhow2(12 * j + 4);

                auto rho_sig_z_i = densityGridAB.rhow1rhow2(12 * j + 5);

                auto rho_lam_xy_r = densityGridAB.rhow1rhow2(12 * j + 6);

                auto rho_lam_xy_i = densityGridAB.rhow1rhow2(12 * j + 7);

                auto rho_lam_xz_r = densityGridAB.rhow1rhow2(12 * j + 8);

                auto rho_lam_xz_i = densityGridAB.rhow1rhow2(12 * j + 9);

                auto rho_lam_yz_r = densityGridAB.rhow1rhow2(12 * j + 10);

                auto rho_lam_yz_i = densityGridAB.rhow1rhow2(12 * j + 11);

                // rxw1rhow part

                auto rxw1rhow2_sig_x_r = densityGridAB.rxw1rhow2(12 * j);

                auto rxw1rhow2_sig_x_i = densityGridAB.rxw1rhow2(12 * j + 1);

                auto rxw1rhow2_sig_y_r = densityGridAB.rxw1rhow2(12 * j + 2);

                auto rxw1rhow2_sig_y_i = densityGridAB.rxw1rhow2(12 * j + 3);

                auto rxw1rhow2_sig_z_r = densityGridAB.rxw1rhow2(12 * j + 4);

                auto rxw1rhow2_sig_z_i = densityGridAB.rxw1rhow2(12 * j + 5);

                auto rxw1rhow2_lam_xy_r = densityGridAB.rxw1rhow2(12 * j + 6);

                auto rxw1rhow2_lam_xy_i = densityGridAB.rxw1rhow2(12 * j + 7);

                auto rxw1rhow2_lam_xz_r = densityGridAB.rxw1rhow2(12 * j + 8);

                auto rxw1rhow2_lam_xz_i = densityGridAB.rxw1rhow2(12 * j + 9);

                auto rxw1rhow2_lam_yz_r = densityGridAB.rxw1rhow2(12 * j + 10);

                auto rxw1rhow2_lam_yz_i = densityGridAB.rxw1rhow2(12 * j + 11);

                // ryw1rhow part

                auto ryw1rhow2_sig_x_r = densityGridAB.ryw1rhow2(12 * j);

                auto ryw1rhow2_sig_x_i = densityGridAB.ryw1rhow2(12 * j + 1);

                auto ryw1rhow2_sig_y_r = densityGridAB.ryw1rhow2(12 * j + 2);

                auto ryw1rhow2_sig_y_i = densityGridAB.ryw1rhow2(12 * j + 3);

                auto ryw1rhow2_sig_z_r = densityGridAB.ryw1rhow2(12 * j + 4);

                auto ryw1rhow2_sig_z_i = densityGridAB.ryw1rhow2(12 * j + 5);

                auto ryw1rhow2_lam_xy_r = densityGridAB.ryw1rhow2(12 * j + 6);

                auto ryw1rhow2_lam_xy_i = densityGridAB.ryw1rhow2(12 * j + 7);

                auto ryw1rhow2_lam_xz_r = densityGridAB.ryw1rhow2(12 * j + 8);

                auto ryw1rhow2_lam_xz_i = densityGridAB.ryw1rhow2(12 * j + 9);

                auto ryw1rhow2_lam_yz_r = densityGridAB.ryw1rhow2(12 * j + 10);

                auto ryw1rhow2_lam_yz_i = densityGridAB.ryw1rhow2(12 * j + 11);

                // rzw1rhow part

                auto rzw1rhow2_sig_x_r = densityGridAB.rzw1rhow2(12 * j);

                auto rzw1rhow2_sig_x_i = densityGridAB.rzw1rhow2(12 * j + 1);

                auto rzw1rhow2_sig_y_r = densityGridAB.rzw1rhow2(12 * j + 2);

                auto rzw1rhow2_sig_y_i = densityGridAB.rzw1rhow2(12 * j + 3);

                auto rzw1rhow2_sig_z_r = densityGridAB.rzw1rhow2(12 * j + 4);

                auto rzw1rhow2_sig_z_i = densityGridAB.rzw1rhow2(12 * j + 5);

                auto rzw1rhow2_lam_xy_r = densityGridAB.rzw1rhow2(12 * j + 6);

                auto rzw1rhow2_lam_xy_i = densityGridAB.rzw1rhow2(12 * j + 7);

                auto rzw1rhow2_lam_xz_r = densityGridAB.rzw1rhow2(12 * j + 8);

                auto rzw1rhow2_lam_xz_i = densityGridAB.rzw1rhow2(12 * j + 9);

                auto rzw1rhow2_lam_yz_r = densityGridAB.rzw1rhow2(12 * j + 10);

                auto rzw1rhow2_lam_yz_i = densityGridAB.rzw1rhow2(12 * j + 11);

                // rxw1rxw2 part

                auto rxw1rxw2_sig_x_r = densityGridAB.rxw1rxw2(12 * j);

                auto rxw1rxw2_sig_x_i = densityGridAB.rxw1rxw2(12 * j + 1);

                auto rxw1rxw2_sig_y_r = densityGridAB.rxw1rxw2(12 * j + 2);

                auto rxw1rxw2_sig_y_i = densityGridAB.rxw1rxw2(12 * j + 3);

                auto rxw1rxw2_sig_z_r = densityGridAB.rxw1rxw2(12 * j + 4);

                auto rxw1rxw2_sig_z_i = densityGridAB.rxw1rxw2(12 * j + 5);

                auto rxw1rxw2_lam_xy_r = densityGridAB.rxw1rxw2(12 * j + 6);

                auto rxw1rxw2_lam_xy_i = densityGridAB.rxw1rxw2(12 * j + 7);

                auto rxw1rxw2_lam_xz_r = densityGridAB.rxw1rxw2(12 * j + 8);

                auto rxw1rxw2_lam_xz_i = densityGridAB.rxw1rxw2(12 * j + 9);

                auto rxw1rxw2_lam_yz_r = densityGridAB.rxw1rxw2(12 * j + 10);

                auto rxw1rxw2_lam_yz_i = densityGridAB.rxw1rxw2(12 * j + 11);

                // rxw1ryw2 part

                auto rxw1ryw2_sig_x_r = densityGridAB.rxw1ryw2(12 * j);

                auto rxw1ryw2_sig_x_i = densityGridAB.rxw1ryw2(12 * j + 1);

                auto rxw1ryw2_sig_y_r = densityGridAB.rxw1ryw2(12 * j + 2);

                auto rxw1ryw2_sig_y_i = densityGridAB.rxw1ryw2(12 * j + 3);

                auto rxw1ryw2_sig_z_r = densityGridAB.rxw1ryw2(12 * j + 4);

                auto rxw1ryw2_sig_z_i = densityGridAB.rxw1ryw2(12 * j + 5);

                auto rxw1ryw2_lam_xy_r = densityGridAB.rxw1ryw2(12 * j + 6);

                auto rxw1ryw2_lam_xy_i = densityGridAB.rxw1ryw2(12 * j + 7);

                auto rxw1ryw2_lam_xz_r = densityGridAB.rxw1ryw2(12 * j + 8);

                auto rxw1ryw2_lam_xz_i = densityGridAB.rxw1ryw2(12 * j + 9);

                auto rxw1ryw2_lam_yz_r = densityGridAB.rxw1ryw2(12 * j + 10);

                auto rxw1ryw2_lam_yz_i = densityGridAB.rxw1ryw2(12 * j + 11);

                // rxw1rzw2 part

                auto rxw1rzw2_sig_x_r = densityGridAB.rxw1rzw2(12 * j);

                auto rxw1rzw2_sig_x_i = densityGridAB.rxw1rzw2(12 * j + 1);

                auto rxw1rzw2_sig_y_r = densityGridAB.rxw1rzw2(12 * j + 2);

                auto rxw1rzw2_sig_y_i = densityGridAB.rxw1rzw2(12 * j + 3);

                auto rxw1rzw2_sig_z_r = densityGridAB.rxw1rzw2(12 * j + 4);

                auto rxw1rzw2_sig_z_i = densityGridAB.rxw1rzw2(12 * j + 5);

                auto rxw1rzw2_lam_xy_r = densityGridAB.rxw1rzw2(12 * j + 6);

                auto rxw1rzw2_lam_xy_i = densityGridAB.rxw1rzw2(12 * j + 7);

                auto rxw1rzw2_lam_xz_r = densityGridAB.rxw1rzw2(12 * j + 8);

                auto rxw1rzw2_lam_xz_i = densityGridAB.rxw1rzw2(12 * j + 9);

                auto rxw1rzw2_lam_yz_r = densityGridAB.rxw1rzw2(12 * j + 10);

                auto rxw1rzw2_lam_yz_i = densityGridAB.rxw1rzw2(12 * j + 11);

                // ryw1rxw2 part

                auto ryw1rxw2_sig_x_r = densityGridAB.ryw1rxw2(12 * j);

                auto ryw1rxw2_sig_x_i = densityGridAB.ryw1rxw2(12 * j + 1);

                auto ryw1rxw2_sig_y_r = densityGridAB.ryw1rxw2(12 * j + 2);

                auto ryw1rxw2_sig_y_i = densityGridAB.ryw1rxw2(12 * j + 3);

                auto ryw1rxw2_sig_z_r = densityGridAB.ryw1rxw2(12 * j + 4);

                auto ryw1rxw2_sig_z_i = densityGridAB.ryw1rxw2(12 * j + 5);

                auto ryw1rxw2_lam_xy_r = densityGridAB.ryw1rxw2(12 * j + 6);

                auto ryw1rxw2_lam_xy_i = densityGridAB.ryw1rxw2(12 * j + 7);

                auto ryw1rxw2_lam_xz_r = densityGridAB.ryw1rxw2(12 * j + 8);

                auto ryw1rxw2_lam_xz_i = densityGridAB.ryw1rxw2(12 * j + 9);

                auto ryw1rxw2_lam_yz_r = densityGridAB.ryw1rxw2(12 * j + 10);

                auto ryw1rxw2_lam_yz_i = densityGridAB.ryw1rxw2(12 * j + 11);

                // ryw1ryw2 part

                auto ryw1ryw2_sig_x_r = densityGridAB.ryw1ryw2(12 * j);

                auto ryw1ryw2_sig_x_i = densityGridAB.ryw1ryw2(12 * j + 1);

                auto ryw1ryw2_sig_y_r = densityGridAB.ryw1ryw2(12 * j + 2);

                auto ryw1ryw2_sig_y_i = densityGridAB.ryw1ryw2(12 * j + 3);

                auto ryw1ryw2_sig_z_r = densityGridAB.ryw1ryw2(12 * j + 4);

                auto ryw1ryw2_sig_z_i = densityGridAB.ryw1ryw2(12 * j + 5);

                auto ryw1ryw2_lam_xy_r = densityGridAB.ryw1ryw2(12 * j + 6);

                auto ryw1ryw2_lam_xy_i = densityGridAB.ryw1ryw2(12 * j + 7);

                auto ryw1ryw2_lam_xz_r = densityGridAB.ryw1ryw2(12 * j + 8);

                auto ryw1ryw2_lam_xz_i = densityGridAB.ryw1ryw2(12 * j + 9);

                auto ryw1ryw2_lam_yz_r = densityGridAB.ryw1ryw2(12 * j + 10);

                auto ryw1ryw2_lam_yz_i = densityGridAB.ryw1ryw2(12 * j + 11);

                // ryw1rzw2 part

                auto ryw1rzw2_sig_x_r = densityGridAB.ryw1rzw2(12 * j);

                auto ryw1rzw2_sig_x_i = densityGridAB.ryw1rzw2(12 * j + 1);

                auto ryw1rzw2_sig_y_r = densityGridAB.ryw1rzw2(12 * j + 2);

                auto ryw1rzw2_sig_y_i = densityGridAB.ryw1rzw2(12 * j + 3);

                auto ryw1rzw2_sig_z_r = densityGridAB.ryw1rzw2(12 * j + 4);

                auto ryw1rzw2_sig_z_i = densityGridAB.ryw1rzw2(12 * j + 5);

                auto ryw1rzw2_lam_xy_r = densityGridAB.ryw1rzw2(12 * j + 6);

                auto ryw1rzw2_lam_xy_i = densityGridAB.ryw1rzw2(12 * j + 7);

                auto ryw1rzw2_lam_xz_r = densityGridAB.ryw1rzw2(12 * j + 8);

                auto ryw1rzw2_lam_xz_i = densityGridAB.ryw1rzw2(12 * j + 9);

                auto ryw1rzw2_lam_yz_r = densityGridAB.ryw1rzw2(12 * j + 10);

                auto ryw1rzw2_lam_yz_i = densityGridAB.ryw1rzw2(12 * j + 11);

                // rzw1rxw2 part

                auto rzw1rxw2_sig_x_r = densityGridAB.rzw1rxw2(12 * j);

                auto rzw1rxw2_sig_x_i = densityGridAB.rzw1rxw2(12 * j + 1);

                auto rzw1rxw2_sig_y_r = densityGridAB.rzw1rxw2(12 * j + 2);

                auto rzw1rxw2_sig_y_i = densityGridAB.rzw1rxw2(12 * j + 3);

                auto rzw1rxw2_sig_z_r = densityGridAB.rzw1rxw2(12 * j + 4);

                auto rzw1rxw2_sig_z_i = densityGridAB.rzw1rxw2(12 * j + 5);

                auto rzw1rxw2_lam_xy_r = densityGridAB.rzw1rxw2(12 * j + 6);

                auto rzw1rxw2_lam_xy_i = densityGridAB.rzw1rxw2(12 * j + 7);

                auto rzw1rxw2_lam_xz_r = densityGridAB.rzw1rxw2(12 * j + 8);

                auto rzw1rxw2_lam_xz_i = densityGridAB.rzw1rxw2(12 * j + 9);

                auto rzw1rxw2_lam_yz_r = densityGridAB.rzw1rxw2(12 * j + 10);

                auto rzw1rxw2_lam_yz_i = densityGridAB.rzw1rxw2(12 * j + 11);

                // rzw1ryw2 part

                auto rzw1ryw2_sig_x_r = densityGridAB.rzw1ryw2(12 * j);

                auto rzw1ryw2_sig_x_i = densityGridAB.rzw1ryw2(12 * j + 1);

                auto rzw1ryw2_sig_y_r = densityGridAB.rzw1ryw2(12 * j + 2);

                auto rzw1ryw2_sig_y_i = densityGridAB.rzw1ryw2(12 * j + 3);

                auto rzw1ryw2_sig_z_r = densityGridAB.rzw1ryw2(12 * j + 4);

                auto rzw1ryw2_sig_z_i = densityGridAB.rzw1ryw2(12 * j + 5);

                auto rzw1ryw2_lam_xy_r = densityGridAB.rzw1ryw2(12 * j + 6);

                auto rzw1ryw2_lam_xy_i = densityGridAB.rzw1ryw2(12 * j + 7);

                auto rzw1ryw2_lam_xz_r = densityGridAB.rzw1ryw2(12 * j + 8);

                auto rzw1ryw2_lam_xz_i = densityGridAB.rzw1ryw2(12 * j + 9);

                auto rzw1ryw2_lam_yz_r = densityGridAB.rzw1ryw2(12 * j + 10);

                auto rzw1ryw2_lam_yz_i = densityGridAB.rzw1ryw2(12 * j + 11);

                // rzw1rzw2 part

                auto rzw1rzw2_sig_x_r = densityGridAB.rzw1rzw2(12 * j);

                auto rzw1rzw2_sig_x_i = densityGridAB.rzw1rzw2(12 * j + 1);

                auto rzw1rzw2_sig_y_r = densityGridAB.rzw1rzw2(12 * j + 2);

                auto rzw1rzw2_sig_y_i = densityGridAB.rzw1rzw2(12 * j + 3);

                auto rzw1rzw2_sig_z_r = densityGridAB.rzw1rzw2(12 * j + 4);

                auto rzw1rzw2_sig_z_i = densityGridAB.rzw1rzw2(12 * j + 5);

                auto rzw1rzw2_lam_xy_r = densityGridAB.rzw1rzw2(12 * j + 6);

                auto rzw1rzw2_lam_xy_i = densityGridAB.rzw1rzw2(12 * j + 7);

                auto rzw1rzw2_lam_xz_r = densityGridAB.rzw1rzw2(12 * j + 8);

                auto rzw1rzw2_lam_xz_i = densityGridAB.rzw1rzw2(12 * j + 9);

                auto rzw1rzw2_lam_yz_r = densityGridAB.rzw1rzw2(12 * j + 10);

                auto rzw1rzw2_lam_yz_i = densityGridAB.rzw1rzw2(12 * j + 11);

                // First-order densities

                // Rho terms

                auto rhow_kx_r = rwdenptr->alphaDensity(6 * j);

                auto rhow_kx_i = rwdenptr->alphaDensity(6 * j + 1);

                auto rhow_ky_r = rwdenptr->alphaDensity(6 * j + 2);

                auto rhow_ky_i = rwdenptr->alphaDensity(6 * j + 3);

                auto rhow_kz_r = rwdenptr->alphaDensity(6 * j + 4);

                auto rhow_kz_i = rwdenptr->alphaDensity(6 * j + 5);

                // Gradx terms

                auto gradxw_kx_r = rwDensityGrid.alphaDensityGradientX(6 * j);

                auto gradxw_kx_i = rwDensityGrid.alphaDensityGradientX(6 * j + 1);

                auto gradxw_ky_r = rwDensityGrid.alphaDensityGradientX(6 * j + 2);

                auto gradxw_ky_i = rwDensityGrid.alphaDensityGradientX(6 * j + 3);

                auto gradxw_kz_r = rwDensityGrid.alphaDensityGradientX(6 * j + 4);

                auto gradxw_kz_i = rwDensityGrid.alphaDensityGradientX(6 * j + 5);

                // Grady terms

                auto gradyw_kx_r = rwDensityGrid.alphaDensityGradientY(6 * j);

                auto gradyw_kx_i = rwDensityGrid.alphaDensityGradientY(6 * j + 1);

                auto gradyw_ky_r = rwDensityGrid.alphaDensityGradientY(6 * j + 2);

                auto gradyw_ky_i = rwDensityGrid.alphaDensityGradientY(6 * j + 3);

                auto gradyw_kz_r = rwDensityGrid.alphaDensityGradientY(6 * j + 4);

                auto gradyw_kz_i = rwDensityGrid.alphaDensityGradientY(6 * j + 5);

                // Gradz terms

                auto gradzw_kx_r = rwDensityGrid.alphaDensityGradientZ(6 * j);

                auto gradzw_kx_i = rwDensityGrid.alphaDensityGradientZ(6 * j + 1);

                auto gradzw_ky_r = rwDensityGrid.alphaDensityGradientZ(6 * j + 2);

                auto gradzw_ky_i = rwDensityGrid.alphaDensityGradientZ(6 * j + 3);

                auto gradzw_kz_r = rwDensityGrid.alphaDensityGradientZ(6 * j + 4);

                auto gradzw_kz_i = rwDensityGrid.alphaDensityGradientZ(6 * j + 5);

                for (int32_t i = 0; i < npoints; i++)
                {
                    // RhoRho real

                    double jj_r = 2.0 * (rhow_kx_r[i] * rhow_kx_r[i] + rhow_ky_r[i] * rhow_ky_r[i] + rhow_kz_r[i] * rhow_kz_r[i])

                                  - 2.0 * (rhow_kx_i[i] * rhow_kx_i[i] + rhow_ky_i[i] * rhow_ky_i[i] + rhow_kz_i[i] * rhow_kz_i[i]);

                    rho_sig_x_r[i] = 4.0 * (rhow_kx_r[i] * rhow_kx_r[i] - rhow_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    rho_sig_y_r[i] = 4.0 * (rhow_ky_r[i] * rhow_ky_r[i] - rhow_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    rho_sig_z_r[i] = 4.0 * (rhow_kz_r[i] * rhow_kz_r[i] - rhow_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    rho_lam_xy_r[i] = (rhow_kx_r[i] * rhow_ky_r[i] - rhow_kx_i[i] * rhow_ky_i[i]

                                       + rhow_ky_r[i] * rhow_kx_r[i] - rhow_ky_i[i] * rhow_kx_i[i]);

                    rho_lam_xz_r[i] = (rhow_kx_r[i] * rhow_kz_r[i] - rhow_kx_i[i] * rhow_kz_i[i]

                                       + rhow_kz_r[i] * rhow_kx_r[i] - rhow_kz_i[i] * rhow_kx_i[i]);

                    rho_lam_yz_r[i] = (rhow_ky_r[i] * rhow_kz_r[i] - rhow_ky_i[i] * rhow_kz_i[i]

                                       + rhow_kz_r[i] * rhow_ky_r[i] - rhow_kz_i[i] * rhow_ky_i[i]);

                    // RhoRho imag

                    double jj_i = 2.0 * (rhow_kx_r[i] * rhow_kx_i[i] + rhow_ky_r[i] * rhow_ky_i[i] + rhow_kz_r[i] * rhow_kz_i[i]

                                         + rhow_kx_i[i] * rhow_kx_r[i] + rhow_ky_i[i] * rhow_ky_r[i] + rhow_kz_i[i] * rhow_kz_r[i]);

                    rho_sig_x_i[i] = 4.0 * (rhow_kx_r[i] * rhow_kx_i[i] + rhow_kx_i[i] * rhow_kx_r[i]) + jj_i;

                    rho_sig_y_i[i] = 4.0 * (rhow_ky_r[i] * rhow_ky_i[i] + rhow_ky_i[i] * rhow_ky_r[i]) + jj_i;

                    rho_sig_z_i[i] = 4.0 * (rhow_kz_r[i] * rhow_kz_i[i] + rhow_kz_i[i] * rhow_kz_r[i]) + jj_i;

                    rho_lam_xy_i[i] = (rhow_kx_r[i] * rhow_ky_i[i] + rhow_kx_i[i] * rhow_ky_r[i]

                                       + rhow_ky_r[i] * rhow_kx_i[i] + rhow_ky_i[i] * rhow_kx_r[i]);

                    rho_lam_xz_i[i] = (rhow_kx_r[i] * rhow_kz_i[i] + rhow_kx_i[i] * rhow_kz_r[i]

                                       + rhow_kz_r[i] * rhow_kx_i[i] + rhow_kz_i[i] * rhow_kx_r[i]);

                    rho_lam_yz_i[i] = (rhow_ky_r[i] * rhow_kz_i[i] + rhow_ky_i[i] * rhow_kz_r[i]

                                       + rhow_kz_r[i] * rhow_ky_i[i] + rhow_kz_i[i] * rhow_ky_r[i]);

                    // rxw1rhow2 real

                    jj_r = 4.0 * (gradxw_kx_r[i] * rhow_kx_r[i] + gradxw_ky_r[i] * rhow_ky_r[i] + gradxw_kz_r[i] * rhow_kz_r[i])

                           - 4.0 * (gradxw_kx_i[i] * rhow_kx_i[i] + gradxw_ky_i[i] * rhow_ky_i[i] + gradxw_kz_i[i] * rhow_kz_i[i]);

                    rxw1rhow2_sig_x_r[i] = 8.0 * (gradxw_kx_r[i] * rhow_kx_r[i] - gradxw_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    rxw1rhow2_sig_y_r[i] = 8.0 * (gradxw_ky_r[i] * rhow_ky_r[i] - gradxw_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    rxw1rhow2_sig_z_r[i] = 8.0 * (gradxw_kz_r[i] * rhow_kz_r[i] - gradxw_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    rxw1rhow2_lam_xy_r[i] = 2.0 * (gradxw_kx_r[i] * rhow_ky_r[i] - gradxw_kx_i[i] * rhow_ky_i[i]

                                                   + gradxw_ky_r[i] * rhow_kx_r[i] - gradxw_ky_i[i] * rhow_kx_i[i]);

                    rxw1rhow2_lam_xz_r[i] = 2.0 * (gradxw_kx_r[i] * rhow_kz_r[i] - gradxw_kx_i[i] * rhow_kz_i[i]

                                                   + gradxw_kz_r[i] * rhow_kx_r[i] - gradxw_kz_i[i] * rhow_kx_i[i]);

                    rxw1rhow2_lam_yz_r[i] = 2.0 * (gradxw_ky_r[i] * rhow_kz_r[i] - gradxw_ky_i[i] * rhow_kz_i[i]

                                                   + gradxw_kz_r[i] * rhow_ky_r[i] - gradxw_kz_i[i] * rhow_ky_i[i]);

                    // rxw1rhow2 imag

                    jj_i = 4.0 * (gradxw_kx_r[i] * rhow_kx_i[i] + gradxw_ky_r[i] * rhow_ky_i[i] + gradxw_kz_r[i] * rhow_kz_i[i]

                                  + gradxw_kx_i[i] * rhow_kx_r[i] + gradxw_ky_i[i] * rhow_ky_r[i] + gradxw_kz_i[i] * rhow_kz_r[i]);

                    rxw1rhow2_sig_x_i[i] = 8.0 * (gradxw_kx_r[i] * rhow_kx_i[i] + gradxw_kx_i[i] * rhow_kx_r[i]) + jj_i;

                    rxw1rhow2_sig_y_i[i] = 8.0 * (gradxw_ky_r[i] * rhow_ky_i[i] + gradxw_ky_i[i] * rhow_ky_r[i]) + jj_i;

                    rxw1rhow2_sig_z_i[i] = 8.0 * (gradxw_kz_r[i] * rhow_kz_i[i] + gradxw_kz_i[i] * rhow_kz_r[i]) + jj_i;

                    rxw1rhow2_lam_xy_i[i] = 2.0 * (gradxw_kx_r[i] * rhow_ky_i[i] + gradxw_kx_i[i] * rhow_ky_r[i]

                                                   + gradxw_ky_r[i] * rhow_kx_i[i] + gradxw_ky_i[i] * rhow_kx_r[i]);

                    rxw1rhow2_lam_xz_i[i] = 2.0 * (gradxw_kx_r[i] * rhow_kz_i[i] + gradxw_kx_i[i] * rhow_kz_r[i]

                                                   + gradxw_kz_r[i] * rhow_kx_i[i] + gradxw_kz_i[i] * rhow_kx_r[i]);

                    rxw1rhow2_lam_yz_i[i] = 2.0 * (gradxw_ky_r[i] * rhow_kz_i[i] + gradxw_ky_i[i] * rhow_kz_r[i]

                                                   + gradxw_kz_r[i] * rhow_ky_i[i] + gradxw_kz_i[i] * rhow_ky_r[i]);

                    // ryw1rhow2 real

                    jj_r = 4.0 * (gradyw_kx_r[i] * rhow_kx_r[i] + gradyw_ky_r[i] * rhow_ky_r[i] + gradyw_kz_r[i] * rhow_kz_r[i])

                           - 4.0 * (gradyw_kx_i[i] * rhow_kx_i[i] + gradyw_ky_i[i] * rhow_ky_i[i] + gradyw_kz_i[i] * rhow_kz_i[i]);

                    ryw1rhow2_sig_x_r[i] = 8.0 * (gradyw_kx_r[i] * rhow_kx_r[i] - gradyw_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    ryw1rhow2_sig_y_r[i] = 8.0 * (gradyw_ky_r[i] * rhow_ky_r[i] - gradyw_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    ryw1rhow2_sig_z_r[i] = 8.0 * (gradyw_kz_r[i] * rhow_kz_r[i] - gradyw_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    ryw1rhow2_lam_xy_r[i] = 2.0 * (gradyw_kx_r[i] * rhow_ky_r[i] - gradyw_kx_i[i] * rhow_ky_i[i]

                                                   + gradyw_ky_r[i] * rhow_kx_r[i] - gradyw_ky_i[i] * rhow_kx_i[i]);

                    ryw1rhow2_lam_xz_r[i] = 2.0 * (gradyw_kx_r[i] * rhow_kz_r[i] - gradyw_kx_i[i] * rhow_kz_i[i]

                                                   + gradyw_kz_r[i] * rhow_kx_r[i] - gradyw_kz_i[i] * rhow_kx_i[i]);

                    ryw1rhow2_lam_yz_r[i] = 2.0 * (gradyw_ky_r[i] * rhow_kz_r[i] - gradyw_ky_i[i] * rhow_kz_i[i]

                                                   + gradyw_kz_r[i] * rhow_ky_r[i] - gradyw_kz_i[i] * rhow_ky_i[i]);

                    // ryw1rhow2 imag

                    jj_i = 4.0 * (gradyw_kx_r[i] * rhow_kx_i[i] + gradyw_ky_r[i] * rhow_ky_i[i] + gradyw_kz_r[i] * rhow_kz_i[i]

                                  + gradyw_kx_i[i] * rhow_kx_r[i] + gradyw_ky_i[i] * rhow_ky_r[i] + gradyw_kz_i[i] * rhow_kz_r[i]);

                    ryw1rhow2_sig_x_i[i] = 8.0 * (gradyw_kx_r[i] * rhow_kx_i[i] + gradyw_kx_i[i] * rhow_kx_r[i]) + jj_i;

                    ryw1rhow2_sig_y_i[i] = 8.0 * (gradyw_ky_r[i] * rhow_ky_i[i] + gradyw_ky_i[i] * rhow_ky_r[i]) + jj_i;

                    ryw1rhow2_sig_z_i[i] = 8.0 * (gradyw_kz_r[i] * rhow_kz_i[i] + gradyw_kz_i[i] * rhow_kz_r[i]) + jj_i;

                    ryw1rhow2_lam_xy_i[i] = 2.0 * (gradyw_kx_r[i] * rhow_ky_i[i] + gradyw_kx_i[i] * rhow_ky_r[i]

                                                   + gradyw_ky_r[i] * rhow_kx_i[i] + gradyw_ky_i[i] * rhow_kx_r[i]);

                    ryw1rhow2_lam_xz_i[i] = 2.0 * (gradyw_kx_r[i] * rhow_kz_i[i] + gradyw_kx_i[i] * rhow_kz_r[i]

                                                   + gradyw_kz_r[i] * rhow_kx_i[i] + gradyw_kz_i[i] * rhow_kx_r[i]);

                    ryw1rhow2_lam_yz_i[i] = 2.0 * (gradyw_ky_r[i] * rhow_kz_i[i] + gradyw_ky_i[i] * rhow_kz_r[i]

                                                   + gradyw_kz_r[i] * rhow_ky_i[i] + gradyw_kz_i[i] * rhow_ky_r[i]);

                    // rzw1rhow2 real

                    jj_r = 4.0 * (gradzw_kx_r[i] * rhow_kx_r[i] + gradzw_ky_r[i] * rhow_ky_r[i] + gradzw_kz_r[i] * rhow_kz_r[i])

                           - 4.0 * (gradzw_kx_i[i] * rhow_kx_i[i] + gradzw_ky_i[i] * rhow_ky_i[i] + gradzw_kz_i[i] * rhow_kz_i[i]);

                    rzw1rhow2_sig_x_r[i] = 8.0 * (gradzw_kx_r[i] * rhow_kx_r[i] - gradzw_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    rzw1rhow2_sig_y_r[i] = 8.0 * (gradzw_ky_r[i] * rhow_ky_r[i] - gradzw_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    rzw1rhow2_sig_z_r[i] = 8.0 * (gradzw_kz_r[i] * rhow_kz_r[i] - gradzw_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    rzw1rhow2_lam_xy_r[i] = 2.0 * (gradzw_kx_r[i] * rhow_ky_r[i] - gradzw_kx_i[i] * rhow_ky_i[i]

                                                   + gradzw_ky_r[i] * rhow_kx_r[i] - gradzw_ky_i[i] * rhow_kx_i[i]);

                    rzw1rhow2_lam_xz_r[i] = 2.0 * (gradzw_kx_r[i] * rhow_kz_r[i] - gradzw_kx_i[i] * rhow_kz_i[i]

                                                   + gradzw_kz_r[i] * rhow_kx_r[i] - gradzw_kz_i[i] * rhow_kx_i[i]);

                    rzw1rhow2_lam_yz_r[i] = 2.0 * (gradzw_ky_r[i] * rhow_kz_r[i] - gradzw_ky_i[i] * rhow_kz_i[i]

                                                   + gradzw_kz_r[i] * rhow_ky_r[i] - gradzw_kz_i[i] * rhow_ky_i[i]);

                    // rzw1rhow2 imag

                    jj_i = 4.0 * (gradzw_kx_r[i] * rhow_kx_i[i] + gradzw_ky_r[i] * rhow_ky_i[i] + gradzw_kz_r[i] * rhow_kz_i[i]

                                  + gradzw_kx_i[i] * rhow_kx_r[i] + gradzw_ky_i[i] * rhow_ky_r[i] + gradzw_kz_i[i] * rhow_kz_r[i]);

                    rzw1rhow2_sig_x_i[i] = 8.0 * (gradzw_kx_r[i] * rhow_kx_i[i] + gradzw_kx_i[i] * rhow_kx_r[i]) + jj_i;

                    rzw1rhow2_sig_y_i[i] = 8.0 * (gradzw_ky_r[i] * rhow_ky_i[i] + gradzw_ky_i[i] * rhow_ky_r[i]) + jj_i;

                    rzw1rhow2_sig_z_i[i] = 8.0 * (gradzw_kz_r[i] * rhow_kz_i[i] + gradzw_kz_i[i] * rhow_kz_r[i]) + jj_i;

                    rzw1rhow2_lam_xy_i[i] = 2.0 * (gradzw_kx_r[i] * rhow_ky_i[i] + gradzw_kx_i[i] * rhow_ky_r[i]

                                                   + gradzw_ky_r[i] * rhow_kx_i[i] + gradzw_ky_i[i] * rhow_kx_r[i]);

                    rzw1rhow2_lam_xz_i[i] = 2.0 * (gradzw_kx_r[i] * rhow_kz_i[i] + gradzw_kx_i[i] * rhow_kz_r[i]

                                                   + gradzw_kz_r[i] * rhow_kx_i[i] + gradzw_kz_i[i] * rhow_kx_r[i]);

                    rzw1rhow2_lam_yz_i[i] = 2.0 * (gradzw_ky_r[i] * rhow_kz_i[i] + gradzw_ky_i[i] * rhow_kz_r[i]

                                                   + gradzw_kz_r[i] * rhow_ky_i[i] + gradzw_kz_i[i] * rhow_ky_r[i]);

                    // rxw1rxw2

                    jj_r = 2.0 * (gradxw_kx_r[i] * gradxw_kx_r[i] + gradxw_ky_r[i] * gradxw_ky_r[i] + gradxw_kz_r[i] * gradxw_kz_r[i])

                           - 2.0 * (gradxw_kx_i[i] * gradxw_kx_i[i] + gradxw_ky_i[i] * gradxw_ky_i[i] + gradxw_kz_i[i] * gradxw_kz_i[i]);

                    rxw1rxw2_sig_x_r[i] = 4.0 * (gradxw_kx_r[i] * gradxw_kx_r[i] - gradxw_kx_i[i] * gradxw_kx_i[i]) + jj_r;

                    rxw1rxw2_sig_y_r[i] = 4.0 * (gradxw_ky_r[i] * gradxw_ky_r[i] - gradxw_ky_i[i] * gradxw_ky_i[i]) + jj_r;

                    rxw1rxw2_sig_z_r[i] = 4.0 * (gradxw_kz_r[i] * gradxw_kz_r[i] - gradxw_kz_i[i] * gradxw_kz_i[i]) + jj_r;

                    rxw1rxw2_lam_xy_r[i] = gradxw_kx_r[i] * gradxw_ky_r[i] - gradxw_kx_i[i] * gradxw_ky_i[i]

                                           + gradxw_ky_r[i] * gradxw_kx_r[i] - gradxw_ky_i[i] * gradxw_kx_i[i];

                    rxw1rxw2_lam_xz_r[i] = gradxw_kx_r[i] * gradxw_kz_r[i] - gradxw_kx_i[i] * gradxw_kz_i[i]

                                           + gradxw_kz_r[i] * gradxw_kx_r[i] - gradxw_kz_i[i] * gradxw_kx_i[i];

                    rxw1rxw2_lam_yz_r[i] = gradxw_ky_r[i] * gradxw_kz_r[i] - gradxw_ky_i[i] * gradxw_kz_i[i]

                                           + gradxw_kz_r[i] * gradxw_ky_r[i] - gradxw_kz_i[i] * gradxw_ky_i[i];

                    // rxw1rxw2 imag

                    jj_i = 2.0 * (gradxw_kx_r[i] * gradxw_kx_i[i] + gradxw_ky_r[i] * gradxw_ky_i[i] + gradxw_kz_r[i] * gradxw_kz_i[i]

                                  + gradxw_kx_i[i] * gradxw_kx_r[i] + gradxw_ky_i[i] * gradxw_ky_r[i] + gradxw_kz_i[i] * gradxw_kz_r[i]);

                    rxw1rxw2_sig_x_i[i] = 4.0 * (gradxw_kx_r[i] * gradxw_kx_i[i] + gradxw_kx_i[i] * gradxw_kx_r[i]) + jj_i;

                    rxw1rxw2_sig_y_i[i] = 4.0 * (gradxw_ky_r[i] * gradxw_ky_i[i] + gradxw_ky_i[i] * gradxw_ky_r[i]) + jj_i;

                    rxw1rxw2_sig_z_i[i] = 4.0 * (gradxw_kz_r[i] * gradxw_kz_i[i] + gradxw_kz_i[i] * gradxw_kz_r[i]) + jj_i;

                    rxw1rxw2_lam_xy_i[i] = gradxw_kx_r[i] * gradxw_ky_i[i] + gradxw_kx_i[i] * gradxw_ky_r[i]

                                           + gradxw_ky_r[i] * gradxw_kx_i[i] + gradxw_ky_i[i] * gradxw_kx_r[i];

                    rxw1rxw2_lam_xz_i[i] = gradxw_kx_r[i] * gradxw_kz_i[i] + gradxw_kx_i[i] * gradxw_kz_r[i]

                                           + gradxw_kz_r[i] * gradxw_kx_i[i] + gradxw_kz_i[i] * gradxw_kx_r[i];

                    rxw1rxw2_lam_yz_i[i] = gradxw_ky_r[i] * gradxw_kz_i[i] + gradxw_ky_i[i] * gradxw_kz_r[i]

                                           + gradxw_kz_r[i] * gradxw_ky_i[i] + gradxw_kz_i[i] * gradxw_ky_r[i];

                    // rxw1ryw2

                    jj_r = 2.0 * (gradxw_kx_r[i] * gradyw_kx_r[i] + gradxw_ky_r[i] * gradyw_ky_r[i] + gradxw_kz_r[i] * gradyw_kz_r[i])

                           - 2.0 * (gradxw_kx_i[i] * gradyw_kx_i[i] + gradxw_ky_i[i] * gradyw_ky_i[i] + gradxw_kz_i[i] * gradyw_kz_i[i]);

                    rxw1ryw2_sig_x_r[i] = 4.0 * (gradxw_kx_r[i] * gradyw_kx_r[i] - gradxw_kx_i[i] * gradyw_kx_i[i]) + jj_r;

                    rxw1ryw2_sig_y_r[i] = 4.0 * (gradxw_ky_r[i] * gradyw_ky_r[i] - gradxw_ky_i[i] * gradyw_ky_i[i]) + jj_r;

                    rxw1ryw2_sig_z_r[i] = 4.0 * (gradxw_kz_r[i] * gradyw_kz_r[i] - gradxw_kz_i[i] * gradyw_kz_i[i]) + jj_r;

                    rxw1ryw2_lam_xy_r[i] = gradxw_kx_r[i] * gradyw_ky_r[i] - gradxw_kx_i[i] * gradyw_ky_i[i]

                                           + gradxw_ky_r[i] * gradyw_kx_r[i] - gradxw_ky_i[i] * gradyw_kx_i[i];

                    rxw1ryw2_lam_xz_r[i] = gradxw_kx_r[i] * gradyw_kz_r[i] - gradxw_kx_i[i] * gradyw_kz_i[i]

                                           + gradxw_kz_r[i] * gradyw_kx_r[i] - gradxw_kz_i[i] * gradyw_kx_i[i];

                    rxw1ryw2_lam_yz_r[i] = gradxw_ky_r[i] * gradyw_kz_r[i] - gradxw_ky_i[i] * gradyw_kz_i[i]

                                           + gradxw_kz_r[i] * gradyw_ky_r[i] - gradxw_kz_i[i] * gradyw_ky_i[i];

                    // rxw1ryw2 imag

                    jj_i = 2.0 * (gradxw_kx_r[i] * gradyw_kx_i[i] + gradxw_ky_r[i] * gradyw_ky_i[i] + gradxw_kz_r[i] * gradyw_kz_i[i]

                                  + gradxw_kx_i[i] * gradyw_kx_r[i] + gradxw_ky_i[i] * gradyw_ky_r[i] + gradxw_kz_i[i] * gradyw_kz_r[i]);

                    rxw1ryw2_sig_x_i[i] = 4.0 * (gradxw_kx_r[i] * gradyw_kx_i[i] + gradxw_kx_i[i] * gradyw_kx_r[i]) + jj_i;

                    rxw1ryw2_sig_y_i[i] = 4.0 * (gradxw_ky_r[i] * gradyw_ky_i[i] + gradxw_ky_i[i] * gradyw_ky_r[i]) + jj_i;

                    rxw1ryw2_sig_z_i[i] = 4.0 * (gradxw_kz_r[i] * gradyw_kz_i[i] + gradxw_kz_i[i] * gradyw_kz_r[i]) + jj_i;

                    rxw1ryw2_lam_xy_i[i] = gradxw_kx_r[i] * gradyw_ky_i[i] + gradxw_kx_i[i] * gradyw_ky_r[i]

                                           + gradxw_ky_r[i] * gradyw_kx_i[i] + gradxw_ky_i[i] * gradyw_kx_r[i];

                    rxw1ryw2_lam_xz_i[i] = gradxw_kx_r[i] * gradyw_kz_i[i] + gradxw_kx_i[i] * gradyw_kz_r[i]

                                           + gradxw_kz_r[i] * gradyw_kx_i[i] + gradxw_kz_i[i] * gradyw_kx_r[i];

                    rxw1ryw2_lam_yz_i[i] = gradxw_ky_r[i] * gradyw_kz_i[i] + gradxw_ky_i[i] * gradyw_kz_r[i]

                                           + gradxw_kz_r[i] * gradyw_ky_i[i] + gradxw_kz_i[i] * gradyw_ky_r[i];

                    // rxw1rzw2

                    jj_r = 2.0 * (gradxw_kx_r[i] * gradzw_kx_r[i] + gradxw_ky_r[i] * gradzw_ky_r[i] + gradxw_kz_r[i] * gradzw_kz_r[i])

                           - 2.0 * (gradxw_kx_i[i] * gradzw_kx_i[i] + gradxw_ky_i[i] * gradzw_ky_i[i] + gradxw_kz_i[i] * gradzw_kz_i[i]);

                    rxw1rzw2_sig_x_r[i] = 4.0 * (gradxw_kx_r[i] * gradzw_kx_r[i] - gradxw_kx_i[i] * gradzw_kx_i[i]) + jj_r;

                    rxw1rzw2_sig_y_r[i] = 4.0 * (gradxw_ky_r[i] * gradzw_ky_r[i] - gradxw_ky_i[i] * gradzw_ky_i[i]) + jj_r;

                    rxw1rzw2_sig_z_r[i] = 4.0 * (gradxw_kz_r[i] * gradzw_kz_r[i] - gradxw_kz_i[i] * gradzw_kz_i[i]) + jj_r;

                    rxw1rzw2_lam_xy_r[i] = gradxw_kx_r[i] * gradzw_ky_r[i] - gradxw_kx_i[i] * gradzw_ky_i[i]

                                           + gradxw_ky_r[i] * gradzw_kx_r[i] - gradxw_ky_i[i] * gradzw_kx_i[i];

                    rxw1rzw2_lam_xz_r[i] = gradxw_kx_r[i] * gradzw_kz_r[i] - gradxw_kx_i[i] * gradzw_kz_i[i]

                                           + gradxw_kz_r[i] * gradzw_kx_r[i] - gradxw_kz_i[i] * gradzw_kx_i[i];

                    rxw1rzw2_lam_yz_r[i] = gradxw_ky_r[i] * gradzw_kz_r[i] - gradxw_ky_i[i] * gradzw_kz_i[i]

                                           + gradxw_kz_r[i] * gradzw_ky_r[i] - gradxw_kz_i[i] * gradzw_ky_i[i];

                    // rxw1rzw2 imag

                    jj_i = 2.0 * (gradxw_kx_r[i] * gradzw_kx_i[i] + gradxw_ky_r[i] * gradzw_ky_i[i] + gradxw_kz_r[i] * gradzw_kz_i[i]

                                  + gradxw_kx_i[i] * gradzw_kx_r[i] + gradxw_ky_i[i] * gradzw_ky_r[i] + gradxw_kz_i[i] * gradzw_kz_r[i]);

                    rxw1rzw2_sig_x_i[i] = 4.0 * (gradxw_kx_r[i] * gradzw_kx_i[i] + gradxw_kx_i[i] * gradzw_kx_r[i]) + jj_i;

                    rxw1rzw2_sig_y_i[i] = 4.0 * (gradxw_ky_r[i] * gradzw_ky_i[i] + gradxw_ky_i[i] * gradzw_ky_r[i]) + jj_i;

                    rxw1rzw2_sig_z_i[i] = 4.0 * (gradxw_kz_r[i] * gradzw_kz_i[i] + gradxw_kz_i[i] * gradzw_kz_r[i]) + jj_i;

                    rxw1rzw2_lam_xy_i[i] = gradxw_kx_r[i] * gradzw_ky_i[i] + gradxw_kx_i[i] * gradzw_ky_r[i]

                                           + gradxw_ky_r[i] * gradzw_kx_i[i] + gradxw_ky_i[i] * gradzw_kx_r[i];

                    rxw1rzw2_lam_xz_i[i] = gradxw_kx_r[i] * gradzw_kz_i[i] + gradxw_kx_i[i] * gradzw_kz_r[i]

                                           + gradxw_kz_r[i] * gradzw_kx_i[i] + gradxw_kz_i[i] * gradzw_kx_r[i];

                    rxw1rzw2_lam_yz_i[i] = gradxw_ky_r[i] * gradzw_kz_i[i] + gradxw_ky_i[i] * gradzw_kz_r[i]

                                           + gradxw_kz_r[i] * gradzw_ky_i[i] + gradxw_kz_i[i] * gradzw_ky_r[i];

                    // ryw1rxw2

                    ryw1rxw2_sig_x_r[i] = rxw1ryw2_sig_x_r[i];

                    ryw1rxw2_sig_y_r[i] = rxw1ryw2_sig_y_r[i];

                    ryw1rxw2_sig_z_r[i] = rxw1ryw2_sig_z_r[i];

                    ryw1rxw2_lam_xy_r[i] = rxw1ryw2_lam_xy_r[i];

                    ryw1rxw2_lam_xz_r[i] = rxw1ryw2_lam_xz_r[i];

                    ryw1rxw2_lam_yz_r[i] = rxw1ryw2_lam_yz_r[i];

                    // ryw1rxw2 imag

                    ryw1rxw2_sig_x_i[i] = rxw1ryw2_sig_x_i[i];

                    ryw1rxw2_sig_y_i[i] = rxw1ryw2_sig_y_i[i];

                    ryw1rxw2_sig_z_i[i] = rxw1ryw2_sig_z_i[i];

                    ryw1rxw2_lam_xy_i[i] = rxw1ryw2_lam_xy_i[i];

                    ryw1rxw2_lam_xz_i[i] = rxw1ryw2_lam_xz_i[i];

                    ryw1rxw2_lam_yz_i[i] = rxw1ryw2_lam_yz_i[i];

                    // ryw1ryw2

                    jj_r = 2.0 * (gradyw_kx_r[i] * gradyw_kx_r[i] + gradyw_ky_r[i] * gradyw_ky_r[i] + gradyw_kz_r[i] * gradyw_kz_r[i])

                           - 2.0 * (gradyw_kx_i[i] * gradyw_kx_i[i] + gradyw_ky_i[i] * gradyw_ky_i[i] + gradyw_kz_i[i] * gradyw_kz_i[i]);

                    ryw1ryw2_sig_x_r[i] = 4.0 * (gradyw_kx_r[i] * gradyw_kx_r[i] - gradyw_kx_i[i] * gradyw_kx_i[i]) + jj_r;

                    ryw1ryw2_sig_y_r[i] = 4.0 * (gradyw_ky_r[i] * gradyw_ky_r[i] - gradyw_ky_i[i] * gradyw_ky_i[i]) + jj_r;

                    ryw1ryw2_sig_z_r[i] = 4.0 * (gradyw_kz_r[i] * gradyw_kz_r[i] - gradyw_kz_i[i] * gradyw_kz_i[i]) + jj_r;

                    ryw1ryw2_lam_xy_r[i] = gradyw_kx_r[i] * gradyw_ky_r[i] - gradyw_kx_i[i] * gradyw_ky_i[i]

                                           + gradyw_ky_r[i] * gradyw_kx_r[i] - gradyw_ky_i[i] * gradyw_kx_i[i];

                    ryw1ryw2_lam_xz_r[i] = gradyw_kx_r[i] * gradyw_kz_r[i] - gradyw_kx_i[i] * gradyw_kz_i[i]

                                           + gradyw_kz_r[i] * gradyw_kx_r[i] - gradyw_kz_i[i] * gradyw_kx_i[i];

                    ryw1ryw2_lam_yz_r[i] = gradyw_ky_r[i] * gradyw_kz_r[i] - gradyw_ky_i[i] * gradyw_kz_i[i]

                                           + gradyw_kz_r[i] * gradyw_ky_r[i] - gradyw_kz_i[i] * gradyw_ky_i[i];

                    // ryw1ryw2 imag

                    jj_i = 2.0 * (gradyw_kx_r[i] * gradyw_kx_i[i] + gradyw_ky_r[i] * gradyw_ky_i[i] + gradyw_kz_r[i] * gradyw_kz_i[i]

                                  + gradyw_kx_i[i] * gradyw_kx_r[i] + gradyw_ky_i[i] * gradyw_ky_r[i] + gradyw_kz_i[i] * gradyw_kz_r[i]);

                    ryw1ryw2_sig_x_i[i] = 4.0 * (gradyw_kx_r[i] * gradyw_kx_i[i] + gradyw_kx_i[i] * gradyw_kx_r[i]) + jj_i;

                    ryw1ryw2_sig_y_i[i] = 4.0 * (gradyw_ky_r[i] * gradyw_ky_i[i] + gradyw_ky_i[i] * gradyw_ky_r[i]) + jj_i;

                    ryw1ryw2_sig_z_i[i] = 4.0 * (gradyw_kz_r[i] * gradyw_kz_i[i] + gradyw_kz_i[i] * gradyw_kz_r[i]) + jj_i;

                    ryw1ryw2_lam_xy_i[i] = gradyw_kx_r[i] * gradyw_ky_i[i] + gradyw_kx_i[i] * gradyw_ky_r[i]

                                           + gradyw_ky_r[i] * gradyw_kx_i[i] + gradyw_ky_i[i] * gradyw_kx_r[i];

                    ryw1ryw2_lam_xz_i[i] = gradyw_kx_r[i] * gradyw_kz_i[i] + gradyw_kx_i[i] * gradyw_kz_r[i]

                                           + gradyw_kz_r[i] * gradyw_kx_i[i] + gradyw_kz_i[i] * gradyw_kx_r[i];

                    ryw1ryw2_lam_yz_i[i] = gradyw_ky_r[i] * gradyw_kz_i[i] + gradyw_ky_i[i] * gradyw_kz_r[i]

                                           + gradyw_kz_r[i] * gradyw_ky_i[i] + gradyw_kz_i[i] * gradyw_ky_r[i];

                    // ryw1rzw2

                    jj_r = 2.0 * (gradyw_kx_r[i] * gradzw_kx_r[i] + gradyw_ky_r[i] * gradzw_ky_r[i] + gradyw_kz_r[i] * gradzw_kz_r[i])

                           - 2.0 * (gradyw_kx_i[i] * gradzw_kx_i[i] + gradyw_ky_i[i] * gradzw_ky_i[i] + gradyw_kz_i[i] * gradzw_kz_i[i]);

                    ryw1rzw2_sig_x_r[i] = 4.0 * (gradyw_kx_r[i] * gradzw_kx_r[i] - gradyw_kx_i[i] * gradzw_kx_i[i]) + jj_r;

                    ryw1rzw2_sig_y_r[i] = 4.0 * (gradyw_ky_r[i] * gradzw_ky_r[i] - gradyw_ky_i[i] * gradzw_ky_i[i]) + jj_r;

                    ryw1rzw2_sig_z_r[i] = 4.0 * (gradyw_kz_r[i] * gradzw_kz_r[i] - gradyw_kz_i[i] * gradzw_kz_i[i]) + jj_r;

                    ryw1rzw2_lam_xy_r[i] = gradyw_kx_r[i] * gradzw_ky_r[i] - gradyw_kx_i[i] * gradzw_ky_i[i]

                                           + gradyw_ky_r[i] * gradzw_kx_r[i] - gradyw_ky_i[i] * gradzw_kx_i[i];

                    ryw1rzw2_lam_xz_r[i] = gradyw_kx_r[i] * gradzw_kz_r[i] - gradyw_kx_i[i] * gradzw_kz_i[i]

                                           + gradyw_kz_r[i] * gradzw_kx_r[i] - gradyw_kz_i[i] * gradzw_kx_i[i];

                    ryw1rzw2_lam_yz_r[i] = gradyw_ky_r[i] * gradzw_kz_r[i] - gradyw_ky_i[i] * gradzw_kz_i[i]

                                           + gradyw_kz_r[i] * gradzw_ky_r[i] - gradyw_kz_i[i] * gradzw_ky_i[i];

                    // ryw1rzw2 imag

                    jj_i = 2.0 * (gradyw_kx_r[i] * gradzw_kx_i[i] + gradyw_ky_r[i] * gradzw_ky_i[i] + gradyw_kz_r[i] * gradzw_kz_i[i]

                                  + gradyw_kx_i[i] * gradzw_kx_r[i] + gradyw_ky_i[i] * gradzw_ky_r[i] + gradyw_kz_i[i] * gradzw_kz_r[i]);

                    ryw1rzw2_sig_x_i[i] = 4.0 * (gradyw_kx_r[i] * gradzw_kx_i[i] + gradyw_kx_i[i] * gradzw_kx_r[i]) + jj_i;

                    ryw1rzw2_sig_y_i[i] = 4.0 * (gradyw_ky_r[i] * gradzw_ky_i[i] + gradyw_ky_i[i] * gradzw_ky_r[i]) + jj_i;

                    ryw1rzw2_sig_z_i[i] = 4.0 * (gradyw_kz_r[i] * gradzw_kz_i[i] + gradyw_kz_i[i] * gradzw_kz_r[i]) + jj_i;

                    ryw1rzw2_lam_xy_i[i] = gradyw_kx_r[i] * gradzw_ky_i[i] + gradyw_kx_i[i] * gradzw_ky_r[i]

                                           + gradyw_ky_r[i] * gradzw_kx_i[i] + gradyw_ky_i[i] * gradzw_kx_r[i];

                    ryw1rzw2_lam_xz_i[i] = gradyw_kx_r[i] * gradzw_kz_i[i] + gradyw_kx_i[i] * gradzw_kz_r[i]

                                           + gradyw_kz_r[i] * gradzw_kx_i[i] + gradyw_kz_i[i] * gradzw_kx_r[i];

                    ryw1rzw2_lam_yz_i[i] = gradyw_ky_r[i] * gradzw_kz_i[i] + gradyw_ky_i[i] * gradzw_kz_r[i]

                                           + gradyw_kz_r[i] * gradzw_ky_i[i] + gradyw_kz_i[i] * gradzw_ky_r[i];

                    // rzw1rxw2

                    rzw1rxw2_sig_x_r[i] = rxw1rzw2_sig_x_r[i];

                    rzw1rxw2_sig_y_r[i] = rxw1rzw2_sig_y_r[i];

                    rzw1rxw2_sig_z_r[i] = rxw1rzw2_sig_z_r[i];

                    rzw1rxw2_lam_xy_r[i] = rxw1rzw2_lam_xy_r[i];

                    rzw1rxw2_lam_xz_r[i] = rxw1rzw2_lam_xz_r[i];

                    rzw1rxw2_lam_yz_r[i] = rxw1rzw2_lam_yz_r[i];

                    // rzw1rxw2 imag

                    rzw1rxw2_sig_x_i[i] = rxw1rzw2_sig_x_i[i];

                    rzw1rxw2_sig_y_i[i] = rxw1rzw2_sig_y_i[i];

                    rzw1rxw2_sig_z_i[i] = rxw1rzw2_sig_z_i[i];

                    rzw1rxw2_lam_xy_i[i] = rxw1rzw2_lam_xy_i[i];

                    rzw1rxw2_lam_xz_i[i] = rxw1rzw2_lam_xz_i[i];

                    rzw1rxw2_lam_yz_i[i] = rxw1rzw2_lam_yz_i[i];

                    // rzw1ryw2

                    rzw1ryw2_sig_x_r[i] = ryw1rzw2_sig_x_r[i];

                    rzw1ryw2_sig_y_r[i] = ryw1rzw2_sig_y_r[i];

                    rzw1ryw2_sig_z_r[i] = ryw1rzw2_sig_z_r[i];

                    rzw1ryw2_lam_xy_r[i] = ryw1rzw2_lam_xy_r[i];

                    rzw1ryw2_lam_xz_r[i] = ryw1rzw2_lam_xz_r[i];

                    rzw1ryw2_lam_yz_r[i] = ryw1rzw2_lam_yz_r[i];

                    // rzw1ryw2 imag

                    rzw1ryw2_sig_x_i[i] = ryw1rzw2_sig_x_i[i];

                    rzw1ryw2_sig_y_i[i] = ryw1rzw2_sig_y_i[i];

                    rzw1ryw2_sig_z_i[i] = ryw1rzw2_sig_z_i[i];

                    rzw1ryw2_lam_xy_i[i] = ryw1rzw2_lam_xy_i[i];

                    rzw1ryw2_lam_xz_i[i] = ryw1rzw2_lam_xz_i[i];

                    rzw1ryw2_lam_yz_i[i] = ryw1rzw2_lam_yz_i[i];

                    // rzw1rzw2

                    jj_r = 2.0 * (gradzw_kx_r[i] * gradzw_kx_r[i] + gradzw_ky_r[i] * gradzw_ky_r[i] + gradzw_kz_r[i] * gradzw_kz_r[i])

                           - 2.0 * (gradzw_kx_i[i] * gradzw_kx_i[i] + gradzw_ky_i[i] * gradzw_ky_i[i] + gradzw_kz_i[i] * gradzw_kz_i[i]);

                    rzw1rzw2_sig_x_r[i] = 4.0 * (gradzw_kx_r[i] * gradzw_kx_r[i] - gradzw_kx_i[i] * gradzw_kx_i[i]) + jj_r;

                    rzw1rzw2_sig_y_r[i] = 4.0 * (gradzw_ky_r[i] * gradzw_ky_r[i] - gradzw_ky_i[i] * gradzw_ky_i[i]) + jj_r;

                    rzw1rzw2_sig_z_r[i] = 4.0 * (gradzw_kz_r[i] * gradzw_kz_r[i] - gradzw_kz_i[i] * gradzw_kz_i[i]) + jj_r;

                    rzw1rzw2_lam_xy_r[i] = gradzw_kx_r[i] * gradzw_ky_r[i] - gradzw_kx_i[i] * gradzw_ky_i[i]

                                           + gradzw_ky_r[i] * gradzw_kx_r[i] - gradzw_ky_i[i] * gradzw_kx_i[i];

                    rzw1rzw2_lam_xz_r[i] = gradzw_kx_r[i] * gradzw_kz_r[i] - gradzw_kx_i[i] * gradzw_kz_i[i]

                                           + gradzw_kz_r[i] * gradzw_kx_r[i] - gradzw_kz_i[i] * gradzw_kx_i[i];

                    rzw1rzw2_lam_yz_r[i] = gradzw_ky_r[i] * gradzw_kz_r[i] - gradzw_ky_i[i] * gradzw_kz_i[i]

                                           + gradzw_kz_r[i] * gradzw_ky_r[i] - gradzw_kz_i[i] * gradzw_ky_i[i];

                    // rzw1rzw2 imag

                    jj_i = 2.0 * (gradzw_kx_r[i] * gradzw_kx_i[i] + gradzw_ky_r[i] * gradzw_ky_i[i] + gradzw_kz_r[i] * gradzw_kz_i[i]

                                  + gradzw_kx_i[i] * gradzw_kx_r[i] + gradzw_ky_i[i] * gradzw_ky_r[i] + gradzw_kz_i[i] * gradzw_kz_r[i]);

                    rzw1rzw2_sig_x_i[i] = 4.0 * (gradzw_kx_r[i] * gradzw_kx_i[i] + gradzw_kx_i[i] * gradzw_kx_r[i]) + jj_i;

                    rzw1rzw2_sig_y_i[i] = 4.0 * (gradzw_ky_r[i] * gradzw_ky_i[i] + gradzw_ky_i[i] * gradzw_ky_r[i]) + jj_i;

                    rzw1rzw2_sig_z_i[i] = 4.0 * (gradzw_kz_r[i] * gradzw_kz_i[i] + gradzw_kz_i[i] * gradzw_kz_r[i]) + jj_i;

                    rzw1rzw2_lam_xy_i[i] = gradzw_kx_r[i] * gradzw_ky_i[i] + gradzw_kx_i[i] * gradzw_ky_r[i]

                                           + gradzw_ky_r[i] * gradzw_kx_i[i] + gradzw_ky_i[i] * gradzw_kx_r[i];

                    rzw1rzw2_lam_xz_i[i] = gradzw_kx_r[i] * gradzw_kz_i[i] + gradzw_kx_i[i] * gradzw_kz_r[i]

                                           + gradzw_kz_r[i] * gradzw_kx_i[i] + gradzw_kz_i[i] * gradzw_kx_r[i];

                    rzw1rzw2_lam_yz_i[i] = gradzw_ky_r[i] * gradzw_kz_i[i] + gradzw_ky_i[i] * gradzw_kz_r[i]

                                           + gradzw_kz_r[i] * gradzw_ky_i[i] + gradzw_kz_i[i] * gradzw_ky_r[i];
                }
            }
        }
        if (fstr::upcase(quadMode) == "SHG_RED")
        {
            for (int32_t j = 0; j < numdens / 6; j++)

            {
                // Density products to be stored

                // RhoRho part

                auto rho_sig_x_r = densityGridAB.rhow1rhow2(6 * j);

                auto rho_sig_y_r = densityGridAB.rhow1rhow2(6 * j + 1);

                auto rho_sig_z_r = densityGridAB.rhow1rhow2(6 * j + 2);

                auto rho_lam_xy_r = densityGridAB.rhow1rhow2(6 * j + 3);

                auto rho_lam_xz_r = densityGridAB.rhow1rhow2(6 * j + 4);

                auto rho_lam_yz_r = densityGridAB.rhow1rhow2(6 * j + 5);

                // rxw1rhow part

                auto rxw1rhow2_sig_x_r = densityGridAB.rxw1rhow2(6 * j);

                auto rxw1rhow2_sig_y_r = densityGridAB.rxw1rhow2(6 * j + 1);

                auto rxw1rhow2_sig_z_r = densityGridAB.rxw1rhow2(6 * j + 2);

                auto rxw1rhow2_lam_xy_r = densityGridAB.rxw1rhow2(6 * j + 3);

                auto rxw1rhow2_lam_xz_r = densityGridAB.rxw1rhow2(6 * j + 4);

                auto rxw1rhow2_lam_yz_r = densityGridAB.rxw1rhow2(6 * j + 5);

                // ryw1rhow part

                auto ryw1rhow2_sig_x_r = densityGridAB.ryw1rhow2(6 * j);

                auto ryw1rhow2_sig_y_r = densityGridAB.ryw1rhow2(6 * j + 1);

                auto ryw1rhow2_sig_z_r = densityGridAB.ryw1rhow2(6 * j + 2);

                auto ryw1rhow2_lam_xy_r = densityGridAB.ryw1rhow2(6 * j + 3);

                auto ryw1rhow2_lam_xz_r = densityGridAB.ryw1rhow2(6 * j + 4);

                auto ryw1rhow2_lam_yz_r = densityGridAB.ryw1rhow2(6 * j + 5);

                // rzw1rhow part

                auto rzw1rhow2_sig_x_r = densityGridAB.rzw1rhow2(6 * j);

                auto rzw1rhow2_sig_y_r = densityGridAB.rzw1rhow2(6 * j + 1);

                auto rzw1rhow2_sig_z_r = densityGridAB.rzw1rhow2(6 * j + 2);

                auto rzw1rhow2_lam_xy_r = densityGridAB.rzw1rhow2(6 * j + 3);

                auto rzw1rhow2_lam_xz_r = densityGridAB.rzw1rhow2(6 * j + 4);

                auto rzw1rhow2_lam_yz_r = densityGridAB.rzw1rhow2(6 * j + 5);

                // rxw1rxw2 part

                auto rxw1rxw2_sig_x_r = densityGridAB.rxw1rxw2(6 * j);

                auto rxw1rxw2_sig_y_r = densityGridAB.rxw1rxw2(6 * j + 1);

                auto rxw1rxw2_sig_z_r = densityGridAB.rxw1rxw2(6 * j + 2);

                auto rxw1rxw2_lam_xy_r = densityGridAB.rxw1rxw2(6 * j + 3);

                auto rxw1rxw2_lam_xz_r = densityGridAB.rxw1rxw2(6 * j + 4);

                auto rxw1rxw2_lam_yz_r = densityGridAB.rxw1rxw2(6 * j + 5);

                // rxw1ryw2 part

                auto rxw1ryw2_sig_x_r = densityGridAB.rxw1ryw2(6 * j);

                auto rxw1ryw2_sig_y_r = densityGridAB.rxw1ryw2(6 * j + 1);

                auto rxw1ryw2_sig_z_r = densityGridAB.rxw1ryw2(6 * j + 2);

                auto rxw1ryw2_lam_xy_r = densityGridAB.rxw1ryw2(6 * j + 3);

                auto rxw1ryw2_lam_xz_r = densityGridAB.rxw1ryw2(6 * j + 4);

                auto rxw1ryw2_lam_yz_r = densityGridAB.rxw1ryw2(6 * j + 5);

                // rxw1rzw2 part

                auto rxw1rzw2_sig_x_r = densityGridAB.rxw1rzw2(6 * j);

                auto rxw1rzw2_sig_y_r = densityGridAB.rxw1rzw2(6 * j + 1);

                auto rxw1rzw2_sig_z_r = densityGridAB.rxw1rzw2(6 * j + 2);

                auto rxw1rzw2_lam_xy_r = densityGridAB.rxw1rzw2(6 * j + 3);

                auto rxw1rzw2_lam_xz_r = densityGridAB.rxw1rzw2(6 * j + 4);

                auto rxw1rzw2_lam_yz_r = densityGridAB.rxw1rzw2(6 * j + 5);

                // ryw1rxw2 part

                auto ryw1rxw2_sig_x_r = densityGridAB.ryw1rxw2(6 * j);

                auto ryw1rxw2_sig_y_r = densityGridAB.ryw1rxw2(6 * j + 1);

                auto ryw1rxw2_sig_z_r = densityGridAB.ryw1rxw2(6 * j + 2);

                auto ryw1rxw2_lam_xy_r = densityGridAB.ryw1rxw2(6 * j + 3);

                auto ryw1rxw2_lam_xz_r = densityGridAB.ryw1rxw2(6 * j + 4);

                auto ryw1rxw2_lam_yz_r = densityGridAB.ryw1rxw2(6 * j + 5);

                // ryw1ryw2 part

                auto ryw1ryw2_sig_x_r = densityGridAB.ryw1ryw2(6 * j);

                auto ryw1ryw2_sig_y_r = densityGridAB.ryw1ryw2(6 * j + 1);

                auto ryw1ryw2_sig_z_r = densityGridAB.ryw1ryw2(6 * j + 2);

                auto ryw1ryw2_lam_xy_r = densityGridAB.ryw1ryw2(6 * j + 3);

                auto ryw1ryw2_lam_xz_r = densityGridAB.ryw1ryw2(6 * j + 4);

                auto ryw1ryw2_lam_yz_r = densityGridAB.ryw1ryw2(6 * j + 5);

                // ryw1rzw2 part

                auto ryw1rzw2_sig_x_r = densityGridAB.ryw1rzw2(6 * j);

                auto ryw1rzw2_sig_y_r = densityGridAB.ryw1rzw2(6 * j + 1);

                auto ryw1rzw2_sig_z_r = densityGridAB.ryw1rzw2(6 * j + 2);

                auto ryw1rzw2_lam_xy_r = densityGridAB.ryw1rzw2(6 * j + 3);

                auto ryw1rzw2_lam_xz_r = densityGridAB.ryw1rzw2(6 * j + 4);

                auto ryw1rzw2_lam_yz_r = densityGridAB.ryw1rzw2(6 * j + 5);

                // rzw1rxw2 part

                auto rzw1rxw2_sig_x_r = densityGridAB.rzw1rxw2(6 * j);

                auto rzw1rxw2_sig_y_r = densityGridAB.rzw1rxw2(6 * j + 1);

                auto rzw1rxw2_sig_z_r = densityGridAB.rzw1rxw2(6 * j + 2);

                auto rzw1rxw2_lam_xy_r = densityGridAB.rzw1rxw2(6 * j + 3);

                auto rzw1rxw2_lam_xz_r = densityGridAB.rzw1rxw2(6 * j + 4);

                auto rzw1rxw2_lam_yz_r = densityGridAB.rzw1rxw2(6 * j + 5);

                // rzw1ryw2 part

                auto rzw1ryw2_sig_x_r = densityGridAB.rzw1ryw2(6 * j);

                auto rzw1ryw2_sig_y_r = densityGridAB.rzw1ryw2(6 * j + 1);

                auto rzw1ryw2_sig_z_r = densityGridAB.rzw1ryw2(6 * j + 2);

                auto rzw1ryw2_lam_xy_r = densityGridAB.rzw1ryw2(6 * j + 3);

                auto rzw1ryw2_lam_xz_r = densityGridAB.rzw1ryw2(6 * j + 4);

                auto rzw1ryw2_lam_yz_r = densityGridAB.rzw1ryw2(6 * j + 5);

                // rzw1rzw2 part

                auto rzw1rzw2_sig_x_r = densityGridAB.rzw1rzw2(6 * j);

                auto rzw1rzw2_sig_y_r = densityGridAB.rzw1rzw2(6 * j + 1);

                auto rzw1rzw2_sig_z_r = densityGridAB.rzw1rzw2(6 * j + 2);

                auto rzw1rzw2_lam_xy_r = densityGridAB.rzw1rzw2(6 * j + 3);

                auto rzw1rzw2_lam_xz_r = densityGridAB.rzw1rzw2(6 * j + 4);

                auto rzw1rzw2_lam_yz_r = densityGridAB.rzw1rzw2(6 * j + 5);

                // First-order densities

                // Rho terms

                auto rhow_kx_r = rwdenptr->alphaDensity(6 * j);

                auto rhow_kx_i = rwdenptr->alphaDensity(6 * j + 1);

                auto rhow_ky_r = rwdenptr->alphaDensity(6 * j + 2);

                auto rhow_ky_i = rwdenptr->alphaDensity(6 * j + 3);

                auto rhow_kz_r = rwdenptr->alphaDensity(6 * j + 4);

                auto rhow_kz_i = rwdenptr->alphaDensity(6 * j + 5);

                // Gradx terms

                auto gradxw_kx_r = rwDensityGrid.alphaDensityGradientX(6 * j);

                auto gradxw_kx_i = rwDensityGrid.alphaDensityGradientX(6 * j + 1);

                auto gradxw_ky_r = rwDensityGrid.alphaDensityGradientX(6 * j + 2);

                auto gradxw_ky_i = rwDensityGrid.alphaDensityGradientX(6 * j + 3);

                auto gradxw_kz_r = rwDensityGrid.alphaDensityGradientX(6 * j + 4);

                auto gradxw_kz_i = rwDensityGrid.alphaDensityGradientX(6 * j + 5);

                // Grady terms

                auto gradyw_kx_r = rwDensityGrid.alphaDensityGradientY(6 * j);

                auto gradyw_kx_i = rwDensityGrid.alphaDensityGradientY(6 * j + 1);

                auto gradyw_ky_r = rwDensityGrid.alphaDensityGradientY(6 * j + 2);

                auto gradyw_ky_i = rwDensityGrid.alphaDensityGradientY(6 * j + 3);

                auto gradyw_kz_r = rwDensityGrid.alphaDensityGradientY(6 * j + 4);

                auto gradyw_kz_i = rwDensityGrid.alphaDensityGradientY(6 * j + 5);

                // Gradz terms

                auto gradzw_kx_r = rwDensityGrid.alphaDensityGradientZ(6 * j);

                auto gradzw_kx_i = rwDensityGrid.alphaDensityGradientZ(6 * j + 1);

                auto gradzw_ky_r = rwDensityGrid.alphaDensityGradientZ(6 * j + 2);

                auto gradzw_ky_i = rwDensityGrid.alphaDensityGradientZ(6 * j + 3);

                auto gradzw_kz_r = rwDensityGrid.alphaDensityGradientZ(6 * j + 4);

                auto gradzw_kz_i = rwDensityGrid.alphaDensityGradientZ(6 * j + 5);

                for (int32_t i = 0; i < npoints; i++)
                {
                    // RhoRho real

                    double jj_r = 2.0 * (rhow_kx_r[i] * rhow_kx_r[i] + rhow_ky_r[i] * rhow_ky_r[i] + rhow_kz_r[i] * rhow_kz_r[i])

                                  - 2.0 * (rhow_kx_i[i] * rhow_kx_i[i] + rhow_ky_i[i] * rhow_ky_i[i] + rhow_kz_i[i] * rhow_kz_i[i]);

                    rho_sig_x_r[i] = 4.0 * (rhow_kx_r[i] * rhow_kx_r[i] - rhow_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    rho_sig_y_r[i] = 4.0 * (rhow_ky_r[i] * rhow_ky_r[i] - rhow_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    rho_sig_z_r[i] = 4.0 * (rhow_kz_r[i] * rhow_kz_r[i] - rhow_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    rho_lam_xy_r[i] = (rhow_kx_r[i] * rhow_ky_r[i] - rhow_kx_i[i] * rhow_ky_i[i]

                                       + rhow_ky_r[i] * rhow_kx_r[i] - rhow_ky_i[i] * rhow_kx_i[i]);

                    rho_lam_xz_r[i] = (rhow_kx_r[i] * rhow_kz_r[i] - rhow_kx_i[i] * rhow_kz_i[i]

                                       + rhow_kz_r[i] * rhow_kx_r[i] - rhow_kz_i[i] * rhow_kx_i[i]);

                    rho_lam_yz_r[i] = (rhow_ky_r[i] * rhow_kz_r[i] - rhow_ky_i[i] * rhow_kz_i[i]

                                       + rhow_kz_r[i] * rhow_ky_r[i] - rhow_kz_i[i] * rhow_ky_i[i]);

                    // RhoRho imag

                    // rxw1rhow2 real

                    jj_r = 4.0 * (gradxw_kx_r[i] * rhow_kx_r[i] + gradxw_ky_r[i] * rhow_ky_r[i] + gradxw_kz_r[i] * rhow_kz_r[i])

                           - 4.0 * (gradxw_kx_i[i] * rhow_kx_i[i] + gradxw_ky_i[i] * rhow_ky_i[i] + gradxw_kz_i[i] * rhow_kz_i[i]);

                    rxw1rhow2_sig_x_r[i] = 8.0 * (gradxw_kx_r[i] * rhow_kx_r[i] - gradxw_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    rxw1rhow2_sig_y_r[i] = 8.0 * (gradxw_ky_r[i] * rhow_ky_r[i] - gradxw_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    rxw1rhow2_sig_z_r[i] = 8.0 * (gradxw_kz_r[i] * rhow_kz_r[i] - gradxw_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    rxw1rhow2_lam_xy_r[i] = 2.0 * (gradxw_kx_r[i] * rhow_ky_r[i] - gradxw_kx_i[i] * rhow_ky_i[i]

                                                   + gradxw_ky_r[i] * rhow_kx_r[i] - gradxw_ky_i[i] * rhow_kx_i[i]);

                    rxw1rhow2_lam_xz_r[i] = 2.0 * (gradxw_kx_r[i] * rhow_kz_r[i] - gradxw_kx_i[i] * rhow_kz_i[i]

                                                   + gradxw_kz_r[i] * rhow_kx_r[i] - gradxw_kz_i[i] * rhow_kx_i[i]);

                    rxw1rhow2_lam_yz_r[i] = 2.0 * (gradxw_ky_r[i] * rhow_kz_r[i] - gradxw_ky_i[i] * rhow_kz_i[i]

                                                   + gradxw_kz_r[i] * rhow_ky_r[i] - gradxw_kz_i[i] * rhow_ky_i[i]);

                    // rxw1rhow2 imag

                    // ryw1rhow2 real

                    jj_r = 4.0 * (gradyw_kx_r[i] * rhow_kx_r[i] + gradyw_ky_r[i] * rhow_ky_r[i] + gradyw_kz_r[i] * rhow_kz_r[i])

                           - 4.0 * (gradyw_kx_i[i] * rhow_kx_i[i] + gradyw_ky_i[i] * rhow_ky_i[i] + gradyw_kz_i[i] * rhow_kz_i[i]);

                    ryw1rhow2_sig_x_r[i] = 8.0 * (gradyw_kx_r[i] * rhow_kx_r[i] - gradyw_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    ryw1rhow2_sig_y_r[i] = 8.0 * (gradyw_ky_r[i] * rhow_ky_r[i] - gradyw_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    ryw1rhow2_sig_z_r[i] = 8.0 * (gradyw_kz_r[i] * rhow_kz_r[i] - gradyw_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    ryw1rhow2_lam_xy_r[i] = 2.0 * (gradyw_kx_r[i] * rhow_ky_r[i] - gradyw_kx_i[i] * rhow_ky_i[i]

                                                   + gradyw_ky_r[i] * rhow_kx_r[i] - gradyw_ky_i[i] * rhow_kx_i[i]);

                    ryw1rhow2_lam_xz_r[i] = 2.0 * (gradyw_kx_r[i] * rhow_kz_r[i] - gradyw_kx_i[i] * rhow_kz_i[i]

                                                   + gradyw_kz_r[i] * rhow_kx_r[i] - gradyw_kz_i[i] * rhow_kx_i[i]);

                    ryw1rhow2_lam_yz_r[i] = 2.0 * (gradyw_ky_r[i] * rhow_kz_r[i] - gradyw_ky_i[i] * rhow_kz_i[i]

                                                   + gradyw_kz_r[i] * rhow_ky_r[i] - gradyw_kz_i[i] * rhow_ky_i[i]);

                    // ryw1rhow2 imag

                    // rzw1rhow2 real

                    jj_r = 4.0 * (gradzw_kx_r[i] * rhow_kx_r[i] + gradzw_ky_r[i] * rhow_ky_r[i] + gradzw_kz_r[i] * rhow_kz_r[i])

                           - 4.0 * (gradzw_kx_i[i] * rhow_kx_i[i] + gradzw_ky_i[i] * rhow_ky_i[i] + gradzw_kz_i[i] * rhow_kz_i[i]);

                    rzw1rhow2_sig_x_r[i] = 8.0 * (gradzw_kx_r[i] * rhow_kx_r[i] - gradzw_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    rzw1rhow2_sig_y_r[i] = 8.0 * (gradzw_ky_r[i] * rhow_ky_r[i] - gradzw_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    rzw1rhow2_sig_z_r[i] = 8.0 * (gradzw_kz_r[i] * rhow_kz_r[i] - gradzw_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    rzw1rhow2_lam_xy_r[i] = 2.0 * (gradzw_kx_r[i] * rhow_ky_r[i] - gradzw_kx_i[i] * rhow_ky_i[i]

                                                   + gradzw_ky_r[i] * rhow_kx_r[i] - gradzw_ky_i[i] * rhow_kx_i[i]);

                    rzw1rhow2_lam_xz_r[i] = 2.0 * (gradzw_kx_r[i] * rhow_kz_r[i] - gradzw_kx_i[i] * rhow_kz_i[i]

                                                   + gradzw_kz_r[i] * rhow_kx_r[i] - gradzw_kz_i[i] * rhow_kx_i[i]);

                    rzw1rhow2_lam_yz_r[i] = 2.0 * (gradzw_ky_r[i] * rhow_kz_r[i] - gradzw_ky_i[i] * rhow_kz_i[i]

                                                   + gradzw_kz_r[i] * rhow_ky_r[i] - gradzw_kz_i[i] * rhow_ky_i[i]);

                    // rzw1rhow2 imag

                    // rxw1rxw2

                    jj_r = 2.0 * (gradxw_kx_r[i] * gradxw_kx_r[i] + gradxw_ky_r[i] * gradxw_ky_r[i] + gradxw_kz_r[i] * gradxw_kz_r[i])

                           - 2.0 * (gradxw_kx_i[i] * gradxw_kx_i[i] + gradxw_ky_i[i] * gradxw_ky_i[i] + gradxw_kz_i[i] * gradxw_kz_i[i]);

                    rxw1rxw2_sig_x_r[i] = 4.0 * (gradxw_kx_r[i] * gradxw_kx_r[i] - gradxw_kx_i[i] * gradxw_kx_i[i]) + jj_r;

                    rxw1rxw2_sig_y_r[i] = 4.0 * (gradxw_ky_r[i] * gradxw_ky_r[i] - gradxw_ky_i[i] * gradxw_ky_i[i]) + jj_r;

                    rxw1rxw2_sig_z_r[i] = 4.0 * (gradxw_kz_r[i] * gradxw_kz_r[i] - gradxw_kz_i[i] * gradxw_kz_i[i]) + jj_r;

                    rxw1rxw2_lam_xy_r[i] = gradxw_kx_r[i] * gradxw_ky_r[i] - gradxw_kx_i[i] * gradxw_ky_i[i]

                                           + gradxw_ky_r[i] * gradxw_kx_r[i] - gradxw_ky_i[i] * gradxw_kx_i[i];

                    rxw1rxw2_lam_xz_r[i] = gradxw_kx_r[i] * gradxw_kz_r[i] - gradxw_kx_i[i] * gradxw_kz_i[i]

                                           + gradxw_kz_r[i] * gradxw_kx_r[i] - gradxw_kz_i[i] * gradxw_kx_i[i];

                    rxw1rxw2_lam_yz_r[i] = gradxw_ky_r[i] * gradxw_kz_r[i] - gradxw_ky_i[i] * gradxw_kz_i[i]

                                           + gradxw_kz_r[i] * gradxw_ky_r[i] - gradxw_kz_i[i] * gradxw_ky_i[i];

                    // rxw1rxw2 imag

                    // rxw1ryw2

                    jj_r = 2.0 * (gradxw_kx_r[i] * gradyw_kx_r[i] + gradxw_ky_r[i] * gradyw_ky_r[i] + gradxw_kz_r[i] * gradyw_kz_r[i])

                           - 2.0 * (gradxw_kx_i[i] * gradyw_kx_i[i] + gradxw_ky_i[i] * gradyw_ky_i[i] + gradxw_kz_i[i] * gradyw_kz_i[i]);

                    rxw1ryw2_sig_x_r[i] = 4.0 * (gradxw_kx_r[i] * gradyw_kx_r[i] - gradxw_kx_i[i] * gradyw_kx_i[i]) + jj_r;

                    rxw1ryw2_sig_y_r[i] = 4.0 * (gradxw_ky_r[i] * gradyw_ky_r[i] - gradxw_ky_i[i] * gradyw_ky_i[i]) + jj_r;

                    rxw1ryw2_sig_z_r[i] = 4.0 * (gradxw_kz_r[i] * gradyw_kz_r[i] - gradxw_kz_i[i] * gradyw_kz_i[i]) + jj_r;

                    rxw1ryw2_lam_xy_r[i] = gradxw_kx_r[i] * gradyw_ky_r[i] - gradxw_kx_i[i] * gradyw_ky_i[i]

                                           + gradxw_ky_r[i] * gradyw_kx_r[i] - gradxw_ky_i[i] * gradyw_kx_i[i];

                    rxw1ryw2_lam_xz_r[i] = gradxw_kx_r[i] * gradyw_kz_r[i] - gradxw_kx_i[i] * gradyw_kz_i[i]

                                           + gradxw_kz_r[i] * gradyw_kx_r[i] - gradxw_kz_i[i] * gradyw_kx_i[i];

                    rxw1ryw2_lam_yz_r[i] = gradxw_ky_r[i] * gradyw_kz_r[i] - gradxw_ky_i[i] * gradyw_kz_i[i]

                                           + gradxw_kz_r[i] * gradyw_ky_r[i] - gradxw_kz_i[i] * gradyw_ky_i[i];

                    // rxw1ryw2 imag

                    // rxw1rzw2

                    jj_r = 2.0 * (gradxw_kx_r[i] * gradzw_kx_r[i] + gradxw_ky_r[i] * gradzw_ky_r[i] + gradxw_kz_r[i] * gradzw_kz_r[i])

                           - 2.0 * (gradxw_kx_i[i] * gradzw_kx_i[i] + gradxw_ky_i[i] * gradzw_ky_i[i] + gradxw_kz_i[i] * gradzw_kz_i[i]);

                    rxw1rzw2_sig_x_r[i] = 4.0 * (gradxw_kx_r[i] * gradzw_kx_r[i] - gradxw_kx_i[i] * gradzw_kx_i[i]) + jj_r;

                    rxw1rzw2_sig_y_r[i] = 4.0 * (gradxw_ky_r[i] * gradzw_ky_r[i] - gradxw_ky_i[i] * gradzw_ky_i[i]) + jj_r;

                    rxw1rzw2_sig_z_r[i] = 4.0 * (gradxw_kz_r[i] * gradzw_kz_r[i] - gradxw_kz_i[i] * gradzw_kz_i[i]) + jj_r;

                    rxw1rzw2_lam_xy_r[i] = gradxw_kx_r[i] * gradzw_ky_r[i] - gradxw_kx_i[i] * gradzw_ky_i[i]

                                           + gradxw_ky_r[i] * gradzw_kx_r[i] - gradxw_ky_i[i] * gradzw_kx_i[i];

                    rxw1rzw2_lam_xz_r[i] = gradxw_kx_r[i] * gradzw_kz_r[i] - gradxw_kx_i[i] * gradzw_kz_i[i]

                                           + gradxw_kz_r[i] * gradzw_kx_r[i] - gradxw_kz_i[i] * gradzw_kx_i[i];

                    rxw1rzw2_lam_yz_r[i] = gradxw_ky_r[i] * gradzw_kz_r[i] - gradxw_ky_i[i] * gradzw_kz_i[i]

                                           + gradxw_kz_r[i] * gradzw_ky_r[i] - gradxw_kz_i[i] * gradzw_ky_i[i];

                    // rxw1rzw2 imag

                    // ryw1rxw2

                    ryw1rxw2_sig_x_r[i] = rxw1ryw2_sig_x_r[i];

                    ryw1rxw2_sig_y_r[i] = rxw1ryw2_sig_y_r[i];

                    ryw1rxw2_sig_z_r[i] = rxw1ryw2_sig_z_r[i];

                    ryw1rxw2_lam_xy_r[i] = rxw1ryw2_lam_xy_r[i];

                    ryw1rxw2_lam_xz_r[i] = rxw1ryw2_lam_xz_r[i];

                    ryw1rxw2_lam_yz_r[i] = rxw1ryw2_lam_yz_r[i];

                    // ryw1rxw2 imag

                    // ryw1ryw2

                    jj_r = 2.0 * (gradyw_kx_r[i] * gradyw_kx_r[i] + gradyw_ky_r[i] * gradyw_ky_r[i] + gradyw_kz_r[i] * gradyw_kz_r[i])

                           - 2.0 * (gradyw_kx_i[i] * gradyw_kx_i[i] + gradyw_ky_i[i] * gradyw_ky_i[i] + gradyw_kz_i[i] * gradyw_kz_i[i]);

                    ryw1ryw2_sig_x_r[i] = 4.0 * (gradyw_kx_r[i] * gradyw_kx_r[i] - gradyw_kx_i[i] * gradyw_kx_i[i]) + jj_r;

                    ryw1ryw2_sig_y_r[i] = 4.0 * (gradyw_ky_r[i] * gradyw_ky_r[i] - gradyw_ky_i[i] * gradyw_ky_i[i]) + jj_r;

                    ryw1ryw2_sig_z_r[i] = 4.0 * (gradyw_kz_r[i] * gradyw_kz_r[i] - gradyw_kz_i[i] * gradyw_kz_i[i]) + jj_r;

                    ryw1ryw2_lam_xy_r[i] = gradyw_kx_r[i] * gradyw_ky_r[i] - gradyw_kx_i[i] * gradyw_ky_i[i]

                                           + gradyw_ky_r[i] * gradyw_kx_r[i] - gradyw_ky_i[i] * gradyw_kx_i[i];

                    ryw1ryw2_lam_xz_r[i] = gradyw_kx_r[i] * gradyw_kz_r[i] - gradyw_kx_i[i] * gradyw_kz_i[i]

                                           + gradyw_kz_r[i] * gradyw_kx_r[i] - gradyw_kz_i[i] * gradyw_kx_i[i];

                    ryw1ryw2_lam_yz_r[i] = gradyw_ky_r[i] * gradyw_kz_r[i] - gradyw_ky_i[i] * gradyw_kz_i[i]

                                           + gradyw_kz_r[i] * gradyw_ky_r[i] - gradyw_kz_i[i] * gradyw_ky_i[i];

                    // ryw1ryw2 imag

                    // ryw1rzw2

                    jj_r = 2.0 * (gradyw_kx_r[i] * gradzw_kx_r[i] + gradyw_ky_r[i] * gradzw_ky_r[i] + gradyw_kz_r[i] * gradzw_kz_r[i])

                           - 2.0 * (gradyw_kx_i[i] * gradzw_kx_i[i] + gradyw_ky_i[i] * gradzw_ky_i[i] + gradyw_kz_i[i] * gradzw_kz_i[i]);

                    ryw1rzw2_sig_x_r[i] = 4.0 * (gradyw_kx_r[i] * gradzw_kx_r[i] - gradyw_kx_i[i] * gradzw_kx_i[i]) + jj_r;

                    ryw1rzw2_sig_y_r[i] = 4.0 * (gradyw_ky_r[i] * gradzw_ky_r[i] - gradyw_ky_i[i] * gradzw_ky_i[i]) + jj_r;

                    ryw1rzw2_sig_z_r[i] = 4.0 * (gradyw_kz_r[i] * gradzw_kz_r[i] - gradyw_kz_i[i] * gradzw_kz_i[i]) + jj_r;

                    ryw1rzw2_lam_xy_r[i] = gradyw_kx_r[i] * gradzw_ky_r[i] - gradyw_kx_i[i] * gradzw_ky_i[i]

                                           + gradyw_ky_r[i] * gradzw_kx_r[i] - gradyw_ky_i[i] * gradzw_kx_i[i];

                    ryw1rzw2_lam_xz_r[i] = gradyw_kx_r[i] * gradzw_kz_r[i] - gradyw_kx_i[i] * gradzw_kz_i[i]

                                           + gradyw_kz_r[i] * gradzw_kx_r[i] - gradyw_kz_i[i] * gradzw_kx_i[i];

                    ryw1rzw2_lam_yz_r[i] = gradyw_ky_r[i] * gradzw_kz_r[i] - gradyw_ky_i[i] * gradzw_kz_i[i]

                                           + gradyw_kz_r[i] * gradzw_ky_r[i] - gradyw_kz_i[i] * gradzw_ky_i[i];

                    // ryw1rzw2 imag

                    // rzw1rxw2

                    rzw1rxw2_sig_x_r[i] = rxw1rzw2_sig_x_r[i];

                    rzw1rxw2_sig_y_r[i] = rxw1rzw2_sig_y_r[i];

                    rzw1rxw2_sig_z_r[i] = rxw1rzw2_sig_z_r[i];

                    rzw1rxw2_lam_xy_r[i] = rxw1rzw2_lam_xy_r[i];

                    rzw1rxw2_lam_xz_r[i] = rxw1rzw2_lam_xz_r[i];

                    rzw1rxw2_lam_yz_r[i] = rxw1rzw2_lam_yz_r[i];

                    // rzw1rxw2 imag

                    // rzw1ryw2

                    rzw1ryw2_sig_x_r[i] = ryw1rzw2_sig_x_r[i];

                    rzw1ryw2_sig_y_r[i] = ryw1rzw2_sig_y_r[i];

                    rzw1ryw2_sig_z_r[i] = ryw1rzw2_sig_z_r[i];

                    rzw1ryw2_lam_xy_r[i] = ryw1rzw2_lam_xy_r[i];

                    rzw1ryw2_lam_xz_r[i] = ryw1rzw2_lam_xz_r[i];

                    rzw1ryw2_lam_yz_r[i] = ryw1rzw2_lam_yz_r[i];

                    // rzw1ryw2 imag

                    // rzw1rzw2

                    jj_r = 2.0 * (gradzw_kx_r[i] * gradzw_kx_r[i] + gradzw_ky_r[i] * gradzw_ky_r[i] + gradzw_kz_r[i] * gradzw_kz_r[i])

                           - 2.0 * (gradzw_kx_i[i] * gradzw_kx_i[i] + gradzw_ky_i[i] * gradzw_ky_i[i] + gradzw_kz_i[i] * gradzw_kz_i[i]);

                    rzw1rzw2_sig_x_r[i] = 4.0 * (gradzw_kx_r[i] * gradzw_kx_r[i] - gradzw_kx_i[i] * gradzw_kx_i[i]) + jj_r;

                    rzw1rzw2_sig_y_r[i] = 4.0 * (gradzw_ky_r[i] * gradzw_ky_r[i] - gradzw_ky_i[i] * gradzw_ky_i[i]) + jj_r;

                    rzw1rzw2_sig_z_r[i] = 4.0 * (gradzw_kz_r[i] * gradzw_kz_r[i] - gradzw_kz_i[i] * gradzw_kz_i[i]) + jj_r;

                    rzw1rzw2_lam_xy_r[i] = gradzw_kx_r[i] * gradzw_ky_r[i] - gradzw_kx_i[i] * gradzw_ky_i[i]

                                           + gradzw_ky_r[i] * gradzw_kx_r[i] - gradzw_ky_i[i] * gradzw_kx_i[i];

                    rzw1rzw2_lam_xz_r[i] = gradzw_kx_r[i] * gradzw_kz_r[i] - gradzw_kx_i[i] * gradzw_kz_i[i]

                                           + gradzw_kz_r[i] * gradzw_kx_r[i] - gradzw_kz_i[i] * gradzw_kx_i[i];

                    rzw1rzw2_lam_yz_r[i] = gradzw_ky_r[i] * gradzw_kz_r[i] - gradzw_ky_i[i] * gradzw_kz_i[i]

                                           + gradzw_kz_r[i] * gradzw_ky_r[i] - gradzw_kz_i[i] * gradzw_ky_i[i];

                    // rzw1rzw2 imag
                }
            }
        }
        if (fstr::upcase(quadMode) == "QRF")
        {
            for (int32_t j = 0; j < numdens / 2; j++)
            {
                // Perturbed densities

                auto rhow1a_r = rwDensityGrid.alphaDensity(4 * j);

                auto gradw1a_x_r = rwDensityGrid.alphaDensityGradientX(4 * j);

                auto gradw1a_y_r = rwDensityGrid.alphaDensityGradientY(4 * j);

                auto gradw1a_z_r = rwDensityGrid.alphaDensityGradientZ(4 * j);

                auto rhow1a_i = rwDensityGrid.alphaDensity(4 * j + 1);

                auto gradw1a_x_i = rwDensityGrid.alphaDensityGradientX(4 * j + 1);

                auto gradw1a_y_i = rwDensityGrid.alphaDensityGradientY(4 * j + 1);

                auto gradw1a_z_i = rwDensityGrid.alphaDensityGradientZ(4 * j + 1);

                auto rhow2a_r = rwDensityGrid.alphaDensity(4 * j + 2);

                auto gradw2a_x_r = rwDensityGrid.alphaDensityGradientX(4 * j + 2);

                auto gradw2a_y_r = rwDensityGrid.alphaDensityGradientY(4 * j + 2);

                auto gradw2a_z_r = rwDensityGrid.alphaDensityGradientZ(4 * j + 2);

                auto rhow2a_i = rwDensityGrid.alphaDensity(4 * j + 3);

                auto gradw2a_x_i = rwDensityGrid.alphaDensityGradientX(4 * j + 3);

                auto gradw2a_y_i = rwDensityGrid.alphaDensityGradientY(4 * j + 3);

                auto gradw2a_z_i = rwDensityGrid.alphaDensityGradientZ(4 * j + 3);

                // Perturbed density products to be stored

                auto rhow1rhow2_r = densityGridAB.rhow1rhow2(2 * j);

                auto rhow1rhow2_i = densityGridAB.rhow1rhow2(2 * j + 1);

                // rxw1rhow2

                auto rxw1rhow2_r = densityGridAB.rxw1rhow2(2 * j);

                auto rxw1rhow2_i = densityGridAB.rxw1rhow2(2 * j + 1);

                // ryw1rhow2

                auto ryw1rhow2_r = densityGridAB.ryw1rhow2(2 * j);

                auto ryw1rhow2_i = densityGridAB.ryw1rhow2(2 * j + 1);

                // rzw1rhow2

                auto rzw1rhow2_r = densityGridAB.rzw1rhow2(2 * j);

                auto rzw1rhow2_i = densityGridAB.rzw1rhow2(2 * j + 1);

                // rAw1rBw2

                auto rxw1rxw2_r = densityGridAB.rxw1rxw2(2 * j);

                auto rxw1rxw2_i = densityGridAB.rxw1rxw2(2 * j + 1);

                auto rxw1ryw2_r = densityGridAB.rxw1ryw2(2 * j);

                auto rxw1ryw2_i = densityGridAB.rxw1ryw2(2 * j + 1);

                auto rxw1rzw2_r = densityGridAB.rxw1rzw2(2 * j);

                auto rxw1rzw2_i = densityGridAB.rxw1rzw2(2 * j + 1);

                auto ryw1rxw2_r = densityGridAB.ryw1rxw2(2 * j);

                auto ryw1rxw2_i = densityGridAB.ryw1rxw2(2 * j + 1);

                auto ryw1ryw2_r = densityGridAB.ryw1ryw2(2 * j);

                auto ryw1ryw2_i = densityGridAB.ryw1ryw2(2 * j + 1);

                auto ryw1rzw2_r = densityGridAB.ryw1rzw2(2 * j);

                auto ryw1rzw2_i = densityGridAB.ryw1rzw2(2 * j + 1);

                auto rzw1rxw2_r = densityGridAB.rzw1rxw2(2 * j);

                auto rzw1rxw2_i = densityGridAB.rzw1rxw2(2 * j + 1);

                auto rzw1ryw2_r = densityGridAB.rzw1ryw2(2 * j);

                auto rzw1ryw2_i = densityGridAB.rzw1ryw2(2 * j + 1);

                auto rzw1rzw2_r = densityGridAB.rzw1rzw2(2 * j);

                auto rzw1rzw2_i = densityGridAB.rzw1rzw2(2 * j + 1);

                for (int32_t i = 0; i < npoints; i++)
                {
                    // RW1 densities

                    double rxw1_r = gradw1a_x_r[i];

                    double ryw1_r = gradw1a_y_r[i];

                    double rzw1_r = gradw1a_z_r[i];

                    double rxw1_i = gradw1a_x_i[i];

                    double ryw1_i = gradw1a_y_i[i];

                    double rzw1_i = gradw1a_z_i[i];

                    // RW2 densities

                    double rxw2_r = gradw2a_x_r[i];

                    double ryw2_r = gradw2a_y_r[i];

                    double rzw2_r = gradw2a_z_r[i];

                    double rxw2_i = gradw2a_x_i[i];

                    double ryw2_i = gradw2a_y_i[i];

                    double rzw2_i = gradw2a_z_i[i];

                    // Densities for terms 1-3

                    rhow1rhow2_r[i] = 2.0 * (rhow1a_r[i] * rhow2a_r[i] - rhow1a_i[i] * rhow2a_i[i]);

                    rhow1rhow2_i[i] = 2.0 * (rhow1a_r[i] * rhow2a_i[i] + rhow1a_i[i] * rhow2a_r[i]);

                    rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2a_r[i] - rxw1_i * rhow2a_i[i]

                                            + rxw2_r * rhow1a_r[i] - rxw2_i * rhow1a_i[i]);

                    rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2a_i[i] + rxw1_i * rhow2a_r[i]

                                            + rxw2_r * rhow1a_i[i] + rxw2_i * rhow1a_r[i]);

                    ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2a_r[i] - ryw1_i * rhow2a_i[i]

                                            + ryw2_r * rhow1a_r[i] - ryw2_i * rhow1a_i[i]);

                    ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2a_i[i] + ryw1_i * rhow2a_r[i]

                                            + ryw2_r * rhow1a_i[i] + ryw2_i * rhow1a_r[i]);

                    rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2a_r[i] - rzw1_i * rhow2a_i[i]

                                            + rzw2_r * rhow1a_r[i] - rzw2_i * rhow1a_i[i]);

                    rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2a_i[i] + rzw1_i * rhow2a_r[i]

                                            + rzw2_r * rhow1a_i[i] + rzw2_i * rhow1a_r[i]);

                    rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;

                    rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;

                    rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;

                    rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;

                    rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;
                }
            }
        }
    }
}

std::ostream&
operator<<(std::ostream& output, const CDensityGridQuad& source)
{
    output << std::endl;

    output << "[CDensityGridQuad (Object):" << &source << "]" << std::endl;

    output << "_gridType: " << to_string(source._gridType) << std::endl;

    output << "_densityValues: " << std::endl;

    output << source._densityValues << std::endl;

    return output;
}
