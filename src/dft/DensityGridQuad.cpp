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

    if (xcFuncType == xcfun::mgga) ncomp = (_gridType == dengrid::ab) ? 24 : 6;

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

double
CDensityGridQuad::prod2_r(double B_r, double B_i, double C_r, double C_i)
{
    double BC  =  (B_r * C_r- B_i * C_i);

    return BC;
}

double
CDensityGridQuad::prod2_i(double B_r, double B_i, double C_r, double C_i)
{
    double BC  =  (B_i * C_r + B_r * C_i);

    return BC;
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
CDensityGridQuad::gam(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::lima) return nullptr;

    return _densityValues.data(iDensityMatrix);
}

double*
CDensityGridQuad::gam(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::lima) return nullptr;

    return _densityValues.data(iDensityMatrix);
}

const double*
CDensityGridQuad::gamX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(1 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::gamX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(1 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::gamY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::gamY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::gamZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::gamZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::gamXX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::gamXX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::gamXY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::gamXY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::gamXZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::gamXZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::gamYX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::gamYX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::gamYY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::gamYY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::gamYZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::gamYZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::gamZX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::gamZX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::gamZY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(11 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::gamZY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(11 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::gamZZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(12 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::gamZZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(12 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridQuad::rt_gam(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(13 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}
    

double*
CDensityGridQuad::rt_gam(const int32_t iDensityMatrix) 
{
    if (_gridType == dengrid::ab) return _densityValues.data(13 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridQuad::rl_gam(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(14 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


double*
CDensityGridQuad::rl_gam(const int32_t iDensityMatrix) 
{
    if (_gridType == dengrid::ab) return _densityValues.data(14 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridQuad::tt_gam(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(15 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


double*
CDensityGridQuad::tt_gam(const int32_t iDensityMatrix) 
{
    if (_gridType == dengrid::ab) return _densityValues.data(15 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridQuad::tl_gam(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(16 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


double*
CDensityGridQuad::tl_gam(const int32_t iDensityMatrix) 
{
    if (_gridType == dengrid::ab) return _densityValues.data(16 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridQuad::ll_gam(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(17 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


double*
CDensityGridQuad::ll_gam(const int32_t iDensityMatrix) 
{
    if (_gridType == dengrid::ab) return _densityValues.data(17 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridQuad::st_gamX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(18 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


double*
CDensityGridQuad::st_gamX(const int32_t iDensityMatrix) 
{
    if (_gridType == dengrid::ab) return _densityValues.data(18 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridQuad::st_gamY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(19 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


double*
CDensityGridQuad::st_gamY(const int32_t iDensityMatrix) 
{
    if (_gridType == dengrid::ab) return _densityValues.data(19 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridQuad::st_gamZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(20 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


double*
CDensityGridQuad::st_gamZ(const int32_t iDensityMatrix) 
{
    if (_gridType == dengrid::ab) return _densityValues.data(20 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridQuad::sl_gamX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(21 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::sl_gamX(const int32_t iDensityMatrix) 
{
    if (_gridType == dengrid::ab) return _densityValues.data(21 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}


const double*
CDensityGridQuad::sl_gamY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(22 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::sl_gamY(const int32_t iDensityMatrix) 
{
    if (_gridType == dengrid::ab) return _densityValues.data(22 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

const double*
CDensityGridQuad::sl_gamZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(23 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

double*
CDensityGridQuad::sl_gamZ(const int32_t iDensityMatrix) 
{
    if (_gridType == dengrid::ab) return _densityValues.data(23 * _nDensityMatrices + iDensityMatrix);

    return nullptr;
}

void
CDensityGridQuad::DensityProd(const CDensityGrid& rwDensityGrid,
                              const xcfun         xcFuncType,
                              const int32_t             numdens,
                              const std::string&  quadMode)

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

                auto rho_sig_x_r = gam(12 * j);

                auto rho_sig_x_i = gam(12 * j + 1);

                auto rho_sig_y_r = gam(12 * j + 2);

                auto rho_sig_y_i = gam(12 * j + 3);

                auto rho_sig_z_r = gam(12 * j + 4);

                auto rho_sig_z_i = gam(12 * j + 5);

                auto rho_lam_xy_r = gam(12 * j + 6);

                auto rho_lam_xy_i = gam(12 * j + 7);

                auto rho_lam_xz_r = gam(12 * j + 8);

                auto rho_lam_xz_i = gam(12 * j + 9);

                auto rho_lam_yz_r = gam(12 * j + 10);

                auto rho_lam_yz_i = gam(12 * j + 11);

                auto rhow_kx_r = rwDensityGrid.alphaDensity(6 * j);

                auto rhow_kx_i = rwDensityGrid.alphaDensity(6 * j + 1);

                auto rhow_ky_r = rwDensityGrid.alphaDensity(6 * j + 2);

                auto rhow_ky_i = rwDensityGrid.alphaDensity(6 * j + 3);

                auto rhow_kz_r = rwDensityGrid.alphaDensity(6 * j + 4);

                auto rhow_kz_i = rwDensityGrid.alphaDensity(6 * j + 5);

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

                auto rho_sig_x_r = gam(6 * j);

                auto rho_sig_y_r = gam(6 * j + 1);

                auto rho_sig_z_r = gam(6 * j + 2);

                auto rho_lam_xy_r = gam(6 * j + 3);

                auto rho_lam_xz_r = gam(6 * j + 4);

                auto rho_lam_yz_r = gam(6 * j + 5);

                auto rhow_kx_r = rwDensityGrid.alphaDensity(6 * j);

                auto rhow_kx_i = rwDensityGrid.alphaDensity(6 * j + 1);

                auto rhow_ky_r = rwDensityGrid.alphaDensity(6 * j + 2);

                auto rhow_ky_i = rwDensityGrid.alphaDensity(6 * j + 3);

                auto rhow_kz_r = rwDensityGrid.alphaDensity(6 * j + 4);

                auto rhow_kz_i = rwDensityGrid.alphaDensity(6 * j + 5);

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
        if (fstr::upcase(quadMode) == "TPA_II") 
        {
            // This code is inteded to compute F_b_sigma fock matrices for the final E3 contraction for tpa calculations.
            for (int32_t j = 0; j < numdens / 6; j++)
            {
                
                auto rr_x_r = gam(6 * j);

                auto rr_x_i = gam(6 * j + 1);

                auto rr_y_r = gam(6 * j + 2);

                auto rr_y_i = gam(6 * j + 3);

                auto rr_z_r = gam(6 * j + 4);

                auto rr_z_i = gam(6 * j + 5);


                auto rhoBx_r = rwDensityGrid.alphaDensity(36 * j);

                auto rhoBx_i = rwDensityGrid.alphaDensity(36 * j + 1);

                auto rhoBy_r = rwDensityGrid.alphaDensity(36 * j + 2);

                auto rhoBy_i = rwDensityGrid.alphaDensity(36 * j + 3);

                auto rhoBz_r = rwDensityGrid.alphaDensity(36 * j + 4);

                auto rhoBz_i = rwDensityGrid.alphaDensity(36 * j + 5);  


                auto rhoCx_r = rwDensityGrid.alphaDensity(36 * j + 6);

                auto rhoCx_i = rwDensityGrid.alphaDensity(36 * j + 7);

                auto rhoCy_r = rwDensityGrid.alphaDensity(36 * j + 8);

                auto rhoCy_i = rwDensityGrid.alphaDensity(36 * j + 9);

                auto rhoCz_r = rwDensityGrid.alphaDensity(36 * j + 10);

                auto rhoCz_i = rwDensityGrid.alphaDensity(36 * j + 11);  


                auto Rho_sig_xx_r = rwDensityGrid.alphaDensity(36 * j + 12);

                auto Rho_sig_xx_i = rwDensityGrid.alphaDensity(36 * j + 13);

                auto Rho_sig_yy_r = rwDensityGrid.alphaDensity(36 * j + 14);

                auto Rho_sig_yy_i = rwDensityGrid.alphaDensity(36 * j + 15);

                auto Rho_sig_zz_r = rwDensityGrid.alphaDensity(36 * j + 16);

                auto Rho_sig_zz_i = rwDensityGrid.alphaDensity(36 * j + 17);


                auto Rho_sig_xy_r = rwDensityGrid.alphaDensity(36 * j + 18);

                auto Rho_sig_xy_i = rwDensityGrid.alphaDensity(36 * j + 19);

                auto Rho_sig_xz_r = rwDensityGrid.alphaDensity(36 * j + 20);

                auto Rho_sig_xz_i = rwDensityGrid.alphaDensity(36 * j + 21);

                auto Rho_sig_yz_r = rwDensityGrid.alphaDensity(36 * j + 22);

                auto Rho_sig_yz_i = rwDensityGrid.alphaDensity(36 * j + 23);


                auto Rho_lamtau_xx_r = rwDensityGrid.alphaDensity(36 * j + 24);

                auto Rho_lamtau_xx_i = rwDensityGrid.alphaDensity(36 * j + 25);

                auto Rho_lamtau_yy_r = rwDensityGrid.alphaDensity(36 * j + 26);

                auto Rho_lamtau_yy_i = rwDensityGrid.alphaDensity(36 * j + 27);

                auto Rho_lamtau_zz_r = rwDensityGrid.alphaDensity(36 * j + 28);

                auto Rho_lamtau_zz_i = rwDensityGrid.alphaDensity(36 * j + 29);


                auto Rho_lamtau_xy_r = rwDensityGrid.alphaDensity(36 * j + 30);

                auto Rho_lamtau_xy_i = rwDensityGrid.alphaDensity(36 * j + 31);

                auto Rho_lamtau_xz_r = rwDensityGrid.alphaDensity(36 * j + 32);

                auto Rho_lamtau_xz_i = rwDensityGrid.alphaDensity(36 * j + 33);

                auto Rho_lamtau_yz_r = rwDensityGrid.alphaDensity(36 * j + 34);

                auto Rho_lamtau_yz_i = rwDensityGrid.alphaDensity(36 * j + 35);
                

                for (int32_t i = 0; i < npoints; i++)
                {

                    rr_x_r[i] = 2.0 * (Rho_sig_xx_r[i] * rhoCx_r[i] - Rho_sig_xx_i[i] * rhoCx_i[i] 
                    
                                    +  Rho_sig_xy_r[i] * rhoCy_r[i] - Rho_sig_xy_i[i] * rhoCy_i[i] 

                                    + Rho_sig_xz_r[i] * rhoCz_r[i] -  Rho_sig_xz_i[i] * rhoCz_i[i] 

                                    + Rho_lamtau_xx_r[i] * rhoBx_r[i] - Rho_lamtau_xx_i[i] * rhoBx_i[i] 

                                    + Rho_lamtau_xy_r[i] * rhoBy_r[i] - Rho_lamtau_xy_i[i] * rhoBy_i[i] 
                                    
                                    + Rho_lamtau_xz_r[i] * rhoBz_r[i] - Rho_lamtau_xz_i[i] * rhoBz_i[i]);


                    rr_y_r[i] = 2.0 * (Rho_sig_xy_r[i]* rhoCx_r[i] - Rho_sig_xy_i[i] * rhoCx_i[i] 
                    
                                     + Rho_sig_yy_r[i]* rhoCy_r[i] - Rho_sig_yy_i[i] * rhoCy_i[i] 
                                                                      
                                     + Rho_sig_yz_r[i]* rhoCz_r[i] - Rho_sig_yz_i[i] * rhoCz_i[i]
     
                                     + Rho_lamtau_xy_r[i]* rhoBx_r[i] - Rho_lamtau_xy_i[i] * rhoBx_i[i] 

                                     + Rho_lamtau_yy_r[i]* rhoBy_r[i] - Rho_lamtau_yy_i[i] * rhoBy_i[i] 
                                     
                                     + Rho_lamtau_yz_r[i]* rhoBz_r[i] - Rho_lamtau_yz_i[i] * rhoBz_i[i]);


                    rr_z_r[i] = 2.0 * (Rho_sig_xz_r[i] * rhoCx_r[i] - Rho_sig_xz_i[i] * rhoCx_i[i] 
                 
                                     + Rho_sig_yz_r[i] * rhoCy_r[i] - Rho_sig_yz_i[i] * rhoCy_i[i] 
                                    
                                    +  Rho_sig_zz_r[i] * rhoCz_r[i] - Rho_sig_zz_i[i] * rhoCz_i[i] 

                                    + Rho_lamtau_xz_r[i]* rhoBx_r[i] - Rho_lamtau_xz_i[i] * rhoBx_i[i]

                                    + Rho_lamtau_yz_r[i]* rhoBy_r[i] - Rho_lamtau_yz_i[i] * rhoBy_i[i] 
                                    
                                    + Rho_lamtau_zz_r[i]* rhoBz_r[i] - Rho_lamtau_zz_i[i] * rhoBz_i[i]);


                    rr_x_i[i] = 2.0 * (Rho_sig_xx_i[i]* rhoCx_r[i] + Rho_sig_xx_r[i] * rhoCx_i[i] 
                    
                                    + Rho_sig_xy_i[i]* rhoCy_r[i] + Rho_sig_xy_r[i] * rhoCy_i[i] 

                                    + Rho_sig_xz_i[i]* rhoCz_r[i] + Rho_sig_xz_r[i] * rhoCz_i[i] 

                                    + Rho_lamtau_xx_i[i]* rhoBx_r[i] + Rho_lamtau_xx_r[i] * rhoBx_i[i] 

                                    + Rho_lamtau_xy_i[i]* rhoBy_r[i] + Rho_lamtau_xy_r[i] * rhoBy_i[i] 

                                    + Rho_lamtau_xz_i[i]* rhoBz_r[i] + Rho_lamtau_xz_r[i] * rhoBz_i[i] );


                    rr_y_i[i] = 2.0 * (Rho_sig_xy_i[i]* rhoCx_r[i] + Rho_sig_xy_r[i] * rhoCx_i[i] 

                                    + Rho_sig_yy_i[i]* rhoCy_r[i] + Rho_sig_yy_r[i] * rhoCy_i[i] 

                                    + Rho_sig_yz_i[i]* rhoCz_r[i] + Rho_sig_yz_r[i] * rhoCz_i[i] 

                                    + Rho_lamtau_xy_i[i]* rhoBx_r[i] + Rho_lamtau_xy_r[i] * rhoBx_i[i] 

                                    + Rho_lamtau_yy_i[i]* rhoBy_r[i] + Rho_lamtau_yy_r[i] * rhoBy_i[i] 
                                
                                    + Rho_lamtau_yz_i[i]* rhoBz_r[i] + Rho_lamtau_yz_r[i] * rhoBz_i[i] );


                    rr_z_i[i] = 2.0 * (Rho_sig_xz_i[i]* rhoCx_r[i] + Rho_sig_xz_r[i] * rhoCx_i[i] 
                                
                                    + Rho_sig_yz_i[i]* rhoCy_r[i] + Rho_sig_yz_r[i] * rhoCy_i[i]

                                    + Rho_sig_zz_i[i]* rhoCz_r[i] + Rho_sig_zz_r[i] * rhoCz_i[i] 

                                    + Rho_lamtau_xz_i[i]* rhoBx_r[i] + Rho_lamtau_xz_r[i] * rhoBx_i[i] 
                                    
                                    + Rho_lamtau_yz_i[i]* rhoBy_r[i] + Rho_lamtau_yz_r[i] * rhoBy_i[i]
                                    
                                    + Rho_lamtau_zz_i[i]* rhoBz_r[i] + Rho_lamtau_zz_r[i] * rhoBz_i[i]);

                }
            }
        }
        if (fstr::upcase(quadMode) == "REDTPA_I")
        {

            // This routine computes the first-order two-times transformed Fock matrices for the E[4] contraction for TPA calculations

            for (int32_t j = 0; j < numdens / 6; j++)
            {
                
                auto gam_sig_xx_r = gam(6 * j);
                auto gam_sig_yy_r = gam(6 * j + 1);
                auto gam_sig_zz_r = gam(6 * j + 2);
                auto gam_sig_xy_r = gam(6 * j + 3);
                auto gam_sig_xz_r = gam(6 * j + 4);
                auto gam_sig_yz_r = gam(6 * j + 5);
                
                auto rhoBx_r = rwDensityGrid.alphaDensity(6 * j);
                auto rhoBx_i = rwDensityGrid.alphaDensity(6 * j + 1);
                auto rhoBy_r = rwDensityGrid.alphaDensity(6 * j + 2);
                auto rhoBy_i = rwDensityGrid.alphaDensity(6 * j + 3);
                auto rhoBz_r = rwDensityGrid.alphaDensity(6 * j + 4);
                auto rhoBz_i = rwDensityGrid.alphaDensity(6 * j + 5);  

                for (int32_t i = 0; i < npoints; i++)
                {

                    double sig_term_r =  + 6.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 6.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 6.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoBz_r[i],rhoBz_i[i]);

                    gam_sig_xx_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBx_r[i],rhoBx_i[i]) + sig_term_r;

                    gam_sig_yy_r[i] = 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoBy_r[i],rhoBy_i[i]) + sig_term_r;

                    gam_sig_zz_r[i] = 12.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoBz_r[i],rhoBz_i[i]) + sig_term_r;

                    gam_sig_xy_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBy_r[i],rhoBy_i[i]);

                    gam_sig_xz_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBz_r[i],rhoBz_i[i]);

                    gam_sig_yz_r[i] = 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoBz_r[i],rhoBz_i[i]);
                    

                }
            }
        }
        if (fstr::upcase(quadMode) == "REDTPA_II")  
        {
            // This code is inteded to compute F_b_cd fock matrices for the final E3 contraction for tpa calculations.

            for (int32_t j = 0; j < numdens / 6; j++)
            {
                
                auto rr_x_r = gam(6 * j);
                auto rr_x_i = gam(6 * j + 1);
                auto rr_y_r = gam(6 * j + 2);
                auto rr_y_i = gam(6 * j + 3);
                auto rr_z_r = gam(6 * j + 4);
                auto rr_z_i = gam(6 * j + 5);
                
                auto rhoCx_r = rwDensityGrid.alphaDensity(18 * j);
                auto rhoCx_i = rwDensityGrid.alphaDensity(18 * j + 1);
                auto rhoCy_r = rwDensityGrid.alphaDensity(18 * j + 2);
                auto rhoCy_i = rwDensityGrid.alphaDensity(18 * j + 3);
                auto rhoCz_r = rwDensityGrid.alphaDensity(18 * j + 4);
                auto rhoCz_i = rwDensityGrid.alphaDensity(18 * j + 5);  
                auto Rho_sig_xx_r = rwDensityGrid.alphaDensity(18 * j + 6);
                auto Rho_sig_xx_i = rwDensityGrid.alphaDensity(18 * j + 7);
                auto Rho_sig_yy_r = rwDensityGrid.alphaDensity(18 * j + 8);
                auto Rho_sig_yy_i = rwDensityGrid.alphaDensity(18 * j + 9);
                auto Rho_sig_zz_r = rwDensityGrid.alphaDensity(18 * j + 10);
                auto Rho_sig_zz_i = rwDensityGrid.alphaDensity(18 * j + 11);
                auto Rho_sig_xy_r = rwDensityGrid.alphaDensity(18 * j + 12);
                auto Rho_sig_xy_i = rwDensityGrid.alphaDensity(18 * j + 13);
                auto Rho_sig_xz_r = rwDensityGrid.alphaDensity(18 * j + 14);
                auto Rho_sig_xz_i = rwDensityGrid.alphaDensity(18 * j + 15);
                auto Rho_sig_yz_r = rwDensityGrid.alphaDensity(18 * j + 16);
                auto Rho_sig_yz_i = rwDensityGrid.alphaDensity(18 * j + 17);

                for (int32_t i = 0; i < npoints; i++)
                {

                    rr_x_r[i] = 2.0 * (Rho_sig_xx_r[i] * rhoCx_r[i] - Rho_sig_xx_i[i] * rhoCx_i[i] 
                    
                                     + Rho_sig_xy_r[i] * rhoCy_r[i] - Rho_sig_xy_i[i] * rhoCy_i[i] 

                                     + Rho_sig_xz_r[i] * rhoCz_r[i] -  Rho_sig_xz_i[i] * rhoCz_i[i]);


                    rr_y_r[i] = 2.0 * (Rho_sig_xy_r[i]* rhoCx_r[i] - Rho_sig_xy_i[i] * rhoCx_i[i] 
                    
                                     + Rho_sig_yy_r[i]* rhoCy_r[i] - Rho_sig_yy_i[i] * rhoCy_i[i] 
                                                                      
                                     + Rho_sig_yz_r[i]* rhoCz_r[i] - Rho_sig_yz_i[i] * rhoCz_i[i]);


                    rr_z_r[i] = 2.0 * (Rho_sig_xz_r[i] * rhoCx_r[i] - Rho_sig_xz_i[i] * rhoCx_i[i] 
                 
                                     + Rho_sig_yz_r[i] * rhoCy_r[i] - Rho_sig_yz_i[i] * rhoCy_i[i] 
                                    
                                    +  Rho_sig_zz_r[i] * rhoCz_r[i] - Rho_sig_zz_i[i] * rhoCz_i[i]);


                    rr_x_i[i] = 2.0 * (Rho_sig_xx_i[i]* rhoCx_r[i] + Rho_sig_xx_r[i] * rhoCx_i[i] 
                    
                                    + Rho_sig_xy_i[i]* rhoCy_r[i] + Rho_sig_xy_r[i] * rhoCy_i[i] 

                                    + Rho_sig_xz_i[i]* rhoCz_r[i] + Rho_sig_xz_r[i] * rhoCz_i[i] );


                    rr_y_i[i] = 2.0 * (Rho_sig_xy_i[i]* rhoCx_r[i] + Rho_sig_xy_r[i] * rhoCx_i[i] 

                                    + Rho_sig_yy_i[i]* rhoCy_r[i] + Rho_sig_yy_r[i] * rhoCy_i[i] 

                                    + Rho_sig_yz_i[i]* rhoCz_r[i] + Rho_sig_yz_r[i] * rhoCz_i[i] );


                    rr_z_i[i] = 2.0 * (Rho_sig_xz_i[i]* rhoCx_r[i] + Rho_sig_xz_r[i] * rhoCx_i[i] 
                                
                                    + Rho_sig_yz_i[i]* rhoCy_r[i] + Rho_sig_yz_r[i] * rhoCy_i[i]

                                    + Rho_sig_zz_i[i]* rhoCz_r[i] + Rho_sig_zz_r[i] * rhoCz_i[i]);
                }
            }
        }
        if (fstr::upcase(quadMode) == "CRF_II")
        {
            // This routine is for computing the second-order fock matrices for the E[3] contraction of the general cubic response function routine

            for (int32_t j = 0; j < numdens / 2; j++)
            {
                auto rr_123_r = gam(2 * j);

                auto rr_123_i = gam(2 * j + 1);

                auto rhowba_r = rwDensityGrid.alphaDensity(12 * j);
                
                auto rhowba_i = rwDensityGrid.alphaDensity(12 * j + 1);

                auto rhowca_r = rwDensityGrid.alphaDensity(12 * j + 2);

                auto rhowca_i = rwDensityGrid.alphaDensity(12 * j + 3);

                auto rhowda_r = rwDensityGrid.alphaDensity(12 * j + 4);

                auto rhowda_i = rwDensityGrid.alphaDensity(12 * j + 5);

                auto rhowbc_r = rwDensityGrid.alphaDensity(12 * j + 6);

                auto rhowbc_i = rwDensityGrid.alphaDensity(12 * j + 7);

                auto rhowbd_r = rwDensityGrid.alphaDensity(12 * j + 8);

                auto rhowbd_i = rwDensityGrid.alphaDensity(12 * j + 9);

                auto rhowcd_r = rwDensityGrid.alphaDensity(12 * j + 10);

                auto rhowcd_i = rwDensityGrid.alphaDensity(12 * j + 11);

                
                for (int32_t i = 0; i < npoints; i++)
                {
                    rr_123_r[i] = 2.0 * (rhowda_r[i]* rhowbc_r[i] - rhowda_i[i] * rhowbc_i[i])

                                + 2.0 * (rhowca_r[i]* rhowbd_r[i] - rhowca_i[i] * rhowbd_i[i])
                                
                                + 2.0 * (rhowba_r[i]* rhowcd_r[i] - rhowba_i[i] * rhowcd_i[i]);


                    rr_123_i[i] = 2.0 * (rhowda_i[i]* rhowbc_r[i] + rhowda_r[i] * rhowbc_i[i])

                                + 2.0 * (rhowca_i[i]* rhowbd_r[i] + rhowca_r[i] * rhowbd_i[i])

                                + 2.0 * (rhowba_i[i]* rhowcd_r[i] + rhowba_r[i] * rhowcd_i[i]);

                }
            }
        }
        if (fstr::upcase(quadMode) == "QRF")
        {
            // This routine computes the Fcb for the general quadratic response function

            for (int32_t j = 0; j < numdens / 2; j++)
            {
                auto rhorho_r = gam(2 * j);

                auto rhorho_i = gam(2 * j + 1);

                auto rhow1a_r = rwDensityGrid.alphaDensity(4 * j);

                auto rhow1a_i = rwDensityGrid.alphaDensity(4 * j + 1);

                auto rhow2a_r = rwDensityGrid.alphaDensity(4 * j + 2);

                auto rhow2a_i = rwDensityGrid.alphaDensity(4 * j + 3);

                for (int32_t i = 0; i < npoints; i++)
                {
                    rhorho_r[i] = 2.0 * (rhow1a_r[i] * rhow2a_r[i] - rhow1a_i[i] * rhow2a_i[i]);

                    rhorho_i[i] = 2.0 * (rhow1a_r[i] * rhow2a_i[i] + rhow1a_i[i] * rhow2a_r[i]);
                }
            }
        }
    }
    if (xcFuncType == xcfun::gga)
    {
        if (fstr::upcase(quadMode) == "REDTPA_I")
        {
           for (int32_t j = 0; j < numdens / 6; j++)
            {

               auto gam_sig_xx_r = gam(6 * j  + 0);
               auto gam_sig_yy_r = gam(6 * j  + 1);
               auto gam_sig_zz_r = gam(6 * j  + 2);
               auto gam_sig_xy_r = gam(6 * j  + 3);
               auto gam_sig_xz_r = gam(6 * j  + 4);
               auto gam_sig_yz_r = gam(6 * j  + 5);
               auto gam_sig_xx_x_r = gamX(6 * j  + 0);
               auto gam_sig_xx_y_r = gamY(6 * j  + 0);
               auto gam_sig_xx_z_r = gamZ(6 * j  + 0);
               auto gam_sig_yy_x_r = gamX(6 * j  + 1);
               auto gam_sig_yy_y_r = gamY(6 * j  + 1);
               auto gam_sig_yy_z_r = gamZ(6 * j  + 1);
               auto gam_sig_zz_x_r = gamX(6 * j  + 2);
               auto gam_sig_zz_y_r = gamY(6 * j  + 2);
               auto gam_sig_zz_z_r = gamZ(6 * j  + 2);
               auto gam_sig_xy_x_r = gamX(6 * j  + 3);
               auto gam_sig_xy_y_r = gamY(6 * j  + 3);
               auto gam_sig_xy_z_r = gamZ(6 * j  + 3);
               auto gam_sig_xz_x_r = gamX(6 * j  + 4);
               auto gam_sig_xz_y_r = gamY(6 * j  + 4);
               auto gam_sig_xz_z_r = gamZ(6 * j  + 4);
               auto gam_sig_yz_x_r = gamX(6 * j  + 5);
               auto gam_sig_yz_y_r = gamY(6 * j  + 5);
               auto gam_sig_yz_z_r = gamZ(6 * j  + 5);
               auto gam_sig_xx_xx_r = gamXX(6 * j  + 0);
               auto gam_sig_xx_xy_r = gamXY(6 * j  + 0);
               auto gam_sig_xx_xz_r = gamXZ(6 * j  + 0);
               auto gam_sig_xx_yx_r = gamYX(6 * j  + 0);
               auto gam_sig_xx_yy_r = gamYY(6 * j  + 0);
               auto gam_sig_xx_yz_r = gamYZ(6 * j  + 0);
               auto gam_sig_xx_zx_r = gamZX(6 * j  + 0);
               auto gam_sig_xx_zy_r = gamZY(6 * j  + 0);
               auto gam_sig_xx_zz_r = gamZZ(6 * j  + 0);
               auto gam_sig_yy_xx_r = gamXX(6 * j  + 1);
               auto gam_sig_yy_xy_r = gamXY(6 * j  + 1);
               auto gam_sig_yy_xz_r = gamXZ(6 * j  + 1);
               auto gam_sig_yy_yx_r = gamYX(6 * j  + 1);
               auto gam_sig_yy_yy_r = gamYY(6 * j  + 1);
               auto gam_sig_yy_yz_r = gamYZ(6 * j  + 1);
               auto gam_sig_yy_zx_r = gamZX(6 * j  + 1);
               auto gam_sig_yy_zy_r = gamZY(6 * j  + 1);
               auto gam_sig_yy_zz_r = gamZZ(6 * j  + 1);
               auto gam_sig_zz_xx_r = gamXX(6 * j  + 2);
               auto gam_sig_zz_xy_r = gamXY(6 * j  + 2);
               auto gam_sig_zz_xz_r = gamXZ(6 * j  + 2);
               auto gam_sig_zz_yx_r = gamYX(6 * j  + 2);
               auto gam_sig_zz_yy_r = gamYY(6 * j  + 2);
               auto gam_sig_zz_yz_r = gamYZ(6 * j  + 2);
               auto gam_sig_zz_zx_r = gamZX(6 * j  + 2);
               auto gam_sig_zz_zy_r = gamZY(6 * j  + 2);
               auto gam_sig_zz_zz_r = gamZZ(6 * j  + 2);
               auto gam_sig_xy_xx_r = gamXX(6 * j  + 3);
               auto gam_sig_xy_xy_r = gamXY(6 * j  + 3);
               auto gam_sig_xy_xz_r = gamXZ(6 * j  + 3);
               auto gam_sig_xy_yx_r = gamYX(6 * j  + 3);
               auto gam_sig_xy_yy_r = gamYY(6 * j  + 3);
               auto gam_sig_xy_yz_r = gamYZ(6 * j  + 3);
               auto gam_sig_xy_zx_r = gamZX(6 * j  + 3);
               auto gam_sig_xy_zy_r = gamZY(6 * j  + 3);
               auto gam_sig_xy_zz_r = gamZZ(6 * j  + 3);
               auto gam_sig_xz_xx_r = gamXX(6 * j  + 4);
               auto gam_sig_xz_xy_r = gamXY(6 * j  + 4);
               auto gam_sig_xz_xz_r = gamXZ(6 * j  + 4);
               auto gam_sig_xz_yx_r = gamYX(6 * j  + 4);
               auto gam_sig_xz_yy_r = gamYY(6 * j  + 4);
               auto gam_sig_xz_yz_r = gamYZ(6 * j  + 4);
               auto gam_sig_xz_zx_r = gamZX(6 * j  + 4);
               auto gam_sig_xz_zy_r = gamZY(6 * j  + 4);
               auto gam_sig_xz_zz_r = gamZZ(6 * j  + 4);
               auto gam_sig_yz_xx_r = gamXX(6 * j  + 5);
               auto gam_sig_yz_xy_r = gamXY(6 * j  + 5);
               auto gam_sig_yz_xz_r = gamXZ(6 * j  + 5);
               auto gam_sig_yz_yx_r = gamYX(6 * j  + 5);
               auto gam_sig_yz_yy_r = gamYY(6 * j  + 5);
               auto gam_sig_yz_yz_r = gamYZ(6 * j  + 5);
               auto gam_sig_yz_zx_r = gamZX(6 * j  + 5);
               auto gam_sig_yz_zy_r = gamZY(6 * j  + 5);
               auto gam_sig_yz_zz_r = gamZZ(6 * j  + 5);

                // First-order densities First variable is for the component of the dipole, the second variable is for the component of the gradient
                // B means positive frequency and C means negative frequency

                auto rhoBx_r = rwDensityGrid.alphaDensity(6 * j + 0);
                auto rhoBx_i = rwDensityGrid.alphaDensity(6 * j + 1);
                auto rhoBy_r = rwDensityGrid.alphaDensity(6 * j + 2);
                auto rhoBy_i = rwDensityGrid.alphaDensity(6 * j + 3);
                auto rhoBz_r = rwDensityGrid.alphaDensity(6 * j + 4);
                auto rhoBz_i = rwDensityGrid.alphaDensity(6 * j + 5);
                auto gradBx_x_r = rwDensityGrid.alphaDensityGradientX(6 * j + 0);
                auto gradBx_x_i = rwDensityGrid.alphaDensityGradientX(6 * j + 1);
                auto gradBx_y_r = rwDensityGrid.alphaDensityGradientY(6 * j + 0);
                auto gradBx_y_i = rwDensityGrid.alphaDensityGradientY(6 * j + 1);
                auto gradBx_z_r = rwDensityGrid.alphaDensityGradientZ(6 * j + 0);
                auto gradBx_z_i = rwDensityGrid.alphaDensityGradientZ(6 * j + 1);
                auto gradBy_x_r = rwDensityGrid.alphaDensityGradientX(6 * j + 2);
                auto gradBy_x_i = rwDensityGrid.alphaDensityGradientX(6 * j + 3);
                auto gradBy_y_r = rwDensityGrid.alphaDensityGradientY(6 * j + 2);
                auto gradBy_y_i = rwDensityGrid.alphaDensityGradientY(6 * j + 3);
                auto gradBy_z_r = rwDensityGrid.alphaDensityGradientZ(6 * j + 2);
                auto gradBy_z_i = rwDensityGrid.alphaDensityGradientZ(6 * j + 3);
                auto gradBz_x_r = rwDensityGrid.alphaDensityGradientX(6 * j + 4);
                auto gradBz_x_i = rwDensityGrid.alphaDensityGradientX(6 * j + 5);
                auto gradBz_y_r = rwDensityGrid.alphaDensityGradientY(6 * j + 4);
                auto gradBz_y_i = rwDensityGrid.alphaDensityGradientY(6 * j + 5);
                auto gradBz_z_r = rwDensityGrid.alphaDensityGradientZ(6 * j + 4);
                auto gradBz_z_i = rwDensityGrid.alphaDensityGradientZ(6 * j + 5);
                

                for (int32_t i = 0; i < npoints; i++)

                {
                    double sig_term_r =  + 6.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBx_r[i],rhoBx_i[i])
                                         + 6.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoBy_r[i],rhoBy_i[i])
                                         + 6.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    double sig_term_x_r =  + 12.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],rhoBx_r[i],rhoBx_i[i])
                                           + 12.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],rhoBy_r[i],rhoBy_i[i])
                                           + 12.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],rhoBz_r[i],rhoBz_i[i]);                 

                    double sig_term_y_r =  + 12.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],rhoBx_r[i],rhoBx_i[i])
                                           + 12.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],rhoBy_r[i],rhoBy_i[i])
                                           + 12.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],rhoBz_r[i],rhoBz_i[i]);                 

                    double sig_term_z_r =  + 12.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],rhoBx_r[i],rhoBx_i[i])
                                           + 12.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],rhoBy_r[i],rhoBy_i[i])
                                           + 12.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],rhoBz_r[i],rhoBz_i[i]);                 
         

                    double sig_term_xx_r =  + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                            + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                            + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]);                   

                    double sig_term_xy_r =  + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                            + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                            + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]);                   

                    double sig_term_xz_r =  + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                            + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                            + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]);                   

                    double sig_term_yx_r =  + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                            + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                            + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]);                   

                    double sig_term_yy_r =  + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                            + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                            + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]);                   

                    double sig_term_yz_r =  + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                            + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                            + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]);                   

                    double sig_term_zx_r =  + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                            + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                            + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]);                   

                    double sig_term_zy_r =  + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                            + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                            + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]);                   

                    double sig_term_zz_r =  + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                            + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                            + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]);                   
                    
                    gam_sig_xx_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBx_r[i],rhoBx_i[i]) + sig_term_r;                 

                    gam_sig_yy_r[i] = 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoBy_r[i],rhoBy_i[i]) + sig_term_r;                 

                    gam_sig_zz_r[i] = 12.0 * prod2_r(rhoBz_r[i],rhoBz_i[i],rhoBz_r[i],rhoBz_i[i]) + sig_term_r;                 

                    gam_sig_xy_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBy_r[i],rhoBy_i[i]);                  

                    gam_sig_xz_r[i] = 12.0 * prod2_r(rhoBx_r[i],rhoBx_i[i],rhoBz_r[i],rhoBz_i[i]);                  

                    gam_sig_yz_r[i] = 12.0 * prod2_r(rhoBy_r[i],rhoBy_i[i],rhoBz_r[i],rhoBz_i[i]);                  
                    
                    gam_sig_xx_x_r[i] = 12.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],rhoBx_r[i],rhoBx_i[i])
                                       +12.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],rhoBx_r[i],rhoBx_i[i]) + sig_term_x_r;                   

                    gam_sig_xx_y_r[i] = 12.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],rhoBx_r[i],rhoBx_i[i])
                                       +12.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],rhoBx_r[i],rhoBx_i[i]) + sig_term_y_r;                   

                    gam_sig_xx_z_r[i] = 12.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],rhoBx_r[i],rhoBx_i[i])
                                       +12.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],rhoBx_r[i],rhoBx_i[i]) + sig_term_z_r;                   

                    gam_sig_yy_x_r[i] = 12.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],rhoBy_r[i],rhoBy_i[i])
                                       +12.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],rhoBy_r[i],rhoBy_i[i]) + sig_term_x_r;                   

                    gam_sig_yy_y_r[i] = 12.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],rhoBy_r[i],rhoBy_i[i])
                                       +12.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],rhoBy_r[i],rhoBy_i[i]) + sig_term_y_r;                   

                    gam_sig_yy_z_r[i] = 12.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],rhoBy_r[i],rhoBy_i[i])
                                       +12.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],rhoBy_r[i],rhoBy_i[i]) + sig_term_z_r;                   

                    gam_sig_zz_x_r[i] = 12.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],rhoBz_r[i],rhoBz_i[i])
                                       +12.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],rhoBz_r[i],rhoBz_i[i]) + sig_term_x_r;                   

                    gam_sig_zz_y_r[i] = 12.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],rhoBz_r[i],rhoBz_i[i])
                                       +12.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],rhoBz_r[i],rhoBz_i[i]) + sig_term_y_r;                   

                    gam_sig_zz_z_r[i] = 12.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],rhoBz_r[i],rhoBz_i[i])
                                       +12.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],rhoBz_r[i],rhoBz_i[i]) + sig_term_z_r;                   

                    gam_sig_xy_x_r[i] = 12.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],rhoBy_r[i],rhoBy_i[i])
                                      + 12.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_xy_y_r[i] = 12.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],rhoBy_r[i],rhoBy_i[i])
                                      + 12.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_xy_z_r[i] = 12.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],rhoBy_r[i],rhoBy_i[i])
                                      + 12.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_xz_x_r[i] = 12.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_xz_y_r[i] = 12.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_xz_z_r[i] = 12.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],rhoBx_r[i],rhoBx_i[i]);
                    gam_sig_yz_x_r[i] = 12.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_sig_yz_y_r[i] = 12.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],rhoBy_r[i],rhoBy_i[i]);
                    gam_sig_yz_z_r[i] = 12.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],rhoBz_r[i],rhoBz_i[i])
                                      + 12.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],rhoBy_r[i],rhoBy_i[i]);
                    
                    gam_sig_xx_xx_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                       + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_x_r[i],gradBx_x_i[i]) + sig_term_xx_r;                    

                    gam_sig_xx_xy_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                       + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_y_r[i],gradBx_y_i[i]) + sig_term_xy_r;                    

                    gam_sig_xx_xz_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                       + 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBx_z_r[i],gradBx_z_i[i]) + sig_term_xz_r;                    

                    gam_sig_xx_yx_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                       + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_x_r[i],gradBx_x_i[i]) + sig_term_yx_r;                    

                    gam_sig_xx_yy_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                       + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_y_r[i],gradBx_y_i[i]) + sig_term_yy_r;                    

                    gam_sig_xx_yz_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                       + 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBx_z_r[i],gradBx_z_i[i]) + sig_term_yz_r;                    

                    gam_sig_xx_zx_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_x_r[i],gradBx_x_i[i])
                                       + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_x_r[i],gradBx_x_i[i]) + sig_term_zx_r;                    

                    gam_sig_xx_zy_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_y_r[i],gradBx_y_i[i])
                                       + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_y_r[i],gradBx_y_i[i]) + sig_term_zy_r;                    

                    gam_sig_xx_zz_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_z_r[i],gradBx_z_i[i])
                                       + 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBx_z_r[i],gradBx_z_i[i]) + sig_term_zz_r;                    

                    gam_sig_yy_xx_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_x_r[i],gradBy_x_i[i]) + sig_term_xx_r;                    

                    gam_sig_yy_xy_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_y_r[i],gradBy_y_i[i]) + sig_term_xy_r;                    

                    gam_sig_yy_xz_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBy_z_r[i],gradBy_z_i[i]) + sig_term_xz_r;                    

                    gam_sig_yy_yx_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_x_r[i],gradBy_x_i[i]) + sig_term_yx_r;                    

                    gam_sig_yy_yy_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_y_r[i],gradBy_y_i[i]) + sig_term_yy_r;                    

                    gam_sig_yy_yz_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBy_z_r[i],gradBy_z_i[i]) + sig_term_yz_r;                    

                    gam_sig_yy_zx_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_x_r[i],gradBy_x_i[i]) + sig_term_zx_r;                    

                    gam_sig_yy_zy_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_y_r[i],gradBy_y_i[i]) + sig_term_zy_r;                    

                    gam_sig_yy_zz_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBy_z_r[i],gradBy_z_i[i]) + sig_term_zz_r;                    

                    gam_sig_zz_xx_r[i] = 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_x_r[i],gradBz_x_i[i]) + sig_term_xx_r;                    

                    gam_sig_zz_xy_r[i] = 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_y_r[i],gradBz_y_i[i]) + sig_term_xy_r;                    

                    gam_sig_zz_xz_r[i] = 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBz_z_r[i],gradBz_z_i[i]) + sig_term_xz_r;                    

                    gam_sig_zz_yx_r[i] = 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_x_r[i],gradBz_x_i[i]) + sig_term_yx_r;                    

                    gam_sig_zz_yy_r[i] = 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_y_r[i],gradBz_y_i[i]) + sig_term_yy_r;                    

                    gam_sig_zz_yz_r[i] = 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBz_z_r[i],gradBz_z_i[i]) + sig_term_yz_r;                    

                    gam_sig_zz_zx_r[i] = 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_x_r[i],gradBz_x_i[i]) + sig_term_zx_r;                    

                    gam_sig_zz_zy_r[i] = 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_y_r[i],gradBz_y_i[i]) + sig_term_zy_r;                    

                    gam_sig_zz_zz_r[i] = 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBz_z_r[i],gradBz_z_i[i]) + sig_term_zz_r;                    

                    gam_sig_xy_xx_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xy_xy_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xy_xz_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_xy_yx_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xy_yy_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xy_yz_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_xy_zx_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBy_x_r[i],gradBy_x_i[i])
                                       + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xy_zy_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBy_y_r[i],gradBy_y_i[i])
                                       + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xy_zz_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBy_z_r[i],gradBy_z_i[i])
                                       + 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_xz_xx_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xz_xy_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xz_xz_r[i] = 6.0 * prod2_r(gradBx_x_r[i],gradBx_x_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_xz_yx_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xz_yy_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xz_yz_r[i] = 6.0 * prod2_r(gradBx_y_r[i],gradBx_y_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_xz_zx_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBx_x_r[i],gradBx_x_i[i]);
                    gam_sig_xz_zy_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBx_y_r[i],gradBx_y_i[i]);
                    gam_sig_xz_zz_r[i] = 6.0 * prod2_r(gradBx_z_r[i],gradBx_z_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBx_z_r[i],gradBx_z_i[i]);
                    gam_sig_yz_xx_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBy_x_r[i],gradBy_x_i[i]);
                    gam_sig_yz_xy_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBy_y_r[i],gradBy_y_i[i]);
                    gam_sig_yz_xz_r[i] = 6.0 * prod2_r(gradBy_x_r[i],gradBy_x_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_x_r[i],gradBz_x_i[i],gradBy_z_r[i],gradBy_z_i[i]);
                    gam_sig_yz_yx_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBy_x_r[i],gradBy_x_i[i]);
                    gam_sig_yz_yy_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBy_y_r[i],gradBy_y_i[i]);
                    gam_sig_yz_yz_r[i] = 6.0 * prod2_r(gradBy_y_r[i],gradBy_y_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_y_r[i],gradBz_y_i[i],gradBy_z_r[i],gradBy_z_i[i]);
                    gam_sig_yz_zx_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBz_x_r[i],gradBz_x_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBy_x_r[i],gradBy_x_i[i]);
                    gam_sig_yz_zy_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBz_y_r[i],gradBz_y_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBy_y_r[i],gradBy_y_i[i]);
                    gam_sig_yz_zz_r[i] = 6.0 * prod2_r(gradBy_z_r[i],gradBy_z_i[i],gradBz_z_r[i],gradBz_z_i[i])
                                       + 6.0 * prod2_r(gradBz_z_r[i],gradBz_z_i[i],gradBy_z_r[i],gradBy_z_i[i]);
                              
                }
            }
        }
        if (fstr::upcase(quadMode) == "REDTPA_II") 
        {
            for (int32_t j = 0; j < numdens / 6; j++)
            {
                
                auto F_x_gam_r = gam(6 * j);

                auto F_x_gam_i = gam(6 * j + 1);

                auto F_x_gamX_r = gamX(6 * j);

                auto F_x_gamX_i = gamX(6 * j + 1);

                auto F_x_gamY_r = gamY(6 * j);

                auto F_x_gamY_i = gamY(6 * j + 1);

                auto F_x_gamZ_r = gamZ(6 * j);

                auto F_x_gamZ_i = gamZ(6 * j + 1);

                auto F_x_gamXX_r = gamXX(6 * j);

                auto F_x_gamXX_i = gamXX(6 * j + 1);

                auto F_x_gamXY_r = gamXY(6 * j);

                auto F_x_gamXY_i = gamXY(6 * j + 1);

                auto F_x_gamXZ_r = gamXZ(6 * j);

                auto F_x_gamXZ_i = gamXZ(6 * j + 1);

                auto F_x_gamYX_r = gamYX(6 * j);

                auto F_x_gamYX_i = gamYX(6 * j + 1);

                auto F_x_gamYY_r = gamYY(6 * j);

                auto F_x_gamYY_i = gamYY(6 * j + 1);

                auto F_x_gamYZ_r = gamYZ(6 * j);

                auto F_x_gamYZ_i = gamYZ(6 * j + 1);

                auto F_x_gamZX_r = gamZX(6 * j);

                auto F_x_gamZX_i = gamZX(6 * j + 1);

                auto F_x_gamZY_r = gamZY(6 * j);

                auto F_x_gamZY_i = gamZY(6 * j + 1);

                auto F_x_gamZZ_r = gamZZ(6 * j);

                auto F_x_gamZZ_i = gamZZ(6 * j + 1);

                // Fy

                auto F_y_gam_r = gam(6 * j  + 2);
 
                auto F_y_gam_i = gam(6 * j  + 3);

                auto F_y_gamX_r = gamX(6 * j  + 2);

                auto F_y_gamX_i = gamX(6 * j  + 3);

                auto F_y_gamY_r = gamY(6 * j  + 2);

                auto F_y_gamY_i = gamY(6 * j  + 3);

                auto F_y_gamZ_r = gamZ(6 * j  + 2);

                auto F_y_gamZ_i = gamZ(6 * j  + 3);

                auto F_y_gamXX_r = gamXX(6 * j  + 2);

                auto F_y_gamXX_i = gamXX(6 * j  + 3);

                auto F_y_gamXY_r = gamXY(6 * j  + 2);

                auto F_y_gamXY_i = gamXY(6 * j  + 3);

                auto F_y_gamXZ_r = gamXZ(6 * j  + 2);

                auto F_y_gamXZ_i = gamXZ(6 * j  + 3);

                auto F_y_gamYX_r = gamYX(6 * j  + 2);

                auto F_y_gamYX_i = gamYX(6 * j  + 3);

                auto F_y_gamYY_r = gamYY(6 * j  + 2);

                auto F_y_gamYY_i = gamYY(6 * j  + 3);

                auto F_y_gamYZ_r = gamYZ(6 * j  + 2);

                auto F_y_gamYZ_i = gamYZ(6 * j  + 3);

                auto F_y_gamZX_r = gamZX(6 * j  + 2);

                auto F_y_gamZX_i = gamZX(6 * j  + 3);

                auto F_y_gamZY_r = gamZY(6 * j  + 2);

                auto F_y_gamZY_i = gamZY(6 * j  + 3);

                auto F_y_gamZZ_r = gamZZ(6 * j  + 2);

                auto F_y_gamZZ_i = gamZZ(6 * j  + 3);

                // Fz

                auto F_z_gam_r = gam(6 * j  + 4);

                auto F_z_gam_i = gam(6 * j  + 5);

                auto F_z_gamX_r = gamX(6 * j  + 4);

                auto F_z_gamX_i = gamX(6 * j  + 5);

                auto F_z_gamY_r = gamY(6 * j  + 4);

                auto F_z_gamY_i = gamY(6 * j  + 5);

                auto F_z_gamZ_r = gamZ(6 * j  + 4);

                auto F_z_gamZ_i = gamZ(6 * j  + 5);

                auto F_z_gamXX_r = gamXX(6 * j  + 4);

                auto F_z_gamXX_i = gamXX(6 * j  + 5);

                auto F_z_gamXY_r = gamXY(6 * j  + 4);

                auto F_z_gamXY_i = gamXY(6 * j  + 5);

                auto F_z_gamXZ_r = gamXZ(6 * j  + 4);

                auto F_z_gamXZ_i = gamXZ(6 * j  + 5);

                auto F_z_gamYX_r = gamYX(6 * j  + 4);

                auto F_z_gamYX_i = gamYX(6 * j  + 5);

                auto F_z_gamYY_r = gamYY(6 * j  + 4);

                auto F_z_gamYY_i = gamYY(6 * j  + 5);

                auto F_z_gamYZ_r = gamYZ(6 * j  + 4);

                auto F_z_gamYZ_i = gamYZ(6 * j  + 5);

                auto F_z_gamZX_r = gamZX(6 * j  + 4);

                auto F_z_gamZX_i = gamZX(6 * j  + 5);

                auto F_z_gamZY_r = gamZY(6 * j  + 4);

                auto F_z_gamZY_i = gamZY(6 * j  + 5);

                auto F_z_gamZZ_r = gamZZ(6 * j  + 4);

                auto F_z_gamZZ_i = gamZZ(6 * j  + 5);


                // Perturbed densities

                auto rhoCx_rho_r = rwDensityGrid.alphaDensity(18 * j + 0);
                auto rhoCx_grad_x_r = rwDensityGrid.alphaDensityGradientX(18 * j + 0);
                auto rhoCx_grad_y_r = rwDensityGrid.alphaDensityGradientY(18 * j + 0);
                auto rhoCx_grad_z_r = rwDensityGrid.alphaDensityGradientZ(18 * j + 0);
                auto rhoCx_rho_i = rwDensityGrid.alphaDensity(18 * j + 1);
                auto rhoCx_grad_x_i = rwDensityGrid.alphaDensityGradientX(18 * j + 1);
                auto rhoCx_grad_y_i = rwDensityGrid.alphaDensityGradientY(18 * j + 1);
                auto rhoCx_grad_z_i = rwDensityGrid.alphaDensityGradientZ(18 * j + 1);

                auto rhoCy_rho_r = rwDensityGrid.alphaDensity(18 * j + 2);
                auto rhoCy_grad_x_r = rwDensityGrid.alphaDensityGradientX(18 * j + 2);
                auto rhoCy_grad_y_r = rwDensityGrid.alphaDensityGradientY(18 * j + 2);
                auto rhoCy_grad_z_r = rwDensityGrid.alphaDensityGradientZ(18 * j + 2);
                auto rhoCy_rho_i = rwDensityGrid.alphaDensity(18 * j + 3);
                auto rhoCy_grad_x_i = rwDensityGrid.alphaDensityGradientX(18 * j + 3);
                auto rhoCy_grad_y_i = rwDensityGrid.alphaDensityGradientY(18 * j + 3);
                auto rhoCy_grad_z_i = rwDensityGrid.alphaDensityGradientZ(18 * j + 3);

                auto rhoCz_rho_r = rwDensityGrid.alphaDensity(18 * j + 4);
                auto rhoCz_grad_x_r = rwDensityGrid.alphaDensityGradientX(18 * j  + 4);
                auto rhoCz_grad_y_r = rwDensityGrid.alphaDensityGradientY(18 * j  + 4);
                auto rhoCz_grad_z_r = rwDensityGrid.alphaDensityGradientZ(18 * j  + 4);
                auto rhoCz_rho_i = rwDensityGrid.alphaDensity(18 * j + 5);
                auto rhoCz_grad_x_i = rwDensityGrid.alphaDensityGradientX(18 * j + 5);
                auto rhoCz_grad_y_i = rwDensityGrid.alphaDensityGradientY(18 * j + 5);
                auto rhoCz_grad_z_i = rwDensityGrid.alphaDensityGradientZ(18 * j + 5);

                auto Rho_sig_xx_rho_r = rwDensityGrid.alphaDensity(18 * j + 6);
                auto Rho_sig_xx_grad_x_r = rwDensityGrid.alphaDensityGradientX(18 * j  + 6);
                auto Rho_sig_xx_grad_y_r = rwDensityGrid.alphaDensityGradientY(18 * j  + 6);
                auto Rho_sig_xx_grad_z_r = rwDensityGrid.alphaDensityGradientZ(18 * j  + 6);
                auto Rho_sig_xx_rho_i = rwDensityGrid.alphaDensity(18 * j + 7);
                auto Rho_sig_xx_grad_x_i = rwDensityGrid.alphaDensityGradientX(18 * j + 7);
                auto Rho_sig_xx_grad_y_i = rwDensityGrid.alphaDensityGradientY(18 * j + 7);
                auto Rho_sig_xx_grad_z_i = rwDensityGrid.alphaDensityGradientZ(18 * j + 7);

                auto Rho_sig_yy_rho_r = rwDensityGrid.alphaDensity(18 * j + 8);
                auto Rho_sig_yy_grad_x_r = rwDensityGrid.alphaDensityGradientX(18 * j  + 8);
                auto Rho_sig_yy_grad_y_r = rwDensityGrid.alphaDensityGradientY(18 * j  + 8);
                auto Rho_sig_yy_grad_z_r = rwDensityGrid.alphaDensityGradientZ(18 * j  + 8);
                auto Rho_sig_yy_rho_i = rwDensityGrid.alphaDensity(18 * j + 9);
                auto Rho_sig_yy_grad_x_i = rwDensityGrid.alphaDensityGradientX(18 * j + 9);
                auto Rho_sig_yy_grad_y_i = rwDensityGrid.alphaDensityGradientY(18 * j + 9);
                auto Rho_sig_yy_grad_z_i = rwDensityGrid.alphaDensityGradientZ(18 * j + 9);

                auto Rho_sig_zz_rho_r = rwDensityGrid.alphaDensity(18 * j + 10);
                auto Rho_sig_zz_grad_x_r = rwDensityGrid.alphaDensityGradientX(18 * j  + 10);
                auto Rho_sig_zz_grad_y_r = rwDensityGrid.alphaDensityGradientY(18 * j  + 10);
                auto Rho_sig_zz_grad_z_r = rwDensityGrid.alphaDensityGradientZ(18 * j  + 10);
                auto Rho_sig_zz_rho_i = rwDensityGrid.alphaDensity(18 * j + 11);
                auto Rho_sig_zz_grad_x_i = rwDensityGrid.alphaDensityGradientX(18 * j + 11);
                auto Rho_sig_zz_grad_y_i = rwDensityGrid.alphaDensityGradientY(18 * j + 11);
                auto Rho_sig_zz_grad_z_i = rwDensityGrid.alphaDensityGradientZ(18 * j + 11);

                auto Rho_sig_xy_rho_r = rwDensityGrid.alphaDensity(18 * j + 12);
                auto Rho_sig_xy_grad_x_r = rwDensityGrid.alphaDensityGradientX(18 * j  + 12);
                auto Rho_sig_xy_grad_y_r = rwDensityGrid.alphaDensityGradientY(18 * j  + 12);
                auto Rho_sig_xy_grad_z_r = rwDensityGrid.alphaDensityGradientZ(18 * j  + 12);
                auto Rho_sig_xy_rho_i = rwDensityGrid.alphaDensity(18 * j + 13);
                auto Rho_sig_xy_grad_x_i = rwDensityGrid.alphaDensityGradientX(18 * j + 13);
                auto Rho_sig_xy_grad_y_i = rwDensityGrid.alphaDensityGradientY(18 * j + 13);
                auto Rho_sig_xy_grad_z_i = rwDensityGrid.alphaDensityGradientZ(18 * j + 13);

                auto Rho_sig_xz_rho_r = rwDensityGrid.alphaDensity(18 * j + 14);
                auto Rho_sig_xz_grad_x_r = rwDensityGrid.alphaDensityGradientX(18 * j  + 14);
                auto Rho_sig_xz_grad_y_r = rwDensityGrid.alphaDensityGradientY(18 * j  + 14);
                auto Rho_sig_xz_grad_z_r = rwDensityGrid.alphaDensityGradientZ(18 * j  + 14);
                auto Rho_sig_xz_rho_i = rwDensityGrid.alphaDensity(18 * j + 15);
                auto Rho_sig_xz_grad_x_i = rwDensityGrid.alphaDensityGradientX(18 * j + 15);
                auto Rho_sig_xz_grad_y_i = rwDensityGrid.alphaDensityGradientY(18 * j + 15);
                auto Rho_sig_xz_grad_z_i = rwDensityGrid.alphaDensityGradientZ(18 * j + 15);

                auto Rho_sig_yz_rho_r = rwDensityGrid.alphaDensity(18 * j + 16);
                auto Rho_sig_yz_grad_x_r = rwDensityGrid.alphaDensityGradientX(18 * j  + 16);
                auto Rho_sig_yz_grad_y_r = rwDensityGrid.alphaDensityGradientY(18 * j  + 16);
                auto Rho_sig_yz_grad_z_r = rwDensityGrid.alphaDensityGradientZ(18 * j  + 16);
                auto Rho_sig_yz_rho_i = rwDensityGrid.alphaDensity(18 * j + 17);
                auto Rho_sig_yz_grad_x_i = rwDensityGrid.alphaDensityGradientX(18 * j + 17);
                auto Rho_sig_yz_grad_y_i = rwDensityGrid.alphaDensityGradientY(18 * j + 17);
                auto Rho_sig_yz_grad_z_i = rwDensityGrid.alphaDensityGradientZ(18 * j + 17);

                for (int32_t i = 0; i < npoints; i++)
                {
                    

                    // Sig_xx rho_x

                    double rhow1_r = Rho_sig_xx_rho_r[i];

                    double rxw1_r = Rho_sig_xx_grad_x_r[i];

                    double ryw1_r = Rho_sig_xx_grad_y_r[i];

                    double rzw1_r = Rho_sig_xx_grad_z_r[i];

                    double rhow1_i = Rho_sig_xx_rho_i[i];

                    double rxw1_i = Rho_sig_xx_grad_x_i[i];

                    double ryw1_i = Rho_sig_xx_grad_y_i[i];

                    double rzw1_i = Rho_sig_xx_grad_z_i[i];


                    double rhow2_r = rhoCx_rho_r[i];

                    double rxw2_r = rhoCx_grad_x_r[i];

                    double ryw2_r = rhoCx_grad_y_r[i];

                    double rzw2_r = rhoCx_grad_z_r[i];

                    double rhow2_i = rhoCx_rho_i[i];

                    double rxw2_i = rhoCx_grad_x_i[i];

                    double ryw2_i = rhoCx_grad_y_i[i];

                    double rzw2_i = rhoCx_grad_z_i[i];
                    

                    F_x_gam_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_x_gam_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_x_gamX_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_x_gamX_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_x_gamY_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_x_gamY_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_x_gamZ_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_x_gamZ_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_x_gamXX_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_x_gamXX_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_x_gamXY_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_x_gamXY_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_x_gamXZ_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_x_gamXZ_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_x_gamYX_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_x_gamYX_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_x_gamYY_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_x_gamYY_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_x_gamYZ_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_x_gamYZ_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_x_gamZX_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_x_gamZX_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_x_gamZY_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_x_gamZY_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_x_gamZZ_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_x_gamZZ_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    // Sig_xy rho_y

                    rhow1_r = Rho_sig_xy_rho_r[i];

                    rxw1_r = Rho_sig_xy_grad_x_r[i];

                    ryw1_r = Rho_sig_xy_grad_y_r[i];

                    rzw1_r = Rho_sig_xy_grad_z_r[i];

                    rhow1_i = Rho_sig_xy_rho_i[i];

                    rxw1_i = Rho_sig_xy_grad_x_i[i];

                    ryw1_i = Rho_sig_xy_grad_y_i[i];

                    rzw1_i = Rho_sig_xy_grad_z_i[i];
                

                    rhow2_r = rhoCy_rho_r[i];

                    rxw2_r = rhoCy_grad_x_r[i];

                    ryw2_r = rhoCy_grad_y_r[i];

                    rzw2_r = rhoCy_grad_z_r[i];

                    rhow2_i = rhoCy_rho_i[i];

                    rxw2_i = rhoCy_grad_x_i[i];

                    ryw2_i = rhoCy_grad_y_i[i];

                    rzw2_i = rhoCy_grad_z_i[i];
                        

                    F_x_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_x_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_x_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_x_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_x_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_x_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_x_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_x_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_x_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_x_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_x_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_x_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_x_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_x_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_x_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_x_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_x_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_x_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_x_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_x_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_x_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_x_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_x_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_x_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_x_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_x_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;



                    // Sig_xz rho_z

                    rhow1_r = Rho_sig_xz_rho_r[i];

                    rxw1_r = Rho_sig_xz_grad_x_r[i];

                    ryw1_r = Rho_sig_xz_grad_y_r[i];

                    rzw1_r = Rho_sig_xz_grad_z_r[i];

                    rhow1_i = Rho_sig_xz_rho_i[i];

                    rxw1_i = Rho_sig_xz_grad_x_i[i];

                    ryw1_i = Rho_sig_xz_grad_y_i[i];

                    rzw1_i = Rho_sig_xz_grad_z_i[i];
                

                    rhow2_r = rhoCz_rho_r[i];

                    rxw2_r = rhoCz_grad_x_r[i];

                    ryw2_r = rhoCz_grad_y_r[i];

                    rzw2_r = rhoCz_grad_z_r[i];

                    rhow2_i = rhoCz_rho_i[i];

                    rxw2_i = rhoCz_grad_x_i[i];

                    ryw2_i = rhoCz_grad_y_i[i];

                    rzw2_i = rhoCz_grad_z_i[i];
                        

                    F_x_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_x_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_x_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_x_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_x_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_x_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_x_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_x_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_x_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_x_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_x_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_x_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_x_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_x_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_x_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_x_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_x_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_x_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_x_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_x_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_x_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_x_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_x_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_x_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_x_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_x_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;

                    

                    // Y case

                    
                    // Sig_yx rho_x

                    rhow1_r = Rho_sig_xy_rho_r[i];

                    rxw1_r = Rho_sig_xy_grad_x_r[i];

                    ryw1_r = Rho_sig_xy_grad_y_r[i];

                    rzw1_r = Rho_sig_xy_grad_z_r[i];

                    rhow1_i = Rho_sig_xy_rho_i[i];

                    rxw1_i = Rho_sig_xy_grad_x_i[i];

                    ryw1_i = Rho_sig_xy_grad_y_i[i];

                    rzw1_i = Rho_sig_xy_grad_z_i[i];


                    rhow2_r = rhoCx_rho_r[i];

                    rxw2_r = rhoCx_grad_x_r[i];

                    ryw2_r = rhoCx_grad_y_r[i];

                    rzw2_r = rhoCx_grad_z_r[i];

                    rhow2_i = rhoCx_rho_i[i];

                    rxw2_i = rhoCx_grad_x_i[i];

                    ryw2_i = rhoCx_grad_y_i[i];

                    rzw2_i = rhoCx_grad_z_i[i];
                    

                    F_y_gam_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_y_gam_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_y_gamX_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_y_gamX_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_y_gamY_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_y_gamY_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_y_gamZ_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_y_gamZ_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_y_gamXX_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_y_gamXX_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_y_gamXY_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_y_gamXY_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_y_gamXZ_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_y_gamXZ_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_y_gamYX_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_y_gamYX_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_y_gamYY_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_y_gamYY_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_y_gamYZ_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_y_gamYZ_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_y_gamZX_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_y_gamZX_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_y_gamZY_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_y_gamZY_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_y_gamZZ_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_y_gamZZ_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    // Sig_yy rho_y

                    rhow1_r = Rho_sig_yy_rho_r[i];

                    rxw1_r = Rho_sig_yy_grad_x_r[i];

                    ryw1_r = Rho_sig_yy_grad_y_r[i];

                    rzw1_r = Rho_sig_yy_grad_z_r[i];

                    rhow1_i = Rho_sig_yy_rho_i[i];

                    rxw1_i = Rho_sig_yy_grad_x_i[i];

                    ryw1_i = Rho_sig_yy_grad_y_i[i];

                    rzw1_i = Rho_sig_yy_grad_z_i[i];
                

                    rhow2_r = rhoCy_rho_r[i];

                    rxw2_r = rhoCy_grad_x_r[i];

                    ryw2_r = rhoCy_grad_y_r[i];

                    rzw2_r = rhoCy_grad_z_r[i];

                    rhow2_i = rhoCy_rho_i[i];

                    rxw2_i = rhoCy_grad_x_i[i];

                    ryw2_i = rhoCy_grad_y_i[i];

                    rzw2_i = rhoCy_grad_z_i[i];
                        

                    F_y_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_y_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_y_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_y_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_y_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_y_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_y_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_y_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_y_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_y_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_y_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_y_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_y_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_y_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_y_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_y_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_y_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_y_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_y_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_y_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_y_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_y_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_y_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_y_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_y_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_y_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    // Sig_yz rho_z

                    rhow1_r = Rho_sig_yz_rho_r[i];

                    rxw1_r = Rho_sig_yz_grad_x_r[i];

                    ryw1_r = Rho_sig_yz_grad_y_r[i];

                    rzw1_r = Rho_sig_yz_grad_z_r[i];

                    rhow1_i = Rho_sig_yz_rho_i[i];

                    rxw1_i = Rho_sig_yz_grad_x_i[i];

                    ryw1_i = Rho_sig_yz_grad_y_i[i];

                    rzw1_i = Rho_sig_yz_grad_z_i[i];
                

                    rhow2_r = rhoCz_rho_r[i];

                    rxw2_r = rhoCz_grad_x_r[i];

                    ryw2_r = rhoCz_grad_y_r[i];

                    rzw2_r = rhoCz_grad_z_r[i];

                    rhow2_i = rhoCz_rho_i[i];

                    rxw2_i = rhoCz_grad_x_i[i];

                    ryw2_i = rhoCz_grad_y_i[i];

                    rzw2_i = rhoCz_grad_z_i[i];
                        

                    F_y_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_y_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_y_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_y_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_y_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_y_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_y_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_y_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_y_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_y_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_y_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_y_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_y_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_y_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_y_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_y_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_y_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_y_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_y_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_y_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_y_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_y_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_y_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_y_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_y_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_y_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;

                  

                    // Z case

                    
                    // Sig_zx rho_x

                    rhow1_r = Rho_sig_xz_rho_r[i];

                    rxw1_r = Rho_sig_xz_grad_x_r[i];

                    ryw1_r = Rho_sig_xz_grad_y_r[i];

                    rzw1_r = Rho_sig_xz_grad_z_r[i];

                    rhow1_i = Rho_sig_xz_rho_i[i];

                    rxw1_i = Rho_sig_xz_grad_x_i[i];

                    ryw1_i = Rho_sig_xz_grad_y_i[i];

                    rzw1_i = Rho_sig_xz_grad_z_i[i];


                    rhow2_r = rhoCx_rho_r[i];

                    rxw2_r = rhoCx_grad_x_r[i];

                    ryw2_r = rhoCx_grad_y_r[i];

                    rzw2_r = rhoCx_grad_z_r[i];

                    rhow2_i = rhoCx_rho_i[i];

                    rxw2_i = rhoCx_grad_x_i[i];

                    ryw2_i = rhoCx_grad_y_i[i];

                    rzw2_i = rhoCx_grad_z_i[i];
                

                    F_z_gam_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_z_gam_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_z_gamX_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_z_gamX_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_z_gamY_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_z_gamY_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_z_gamZ_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_z_gamZ_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_z_gamXX_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_z_gamXX_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_z_gamXY_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_z_gamXY_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_z_gamXZ_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_z_gamXZ_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_z_gamYX_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_z_gamYX_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_z_gamYY_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_z_gamYY_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_z_gamYZ_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_z_gamYZ_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_z_gamZX_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_z_gamZX_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_z_gamZY_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_z_gamZY_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_z_gamZZ_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_z_gamZZ_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    // Sig_zy rho_y

                    rhow1_r = Rho_sig_yz_rho_r[i];

                    rxw1_r = Rho_sig_yz_grad_x_r[i];

                    ryw1_r = Rho_sig_yz_grad_y_r[i];

                    rzw1_r = Rho_sig_yz_grad_z_r[i];

                    rhow1_i = Rho_sig_yz_rho_i[i];

                    rxw1_i = Rho_sig_yz_grad_x_i[i];

                    ryw1_i = Rho_sig_yz_grad_y_i[i];

                    rzw1_i = Rho_sig_yz_grad_z_i[i];
                

                    rhow2_r = rhoCy_rho_r[i];

                    rxw2_r = rhoCy_grad_x_r[i];

                    ryw2_r = rhoCy_grad_y_r[i];

                    rzw2_r = rhoCy_grad_z_r[i];

                    rhow2_i = rhoCy_rho_i[i];

                    rxw2_i = rhoCy_grad_x_i[i];

                    ryw2_i = rhoCy_grad_y_i[i];

                    rzw2_i = rhoCy_grad_z_i[i];
                        

                    F_z_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_z_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_z_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_z_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_z_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_z_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_z_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_z_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_z_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_z_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_z_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_z_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_z_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_z_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_z_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_z_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_z_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_z_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_z_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_z_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_z_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_z_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_z_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_z_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_z_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_z_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;



                    // Sig_zx rho_z

                    rhow1_r = Rho_sig_zz_rho_r[i];

                    rxw1_r = Rho_sig_zz_grad_x_r[i];

                    ryw1_r = Rho_sig_zz_grad_y_r[i];

                    rzw1_r = Rho_sig_zz_grad_z_r[i];

                    rhow1_i = Rho_sig_zz_rho_i[i];

                    rxw1_i = Rho_sig_zz_grad_x_i[i];

                    ryw1_i = Rho_sig_zz_grad_y_i[i];

                    rzw1_i = Rho_sig_zz_grad_z_i[i];
                

                    rhow2_r = rhoCz_rho_r[i];

                    rxw2_r = rhoCz_grad_x_r[i];

                    ryw2_r = rhoCz_grad_y_r[i];

                    rzw2_r = rhoCz_grad_z_r[i];

                    rhow2_i = rhoCz_rho_i[i];

                    rxw2_i = rhoCz_grad_x_i[i];

                    ryw2_i = rhoCz_grad_y_i[i];

                    rzw2_i = rhoCz_grad_z_i[i];
                        

                    F_z_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_z_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_z_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_z_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_z_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_z_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_z_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_z_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_z_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_z_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_z_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_z_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_z_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_z_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_z_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_z_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_z_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_z_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_z_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_z_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_z_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_z_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_z_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_z_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_z_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_z_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;
                    


                }
        }
        }
        if (fstr::upcase(quadMode) == "TPA_II") 
        {
            // This code is inteded to compute F_b_sigma fock matrices for the final E3 contraction for tpa calculations.

            for (int32_t j = 0; j < numdens / 6; j++)
            {
                
                auto F_x_gam_r = gam(6 * j);

                auto F_x_gam_i = gam(6 * j + 1);

                auto F_x_gamX_r = gamX(6 * j);

                auto F_x_gamX_i = gamX(6 * j + 1);

                auto F_x_gamY_r = gamY(6 * j);

                auto F_x_gamY_i = gamY(6 * j + 1);

                auto F_x_gamZ_r = gamZ(6 * j);

                auto F_x_gamZ_i = gamZ(6 * j + 1);

                auto F_x_gamXX_r = gamXX(6 * j);

                auto F_x_gamXX_i = gamXX(6 * j + 1);

                auto F_x_gamXY_r = gamXY(6 * j);

                auto F_x_gamXY_i = gamXY(6 * j + 1);

                auto F_x_gamXZ_r = gamXZ(6 * j);

                auto F_x_gamXZ_i = gamXZ(6 * j + 1);

                auto F_x_gamYX_r = gamYX(6 * j);

                auto F_x_gamYX_i = gamYX(6 * j + 1);

                auto F_x_gamYY_r = gamYY(6 * j);

                auto F_x_gamYY_i = gamYY(6 * j + 1);

                auto F_x_gamYZ_r = gamYZ(6 * j);

                auto F_x_gamYZ_i = gamYZ(6 * j + 1);

                auto F_x_gamZX_r = gamZX(6 * j);

                auto F_x_gamZX_i = gamZX(6 * j + 1);

                auto F_x_gamZY_r = gamZY(6 * j);

                auto F_x_gamZY_i = gamZY(6 * j + 1);

                auto F_x_gamZZ_r = gamZZ(6 * j);

                auto F_x_gamZZ_i = gamZZ(6 * j + 1);

                // Fy

                auto F_y_gam_r = gam(6 * j  + 2);
 
                auto F_y_gam_i = gam(6 * j  + 3);

                auto F_y_gamX_r = gamX(6 * j  + 2);

                auto F_y_gamX_i = gamX(6 * j  + 3);

                auto F_y_gamY_r = gamY(6 * j  + 2);

                auto F_y_gamY_i = gamY(6 * j  + 3);

                auto F_y_gamZ_r = gamZ(6 * j  + 2);

                auto F_y_gamZ_i = gamZ(6 * j  + 3);

                auto F_y_gamXX_r = gamXX(6 * j  + 2);

                auto F_y_gamXX_i = gamXX(6 * j  + 3);

                auto F_y_gamXY_r = gamXY(6 * j  + 2);

                auto F_y_gamXY_i = gamXY(6 * j  + 3);

                auto F_y_gamXZ_r = gamXZ(6 * j  + 2);

                auto F_y_gamXZ_i = gamXZ(6 * j  + 3);

                auto F_y_gamYX_r = gamYX(6 * j  + 2);

                auto F_y_gamYX_i = gamYX(6 * j  + 3);

                auto F_y_gamYY_r = gamYY(6 * j  + 2);

                auto F_y_gamYY_i = gamYY(6 * j  + 3);

                auto F_y_gamYZ_r = gamYZ(6 * j  + 2);

                auto F_y_gamYZ_i = gamYZ(6 * j  + 3);

                auto F_y_gamZX_r = gamZX(6 * j  + 2);

                auto F_y_gamZX_i = gamZX(6 * j  + 3);

                auto F_y_gamZY_r = gamZY(6 * j  + 2);

                auto F_y_gamZY_i = gamZY(6 * j  + 3);

                auto F_y_gamZZ_r = gamZZ(6 * j  + 2);

                auto F_y_gamZZ_i = gamZZ(6 * j  + 3);

                // Fz

               auto F_z_gam_r = gam(6 * j  + 4);

                auto F_z_gam_i = gam(6 * j  + 5);

                auto F_z_gamX_r = gamX(6 * j  + 4);

                auto F_z_gamX_i = gamX(6 * j  + 5);

                auto F_z_gamY_r = gamY(6 * j  + 4);

                auto F_z_gamY_i = gamY(6 * j  + 5);

                auto F_z_gamZ_r = gamZ(6 * j  + 4);

                auto F_z_gamZ_i = gamZ(6 * j  + 5);

                auto F_z_gamXX_r = gamXX(6 * j  + 4);

                auto F_z_gamXX_i = gamXX(6 * j  + 5);

                auto F_z_gamXY_r = gamXY(6 * j  + 4);

                auto F_z_gamXY_i = gamXY(6 * j  + 5);

                auto F_z_gamXZ_r = gamXZ(6 * j  + 4);

                auto F_z_gamXZ_i = gamXZ(6 * j  + 5);

                auto F_z_gamYX_r = gamYX(6 * j  + 4);

                auto F_z_gamYX_i = gamYX(6 * j  + 5);

                auto F_z_gamYY_r = gamYY(6 * j  + 4);

                auto F_z_gamYY_i = gamYY(6 * j  + 5);

                auto F_z_gamYZ_r = gamYZ(6 * j  + 4);

                auto F_z_gamYZ_i = gamYZ(6 * j  + 5);

                auto F_z_gamZX_r = gamZX(6 * j  + 4);

                auto F_z_gamZX_i = gamZX(6 * j  + 5);

                auto F_z_gamZY_r = gamZY(6 * j  + 4);

                auto F_z_gamZY_i = gamZY(6 * j  + 5);

                auto F_z_gamZZ_r = gamZZ(6 * j  + 4);

                auto F_z_gamZZ_i = gamZZ(6 * j  + 5);


                // Perturbed densities

                auto rhoBx_rho_r = rwDensityGrid.alphaDensity(36 * j);

                auto rhoBx_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j);

                auto rhoBx_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j);

                auto rhoBx_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j);


                auto rhoBx_rho_i = rwDensityGrid.alphaDensity(36 * j + 1);

                auto rhoBx_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 1);

                auto rhoBx_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 1);

                auto rhoBx_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 1);


                auto rhoBy_rho_r = rwDensityGrid.alphaDensity(36 * j + 2);

                auto rhoBy_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j + 2);

                auto rhoBy_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j + 2 );

                auto rhoBy_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j + 2);


                auto rhoBy_rho_i = rwDensityGrid.alphaDensity(36 * j + 3);

                auto rhoBy_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 3);

                auto rhoBy_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 3);

                auto rhoBy_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 3);


                auto rhoBz_rho_r = rwDensityGrid.alphaDensity(36 * j + 4);

                auto rhoBz_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j+ 4);

                auto rhoBz_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j+ 4);

                auto rhoBz_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j+ 4);


                auto rhoBz_rho_i = rwDensityGrid.alphaDensity(36 * j + 5);

                auto rhoBz_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 5);

                auto rhoBz_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 5);

                auto rhoBz_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 5);

                // C 

                auto rhoCx_rho_r = rwDensityGrid.alphaDensity(36 * j + 6);

                auto rhoCx_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j + 6);

                auto rhoCx_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j + 6);

                auto rhoCx_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j + 6);

                auto rhoCx_rho_i = rwDensityGrid.alphaDensity(36 * j + 7);

                auto rhoCx_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 7);

                auto rhoCx_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 7);

                auto rhoCx_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 7);


                auto rhoCy_rho_r = rwDensityGrid.alphaDensity(36 * j + 8);

                auto rhoCy_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j + 8);

                auto rhoCy_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j + 8);

                auto rhoCy_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j + 8);

                auto rhoCy_rho_i = rwDensityGrid.alphaDensity(36 * j + 9);

                auto rhoCy_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 9);

                auto rhoCy_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 9);

                auto rhoCy_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 9);


                auto rhoCz_rho_r = rwDensityGrid.alphaDensity(36 * j + 10);

                auto rhoCz_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j  + 10);

                auto rhoCz_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j  + 10);

                auto rhoCz_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j  + 10);

                auto rhoCz_rho_i = rwDensityGrid.alphaDensity(36 * j + 11);

                auto rhoCz_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 11);

                auto rhoCz_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 11);

                auto rhoCz_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 11);


                auto Rho_sig_xx_rho_r = rwDensityGrid.alphaDensity(36 * j + 12);

                auto Rho_sig_xx_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j  + 12);

                auto Rho_sig_xx_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j  + 12);

                auto Rho_sig_xx_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j  + 12);

                auto Rho_sig_xx_rho_i = rwDensityGrid.alphaDensity(36 * j + 13);

                auto Rho_sig_xx_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 13);

                auto Rho_sig_xx_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 13);

                auto Rho_sig_xx_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 13);



                auto Rho_sig_yy_rho_r = rwDensityGrid.alphaDensity(36 * j + 14);

                auto Rho_sig_yy_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j  + 14);

                auto Rho_sig_yy_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j  + 14);

                auto Rho_sig_yy_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j  + 14);

                auto Rho_sig_yy_rho_i = rwDensityGrid.alphaDensity(36 * j + 15);

                auto Rho_sig_yy_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 15);

                auto Rho_sig_yy_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 15);

                auto Rho_sig_yy_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 15);



                auto Rho_sig_zz_rho_r = rwDensityGrid.alphaDensity(36 * j + 16);

                auto Rho_sig_zz_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j  + 16);

                auto Rho_sig_zz_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j  + 16);

                auto Rho_sig_zz_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j  + 16);

                auto Rho_sig_zz_rho_i = rwDensityGrid.alphaDensity(36 * j + 17);

                auto Rho_sig_zz_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 17);

                auto Rho_sig_zz_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 17);

                auto Rho_sig_zz_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 17);


                auto Rho_sig_xy_rho_r = rwDensityGrid.alphaDensity(36 * j + 18);

                auto Rho_sig_xy_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j  + 18);

                auto Rho_sig_xy_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j  + 18);

                auto Rho_sig_xy_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j  + 18);

                auto Rho_sig_xy_rho_i = rwDensityGrid.alphaDensity(36 * j + 19);

                auto Rho_sig_xy_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 19);

                auto Rho_sig_xy_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 19);

                auto Rho_sig_xy_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 19);


                auto Rho_sig_xz_rho_r = rwDensityGrid.alphaDensity(36 * j + 20);

                auto Rho_sig_xz_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j  + 20);

                auto Rho_sig_xz_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j  + 20);

                auto Rho_sig_xz_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j  + 20);

                auto Rho_sig_xz_rho_i = rwDensityGrid.alphaDensity(36 * j + 21);

                auto Rho_sig_xz_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 21);

                auto Rho_sig_xz_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 21);

                auto Rho_sig_xz_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 21);



                auto Rho_sig_yz_rho_r = rwDensityGrid.alphaDensity(36 * j + 22);

                auto Rho_sig_yz_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j  + 22);

                auto Rho_sig_yz_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j  + 22);

                auto Rho_sig_yz_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j  + 22);

                auto Rho_sig_yz_rho_i = rwDensityGrid.alphaDensity(36 * j + 23);

                auto Rho_sig_yz_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 23);

                auto Rho_sig_yz_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 23);

                auto Rho_sig_yz_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 23);


                auto Rho_lamtau_xx_rho_r = rwDensityGrid.alphaDensity(36 * j + 24);

                auto Rho_lamtau_xx_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j  + 24);

                auto Rho_lamtau_xx_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j  + 24);

                auto Rho_lamtau_xx_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j  + 24);

                auto Rho_lamtau_xx_rho_i = rwDensityGrid.alphaDensity(36 * j + 25);

                auto Rho_lamtau_xx_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 25);

                auto Rho_lamtau_xx_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 25);

                auto Rho_lamtau_xx_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 25);


                auto Rho_lamtau_yy_rho_r = rwDensityGrid.alphaDensity(36 * j + 26);

                auto Rho_lamtau_yy_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j  + 26);

                auto Rho_lamtau_yy_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j  + 26);

                auto Rho_lamtau_yy_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j  + 26);

                auto Rho_lamtau_yy_rho_i = rwDensityGrid.alphaDensity(36 * j + 27);

                auto Rho_lamtau_yy_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 27);

                auto Rho_lamtau_yy_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 27);

                auto Rho_lamtau_yy_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 27);


                auto Rho_lamtau_zz_rho_r = rwDensityGrid.alphaDensity(36 * j + 28);

                auto Rho_lamtau_zz_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j  + 28);

                auto Rho_lamtau_zz_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j  + 28);

                auto Rho_lamtau_zz_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j  + 28);

                auto Rho_lamtau_zz_rho_i = rwDensityGrid.alphaDensity(36 * j + 29);

                auto Rho_lamtau_zz_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 29);

                auto Rho_lamtau_zz_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 29);

                auto Rho_lamtau_zz_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 29);


                auto Rho_lamtau_xy_rho_r = rwDensityGrid.alphaDensity(36 * j + 30);

                auto Rho_lamtau_xy_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j  + 30);

                auto Rho_lamtau_xy_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j  + 30);

                auto Rho_lamtau_xy_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j  + 30);

                auto Rho_lamtau_xy_rho_i = rwDensityGrid.alphaDensity(36 * j + 31);

                auto Rho_lamtau_xy_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 31);

                auto Rho_lamtau_xy_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 31);

                auto Rho_lamtau_xy_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 31);


                auto Rho_lamtau_xz_rho_r = rwDensityGrid.alphaDensity(36 * j + 32);

                auto Rho_lamtau_xz_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j  + 32);

                auto Rho_lamtau_xz_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j  + 32);

                auto Rho_lamtau_xz_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j  + 32);

                auto Rho_lamtau_xz_rho_i = rwDensityGrid.alphaDensity(36 * j + 33);

                auto Rho_lamtau_xz_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 33);

                auto Rho_lamtau_xz_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 33);

                auto Rho_lamtau_xz_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 33);
                

                auto Rho_lamtau_yz_rho_r = rwDensityGrid.alphaDensity(36 * j + 34);

                auto Rho_lamtau_yz_grad_x_r = rwDensityGrid.alphaDensityGradientX(36 * j  + 34);

                auto Rho_lamtau_yz_grad_y_r = rwDensityGrid.alphaDensityGradientY(36 * j  + 34);

                auto Rho_lamtau_yz_grad_z_r = rwDensityGrid.alphaDensityGradientZ(36 * j  + 34);

                auto Rho_lamtau_yz_rho_i = rwDensityGrid.alphaDensity(36 * j + 35);

                auto Rho_lamtau_yz_grad_x_i = rwDensityGrid.alphaDensityGradientX(36 * j + 35);

                auto Rho_lamtau_yz_grad_y_i = rwDensityGrid.alphaDensityGradientY(36 * j + 35);

                auto Rho_lamtau_yz_grad_z_i = rwDensityGrid.alphaDensityGradientZ(36 * j + 35);


                for (int32_t i = 0; i < npoints; i++)
                {
                    

                    // Sig_xx rho_x

                    double rhow1_r = Rho_sig_xx_rho_r[i];

                    double rxw1_r = Rho_sig_xx_grad_x_r[i];

                    double ryw1_r = Rho_sig_xx_grad_y_r[i];

                    double rzw1_r = Rho_sig_xx_grad_z_r[i];

                    double rhow1_i = Rho_sig_xx_rho_i[i];

                    double rxw1_i = Rho_sig_xx_grad_x_i[i];

                    double ryw1_i = Rho_sig_xx_grad_y_i[i];

                    double rzw1_i = Rho_sig_xx_grad_z_i[i];


                    double rhow2_r = rhoCx_rho_r[i];

                    double rxw2_r = rhoCx_grad_x_r[i];

                    double ryw2_r = rhoCx_grad_y_r[i];

                    double rzw2_r = rhoCx_grad_z_r[i];

                    double rhow2_i = rhoCx_rho_i[i];

                    double rxw2_i = rhoCx_grad_x_i[i];

                    double ryw2_i = rhoCx_grad_y_i[i];

                    double rzw2_i = rhoCx_grad_z_i[i];
                    

                    F_x_gam_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_x_gam_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_x_gamX_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_x_gamX_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_x_gamY_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_x_gamY_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_x_gamZ_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_x_gamZ_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_x_gamXX_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_x_gamXX_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_x_gamXY_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_x_gamXY_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_x_gamXZ_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_x_gamXZ_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_x_gamYX_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_x_gamYX_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_x_gamYY_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_x_gamYY_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_x_gamYZ_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_x_gamYZ_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_x_gamZX_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_x_gamZX_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_x_gamZY_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_x_gamZY_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_x_gamZZ_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_x_gamZZ_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    // Sig_xy rho_y

                    rhow1_r = Rho_sig_xy_rho_r[i];

                    rxw1_r = Rho_sig_xy_grad_x_r[i];

                    ryw1_r = Rho_sig_xy_grad_y_r[i];

                    rzw1_r = Rho_sig_xy_grad_z_r[i];

                    rhow1_i = Rho_sig_xy_rho_i[i];

                    rxw1_i = Rho_sig_xy_grad_x_i[i];

                    ryw1_i = Rho_sig_xy_grad_y_i[i];

                    rzw1_i = Rho_sig_xy_grad_z_i[i];
                

                    rhow2_r = rhoCy_rho_r[i];

                    rxw2_r = rhoCy_grad_x_r[i];

                    ryw2_r = rhoCy_grad_y_r[i];

                    rzw2_r = rhoCy_grad_z_r[i];

                    rhow2_i = rhoCy_rho_i[i];

                    rxw2_i = rhoCy_grad_x_i[i];

                    ryw2_i = rhoCy_grad_y_i[i];

                    rzw2_i = rhoCy_grad_z_i[i];
                        

                    F_x_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_x_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_x_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_x_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_x_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_x_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_x_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_x_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_x_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_x_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_x_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_x_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_x_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_x_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_x_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_x_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_x_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_x_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_x_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_x_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_x_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_x_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_x_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_x_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_x_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_x_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;



                    // Sig_xz rho_z

                    rhow1_r = Rho_sig_xz_rho_r[i];

                    rxw1_r = Rho_sig_xz_grad_x_r[i];

                    ryw1_r = Rho_sig_xz_grad_y_r[i];

                    rzw1_r = Rho_sig_xz_grad_z_r[i];

                    rhow1_i = Rho_sig_xz_rho_i[i];

                    rxw1_i = Rho_sig_xz_grad_x_i[i];

                    ryw1_i = Rho_sig_xz_grad_y_i[i];

                    rzw1_i = Rho_sig_xz_grad_z_i[i];
                

                    rhow2_r = rhoCz_rho_r[i];

                    rxw2_r = rhoCz_grad_x_r[i];

                    ryw2_r = rhoCz_grad_y_r[i];

                    rzw2_r = rhoCz_grad_z_r[i];

                    rhow2_i = rhoCz_rho_i[i];

                    rxw2_i = rhoCz_grad_x_i[i];

                    ryw2_i = rhoCz_grad_y_i[i];

                    rzw2_i = rhoCz_grad_z_i[i];
                        

                    F_x_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_x_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_x_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_x_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_x_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_x_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_x_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_x_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_x_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_x_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_x_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_x_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_x_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_x_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_x_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_x_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_x_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_x_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_x_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_x_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_x_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_x_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_x_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_x_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_x_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_x_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;



                    // lamtau_xx rho_x

                    rhow1_r = Rho_lamtau_xx_rho_r[i];

                    rxw1_r = Rho_lamtau_xx_grad_x_r[i];

                    ryw1_r = Rho_lamtau_xx_grad_y_r[i];

                    rzw1_r = Rho_lamtau_xx_grad_z_r[i];

                    rhow1_i = Rho_lamtau_xx_rho_i[i];

                    rxw1_i = Rho_lamtau_xx_grad_x_i[i];

                    ryw1_i = Rho_lamtau_xx_grad_y_i[i];

                    rzw1_i = Rho_lamtau_xx_grad_z_i[i];


                    rhow2_r = rhoBx_rho_r[i];

                    rxw2_r = rhoBx_grad_x_r[i];

                    ryw2_r = rhoBx_grad_y_r[i];

                    rzw2_r = rhoBx_grad_z_r[i];

                    rhow2_i = rhoBx_rho_i[i];

                    rxw2_i = rhoBx_grad_x_i[i];

                    ryw2_i = rhoBx_grad_y_i[i];

                    rzw2_i = rhoBx_grad_z_i[i];
                    

                    F_x_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_x_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_x_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_x_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_x_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_x_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_x_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_x_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_x_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_x_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_x_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_x_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_x_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_x_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_x_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_x_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_x_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_x_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_x_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_x_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_x_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_x_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_x_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_x_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_x_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_x_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    // lamtau_xy rho_y

                    rhow1_r = Rho_lamtau_xy_rho_r[i];

                    rxw1_r = Rho_lamtau_xy_grad_x_r[i];

                    ryw1_r = Rho_lamtau_xy_grad_y_r[i];

                    rzw1_r = Rho_lamtau_xy_grad_z_r[i];

                    rhow1_i = Rho_lamtau_xy_rho_i[i];

                    rxw1_i = Rho_lamtau_xy_grad_x_i[i];

                    ryw1_i = Rho_lamtau_xy_grad_y_i[i];

                    rzw1_i = Rho_lamtau_xy_grad_z_i[i];
                

                    rhow2_r = rhoBy_rho_r[i];

                    rxw2_r = rhoBy_grad_x_r[i];

                    ryw2_r = rhoBy_grad_y_r[i];

                    rzw2_r = rhoBy_grad_z_r[i];

                    rhow2_i = rhoBy_rho_i[i];

                    rxw2_i = rhoBy_grad_x_i[i];

                    ryw2_i = rhoBy_grad_y_i[i];

                    rzw2_i = rhoBy_grad_z_i[i];
                        

                    F_x_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_x_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_x_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_x_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_x_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_x_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_x_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_x_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_x_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_x_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_x_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_x_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_x_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_x_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_x_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_x_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_x_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_x_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_x_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_x_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_x_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_x_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_x_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_x_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_x_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_x_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;



                    // lamtau_xz rho_z

                    rhow1_r = Rho_lamtau_xz_rho_r[i];

                    rxw1_r = Rho_lamtau_xz_grad_x_r[i];

                    ryw1_r = Rho_lamtau_xz_grad_y_r[i];

                    rzw1_r = Rho_lamtau_xz_grad_z_r[i];

                    rhow1_i = Rho_lamtau_xz_rho_i[i];

                    rxw1_i = Rho_lamtau_xz_grad_x_i[i];

                    ryw1_i = Rho_lamtau_xz_grad_y_i[i];

                    rzw1_i = Rho_lamtau_xz_grad_z_i[i];
                

                    rhow2_r = rhoBz_rho_r[i];

                    rxw2_r = rhoBz_grad_x_r[i];

                    ryw2_r = rhoBz_grad_y_r[i];

                    rzw2_r = rhoBz_grad_z_r[i];

                    rhow2_i = rhoBz_rho_i[i];

                    rxw2_i = rhoBz_grad_x_i[i];

                    ryw2_i = rhoBz_grad_y_i[i];

                    rzw2_i = rhoBz_grad_z_i[i];
                        

                    F_x_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_x_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_x_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_x_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_x_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_x_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_x_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_x_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_x_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_x_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_x_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_x_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_x_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_x_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_x_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_x_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_x_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_x_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_x_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_x_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_x_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_x_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_x_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_x_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_x_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_x_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;




                    // Y case

                    
                    // Sig_yx rho_x

                    rhow1_r = Rho_sig_xy_rho_r[i];

                    rxw1_r = Rho_sig_xy_grad_x_r[i];

                    ryw1_r = Rho_sig_xy_grad_y_r[i];

                    rzw1_r = Rho_sig_xy_grad_z_r[i];

                    rhow1_i = Rho_sig_xy_rho_i[i];

                    rxw1_i = Rho_sig_xy_grad_x_i[i];

                    ryw1_i = Rho_sig_xy_grad_y_i[i];

                    rzw1_i = Rho_sig_xy_grad_z_i[i];


                    rhow2_r = rhoCx_rho_r[i];

                    rxw2_r = rhoCx_grad_x_r[i];

                    ryw2_r = rhoCx_grad_y_r[i];

                    rzw2_r = rhoCx_grad_z_r[i];

                    rhow2_i = rhoCx_rho_i[i];

                    rxw2_i = rhoCx_grad_x_i[i];

                    ryw2_i = rhoCx_grad_y_i[i];

                    rzw2_i = rhoCx_grad_z_i[i];
                    

                    F_y_gam_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_y_gam_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_y_gamX_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_y_gamX_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_y_gamY_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_y_gamY_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_y_gamZ_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_y_gamZ_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_y_gamXX_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_y_gamXX_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_y_gamXY_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_y_gamXY_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_y_gamXZ_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_y_gamXZ_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_y_gamYX_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_y_gamYX_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_y_gamYY_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_y_gamYY_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_y_gamYZ_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_y_gamYZ_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_y_gamZX_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_y_gamZX_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_y_gamZY_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_y_gamZY_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_y_gamZZ_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_y_gamZZ_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    // Sig_yy rho_y

                    rhow1_r = Rho_sig_yy_rho_r[i];

                    rxw1_r = Rho_sig_yy_grad_x_r[i];

                    ryw1_r = Rho_sig_yy_grad_y_r[i];

                    rzw1_r = Rho_sig_yy_grad_z_r[i];

                    rhow1_i = Rho_sig_yy_rho_i[i];

                    rxw1_i = Rho_sig_yy_grad_x_i[i];

                    ryw1_i = Rho_sig_yy_grad_y_i[i];

                    rzw1_i = Rho_sig_yy_grad_z_i[i];
                

                    rhow2_r = rhoCy_rho_r[i];

                    rxw2_r = rhoCy_grad_x_r[i];

                    ryw2_r = rhoCy_grad_y_r[i];

                    rzw2_r = rhoCy_grad_z_r[i];

                    rhow2_i = rhoCy_rho_i[i];

                    rxw2_i = rhoCy_grad_x_i[i];

                    ryw2_i = rhoCy_grad_y_i[i];

                    rzw2_i = rhoCy_grad_z_i[i];
                        

                    F_y_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_y_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_y_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_y_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_y_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_y_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_y_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_y_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_y_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_y_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_y_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_y_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_y_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_y_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_y_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_y_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_y_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_y_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_y_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_y_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_y_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_y_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_y_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_y_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_y_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_y_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    // Sig_yz rho_z

                    rhow1_r = Rho_sig_yz_rho_r[i];

                    rxw1_r = Rho_sig_yz_grad_x_r[i];

                    ryw1_r = Rho_sig_yz_grad_y_r[i];

                    rzw1_r = Rho_sig_yz_grad_z_r[i];

                    rhow1_i = Rho_sig_yz_rho_i[i];

                    rxw1_i = Rho_sig_yz_grad_x_i[i];

                    ryw1_i = Rho_sig_yz_grad_y_i[i];

                    rzw1_i = Rho_sig_yz_grad_z_i[i];
                

                    rhow2_r = rhoCz_rho_r[i];

                    rxw2_r = rhoCz_grad_x_r[i];

                    ryw2_r = rhoCz_grad_y_r[i];

                    rzw2_r = rhoCz_grad_z_r[i];

                    rhow2_i = rhoCz_rho_i[i];

                    rxw2_i = rhoCz_grad_x_i[i];

                    ryw2_i = rhoCz_grad_y_i[i];

                    rzw2_i = rhoCz_grad_z_i[i];
                        

                    F_y_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_y_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_y_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_y_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_y_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_y_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_y_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_y_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_y_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_y_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_y_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_y_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_y_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_y_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_y_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_y_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_y_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_y_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_y_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_y_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_y_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_y_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_y_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_y_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_y_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_y_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;

                    // lamtau_yx rho_x

                    rhow1_r = Rho_lamtau_xy_rho_r[i];

                    rxw1_r = Rho_lamtau_xy_grad_x_r[i];

                    ryw1_r = Rho_lamtau_xy_grad_y_r[i];

                    rzw1_r = Rho_lamtau_xy_grad_z_r[i];

                    rhow1_i = Rho_lamtau_xy_rho_i[i];

                    rxw1_i = Rho_lamtau_xy_grad_x_i[i];

                    ryw1_i = Rho_lamtau_xy_grad_y_i[i];

                    rzw1_i = Rho_lamtau_xy_grad_z_i[i];


                    rhow2_r = rhoBx_rho_r[i];

                    rxw2_r = rhoBx_grad_x_r[i];

                    ryw2_r = rhoBx_grad_y_r[i];

                    rzw2_r = rhoBx_grad_z_r[i];

                    rhow2_i = rhoBx_rho_i[i];

                    rxw2_i = rhoBx_grad_x_i[i];

                    ryw2_i = rhoBx_grad_y_i[i];

                    rzw2_i = rhoBx_grad_z_i[i];
                    

                    F_y_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_y_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_y_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_y_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_y_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_y_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_y_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_y_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_y_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_y_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_y_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_y_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_y_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_y_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_y_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_y_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_y_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_y_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_y_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_y_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_y_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_y_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_y_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_y_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_y_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_y_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    // lamtau_yy rho_y

                    rhow1_r = Rho_lamtau_yy_rho_r[i];

                    rxw1_r = Rho_lamtau_yy_grad_x_r[i];

                    ryw1_r = Rho_lamtau_yy_grad_y_r[i];

                    rzw1_r = Rho_lamtau_yy_grad_z_r[i];

                    rhow1_i = Rho_lamtau_yy_rho_i[i];

                    rxw1_i = Rho_lamtau_yy_grad_x_i[i];

                    ryw1_i = Rho_lamtau_yy_grad_y_i[i];

                    rzw1_i = Rho_lamtau_yy_grad_z_i[i];
                

                    rhow2_r = rhoBy_rho_r[i];

                    rxw2_r = rhoBy_grad_x_r[i];

                    ryw2_r = rhoBy_grad_y_r[i];

                    rzw2_r = rhoBy_grad_z_r[i];

                    rhow2_i = rhoBy_rho_i[i];

                    rxw2_i = rhoBy_grad_x_i[i];

                    ryw2_i = rhoBy_grad_y_i[i];

                    rzw2_i = rhoBy_grad_z_i[i];
                        

                    F_y_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_y_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_y_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_y_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_y_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_y_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_y_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_y_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_y_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_y_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_y_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_y_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_y_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_y_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_y_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_y_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_y_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_y_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_y_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_y_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_y_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_y_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_y_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_y_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_y_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_y_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;



                    // lamtau_yz rho_z

                    rhow1_r = Rho_lamtau_yz_rho_r[i];

                    rxw1_r = Rho_lamtau_yz_grad_x_r[i];

                    ryw1_r = Rho_lamtau_yz_grad_y_r[i];

                    rzw1_r = Rho_lamtau_yz_grad_z_r[i];

                    rhow1_i = Rho_lamtau_yz_rho_i[i];

                    rxw1_i = Rho_lamtau_yz_grad_x_i[i];

                    ryw1_i = Rho_lamtau_yz_grad_y_i[i];

                    rzw1_i = Rho_lamtau_yz_grad_z_i[i];
                

                    rhow2_r = rhoBz_rho_r[i];

                    rxw2_r = rhoBz_grad_x_r[i];

                    ryw2_r = rhoBz_grad_y_r[i];

                    rzw2_r = rhoBz_grad_z_r[i];

                    rhow2_i = rhoBz_rho_i[i];

                    rxw2_i = rhoBz_grad_x_i[i];

                    ryw2_i = rhoBz_grad_y_i[i];

                    rzw2_i = rhoBz_grad_z_i[i];
                        

                    F_y_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_y_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_y_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_y_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_y_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_y_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_y_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_y_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_y_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_y_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_y_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_y_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_y_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_y_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_y_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_y_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_y_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_y_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_y_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_y_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_y_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_y_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_y_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_y_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_y_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_y_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    // Z case


                    
                    // Sig_zx rho_x

                    rhow1_r = Rho_sig_xz_rho_r[i];

                    rxw1_r = Rho_sig_xz_grad_x_r[i];

                    ryw1_r = Rho_sig_xz_grad_y_r[i];

                    rzw1_r = Rho_sig_xz_grad_z_r[i];

                    rhow1_i = Rho_sig_xz_rho_i[i];

                    rxw1_i = Rho_sig_xz_grad_x_i[i];

                    ryw1_i = Rho_sig_xz_grad_y_i[i];

                    rzw1_i = Rho_sig_xz_grad_z_i[i];


                    rhow2_r = rhoCx_rho_r[i];

                    rxw2_r = rhoCx_grad_x_r[i];

                    ryw2_r = rhoCx_grad_y_r[i];

                    rzw2_r = rhoCx_grad_z_r[i];

                    rhow2_i = rhoCx_rho_i[i];

                    rxw2_i = rhoCx_grad_x_i[i];

                    ryw2_i = rhoCx_grad_y_i[i];

                    rzw2_i = rhoCx_grad_z_i[i];
                

                    F_z_gam_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_z_gam_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_z_gamX_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_z_gamX_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_z_gamY_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_z_gamY_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_z_gamZ_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_z_gamZ_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_z_gamXX_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_z_gamXX_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_z_gamXY_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_z_gamXY_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_z_gamXZ_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_z_gamXZ_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_z_gamYX_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_z_gamYX_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_z_gamYY_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_z_gamYY_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_z_gamYZ_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_z_gamYZ_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_z_gamZX_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_z_gamZX_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_z_gamZY_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_z_gamZY_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_z_gamZZ_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_z_gamZZ_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    // Sig_zy rho_y

                    rhow1_r = Rho_sig_yz_rho_r[i];

                    rxw1_r = Rho_sig_yz_grad_x_r[i];

                    ryw1_r = Rho_sig_yz_grad_y_r[i];

                    rzw1_r = Rho_sig_yz_grad_z_r[i];

                    rhow1_i = Rho_sig_yz_rho_i[i];

                    rxw1_i = Rho_sig_yz_grad_x_i[i];

                    ryw1_i = Rho_sig_yz_grad_y_i[i];

                    rzw1_i = Rho_sig_yz_grad_z_i[i];
                

                    rhow2_r = rhoCy_rho_r[i];

                    rxw2_r = rhoCy_grad_x_r[i];

                    ryw2_r = rhoCy_grad_y_r[i];

                    rzw2_r = rhoCy_grad_z_r[i];

                    rhow2_i = rhoCy_rho_i[i];

                    rxw2_i = rhoCy_grad_x_i[i];

                    ryw2_i = rhoCy_grad_y_i[i];

                    rzw2_i = rhoCy_grad_z_i[i];
                        

                    F_z_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_z_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_z_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_z_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_z_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_z_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_z_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_z_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_z_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_z_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_z_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_z_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_z_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_z_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_z_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_z_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_z_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_z_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_z_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_z_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_z_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_z_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_z_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_z_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_z_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_z_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;



                    // Sig_zx rho_z

                    rhow1_r = Rho_sig_zz_rho_r[i];

                    rxw1_r = Rho_sig_zz_grad_x_r[i];

                    ryw1_r = Rho_sig_zz_grad_y_r[i];

                    rzw1_r = Rho_sig_zz_grad_z_r[i];

                    rhow1_i = Rho_sig_zz_rho_i[i];

                    rxw1_i = Rho_sig_zz_grad_x_i[i];

                    ryw1_i = Rho_sig_zz_grad_y_i[i];

                    rzw1_i = Rho_sig_zz_grad_z_i[i];
                

                    rhow2_r = rhoCz_rho_r[i];

                    rxw2_r = rhoCz_grad_x_r[i];

                    ryw2_r = rhoCz_grad_y_r[i];

                    rzw2_r = rhoCz_grad_z_r[i];

                    rhow2_i = rhoCz_rho_i[i];

                    rxw2_i = rhoCz_grad_x_i[i];

                    ryw2_i = rhoCz_grad_y_i[i];

                    rzw2_i = rhoCz_grad_z_i[i];
                        

                    F_z_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_z_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_z_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_z_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_z_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_z_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_z_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_z_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_z_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_z_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_z_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_z_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_z_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_z_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_z_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_z_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_z_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_z_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_z_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_z_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_z_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_z_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_z_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_z_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_z_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_z_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;



                    // lamtau_zx rho_x

                    rhow1_r = Rho_lamtau_xz_rho_r[i];

                    rxw1_r = Rho_lamtau_xz_grad_x_r[i];

                    ryw1_r = Rho_lamtau_xz_grad_y_r[i];

                    rzw1_r = Rho_lamtau_xz_grad_z_r[i];

                    rhow1_i = Rho_lamtau_xz_rho_i[i];

                    rxw1_i = Rho_lamtau_xz_grad_x_i[i];

                    ryw1_i = Rho_lamtau_xz_grad_y_i[i];

                    rzw1_i = Rho_lamtau_xz_grad_z_i[i];


                    rhow2_r = rhoBx_rho_r[i];

                    rxw2_r = rhoBx_grad_x_r[i];

                    ryw2_r = rhoBx_grad_y_r[i];

                    rzw2_r = rhoBx_grad_z_r[i];

                    rhow2_i = rhoBx_rho_i[i];

                    rxw2_i = rhoBx_grad_x_i[i];

                    ryw2_i = rhoBx_grad_y_i[i];

                    rzw2_i = rhoBx_grad_z_i[i];
                    

                    F_z_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_z_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_z_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_z_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_z_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_z_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_z_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_z_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_z_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_z_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_z_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_z_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_z_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_z_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_z_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_z_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_z_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_z_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_z_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_z_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_z_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_z_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_z_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_z_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_z_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_z_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    // lamtau_yz rho_y

                    rhow1_r = Rho_lamtau_yz_rho_r[i];

                    rxw1_r = Rho_lamtau_yz_grad_x_r[i];

                    ryw1_r = Rho_lamtau_yz_grad_y_r[i];

                    rzw1_r = Rho_lamtau_yz_grad_z_r[i];

                    rhow1_i = Rho_lamtau_yz_rho_i[i];

                    rxw1_i = Rho_lamtau_yz_grad_x_i[i];

                    ryw1_i = Rho_lamtau_yz_grad_y_i[i];

                    rzw1_i = Rho_lamtau_yz_grad_z_i[i];
                

                    rhow2_r = rhoBy_rho_r[i];

                    rxw2_r = rhoBy_grad_x_r[i];

                    ryw2_r = rhoBy_grad_y_r[i];

                    rzw2_r = rhoBy_grad_z_r[i];

                    rhow2_i = rhoBy_rho_i[i];

                    rxw2_i = rhoBy_grad_x_i[i];

                    ryw2_i = rhoBy_grad_y_i[i];

                    rzw2_i = rhoBy_grad_z_i[i];
                        

                    F_z_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_z_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_z_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_z_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_z_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_z_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_z_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_z_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_z_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_z_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_z_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_z_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_z_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_z_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_z_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_z_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_z_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_z_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_z_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_z_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_z_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_z_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_z_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_z_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_z_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_z_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;



                    // lamtau_zz rho_z

                    rhow1_r = Rho_lamtau_zz_rho_r[i];

                    rxw1_r = Rho_lamtau_zz_grad_x_r[i];

                    ryw1_r = Rho_lamtau_zz_grad_y_r[i];

                    rzw1_r = Rho_lamtau_zz_grad_z_r[i];

                    rhow1_i = Rho_lamtau_zz_rho_i[i];

                    rxw1_i = Rho_lamtau_zz_grad_x_i[i];

                    ryw1_i = Rho_lamtau_zz_grad_y_i[i];

                    rzw1_i = Rho_lamtau_zz_grad_z_i[i];
                

                    rhow2_r = rhoBz_rho_r[i];

                    rxw2_r = rhoBz_grad_x_r[i];

                    ryw2_r = rhoBz_grad_y_r[i];

                    rzw2_r = rhoBz_grad_z_r[i];

                    rhow2_i = rhoBz_rho_i[i];

                    rxw2_i = rhoBz_grad_x_i[i];

                    ryw2_i = rhoBz_grad_y_i[i];

                    rzw2_i = rhoBz_grad_z_i[i];
                        

                    F_z_gam_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_z_gam_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_z_gamX_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_z_gamX_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_z_gamY_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_z_gamY_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_z_gamZ_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_z_gamZ_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_z_gamXX_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_z_gamXX_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_z_gamXY_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_z_gamXY_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_z_gamXZ_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_z_gamXZ_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_z_gamYX_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_z_gamYX_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_z_gamYY_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_z_gamYY_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_z_gamYZ_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_z_gamYZ_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_z_gamZX_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_z_gamZX_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_z_gamZY_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_z_gamZY_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_z_gamZZ_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_z_gamZZ_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                }
            }
    }    
        if (fstr::upcase(quadMode) == "CRF_II")
        {
            for (int32_t j = 0; j < numdens / 2; j++)
            {
                // Perturbed densities

                // Perturbed density products to be stored

                auto gam_r = gam(2 * j);
                auto gam_i = gam(2 * j + 1);

                auto gam_X_r = gamX(2 * j);
                auto gam_X_i = gamX(2 * j + 1);
                auto gam_Y_r = gamY(2 * j);
                auto gam_Y_i = gamY(2 * j + 1);
                auto gam_Z_r = gamZ(2 * j);
                auto gam_Z_i = gamZ(2 * j + 1);

                auto gam_XX_r = gamXX(2 * j);
                auto gam_XX_i = gamXX(2 * j + 1);
                auto gam_XY_r = gamXY(2 * j);
                auto gam_XY_i = gamXY(2 * j + 1);
                auto gam_XZ_r = gamXZ(2 * j);
                auto gam_XZ_i = gamXZ(2 * j + 1);
                auto gam_YX_r = gamYX(2 * j);
                auto gam_YX_i = gamYX(2 * j + 1);
                auto gam_YY_r = gamYY(2 * j);
                auto gam_YY_i = gamYY(2 * j + 1);
                auto gam_YZ_r = gamYZ(2 * j);
                auto gam_YZ_i = gamYZ(2 * j + 1);
                auto gam_ZX_r = gamZX(2 * j);
                auto gam_ZX_i = gamZX(2 * j + 1);
                auto gam_ZY_r = gamZY(2 * j);
                auto gam_ZY_i = gamZY(2 * j + 1);
                auto gam_ZZ_r = gamZZ(2 * j);
                auto gam_ZZ_i = gamZZ(2 * j + 1);

                // Input densities 

                auto rhoB_r = rwDensityGrid.alphaDensity(12 * j);
                auto gradB_x_r = rwDensityGrid.alphaDensityGradientX(12 * j);
                auto gradB_y_r = rwDensityGrid.alphaDensityGradientY(12 * j);
                auto gradB_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j);

                auto rhoB_i = rwDensityGrid.alphaDensity(12 * j + 1);
                auto gradB_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 1);
                auto gradB_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 1);
                auto gradB_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 1);

                auto rhoC_r = rwDensityGrid.alphaDensity(12 * j + 2);
                auto gradC_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 2);
                auto gradC_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 2);
                auto gradC_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 2);

                auto rhoC_i = rwDensityGrid.alphaDensity(12 * j + 3);
                auto gradC_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 3);
                auto gradC_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 3);
                auto gradC_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 3);

                auto rhoD_r = rwDensityGrid.alphaDensity(12 * j + 4);
                auto gradD_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 4);
                auto gradD_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 4);
                auto gradD_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 4);

                auto rhoD_i = rwDensityGrid.alphaDensity(12 * j + 5);
                auto gradD_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 5);
                auto gradD_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 5);
                auto gradD_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 5);

                auto rhoBC_r = rwDensityGrid.alphaDensity(12 * j + 6);
                auto gradBC_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 6);
                auto gradBC_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 6);
                auto gradBC_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 6);

                auto rhoBC_i = rwDensityGrid.alphaDensity(12 * j + 7);
                auto gradBC_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 7);
                auto gradBC_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 7);
                auto gradBC_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 7);

                auto rhoBD_r = rwDensityGrid.alphaDensity(12 * j + 8);
                auto gradBD_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 8);
                auto gradBD_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 8);
                auto gradBD_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 8);

                auto rhoBD_i = rwDensityGrid.alphaDensity(12 * j + 9);
                auto gradBD_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 9);
                auto gradBD_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 9);
                auto gradBD_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 9);
        
                auto rhoCD_r = rwDensityGrid.alphaDensity(12 * j + 10);
                auto gradCD_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 10);
                auto gradCD_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 10);
                auto gradCD_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 10);

                auto rhoCD_i = rwDensityGrid.alphaDensity(12 * j + 11);
                auto gradCD_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 11);
                auto gradCD_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 11);
                auto gradCD_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 11);


                for (int32_t i = 0; i < npoints; i++)
                {

                    gam_r[i] += 2.0 * prod2_r(rhoBC_r[i],rhoBC_i[i],rhoD_r[i],rhoD_i[i])
                              + 2.0 * prod2_r(rhoBD_r[i],rhoBD_i[i],rhoC_r[i],rhoC_i[i])
                              + 2.0 * prod2_r(rhoCD_r[i],rhoCD_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_i[i] += 2.0 * prod2_i(rhoBC_r[i],rhoBC_i[i],rhoD_r[i],rhoD_i[i])
                              + 2.0 * prod2_i(rhoBD_r[i],rhoBD_i[i],rhoC_r[i],rhoC_i[i])
                              + 2.0 * prod2_i(rhoCD_r[i],rhoCD_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_X_r[i] += 2.0 * prod2_r(gradBC_x_r[i],gradBC_x_i[i],rhoD_r[i],rhoD_i[i])
                                + 2.0 * prod2_r(gradBD_x_r[i],gradBD_x_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_r(gradCD_x_r[i],gradCD_x_i[i],rhoB_r[i],rhoB_i[i])
                                + 2.0 * prod2_r(gradD_x_r[i],gradD_x_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 2.0 * prod2_r(gradC_x_r[i],gradC_x_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 2.0 * prod2_r(gradB_x_r[i],gradB_x_i[i],rhoCD_r[i],rhoCD_i[i]);                   

                    gam_Y_r[i] += 2.0 * prod2_r(gradBC_y_r[i],gradBC_y_i[i],rhoD_r[i],rhoD_i[i])
                                + 2.0 * prod2_r(gradBD_y_r[i],gradBD_y_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_r(gradCD_y_r[i],gradCD_y_i[i],rhoB_r[i],rhoB_i[i])
                                + 2.0 * prod2_r(gradD_y_r[i],gradD_y_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 2.0 * prod2_r(gradC_y_r[i],gradC_y_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 2.0 * prod2_r(gradB_y_r[i],gradB_y_i[i],rhoCD_r[i],rhoCD_i[i]);                   

                    gam_Z_r[i] += 2.0 * prod2_r(gradBC_z_r[i],gradBC_z_i[i],rhoD_r[i],rhoD_i[i])
                                + 2.0 * prod2_r(gradBD_z_r[i],gradBD_z_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_r(gradCD_z_r[i],gradCD_z_i[i],rhoB_r[i],rhoB_i[i])
                                + 2.0 * prod2_r(gradD_z_r[i],gradD_z_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 2.0 * prod2_r(gradC_z_r[i],gradC_z_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 2.0 * prod2_r(gradB_z_r[i],gradB_z_i[i],rhoCD_r[i],rhoCD_i[i]);                   

                    gam_X_i[i] += 2.0 * prod2_i(gradBC_x_r[i],gradBC_x_i[i],rhoD_r[i],rhoD_i[i])
                                + 2.0 * prod2_i(gradBD_x_r[i],gradBD_x_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_i(gradCD_x_r[i],gradCD_x_i[i],rhoB_r[i],rhoB_i[i])
                                + 2.0 * prod2_i(gradD_x_r[i],gradD_x_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 2.0 * prod2_i(gradC_x_r[i],gradC_x_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 2.0 * prod2_i(gradB_x_r[i],gradB_x_i[i],rhoCD_r[i],rhoCD_i[i]);                   

                    gam_Y_i[i] += 2.0 * prod2_i(gradBC_y_r[i],gradBC_y_i[i],rhoD_r[i],rhoD_i[i])
                                + 2.0 * prod2_i(gradBD_y_r[i],gradBD_y_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_i(gradCD_y_r[i],gradCD_y_i[i],rhoB_r[i],rhoB_i[i])
                                + 2.0 * prod2_i(gradD_y_r[i],gradD_y_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 2.0 * prod2_i(gradC_y_r[i],gradC_y_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 2.0 * prod2_i(gradB_y_r[i],gradB_y_i[i],rhoCD_r[i],rhoCD_i[i]);                   

                    gam_Z_i[i] += 2.0 * prod2_i(gradBC_z_r[i],gradBC_z_i[i],rhoD_r[i],rhoD_i[i])
                                + 2.0 * prod2_i(gradBD_z_r[i],gradBD_z_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_i(gradCD_z_r[i],gradCD_z_i[i],rhoB_r[i],rhoB_i[i])
                                + 2.0 * prod2_i(gradD_z_r[i],gradD_z_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 2.0 * prod2_i(gradC_z_r[i],gradC_z_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 2.0 * prod2_i(gradB_z_r[i],gradB_z_i[i],rhoCD_r[i],rhoCD_i[i]);                   

                    gam_XX_r[i] += prod2_r(gradBC_x_r[i],gradBC_x_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_r(gradBD_x_r[i],gradBD_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gradCD_x_r[i],gradCD_x_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_r(gradBC_x_r[i],gradBC_x_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_r(gradBD_x_r[i],gradBD_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gradCD_x_r[i],gradCD_x_i[i],gradB_x_r[i],gradB_x_i[i]);
                    
                    gam_XY_r[i] += prod2_r(gradBC_x_r[i],gradBC_x_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_r(gradBD_x_r[i],gradBD_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gradCD_x_r[i],gradCD_x_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_r(gradBC_y_r[i],gradBC_y_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_r(gradBD_y_r[i],gradBD_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gradCD_y_r[i],gradCD_y_i[i],gradB_x_r[i],gradB_x_i[i]);
                    
                    gam_XZ_r[i] += prod2_r(gradBC_x_r[i],gradBC_x_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_r(gradBD_x_r[i],gradBD_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gradCD_x_r[i],gradCD_x_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_r(gradBC_z_r[i],gradBC_z_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_r(gradBD_z_r[i],gradBD_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gradCD_z_r[i],gradCD_z_i[i],gradB_x_r[i],gradB_x_i[i]);
                    
                    gam_YX_r[i] += prod2_r(gradBC_y_r[i],gradBC_y_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_r(gradBD_y_r[i],gradBD_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gradCD_y_r[i],gradCD_y_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_r(gradBC_x_r[i],gradBC_x_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_r(gradBD_x_r[i],gradBD_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gradCD_x_r[i],gradCD_x_i[i],gradB_y_r[i],gradB_y_i[i]);
                    
                    gam_YY_r[i] += prod2_r(gradBC_y_r[i],gradBC_y_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_r(gradBD_y_r[i],gradBD_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gradCD_y_r[i],gradCD_y_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_r(gradBC_y_r[i],gradBC_y_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_r(gradBD_y_r[i],gradBD_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gradCD_y_r[i],gradCD_y_i[i],gradB_y_r[i],gradB_y_i[i]);
                    
                    gam_YZ_r[i] += prod2_r(gradBC_y_r[i],gradBC_y_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_r(gradBD_y_r[i],gradBD_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gradCD_y_r[i],gradCD_y_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_r(gradBC_z_r[i],gradBC_z_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_r(gradBD_z_r[i],gradBD_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gradCD_z_r[i],gradCD_z_i[i],gradB_y_r[i],gradB_y_i[i]);
                    
                    gam_ZX_r[i] += prod2_r(gradBC_z_r[i],gradBC_z_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_r(gradBD_z_r[i],gradBD_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gradCD_z_r[i],gradCD_z_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_r(gradBC_x_r[i],gradBC_x_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_r(gradBD_x_r[i],gradBD_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gradCD_x_r[i],gradCD_x_i[i],gradB_z_r[i],gradB_z_i[i]);
                    
                    gam_ZY_r[i] += prod2_r(gradBC_z_r[i],gradBC_z_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_r(gradBD_z_r[i],gradBD_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gradCD_z_r[i],gradCD_z_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_r(gradBC_y_r[i],gradBC_y_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_r(gradBD_y_r[i],gradBD_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gradCD_y_r[i],gradCD_y_i[i],gradB_z_r[i],gradB_z_i[i]);
                    
                    gam_ZZ_r[i] += prod2_r(gradBC_z_r[i],gradBC_z_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_r(gradBD_z_r[i],gradBD_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gradCD_z_r[i],gradCD_z_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_r(gradBC_z_r[i],gradBC_z_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_r(gradBD_z_r[i],gradBD_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gradCD_z_r[i],gradCD_z_i[i],gradB_z_r[i],gradB_z_i[i]);
                    
                    gam_XX_i[i] += prod2_i(gradBC_x_r[i],gradBC_x_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_i(gradBD_x_r[i],gradBD_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gradCD_x_r[i],gradCD_x_i[i],gradB_x_r[i],gradB_x_i[i])
                                  + prod2_i(gradBC_x_r[i],gradBC_x_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_i(gradBD_x_r[i],gradBD_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gradCD_x_r[i],gradCD_x_i[i],gradB_x_r[i],gradB_x_i[i]);
                    
                    gam_XY_i[i] += prod2_i(gradBC_x_r[i],gradBC_x_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_i(gradBD_x_r[i],gradBD_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gradCD_x_r[i],gradCD_x_i[i],gradB_y_r[i],gradB_y_i[i])
                                  + prod2_i(gradBC_y_r[i],gradBC_y_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_i(gradBD_y_r[i],gradBD_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gradCD_y_r[i],gradCD_y_i[i],gradB_x_r[i],gradB_x_i[i]);
                    
                    gam_XZ_i[i] += prod2_i(gradBC_x_r[i],gradBC_x_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_i(gradBD_x_r[i],gradBD_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gradCD_x_r[i],gradCD_x_i[i],gradB_z_r[i],gradB_z_i[i])
                                  + prod2_i(gradBC_z_r[i],gradBC_z_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_i(gradBD_z_r[i],gradBD_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gradCD_z_r[i],gradCD_z_i[i],gradB_x_r[i],gradB_x_i[i]);
                    
                    gam_YX_i[i] += prod2_i(gradBC_y_r[i],gradBC_y_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_i(gradBD_y_r[i],gradBD_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gradCD_y_r[i],gradCD_y_i[i],gradB_x_r[i],gradB_x_i[i])
                                  + prod2_i(gradBC_x_r[i],gradBC_x_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_i(gradBD_x_r[i],gradBD_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gradCD_x_r[i],gradCD_x_i[i],gradB_y_r[i],gradB_y_i[i]);
                    
                    gam_YY_i[i] += prod2_i(gradBC_y_r[i],gradBC_y_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_i(gradBD_y_r[i],gradBD_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gradCD_y_r[i],gradCD_y_i[i],gradB_y_r[i],gradB_y_i[i])
                                  + prod2_i(gradBC_y_r[i],gradBC_y_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_i(gradBD_y_r[i],gradBD_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gradCD_y_r[i],gradCD_y_i[i],gradB_y_r[i],gradB_y_i[i]);
                    
                    gam_YZ_i[i] += prod2_i(gradBC_y_r[i],gradBC_y_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_i(gradBD_y_r[i],gradBD_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gradCD_y_r[i],gradCD_y_i[i],gradB_z_r[i],gradB_z_i[i])
                                  + prod2_i(gradBC_z_r[i],gradBC_z_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_i(gradBD_z_r[i],gradBD_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gradCD_z_r[i],gradCD_z_i[i],gradB_y_r[i],gradB_y_i[i]);
                    
                    gam_ZX_i[i] += prod2_i(gradBC_z_r[i],gradBC_z_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_i(gradBD_z_r[i],gradBD_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gradCD_z_r[i],gradCD_z_i[i],gradB_x_r[i],gradB_x_i[i])
                                  + prod2_i(gradBC_x_r[i],gradBC_x_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_i(gradBD_x_r[i],gradBD_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gradCD_x_r[i],gradCD_x_i[i],gradB_z_r[i],gradB_z_i[i]);
                    
                    gam_ZY_i[i] += prod2_i(gradBC_z_r[i],gradBC_z_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_i(gradBD_z_r[i],gradBD_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gradCD_z_r[i],gradCD_z_i[i],gradB_y_r[i],gradB_y_i[i])
                                  + prod2_i(gradBC_y_r[i],gradBC_y_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_i(gradBD_y_r[i],gradBD_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gradCD_y_r[i],gradCD_y_i[i],gradB_z_r[i],gradB_z_i[i]);
                    
                    gam_ZZ_i[i] += prod2_i(gradBC_z_r[i],gradBC_z_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_i(gradBD_z_r[i],gradBD_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gradCD_z_r[i],gradCD_z_i[i],gradB_z_r[i],gradB_z_i[i])
                                  + prod2_i(gradBC_z_r[i],gradBC_z_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_i(gradBD_z_r[i],gradBD_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gradCD_z_r[i],gradCD_z_i[i],gradB_z_r[i],gradB_z_i[i]);
                   
                }
            }
        }
        if (fstr::upcase(quadMode) == "SHG")
        {
            for (int32_t j = 0; j < numdens / 12; j++)

            {
                // Density products to be stored

                // RhoRho part

                auto rho_sig_x_r = gam(12 * j);

                auto rho_sig_x_i = gam(12 * j + 1);

                auto rho_sig_y_r = gam(12 * j + 2);

                auto rho_sig_y_i = gam(12 * j + 3);

                auto rho_sig_z_r = gam(12 * j + 4);

                auto rho_sig_z_i = gam(12 * j + 5);

                auto rho_lam_xy_r = gam(12 * j + 6);

                auto rho_lam_xy_i = gam(12 * j + 7);

                auto rho_lam_xz_r = gam(12 * j + 8);

                auto rho_lam_xz_i = gam(12 * j + 9);

                auto rho_lam_yz_r = gam(12 * j + 10);

                auto rho_lam_yz_i = gam(12 * j + 11);

                // rxw1rhow part

                auto gamX_sig_x_r = gamX(12 * j);

                auto gamX_sig_x_i = gamX(12 * j + 1);

                auto gamX_sig_y_r = gamX(12 * j + 2);

                auto gamX_sig_y_i = gamX(12 * j + 3);

                auto gamX_sig_z_r = gamX(12 * j + 4);

                auto gamX_sig_z_i = gamX(12 * j + 5);

                auto gamX_lam_xy_r = gamX(12 * j + 6);

                auto gamX_lam_xy_i = gamX(12 * j + 7);

                auto gamX_lam_xz_r = gamX(12 * j + 8);

                auto gamX_lam_xz_i = gamX(12 * j + 9);

                auto gamX_lam_yz_r = gamX(12 * j + 10);

                auto gamX_lam_yz_i = gamX(12 * j + 11);

                // ryw1rhow part

                auto gamY_sig_x_r = gamY(12 * j);

                auto gamY_sig_x_i = gamY(12 * j + 1);

                auto gamY_sig_y_r = gamY(12 * j + 2);

                auto gamY_sig_y_i = gamY(12 * j + 3);

                auto gamY_sig_z_r = gamY(12 * j + 4);

                auto gamY_sig_z_i = gamY(12 * j + 5);

                auto gamY_lam_xy_r = gamY(12 * j + 6);

                auto gamY_lam_xy_i = gamY(12 * j + 7);

                auto gamY_lam_xz_r = gamY(12 * j + 8);

                auto gamY_lam_xz_i = gamY(12 * j + 9);

                auto gamY_lam_yz_r = gamY(12 * j + 10);

                auto gamY_lam_yz_i = gamY(12 * j + 11);

                // rzw1rhow part

                auto gamZ_sig_x_r = gamZ(12 * j);

                auto gamZ_sig_x_i = gamZ(12 * j + 1);

                auto gamZ_sig_y_r = gamZ(12 * j + 2);

                auto gamZ_sig_y_i = gamZ(12 * j + 3);

                auto gamZ_sig_z_r = gamZ(12 * j + 4);

                auto gamZ_sig_z_i = gamZ(12 * j + 5);

                auto gamZ_lam_xy_r = gamZ(12 * j + 6);

                auto gamZ_lam_xy_i = gamZ(12 * j + 7);

                auto gamZ_lam_xz_r = gamZ(12 * j + 8);

                auto gamZ_lam_xz_i = gamZ(12 * j + 9);

                auto gamZ_lam_yz_r = gamZ(12 * j + 10);

                auto gamZ_lam_yz_i = gamZ(12 * j + 11);

                // gamXX part

                auto gamXX_sig_x_r = gamXX(12 * j);

                auto gamXX_sig_x_i = gamXX(12 * j + 1);

                auto gamXX_sig_y_r = gamXX(12 * j + 2);

                auto gamXX_sig_y_i = gamXX(12 * j + 3);

                auto gamXX_sig_z_r = gamXX(12 * j + 4);

                auto gamXX_sig_z_i = gamXX(12 * j + 5);

                auto gamXX_lam_xy_r = gamXX(12 * j + 6);

                auto gamXX_lam_xy_i = gamXX(12 * j + 7);

                auto gamXX_lam_xz_r = gamXX(12 * j + 8);

                auto gamXX_lam_xz_i = gamXX(12 * j + 9);

                auto gamXX_lam_yz_r = gamXX(12 * j + 10);

                auto gamXX_lam_yz_i = gamXX(12 * j + 11);

                // gamXY part

                auto gamXY_sig_x_r = gamXY(12 * j);

                auto gamXY_sig_x_i = gamXY(12 * j + 1);

                auto gamXY_sig_y_r = gamXY(12 * j + 2);

                auto gamXY_sig_y_i = gamXY(12 * j + 3);

                auto gamXY_sig_z_r = gamXY(12 * j + 4);

                auto gamXY_sig_z_i = gamXY(12 * j + 5);

                auto gamXY_lam_xy_r = gamXY(12 * j + 6);

                auto gamXY_lam_xy_i = gamXY(12 * j + 7);

                auto gamXY_lam_xz_r = gamXY(12 * j + 8);

                auto gamXY_lam_xz_i = gamXY(12 * j + 9);

                auto gamXY_lam_yz_r = gamXY(12 * j + 10);

                auto gamXY_lam_yz_i = gamXY(12 * j + 11);

                // gamXZ part

                auto gamXZ_sig_x_r = gamXZ(12 * j);

                auto gamXZ_sig_x_i = gamXZ(12 * j + 1);

                auto gamXZ_sig_y_r = gamXZ(12 * j + 2);

                auto gamXZ_sig_y_i = gamXZ(12 * j + 3);

                auto gamXZ_sig_z_r = gamXZ(12 * j + 4);

                auto gamXZ_sig_z_i = gamXZ(12 * j + 5);

                auto gamXZ_lam_xy_r = gamXZ(12 * j + 6);

                auto gamXZ_lam_xy_i = gamXZ(12 * j + 7);

                auto gamXZ_lam_xz_r = gamXZ(12 * j + 8);

                auto gamXZ_lam_xz_i = gamXZ(12 * j + 9);

                auto gamXZ_lam_yz_r = gamXZ(12 * j + 10);

                auto gamXZ_lam_yz_i = gamXZ(12 * j + 11);

                // gamYX part

                auto gamYX_sig_x_r = gamYX(12 * j);

                auto gamYX_sig_x_i = gamYX(12 * j + 1);

                auto gamYX_sig_y_r = gamYX(12 * j + 2);

                auto gamYX_sig_y_i = gamYX(12 * j + 3);

                auto gamYX_sig_z_r = gamYX(12 * j + 4);

                auto gamYX_sig_z_i = gamYX(12 * j + 5);

                auto gamYX_lam_xy_r = gamYX(12 * j + 6);

                auto gamYX_lam_xy_i = gamYX(12 * j + 7);

                auto gamYX_lam_xz_r = gamYX(12 * j + 8);

                auto gamYX_lam_xz_i = gamYX(12 * j + 9);

                auto gamYX_lam_yz_r = gamYX(12 * j + 10);

                auto gamYX_lam_yz_i = gamYX(12 * j + 11);

                // gamYY part

                auto gamYY_sig_x_r = gamYY(12 * j);

                auto gamYY_sig_x_i = gamYY(12 * j + 1);

                auto gamYY_sig_y_r = gamYY(12 * j + 2);

                auto gamYY_sig_y_i = gamYY(12 * j + 3);

                auto gamYY_sig_z_r = gamYY(12 * j + 4);

                auto gamYY_sig_z_i = gamYY(12 * j + 5);

                auto gamYY_lam_xy_r = gamYY(12 * j + 6);

                auto gamYY_lam_xy_i = gamYY(12 * j + 7);

                auto gamYY_lam_xz_r = gamYY(12 * j + 8);

                auto gamYY_lam_xz_i = gamYY(12 * j + 9);

                auto gamYY_lam_yz_r = gamYY(12 * j + 10);

                auto gamYY_lam_yz_i = gamYY(12 * j + 11);

                // gamYZ part

                auto gamYZ_sig_x_r = gamYZ(12 * j);

                auto gamYZ_sig_x_i = gamYZ(12 * j + 1);

                auto gamYZ_sig_y_r = gamYZ(12 * j + 2);

                auto gamYZ_sig_y_i = gamYZ(12 * j + 3);

                auto gamYZ_sig_z_r = gamYZ(12 * j + 4);

                auto gamYZ_sig_z_i = gamYZ(12 * j + 5);

                auto gamYZ_lam_xy_r = gamYZ(12 * j + 6);

                auto gamYZ_lam_xy_i = gamYZ(12 * j + 7);

                auto gamYZ_lam_xz_r = gamYZ(12 * j + 8);

                auto gamYZ_lam_xz_i = gamYZ(12 * j + 9);

                auto gamYZ_lam_yz_r = gamYZ(12 * j + 10);

                auto gamYZ_lam_yz_i = gamYZ(12 * j + 11);

                // gamZX part

                auto gamZX_sig_x_r = gamZX(12 * j);

                auto gamZX_sig_x_i = gamZX(12 * j + 1);

                auto gamZX_sig_y_r = gamZX(12 * j + 2);

                auto gamZX_sig_y_i = gamZX(12 * j + 3);

                auto gamZX_sig_z_r = gamZX(12 * j + 4);

                auto gamZX_sig_z_i = gamZX(12 * j + 5);

                auto gamZX_lam_xy_r = gamZX(12 * j + 6);

                auto gamZX_lam_xy_i = gamZX(12 * j + 7);

                auto gamZX_lam_xz_r = gamZX(12 * j + 8);

                auto gamZX_lam_xz_i = gamZX(12 * j + 9);

                auto gamZX_lam_yz_r = gamZX(12 * j + 10);

                auto gamZX_lam_yz_i = gamZX(12 * j + 11);

                // gamZY part

                auto gamZY_sig_x_r = gamZY(12 * j);

                auto gamZY_sig_x_i = gamZY(12 * j + 1);

                auto gamZY_sig_y_r = gamZY(12 * j + 2);

                auto gamZY_sig_y_i = gamZY(12 * j + 3);

                auto gamZY_sig_z_r = gamZY(12 * j + 4);

                auto gamZY_sig_z_i = gamZY(12 * j + 5);

                auto gamZY_lam_xy_r = gamZY(12 * j + 6);

                auto gamZY_lam_xy_i = gamZY(12 * j + 7);

                auto gamZY_lam_xz_r = gamZY(12 * j + 8);

                auto gamZY_lam_xz_i = gamZY(12 * j + 9);

                auto gamZY_lam_yz_r = gamZY(12 * j + 10);

                auto gamZY_lam_yz_i = gamZY(12 * j + 11);

                // gamZZ part

                auto gamZZ_sig_x_r = gamZZ(12 * j);

                auto gamZZ_sig_x_i = gamZZ(12 * j + 1);

                auto gamZZ_sig_y_r = gamZZ(12 * j + 2);

                auto gamZZ_sig_y_i = gamZZ(12 * j + 3);

                auto gamZZ_sig_z_r = gamZZ(12 * j + 4);

                auto gamZZ_sig_z_i = gamZZ(12 * j + 5);

                auto gamZZ_lam_xy_r = gamZZ(12 * j + 6);

                auto gamZZ_lam_xy_i = gamZZ(12 * j + 7);

                auto gamZZ_lam_xz_r = gamZZ(12 * j + 8);

                auto gamZZ_lam_xz_i = gamZZ(12 * j + 9);

                auto gamZZ_lam_yz_r = gamZZ(12 * j + 10);

                auto gamZZ_lam_yz_i = gamZZ(12 * j + 11);

                // First-order densities

                // Rho terms

                auto rhow_kx_r = rwDensityGrid.alphaDensity(6 * j);

                auto rhow_kx_i = rwDensityGrid.alphaDensity(6 * j + 1);

                auto rhow_ky_r = rwDensityGrid.alphaDensity(6 * j + 2);

                auto rhow_ky_i = rwDensityGrid.alphaDensity(6 * j + 3);

                auto rhow_kz_r = rwDensityGrid.alphaDensity(6 * j + 4);

                auto rhow_kz_i = rwDensityGrid.alphaDensity(6 * j + 5);

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

                    // gamX real

                    jj_r = 4.0 * (gradxw_kx_r[i] * rhow_kx_r[i] + gradxw_ky_r[i] * rhow_ky_r[i] + gradxw_kz_r[i] * rhow_kz_r[i])

                           - 4.0 * (gradxw_kx_i[i] * rhow_kx_i[i] + gradxw_ky_i[i] * rhow_ky_i[i] + gradxw_kz_i[i] * rhow_kz_i[i]);

                    gamX_sig_x_r[i] = 8.0 * (gradxw_kx_r[i] * rhow_kx_r[i] - gradxw_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    gamX_sig_y_r[i] = 8.0 * (gradxw_ky_r[i] * rhow_ky_r[i] - gradxw_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    gamX_sig_z_r[i] = 8.0 * (gradxw_kz_r[i] * rhow_kz_r[i] - gradxw_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    gamX_lam_xy_r[i] = 2.0 * (gradxw_kx_r[i] * rhow_ky_r[i] - gradxw_kx_i[i] * rhow_ky_i[i]

                                                   + gradxw_ky_r[i] * rhow_kx_r[i] - gradxw_ky_i[i] * rhow_kx_i[i]);

                    gamX_lam_xz_r[i] = 2.0 * (gradxw_kx_r[i] * rhow_kz_r[i] - gradxw_kx_i[i] * rhow_kz_i[i]

                                                   + gradxw_kz_r[i] * rhow_kx_r[i] - gradxw_kz_i[i] * rhow_kx_i[i]);

                    gamX_lam_yz_r[i] = 2.0 * (gradxw_ky_r[i] * rhow_kz_r[i] - gradxw_ky_i[i] * rhow_kz_i[i]

                                                   + gradxw_kz_r[i] * rhow_ky_r[i] - gradxw_kz_i[i] * rhow_ky_i[i]);

                    // gamX imag

                    jj_i = 4.0 * (gradxw_kx_r[i] * rhow_kx_i[i] + gradxw_ky_r[i] * rhow_ky_i[i] + gradxw_kz_r[i] * rhow_kz_i[i]

                                  + gradxw_kx_i[i] * rhow_kx_r[i] + gradxw_ky_i[i] * rhow_ky_r[i] + gradxw_kz_i[i] * rhow_kz_r[i]);

                    gamX_sig_x_i[i] = 8.0 * (gradxw_kx_r[i] * rhow_kx_i[i] + gradxw_kx_i[i] * rhow_kx_r[i]) + jj_i;

                    gamX_sig_y_i[i] = 8.0 * (gradxw_ky_r[i] * rhow_ky_i[i] + gradxw_ky_i[i] * rhow_ky_r[i]) + jj_i;

                    gamX_sig_z_i[i] = 8.0 * (gradxw_kz_r[i] * rhow_kz_i[i] + gradxw_kz_i[i] * rhow_kz_r[i]) + jj_i;

                    gamX_lam_xy_i[i] = 2.0 * (gradxw_kx_r[i] * rhow_ky_i[i] + gradxw_kx_i[i] * rhow_ky_r[i]

                                                   + gradxw_ky_r[i] * rhow_kx_i[i] + gradxw_ky_i[i] * rhow_kx_r[i]);

                    gamX_lam_xz_i[i] = 2.0 * (gradxw_kx_r[i] * rhow_kz_i[i] + gradxw_kx_i[i] * rhow_kz_r[i]

                                                   + gradxw_kz_r[i] * rhow_kx_i[i] + gradxw_kz_i[i] * rhow_kx_r[i]);

                    gamX_lam_yz_i[i] = 2.0 * (gradxw_ky_r[i] * rhow_kz_i[i] + gradxw_ky_i[i] * rhow_kz_r[i]

                                                   + gradxw_kz_r[i] * rhow_ky_i[i] + gradxw_kz_i[i] * rhow_ky_r[i]);

                    // gamY real

                    jj_r = 4.0 * (gradyw_kx_r[i] * rhow_kx_r[i] + gradyw_ky_r[i] * rhow_ky_r[i] + gradyw_kz_r[i] * rhow_kz_r[i])

                           - 4.0 * (gradyw_kx_i[i] * rhow_kx_i[i] + gradyw_ky_i[i] * rhow_ky_i[i] + gradyw_kz_i[i] * rhow_kz_i[i]);

                    gamY_sig_x_r[i] = 8.0 * (gradyw_kx_r[i] * rhow_kx_r[i] - gradyw_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    gamY_sig_y_r[i] = 8.0 * (gradyw_ky_r[i] * rhow_ky_r[i] - gradyw_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    gamY_sig_z_r[i] = 8.0 * (gradyw_kz_r[i] * rhow_kz_r[i] - gradyw_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    gamY_lam_xy_r[i] = 2.0 * (gradyw_kx_r[i] * rhow_ky_r[i] - gradyw_kx_i[i] * rhow_ky_i[i]

                                                   + gradyw_ky_r[i] * rhow_kx_r[i] - gradyw_ky_i[i] * rhow_kx_i[i]);

                    gamY_lam_xz_r[i] = 2.0 * (gradyw_kx_r[i] * rhow_kz_r[i] - gradyw_kx_i[i] * rhow_kz_i[i]

                                                   + gradyw_kz_r[i] * rhow_kx_r[i] - gradyw_kz_i[i] * rhow_kx_i[i]);

                    gamY_lam_yz_r[i] = 2.0 * (gradyw_ky_r[i] * rhow_kz_r[i] - gradyw_ky_i[i] * rhow_kz_i[i]

                                                   + gradyw_kz_r[i] * rhow_ky_r[i] - gradyw_kz_i[i] * rhow_ky_i[i]);

                    // gamY imag

                    jj_i = 4.0 * (gradyw_kx_r[i] * rhow_kx_i[i] + gradyw_ky_r[i] * rhow_ky_i[i] + gradyw_kz_r[i] * rhow_kz_i[i]

                                  + gradyw_kx_i[i] * rhow_kx_r[i] + gradyw_ky_i[i] * rhow_ky_r[i] + gradyw_kz_i[i] * rhow_kz_r[i]);

                    gamY_sig_x_i[i] = 8.0 * (gradyw_kx_r[i] * rhow_kx_i[i] + gradyw_kx_i[i] * rhow_kx_r[i]) + jj_i;

                    gamY_sig_y_i[i] = 8.0 * (gradyw_ky_r[i] * rhow_ky_i[i] + gradyw_ky_i[i] * rhow_ky_r[i]) + jj_i;

                    gamY_sig_z_i[i] = 8.0 * (gradyw_kz_r[i] * rhow_kz_i[i] + gradyw_kz_i[i] * rhow_kz_r[i]) + jj_i;

                    gamY_lam_xy_i[i] = 2.0 * (gradyw_kx_r[i] * rhow_ky_i[i] + gradyw_kx_i[i] * rhow_ky_r[i]

                                                   + gradyw_ky_r[i] * rhow_kx_i[i] + gradyw_ky_i[i] * rhow_kx_r[i]);

                    gamY_lam_xz_i[i] = 2.0 * (gradyw_kx_r[i] * rhow_kz_i[i] + gradyw_kx_i[i] * rhow_kz_r[i]

                                                   + gradyw_kz_r[i] * rhow_kx_i[i] + gradyw_kz_i[i] * rhow_kx_r[i]);

                    gamY_lam_yz_i[i] = 2.0 * (gradyw_ky_r[i] * rhow_kz_i[i] + gradyw_ky_i[i] * rhow_kz_r[i]

                                                   + gradyw_kz_r[i] * rhow_ky_i[i] + gradyw_kz_i[i] * rhow_ky_r[i]);

                    // gamZ real

                    jj_r = 4.0 * (gradzw_kx_r[i] * rhow_kx_r[i] + gradzw_ky_r[i] * rhow_ky_r[i] + gradzw_kz_r[i] * rhow_kz_r[i])

                           - 4.0 * (gradzw_kx_i[i] * rhow_kx_i[i] + gradzw_ky_i[i] * rhow_ky_i[i] + gradzw_kz_i[i] * rhow_kz_i[i]);

                    gamZ_sig_x_r[i] = 8.0 * (gradzw_kx_r[i] * rhow_kx_r[i] - gradzw_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    gamZ_sig_y_r[i] = 8.0 * (gradzw_ky_r[i] * rhow_ky_r[i] - gradzw_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    gamZ_sig_z_r[i] = 8.0 * (gradzw_kz_r[i] * rhow_kz_r[i] - gradzw_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    gamZ_lam_xy_r[i] = 2.0 * (gradzw_kx_r[i] * rhow_ky_r[i] - gradzw_kx_i[i] * rhow_ky_i[i]

                                                   + gradzw_ky_r[i] * rhow_kx_r[i] - gradzw_ky_i[i] * rhow_kx_i[i]);

                    gamZ_lam_xz_r[i] = 2.0 * (gradzw_kx_r[i] * rhow_kz_r[i] - gradzw_kx_i[i] * rhow_kz_i[i]

                                                   + gradzw_kz_r[i] * rhow_kx_r[i] - gradzw_kz_i[i] * rhow_kx_i[i]);

                    gamZ_lam_yz_r[i] = 2.0 * (gradzw_ky_r[i] * rhow_kz_r[i] - gradzw_ky_i[i] * rhow_kz_i[i]

                                                   + gradzw_kz_r[i] * rhow_ky_r[i] - gradzw_kz_i[i] * rhow_ky_i[i]);

                    // gamZ imag

                    jj_i = 4.0 * (gradzw_kx_r[i] * rhow_kx_i[i] + gradzw_ky_r[i] * rhow_ky_i[i] + gradzw_kz_r[i] * rhow_kz_i[i]

                                  + gradzw_kx_i[i] * rhow_kx_r[i] + gradzw_ky_i[i] * rhow_ky_r[i] + gradzw_kz_i[i] * rhow_kz_r[i]);

                    gamZ_sig_x_i[i] = 8.0 * (gradzw_kx_r[i] * rhow_kx_i[i] + gradzw_kx_i[i] * rhow_kx_r[i]) + jj_i;

                    gamZ_sig_y_i[i] = 8.0 * (gradzw_ky_r[i] * rhow_ky_i[i] + gradzw_ky_i[i] * rhow_ky_r[i]) + jj_i;

                    gamZ_sig_z_i[i] = 8.0 * (gradzw_kz_r[i] * rhow_kz_i[i] + gradzw_kz_i[i] * rhow_kz_r[i]) + jj_i;

                    gamZ_lam_xy_i[i] = 2.0 * (gradzw_kx_r[i] * rhow_ky_i[i] + gradzw_kx_i[i] * rhow_ky_r[i]

                                                   + gradzw_ky_r[i] * rhow_kx_i[i] + gradzw_ky_i[i] * rhow_kx_r[i]);

                    gamZ_lam_xz_i[i] = 2.0 * (gradzw_kx_r[i] * rhow_kz_i[i] + gradzw_kx_i[i] * rhow_kz_r[i]

                                                   + gradzw_kz_r[i] * rhow_kx_i[i] + gradzw_kz_i[i] * rhow_kx_r[i]);

                    gamZ_lam_yz_i[i] = 2.0 * (gradzw_ky_r[i] * rhow_kz_i[i] + gradzw_ky_i[i] * rhow_kz_r[i]

                                                   + gradzw_kz_r[i] * rhow_ky_i[i] + gradzw_kz_i[i] * rhow_ky_r[i]);

                    // gamXX

                    jj_r = 2.0 * (gradxw_kx_r[i] * gradxw_kx_r[i] + gradxw_ky_r[i] * gradxw_ky_r[i] + gradxw_kz_r[i] * gradxw_kz_r[i])

                           - 2.0 * (gradxw_kx_i[i] * gradxw_kx_i[i] + gradxw_ky_i[i] * gradxw_ky_i[i] + gradxw_kz_i[i] * gradxw_kz_i[i]);

                    gamXX_sig_x_r[i] = 4.0 * (gradxw_kx_r[i] * gradxw_kx_r[i] - gradxw_kx_i[i] * gradxw_kx_i[i]) + jj_r;

                    gamXX_sig_y_r[i] = 4.0 * (gradxw_ky_r[i] * gradxw_ky_r[i] - gradxw_ky_i[i] * gradxw_ky_i[i]) + jj_r;

                    gamXX_sig_z_r[i] = 4.0 * (gradxw_kz_r[i] * gradxw_kz_r[i] - gradxw_kz_i[i] * gradxw_kz_i[i]) + jj_r;

                    gamXX_lam_xy_r[i] = gradxw_kx_r[i] * gradxw_ky_r[i] - gradxw_kx_i[i] * gradxw_ky_i[i]

                                           + gradxw_ky_r[i] * gradxw_kx_r[i] - gradxw_ky_i[i] * gradxw_kx_i[i];

                    gamXX_lam_xz_r[i] = gradxw_kx_r[i] * gradxw_kz_r[i] - gradxw_kx_i[i] * gradxw_kz_i[i]

                                           + gradxw_kz_r[i] * gradxw_kx_r[i] - gradxw_kz_i[i] * gradxw_kx_i[i];

                    gamXX_lam_yz_r[i] = gradxw_ky_r[i] * gradxw_kz_r[i] - gradxw_ky_i[i] * gradxw_kz_i[i]

                                           + gradxw_kz_r[i] * gradxw_ky_r[i] - gradxw_kz_i[i] * gradxw_ky_i[i];

                    // gamXX imag

                    jj_i = 2.0 * (gradxw_kx_r[i] * gradxw_kx_i[i] + gradxw_ky_r[i] * gradxw_ky_i[i] + gradxw_kz_r[i] * gradxw_kz_i[i]

                                  + gradxw_kx_i[i] * gradxw_kx_r[i] + gradxw_ky_i[i] * gradxw_ky_r[i] + gradxw_kz_i[i] * gradxw_kz_r[i]);

                    gamXX_sig_x_i[i] = 4.0 * (gradxw_kx_r[i] * gradxw_kx_i[i] + gradxw_kx_i[i] * gradxw_kx_r[i]) + jj_i;

                    gamXX_sig_y_i[i] = 4.0 * (gradxw_ky_r[i] * gradxw_ky_i[i] + gradxw_ky_i[i] * gradxw_ky_r[i]) + jj_i;

                    gamXX_sig_z_i[i] = 4.0 * (gradxw_kz_r[i] * gradxw_kz_i[i] + gradxw_kz_i[i] * gradxw_kz_r[i]) + jj_i;

                    gamXX_lam_xy_i[i] = gradxw_kx_r[i] * gradxw_ky_i[i] + gradxw_kx_i[i] * gradxw_ky_r[i]

                                           + gradxw_ky_r[i] * gradxw_kx_i[i] + gradxw_ky_i[i] * gradxw_kx_r[i];

                    gamXX_lam_xz_i[i] = gradxw_kx_r[i] * gradxw_kz_i[i] + gradxw_kx_i[i] * gradxw_kz_r[i]

                                           + gradxw_kz_r[i] * gradxw_kx_i[i] + gradxw_kz_i[i] * gradxw_kx_r[i];

                    gamXX_lam_yz_i[i] = gradxw_ky_r[i] * gradxw_kz_i[i] + gradxw_ky_i[i] * gradxw_kz_r[i]

                                           + gradxw_kz_r[i] * gradxw_ky_i[i] + gradxw_kz_i[i] * gradxw_ky_r[i];

                    // gamXY

                    jj_r = 2.0 * (gradxw_kx_r[i] * gradyw_kx_r[i] + gradxw_ky_r[i] * gradyw_ky_r[i] + gradxw_kz_r[i] * gradyw_kz_r[i])

                           - 2.0 * (gradxw_kx_i[i] * gradyw_kx_i[i] + gradxw_ky_i[i] * gradyw_ky_i[i] + gradxw_kz_i[i] * gradyw_kz_i[i]);

                    gamXY_sig_x_r[i] = 4.0 * (gradxw_kx_r[i] * gradyw_kx_r[i] - gradxw_kx_i[i] * gradyw_kx_i[i]) + jj_r;

                    gamXY_sig_y_r[i] = 4.0 * (gradxw_ky_r[i] * gradyw_ky_r[i] - gradxw_ky_i[i] * gradyw_ky_i[i]) + jj_r;

                    gamXY_sig_z_r[i] = 4.0 * (gradxw_kz_r[i] * gradyw_kz_r[i] - gradxw_kz_i[i] * gradyw_kz_i[i]) + jj_r;

                    gamXY_lam_xy_r[i] = gradxw_kx_r[i] * gradyw_ky_r[i] - gradxw_kx_i[i] * gradyw_ky_i[i]

                                           + gradxw_ky_r[i] * gradyw_kx_r[i] - gradxw_ky_i[i] * gradyw_kx_i[i];

                    gamXY_lam_xz_r[i] = gradxw_kx_r[i] * gradyw_kz_r[i] - gradxw_kx_i[i] * gradyw_kz_i[i]

                                           + gradxw_kz_r[i] * gradyw_kx_r[i] - gradxw_kz_i[i] * gradyw_kx_i[i];

                    gamXY_lam_yz_r[i] = gradxw_ky_r[i] * gradyw_kz_r[i] - gradxw_ky_i[i] * gradyw_kz_i[i]

                                           + gradxw_kz_r[i] * gradyw_ky_r[i] - gradxw_kz_i[i] * gradyw_ky_i[i];

                    // gamXY imag

                    jj_i = 2.0 * (gradxw_kx_r[i] * gradyw_kx_i[i] + gradxw_ky_r[i] * gradyw_ky_i[i] + gradxw_kz_r[i] * gradyw_kz_i[i]

                                  + gradxw_kx_i[i] * gradyw_kx_r[i] + gradxw_ky_i[i] * gradyw_ky_r[i] + gradxw_kz_i[i] * gradyw_kz_r[i]);

                    gamXY_sig_x_i[i] = 4.0 * (gradxw_kx_r[i] * gradyw_kx_i[i] + gradxw_kx_i[i] * gradyw_kx_r[i]) + jj_i;

                    gamXY_sig_y_i[i] = 4.0 * (gradxw_ky_r[i] * gradyw_ky_i[i] + gradxw_ky_i[i] * gradyw_ky_r[i]) + jj_i;

                    gamXY_sig_z_i[i] = 4.0 * (gradxw_kz_r[i] * gradyw_kz_i[i] + gradxw_kz_i[i] * gradyw_kz_r[i]) + jj_i;

                    gamXY_lam_xy_i[i] = gradxw_kx_r[i] * gradyw_ky_i[i] + gradxw_kx_i[i] * gradyw_ky_r[i]

                                           + gradxw_ky_r[i] * gradyw_kx_i[i] + gradxw_ky_i[i] * gradyw_kx_r[i];

                    gamXY_lam_xz_i[i] = gradxw_kx_r[i] * gradyw_kz_i[i] + gradxw_kx_i[i] * gradyw_kz_r[i]

                                           + gradxw_kz_r[i] * gradyw_kx_i[i] + gradxw_kz_i[i] * gradyw_kx_r[i];

                    gamXY_lam_yz_i[i] = gradxw_ky_r[i] * gradyw_kz_i[i] + gradxw_ky_i[i] * gradyw_kz_r[i]

                                           + gradxw_kz_r[i] * gradyw_ky_i[i] + gradxw_kz_i[i] * gradyw_ky_r[i];

                    // gamXZ

                    jj_r = 2.0 * (gradxw_kx_r[i] * gradzw_kx_r[i] + gradxw_ky_r[i] * gradzw_ky_r[i] + gradxw_kz_r[i] * gradzw_kz_r[i])

                           - 2.0 * (gradxw_kx_i[i] * gradzw_kx_i[i] + gradxw_ky_i[i] * gradzw_ky_i[i] + gradxw_kz_i[i] * gradzw_kz_i[i]);

                    gamXZ_sig_x_r[i] = 4.0 * (gradxw_kx_r[i] * gradzw_kx_r[i] - gradxw_kx_i[i] * gradzw_kx_i[i]) + jj_r;

                    gamXZ_sig_y_r[i] = 4.0 * (gradxw_ky_r[i] * gradzw_ky_r[i] - gradxw_ky_i[i] * gradzw_ky_i[i]) + jj_r;

                    gamXZ_sig_z_r[i] = 4.0 * (gradxw_kz_r[i] * gradzw_kz_r[i] - gradxw_kz_i[i] * gradzw_kz_i[i]) + jj_r;

                    gamXZ_lam_xy_r[i] = gradxw_kx_r[i] * gradzw_ky_r[i] - gradxw_kx_i[i] * gradzw_ky_i[i]

                                           + gradxw_ky_r[i] * gradzw_kx_r[i] - gradxw_ky_i[i] * gradzw_kx_i[i];

                    gamXZ_lam_xz_r[i] = gradxw_kx_r[i] * gradzw_kz_r[i] - gradxw_kx_i[i] * gradzw_kz_i[i]

                                           + gradxw_kz_r[i] * gradzw_kx_r[i] - gradxw_kz_i[i] * gradzw_kx_i[i];

                    gamXZ_lam_yz_r[i] = gradxw_ky_r[i] * gradzw_kz_r[i] - gradxw_ky_i[i] * gradzw_kz_i[i]

                                           + gradxw_kz_r[i] * gradzw_ky_r[i] - gradxw_kz_i[i] * gradzw_ky_i[i];

                    // gamXZ imag

                    jj_i = 2.0 * (gradxw_kx_r[i] * gradzw_kx_i[i] + gradxw_ky_r[i] * gradzw_ky_i[i] + gradxw_kz_r[i] * gradzw_kz_i[i]

                                  + gradxw_kx_i[i] * gradzw_kx_r[i] + gradxw_ky_i[i] * gradzw_ky_r[i] + gradxw_kz_i[i] * gradzw_kz_r[i]);

                    gamXZ_sig_x_i[i] = 4.0 * (gradxw_kx_r[i] * gradzw_kx_i[i] + gradxw_kx_i[i] * gradzw_kx_r[i]) + jj_i;

                    gamXZ_sig_y_i[i] = 4.0 * (gradxw_ky_r[i] * gradzw_ky_i[i] + gradxw_ky_i[i] * gradzw_ky_r[i]) + jj_i;

                    gamXZ_sig_z_i[i] = 4.0 * (gradxw_kz_r[i] * gradzw_kz_i[i] + gradxw_kz_i[i] * gradzw_kz_r[i]) + jj_i;

                    gamXZ_lam_xy_i[i] = gradxw_kx_r[i] * gradzw_ky_i[i] + gradxw_kx_i[i] * gradzw_ky_r[i]

                                           + gradxw_ky_r[i] * gradzw_kx_i[i] + gradxw_ky_i[i] * gradzw_kx_r[i];

                    gamXZ_lam_xz_i[i] = gradxw_kx_r[i] * gradzw_kz_i[i] + gradxw_kx_i[i] * gradzw_kz_r[i]

                                           + gradxw_kz_r[i] * gradzw_kx_i[i] + gradxw_kz_i[i] * gradzw_kx_r[i];

                    gamXZ_lam_yz_i[i] = gradxw_ky_r[i] * gradzw_kz_i[i] + gradxw_ky_i[i] * gradzw_kz_r[i]

                                           + gradxw_kz_r[i] * gradzw_ky_i[i] + gradxw_kz_i[i] * gradzw_ky_r[i];

                    // gamYX

                    gamYX_sig_x_r[i] = gamXY_sig_x_r[i];

                    gamYX_sig_y_r[i] = gamXY_sig_y_r[i];

                    gamYX_sig_z_r[i] = gamXY_sig_z_r[i];

                    gamYX_lam_xy_r[i] = gamXY_lam_xy_r[i];

                    gamYX_lam_xz_r[i] = gamXY_lam_xz_r[i];

                    gamYX_lam_yz_r[i] = gamXY_lam_yz_r[i];

                    // gamYX imag

                    gamYX_sig_x_i[i] = gamXY_sig_x_i[i];

                    gamYX_sig_y_i[i] = gamXY_sig_y_i[i];

                    gamYX_sig_z_i[i] = gamXY_sig_z_i[i];

                    gamYX_lam_xy_i[i] = gamXY_lam_xy_i[i];

                    gamYX_lam_xz_i[i] = gamXY_lam_xz_i[i];

                    gamYX_lam_yz_i[i] = gamXY_lam_yz_i[i];

                    // gamYY

                    jj_r = 2.0 * (gradyw_kx_r[i] * gradyw_kx_r[i] + gradyw_ky_r[i] * gradyw_ky_r[i] + gradyw_kz_r[i] * gradyw_kz_r[i])

                           - 2.0 * (gradyw_kx_i[i] * gradyw_kx_i[i] + gradyw_ky_i[i] * gradyw_ky_i[i] + gradyw_kz_i[i] * gradyw_kz_i[i]);

                    gamYY_sig_x_r[i] = 4.0 * (gradyw_kx_r[i] * gradyw_kx_r[i] - gradyw_kx_i[i] * gradyw_kx_i[i]) + jj_r;

                    gamYY_sig_y_r[i] = 4.0 * (gradyw_ky_r[i] * gradyw_ky_r[i] - gradyw_ky_i[i] * gradyw_ky_i[i]) + jj_r;

                    gamYY_sig_z_r[i] = 4.0 * (gradyw_kz_r[i] * gradyw_kz_r[i] - gradyw_kz_i[i] * gradyw_kz_i[i]) + jj_r;

                    gamYY_lam_xy_r[i] = gradyw_kx_r[i] * gradyw_ky_r[i] - gradyw_kx_i[i] * gradyw_ky_i[i]

                                           + gradyw_ky_r[i] * gradyw_kx_r[i] - gradyw_ky_i[i] * gradyw_kx_i[i];

                    gamYY_lam_xz_r[i] = gradyw_kx_r[i] * gradyw_kz_r[i] - gradyw_kx_i[i] * gradyw_kz_i[i]

                                           + gradyw_kz_r[i] * gradyw_kx_r[i] - gradyw_kz_i[i] * gradyw_kx_i[i];

                    gamYY_lam_yz_r[i] = gradyw_ky_r[i] * gradyw_kz_r[i] - gradyw_ky_i[i] * gradyw_kz_i[i]

                                           + gradyw_kz_r[i] * gradyw_ky_r[i] - gradyw_kz_i[i] * gradyw_ky_i[i];

                    // gamYY imag

                    jj_i = 2.0 * (gradyw_kx_r[i] * gradyw_kx_i[i] + gradyw_ky_r[i] * gradyw_ky_i[i] + gradyw_kz_r[i] * gradyw_kz_i[i]

                                  + gradyw_kx_i[i] * gradyw_kx_r[i] + gradyw_ky_i[i] * gradyw_ky_r[i] + gradyw_kz_i[i] * gradyw_kz_r[i]);

                    gamYY_sig_x_i[i] = 4.0 * (gradyw_kx_r[i] * gradyw_kx_i[i] + gradyw_kx_i[i] * gradyw_kx_r[i]) + jj_i;

                    gamYY_sig_y_i[i] = 4.0 * (gradyw_ky_r[i] * gradyw_ky_i[i] + gradyw_ky_i[i] * gradyw_ky_r[i]) + jj_i;

                    gamYY_sig_z_i[i] = 4.0 * (gradyw_kz_r[i] * gradyw_kz_i[i] + gradyw_kz_i[i] * gradyw_kz_r[i]) + jj_i;

                    gamYY_lam_xy_i[i] = gradyw_kx_r[i] * gradyw_ky_i[i] + gradyw_kx_i[i] * gradyw_ky_r[i]

                                           + gradyw_ky_r[i] * gradyw_kx_i[i] + gradyw_ky_i[i] * gradyw_kx_r[i];

                    gamYY_lam_xz_i[i] = gradyw_kx_r[i] * gradyw_kz_i[i] + gradyw_kx_i[i] * gradyw_kz_r[i]

                                           + gradyw_kz_r[i] * gradyw_kx_i[i] + gradyw_kz_i[i] * gradyw_kx_r[i];

                    gamYY_lam_yz_i[i] = gradyw_ky_r[i] * gradyw_kz_i[i] + gradyw_ky_i[i] * gradyw_kz_r[i]

                                           + gradyw_kz_r[i] * gradyw_ky_i[i] + gradyw_kz_i[i] * gradyw_ky_r[i];

                    // gamYZ

                    jj_r = 2.0 * (gradyw_kx_r[i] * gradzw_kx_r[i] + gradyw_ky_r[i] * gradzw_ky_r[i] + gradyw_kz_r[i] * gradzw_kz_r[i])

                           - 2.0 * (gradyw_kx_i[i] * gradzw_kx_i[i] + gradyw_ky_i[i] * gradzw_ky_i[i] + gradyw_kz_i[i] * gradzw_kz_i[i]);

                    gamYZ_sig_x_r[i] = 4.0 * (gradyw_kx_r[i] * gradzw_kx_r[i] - gradyw_kx_i[i] * gradzw_kx_i[i]) + jj_r;

                    gamYZ_sig_y_r[i] = 4.0 * (gradyw_ky_r[i] * gradzw_ky_r[i] - gradyw_ky_i[i] * gradzw_ky_i[i]) + jj_r;

                    gamYZ_sig_z_r[i] = 4.0 * (gradyw_kz_r[i] * gradzw_kz_r[i] - gradyw_kz_i[i] * gradzw_kz_i[i]) + jj_r;

                    gamYZ_lam_xy_r[i] = gradyw_kx_r[i] * gradzw_ky_r[i] - gradyw_kx_i[i] * gradzw_ky_i[i]

                                           + gradyw_ky_r[i] * gradzw_kx_r[i] - gradyw_ky_i[i] * gradzw_kx_i[i];

                    gamYZ_lam_xz_r[i] = gradyw_kx_r[i] * gradzw_kz_r[i] - gradyw_kx_i[i] * gradzw_kz_i[i]

                                           + gradyw_kz_r[i] * gradzw_kx_r[i] - gradyw_kz_i[i] * gradzw_kx_i[i];

                    gamYZ_lam_yz_r[i] = gradyw_ky_r[i] * gradzw_kz_r[i] - gradyw_ky_i[i] * gradzw_kz_i[i]

                                           + gradyw_kz_r[i] * gradzw_ky_r[i] - gradyw_kz_i[i] * gradzw_ky_i[i];

                    // gamYZ imag

                    jj_i = 2.0 * (gradyw_kx_r[i] * gradzw_kx_i[i] + gradyw_ky_r[i] * gradzw_ky_i[i] + gradyw_kz_r[i] * gradzw_kz_i[i]

                                  + gradyw_kx_i[i] * gradzw_kx_r[i] + gradyw_ky_i[i] * gradzw_ky_r[i] + gradyw_kz_i[i] * gradzw_kz_r[i]);

                    gamYZ_sig_x_i[i] = 4.0 * (gradyw_kx_r[i] * gradzw_kx_i[i] + gradyw_kx_i[i] * gradzw_kx_r[i]) + jj_i;

                    gamYZ_sig_y_i[i] = 4.0 * (gradyw_ky_r[i] * gradzw_ky_i[i] + gradyw_ky_i[i] * gradzw_ky_r[i]) + jj_i;

                    gamYZ_sig_z_i[i] = 4.0 * (gradyw_kz_r[i] * gradzw_kz_i[i] + gradyw_kz_i[i] * gradzw_kz_r[i]) + jj_i;

                    gamYZ_lam_xy_i[i] = gradyw_kx_r[i] * gradzw_ky_i[i] + gradyw_kx_i[i] * gradzw_ky_r[i]

                                           + gradyw_ky_r[i] * gradzw_kx_i[i] + gradyw_ky_i[i] * gradzw_kx_r[i];

                    gamYZ_lam_xz_i[i] = gradyw_kx_r[i] * gradzw_kz_i[i] + gradyw_kx_i[i] * gradzw_kz_r[i]

                                           + gradyw_kz_r[i] * gradzw_kx_i[i] + gradyw_kz_i[i] * gradzw_kx_r[i];

                    gamYZ_lam_yz_i[i] = gradyw_ky_r[i] * gradzw_kz_i[i] + gradyw_ky_i[i] * gradzw_kz_r[i]

                                           + gradyw_kz_r[i] * gradzw_ky_i[i] + gradyw_kz_i[i] * gradzw_ky_r[i];

                    // gamZX

                    gamZX_sig_x_r[i] = gamXZ_sig_x_r[i];

                    gamZX_sig_y_r[i] = gamXZ_sig_y_r[i];

                    gamZX_sig_z_r[i] = gamXZ_sig_z_r[i];

                    gamZX_lam_xy_r[i] = gamXZ_lam_xy_r[i];

                    gamZX_lam_xz_r[i] = gamXZ_lam_xz_r[i];

                    gamZX_lam_yz_r[i] = gamXZ_lam_yz_r[i];

                    // gamZX imag

                    gamZX_sig_x_i[i] = gamXZ_sig_x_i[i];

                    gamZX_sig_y_i[i] = gamXZ_sig_y_i[i];

                    gamZX_sig_z_i[i] = gamXZ_sig_z_i[i];

                    gamZX_lam_xy_i[i] = gamXZ_lam_xy_i[i];

                    gamZX_lam_xz_i[i] = gamXZ_lam_xz_i[i];

                    gamZX_lam_yz_i[i] = gamXZ_lam_yz_i[i];

                    // gamZY

                    gamZY_sig_x_r[i] = gamYZ_sig_x_r[i];

                    gamZY_sig_y_r[i] = gamYZ_sig_y_r[i];

                    gamZY_sig_z_r[i] = gamYZ_sig_z_r[i];

                    gamZY_lam_xy_r[i] = gamYZ_lam_xy_r[i];

                    gamZY_lam_xz_r[i] = gamYZ_lam_xz_r[i];

                    gamZY_lam_yz_r[i] = gamYZ_lam_yz_r[i];

                    // gamZY imag

                    gamZY_sig_x_i[i] = gamYZ_sig_x_i[i];

                    gamZY_sig_y_i[i] = gamYZ_sig_y_i[i];

                    gamZY_sig_z_i[i] = gamYZ_sig_z_i[i];

                    gamZY_lam_xy_i[i] = gamYZ_lam_xy_i[i];

                    gamZY_lam_xz_i[i] = gamYZ_lam_xz_i[i];

                    gamZY_lam_yz_i[i] = gamYZ_lam_yz_i[i];

                    // gamZZ

                    jj_r = 2.0 * (gradzw_kx_r[i] * gradzw_kx_r[i] + gradzw_ky_r[i] * gradzw_ky_r[i] + gradzw_kz_r[i] * gradzw_kz_r[i])

                           - 2.0 * (gradzw_kx_i[i] * gradzw_kx_i[i] + gradzw_ky_i[i] * gradzw_ky_i[i] + gradzw_kz_i[i] * gradzw_kz_i[i]);

                    gamZZ_sig_x_r[i] = 4.0 * (gradzw_kx_r[i] * gradzw_kx_r[i] - gradzw_kx_i[i] * gradzw_kx_i[i]) + jj_r;

                    gamZZ_sig_y_r[i] = 4.0 * (gradzw_ky_r[i] * gradzw_ky_r[i] - gradzw_ky_i[i] * gradzw_ky_i[i]) + jj_r;

                    gamZZ_sig_z_r[i] = 4.0 * (gradzw_kz_r[i] * gradzw_kz_r[i] - gradzw_kz_i[i] * gradzw_kz_i[i]) + jj_r;

                    gamZZ_lam_xy_r[i] = gradzw_kx_r[i] * gradzw_ky_r[i] - gradzw_kx_i[i] * gradzw_ky_i[i]

                                           + gradzw_ky_r[i] * gradzw_kx_r[i] - gradzw_ky_i[i] * gradzw_kx_i[i];

                    gamZZ_lam_xz_r[i] = gradzw_kx_r[i] * gradzw_kz_r[i] - gradzw_kx_i[i] * gradzw_kz_i[i]

                                           + gradzw_kz_r[i] * gradzw_kx_r[i] - gradzw_kz_i[i] * gradzw_kx_i[i];

                    gamZZ_lam_yz_r[i] = gradzw_ky_r[i] * gradzw_kz_r[i] - gradzw_ky_i[i] * gradzw_kz_i[i]

                                           + gradzw_kz_r[i] * gradzw_ky_r[i] - gradzw_kz_i[i] * gradzw_ky_i[i];

                    // gamZZ imag

                    jj_i = 2.0 * (gradzw_kx_r[i] * gradzw_kx_i[i] + gradzw_ky_r[i] * gradzw_ky_i[i] + gradzw_kz_r[i] * gradzw_kz_i[i]

                                  + gradzw_kx_i[i] * gradzw_kx_r[i] + gradzw_ky_i[i] * gradzw_ky_r[i] + gradzw_kz_i[i] * gradzw_kz_r[i]);

                    gamZZ_sig_x_i[i] = 4.0 * (gradzw_kx_r[i] * gradzw_kx_i[i] + gradzw_kx_i[i] * gradzw_kx_r[i]) + jj_i;

                    gamZZ_sig_y_i[i] = 4.0 * (gradzw_ky_r[i] * gradzw_ky_i[i] + gradzw_ky_i[i] * gradzw_ky_r[i]) + jj_i;

                    gamZZ_sig_z_i[i] = 4.0 * (gradzw_kz_r[i] * gradzw_kz_i[i] + gradzw_kz_i[i] * gradzw_kz_r[i]) + jj_i;

                    gamZZ_lam_xy_i[i] = gradzw_kx_r[i] * gradzw_ky_i[i] + gradzw_kx_i[i] * gradzw_ky_r[i]

                                           + gradzw_ky_r[i] * gradzw_kx_i[i] + gradzw_ky_i[i] * gradzw_kx_r[i];

                    gamZZ_lam_xz_i[i] = gradzw_kx_r[i] * gradzw_kz_i[i] + gradzw_kx_i[i] * gradzw_kz_r[i]

                                           + gradzw_kz_r[i] * gradzw_kx_i[i] + gradzw_kz_i[i] * gradzw_kx_r[i];

                    gamZZ_lam_yz_i[i] = gradzw_ky_r[i] * gradzw_kz_i[i] + gradzw_ky_i[i] * gradzw_kz_r[i]

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

                auto rho_sig_x_r = gam(6 * j);

                auto rho_sig_y_r = gam(6 * j + 1);

                auto rho_sig_z_r = gam(6 * j + 2);

                auto rho_lam_xy_r = gam(6 * j + 3);

                auto rho_lam_xz_r = gam(6 * j + 4);

                auto rho_lam_yz_r = gam(6 * j + 5);

                // rxw1rhow part

                auto gamX_sig_x_r = gamX(6 * j);

                auto gamX_sig_y_r = gamX(6 * j + 1);

                auto gamX_sig_z_r = gamX(6 * j + 2);

                auto gamX_lam_xy_r = gamX(6 * j + 3);

                auto gamX_lam_xz_r = gamX(6 * j + 4);

                auto gamX_lam_yz_r = gamX(6 * j + 5);

                // ryw1rhow part

                auto gamY_sig_x_r = gamY(6 * j);

                auto gamY_sig_y_r = gamY(6 * j + 1);

                auto gamY_sig_z_r = gamY(6 * j + 2);

                auto gamY_lam_xy_r = gamY(6 * j + 3);

                auto gamY_lam_xz_r = gamY(6 * j + 4);

                auto gamY_lam_yz_r = gamY(6 * j + 5);

                // rzw1rhow part

                auto gamZ_sig_x_r = gamZ(6 * j);

                auto gamZ_sig_y_r = gamZ(6 * j + 1);

                auto gamZ_sig_z_r = gamZ(6 * j + 2);

                auto gamZ_lam_xy_r = gamZ(6 * j + 3);

                auto gamZ_lam_xz_r = gamZ(6 * j + 4);

                auto gamZ_lam_yz_r = gamZ(6 * j + 5);

                // gamXX part

                auto gamXX_sig_x_r = gamXX(6 * j);

                auto gamXX_sig_y_r = gamXX(6 * j + 1);

                auto gamXX_sig_z_r = gamXX(6 * j + 2);

                auto gamXX_lam_xy_r = gamXX(6 * j + 3);

                auto gamXX_lam_xz_r = gamXX(6 * j + 4);

                auto gamXX_lam_yz_r = gamXX(6 * j + 5);

                // gamXY part

                auto gamXY_sig_x_r = gamXY(6 * j);

                auto gamXY_sig_y_r = gamXY(6 * j + 1);

                auto gamXY_sig_z_r = gamXY(6 * j + 2);

                auto gamXY_lam_xy_r = gamXY(6 * j + 3);

                auto gamXY_lam_xz_r = gamXY(6 * j + 4);

                auto gamXY_lam_yz_r = gamXY(6 * j + 5);

                // gamXZ part

                auto gamXZ_sig_x_r = gamXZ(6 * j);

                auto gamXZ_sig_y_r = gamXZ(6 * j + 1);

                auto gamXZ_sig_z_r = gamXZ(6 * j + 2);

                auto gamXZ_lam_xy_r = gamXZ(6 * j + 3);

                auto gamXZ_lam_xz_r = gamXZ(6 * j + 4);

                auto gamXZ_lam_yz_r = gamXZ(6 * j + 5);

                // gamYX part

                auto gamYX_sig_x_r = gamYX(6 * j);

                auto gamYX_sig_y_r = gamYX(6 * j + 1);

                auto gamYX_sig_z_r = gamYX(6 * j + 2);

                auto gamYX_lam_xy_r = gamYX(6 * j + 3);

                auto gamYX_lam_xz_r = gamYX(6 * j + 4);

                auto gamYX_lam_yz_r = gamYX(6 * j + 5);

                // gamYY part

                auto gamYY_sig_x_r = gamYY(6 * j);

                auto gamYY_sig_y_r = gamYY(6 * j + 1);

                auto gamYY_sig_z_r = gamYY(6 * j + 2);

                auto gamYY_lam_xy_r = gamYY(6 * j + 3);

                auto gamYY_lam_xz_r = gamYY(6 * j + 4);

                auto gamYY_lam_yz_r = gamYY(6 * j + 5);

                // gamYZ part

                auto gamYZ_sig_x_r = gamYZ(6 * j);

                auto gamYZ_sig_y_r = gamYZ(6 * j + 1);

                auto gamYZ_sig_z_r = gamYZ(6 * j + 2);

                auto gamYZ_lam_xy_r = gamYZ(6 * j + 3);

                auto gamYZ_lam_xz_r = gamYZ(6 * j + 4);

                auto gamYZ_lam_yz_r = gamYZ(6 * j + 5);

                // gamZX part

                auto gamZX_sig_x_r = gamZX(6 * j);

                auto gamZX_sig_y_r = gamZX(6 * j + 1);

                auto gamZX_sig_z_r = gamZX(6 * j + 2);

                auto gamZX_lam_xy_r = gamZX(6 * j + 3);

                auto gamZX_lam_xz_r = gamZX(6 * j + 4);

                auto gamZX_lam_yz_r = gamZX(6 * j + 5);

                // gamZY part

                auto gamZY_sig_x_r = gamZY(6 * j);

                auto gamZY_sig_y_r = gamZY(6 * j + 1);

                auto gamZY_sig_z_r = gamZY(6 * j + 2);

                auto gamZY_lam_xy_r = gamZY(6 * j + 3);

                auto gamZY_lam_xz_r = gamZY(6 * j + 4);

                auto gamZY_lam_yz_r = gamZY(6 * j + 5);

                // gamZZ part

                auto gamZZ_sig_x_r = gamZZ(6 * j);

                auto gamZZ_sig_y_r = gamZZ(6 * j + 1);

                auto gamZZ_sig_z_r = gamZZ(6 * j + 2);

                auto gamZZ_lam_xy_r = gamZZ(6 * j + 3);

                auto gamZZ_lam_xz_r = gamZZ(6 * j + 4);

                auto gamZZ_lam_yz_r = gamZZ(6 * j + 5);

                // First-order densities

                // Rho terms

                auto rhow_kx_r = rwDensityGrid.alphaDensity(6 * j);

                auto rhow_kx_i = rwDensityGrid.alphaDensity(6 * j + 1);

                auto rhow_ky_r = rwDensityGrid.alphaDensity(6 * j + 2);

                auto rhow_ky_i = rwDensityGrid.alphaDensity(6 * j + 3);

                auto rhow_kz_r = rwDensityGrid.alphaDensity(6 * j + 4);

                auto rhow_kz_i = rwDensityGrid.alphaDensity(6 * j + 5);

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

                    // gamX real

                    jj_r = 4.0 * (gradxw_kx_r[i] * rhow_kx_r[i] + gradxw_ky_r[i] * rhow_ky_r[i] + gradxw_kz_r[i] * rhow_kz_r[i])

                           - 4.0 * (gradxw_kx_i[i] * rhow_kx_i[i] + gradxw_ky_i[i] * rhow_ky_i[i] + gradxw_kz_i[i] * rhow_kz_i[i]);

                    gamX_sig_x_r[i] = 8.0 * (gradxw_kx_r[i] * rhow_kx_r[i] - gradxw_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    gamX_sig_y_r[i] = 8.0 * (gradxw_ky_r[i] * rhow_ky_r[i] - gradxw_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    gamX_sig_z_r[i] = 8.0 * (gradxw_kz_r[i] * rhow_kz_r[i] - gradxw_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    gamX_lam_xy_r[i] = 2.0 * (gradxw_kx_r[i] * rhow_ky_r[i] - gradxw_kx_i[i] * rhow_ky_i[i]

                                                   + gradxw_ky_r[i] * rhow_kx_r[i] - gradxw_ky_i[i] * rhow_kx_i[i]);

                    gamX_lam_xz_r[i] = 2.0 * (gradxw_kx_r[i] * rhow_kz_r[i] - gradxw_kx_i[i] * rhow_kz_i[i]

                                                   + gradxw_kz_r[i] * rhow_kx_r[i] - gradxw_kz_i[i] * rhow_kx_i[i]);

                    gamX_lam_yz_r[i] = 2.0 * (gradxw_ky_r[i] * rhow_kz_r[i] - gradxw_ky_i[i] * rhow_kz_i[i]

                                                   + gradxw_kz_r[i] * rhow_ky_r[i] - gradxw_kz_i[i] * rhow_ky_i[i]);

                    // gamX imag

                    // gamY real

                    jj_r = 4.0 * (gradyw_kx_r[i] * rhow_kx_r[i] + gradyw_ky_r[i] * rhow_ky_r[i] + gradyw_kz_r[i] * rhow_kz_r[i])

                           - 4.0 * (gradyw_kx_i[i] * rhow_kx_i[i] + gradyw_ky_i[i] * rhow_ky_i[i] + gradyw_kz_i[i] * rhow_kz_i[i]);

                    gamY_sig_x_r[i] = 8.0 * (gradyw_kx_r[i] * rhow_kx_r[i] - gradyw_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    gamY_sig_y_r[i] = 8.0 * (gradyw_ky_r[i] * rhow_ky_r[i] - gradyw_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    gamY_sig_z_r[i] = 8.0 * (gradyw_kz_r[i] * rhow_kz_r[i] - gradyw_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    gamY_lam_xy_r[i] = 2.0 * (gradyw_kx_r[i] * rhow_ky_r[i] - gradyw_kx_i[i] * rhow_ky_i[i]

                                                   + gradyw_ky_r[i] * rhow_kx_r[i] - gradyw_ky_i[i] * rhow_kx_i[i]);

                    gamY_lam_xz_r[i] = 2.0 * (gradyw_kx_r[i] * rhow_kz_r[i] - gradyw_kx_i[i] * rhow_kz_i[i]

                                                   + gradyw_kz_r[i] * rhow_kx_r[i] - gradyw_kz_i[i] * rhow_kx_i[i]);

                    gamY_lam_yz_r[i] = 2.0 * (gradyw_ky_r[i] * rhow_kz_r[i] - gradyw_ky_i[i] * rhow_kz_i[i]

                                                   + gradyw_kz_r[i] * rhow_ky_r[i] - gradyw_kz_i[i] * rhow_ky_i[i]);

                    // gamY imag

                    // gamZ real

                    jj_r = 4.0 * (gradzw_kx_r[i] * rhow_kx_r[i] + gradzw_ky_r[i] * rhow_ky_r[i] + gradzw_kz_r[i] * rhow_kz_r[i])

                           - 4.0 * (gradzw_kx_i[i] * rhow_kx_i[i] + gradzw_ky_i[i] * rhow_ky_i[i] + gradzw_kz_i[i] * rhow_kz_i[i]);

                    gamZ_sig_x_r[i] = 8.0 * (gradzw_kx_r[i] * rhow_kx_r[i] - gradzw_kx_i[i] * rhow_kx_i[i]) + jj_r;

                    gamZ_sig_y_r[i] = 8.0 * (gradzw_ky_r[i] * rhow_ky_r[i] - gradzw_ky_i[i] * rhow_ky_i[i]) + jj_r;

                    gamZ_sig_z_r[i] = 8.0 * (gradzw_kz_r[i] * rhow_kz_r[i] - gradzw_kz_i[i] * rhow_kz_i[i]) + jj_r;

                    gamZ_lam_xy_r[i] = 2.0 * (gradzw_kx_r[i] * rhow_ky_r[i] - gradzw_kx_i[i] * rhow_ky_i[i]

                                                   + gradzw_ky_r[i] * rhow_kx_r[i] - gradzw_ky_i[i] * rhow_kx_i[i]);

                    gamZ_lam_xz_r[i] = 2.0 * (gradzw_kx_r[i] * rhow_kz_r[i] - gradzw_kx_i[i] * rhow_kz_i[i]

                                                   + gradzw_kz_r[i] * rhow_kx_r[i] - gradzw_kz_i[i] * rhow_kx_i[i]);

                    gamZ_lam_yz_r[i] = 2.0 * (gradzw_ky_r[i] * rhow_kz_r[i] - gradzw_ky_i[i] * rhow_kz_i[i]

                                                   + gradzw_kz_r[i] * rhow_ky_r[i] - gradzw_kz_i[i] * rhow_ky_i[i]);

                    // gamZ imag

                    // gamXX

                    jj_r = 2.0 * (gradxw_kx_r[i] * gradxw_kx_r[i] + gradxw_ky_r[i] * gradxw_ky_r[i] + gradxw_kz_r[i] * gradxw_kz_r[i])

                           - 2.0 * (gradxw_kx_i[i] * gradxw_kx_i[i] + gradxw_ky_i[i] * gradxw_ky_i[i] + gradxw_kz_i[i] * gradxw_kz_i[i]);

                    gamXX_sig_x_r[i] = 4.0 * (gradxw_kx_r[i] * gradxw_kx_r[i] - gradxw_kx_i[i] * gradxw_kx_i[i]) + jj_r;

                    gamXX_sig_y_r[i] = 4.0 * (gradxw_ky_r[i] * gradxw_ky_r[i] - gradxw_ky_i[i] * gradxw_ky_i[i]) + jj_r;

                    gamXX_sig_z_r[i] = 4.0 * (gradxw_kz_r[i] * gradxw_kz_r[i] - gradxw_kz_i[i] * gradxw_kz_i[i]) + jj_r;

                    gamXX_lam_xy_r[i] = gradxw_kx_r[i] * gradxw_ky_r[i] - gradxw_kx_i[i] * gradxw_ky_i[i]

                                           + gradxw_ky_r[i] * gradxw_kx_r[i] - gradxw_ky_i[i] * gradxw_kx_i[i];

                    gamXX_lam_xz_r[i] = gradxw_kx_r[i] * gradxw_kz_r[i] - gradxw_kx_i[i] * gradxw_kz_i[i]

                                           + gradxw_kz_r[i] * gradxw_kx_r[i] - gradxw_kz_i[i] * gradxw_kx_i[i];

                    gamXX_lam_yz_r[i] = gradxw_ky_r[i] * gradxw_kz_r[i] - gradxw_ky_i[i] * gradxw_kz_i[i]

                                           + gradxw_kz_r[i] * gradxw_ky_r[i] - gradxw_kz_i[i] * gradxw_ky_i[i];

                    // gamXX imag

                    // gamXY

                    jj_r = 2.0 * (gradxw_kx_r[i] * gradyw_kx_r[i] + gradxw_ky_r[i] * gradyw_ky_r[i] + gradxw_kz_r[i] * gradyw_kz_r[i])

                           - 2.0 * (gradxw_kx_i[i] * gradyw_kx_i[i] + gradxw_ky_i[i] * gradyw_ky_i[i] + gradxw_kz_i[i] * gradyw_kz_i[i]);

                    gamXY_sig_x_r[i] = 4.0 * (gradxw_kx_r[i] * gradyw_kx_r[i] - gradxw_kx_i[i] * gradyw_kx_i[i]) + jj_r;

                    gamXY_sig_y_r[i] = 4.0 * (gradxw_ky_r[i] * gradyw_ky_r[i] - gradxw_ky_i[i] * gradyw_ky_i[i]) + jj_r;

                    gamXY_sig_z_r[i] = 4.0 * (gradxw_kz_r[i] * gradyw_kz_r[i] - gradxw_kz_i[i] * gradyw_kz_i[i]) + jj_r;

                    gamXY_lam_xy_r[i] = gradxw_kx_r[i] * gradyw_ky_r[i] - gradxw_kx_i[i] * gradyw_ky_i[i]

                                           + gradxw_ky_r[i] * gradyw_kx_r[i] - gradxw_ky_i[i] * gradyw_kx_i[i];

                    gamXY_lam_xz_r[i] = gradxw_kx_r[i] * gradyw_kz_r[i] - gradxw_kx_i[i] * gradyw_kz_i[i]

                                           + gradxw_kz_r[i] * gradyw_kx_r[i] - gradxw_kz_i[i] * gradyw_kx_i[i];

                    gamXY_lam_yz_r[i] = gradxw_ky_r[i] * gradyw_kz_r[i] - gradxw_ky_i[i] * gradyw_kz_i[i]

                                           + gradxw_kz_r[i] * gradyw_ky_r[i] - gradxw_kz_i[i] * gradyw_ky_i[i];

                    // gamXY imag

                    // gamXZ

                    jj_r = 2.0 * (gradxw_kx_r[i] * gradzw_kx_r[i] + gradxw_ky_r[i] * gradzw_ky_r[i] + gradxw_kz_r[i] * gradzw_kz_r[i])

                           - 2.0 * (gradxw_kx_i[i] * gradzw_kx_i[i] + gradxw_ky_i[i] * gradzw_ky_i[i] + gradxw_kz_i[i] * gradzw_kz_i[i]);

                    gamXZ_sig_x_r[i] = 4.0 * (gradxw_kx_r[i] * gradzw_kx_r[i] - gradxw_kx_i[i] * gradzw_kx_i[i]) + jj_r;

                    gamXZ_sig_y_r[i] = 4.0 * (gradxw_ky_r[i] * gradzw_ky_r[i] - gradxw_ky_i[i] * gradzw_ky_i[i]) + jj_r;

                    gamXZ_sig_z_r[i] = 4.0 * (gradxw_kz_r[i] * gradzw_kz_r[i] - gradxw_kz_i[i] * gradzw_kz_i[i]) + jj_r;

                    gamXZ_lam_xy_r[i] = gradxw_kx_r[i] * gradzw_ky_r[i] - gradxw_kx_i[i] * gradzw_ky_i[i]

                                           + gradxw_ky_r[i] * gradzw_kx_r[i] - gradxw_ky_i[i] * gradzw_kx_i[i];

                    gamXZ_lam_xz_r[i] = gradxw_kx_r[i] * gradzw_kz_r[i] - gradxw_kx_i[i] * gradzw_kz_i[i]

                                           + gradxw_kz_r[i] * gradzw_kx_r[i] - gradxw_kz_i[i] * gradzw_kx_i[i];

                    gamXZ_lam_yz_r[i] = gradxw_ky_r[i] * gradzw_kz_r[i] - gradxw_ky_i[i] * gradzw_kz_i[i]

                                           + gradxw_kz_r[i] * gradzw_ky_r[i] - gradxw_kz_i[i] * gradzw_ky_i[i];

                    // gamXZ imag

                    // gamYX

                    gamYX_sig_x_r[i] = gamXY_sig_x_r[i];

                    gamYX_sig_y_r[i] = gamXY_sig_y_r[i];

                    gamYX_sig_z_r[i] = gamXY_sig_z_r[i];

                    gamYX_lam_xy_r[i] = gamXY_lam_xy_r[i];

                    gamYX_lam_xz_r[i] = gamXY_lam_xz_r[i];

                    gamYX_lam_yz_r[i] = gamXY_lam_yz_r[i];

                    // gamYX imag

                    // gamYY

                    jj_r = 2.0 * (gradyw_kx_r[i] * gradyw_kx_r[i] + gradyw_ky_r[i] * gradyw_ky_r[i] + gradyw_kz_r[i] * gradyw_kz_r[i])

                           - 2.0 * (gradyw_kx_i[i] * gradyw_kx_i[i] + gradyw_ky_i[i] * gradyw_ky_i[i] + gradyw_kz_i[i] * gradyw_kz_i[i]);

                    gamYY_sig_x_r[i] = 4.0 * (gradyw_kx_r[i] * gradyw_kx_r[i] - gradyw_kx_i[i] * gradyw_kx_i[i]) + jj_r;

                    gamYY_sig_y_r[i] = 4.0 * (gradyw_ky_r[i] * gradyw_ky_r[i] - gradyw_ky_i[i] * gradyw_ky_i[i]) + jj_r;

                    gamYY_sig_z_r[i] = 4.0 * (gradyw_kz_r[i] * gradyw_kz_r[i] - gradyw_kz_i[i] * gradyw_kz_i[i]) + jj_r;

                    gamYY_lam_xy_r[i] = gradyw_kx_r[i] * gradyw_ky_r[i] - gradyw_kx_i[i] * gradyw_ky_i[i]

                                           + gradyw_ky_r[i] * gradyw_kx_r[i] - gradyw_ky_i[i] * gradyw_kx_i[i];

                    gamYY_lam_xz_r[i] = gradyw_kx_r[i] * gradyw_kz_r[i] - gradyw_kx_i[i] * gradyw_kz_i[i]

                                           + gradyw_kz_r[i] * gradyw_kx_r[i] - gradyw_kz_i[i] * gradyw_kx_i[i];

                    gamYY_lam_yz_r[i] = gradyw_ky_r[i] * gradyw_kz_r[i] - gradyw_ky_i[i] * gradyw_kz_i[i]

                                           + gradyw_kz_r[i] * gradyw_ky_r[i] - gradyw_kz_i[i] * gradyw_ky_i[i];

                    // gamYY imag

                    // gamYZ

                    jj_r = 2.0 * (gradyw_kx_r[i] * gradzw_kx_r[i] + gradyw_ky_r[i] * gradzw_ky_r[i] + gradyw_kz_r[i] * gradzw_kz_r[i])

                           - 2.0 * (gradyw_kx_i[i] * gradzw_kx_i[i] + gradyw_ky_i[i] * gradzw_ky_i[i] + gradyw_kz_i[i] * gradzw_kz_i[i]);

                    gamYZ_sig_x_r[i] = 4.0 * (gradyw_kx_r[i] * gradzw_kx_r[i] - gradyw_kx_i[i] * gradzw_kx_i[i]) + jj_r;

                    gamYZ_sig_y_r[i] = 4.0 * (gradyw_ky_r[i] * gradzw_ky_r[i] - gradyw_ky_i[i] * gradzw_ky_i[i]) + jj_r;

                    gamYZ_sig_z_r[i] = 4.0 * (gradyw_kz_r[i] * gradzw_kz_r[i] - gradyw_kz_i[i] * gradzw_kz_i[i]) + jj_r;

                    gamYZ_lam_xy_r[i] = gradyw_kx_r[i] * gradzw_ky_r[i] - gradyw_kx_i[i] * gradzw_ky_i[i]

                                           + gradyw_ky_r[i] * gradzw_kx_r[i] - gradyw_ky_i[i] * gradzw_kx_i[i];

                    gamYZ_lam_xz_r[i] = gradyw_kx_r[i] * gradzw_kz_r[i] - gradyw_kx_i[i] * gradzw_kz_i[i]

                                           + gradyw_kz_r[i] * gradzw_kx_r[i] - gradyw_kz_i[i] * gradzw_kx_i[i];

                    gamYZ_lam_yz_r[i] = gradyw_ky_r[i] * gradzw_kz_r[i] - gradyw_ky_i[i] * gradzw_kz_i[i]

                                           + gradyw_kz_r[i] * gradzw_ky_r[i] - gradyw_kz_i[i] * gradzw_ky_i[i];

                    // gamYZ imag

                    // gamZX

                    gamZX_sig_x_r[i] = gamXZ_sig_x_r[i];

                    gamZX_sig_y_r[i] = gamXZ_sig_y_r[i];

                    gamZX_sig_z_r[i] = gamXZ_sig_z_r[i];

                    gamZX_lam_xy_r[i] = gamXZ_lam_xy_r[i];

                    gamZX_lam_xz_r[i] = gamXZ_lam_xz_r[i];

                    gamZX_lam_yz_r[i] = gamXZ_lam_yz_r[i];

                    // gamZX imag

                    // gamZY

                    gamZY_sig_x_r[i] = gamYZ_sig_x_r[i];

                    gamZY_sig_y_r[i] = gamYZ_sig_y_r[i];

                    gamZY_sig_z_r[i] = gamYZ_sig_z_r[i];

                    gamZY_lam_xy_r[i] = gamYZ_lam_xy_r[i];

                    gamZY_lam_xz_r[i] = gamYZ_lam_xz_r[i];

                    gamZY_lam_yz_r[i] = gamYZ_lam_yz_r[i];

                    // gamZY imag

                    // gamZZ

                    jj_r = 2.0 * (gradzw_kx_r[i] * gradzw_kx_r[i] + gradzw_ky_r[i] * gradzw_ky_r[i] + gradzw_kz_r[i] * gradzw_kz_r[i])

                           - 2.0 * (gradzw_kx_i[i] * gradzw_kx_i[i] + gradzw_ky_i[i] * gradzw_ky_i[i] + gradzw_kz_i[i] * gradzw_kz_i[i]);

                    gamZZ_sig_x_r[i] = 4.0 * (gradzw_kx_r[i] * gradzw_kx_r[i] - gradzw_kx_i[i] * gradzw_kx_i[i]) + jj_r;

                    gamZZ_sig_y_r[i] = 4.0 * (gradzw_ky_r[i] * gradzw_ky_r[i] - gradzw_ky_i[i] * gradzw_ky_i[i]) + jj_r;

                    gamZZ_sig_z_r[i] = 4.0 * (gradzw_kz_r[i] * gradzw_kz_r[i] - gradzw_kz_i[i] * gradzw_kz_i[i]) + jj_r;

                    gamZZ_lam_xy_r[i] = gradzw_kx_r[i] * gradzw_ky_r[i] - gradzw_kx_i[i] * gradzw_ky_i[i]

                                           + gradzw_ky_r[i] * gradzw_kx_r[i] - gradzw_ky_i[i] * gradzw_kx_i[i];

                    gamZZ_lam_xz_r[i] = gradzw_kx_r[i] * gradzw_kz_r[i] - gradzw_kx_i[i] * gradzw_kz_i[i]

                                           + gradzw_kz_r[i] * gradzw_kx_r[i] - gradzw_kz_i[i] * gradzw_kx_i[i];

                    gamZZ_lam_yz_r[i] = gradzw_ky_r[i] * gradzw_kz_r[i] - gradzw_ky_i[i] * gradzw_kz_i[i]

                                           + gradzw_kz_r[i] * gradzw_ky_r[i] - gradzw_kz_i[i] * gradzw_ky_i[i];

                    // gamZZ imag
                }
            }
        }
        if (fstr::upcase(quadMode) == "QRF")
        {
            for (int32_t j = 0; j < numdens / 2; j++)
            {
               // Perturbed densities

                auto rhoB_r = rwDensityGrid.alphaDensity(4 * j);
                auto gradB_x_r = rwDensityGrid.alphaDensityGradientX(4 * j);
                auto gradB_y_r = rwDensityGrid.alphaDensityGradientY(4 * j);
                auto gradB_z_r = rwDensityGrid.alphaDensityGradientZ(4 * j);

                auto rhoB_i = rwDensityGrid.alphaDensity(4 * j + 1);
                auto gradB_x_i = rwDensityGrid.alphaDensityGradientX(4 * j + 1);
                auto gradB_y_i = rwDensityGrid.alphaDensityGradientY(4 * j + 1);
                auto gradB_z_i = rwDensityGrid.alphaDensityGradientZ(4 * j + 1);

                auto rhoC_r = rwDensityGrid.alphaDensity(4 * j + 2);
                auto gradC_x_r = rwDensityGrid.alphaDensityGradientX(4 * j + 2);
                auto gradC_y_r = rwDensityGrid.alphaDensityGradientY(4 * j + 2);
                auto gradC_z_r = rwDensityGrid.alphaDensityGradientZ(4 * j + 2);

                auto rhoC_i = rwDensityGrid.alphaDensity(4 * j + 3);
                auto gradC_x_i = rwDensityGrid.alphaDensityGradientX(4 * j + 3);
                auto gradC_y_i = rwDensityGrid.alphaDensityGradientY(4 * j + 3);
                auto gradC_z_i = rwDensityGrid.alphaDensityGradientZ(4 * j + 3);

                // Perturbed density products to be stored

                auto gam_r = gam(2 * j);
                auto gam_i = gam(2 * j + 1);

                auto gam_x_r = gamX(2 * j);
                auto gam_x_i = gamX(2 * j + 1);
                auto gam_y_r = gamY(2 * j);
                auto gam_y_i = gamY(2 * j + 1);
                auto gam_z_r = gamZ(2 * j);
                auto gam_z_i = gamZ(2 * j + 1);

                auto gam_xx_r = gamXX(2 * j);
                auto gam_xx_i = gamXX(2 * j + 1);
                auto gam_xy_r = gamXY(2 * j);
                auto gam_xy_i = gamXY(2 * j + 1);
                auto gam_xz_r = gamXZ(2 * j);
                auto gam_xz_i = gamXZ(2 * j + 1);
                auto gam_yx_r = gamYX(2 * j);
                auto gam_yx_i = gamYX(2 * j + 1);
                auto gam_yy_r = gamYY(2 * j);
                auto gam_yy_i = gamYY(2 * j + 1);
                auto gam_yz_r = gamYZ(2 * j);
                auto gam_yz_i = gamYZ(2 * j + 1);
                auto gam_zx_r = gamZX(2 * j);
                auto gam_zx_i = gamZX(2 * j + 1);
                auto gam_zy_r = gamZY(2 * j);
                auto gam_zy_i = gamZY(2 * j + 1);
                auto gam_zz_r = gamZZ(2 * j);
                auto gam_zz_i = gamZZ(2 * j + 1);

                for (int32_t i = 0; i < npoints; i++)
                {
                    gam_r[i] += prod2_r(rhoB_r[i],rhoB_i[i],rhoC_r[i],rhoC_i[i])
                              + prod2_r(rhoC_r[i],rhoC_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_i[i] += prod2_i(rhoB_r[i],rhoB_i[i],rhoC_r[i],rhoC_i[i])
                              + prod2_i(rhoC_r[i],rhoC_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_x_r[i] += 2.0 * prod2_r(gradB_x_r[i],gradB_x_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_r(gradC_x_r[i],gradC_x_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_y_r[i] += 2.0 * prod2_r(gradB_y_r[i],gradB_y_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_r(gradC_y_r[i],gradC_y_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_z_r[i] += 2.0 * prod2_r(gradB_z_r[i],gradB_z_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_r(gradC_z_r[i],gradC_z_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_x_i[i] += 2.0 * prod2_i(gradB_x_r[i],gradB_x_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_i(gradC_x_r[i],gradC_x_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_y_i[i] += 2.0 * prod2_i(gradB_y_r[i],gradB_y_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_i(gradC_y_r[i],gradC_y_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_z_i[i] += 2.0 * prod2_i(gradB_z_r[i],gradB_z_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_i(gradC_z_r[i],gradC_z_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_xx_r[i] += prod2_r(gradB_x_r[i],gradB_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                +  prod2_r(gradC_x_r[i],gradC_x_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_xy_r[i] += prod2_r(gradB_x_r[i],gradB_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                +  prod2_r(gradC_x_r[i],gradC_x_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_xz_r[i] += prod2_r(gradB_x_r[i],gradB_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                +  prod2_r(gradC_x_r[i],gradC_x_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_yx_r[i] += prod2_r(gradB_y_r[i],gradB_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                +  prod2_r(gradC_y_r[i],gradC_y_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_yy_r[i] += prod2_r(gradB_y_r[i],gradB_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                +  prod2_r(gradC_y_r[i],gradC_y_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_yz_r[i] += prod2_r(gradB_y_r[i],gradB_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                +  prod2_r(gradC_y_r[i],gradC_y_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_zx_r[i] += prod2_r(gradB_z_r[i],gradB_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                +  prod2_r(gradC_z_r[i],gradC_z_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_zy_r[i] += prod2_r(gradB_z_r[i],gradB_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                +  prod2_r(gradC_z_r[i],gradC_z_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_zz_r[i] += prod2_r(gradB_z_r[i],gradB_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                +  prod2_r(gradC_z_r[i],gradC_z_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_xx_i[i] += prod2_i(gradB_x_r[i],gradB_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                +  prod2_i(gradC_x_r[i],gradC_x_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_xy_i[i] += prod2_i(gradB_x_r[i],gradB_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                +  prod2_i(gradC_x_r[i],gradC_x_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_xz_i[i] += prod2_i(gradB_x_r[i],gradB_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                +  prod2_i(gradC_x_r[i],gradC_x_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_yx_i[i] += prod2_i(gradB_y_r[i],gradB_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                +  prod2_i(gradC_y_r[i],gradC_y_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_yy_i[i] += prod2_i(gradB_y_r[i],gradB_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                +  prod2_i(gradC_y_r[i],gradC_y_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_yz_i[i] += prod2_i(gradB_y_r[i],gradB_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                +  prod2_i(gradC_y_r[i],gradC_y_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_zx_i[i] += prod2_i(gradB_z_r[i],gradB_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                +  prod2_i(gradC_z_r[i],gradC_z_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_zy_i[i] += prod2_i(gradB_z_r[i],gradB_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                +  prod2_i(gradC_z_r[i],gradC_z_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_zz_i[i] += prod2_i(gradB_z_r[i],gradB_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                +  prod2_i(gradC_z_r[i],gradC_z_i[i],gradB_z_r[i],gradB_z_i[i]);

                }
            }
        }
    }
    if (xcFuncType == xcfun::mgga)
    {
        if (fstr::upcase(quadMode) == "QRF")
        {
            for (int32_t j = 0; j < numdens / 2; j++)
            {
                // Perturbed densities
                
                auto rhoB_r = rwDensityGrid.alphaDensity(4 * j);
                auto tauB_r = rwDensityGrid.alphaDensitytau(4 * j);
                auto laplB_r = rwDensityGrid.alphaDensitylapl(4 * j);
                auto gradB_x_r = rwDensityGrid.alphaDensityGradientX(4 * j);
                auto gradB_y_r = rwDensityGrid.alphaDensityGradientY(4 * j);
                auto gradB_z_r = rwDensityGrid.alphaDensityGradientZ(4 * j);

                auto rhoB_i = rwDensityGrid.alphaDensity(4 * j + 1);
                auto tauB_i = rwDensityGrid.alphaDensitytau(4 * j + 1);
                auto laplB_i = rwDensityGrid.alphaDensitylapl(4 * j + 1);
                auto gradB_x_i = rwDensityGrid.alphaDensityGradientX(4 * j + 1);
                auto gradB_y_i = rwDensityGrid.alphaDensityGradientY(4 * j + 1);
                auto gradB_z_i = rwDensityGrid.alphaDensityGradientZ(4 * j + 1);

                auto rhoC_r = rwDensityGrid.alphaDensity(4 * j + 2);
                auto tauC_r = rwDensityGrid.alphaDensitytau(4 * j + 2);
                auto laplC_r = rwDensityGrid.alphaDensitylapl(4 * j + 2);
                auto gradC_x_r = rwDensityGrid.alphaDensityGradientX(4 * j + 2);
                auto gradC_y_r = rwDensityGrid.alphaDensityGradientY(4 * j + 2);
                auto gradC_z_r = rwDensityGrid.alphaDensityGradientZ(4 * j + 2);

                auto rhoC_i = rwDensityGrid.alphaDensity(4 * j + 3);
                auto tauC_i = rwDensityGrid.alphaDensitytau(4 * j + 3);
                auto laplC_i = rwDensityGrid.alphaDensitylapl(4 * j + 3);
                auto gradC_x_i = rwDensityGrid.alphaDensityGradientX(4 * j + 3);
                auto gradC_y_i = rwDensityGrid.alphaDensityGradientY(4 * j + 3);
                auto gradC_z_i = rwDensityGrid.alphaDensityGradientZ(4 * j + 3);

                // Perturbed density products to be stored

                auto gam_r = gam(2 * j);
                auto rt_gam_r = rt_gam(2 * j);
                auto rl_gam_r = rl_gam(2 * j);
                auto tt_gam_r = tt_gam(2 * j);
                auto tl_gam_r = tl_gam(2 * j);
                auto ll_gam_r = ll_gam(2 * j);
                
                auto gam_i = gam(2 * j + 1);
                auto rt_gam_i = rt_gam(2 * j + 1);
                auto rl_gam_i = rl_gam(2 * j + 1);
                auto tt_gam_i = tt_gam(2 * j + 1);
                auto tl_gam_i = tl_gam(2 * j + 1);
                auto ll_gam_i = ll_gam(2 * j + 1);

                auto gam_x_r = gamX(2 * j);
                auto gam_y_r = gamY(2 * j);
                auto gam_z_r = gamZ(2 * j);

                auto gam_x_i = gamX(2 * j + 1);
                auto gam_y_i = gamY(2 * j + 1);
                auto gam_z_i = gamZ(2 * j + 1);

                auto st_gam_x_r = st_gamX(2 * j);
                auto st_gam_y_r = st_gamY(2 * j);
                auto st_gam_z_r = st_gamZ(2 * j);

                auto sl_gam_x_r = sl_gamX(2 * j);
                auto sl_gam_y_r = sl_gamY(2 * j);
                auto sl_gam_z_r = sl_gamZ(2 * j);

                auto st_gam_x_i = st_gamX(2 * j + 1);
                auto st_gam_y_i = st_gamY(2 * j + 1);
                auto st_gam_z_i = st_gamZ(2 * j + 1);

                auto sl_gam_x_i = sl_gamX(2 * j + 1);
                auto sl_gam_y_i = sl_gamY(2 * j + 1);
                auto sl_gam_z_i = sl_gamZ(2 * j + 1);

                auto gam_xx_r = gamXX(2 * j);
                auto gam_xx_i = gamXX(2 * j + 1);
                auto gam_xy_r = gamXY(2 * j);
                auto gam_xy_i = gamXY(2 * j + 1);
                auto gam_xz_r = gamXZ(2 * j);
                auto gam_xz_i = gamXZ(2 * j + 1);
                auto gam_yx_r = gamYX(2 * j);
                auto gam_yx_i = gamYX(2 * j + 1);
                auto gam_yy_r = gamYY(2 * j);
                auto gam_yy_i = gamYY(2 * j + 1);
                auto gam_yz_r = gamYZ(2 * j);
                auto gam_yz_i = gamYZ(2 * j + 1);
                auto gam_zx_r = gamZX(2 * j);
                auto gam_zx_i = gamZX(2 * j + 1);
                auto gam_zy_r = gamZY(2 * j);
                auto gam_zy_i = gamZY(2 * j + 1);
                auto gam_zz_r = gamZZ(2 * j);
                auto gam_zz_i = gamZZ(2 * j + 1);

                for (int32_t i = 0; i < npoints; i++)
                {

                    gam_r[i] = prod2_r(rhoB_r[i],rhoB_i[i],rhoC_r[i],rhoC_i[i])
                              + prod2_r(rhoC_r[i],rhoC_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_i[i] = prod2_i(rhoB_r[i],rhoB_i[i],rhoC_r[i],rhoC_i[i])
                              + prod2_i(rhoC_r[i],rhoC_i[i],rhoB_r[i],rhoB_i[i]);

                    ll_gam_r[i] = prod2_r(laplB_r[i],laplB_i[i],laplC_r[i],laplC_i[i])
                                + prod2_r(laplC_r[i],laplC_i[i],laplB_r[i],laplB_i[i]);

                    ll_gam_i[i] = prod2_i(laplB_r[i],laplB_i[i],laplC_r[i],laplC_i[i])
                                 + prod2_i(laplC_r[i],laplC_i[i],laplB_r[i],laplB_i[i]);

                    tt_gam_r[i] = prod2_r(tauB_r[i],tauB_i[i],tauC_r[i],tauC_i[i])
                                 + prod2_r(tauC_r[i],tauC_i[i],tauB_r[i],tauB_i[i]);

                    tt_gam_i[i] = prod2_i(tauB_r[i],tauB_i[i],tauC_r[i],tauC_i[i])
                                 + prod2_i(tauC_r[i],tauC_i[i],tauB_r[i],tauB_i[i]);

                    rt_gam_r[i] = 2.0 * (prod2_r(rhoB_r[i],rhoB_i[i],tauC_r[i],tauC_i[i])
                                       + prod2_r(rhoC_r[i],rhoC_i[i],tauB_r[i],tauB_i[i]));

                    rt_gam_i[i] = 2.0 * (prod2_i(rhoB_r[i],rhoB_i[i],tauC_r[i],tauC_i[i])
                                      + prod2_i(rhoC_r[i],rhoC_i[i],tauB_r[i],tauB_i[i]));

                    rl_gam_r[i] = 2.0 * (prod2_r(rhoB_r[i],rhoB_i[i],laplC_r[i],laplC_i[i])
                                       + prod2_r(rhoC_r[i],rhoC_i[i],laplB_r[i],laplB_i[i]));

                    rl_gam_i[i] = 2.0 * (prod2_i(rhoB_r[i],rhoB_i[i],laplC_r[i],laplC_i[i])
                                       + prod2_i(rhoC_r[i],rhoC_i[i],laplB_r[i],laplB_i[i]));

                    tl_gam_r[i] = 2.0 * (prod2_r(tauB_r[i],tauB_i[i],laplC_r[i],laplC_i[i])
                                      + prod2_r(tauC_r[i],tauC_i[i],laplB_r[i],laplB_i[i]));

                    tl_gam_i[i] = 2.0 * (prod2_i(tauB_r[i],tauB_i[i],laplC_r[i],laplC_i[i])
                                      + prod2_i(tauC_r[i],tauC_i[i],laplB_r[i],laplB_i[i]));

                    gam_x_r[i] = 2.0 * prod2_r(gradB_x_r[i],gradB_x_i[i],rhoC_r[i],rhoC_i[i])
                               + 2.0 * prod2_r(gradC_x_r[i],gradC_x_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_y_r[i] = 2.0 * prod2_r(gradB_y_r[i],gradB_y_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_r(gradC_y_r[i],gradC_y_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_z_r[i] = 2.0 * prod2_r(gradB_z_r[i],gradB_z_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_r(gradC_z_r[i],gradC_z_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_x_i[i] = 2.0 * prod2_i(gradB_x_r[i],gradB_x_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_i(gradC_x_r[i],gradC_x_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_y_i[i] = 2.0 * prod2_i(gradB_y_r[i],gradB_y_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_i(gradC_y_r[i],gradC_y_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_z_i[i] = 2.0 * prod2_i(gradB_z_r[i],gradB_z_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_i(gradC_z_r[i],gradC_z_i[i],rhoB_r[i],rhoB_i[i]);
                                

                    st_gam_x_r[i] = 2.0 * prod2_r(gradB_x_r[i],gradB_x_i[i],tauC_r[i],tauC_i[i])
                                   + 2.0 * prod2_r(gradC_x_r[i],gradC_x_i[i],tauB_r[i],tauB_i[i]);                      

                    st_gam_y_r[i] = 2.0 * prod2_r(gradB_y_r[i],gradB_y_i[i],tauC_r[i],tauC_i[i])
                                   + 2.0 * prod2_r(gradC_y_r[i],gradC_y_i[i],tauB_r[i],tauB_i[i]);                      

                    st_gam_z_r[i] = 2.0 * prod2_r(gradB_z_r[i],gradB_z_i[i],tauC_r[i],tauC_i[i])
                                   + 2.0 * prod2_r(gradC_z_r[i],gradC_z_i[i],tauB_r[i],tauB_i[i]);                      

                    sl_gam_x_r[i] = 2.0 * prod2_r(gradB_x_r[i],gradB_x_i[i],laplC_r[i],laplC_i[i])
                                   + 2.0 * prod2_r(gradC_x_r[i],gradC_x_i[i],laplB_r[i],laplB_i[i]);                    

                    sl_gam_y_r[i] = 2.0 * prod2_r(gradB_y_r[i],gradB_y_i[i],laplC_r[i],laplC_i[i])
                                   + 2.0 * prod2_r(gradC_y_r[i],gradC_y_i[i],laplB_r[i],laplB_i[i]);                    

                    sl_gam_z_r[i] = 2.0 * prod2_r(gradB_z_r[i],gradB_z_i[i],laplC_r[i],laplC_i[i])
                                   + 2.0 * prod2_r(gradC_z_r[i],gradC_z_i[i],laplB_r[i],laplB_i[i]);                    

                    st_gam_x_i[i] = 2.0 * prod2_i(gradB_x_r[i],gradB_x_i[i],tauC_r[i],tauC_i[i])
                                   + 2.0 * prod2_i(gradC_x_r[i],gradC_x_i[i],tauB_r[i],tauB_i[i]);                      

                    st_gam_y_i[i] = 2.0 * prod2_i(gradB_y_r[i],gradB_y_i[i],tauC_r[i],tauC_i[i])
                                   + 2.0 * prod2_i(gradC_y_r[i],gradC_y_i[i],tauB_r[i],tauB_i[i]);                      

                    st_gam_z_i[i] = 2.0 * prod2_i(gradB_z_r[i],gradB_z_i[i],tauC_r[i],tauC_i[i])
                                   + 2.0 * prod2_i(gradC_z_r[i],gradC_z_i[i],tauB_r[i],tauB_i[i]);                      

                    sl_gam_x_i[i] = 2.0 * prod2_i(gradB_x_r[i],gradB_x_i[i],laplC_r[i],laplC_i[i])
                                   + 2.0 * prod2_i(gradC_x_r[i],gradC_x_i[i],laplB_r[i],laplB_i[i]);                    

                    sl_gam_y_i[i] = 2.0 * prod2_i(gradB_y_r[i],gradB_y_i[i],laplC_r[i],laplC_i[i])
                                   + 2.0 * prod2_i(gradC_y_r[i],gradC_y_i[i],laplB_r[i],laplB_i[i]);                    

                    sl_gam_z_i[i] = 2.0 * prod2_i(gradB_z_r[i],gradB_z_i[i],laplC_r[i],laplC_i[i])
                                   + 2.0 * prod2_i(gradC_z_r[i],gradC_z_i[i],laplB_r[i],laplB_i[i]);                    

                    gam_xx_r[i] = prod2_r(gradB_x_r[i],gradB_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                +  prod2_r(gradC_x_r[i],gradC_x_i[i],gradB_x_r[i],gradB_x_i[i]);                    

                    gam_xy_r[i] = prod2_r(gradB_x_r[i],gradB_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                +  prod2_r(gradC_x_r[i],gradC_x_i[i],gradB_y_r[i],gradB_y_i[i]);                    

                    gam_xz_r[i] = prod2_r(gradB_x_r[i],gradB_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                +  prod2_r(gradC_x_r[i],gradC_x_i[i],gradB_z_r[i],gradB_z_i[i]);                    

                    gam_yx_r[i] = prod2_r(gradB_y_r[i],gradB_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                +  prod2_r(gradC_y_r[i],gradC_y_i[i],gradB_x_r[i],gradB_x_i[i]);                    

                    gam_yy_r[i] = prod2_r(gradB_y_r[i],gradB_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                +  prod2_r(gradC_y_r[i],gradC_y_i[i],gradB_y_r[i],gradB_y_i[i]);                    

                    gam_yz_r[i] = prod2_r(gradB_y_r[i],gradB_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                +  prod2_r(gradC_y_r[i],gradC_y_i[i],gradB_z_r[i],gradB_z_i[i]);                    

                    gam_zx_r[i] = prod2_r(gradB_z_r[i],gradB_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                +  prod2_r(gradC_z_r[i],gradC_z_i[i],gradB_x_r[i],gradB_x_i[i]);                    

                    gam_zy_r[i] = prod2_r(gradB_z_r[i],gradB_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                +  prod2_r(gradC_z_r[i],gradC_z_i[i],gradB_y_r[i],gradB_y_i[i]);                    

                    gam_zz_r[i] = prod2_r(gradB_z_r[i],gradB_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                +  prod2_r(gradC_z_r[i],gradC_z_i[i],gradB_z_r[i],gradB_z_i[i]);                    

                    gam_xx_i[i] = prod2_i(gradB_x_r[i],gradB_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                +  prod2_i(gradC_x_r[i],gradC_x_i[i],gradB_x_r[i],gradB_x_i[i]);                    

                    gam_xy_i[i] = prod2_i(gradB_x_r[i],gradB_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                +  prod2_i(gradC_x_r[i],gradC_x_i[i],gradB_y_r[i],gradB_y_i[i]);                    

                    gam_xz_i[i] = prod2_i(gradB_x_r[i],gradB_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                +  prod2_i(gradC_x_r[i],gradC_x_i[i],gradB_z_r[i],gradB_z_i[i]);                    

                    gam_yx_i[i] = prod2_i(gradB_y_r[i],gradB_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                +  prod2_i(gradC_y_r[i],gradC_y_i[i],gradB_x_r[i],gradB_x_i[i]);                    

                    gam_yy_i[i] = prod2_i(gradB_y_r[i],gradB_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                +  prod2_i(gradC_y_r[i],gradC_y_i[i],gradB_y_r[i],gradB_y_i[i]);                    

                    gam_yz_i[i] = prod2_i(gradB_y_r[i],gradB_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                +  prod2_i(gradC_y_r[i],gradC_y_i[i],gradB_z_r[i],gradB_z_i[i]);                    

                    gam_zx_i[i] = prod2_i(gradB_z_r[i],gradB_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                +  prod2_i(gradC_z_r[i],gradC_z_i[i],gradB_x_r[i],gradB_x_i[i]);                    

                    gam_zy_i[i] = prod2_i(gradB_z_r[i],gradB_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                +  prod2_i(gradC_z_r[i],gradC_z_i[i],gradB_y_r[i],gradB_y_i[i]);                    

                    gam_zz_i[i] = prod2_i(gradB_z_r[i],gradB_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                +  prod2_i(gradC_z_r[i],gradC_z_i[i],gradB_z_r[i],gradB_z_i[i]);

                }
            }
        }
        if (fstr::upcase(quadMode) == "CRF_II")
        {
            for (int32_t j = 0; j < numdens / 2; j++)
            {
                // Perturbed densities

                // Perturbed density products to be stored

                auto gam_r = gam(2 * j);
                auto rt_gam_r = rt_gam(2 * j);
                auto rl_gam_r = rl_gam(2 * j);
                auto tt_gam_r = tt_gam(2 * j);
                auto tl_gam_r = tl_gam(2 * j);
                auto ll_gam_r = ll_gam(2 * j);
                
                auto gam_i = gam(2 * j + 1);
                auto rt_gam_i = rt_gam(2 * j + 1);
                auto rl_gam_i = rl_gam(2 * j + 1);
                auto tt_gam_i = tt_gam(2 * j + 1);
                auto tl_gam_i = tl_gam(2 * j + 1);
                auto ll_gam_i = ll_gam(2 * j + 1);

                auto gam_X_r = gamX(2 * j);
                auto gam_Y_r = gamY(2 * j);
                auto gam_Z_r = gamZ(2 * j);

                auto gam_X_i = gamX(2 * j + 1);
                auto gam_Y_i = gamY(2 * j + 1);
                auto gam_Z_i = gamZ(2 * j + 1);

                auto st_gam_X_r = st_gamX(2 * j);
                auto st_gam_Y_r = st_gamY(2 * j);
                auto st_gam_Z_r = st_gamZ(2 * j);

                auto sl_gam_X_r = sl_gamX(2 * j);
                auto sl_gam_Y_r = sl_gamY(2 * j);
                auto sl_gam_Z_r = sl_gamZ(2 * j);

                auto st_gam_X_i = st_gamX(2 * j + 1);
                auto st_gam_Y_i = st_gamY(2 * j + 1);
                auto st_gam_Z_i = st_gamZ(2 * j + 1);

                auto sl_gam_X_i = sl_gamX(2 * j + 1);
                auto sl_gam_Y_i = sl_gamY(2 * j + 1);
                auto sl_gam_Z_i = sl_gamZ(2 * j + 1);

                auto gam_XX_r = gamXX(2 * j);
                auto gam_XX_i = gamXX(2 * j + 1);
                auto gam_XY_r = gamXY(2 * j);
                auto gam_XY_i = gamXY(2 * j + 1);
                auto gam_XZ_r = gamXZ(2 * j);
                auto gam_XZ_i = gamXZ(2 * j + 1);
                auto gam_YX_r = gamYX(2 * j);
                auto gam_YX_i = gamYX(2 * j + 1);
                auto gam_YY_r = gamYY(2 * j);
                auto gam_YY_i = gamYY(2 * j + 1);
                auto gam_YZ_r = gamYZ(2 * j);
                auto gam_YZ_i = gamYZ(2 * j + 1);
                auto gam_ZX_r = gamZX(2 * j);
                auto gam_ZX_i = gamZX(2 * j + 1);
                auto gam_ZY_r = gamZY(2 * j);
                auto gam_ZY_i = gamZY(2 * j + 1);
                auto gam_ZZ_r = gamZZ(2 * j);
                auto gam_ZZ_i = gamZZ(2 * j + 1);

                // Input densities 

                auto rhoB_r = rwDensityGrid.alphaDensity(12 * j);
                auto tauB_r = rwDensityGrid.alphaDensitytau(12 * j);
                auto laplB_r = rwDensityGrid.alphaDensitylapl(12 * j);
                auto gradB_x_r = rwDensityGrid.alphaDensityGradientX(12 * j);
                auto gradB_y_r = rwDensityGrid.alphaDensityGradientY(12 * j);
                auto gradB_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j);

                auto rhoB_i = rwDensityGrid.alphaDensity(12 * j + 1);
                auto tauB_i = rwDensityGrid.alphaDensitytau(12 * j + 1);
                auto laplB_i = rwDensityGrid.alphaDensitylapl(12 * j + 1);
                auto gradB_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 1);
                auto gradB_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 1);
                auto gradB_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 1);

                auto rhoC_r = rwDensityGrid.alphaDensity(12 * j + 2);
                auto tauC_r = rwDensityGrid.alphaDensitytau(12 * j + 2);
                auto laplC_r = rwDensityGrid.alphaDensitylapl(12 * j + 2);
                auto gradC_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 2);
                auto gradC_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 2);
                auto gradC_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 2);

                auto rhoC_i = rwDensityGrid.alphaDensity(12 * j + 3);
                auto tauC_i = rwDensityGrid.alphaDensitytau(12 * j + 3);
                auto laplC_i = rwDensityGrid.alphaDensitylapl(12 * j + 3);
                auto gradC_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 3);
                auto gradC_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 3);
                auto gradC_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 3);

                auto rhoD_r = rwDensityGrid.alphaDensity(12 * j + 4);
                auto tauD_r = rwDensityGrid.alphaDensitytau(12 * j + 4);
                auto laplD_r = rwDensityGrid.alphaDensitylapl(12 * j + 4);
                auto gradD_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 4);
                auto gradD_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 4);
                auto gradD_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 4);

                auto rhoD_i = rwDensityGrid.alphaDensity(12 * j + 5);
                auto tauD_i = rwDensityGrid.alphaDensitytau(12 * j + 5);
                auto laplD_i = rwDensityGrid.alphaDensitylapl(12 * j + 5);
                auto gradD_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 5);
                auto gradD_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 5);
                auto gradD_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 5);

                auto rhoBC_r = rwDensityGrid.alphaDensity(12 * j + 6);
                auto tauBC_r = rwDensityGrid.alphaDensitytau(12 * j + 6);
                auto laplBC_r = rwDensityGrid.alphaDensitylapl(12 * j + 6);
                auto gradBC_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 6);
                auto gradBC_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 6);
                auto gradBC_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 6);

                auto rhoBC_i = rwDensityGrid.alphaDensity(12 * j + 7);
                auto tauBC_i = rwDensityGrid.alphaDensitytau(12 * j + 7);
                auto laplBC_i = rwDensityGrid.alphaDensitylapl(12 * j + 7);
                auto gradBC_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 7);
                auto gradBC_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 7);
                auto gradBC_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 7);

                auto rhoBD_r = rwDensityGrid.alphaDensity(12 * j + 8);
                auto tauBD_r = rwDensityGrid.alphaDensitytau(12 * j + 8);
                auto laplBD_r = rwDensityGrid.alphaDensitylapl(12 * j + 8);
                auto gradBD_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 8);
                auto gradBD_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 8);
                auto gradBD_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 8);

                auto rhoBD_i = rwDensityGrid.alphaDensity(12 * j + 9);
                auto tauBD_i = rwDensityGrid.alphaDensitytau(12 * j + 9);
                auto laplBD_i = rwDensityGrid.alphaDensitylapl(12 * j + 9);
                auto gradBD_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 9);
                auto gradBD_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 9);
                auto gradBD_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 9);
        
                auto rhoCD_r = rwDensityGrid.alphaDensity(12 * j + 10);
                auto tauCD_r = rwDensityGrid.alphaDensitytau(12 * j + 10);
                auto laplCD_r = rwDensityGrid.alphaDensitylapl(12 * j + 10);
                auto gradCD_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 10);
                auto gradCD_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 10);
                auto gradCD_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 10);

                auto rhoCD_i = rwDensityGrid.alphaDensity(12 * j + 11);
                auto tauCD_i = rwDensityGrid.alphaDensitytau(12 * j + 11);
                auto laplCD_i = rwDensityGrid.alphaDensitylapl(12 * j + 11);
                auto gradCD_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 11);
                auto gradCD_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 11);
                auto gradCD_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 11);


                for (int32_t i = 0; i < npoints; i++)
                {

                    gam_r[i] += 2.0 * prod2_r(rhoBC_r[i],rhoBC_i[i],rhoD_r[i],rhoD_i[i])
                              + 2.0 * prod2_r(rhoBD_r[i],rhoBD_i[i],rhoC_r[i],rhoC_i[i])
                              + 2.0 * prod2_r(rhoCD_r[i],rhoCD_i[i],rhoB_r[i],rhoB_i[i]);

                    gam_i[i] += 2.0 * prod2_i(rhoBC_r[i],rhoBC_i[i],rhoD_r[i],rhoD_i[i])
                              + 2.0 * prod2_i(rhoBD_r[i],rhoBD_i[i],rhoC_r[i],rhoC_i[i])
                              + 2.0 * prod2_i(rhoCD_r[i],rhoCD_i[i],rhoB_r[i],rhoB_i[i]);

                    // tt term

                    tt_gam_r[i] += 2.0 * prod2_r(tauBC_r[i],tauBC_i[i],tauD_r[i],tauD_i[i])
                                 + 2.0 * prod2_r(tauBD_r[i],tauBD_i[i],tauC_r[i],tauC_i[i])
                                 + 2.0 * prod2_r(tauCD_r[i],tauCD_i[i],tauB_r[i],tauB_i[i]);

                    tt_gam_i[i] += 2.0 * prod2_i(tauBC_r[i],tauBC_i[i],tauD_r[i],tauD_i[i])
                                 + 2.0 * prod2_i(tauBD_r[i],tauBD_i[i],tauC_r[i],tauC_i[i])
                                 + 2.0 * prod2_i(tauCD_r[i],tauCD_i[i],tauB_r[i],tauB_i[i]);

                    // ll term
                    ll_gam_r[i] += 2.0 * prod2_r(laplBC_r[i],laplBC_i[i],laplD_r[i],laplD_i[i])
                                 + 2.0 * prod2_r(laplBD_r[i],laplBD_i[i],laplC_r[i],laplC_i[i])
                                 + 2.0 * prod2_r(laplCD_r[i],laplCD_i[i],laplB_r[i],laplB_i[i]);

                    ll_gam_i[i] += 2.0 * prod2_i(laplBC_r[i],laplBC_i[i],laplD_r[i],laplD_i[i])
                                 + 2.0 * prod2_i(laplBD_r[i],laplBD_i[i],laplC_r[i],laplC_i[i])
                                 + 2.0 * prod2_i(laplCD_r[i],laplCD_i[i],laplB_r[i],laplB_i[i]);

                    // rt term
                    rt_gam_r[i] += 2.0 * prod2_r(tauBC_r[i],tauBC_i[i],rhoD_r[i],rhoD_i[i])
                                 + 2.0 * prod2_r(tauBD_r[i],tauBD_i[i],rhoC_r[i],rhoC_i[i])
                                 + 2.0 * prod2_r(tauCD_r[i],tauCD_i[i],rhoB_r[i],rhoB_i[i])
                                 + 2.0 * prod2_r(rhoBC_r[i],rhoBC_i[i],tauD_r[i],tauD_i[i])
                                 + 2.0 * prod2_r(rhoBD_r[i],rhoBD_i[i],tauC_r[i],tauC_i[i])
                                 + 2.0 * prod2_r(rhoCD_r[i],rhoCD_i[i],tauB_r[i],tauB_i[i]);

                    rt_gam_i[i] += 2.0 * prod2_i(tauBC_r[i],tauBC_i[i],rhoD_r[i],rhoD_i[i])
                                 + 2.0 * prod2_i(tauBD_r[i],tauBD_i[i],rhoC_r[i],rhoC_i[i])
                                 + 2.0 * prod2_i(tauCD_r[i],tauCD_i[i],rhoB_r[i],rhoB_i[i])
                                 + 2.0 * prod2_i(rhoBC_r[i],rhoBC_i[i],tauD_r[i],tauD_i[i])
                                 + 2.0 * prod2_i(rhoBD_r[i],rhoBD_i[i],tauC_r[i],tauC_i[i])
                                 + 2.0 * prod2_i(rhoCD_r[i],rhoCD_i[i],tauB_r[i],tauB_i[i]);

                    gam_X_r[i] += 2.0 * prod2_r(gradBC_x_r[i],gradBC_x_i[i],rhoD_r[i],rhoD_i[i])
                                + 2.0 * prod2_r(gradBD_x_r[i],gradBD_x_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_r(gradCD_x_r[i],gradCD_x_i[i],rhoB_r[i],rhoB_i[i])
                                + 2.0 * prod2_r(gradD_x_r[i],gradD_x_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 2.0 * prod2_r(gradC_x_r[i],gradC_x_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 2.0 * prod2_r(gradB_x_r[i],gradB_x_i[i],rhoCD_r[i],rhoCD_i[i]);                   

                    gam_Y_r[i] += 2.0 * prod2_r(gradBC_y_r[i],gradBC_y_i[i],rhoD_r[i],rhoD_i[i])
                                + 2.0 * prod2_r(gradBD_y_r[i],gradBD_y_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_r(gradCD_y_r[i],gradCD_y_i[i],rhoB_r[i],rhoB_i[i])
                                + 2.0 * prod2_r(gradD_y_r[i],gradD_y_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 2.0 * prod2_r(gradC_y_r[i],gradC_y_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 2.0 * prod2_r(gradB_y_r[i],gradB_y_i[i],rhoCD_r[i],rhoCD_i[i]);                   

                    gam_Z_r[i] += 2.0 * prod2_r(gradBC_z_r[i],gradBC_z_i[i],rhoD_r[i],rhoD_i[i])
                                + 2.0 * prod2_r(gradBD_z_r[i],gradBD_z_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_r(gradCD_z_r[i],gradCD_z_i[i],rhoB_r[i],rhoB_i[i])
                                + 2.0 * prod2_r(gradD_z_r[i],gradD_z_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 2.0 * prod2_r(gradC_z_r[i],gradC_z_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 2.0 * prod2_r(gradB_z_r[i],gradB_z_i[i],rhoCD_r[i],rhoCD_i[i]);                   

                    gam_X_i[i] += 2.0 * prod2_i(gradBC_x_r[i],gradBC_x_i[i],rhoD_r[i],rhoD_i[i])
                                + 2.0 * prod2_i(gradBD_x_r[i],gradBD_x_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_i(gradCD_x_r[i],gradCD_x_i[i],rhoB_r[i],rhoB_i[i])
                                + 2.0 * prod2_i(gradD_x_r[i],gradD_x_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 2.0 * prod2_i(gradC_x_r[i],gradC_x_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 2.0 * prod2_i(gradB_x_r[i],gradB_x_i[i],rhoCD_r[i],rhoCD_i[i]);                   

                    gam_Y_i[i] += 2.0 * prod2_i(gradBC_y_r[i],gradBC_y_i[i],rhoD_r[i],rhoD_i[i])
                                + 2.0 * prod2_i(gradBD_y_r[i],gradBD_y_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_i(gradCD_y_r[i],gradCD_y_i[i],rhoB_r[i],rhoB_i[i])
                                + 2.0 * prod2_i(gradD_y_r[i],gradD_y_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 2.0 * prod2_i(gradC_y_r[i],gradC_y_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 2.0 * prod2_i(gradB_y_r[i],gradB_y_i[i],rhoCD_r[i],rhoCD_i[i]);                   

                    gam_Z_i[i] += 2.0 * prod2_i(gradBC_z_r[i],gradBC_z_i[i],rhoD_r[i],rhoD_i[i])
                                + 2.0 * prod2_i(gradBD_z_r[i],gradBD_z_i[i],rhoC_r[i],rhoC_i[i])
                                + 2.0 * prod2_i(gradCD_z_r[i],gradCD_z_i[i],rhoB_r[i],rhoB_i[i])
                                + 2.0 * prod2_i(gradD_z_r[i],gradD_z_i[i],rhoBC_r[i],rhoBC_i[i])
                                + 2.0 * prod2_i(gradC_z_r[i],gradC_z_i[i],rhoBD_r[i],rhoBD_i[i])
                                + 2.0 * prod2_i(gradB_z_r[i],gradB_z_i[i],rhoCD_r[i],rhoCD_i[i]);   
                                

                    // ST 
                    
                    st_gam_X_r[i] += 2.0 * prod2_r(gradBC_x_r[i],gradBC_x_i[i],tauD_r[i],tauD_i[i])
                                   + 2.0 * prod2_r(gradBD_x_r[i],gradBD_x_i[i],tauC_r[i],tauC_i[i])
                                   + 2.0 * prod2_r(gradCD_x_r[i],gradCD_x_i[i],tauB_r[i],tauB_i[i])
                                   + 2.0 * prod2_r(gradD_x_r[i],gradD_x_i[i],tauBC_r[i],tauBC_i[i])
                                   + 2.0 * prod2_r(gradC_x_r[i],gradC_x_i[i],tauBD_r[i],tauBD_i[i])
                                   + 2.0 * prod2_r(gradB_x_r[i],gradB_x_i[i],tauCD_r[i],tauCD_i[i]);                   

                    st_gam_Y_r[i] += 2.0 * prod2_r(gradBC_y_r[i],gradBC_y_i[i],tauD_r[i],tauD_i[i])
                                   + 2.0 * prod2_r(gradBD_y_r[i],gradBD_y_i[i],tauC_r[i],tauC_i[i])
                                   + 2.0 * prod2_r(gradCD_y_r[i],gradCD_y_i[i],tauB_r[i],tauB_i[i])
                                   + 2.0 * prod2_r(gradD_y_r[i],gradD_y_i[i],tauBC_r[i],tauBC_i[i])
                                   + 2.0 * prod2_r(gradC_y_r[i],gradC_y_i[i],tauBD_r[i],tauBD_i[i])
                                   + 2.0 * prod2_r(gradB_y_r[i],gradB_y_i[i],tauCD_r[i],tauCD_i[i]);                   

                    st_gam_Z_r[i] += 2.0 * prod2_r(gradBC_z_r[i],gradBC_z_i[i],tauD_r[i],tauD_i[i])
                                   + 2.0 * prod2_r(gradBD_z_r[i],gradBD_z_i[i],tauC_r[i],tauC_i[i])
                                   + 2.0 * prod2_r(gradCD_z_r[i],gradCD_z_i[i],tauB_r[i],tauB_i[i])
                                   + 2.0 * prod2_r(gradD_z_r[i],gradD_z_i[i],tauBC_r[i],tauBC_i[i])
                                   + 2.0 * prod2_r(gradC_z_r[i],gradC_z_i[i],tauBD_r[i],tauBD_i[i])
                                   + 2.0 * prod2_r(gradB_z_r[i],gradB_z_i[i],tauCD_r[i],tauCD_i[i]);                   

                    st_gam_X_i[i] += 2.0 * prod2_i(gradBC_x_r[i],gradBC_x_i[i],tauD_r[i],tauD_i[i])
                                   + 2.0 * prod2_i(gradBD_x_r[i],gradBD_x_i[i],tauC_r[i],tauC_i[i])
                                   + 2.0 * prod2_i(gradCD_x_r[i],gradCD_x_i[i],tauB_r[i],tauB_i[i])
                                   + 2.0 * prod2_i(gradD_x_r[i],gradD_x_i[i],tauBC_r[i],tauBC_i[i])
                                   + 2.0 * prod2_i(gradC_x_r[i],gradC_x_i[i],tauBD_r[i],tauBD_i[i])
                                   + 2.0 * prod2_i(gradB_x_r[i],gradB_x_i[i],tauCD_r[i],tauCD_i[i]);                   

                    st_gam_Y_i[i] += 2.0 * prod2_i(gradBC_y_r[i],gradBC_y_i[i],tauD_r[i],tauD_i[i])
                                   + 2.0 * prod2_i(gradBD_y_r[i],gradBD_y_i[i],tauC_r[i],tauC_i[i])
                                   + 2.0 * prod2_i(gradCD_y_r[i],gradCD_y_i[i],tauB_r[i],tauB_i[i])
                                   + 2.0 * prod2_i(gradD_y_r[i],gradD_y_i[i],tauBC_r[i],tauBC_i[i])
                                   + 2.0 * prod2_i(gradC_y_r[i],gradC_y_i[i],tauBD_r[i],tauBD_i[i])
                                   + 2.0 * prod2_i(gradB_y_r[i],gradB_y_i[i],tauCD_r[i],tauCD_i[i]);                   

                    st_gam_Z_i[i] += 2.0 * prod2_i(gradBC_z_r[i],gradBC_z_i[i],tauD_r[i],tauD_i[i])
                                   + 2.0 * prod2_i(gradBD_z_r[i],gradBD_z_i[i],tauC_r[i],tauC_i[i])
                                   + 2.0 * prod2_i(gradCD_z_r[i],gradCD_z_i[i],tauB_r[i],tauB_i[i])
                                   + 2.0 * prod2_i(gradD_z_r[i],gradD_z_i[i],tauBC_r[i],tauBC_i[i])
                                   + 2.0 * prod2_i(gradC_z_r[i],gradC_z_i[i],tauBD_r[i],tauBD_i[i])
                                   + 2.0 * prod2_i(gradB_z_r[i],gradB_z_i[i],tauCD_r[i],tauCD_i[i]);   


                    // SL 
                    sl_gam_X_r[i] += 2.0 * prod2_r(gradBC_x_r[i],gradBC_x_i[i],laplD_r[i],laplD_i[i])
                                   + 2.0 * prod2_r(gradBD_x_r[i],gradBD_x_i[i],laplC_r[i],laplC_i[i])
                                   + 2.0 * prod2_r(gradCD_x_r[i],gradCD_x_i[i],laplB_r[i],laplB_i[i])
                                   + 2.0 * prod2_r(gradD_x_r[i],gradD_x_i[i],laplBC_r[i],laplBC_i[i])
                                   + 2.0 * prod2_r(gradC_x_r[i],gradC_x_i[i],laplBD_r[i],laplBD_i[i])
                                   + 2.0 * prod2_r(gradB_x_r[i],gradB_x_i[i],laplCD_r[i],laplCD_i[i]);                   

                    sl_gam_Y_r[i] += 2.0 * prod2_r(gradBC_y_r[i],gradBC_y_i[i],laplD_r[i],laplD_i[i])
                                   + 2.0 * prod2_r(gradBD_y_r[i],gradBD_y_i[i],laplC_r[i],laplC_i[i])
                                   + 2.0 * prod2_r(gradCD_y_r[i],gradCD_y_i[i],laplB_r[i],laplB_i[i])
                                   + 2.0 * prod2_r(gradD_y_r[i],gradD_y_i[i],laplBC_r[i],laplBC_i[i])
                                   + 2.0 * prod2_r(gradC_y_r[i],gradC_y_i[i],laplBD_r[i],laplBD_i[i])
                                   + 2.0 * prod2_r(gradB_y_r[i],gradB_y_i[i],laplCD_r[i],laplCD_i[i]);                   

                    sl_gam_Z_r[i] += 2.0 * prod2_r(gradBC_z_r[i],gradBC_z_i[i],laplD_r[i],laplD_i[i])
                                   + 2.0 * prod2_r(gradBD_z_r[i],gradBD_z_i[i],laplC_r[i],laplC_i[i])
                                   + 2.0 * prod2_r(gradCD_z_r[i],gradCD_z_i[i],laplB_r[i],laplB_i[i])
                                   + 2.0 * prod2_r(gradD_z_r[i],gradD_z_i[i],laplBC_r[i],laplBC_i[i])
                                   + 2.0 * prod2_r(gradC_z_r[i],gradC_z_i[i],laplBD_r[i],laplBD_i[i])
                                   + 2.0 * prod2_r(gradB_z_r[i],gradB_z_i[i],laplCD_r[i],laplCD_i[i]);                   

                    sl_gam_X_i[i] += 2.0 * prod2_i(gradBC_x_r[i],gradBC_x_i[i],laplD_r[i],laplD_i[i])
                                   + 2.0 * prod2_i(gradBD_x_r[i],gradBD_x_i[i],laplC_r[i],laplC_i[i])
                                   + 2.0 * prod2_i(gradCD_x_r[i],gradCD_x_i[i],laplB_r[i],laplB_i[i])
                                   + 2.0 * prod2_i(gradD_x_r[i],gradD_x_i[i],laplBC_r[i],laplBC_i[i])
                                   + 2.0 * prod2_i(gradC_x_r[i],gradC_x_i[i],laplBD_r[i],laplBD_i[i])
                                   + 2.0 * prod2_i(gradB_x_r[i],gradB_x_i[i],laplCD_r[i],laplCD_i[i]);                   

                    sl_gam_Y_i[i] += 2.0 * prod2_i(gradBC_y_r[i],gradBC_y_i[i],laplD_r[i],laplD_i[i])
                                   + 2.0 * prod2_i(gradBD_y_r[i],gradBD_y_i[i],laplC_r[i],laplC_i[i])
                                   + 2.0 * prod2_i(gradCD_y_r[i],gradCD_y_i[i],laplB_r[i],laplB_i[i])
                                   + 2.0 * prod2_i(gradD_y_r[i],gradD_y_i[i],laplBC_r[i],laplBC_i[i])
                                   + 2.0 * prod2_i(gradC_y_r[i],gradC_y_i[i],laplBD_r[i],laplBD_i[i])
                                   + 2.0 * prod2_i(gradB_y_r[i],gradB_y_i[i],laplCD_r[i],laplCD_i[i]);                   

                    sl_gam_Z_i[i] += 2.0 * prod2_i(gradBC_z_r[i],gradBC_z_i[i],laplD_r[i],laplD_i[i])
                                   + 2.0 * prod2_i(gradBD_z_r[i],gradBD_z_i[i],laplC_r[i],laplC_i[i])
                                   + 2.0 * prod2_i(gradCD_z_r[i],gradCD_z_i[i],laplB_r[i],laplB_i[i])
                                   + 2.0 * prod2_i(gradD_z_r[i],gradD_z_i[i],laplBC_r[i],laplBC_i[i])
                                   + 2.0 * prod2_i(gradC_z_r[i],gradC_z_i[i],laplBD_r[i],laplBD_i[i])
                                   + 2.0 * prod2_i(gradB_z_r[i],gradB_z_i[i],laplCD_r[i],laplCD_i[i]);   


                    gam_XX_r[i] += prod2_r(gradBC_x_r[i],gradBC_x_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_r(gradBD_x_r[i],gradBD_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gradCD_x_r[i],gradCD_x_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_r(gradBC_x_r[i],gradBC_x_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_r(gradBD_x_r[i],gradBD_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gradCD_x_r[i],gradCD_x_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_XY_r[i] += prod2_r(gradBC_x_r[i],gradBC_x_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_r(gradBD_x_r[i],gradBD_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gradCD_x_r[i],gradCD_x_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_r(gradBC_y_r[i],gradBC_y_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_r(gradBD_y_r[i],gradBD_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gradCD_y_r[i],gradCD_y_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_XZ_r[i] += prod2_r(gradBC_x_r[i],gradBC_x_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_r(gradBD_x_r[i],gradBD_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gradCD_x_r[i],gradCD_x_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_r(gradBC_z_r[i],gradBC_z_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_r(gradBD_z_r[i],gradBD_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gradCD_z_r[i],gradCD_z_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_YX_r[i] += prod2_r(gradBC_y_r[i],gradBC_y_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_r(gradBD_y_r[i],gradBD_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gradCD_y_r[i],gradCD_y_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_r(gradBC_x_r[i],gradBC_x_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_r(gradBD_x_r[i],gradBD_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gradCD_x_r[i],gradCD_x_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_YY_r[i] += prod2_r(gradBC_y_r[i],gradBC_y_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_r(gradBD_y_r[i],gradBD_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gradCD_y_r[i],gradCD_y_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_r(gradBC_y_r[i],gradBC_y_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_r(gradBD_y_r[i],gradBD_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gradCD_y_r[i],gradCD_y_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_YZ_r[i] += prod2_r(gradBC_y_r[i],gradBC_y_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_r(gradBD_y_r[i],gradBD_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gradCD_y_r[i],gradCD_y_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_r(gradBC_z_r[i],gradBC_z_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_r(gradBD_z_r[i],gradBD_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gradCD_z_r[i],gradCD_z_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_ZX_r[i] += prod2_r(gradBC_z_r[i],gradBC_z_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_r(gradBD_z_r[i],gradBD_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_r(gradCD_z_r[i],gradCD_z_i[i],gradB_x_r[i],gradB_x_i[i])
                                 + prod2_r(gradBC_x_r[i],gradBC_x_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_r(gradBD_x_r[i],gradBD_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gradCD_x_r[i],gradCD_x_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_ZY_r[i] += prod2_r(gradBC_z_r[i],gradBC_z_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_r(gradBD_z_r[i],gradBD_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_r(gradCD_z_r[i],gradCD_z_i[i],gradB_y_r[i],gradB_y_i[i])
                                 + prod2_r(gradBC_y_r[i],gradBC_y_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_r(gradBD_y_r[i],gradBD_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gradCD_y_r[i],gradCD_y_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_ZZ_r[i] += prod2_r(gradBC_z_r[i],gradBC_z_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_r(gradBD_z_r[i],gradBD_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gradCD_z_r[i],gradCD_z_i[i],gradB_z_r[i],gradB_z_i[i])
                                 + prod2_r(gradBC_z_r[i],gradBC_z_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_r(gradBD_z_r[i],gradBD_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_r(gradCD_z_r[i],gradCD_z_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_XX_i[i] += prod2_i(gradBC_x_r[i],gradBC_x_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_i(gradBD_x_r[i],gradBD_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gradCD_x_r[i],gradCD_x_i[i],gradB_x_r[i],gradB_x_i[i])
                                  + prod2_i(gradBC_x_r[i],gradBC_x_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_i(gradBD_x_r[i],gradBD_x_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gradCD_x_r[i],gradCD_x_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_XY_i[i] += prod2_i(gradBC_x_r[i],gradBC_x_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_i(gradBD_x_r[i],gradBD_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gradCD_x_r[i],gradCD_x_i[i],gradB_y_r[i],gradB_y_i[i])
                                  + prod2_i(gradBC_y_r[i],gradBC_y_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_i(gradBD_y_r[i],gradBD_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gradCD_y_r[i],gradCD_y_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_XZ_i[i] += prod2_i(gradBC_x_r[i],gradBC_x_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_i(gradBD_x_r[i],gradBD_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gradCD_x_r[i],gradCD_x_i[i],gradB_z_r[i],gradB_z_i[i])
                                  + prod2_i(gradBC_z_r[i],gradBC_z_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_i(gradBD_z_r[i],gradBD_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gradCD_z_r[i],gradCD_z_i[i],gradB_x_r[i],gradB_x_i[i]);

                    gam_YX_i[i] += prod2_i(gradBC_y_r[i],gradBC_y_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_i(gradBD_y_r[i],gradBD_y_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gradCD_y_r[i],gradCD_y_i[i],gradB_x_r[i],gradB_x_i[i])
                                  + prod2_i(gradBC_x_r[i],gradBC_x_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_i(gradBD_x_r[i],gradBD_x_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gradCD_x_r[i],gradCD_x_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_YY_i[i] += prod2_i(gradBC_y_r[i],gradBC_y_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_i(gradBD_y_r[i],gradBD_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gradCD_y_r[i],gradCD_y_i[i],gradB_y_r[i],gradB_y_i[i])
                                  + prod2_i(gradBC_y_r[i],gradBC_y_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_i(gradBD_y_r[i],gradBD_y_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gradCD_y_r[i],gradCD_y_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_YZ_i[i] += prod2_i(gradBC_y_r[i],gradBC_y_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_i(gradBD_y_r[i],gradBD_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gradCD_y_r[i],gradCD_y_i[i],gradB_z_r[i],gradB_z_i[i])
                                  + prod2_i(gradBC_z_r[i],gradBC_z_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_i(gradBD_z_r[i],gradBD_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gradCD_z_r[i],gradCD_z_i[i],gradB_y_r[i],gradB_y_i[i]);

                    gam_ZX_i[i] += prod2_i(gradBC_z_r[i],gradBC_z_i[i],gradD_x_r[i],gradD_x_i[i])
                                 + prod2_i(gradBD_z_r[i],gradBD_z_i[i],gradC_x_r[i],gradC_x_i[i])
                                 + prod2_i(gradCD_z_r[i],gradCD_z_i[i],gradB_x_r[i],gradB_x_i[i])
                                  + prod2_i(gradBC_x_r[i],gradBC_x_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_i(gradBD_x_r[i],gradBD_x_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gradCD_x_r[i],gradCD_x_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_ZY_i[i] += prod2_i(gradBC_z_r[i],gradBC_z_i[i],gradD_y_r[i],gradD_y_i[i])
                                 + prod2_i(gradBD_z_r[i],gradBD_z_i[i],gradC_y_r[i],gradC_y_i[i])
                                 + prod2_i(gradCD_z_r[i],gradCD_z_i[i],gradB_y_r[i],gradB_y_i[i])
                                  + prod2_i(gradBC_y_r[i],gradBC_y_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_i(gradBD_y_r[i],gradBD_y_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gradCD_y_r[i],gradCD_y_i[i],gradB_z_r[i],gradB_z_i[i]);

                    gam_ZZ_i[i] += prod2_i(gradBC_z_r[i],gradBC_z_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_i(gradBD_z_r[i],gradBD_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gradCD_z_r[i],gradCD_z_i[i],gradB_z_r[i],gradB_z_i[i])
                                  + prod2_i(gradBC_z_r[i],gradBC_z_i[i],gradD_z_r[i],gradD_z_i[i])
                                 + prod2_i(gradBD_z_r[i],gradBD_z_i[i],gradC_z_r[i],gradC_z_i[i])
                                 + prod2_i(gradCD_z_r[i],gradCD_z_i[i],gradB_z_r[i],gradB_z_i[i]);
               
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
