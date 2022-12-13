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

                auto rho_sig_x_r = rhow1rhow2(12 * j);

                auto rho_sig_x_i = rhow1rhow2(12 * j + 1);

                auto rho_sig_y_r = rhow1rhow2(12 * j + 2);

                auto rho_sig_y_i = rhow1rhow2(12 * j + 3);

                auto rho_sig_z_r = rhow1rhow2(12 * j + 4);

                auto rho_sig_z_i = rhow1rhow2(12 * j + 5);

                auto rho_lam_xy_r = rhow1rhow2(12 * j + 6);

                auto rho_lam_xy_i = rhow1rhow2(12 * j + 7);

                auto rho_lam_xz_r = rhow1rhow2(12 * j + 8);

                auto rho_lam_xz_i = rhow1rhow2(12 * j + 9);

                auto rho_lam_yz_r = rhow1rhow2(12 * j + 10);

                auto rho_lam_yz_i = rhow1rhow2(12 * j + 11);

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

                auto rho_sig_x_r = rhow1rhow2(6 * j);

                auto rho_sig_y_r = rhow1rhow2(6 * j + 1);

                auto rho_sig_z_r = rhow1rhow2(6 * j + 2);

                auto rho_lam_xy_r = rhow1rhow2(6 * j + 3);

                auto rho_lam_xz_r = rhow1rhow2(6 * j + 4);

                auto rho_lam_yz_r = rhow1rhow2(6 * j + 5);

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

        if (fstr::upcase(quadMode) == "TPA_I")
        {

            // This routine computes the first-order two-times transformed Fock matrices for the E[4] contraction for TPA calculations

            for (int32_t j = 0; j < numdens / 24; j++)
            {
                
                auto rr_sig_xx_r = rhow1rhow2(24 * j);

                auto rr_sig_xx_i = rhow1rhow2(24 * j + 1);

                auto rr_sig_yy_r = rhow1rhow2(24 * j + 2);

                auto rr_sig_yy_i = rhow1rhow2(24 * j + 3);

                auto rr_sig_zz_r = rhow1rhow2(24 * j + 4);

                auto rr_sig_zz_i = rhow1rhow2(24 * j + 5);

                auto rr_sig_xy_r = rhow1rhow2(24 * j + 6);

                auto rr_sig_xy_i = rhow1rhow2(24 * j + 7);

                auto rr_sig_xz_r = rhow1rhow2(24 * j + 8);

                auto rr_sig_xz_i = rhow1rhow2(24 * j + 9);

                auto rr_sig_yz_r = rhow1rhow2(24 * j + 10);

                auto rr_sig_yz_i = rhow1rhow2(24 * j + 11);

                auto rr_lamtau_xx_r = rhow1rhow2(24 * j + 12);

                auto rr_lamtau_xx_i = rhow1rhow2(24 * j + 13);

                auto rr_lamtau_yy_r = rhow1rhow2(24 * j + 14);

                auto rr_lamtau_yy_i = rhow1rhow2(24 * j + 15);

                auto rr_lamtau_zz_r = rhow1rhow2(24 * j + 16);

                auto rr_lamtau_zz_i = rhow1rhow2(24 * j + 17);

                auto rr_lamtau_xy_r = rhow1rhow2(24 * j + 18);

                auto rr_lamtau_xy_i = rhow1rhow2(24 * j + 19);

                auto rr_lamtau_xz_r = rhow1rhow2(24 * j + 20);

                auto rr_lamtau_xz_i = rhow1rhow2(24 * j + 21);

                auto rr_lamtau_yz_r = rhow1rhow2(24 * j + 22);

                auto rr_lamtau_yz_i = rhow1rhow2(24 * j + 23);


                auto rhoBx_r = rwDensityGrid.alphaDensity(12 * j);

                auto rhoBx_i = rwDensityGrid.alphaDensity(12 * j + 1);

                auto rhoBy_r = rwDensityGrid.alphaDensity(12 * j + 2);

                auto rhoBy_i = rwDensityGrid.alphaDensity(12 * j + 3);

                auto rhoBz_r = rwDensityGrid.alphaDensity(12 * j + 4);

                auto rhoBz_i = rwDensityGrid.alphaDensity(12 * j + 5);  

                auto rhoCx_r = rwDensityGrid.alphaDensity(12 * j + 6);

                auto rhoCx_i = rwDensityGrid.alphaDensity(12 * j + 7);

                auto rhoCy_r = rwDensityGrid.alphaDensity(12 * j + 8);

                auto rhoCy_i = rwDensityGrid.alphaDensity(12 * j + 9);

                auto rhoCz_r = rwDensityGrid.alphaDensity(12 * j + 10);

                auto rhoCz_i = rwDensityGrid.alphaDensity(12 * j + 11);  


                for (int32_t i = 0; i < npoints; i++)
                {

                    double rhobx_r = rhoBx_r[i]; double rhobx_i = rhoBx_i[i]; double rhoby_r = rhoBy_r[i];

                    double rhoby_i = rhoBy_i[i]; double rhobz_r = rhoBz_r[i]; double rhobz_i = rhoBz_i[i];

                    double rhodx_r = rhoBx_r[i]; double rhodx_i = rhoBx_i[i]; double rhody_r = rhoBy_r[i];

                    double rhody_i = rhoBy_i[i]; double rhodz_r = rhoBz_r[i]; double rhodz_i = rhoBz_i[i];

                    double rhocx_r = rhoCx_r[i]; double rhocx_i = rhoCx_i[i]; double rhocy_r = rhoCy_r[i];

                    double rhocy_i = rhoCy_i[i]; double rhocz_r = rhoCz_r[i]; double rhocz_i = rhoCz_i[i];


                    double sig_fac_r = 6.0 * ( ( rhobx_r* rhodx_r- rhobx_i * rhodx_i) +  ( rhoby_r* rhody_r- rhoby_i * rhody_i) +  ( rhobz_r* rhodz_r- rhobz_i * rhodz_i) );

                    double sig_fac_i =  6.0 * ( ( rhobx_i* rhodx_r + rhobx_r * rhodx_i) + ( rhoby_i* rhody_r + rhoby_r * rhody_i) + ( rhobz_i* rhodz_r + rhobz_r * rhodz_i) );

                    double lamtau_fac_r = 12.0 *( ( rhobx_r* rhocx_r- rhobx_i * rhocx_i) + ( rhoby_r* rhocy_r- rhoby_i * rhocy_i)  + ( rhobz_r* rhocz_r- rhobz_i * rhocz_i) );

                    double lamtau_fac_i = 12.0 *( ( rhobx_i* rhocx_r + rhobx_r * rhocx_i) + ( rhoby_i* rhocy_r + rhoby_r * rhocy_i) + ( rhobz_i* rhocz_r + rhobz_r * rhocz_i) );


                    rr_sig_xx_r[i] = 12.0 *  ( rhobx_r* rhodx_r- rhobx_i * rhodx_i) + sig_fac_r ;

                    rr_sig_xy_r[i] = 12.0 *  ( rhobx_r* rhody_r- rhobx_i * rhody_i) ;

                    rr_sig_xz_r[i] = 12.0 *  ( rhobx_r* rhodz_r- rhobx_i * rhodz_i) ;

                    rr_sig_yy_r[i] = 12.0 *  ( rhoby_r* rhody_r- rhoby_i * rhody_i) + sig_fac_r ;

                    rr_sig_yz_r[i] = 12.0 *  ( rhoby_r* rhodz_r- rhoby_i * rhodz_i) ;

                    rr_sig_zz_r[i] = 12.0 *  ( rhobz_r* rhodz_r- rhobz_i * rhodz_i) + sig_fac_r ;


                    rr_sig_xx_i[i] = 12.0 *  ( rhobx_i* rhodx_r + rhobx_r * rhodx_i) + sig_fac_i ;

                    rr_sig_xy_i[i] = 12.0 *  ( rhobx_i* rhody_r + rhobx_r * rhody_i) ;

                    rr_sig_xz_i[i] =  12.0 *  ( rhobx_i* rhodz_r + rhobx_r * rhodz_i) ;

                    rr_sig_yy_i[i] = 12.0 *  ( rhoby_i* rhody_r + rhoby_r * rhody_i) + sig_fac_i ;

                    rr_sig_yz_i[i] = 12.0 *  ( rhoby_i* rhodz_r + rhoby_r * rhodz_i) ;

                    rr_sig_zz_i[i] = 12.0 *  ( rhobz_i* rhodz_r + rhobz_r * rhodz_i) + sig_fac_i ;


                    rr_lamtau_xx_r[i] = 12.0 *  ( rhobx_r* rhocx_r- rhobx_i * rhocx_i)

                                        + 12.0 *  ( rhobx_r* rhocx_r- rhobx_i * rhocx_i) + lamtau_fac_r ;


                    rr_lamtau_xy_r[i] = 12.0 *  ( rhobx_r* rhocy_r- rhobx_i * rhocy_i) 

                                          + 12.0 *  ( rhoby_r* rhocx_r- rhoby_i * rhocx_i);

                    rr_lamtau_xz_r[i] = 12.0 *  ( rhobx_r* rhocz_r- rhobx_i * rhocz_i) 

                                        + 12.0 *  ( rhobz_r* rhocx_r- rhobz_i * rhocx_i);


                    rr_lamtau_yy_r[i] = 12.0 *  ( rhoby_r* rhocy_r- rhoby_i * rhocy_i)

                                        + 12.0 *  ( rhoby_r* rhocy_r- rhoby_i * rhocy_i) + lamtau_fac_r ;

                    rr_lamtau_yz_r[i] = 12.0 *  ( rhoby_r* rhocz_r- rhoby_i * rhocz_i) 

                                         + 12.0 *  ( rhobz_r* rhocy_r- rhobz_i * rhocy_i);


                    rr_lamtau_zz_r[i] = 12.0 *  ( rhobz_r* rhocz_r- rhobz_i * rhocz_i)
                    
                                        + 12.0 *  ( rhobz_r* rhocz_r- rhobz_i * rhocz_i) + lamtau_fac_r ;


                    rr_lamtau_xx_i[i] = 12.0 *  ( rhobx_i* rhocx_r + rhobx_r * rhocx_i)
                    
                                    + 12.0 *  ( rhobx_i* rhocx_r + rhobx_r * rhocx_i) + lamtau_fac_i ;


                    rr_lamtau_xy_i[i] =  12.0 *  ( rhobx_i* rhocy_r + rhobx_r * rhocy_i) 

                                        + 12.0 *  ( rhoby_i* rhocx_r + rhoby_r * rhocx_i);


                    rr_lamtau_xz_i[i] =  12.0 *  ( rhobx_i* rhocz_r + rhobx_r * rhocz_i) 

                                        + 12.0 *  ( rhobz_i* rhocx_r + rhobz_r * rhocx_i);


                    rr_lamtau_yy_i[i] = 12.0 *  ( rhoby_i* rhocy_r + rhoby_r * rhocy_i)
                    
                                            + 12.0 *  ( rhoby_i* rhocy_r + rhoby_r * rhocy_i) + lamtau_fac_i ;


                    rr_lamtau_yz_i[i] =  12.0 *  ( rhoby_i* rhocz_r + rhoby_r * rhocz_i) 

                                            + 12.0 *  ( rhobz_i* rhocy_r + rhobz_r * rhocy_i);


                    rr_lamtau_zz_i[i] = 12.0 *  ( rhobz_i* rhocz_r + rhobz_r * rhocz_i)
                    
                                        + 12.0 *  ( rhobz_i* rhocz_r + rhobz_r * rhocz_i) + lamtau_fac_i ;

                }
            }
        }
        if (fstr::upcase(quadMode) == "TPA_II") 
        {
            // This code is inteded to compute F_b_sigma fock matrices for the final E3 contraction for tpa calculations.
            for (int32_t j = 0; j < numdens / 6; j++)
            {
                
                auto rr_x_r = rhow1rhow2(6 * j);

                auto rr_x_i = rhow1rhow2(6 * j + 1);

                auto rr_y_r = rhow1rhow2(6 * j + 2);

                auto rr_y_i = rhow1rhow2(6 * j + 3);

                auto rr_z_r = rhow1rhow2(6 * j + 4);

                auto rr_z_i = rhow1rhow2(6 * j + 5);


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

            for (int32_t j = 0; j < numdens / 12; j++)
            {
                
                auto rr_sig_xx_r = rhow1rhow2(12 * j);

                auto rr_sig_xx_i = rhow1rhow2(12 * j + 1);

                auto rr_sig_yy_r = rhow1rhow2(12 * j + 2);

                auto rr_sig_yy_i = rhow1rhow2(12 * j + 3);

                auto rr_sig_zz_r = rhow1rhow2(12 * j + 4);

                auto rr_sig_zz_i = rhow1rhow2(12 * j + 5);

                auto rr_sig_xy_r = rhow1rhow2(12 * j + 6);

                auto rr_sig_xy_i = rhow1rhow2(12 * j + 7);

                auto rr_sig_xz_r = rhow1rhow2(12 * j + 8);

                auto rr_sig_xz_i = rhow1rhow2(12 * j + 9);

                auto rr_sig_yz_r = rhow1rhow2(12 * j + 10);

                auto rr_sig_yz_i = rhow1rhow2(12 * j + 11);


                auto rhoBx_r = rwDensityGrid.alphaDensity(6 * j);

                auto rhoBx_i = rwDensityGrid.alphaDensity(6 * j + 1);

                auto rhoBy_r = rwDensityGrid.alphaDensity(6 * j + 2);

                auto rhoBy_i = rwDensityGrid.alphaDensity(6 * j + 3);

                auto rhoBz_r = rwDensityGrid.alphaDensity(6 * j + 4);

                auto rhoBz_i = rwDensityGrid.alphaDensity(6 * j + 5);  



                for (int32_t i = 0; i < npoints; i++)
                {

                    double rhobx_r = rhoBx_r[i]; double rhobx_i = rhoBx_i[i]; double rhoby_r = rhoBy_r[i];

                    double rhoby_i = rhoBy_i[i]; double rhobz_r = rhoBz_r[i]; double rhobz_i = rhoBz_i[i];

                    double rhodx_r = rhoBx_r[i]; double rhodx_i = rhoBx_i[i]; double rhody_r = rhoBy_r[i];

                    double rhody_i = rhoBy_i[i]; double rhodz_r = rhoBz_r[i]; double rhodz_i = rhoBz_i[i];



                    double sig_fac_r = 6.0 * ( ( rhobx_r* rhodx_r- rhobx_i * rhodx_i) +  ( rhoby_r* rhody_r- rhoby_i * rhody_i) +  ( rhobz_r* rhodz_r- rhobz_i * rhodz_i) );

                    double sig_fac_i =  6.0 * ( ( rhobx_i* rhodx_r + rhobx_r * rhodx_i) + ( rhoby_i* rhody_r + rhoby_r * rhody_i) + ( rhobz_i* rhodz_r + rhobz_r * rhodz_i) );


                    rr_sig_xx_r[i] = 12.0 *  ( rhobx_r* rhodx_r- rhobx_i * rhodx_i) + sig_fac_r ;

                    rr_sig_xy_r[i] = 12.0 *  ( rhobx_r* rhody_r- rhobx_i * rhody_i) ;

                    rr_sig_xz_r[i] = 12.0 *  ( rhobx_r* rhodz_r- rhobx_i * rhodz_i) ;

                    rr_sig_yy_r[i] = 12.0 *  ( rhoby_r* rhody_r- rhoby_i * rhody_i) + sig_fac_r ;

                    rr_sig_yz_r[i] = 12.0 *  ( rhoby_r* rhodz_r- rhoby_i * rhodz_i) ;

                    rr_sig_zz_r[i] = 12.0 *  ( rhobz_r* rhodz_r- rhobz_i * rhodz_i) + sig_fac_r ;


                    rr_sig_xx_i[i] = 12.0 *  ( rhobx_i* rhodx_r + rhobx_r * rhodx_i) + sig_fac_i ;

                    rr_sig_xy_i[i] = 12.0 *  ( rhobx_i* rhody_r + rhobx_r * rhody_i) ;

                    rr_sig_xz_i[i] =  12.0 *  ( rhobx_i* rhodz_r + rhobx_r * rhodz_i) ;

                    rr_sig_yy_i[i] = 12.0 *  ( rhoby_i* rhody_r + rhoby_r * rhody_i) + sig_fac_i ;

                    rr_sig_yz_i[i] = 12.0 *  ( rhoby_i* rhodz_r + rhoby_r * rhodz_i) ;

                    rr_sig_zz_i[i] = 12.0 *  ( rhobz_i* rhodz_r + rhobz_r * rhodz_i) + sig_fac_i ;


                }
            }
        }
        if (fstr::upcase(quadMode) == "REDTPA_II")  
        {
            // This code is inteded to compute F_b_cd fock matrices for the final E3 contraction for tpa calculations.

            for (int32_t j = 0; j < numdens / 6; j++)
            {
                
                auto rr_x_r = rhow1rhow2(6 * j);

                auto rr_x_i = rhow1rhow2(6 * j + 1);

                auto rr_y_r = rhow1rhow2(6 * j + 2);

                auto rr_y_i = rhow1rhow2(6 * j + 3);

                auto rr_z_r = rhow1rhow2(6 * j + 4);

                auto rr_z_i = rhow1rhow2(6 * j + 5);
                
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

                    double rhocx_r = rhoCx_r[i]; double rhocx_i = rhoCx_i[i]; double rhocy_r = rhoCy_r[i];

                    double rhocy_i = rhoCy_i[i]; double rhocz_r = rhoCz_r[i]; double rhocz_i = rhoCz_i[i];

                    double rho_sig_xx_r = Rho_sig_xx_r[i];

                    double rho_sig_xx_i = Rho_sig_xx_i[i];

                    double rho_sig_yy_r = Rho_sig_yy_r[i];

                    double rho_sig_yy_i = Rho_sig_yy_i[i] ;

                    double rho_sig_zz_r = Rho_sig_zz_r[i];

                    double rho_sig_zz_i = Rho_sig_zz_i[i];

                    double rho_sig_xy_r = Rho_sig_xy_r[i];

                    double rho_sig_xy_i = Rho_sig_xy_i[i];

                    double rho_sig_yx_r = Rho_sig_xy_r[i];

                    double rho_sig_yx_i = Rho_sig_xy_i[i];

                    double rho_sig_xz_r = Rho_sig_xz_r[i];

                    double rho_sig_xz_i = Rho_sig_xz_i[i];

                    double rho_sig_zx_r = Rho_sig_xz_r[i];

                    double rho_sig_zx_i = Rho_sig_xz_i[i];

                    double rho_sig_yz_r = Rho_sig_yz_r[i];

                    double rho_sig_yz_i = Rho_sig_yz_i[i];

                    double rho_sig_zy_r = Rho_sig_yz_r[i];

                    double rho_sig_zy_i = Rho_sig_yz_i[i];

            
                    rr_x_r[i] = rhocx_r* rho_sig_xx_r - rhocx_i * rho_sig_xx_i + rhocy_r* rho_sig_xy_r - rhocy_i * rho_sig_xy_i + rhocz_r* rho_sig_xz_r - rhocz_i * rho_sig_xz_i;
                    rr_y_r[i] = rhocx_r* rho_sig_yx_r - rhocx_i * rho_sig_yx_i + rhocy_r* rho_sig_yy_r - rhocy_i * rho_sig_yy_i + rhocz_r* rho_sig_yz_r - rhocz_i * rho_sig_yz_i;
                    rr_z_r[i] = rhocx_r* rho_sig_zx_r - rhocx_i * rho_sig_zx_i + rhocy_r* rho_sig_zy_r - rhocy_i * rho_sig_zy_i + rhocz_r* rho_sig_zz_r - rhocz_i * rho_sig_zz_i;

                    rr_x_i[i] = rhocx_i* rho_sig_xx_r + rhocx_r * rho_sig_xx_i + rhocy_i* rho_sig_xy_r + rhocy_r * rho_sig_xy_i + rhocz_i* rho_sig_xz_r + rhocz_r * rho_sig_xz_i;
                    rr_y_i[i] = rhocx_i* rho_sig_yx_r + rhocx_r * rho_sig_yx_i + rhocy_i* rho_sig_yy_r + rhocy_r * rho_sig_yy_i + rhocz_i* rho_sig_yz_r + rhocz_r * rho_sig_yz_i;
                    rr_z_i[i] = rhocx_i* rho_sig_zx_r + rhocx_r * rho_sig_zx_i + rhocy_i* rho_sig_zy_r + rhocy_r * rho_sig_zy_i + rhocz_i* rho_sig_zz_r + rhocz_r * rho_sig_zz_i;

                }
            }
        }
        if (fstr::upcase(quadMode) == "CRF_I")
        {
            // This routine is for computing the Fbc, Fbd, Fcd first-order fock matrices for the general cubic response function
            
            std::cout << "crf_i" << std::endl;

            for (int32_t j = 0; j < numdens / 6; j++)
            {
                auto rho_bc_r = rhow1rhow2(6 * j);

                auto rho_bc_i = rhow1rhow2(6 * j + 1);

                auto rho_bd_r = rhow1rhow2(6 * j + 2);

                auto rho_bd_i = rhow1rhow2(6 * j + 3);

                auto rho_cd_r = rhow1rhow2(6 * j + 4);

                auto rho_cd_i = rhow1rhow2(6 * j + 5);


                auto rhowba_r = rwDensityGrid.alphaDensity(6 * j);

                auto rhowba_i = rwDensityGrid.alphaDensity(6 * j + 1);

                auto rhowca_r = rwDensityGrid.alphaDensity(6 * j + 2);

                auto rhowca_i = rwDensityGrid.alphaDensity(6 * j + 3);

                auto rhowda_r = rwDensityGrid.alphaDensity(6 * j + 4);

                auto rhowda_i = rwDensityGrid.alphaDensity(6 * j + 5);

                for (int32_t i = 0; i < npoints; i++)
                {
                    rho_bc_r[i] =  2.0 * (rhowba_r[i] * rhowca_r[i] - rhowba_i[i] * rhowca_i[i]);

                    rho_bc_i[i] =  2.0 * (rhowba_r[i] * rhowca_i[i] + rhowba_i[i] * rhowca_r[i]);

                    rho_bd_r[i] =  2.0 * (rhowba_r[i] * rhowda_r[i] - rhowba_i[i] * rhowda_i[i]);

                    rho_bd_i[i] =  2.0 * (rhowba_r[i] * rhowda_i[i] + rhowba_i[i] * rhowda_r[i]);

                    rho_cd_r[i] =  2.0 * (rhowca_r[i] * rhowda_r[i] - rhowca_i[i] * rhowda_i[i]);

                    rho_cd_i[i] =  2.0 * (rhowca_r[i] * rhowda_i[i] + rhowca_i[i] * rhowda_r[i]);

                }
            }
        }
        if (fstr::upcase(quadMode) == "CRF_II")
        {
            // This routine is for computing the second-order fock matrices for the E[3] contraction of the general cubic response function routine

            for (int32_t j = 0; j < numdens / 2; j++)
            {
                auto rr_123_r = rhow1rhow2(2 * j);

                auto rr_123_i = rhow1rhow2(2 * j + 1);

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

                                +  2.0 * (rhowca_r[i]* rhowbd_r[i] - rhowca_i[i] * rhowbd_i[i])
                                
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
                auto rhorho_r = rhow1rhow2(2 * j);

                auto rhorho_i = rhow1rhow2(2 * j + 1);

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
        if (fstr::upcase(quadMode) == "TPA_II") 
        {
            // This code is inteded to compute F_b_sigma fock matrices for the final E3 contraction for tpa calculations.

            for (int32_t j = 0; j < numdens / 6; j++)
            {
                
                auto F_x_rhow1rhow2_r = rhow1rhow2(6 * j);

                auto F_x_rhow1rhow2_i = rhow1rhow2(6 * j + 1);

                auto F_x_rxw1rhow2_r = rxw1rhow2(6 * j);

                auto F_x_rxw1rhow2_i = rxw1rhow2(6 * j + 1);

                auto F_x_ryw1rhow2_r = ryw1rhow2(6 * j);

                auto F_x_ryw1rhow2_i = ryw1rhow2(6 * j + 1);

                auto F_x_rzw1rhow2_r = rzw1rhow2(6 * j);

                auto F_x_rzw1rhow2_i = rzw1rhow2(6 * j + 1);

                auto F_x_rxw1rxw2_r = rxw1rxw2(6 * j);

                auto F_x_rxw1rxw2_i = rxw1rxw2(6 * j + 1);

                auto F_x_rxw1ryw2_r = rxw1ryw2(6 * j);

                auto F_x_rxw1ryw2_i = rxw1ryw2(6 * j + 1);

                auto F_x_rxw1rzw2_r = rxw1rzw2(6 * j);

                auto F_x_rxw1rzw2_i = rxw1rzw2(6 * j + 1);

                auto F_x_ryw1rxw2_r = ryw1rxw2(6 * j);

                auto F_x_ryw1rxw2_i = ryw1rxw2(6 * j + 1);

                auto F_x_ryw1ryw2_r = ryw1ryw2(6 * j);

                auto F_x_ryw1ryw2_i = ryw1ryw2(6 * j + 1);

                auto F_x_ryw1rzw2_r = ryw1rzw2(6 * j);

                auto F_x_ryw1rzw2_i = ryw1rzw2(6 * j + 1);

                auto F_x_rzw1rxw2_r = rzw1rxw2(6 * j);

                auto F_x_rzw1rxw2_i = rzw1rxw2(6 * j + 1);

                auto F_x_rzw1ryw2_r = rzw1ryw2(6 * j);

                auto F_x_rzw1ryw2_i = rzw1ryw2(6 * j + 1);

                auto F_x_rzw1rzw2_r = rzw1rzw2(6 * j);

                auto F_x_rzw1rzw2_i = rzw1rzw2(6 * j + 1);

                // Fy

                auto F_y_rhow1rhow2_r = rhow1rhow2(6 * j  + 2);
 
                auto F_y_rhow1rhow2_i = rhow1rhow2(6 * j  + 3);

                auto F_y_rxw1rhow2_r = rxw1rhow2(6 * j  + 2);

                auto F_y_rxw1rhow2_i = rxw1rhow2(6 * j  + 3);

                auto F_y_ryw1rhow2_r = ryw1rhow2(6 * j  + 2);

                auto F_y_ryw1rhow2_i = ryw1rhow2(6 * j  + 3);

                auto F_y_rzw1rhow2_r = rzw1rhow2(6 * j  + 2);

                auto F_y_rzw1rhow2_i = rzw1rhow2(6 * j  + 3);

                auto F_y_rxw1rxw2_r = rxw1rxw2(6 * j  + 2);

                auto F_y_rxw1rxw2_i = rxw1rxw2(6 * j  + 3);

                auto F_y_rxw1ryw2_r = rxw1ryw2(6 * j  + 2);

                auto F_y_rxw1ryw2_i = rxw1ryw2(6 * j  + 3);

                auto F_y_rxw1rzw2_r = rxw1rzw2(6 * j  + 2);

                auto F_y_rxw1rzw2_i = rxw1rzw2(6 * j  + 3);

                auto F_y_ryw1rxw2_r = ryw1rxw2(6 * j  + 2);

                auto F_y_ryw1rxw2_i = ryw1rxw2(6 * j  + 3);

                auto F_y_ryw1ryw2_r = ryw1ryw2(6 * j  + 2);

                auto F_y_ryw1ryw2_i = ryw1ryw2(6 * j  + 3);

                auto F_y_ryw1rzw2_r = ryw1rzw2(6 * j  + 2);

                auto F_y_ryw1rzw2_i = ryw1rzw2(6 * j  + 3);

                auto F_y_rzw1rxw2_r = rzw1rxw2(6 * j  + 2);

                auto F_y_rzw1rxw2_i = rzw1rxw2(6 * j  + 3);

                auto F_y_rzw1ryw2_r = rzw1ryw2(6 * j  + 2);

                auto F_y_rzw1ryw2_i = rzw1ryw2(6 * j  + 3);

                auto F_y_rzw1rzw2_r = rzw1rzw2(6 * j  + 2);

                auto F_y_rzw1rzw2_i = rzw1rzw2(6 * j  + 3);

                // Fz

               auto F_z_rhow1rhow2_r = rhow1rhow2(6 * j  + 4);

                auto F_z_rhow1rhow2_i = rhow1rhow2(6 * j  + 5);

                auto F_z_rxw1rhow2_r = rxw1rhow2(6 * j  + 4);

                auto F_z_rxw1rhow2_i = rxw1rhow2(6 * j  + 5);

                auto F_z_ryw1rhow2_r = ryw1rhow2(6 * j  + 4);

                auto F_z_ryw1rhow2_i = ryw1rhow2(6 * j  + 5);

                auto F_z_rzw1rhow2_r = rzw1rhow2(6 * j  + 4);

                auto F_z_rzw1rhow2_i = rzw1rhow2(6 * j  + 5);

                auto F_z_rxw1rxw2_r = rxw1rxw2(6 * j  + 4);

                auto F_z_rxw1rxw2_i = rxw1rxw2(6 * j  + 5);

                auto F_z_rxw1ryw2_r = rxw1ryw2(6 * j  + 4);

                auto F_z_rxw1ryw2_i = rxw1ryw2(6 * j  + 5);

                auto F_z_rxw1rzw2_r = rxw1rzw2(6 * j  + 4);

                auto F_z_rxw1rzw2_i = rxw1rzw2(6 * j  + 5);

                auto F_z_ryw1rxw2_r = ryw1rxw2(6 * j  + 4);

                auto F_z_ryw1rxw2_i = ryw1rxw2(6 * j  + 5);

                auto F_z_ryw1ryw2_r = ryw1ryw2(6 * j  + 4);

                auto F_z_ryw1ryw2_i = ryw1ryw2(6 * j  + 5);

                auto F_z_ryw1rzw2_r = ryw1rzw2(6 * j  + 4);

                auto F_z_ryw1rzw2_i = ryw1rzw2(6 * j  + 5);

                auto F_z_rzw1rxw2_r = rzw1rxw2(6 * j  + 4);

                auto F_z_rzw1rxw2_i = rzw1rxw2(6 * j  + 5);

                auto F_z_rzw1ryw2_r = rzw1ryw2(6 * j  + 4);

                auto F_z_rzw1ryw2_i = rzw1ryw2(6 * j  + 5);

                auto F_z_rzw1rzw2_r = rzw1rzw2(6 * j  + 4);

                auto F_z_rzw1rzw2_i = rzw1rzw2(6 * j  + 5);


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
                    

                    F_x_rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_x_rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_x_rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_x_rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_x_ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_x_ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_x_rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_x_rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_x_rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_x_rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_x_rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_x_rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_x_rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_x_rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_x_ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_x_ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_x_ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_x_ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_x_ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_x_ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_x_rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_x_rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_x_rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_x_rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_x_rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_x_rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                        

                    F_x_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_x_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_x_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_x_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_x_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_x_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_x_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_x_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_x_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_x_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_x_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_x_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_x_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_x_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_x_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_x_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_x_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_x_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_x_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_x_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_x_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_x_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_x_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_x_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_x_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_x_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                        

                    F_x_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_x_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_x_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_x_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_x_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_x_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_x_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_x_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_x_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_x_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_x_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_x_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_x_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_x_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_x_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_x_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_x_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_x_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_x_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_x_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_x_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_x_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_x_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_x_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_x_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_x_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                    

                    F_x_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_x_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_x_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_x_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_x_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_x_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_x_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_x_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_x_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_x_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_x_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_x_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_x_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_x_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_x_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_x_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_x_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_x_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_x_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_x_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_x_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_x_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_x_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_x_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_x_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_x_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                        

                    F_x_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_x_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_x_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_x_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_x_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_x_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_x_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_x_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_x_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_x_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_x_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_x_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_x_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_x_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_x_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_x_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_x_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_x_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_x_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_x_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_x_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_x_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_x_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_x_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_x_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_x_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                        

                    F_x_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_x_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_x_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_x_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_x_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_x_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_x_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_x_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_x_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_x_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_x_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_x_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_x_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_x_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_x_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_x_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_x_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_x_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_x_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_x_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_x_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_x_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_x_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_x_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_x_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_x_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                    

                    F_y_rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_y_rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_y_rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_y_rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_y_ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_y_ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_y_rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_y_rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_y_rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_y_rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_y_rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_y_rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_y_rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_y_rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_y_ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_y_ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_y_ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_y_ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_y_ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_y_ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_y_rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_y_rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_y_rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_y_rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_y_rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_y_rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                        

                    F_y_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_y_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_y_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_y_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_y_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_y_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_y_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_y_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_y_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_y_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_y_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_y_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_y_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_y_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_y_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_y_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_y_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_y_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_y_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_y_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_y_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_y_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_y_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_y_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_y_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_y_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                        

                    F_y_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_y_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_y_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_y_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_y_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_y_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_y_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_y_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_y_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_y_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_y_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_y_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_y_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_y_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_y_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_y_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_y_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_y_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_y_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_y_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_y_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_y_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_y_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_y_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_y_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_y_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                    

                    F_y_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_y_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_y_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_y_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_y_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_y_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_y_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_y_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_y_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_y_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_y_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_y_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_y_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_y_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_y_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_y_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_y_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_y_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_y_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_y_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_y_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_y_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_y_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_y_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_y_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_y_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                        

                    F_y_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_y_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_y_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_y_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_y_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_y_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_y_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_y_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_y_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_y_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_y_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_y_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_y_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_y_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_y_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_y_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_y_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_y_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_y_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_y_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_y_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_y_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_y_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_y_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_y_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_y_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                        

                    F_y_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_y_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_y_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_y_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_y_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_y_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_y_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_y_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_y_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_y_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_y_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_y_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_y_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_y_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_y_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_y_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_y_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_y_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_y_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_y_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_y_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_y_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_y_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_y_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_y_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_y_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                

                    F_z_rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_z_rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_z_rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_z_rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_z_ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_z_ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_z_rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_z_rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_z_rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_z_rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_z_rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_z_rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_z_rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_z_rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_z_ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_z_ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_z_ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_z_ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_z_ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_z_ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_z_rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_z_rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_z_rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_z_rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_z_rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_z_rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                        

                    F_z_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_z_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_z_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_z_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_z_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_z_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_z_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_z_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_z_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_z_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_z_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_z_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_z_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_z_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_z_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_z_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_z_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_z_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_z_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_z_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_z_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_z_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_z_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_z_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_z_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_z_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                        

                    F_z_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_z_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_z_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_z_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_z_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_z_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_z_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_z_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_z_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_z_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_z_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_z_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_z_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_z_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_z_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_z_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_z_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_z_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_z_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_z_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_z_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_z_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_z_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_z_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_z_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_z_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                    

                    F_z_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_z_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_z_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_z_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_z_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_z_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_z_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_z_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_z_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_z_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_z_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_z_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_z_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_z_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_z_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_z_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_z_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_z_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_z_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_z_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_z_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_z_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_z_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_z_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_z_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_z_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                        

                    F_z_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_z_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_z_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_z_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_z_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_z_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_z_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_z_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_z_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_z_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_z_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_z_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_z_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_z_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_z_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_z_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_z_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_z_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_z_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_z_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_z_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_z_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_z_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_z_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_z_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_z_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

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
                        

                    F_z_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    F_z_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    F_z_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    F_z_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    F_z_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    F_z_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    F_z_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    F_z_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    F_z_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    F_z_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    F_z_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    F_z_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    F_z_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    F_z_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    F_z_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    F_z_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    F_z_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    F_z_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    F_z_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    F_z_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    F_z_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    F_z_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    F_z_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    F_z_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    F_z_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    F_z_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                }
            }
    }    
    if (fstr::upcase(quadMode) == "TPA_I")
        {
        for (int32_t j = 0; j < numdens / 24; j++)
            {
                
                // Sig_xx
                
                auto sig_xx_rhow1rhow2_r = rhow1rhow2(24 * j);

                auto sig_xx_rhow1rhow2_i = rhow1rhow2(24 * j + 1);

                auto sig_xx_rxw1rhow2_r = rxw1rhow2(24 * j);

                auto sig_xx_rxw1rhow2_i = rxw1rhow2(24 * j + 1);

                auto sig_xx_ryw1rhow2_r = ryw1rhow2(24 * j);

                auto sig_xx_ryw1rhow2_i = ryw1rhow2(24 * j + 1);

                auto sig_xx_rzw1rhow2_r = rzw1rhow2(24 * j);

                auto sig_xx_rzw1rhow2_i = rzw1rhow2(24 * j + 1);

                auto sig_xx_rxw1rxw2_r = rxw1rxw2(24 * j);

                auto sig_xx_rxw1rxw2_i = rxw1rxw2(24 * j + 1);

                auto sig_xx_rxw1ryw2_r = rxw1ryw2(24 * j);

                auto sig_xx_rxw1ryw2_i = rxw1ryw2(24 * j + 1);

                auto sig_xx_rxw1rzw2_r = rxw1rzw2(24 * j);

                auto sig_xx_rxw1rzw2_i = rxw1rzw2(24 * j + 1);

                auto sig_xx_ryw1rxw2_r = ryw1rxw2(24 * j);

                auto sig_xx_ryw1rxw2_i = ryw1rxw2(24 * j + 1);

                auto sig_xx_ryw1ryw2_r = ryw1ryw2(24 * j);

                auto sig_xx_ryw1ryw2_i = ryw1ryw2(24 * j + 1);

                auto sig_xx_ryw1rzw2_r = ryw1rzw2(24 * j);

                auto sig_xx_ryw1rzw2_i = ryw1rzw2(24 * j + 1);

                auto sig_xx_rzw1rxw2_r = rzw1rxw2(24 * j);

                auto sig_xx_rzw1rxw2_i = rzw1rxw2(24 * j + 1);

                auto sig_xx_rzw1ryw2_r = rzw1ryw2(24 * j);

                auto sig_xx_rzw1ryw2_i = rzw1ryw2(24 * j + 1);

                auto sig_xx_rzw1rzw2_r = rzw1rzw2(24 * j);

                auto sig_xx_rzw1rzw2_i = rzw1rzw2(24 * j + 1);


                // Sig_yy

                auto sig_yy_rhow1rhow2_r = rhow1rhow2(24 * j  + 2);
 
                auto sig_yy_rhow1rhow2_i = rhow1rhow2(24 * j  + 3);

                auto sig_yy_rxw1rhow2_r = rxw1rhow2(24 * j  + 2);

                auto sig_yy_rxw1rhow2_i = rxw1rhow2(24 * j  + 3);

                auto sig_yy_ryw1rhow2_r = ryw1rhow2(24 * j  + 2);

                auto sig_yy_ryw1rhow2_i = ryw1rhow2(24 * j  + 3);

                auto sig_yy_rzw1rhow2_r = rzw1rhow2(24 * j  + 2);

                auto sig_yy_rzw1rhow2_i = rzw1rhow2(24 * j  + 3);

                auto sig_yy_rxw1rxw2_r = rxw1rxw2(24 * j  + 2);

                auto sig_yy_rxw1rxw2_i = rxw1rxw2(24 * j  + 3);

                auto sig_yy_rxw1ryw2_r = rxw1ryw2(24 * j  + 2);

                auto sig_yy_rxw1ryw2_i = rxw1ryw2(24 * j  + 3);

                auto sig_yy_rxw1rzw2_r = rxw1rzw2(24 * j  + 2);

                auto sig_yy_rxw1rzw2_i = rxw1rzw2(24 * j  + 3);

                auto sig_yy_ryw1rxw2_r = ryw1rxw2(24 * j  + 2);

                auto sig_yy_ryw1rxw2_i = ryw1rxw2(24 * j  + 3);

                auto sig_yy_ryw1ryw2_r = ryw1ryw2(24 * j  + 2);

                auto sig_yy_ryw1ryw2_i = ryw1ryw2(24 * j  + 3);

                auto sig_yy_ryw1rzw2_r = ryw1rzw2(24 * j  + 2);

                auto sig_yy_ryw1rzw2_i = ryw1rzw2(24 * j  + 3);

                auto sig_yy_rzw1rxw2_r = rzw1rxw2(24 * j  + 2);

                auto sig_yy_rzw1rxw2_i = rzw1rxw2(24 * j  + 3);

                auto sig_yy_rzw1ryw2_r = rzw1ryw2(24 * j  + 2);

                auto sig_yy_rzw1ryw2_i = rzw1ryw2(24 * j  + 3);

                auto sig_yy_rzw1rzw2_r = rzw1rzw2(24 * j  + 2);

                auto sig_yy_rzw1rzw2_i = rzw1rzw2(24 * j  + 3);

                // Sig_zz

                auto sig_zz_rhow1rhow2_r = rhow1rhow2(24 * j  + 4);

                auto sig_zz_rhow1rhow2_i = rhow1rhow2(24 * j  + 5);

                auto sig_zz_rxw1rhow2_r = rxw1rhow2(24 * j  + 4);

                auto sig_zz_rxw1rhow2_i = rxw1rhow2(24 * j  + 5);

                auto sig_zz_ryw1rhow2_r = ryw1rhow2(24 * j  + 4);

                auto sig_zz_ryw1rhow2_i = ryw1rhow2(24 * j  + 5);

                auto sig_zz_rzw1rhow2_r = rzw1rhow2(24 * j  + 4);

                auto sig_zz_rzw1rhow2_i = rzw1rhow2(24 * j  + 5);

                auto sig_zz_rxw1rxw2_r = rxw1rxw2(24 * j  + 4);

                auto sig_zz_rxw1rxw2_i = rxw1rxw2(24 * j  + 5);

                auto sig_zz_rxw1ryw2_r = rxw1ryw2(24 * j  + 4);

                auto sig_zz_rxw1ryw2_i = rxw1ryw2(24 * j  + 5);

                auto sig_zz_rxw1rzw2_r = rxw1rzw2(24 * j  + 4);

                auto sig_zz_rxw1rzw2_i = rxw1rzw2(24 * j  + 5);

                auto sig_zz_ryw1rxw2_r = ryw1rxw2(24 * j  + 4);

                auto sig_zz_ryw1rxw2_i = ryw1rxw2(24 * j  + 5);

                auto sig_zz_ryw1ryw2_r = ryw1ryw2(24 * j  + 4);

                auto sig_zz_ryw1ryw2_i = ryw1ryw2(24 * j  + 5);

                auto sig_zz_ryw1rzw2_r = ryw1rzw2(24 * j  + 4);

                auto sig_zz_ryw1rzw2_i = ryw1rzw2(24 * j  + 5);

                auto sig_zz_rzw1rxw2_r = rzw1rxw2(24 * j  + 4);

                auto sig_zz_rzw1rxw2_i = rzw1rxw2(24 * j  + 5);

                auto sig_zz_rzw1ryw2_r = rzw1ryw2(24 * j  + 4);

                auto sig_zz_rzw1ryw2_i = rzw1ryw2(24 * j  + 5);

                auto sig_zz_rzw1rzw2_r = rzw1rzw2(24 * j  + 4);

                auto sig_zz_rzw1rzw2_i = rzw1rzw2(24 * j  + 5);

                // Sig_xy

                auto sig_xy_rhow1rhow2_r = rhow1rhow2(24 * j  + 6);

                auto sig_xy_rhow1rhow2_i = rhow1rhow2(24 * j  + 7);

                auto sig_xy_rxw1rhow2_r = rxw1rhow2(24 * j  + 6);

                auto sig_xy_rxw1rhow2_i = rxw1rhow2(24 * j  + 7);

                auto sig_xy_ryw1rhow2_r = ryw1rhow2(24 * j  + 6);

                auto sig_xy_ryw1rhow2_i = ryw1rhow2(24 * j  + 7);

                auto sig_xy_rzw1rhow2_r = rzw1rhow2(24 * j  + 6);

                auto sig_xy_rzw1rhow2_i = rzw1rhow2(24 * j  + 7);

                auto sig_xy_rxw1rxw2_r = rxw1rxw2(24 * j  + 6);

                auto sig_xy_rxw1rxw2_i = rxw1rxw2(24 * j  + 7);

                auto sig_xy_rxw1ryw2_r = rxw1ryw2(24 * j  + 6);

                auto sig_xy_rxw1ryw2_i = rxw1ryw2(24 * j  + 7);

                auto sig_xy_rxw1rzw2_r = rxw1rzw2(24 * j  + 6);

                auto sig_xy_rxw1rzw2_i = rxw1rzw2(24 * j  + 7);

                auto sig_xy_ryw1rxw2_r = ryw1rxw2(24 * j  + 6);

                auto sig_xy_ryw1rxw2_i = ryw1rxw2(24 * j  + 7);

                auto sig_xy_ryw1ryw2_r = ryw1ryw2(24 * j  + 6);

                auto sig_xy_ryw1ryw2_i = ryw1ryw2(24 * j  + 7);

                auto sig_xy_ryw1rzw2_r = ryw1rzw2(24 * j  + 6);

                auto sig_xy_ryw1rzw2_i = ryw1rzw2(24 * j  + 7);

                auto sig_xy_rzw1rxw2_r = rzw1rxw2(24 * j  + 6);

                auto sig_xy_rzw1rxw2_i = rzw1rxw2(24 * j  + 7);

                auto sig_xy_rzw1ryw2_r = rzw1ryw2(24 * j  + 6);

                auto sig_xy_rzw1ryw2_i = rzw1ryw2(24 * j  + 7);

                auto sig_xy_rzw1rzw2_r = rzw1rzw2(24 * j  + 6);

                auto sig_xy_rzw1rzw2_i = rzw1rzw2(24 * j  + 7);

                // Sig_xz

                auto sig_xz_rhow1rhow2_r = rhow1rhow2(24 * j  + 8);
 
                auto sig_xz_rhow1rhow2_i = rhow1rhow2(24 * j  + 9);

                auto sig_xz_rxw1rhow2_r = rxw1rhow2(24 * j  + 8);

                auto sig_xz_rxw1rhow2_i = rxw1rhow2(24 * j  + 9);

                auto sig_xz_ryw1rhow2_r = ryw1rhow2(24 * j  + 8);

                auto sig_xz_ryw1rhow2_i = ryw1rhow2(24 * j  + 9);

                auto sig_xz_rzw1rhow2_r = rzw1rhow2(24 * j  + 8);

                auto sig_xz_rzw1rhow2_i = rzw1rhow2(24 * j  + 9);

                auto sig_xz_rxw1rxw2_r = rxw1rxw2(24 * j  + 8);

                auto sig_xz_rxw1rxw2_i = rxw1rxw2(24 * j  + 9);

                auto sig_xz_rxw1ryw2_r = rxw1ryw2(24 * j  + 8);

                auto sig_xz_rxw1ryw2_i = rxw1ryw2(24 * j  + 9);

                auto sig_xz_rxw1rzw2_r = rxw1rzw2(24 * j  + 8);

                auto sig_xz_rxw1rzw2_i = rxw1rzw2(24 * j  + 9);

                auto sig_xz_ryw1rxw2_r = ryw1rxw2(24 * j  + 8);

                auto sig_xz_ryw1rxw2_i = ryw1rxw2(24 * j  + 9);

                auto sig_xz_ryw1ryw2_r = ryw1ryw2(24 * j  + 8);

                auto sig_xz_ryw1ryw2_i = ryw1ryw2(24 * j  + 9);

                auto sig_xz_ryw1rzw2_r = ryw1rzw2(24 * j  + 8);

                auto sig_xz_ryw1rzw2_i = ryw1rzw2(24 * j  + 9);

                auto sig_xz_rzw1rxw2_r = rzw1rxw2(24 * j  + 8);

                auto sig_xz_rzw1rxw2_i = rzw1rxw2(24 * j  + 9);

                auto sig_xz_rzw1ryw2_r = rzw1ryw2(24 * j  + 8);

                auto sig_xz_rzw1ryw2_i = rzw1ryw2(24 * j  + 9);

                auto sig_xz_rzw1rzw2_r = rzw1rzw2(24 * j  + 8);

                auto sig_xz_rzw1rzw2_i = rzw1rzw2(24 * j  + 9);

                // Sig_yz

                auto sig_yz_rhow1rhow2_r = rhow1rhow2(24 * j  + 10);

                auto sig_yz_rhow1rhow2_i = rhow1rhow2(24 * j  + 11);

                auto sig_yz_rxw1rhow2_r = rxw1rhow2(24 * j  + 10);

                auto sig_yz_rxw1rhow2_i = rxw1rhow2(24 * j  + 11);

                auto sig_yz_ryw1rhow2_r = ryw1rhow2(24 * j  + 10);

                auto sig_yz_ryw1rhow2_i = ryw1rhow2(24 * j  + 11);

                auto sig_yz_rzw1rhow2_r = rzw1rhow2(24 * j  + 10);

                auto sig_yz_rzw1rhow2_i = rzw1rhow2(24 * j  + 11);

                auto sig_yz_rxw1rxw2_r = rxw1rxw2(24 * j  + 10);

                auto sig_yz_rxw1rxw2_i = rxw1rxw2(24 * j  + 11);

                auto sig_yz_rxw1ryw2_r = rxw1ryw2(24 * j  + 10);

                auto sig_yz_rxw1ryw2_i = rxw1ryw2(24 * j  + 11);

                auto sig_yz_rxw1rzw2_r = rxw1rzw2(24 * j  + 10);

                auto sig_yz_rxw1rzw2_i = rxw1rzw2(24 * j  + 11);

                auto sig_yz_ryw1rxw2_r = ryw1rxw2(24 * j  + 10);

                auto sig_yz_ryw1rxw2_i = ryw1rxw2(24 * j  + 11);

                auto sig_yz_ryw1ryw2_r = ryw1ryw2(24 * j  + 10);

                auto sig_yz_ryw1ryw2_i = ryw1ryw2(24 * j  + 11);

                auto sig_yz_ryw1rzw2_r = ryw1rzw2(24 * j  + 10);

                auto sig_yz_ryw1rzw2_i = ryw1rzw2(24 * j  + 11);

                auto sig_yz_rzw1rxw2_r = rzw1rxw2(24 * j  + 10);

                auto sig_yz_rzw1rxw2_i = rzw1rxw2(24 * j  + 11);

                auto sig_yz_rzw1ryw2_r = rzw1ryw2(24 * j  + 10);

                auto sig_yz_rzw1ryw2_i = rzw1ryw2(24 * j  + 11);

                auto sig_yz_rzw1rzw2_r = rzw1rzw2(24 * j  + 10);

                auto sig_yz_rzw1rzw2_i = rzw1rzw2(24 * j  + 11);


                // lamtau_xx
                
                auto lamtau_xx_rhow1rhow2_r = rhow1rhow2(24 * j  + 12);
            
                auto lamtau_xx_rhow1rhow2_i = rhow1rhow2(24 * j  + 13);

                auto lamtau_xx_rxw1rhow2_r = rxw1rhow2(24 * j  + 12);

                auto lamtau_xx_rxw1rhow2_i = rxw1rhow2(24 * j  + 13);

                auto lamtau_xx_ryw1rhow2_r = ryw1rhow2(24 * j  + 12);

                auto lamtau_xx_ryw1rhow2_i = ryw1rhow2(24 * j  + 13);

                auto lamtau_xx_rzw1rhow2_r = rzw1rhow2(24 * j  + 12);

                auto lamtau_xx_rzw1rhow2_i = rzw1rhow2(24 * j  + 13);

                auto lamtau_xx_rxw1rxw2_r = rxw1rxw2(24 * j  + 12);

                auto lamtau_xx_rxw1rxw2_i = rxw1rxw2(24 * j  + 13);

                auto lamtau_xx_rxw1ryw2_r = rxw1ryw2(24 * j  + 12);

                auto lamtau_xx_rxw1ryw2_i = rxw1ryw2(24 * j  + 13);

                auto lamtau_xx_rxw1rzw2_r = rxw1rzw2(24 * j  + 12);

                auto lamtau_xx_rxw1rzw2_i = rxw1rzw2(24 * j  + 13);

                auto lamtau_xx_ryw1rxw2_r = ryw1rxw2(24 * j  + 12);

                auto lamtau_xx_ryw1rxw2_i = ryw1rxw2(24 * j  + 13);

                auto lamtau_xx_ryw1ryw2_r = ryw1ryw2(24 * j  + 12);

                auto lamtau_xx_ryw1ryw2_i = ryw1ryw2(24 * j  + 13);

                auto lamtau_xx_ryw1rzw2_r = ryw1rzw2(24 * j  + 12);

                auto lamtau_xx_ryw1rzw2_i = ryw1rzw2(24 * j  + 13);

                auto lamtau_xx_rzw1rxw2_r = rzw1rxw2(24 * j  + 12);

                auto lamtau_xx_rzw1rxw2_i = rzw1rxw2(24 * j  + 13);

                auto lamtau_xx_rzw1ryw2_r = rzw1ryw2(24 * j  + 12);

                auto lamtau_xx_rzw1ryw2_i = rzw1ryw2(24 * j  + 13);

                auto lamtau_xx_rzw1rzw2_r = rzw1rzw2(24 * j  + 12);

                auto lamtau_xx_rzw1rzw2_i = rzw1rzw2(24 * j  + 13);


                // lamtau_yy

                auto lamtau_yy_rhow1rhow2_r = rhow1rhow2(24 * j  + 14);
 
                auto lamtau_yy_rhow1rhow2_i = rhow1rhow2(24 * j  + 15);

                auto lamtau_yy_rxw1rhow2_r = rxw1rhow2(24 * j  + 14);

                auto lamtau_yy_rxw1rhow2_i = rxw1rhow2(24 * j  + 15);

                auto lamtau_yy_ryw1rhow2_r = ryw1rhow2(24 * j  + 14);

                auto lamtau_yy_ryw1rhow2_i = ryw1rhow2(24 * j  + 15);

                auto lamtau_yy_rzw1rhow2_r = rzw1rhow2(24 * j  + 14);

                auto lamtau_yy_rzw1rhow2_i = rzw1rhow2(24 * j  + 15);

                auto lamtau_yy_rxw1rxw2_r = rxw1rxw2(24 * j  + 14);

                auto lamtau_yy_rxw1rxw2_i = rxw1rxw2(24 * j  + 15);

                auto lamtau_yy_rxw1ryw2_r = rxw1ryw2(24 * j  + 14);

                auto lamtau_yy_rxw1ryw2_i = rxw1ryw2(24 * j  + 15);

                auto lamtau_yy_rxw1rzw2_r = rxw1rzw2(24 * j  + 14);

                auto lamtau_yy_rxw1rzw2_i = rxw1rzw2(24 * j  + 15);

                auto lamtau_yy_ryw1rxw2_r = ryw1rxw2(24 * j  + 14);

                auto lamtau_yy_ryw1rxw2_i = ryw1rxw2(24 * j  + 15);

                auto lamtau_yy_ryw1ryw2_r = ryw1ryw2(24 * j  + 14);

                auto lamtau_yy_ryw1ryw2_i = ryw1ryw2(24 * j  + 15);

                auto lamtau_yy_ryw1rzw2_r = ryw1rzw2(24 * j  + 14);

                auto lamtau_yy_ryw1rzw2_i = ryw1rzw2(24 * j  + 15);

                auto lamtau_yy_rzw1rxw2_r = rzw1rxw2(24 * j  + 14);

                auto lamtau_yy_rzw1rxw2_i = rzw1rxw2(24 * j  + 15);

                auto lamtau_yy_rzw1ryw2_r = rzw1ryw2(24 * j  + 14);

                auto lamtau_yy_rzw1ryw2_i = rzw1ryw2(24 * j  + 15);

                auto lamtau_yy_rzw1rzw2_r = rzw1rzw2(24 * j  + 14);

                auto lamtau_yy_rzw1rzw2_i = rzw1rzw2(24 * j  + 15);

                // lamtau_zz

                auto lamtau_zz_rhow1rhow2_r = rhow1rhow2(24 * j  + 16);
 
                auto lamtau_zz_rhow1rhow2_i = rhow1rhow2(24 * j  + 17);

                auto lamtau_zz_rxw1rhow2_r = rxw1rhow2(24 * j  + 16);

                auto lamtau_zz_rxw1rhow2_i = rxw1rhow2(24 * j  + 17);

                auto lamtau_zz_ryw1rhow2_r = ryw1rhow2(24 * j  + 16);

                auto lamtau_zz_ryw1rhow2_i = ryw1rhow2(24 * j  + 17);

                auto lamtau_zz_rzw1rhow2_r = rzw1rhow2(24 * j  + 16);

                auto lamtau_zz_rzw1rhow2_i = rzw1rhow2(24 * j  + 17);

                auto lamtau_zz_rxw1rxw2_r = rxw1rxw2(24 * j  + 16);

                auto lamtau_zz_rxw1rxw2_i = rxw1rxw2(24 * j  + 17);

                auto lamtau_zz_rxw1ryw2_r = rxw1ryw2(24 * j  + 16);

                auto lamtau_zz_rxw1ryw2_i = rxw1ryw2(24 * j  + 17);

                auto lamtau_zz_rxw1rzw2_r = rxw1rzw2(24 * j  + 16);

                auto lamtau_zz_rxw1rzw2_i = rxw1rzw2(24 * j  + 17);

                auto lamtau_zz_ryw1rxw2_r = ryw1rxw2(24 * j  + 16);

                auto lamtau_zz_ryw1rxw2_i = ryw1rxw2(24 * j  + 17);

                auto lamtau_zz_ryw1ryw2_r = ryw1ryw2(24 * j  + 16);

                auto lamtau_zz_ryw1ryw2_i = ryw1ryw2(24 * j  + 17);

                auto lamtau_zz_ryw1rzw2_r = ryw1rzw2(24 * j  + 16);

                auto lamtau_zz_ryw1rzw2_i = ryw1rzw2(24 * j  + 17);

                auto lamtau_zz_rzw1rxw2_r = rzw1rxw2(24 * j  + 16);

                auto lamtau_zz_rzw1rxw2_i = rzw1rxw2(24 * j  + 17);

                auto lamtau_zz_rzw1ryw2_r = rzw1ryw2(24 * j  + 16);

                auto lamtau_zz_rzw1ryw2_i = rzw1ryw2(24 * j  + 17);

                auto lamtau_zz_rzw1rzw2_r = rzw1rzw2(24 * j  + 16);

                auto lamtau_zz_rzw1rzw2_i = rzw1rzw2(24 * j  + 17);

                // lamtau_xy

                auto lamtau_xy_rhow1rhow2_r = rhow1rhow2(24 * j  + 18);

                auto lamtau_xy_rhow1rhow2_i = rhow1rhow2(24 * j  + 19);
 
                auto lamtau_xy_rxw1rhow2_r = rxw1rhow2(24 * j  + 18);

                auto lamtau_xy_rxw1rhow2_i = rxw1rhow2(24 * j  + 19);

                auto lamtau_xy_ryw1rhow2_r = ryw1rhow2(24 * j  + 18);

                auto lamtau_xy_ryw1rhow2_i = ryw1rhow2(24 * j  + 19);

                auto lamtau_xy_rzw1rhow2_r = rzw1rhow2(24 * j  + 18);

                auto lamtau_xy_rzw1rhow2_i = rzw1rhow2(24 * j  + 19);

                auto lamtau_xy_rxw1rxw2_r = rxw1rxw2(24 * j  + 18);

                auto lamtau_xy_rxw1rxw2_i = rxw1rxw2(24 * j  + 19);

                auto lamtau_xy_rxw1ryw2_r = rxw1ryw2(24 * j  + 18);

                auto lamtau_xy_rxw1ryw2_i = rxw1ryw2(24 * j  + 19);

                auto lamtau_xy_rxw1rzw2_r = rxw1rzw2(24 * j  + 18);

                auto lamtau_xy_rxw1rzw2_i = rxw1rzw2(24 * j  + 19);

                auto lamtau_xy_ryw1rxw2_r = ryw1rxw2(24 * j  + 18);

                auto lamtau_xy_ryw1rxw2_i = ryw1rxw2(24 * j  + 19);

                auto lamtau_xy_ryw1ryw2_r = ryw1ryw2(24 * j  + 18);

                auto lamtau_xy_ryw1ryw2_i = ryw1ryw2(24 * j  + 19);

                auto lamtau_xy_ryw1rzw2_r = ryw1rzw2(24 * j  + 18);

                auto lamtau_xy_ryw1rzw2_i = ryw1rzw2(24 * j  + 19);

                auto lamtau_xy_rzw1rxw2_r = rzw1rxw2(24 * j  + 18);

                auto lamtau_xy_rzw1rxw2_i = rzw1rxw2(24 * j  + 19);

                auto lamtau_xy_rzw1ryw2_r = rzw1ryw2(24 * j  + 18);

                auto lamtau_xy_rzw1ryw2_i = rzw1ryw2(24 * j  + 19);

                auto lamtau_xy_rzw1rzw2_r = rzw1rzw2(24 * j  + 18);

                auto lamtau_xy_rzw1rzw2_i = rzw1rzw2(24 * j  + 19);

                // lamtau_xz

                auto lamtau_xz_rhow1rhow2_r = rhow1rhow2(24 * j  + 20);

                auto lamtau_xz_rhow1rhow2_i = rhow1rhow2(24 * j  + 21);

                auto lamtau_xz_rxw1rhow2_r = rxw1rhow2(24 * j  + 20);

                auto lamtau_xz_rxw1rhow2_i = rxw1rhow2(24 * j  + 21);

                auto lamtau_xz_ryw1rhow2_r = ryw1rhow2(24 * j  + 20);

                auto lamtau_xz_ryw1rhow2_i = ryw1rhow2(24 * j  + 21);

                auto lamtau_xz_rzw1rhow2_r = rzw1rhow2(24 * j  + 20);

                auto lamtau_xz_rzw1rhow2_i = rzw1rhow2(24 * j  + 21);

                auto lamtau_xz_rxw1rxw2_r = rxw1rxw2(24 * j  + 20);

                auto lamtau_xz_rxw1rxw2_i = rxw1rxw2(24 * j  + 21);

                auto lamtau_xz_rxw1ryw2_r = rxw1ryw2(24 * j  + 20);

                auto lamtau_xz_rxw1ryw2_i = rxw1ryw2(24 * j  + 21);

                auto lamtau_xz_rxw1rzw2_r = rxw1rzw2(24 * j  + 20);

                auto lamtau_xz_rxw1rzw2_i = rxw1rzw2(24 * j  + 21);

                auto lamtau_xz_ryw1rxw2_r = ryw1rxw2(24 * j  + 20);

                auto lamtau_xz_ryw1rxw2_i = ryw1rxw2(24 * j  + 21);

                auto lamtau_xz_ryw1ryw2_r = ryw1ryw2(24 * j  + 20);

                auto lamtau_xz_ryw1ryw2_i = ryw1ryw2(24 * j  + 21);

                auto lamtau_xz_ryw1rzw2_r = ryw1rzw2(24 * j  + 20);

                auto lamtau_xz_ryw1rzw2_i = ryw1rzw2(24 * j  + 21);

                auto lamtau_xz_rzw1rxw2_r = rzw1rxw2(24 * j  + 20);

                auto lamtau_xz_rzw1rxw2_i = rzw1rxw2(24 * j  + 21);

                auto lamtau_xz_rzw1ryw2_r = rzw1ryw2(24 * j  + 20);

                auto lamtau_xz_rzw1ryw2_i = rzw1ryw2(24 * j  + 21);

                auto lamtau_xz_rzw1rzw2_r = rzw1rzw2(24 * j  + 20);

                auto lamtau_xz_rzw1rzw2_i = rzw1rzw2(24 * j  + 21);

                // lamtau_yz

                auto lamtau_yz_rhow1rhow2_r = rhow1rhow2(24 * j  + 22);
                
                auto lamtau_yz_rhow1rhow2_i = rhow1rhow2(24 * j  + 23);
 
                auto lamtau_yz_rxw1rhow2_r = rxw1rhow2(24 * j  + 22);

                auto lamtau_yz_rxw1rhow2_i = rxw1rhow2(24 * j  + 23);

                auto lamtau_yz_ryw1rhow2_r = ryw1rhow2(24 * j  + 22);

                auto lamtau_yz_ryw1rhow2_i = ryw1rhow2(24 * j  + 23);

                auto lamtau_yz_rzw1rhow2_r = rzw1rhow2(24 * j  + 22);

                auto lamtau_yz_rzw1rhow2_i = rzw1rhow2(24 * j  + 23);

                auto lamtau_yz_rxw1rxw2_r = rxw1rxw2(24 * j  + 22);

                auto lamtau_yz_rxw1rxw2_i = rxw1rxw2(24 * j  + 23);

                auto lamtau_yz_rxw1ryw2_r = rxw1ryw2(24 * j  + 22);

                auto lamtau_yz_rxw1ryw2_i = rxw1ryw2(24 * j  + 23);

                auto lamtau_yz_rxw1rzw2_r = rxw1rzw2(24 * j  + 22);

                auto lamtau_yz_rxw1rzw2_i = rxw1rzw2(24 * j  + 23);

                auto lamtau_yz_ryw1rxw2_r = ryw1rxw2(24 * j  + 22);

                auto lamtau_yz_ryw1rxw2_i = ryw1rxw2(24 * j  + 23);

                auto lamtau_yz_ryw1ryw2_r = ryw1ryw2(24 * j  + 22);

                auto lamtau_yz_ryw1ryw2_i = ryw1ryw2(24 * j  + 23);

                auto lamtau_yz_ryw1rzw2_r = ryw1rzw2(24 * j  + 22);

                auto lamtau_yz_ryw1rzw2_i = ryw1rzw2(24 * j  + 23);

                auto lamtau_yz_rzw1rxw2_r = rzw1rxw2(24 * j  + 22);

                auto lamtau_yz_rzw1rxw2_i = rzw1rxw2(24 * j  + 23);

                auto lamtau_yz_rzw1ryw2_r = rzw1ryw2(24 * j  + 22);

                auto lamtau_yz_rzw1ryw2_i = rzw1ryw2(24 * j  + 23);

                auto lamtau_yz_rzw1rzw2_r = rzw1rzw2(24 * j  + 22);

                auto lamtau_yz_rzw1rzw2_i = rzw1rzw2(24 * j  + 23);


                // Perturbed densities

                auto rhoBx_rho_r = rwDensityGrid.alphaDensity(12 * j);

                auto rhoBx_grad_x_r = rwDensityGrid.alphaDensityGradientX(12 * j);

                auto rhoBx_grad_y_r = rwDensityGrid.alphaDensityGradientY(12 * j);

                auto rhoBx_grad_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j);


                auto rhoBx_rho_i = rwDensityGrid.alphaDensity(12 * j + 1);

                auto rhoBx_grad_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 1);

                auto rhoBx_grad_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 1);

                auto rhoBx_grad_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 1);


                auto rhoBy_rho_r = rwDensityGrid.alphaDensity(12 * j + 2);

                auto rhoBy_grad_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 2);

                auto rhoBy_grad_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 2 );

                auto rhoBy_grad_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 2);


                auto rhoBy_rho_i = rwDensityGrid.alphaDensity(12 * j + 3);

                auto rhoBy_grad_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 3);

                auto rhoBy_grad_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 3);

                auto rhoBy_grad_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 3);


                auto rhoBz_rho_r = rwDensityGrid.alphaDensity(12 * j + 4);

                auto rhoBz_grad_x_r = rwDensityGrid.alphaDensityGradientX(12 * j+ 4);

                auto rhoBz_grad_y_r = rwDensityGrid.alphaDensityGradientY(12 * j+ 4);

                auto rhoBz_grad_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j+ 4);


                auto rhoBz_rho_i = rwDensityGrid.alphaDensity(12 * j + 5);

                auto rhoBz_grad_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 5);

                auto rhoBz_grad_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 5);

                auto rhoBz_grad_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 5);

                // C 

                auto rhoCx_rho_r = rwDensityGrid.alphaDensity(12 * j + 6);

                auto rhoCx_grad_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 6);

                auto rhoCx_grad_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 6);

                auto rhoCx_grad_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 6);

                auto rhoCx_rho_i = rwDensityGrid.alphaDensity(12 * j + 7);

                auto rhoCx_grad_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 7);

                auto rhoCx_grad_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 7);

                auto rhoCx_grad_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 7);


                auto rhoCy_rho_r = rwDensityGrid.alphaDensity(12 * j + 8);

                auto rhoCy_grad_x_r = rwDensityGrid.alphaDensityGradientX(12 * j + 8);

                auto rhoCy_grad_y_r = rwDensityGrid.alphaDensityGradientY(12 * j + 8);

                auto rhoCy_grad_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 8);

                auto rhoCy_rho_i = rwDensityGrid.alphaDensity(12 * j + 9);

                auto rhoCy_grad_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 9);

                auto rhoCy_grad_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 9);

                auto rhoCy_grad_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 9);


                auto rhoCz_rho_r = rwDensityGrid.alphaDensity(12 * j + 10);

                auto rhoCz_grad_x_r = rwDensityGrid.alphaDensityGradientX(12 * j  + 10);

                auto rhoCz_grad_y_r = rwDensityGrid.alphaDensityGradientY(12 * j  + 10);

                auto rhoCz_grad_z_r = rwDensityGrid.alphaDensityGradientZ(12 * j  + 10);

                auto rhoCz_rho_i = rwDensityGrid.alphaDensity(12 * j + 11);

                auto rhoCz_grad_x_i = rwDensityGrid.alphaDensityGradientX(12 * j + 11);

                auto rhoCz_grad_y_i = rwDensityGrid.alphaDensityGradientY(12 * j + 11);

                auto rhoCz_grad_z_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 11);
                

                double fac_sig_rhow1rhow2_r ; 


                double fac_sig_rhow1rhow2_i ; 


                double fac_sig_rxw1rhow2_r ;

                                        
                double fac_sig_rxw1rhow2_i ;
                                        

                double fac_sig_ryw1rhow2_r ;


                double fac_sig_ryw1rhow2_i ;


                double fac_sig_rzw1rhow2_r ;


                double fac_sig_rzw1rhow2_i ;
                                        

                double fac_sig_rxw1rxw2_r ;


                double fac_sig_rxw1rxw2_i ;


                double fac_sig_rxw1ryw2_r ;


                double fac_sig_rxw1ryw2_i ;


                double fac_sig_rxw1rzw2_r ;


                double fac_sig_rxw1rzw2_i ;

                            
                double fac_sig_ryw1rxw2_r ;
                                

                double fac_sig_ryw1rxw2_i ;


                double fac_sig_ryw1ryw2_r ;


                double fac_sig_ryw1ryw2_i ;


                double fac_sig_ryw1rzw2_r ;


                double fac_sig_ryw1rzw2_i ;


                double fac_sig_rzw1rxw2_r ;


                double fac_sig_rzw1rxw2_i ;


                double fac_sig_rzw1ryw2_r ;


                double fac_sig_rzw1ryw2_i ;


                double fac_sig_rzw1rzw2_r ;


                double fac_sig_rzw1rzw2_i ;


                double fac_lamtau_rhow1rhow2_r ; 


                double fac_lamtau_rhow1rhow2_i ; 


                double fac_lamtau_rxw1rhow2_r ;

                                        
                double fac_lamtau_rxw1rhow2_i ;
                                        

                double fac_lamtau_ryw1rhow2_r ;


                double fac_lamtau_ryw1rhow2_i ;


                double fac_lamtau_rzw1rhow2_r ;


                double fac_lamtau_rzw1rhow2_i ;
                                        

                double fac_lamtau_rxw1rxw2_r ;


                double fac_lamtau_rxw1rxw2_i ;


                double fac_lamtau_rxw1ryw2_r ;


                double fac_lamtau_rxw1ryw2_i ;


                double fac_lamtau_rxw1rzw2_r ;


                double fac_lamtau_rxw1rzw2_i ;

                            
                double fac_lamtau_ryw1rxw2_r ;
                                

                double fac_lamtau_ryw1rxw2_i ;


                double fac_lamtau_ryw1ryw2_r ;


                double fac_lamtau_ryw1ryw2_i ;


                double fac_lamtau_ryw1rzw2_r ;


                double fac_lamtau_ryw1rzw2_i ;


                double fac_lamtau_rzw1rxw2_r ;


                double fac_lamtau_rzw1rxw2_i ;


                double fac_lamtau_rzw1ryw2_r ;


                double fac_lamtau_rzw1ryw2_i ;


                double fac_lamtau_rzw1rzw2_r ;


                double fac_lamtau_rzw1rzw2_i ;

                for (int32_t i = 0; i < npoints; i++)
                {
                    

                    //  Sig_xx

                    double rhow1_r = rhoBx_rho_r[i];

                    double rxw1_r = rhoBx_grad_x_r[i];

                    double ryw1_r = rhoBx_grad_y_r[i];

                    double rzw1_r = rhoBx_grad_z_r[i];

                    double rhow1_i = rhoBx_rho_i[i];

                    double rxw1_i = rhoBx_grad_x_i[i];

                    double ryw1_i = rhoBx_grad_y_i[i];

                    double rzw1_i = rhoBx_grad_z_i[i];


                    double rhow2_r = rhoBx_rho_r[i];

                    double rxw2_r = rhoBx_grad_x_r[i];

                    double ryw2_r = rhoBx_grad_y_r[i];

                    double rzw2_r = rhoBx_grad_z_r[i];

                    double rhow2_i = rhoBx_rho_i[i];

                    double rxw2_i = rhoBx_grad_x_i[i];

                    double ryw2_i = rhoBx_grad_y_i[i];

                    double rzw2_i = rhoBx_grad_z_i[i];
                    

                    sig_xx_rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    sig_xx_rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    sig_xx_rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    sig_xx_rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    sig_xx_ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    sig_xx_ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    sig_xx_rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    sig_xx_rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    sig_xx_rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    sig_xx_rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    sig_xx_rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    sig_xx_rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    sig_xx_rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    sig_xx_rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    sig_xx_ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    sig_xx_ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    sig_xx_ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    sig_xx_ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    sig_xx_ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    sig_xx_ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    sig_xx_rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    sig_xx_rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    sig_xx_rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    sig_xx_rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    sig_xx_rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    sig_xx_rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;

                    //  sig_yy

                    rhow1_r = rhoBy_rho_r[i];

                    rxw1_r = rhoBy_grad_x_r[i];

                    ryw1_r = rhoBy_grad_y_r[i];

                    rzw1_r = rhoBy_grad_z_r[i];

                    rhow1_i = rhoBy_rho_i[i];

                    rxw1_i = rhoBy_grad_x_i[i];

                    ryw1_i = rhoBy_grad_y_i[i];

                    rzw1_i = rhoBy_grad_z_i[i];


                    rhow2_r = rhoBy_rho_r[i];

                    rxw2_r = rhoBy_grad_x_r[i];

                    ryw2_r = rhoBy_grad_y_r[i];

                    rzw2_r = rhoBy_grad_z_r[i];

                    rhow2_i = rhoBy_rho_i[i];

                    rxw2_i = rhoBy_grad_x_i[i];

                    ryw2_i = rhoBy_grad_y_i[i];

                    rzw2_i = rhoBy_grad_z_i[i];




                    sig_yy_rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);


                    sig_yy_rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    sig_yy_rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    sig_yy_rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    sig_yy_ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    sig_yy_ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    sig_yy_rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    sig_yy_rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    sig_yy_rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    sig_yy_rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    sig_yy_rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    sig_yy_rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    sig_yy_rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    sig_yy_rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    sig_yy_ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    sig_yy_ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    sig_yy_ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    sig_yy_ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    sig_yy_ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    sig_yy_ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    sig_yy_rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    sig_yy_rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    sig_yy_rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    sig_yy_rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    sig_yy_rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    sig_yy_rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    //  sig_zz

                    rhow1_r = rhoBz_rho_r[i];

                    rxw1_r = rhoBz_grad_x_r[i];

                    ryw1_r = rhoBz_grad_y_r[i];

                    rzw1_r = rhoBz_grad_z_r[i];

                    rhow1_i = rhoBz_rho_i[i];

                    rxw1_i = rhoBz_grad_x_i[i];

                    ryw1_i = rhoBz_grad_y_i[i];

                    rzw1_i = rhoBz_grad_z_i[i];


                    rhow2_r = rhoBz_rho_r[i];

                    rxw2_r = rhoBz_grad_x_r[i];

                    ryw2_r = rhoBz_grad_y_r[i];

                    rzw2_r = rhoBz_grad_z_r[i];

                    rhow2_i = rhoBz_rho_i[i];

                    rxw2_i = rhoBz_grad_x_i[i];

                    ryw2_i = rhoBz_grad_y_i[i];

                    rzw2_i = rhoBz_grad_z_i[i];



                    sig_zz_rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);


                    sig_zz_rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    sig_zz_rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    sig_zz_rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    sig_zz_ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    sig_zz_ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    sig_zz_rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    sig_zz_rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    sig_zz_rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    sig_zz_rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    sig_zz_rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    sig_zz_rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    sig_zz_rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    sig_zz_rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    sig_zz_ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    sig_zz_ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    sig_zz_ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    sig_zz_ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    sig_zz_ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    sig_zz_ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    sig_zz_rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    sig_zz_rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    sig_zz_rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    sig_zz_rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    sig_zz_rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    sig_zz_rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;



                    // Sigma factors

                    fac_sig_rhow1rhow2_r  =  (sig_xx_rhow1rhow2_r[i] +  sig_yy_rhow1rhow2_r[i] +  sig_zz_rhow1rhow2_r[i]) ;  
  
  
                    fac_sig_rhow1rhow2_i  =  (sig_xx_rhow1rhow2_i[i] +  sig_yy_rhow1rhow2_i[i] +  sig_zz_rhow1rhow2_i[i]) ;  
  
  
                    fac_sig_rxw1rhow2_r  = (sig_xx_rxw1rhow2_r[i] + sig_yy_rxw1rhow2_r[i] + sig_zz_rxw1rhow2_r[i] );
  
                                             
                    fac_sig_rxw1rhow2_i  =  (sig_xx_rxw1rhow2_i[i] + sig_yy_rxw1rhow2_i[i] + sig_zz_rxw1rhow2_i[i] );
                                           
  
                    fac_sig_ryw1rhow2_r  = (sig_xx_ryw1rhow2_r[i] + sig_yy_ryw1rhow2_r[i] + sig_zz_ryw1rhow2_r[i] );
  
  
                    fac_sig_ryw1rhow2_i  =  (sig_xx_ryw1rhow2_i[i] + sig_yy_ryw1rhow2_i[i] + sig_zz_ryw1rhow2_i[i] );
   
  
                    fac_sig_rzw1rhow2_r  = (sig_xx_rzw1rhow2_r[i] + sig_yy_rzw1rhow2_r[i] + sig_zz_rzw1rhow2_r[i] );
  
  
                    fac_sig_rzw1rhow2_i  =  (sig_xx_rzw1rhow2_i[i] + sig_yy_rzw1rhow2_i[i] + sig_zz_rzw1rhow2_i[i] );
                                             
  
                    fac_sig_rxw1rxw2_r  = ( sig_xx_rxw1rxw2_r[i] + sig_yy_rxw1rxw2_r[i] + sig_zz_rxw1rxw2_r[i] );
  
    
                    fac_sig_rxw1rxw2_i  = ( sig_xx_rxw1rxw2_i[i] + sig_yy_rxw1rxw2_i[i] + sig_zz_rxw1rxw2_i[i] );
                              
  
                    fac_sig_rxw1ryw2_r  = ( sig_xx_rxw1ryw2_r[i] + sig_yy_rxw1ryw2_r[i] + sig_zz_rxw1ryw2_r[i] );
  
  
                    fac_sig_rxw1ryw2_i  = ( sig_xx_rxw1ryw2_i[i] + sig_yy_rxw1ryw2_i[i] + sig_zz_rxw1ryw2_i[i] );
  
  
                    fac_sig_rxw1rzw2_r  = (sig_xx_rxw1rzw2_r[i] + sig_yy_rxw1rzw2_r[i] + sig_zz_rxw1rzw2_r[i] );
  
  
                    fac_sig_rxw1rzw2_i  = (sig_xx_rxw1rzw2_i[i] + sig_yy_rxw1rzw2_i[i] + sig_zz_rxw1rzw2_i[i] );
  
                                  
                    fac_sig_ryw1rxw2_r  = (sig_xx_ryw1rxw2_r[i] + sig_yy_ryw1rxw2_r[i]  + sig_zz_ryw1rxw2_r[i]);
                                      
  
                    fac_sig_ryw1rxw2_i  = (sig_xx_ryw1rxw2_i[i] + sig_yy_ryw1rxw2_i[i]  + sig_zz_ryw1rxw2_i[i]);
  
  
                    fac_sig_ryw1ryw2_r  = (sig_xx_ryw1ryw2_r[i] + sig_yy_ryw1ryw2_r[i]  + sig_zz_ryw1ryw2_r[i]);

                                        
                    fac_sig_ryw1ryw2_i  = (sig_xx_ryw1ryw2_i[i] + sig_yy_ryw1ryw2_i[i]  + sig_zz_ryw1ryw2_i[i] );
  

                    fac_sig_ryw1rzw2_r  = (sig_xx_ryw1rzw2_r[i] + sig_yy_ryw1rzw2_r[i]  + sig_zz_ryw1rzw2_r[i]);
  

                    fac_sig_ryw1rzw2_i  = (sig_xx_ryw1rzw2_i[i] + sig_yy_ryw1rzw2_i[i]  + sig_zz_ryw1rzw2_i[i]);
  
  
                    fac_sig_rzw1rxw2_r  = (sig_xx_rzw1rxw2_r[i] + sig_yy_rzw1rxw2_r[i]  + sig_zz_rzw1rxw2_r[i]);
  
  
                    fac_sig_rzw1rxw2_i  = (sig_xx_rzw1rxw2_i[i] + sig_yy_rzw1rxw2_i[i]  + sig_zz_rzw1rxw2_i[i]);
  
  
                    fac_sig_rzw1ryw2_r  = (sig_xx_rzw1ryw2_r[i] + sig_yy_rzw1ryw2_r[i]  + sig_zz_rzw1ryw2_r[i]);
  
  
                    fac_sig_rzw1ryw2_i  = (sig_xx_rzw1ryw2_i[i] + sig_yy_rzw1ryw2_i[i]  + sig_zz_rzw1ryw2_i[i]);
  
  
                    fac_sig_rzw1rzw2_r = (sig_xx_rzw1rzw2_r[i] + sig_yy_rzw1rzw2_r[i]  + sig_zz_rzw1rzw2_r[i]);


                    fac_sig_rzw1rzw2_i = (sig_xx_rzw1rzw2_i[i] + sig_yy_rzw1rzw2_i[i]  + sig_zz_rzw1rzw2_i[i]);



                    sig_xx_rhow1rhow2_r[i] *= 6.0 ; 

                    sig_xx_rhow1rhow2_i[i] *= 6.0 ;

                    sig_xx_rxw1rhow2_r[i] *= 6.0 ; 

                    sig_xx_rxw1rhow2_i[i] *= 6.0 ; 

                    sig_xx_ryw1rhow2_r[i] *= 6.0 ; 

                    sig_xx_ryw1rhow2_i[i] *= 6.0 ; 

                    sig_xx_rzw1rhow2_r[i] *= 6.0 ; 

                    sig_xx_rzw1rhow2_i[i] *= 6.0 ; 

                    sig_xx_rxw1rxw2_r[i] *= 6.0 ; 

                    sig_xx_rxw1rxw2_i[i] *= 6.0 ; 

                    sig_xx_rxw1ryw2_r[i] *= 6.0 ;

                    sig_xx_rxw1ryw2_i[i] *= 6.0 ;

                    sig_xx_rxw1rzw2_r[i] *= 6.0 ; 

                    sig_xx_rxw1rzw2_i[i] *= 6.0 ; 

                    sig_xx_ryw1rxw2_r[i] *= 6.0 ; 

                    sig_xx_ryw1rxw2_i[i] *= 6.0 ; 

                    sig_xx_ryw1ryw2_r[i] *= 6.0 ; 

                    sig_xx_ryw1ryw2_i[i] *= 6.0 ; 

                    sig_xx_ryw1rzw2_r[i] *= 6.0 ; 

                    sig_xx_ryw1rzw2_i[i] *= 6.0 ; 

                    sig_xx_rzw1rxw2_r[i] *= 6.0 ; 

                    sig_xx_rzw1rxw2_i[i] *= 6.0 ; 

                    sig_xx_rzw1ryw2_r[i] *= 6.0 ; 

                    sig_xx_rzw1ryw2_i[i] *= 6.0 ; 

                    sig_xx_rzw1rzw2_r[i] *= 6.0 ; 

                    sig_xx_rzw1rzw2_i[i] *= 6.0 ; 


                    sig_yy_rhow1rhow2_r[i] *= 6.0 ; 

                    sig_yy_rhow1rhow2_i[i] *= 6.0 ;

                    sig_yy_rxw1rhow2_r[i] *= 6.0 ; 

                    sig_yy_rxw1rhow2_i[i] *= 6.0 ; 

                    sig_yy_ryw1rhow2_r[i] *= 6.0 ; 

                    sig_yy_ryw1rhow2_i[i] *= 6.0 ; 

                    sig_yy_rzw1rhow2_r[i] *= 6.0 ; 

                    sig_yy_rzw1rhow2_i[i] *= 6.0 ; 

                    sig_yy_rxw1rxw2_r[i] *= 6.0 ; 

                    sig_yy_rxw1rxw2_i[i] *= 6.0 ; 

                    sig_yy_rxw1ryw2_r[i] *= 6.0 ;

                    sig_yy_rxw1ryw2_i[i] *= 6.0 ;

                    sig_yy_rxw1rzw2_r[i] *= 6.0 ; 

                    sig_yy_rxw1rzw2_i[i] *= 6.0 ; 

                    sig_yy_ryw1rxw2_r[i] *= 6.0 ; 

                    sig_yy_ryw1rxw2_i[i] *= 6.0 ; 

                    sig_yy_ryw1ryw2_r[i] *= 6.0 ; 

                    sig_yy_ryw1ryw2_i[i] *= 6.0 ; 

                    sig_yy_ryw1rzw2_r[i] *= 6.0 ; 

                    sig_yy_ryw1rzw2_i[i] *= 6.0 ; 

                    sig_yy_rzw1rxw2_r[i] *= 6.0 ; 

                    sig_yy_rzw1rxw2_i[i] *= 6.0 ; 

                    sig_yy_rzw1ryw2_r[i] *= 6.0 ; 

                    sig_yy_rzw1ryw2_i[i] *= 6.0 ; 

                    sig_yy_rzw1rzw2_r[i] *= 6.0 ; 

                    sig_yy_rzw1rzw2_i[i] *= 6.0 ; 



                    sig_zz_rhow1rhow2_r[i] *= 6.0 ; 

                    sig_zz_rhow1rhow2_i[i] *= 6.0 ;

                    sig_zz_rxw1rhow2_r[i] *= 6.0 ; 

                    sig_zz_rxw1rhow2_i[i] *= 6.0 ; 

                    sig_zz_ryw1rhow2_r[i] *= 6.0 ; 

                    sig_zz_ryw1rhow2_i[i] *= 6.0 ; 

                    sig_zz_rzw1rhow2_r[i] *= 6.0 ; 

                    sig_zz_rzw1rhow2_i[i] *= 6.0 ; 

                    sig_zz_rxw1rxw2_r[i] *= 6.0 ; 

                    sig_zz_rxw1rxw2_i[i] *= 6.0 ; 

                    sig_zz_rxw1ryw2_r[i] *= 6.0 ;

                    sig_zz_rxw1ryw2_i[i] *= 6.0 ;

                    sig_zz_rxw1rzw2_r[i] *= 6.0 ; 

                    sig_zz_rxw1rzw2_i[i] *= 6.0 ; 

                    sig_zz_ryw1rxw2_r[i] *= 6.0 ; 

                    sig_zz_ryw1rxw2_i[i] *= 6.0 ; 

                    sig_zz_ryw1ryw2_r[i] *= 6.0 ; 

                    sig_zz_ryw1ryw2_i[i] *= 6.0 ; 

                    sig_zz_ryw1rzw2_r[i] *= 6.0 ; 

                    sig_zz_ryw1rzw2_i[i] *= 6.0 ; 

                    sig_zz_rzw1rxw2_r[i] *= 6.0 ; 

                    sig_zz_rzw1rxw2_i[i] *= 6.0 ; 

                    sig_zz_rzw1ryw2_r[i] *= 6.0 ; 

                    sig_zz_rzw1ryw2_i[i] *= 6.0 ; 

                    sig_zz_rzw1rzw2_r[i] *= 6.0 ; 

                    sig_zz_rzw1rzw2_i[i] *= 6.0 ; 


                    // Add factors 


                    sig_xx_rhow1rhow2_r[i] += 3.0 * fac_sig_rhow1rhow2_r;

                    sig_xx_rhow1rhow2_i[i] += 3.0 *  fac_sig_rhow1rhow2_i;

                    sig_xx_rxw1rhow2_r[i] += 3.0 * fac_sig_rxw1rhow2_r;

                    sig_xx_rxw1rhow2_i[i] += 3.0 * fac_sig_rxw1rhow2_i;

                    sig_xx_ryw1rhow2_r[i] += 3.0 * fac_sig_ryw1rhow2_r;

                    sig_xx_ryw1rhow2_i[i] += 3.0 * fac_sig_ryw1rhow2_i;

                    sig_xx_rzw1rhow2_r[i] += 3.0 * fac_sig_rzw1rhow2_r;

                    sig_xx_rzw1rhow2_i[i] += 3.0 * fac_sig_rzw1rhow2_i;

                    sig_xx_rxw1rxw2_r[i] += 3.0 * fac_sig_rxw1rxw2_r;

                    sig_xx_rxw1rxw2_i[i] += 3.0 * fac_sig_rxw1rxw2_i;

                    sig_xx_rxw1ryw2_r[i] += 3.0 * fac_sig_rxw1ryw2_r;

                    sig_xx_rxw1ryw2_i[i] += 3.0 * fac_sig_rxw1ryw2_i;

                    sig_xx_rxw1rzw2_r[i] += 3.0 * fac_sig_rxw1rzw2_r;

                    sig_xx_rxw1rzw2_i[i] += 3.0 * fac_sig_rxw1rzw2_i;

                    sig_xx_ryw1rxw2_r[i] += 3.0 * fac_sig_ryw1rxw2_r;

                    sig_xx_ryw1rxw2_i[i] += 3.0 * fac_sig_ryw1rxw2_i;

                    sig_xx_ryw1ryw2_r[i] += 3.0 * fac_sig_ryw1ryw2_r;

                    sig_xx_ryw1ryw2_i[i] += 3.0 * fac_sig_ryw1ryw2_i;

                    sig_xx_ryw1rzw2_r[i] += 3.0 * fac_sig_ryw1rzw2_r;

                    sig_xx_ryw1rzw2_i[i] += 3.0 * fac_sig_ryw1rzw2_i;

                    sig_xx_rzw1rxw2_r[i] += 3.0 * fac_sig_rzw1rxw2_r;

                    sig_xx_rzw1rxw2_i[i] += 3.0 * fac_sig_rzw1rxw2_i;

                    sig_xx_rzw1ryw2_r[i] += 3.0 * fac_sig_rzw1ryw2_r;

                    sig_xx_rzw1ryw2_i[i] += 3.0 * fac_sig_rzw1ryw2_i;

                    sig_xx_rzw1rzw2_r[i] += 3.0 * fac_sig_rzw1rzw2_r;

                    sig_xx_rzw1rzw2_i[i] += 3.0 * fac_sig_rzw1rzw2_i;


                    sig_yy_rhow1rhow2_r[i] += 3.0 * fac_sig_rhow1rhow2_r;

                    sig_yy_rhow1rhow2_i[i] += 3.0 * fac_sig_rhow1rhow2_i;

                    sig_yy_rxw1rhow2_r[i] += 3.0 * fac_sig_rxw1rhow2_r;

                    sig_yy_rxw1rhow2_i[i] += 3.0 * fac_sig_rxw1rhow2_i;

                    sig_yy_ryw1rhow2_r[i] += 3.0 * fac_sig_ryw1rhow2_r;

                    sig_yy_ryw1rhow2_i[i] += 3.0 * fac_sig_ryw1rhow2_i;

                    sig_yy_rzw1rhow2_r[i] += 3.0 * fac_sig_rzw1rhow2_r;

                    sig_yy_rzw1rhow2_i[i] += 3.0 * fac_sig_rzw1rhow2_i;

                    sig_yy_rxw1rxw2_r[i] += 3.0 * fac_sig_rxw1rxw2_r;

                    sig_yy_rxw1rxw2_i[i] += 3.0 * fac_sig_rxw1rxw2_i;

                    sig_yy_rxw1ryw2_r[i] += 3.0 * fac_sig_rxw1ryw2_r;

                    sig_yy_rxw1ryw2_i[i] += 3.0 * fac_sig_rxw1ryw2_i;

                    sig_yy_rxw1rzw2_r[i] += 3.0 * fac_sig_rxw1rzw2_r;

                    sig_yy_rxw1rzw2_i[i] += 3.0 * fac_sig_rxw1rzw2_i;

                    sig_yy_ryw1rxw2_r[i] += 3.0 * fac_sig_ryw1rxw2_r;

                    sig_yy_ryw1rxw2_i[i] += 3.0 * fac_sig_ryw1rxw2_i;

                    sig_yy_ryw1ryw2_r[i] += 3.0 * fac_sig_ryw1ryw2_r;

                    sig_yy_ryw1ryw2_i[i] += 3.0 * fac_sig_ryw1ryw2_i;

                    sig_yy_ryw1rzw2_r[i] += 3.0 * fac_sig_ryw1rzw2_r;

                    sig_yy_ryw1rzw2_i[i] += 3.0 * fac_sig_ryw1rzw2_i;

                    sig_yy_rzw1rxw2_r[i] += 3.0 * fac_sig_rzw1rxw2_r;

                    sig_yy_rzw1rxw2_i[i] += 3.0 * fac_sig_rzw1rxw2_i;

                    sig_yy_rzw1ryw2_r[i] += 3.0 * fac_sig_rzw1ryw2_r;

                    sig_yy_rzw1ryw2_i[i] += 3.0 * fac_sig_rzw1ryw2_i;

                    sig_yy_rzw1rzw2_r[i] += 3.0 * fac_sig_rzw1rzw2_r;

                    sig_yy_rzw1rzw2_i[i] += 3.0 * fac_sig_rzw1rzw2_i;


                    sig_zz_rhow1rhow2_r[i] += 3.0 * fac_sig_rhow1rhow2_r;

                    sig_zz_rhow1rhow2_i[i] += 3.0 * fac_sig_rhow1rhow2_i;

                    sig_zz_rxw1rhow2_r[i] += 3.0 * fac_sig_rxw1rhow2_r;

                    sig_zz_rxw1rhow2_i[i] += 3.0 * fac_sig_rxw1rhow2_i;

                    sig_zz_ryw1rhow2_r[i] += 3.0 * fac_sig_ryw1rhow2_r;

                    sig_zz_ryw1rhow2_i[i] += 3.0 * fac_sig_ryw1rhow2_i;

                    sig_zz_rzw1rhow2_r[i] += 3.0 * fac_sig_rzw1rhow2_r;

                    sig_zz_rzw1rhow2_i[i] += 3.0 * fac_sig_rzw1rhow2_i;

                    sig_zz_rxw1rxw2_r[i] += 3.0 * fac_sig_rxw1rxw2_r;

                    sig_zz_rxw1rxw2_i[i] += 3.0 * fac_sig_rxw1rxw2_i;

                    sig_zz_rxw1ryw2_r[i] += 3.0 * fac_sig_rxw1ryw2_r;

                    sig_zz_rxw1ryw2_i[i] += 3.0 * fac_sig_rxw1ryw2_i;

                    sig_zz_rxw1rzw2_r[i] += 3.0 * fac_sig_rxw1rzw2_r;

                    sig_zz_rxw1rzw2_i[i] += 3.0 * fac_sig_rxw1rzw2_i;

                    sig_zz_ryw1rxw2_r[i] += 3.0 * fac_sig_ryw1rxw2_r;

                    sig_zz_ryw1rxw2_i[i] += 3.0 * fac_sig_ryw1rxw2_i;

                    sig_zz_ryw1ryw2_r[i] += 3.0 * fac_sig_ryw1ryw2_r;

                    sig_zz_ryw1ryw2_i[i] += 3.0 * fac_sig_ryw1ryw2_i;

                    sig_zz_ryw1rzw2_r[i] += 3.0 * fac_sig_ryw1rzw2_r;

                    sig_zz_ryw1rzw2_i[i] += 3.0 * fac_sig_ryw1rzw2_i;

                    sig_zz_rzw1rxw2_r[i] += 3.0 * fac_sig_rzw1rxw2_r;

                    sig_zz_rzw1rxw2_i[i] += 3.0 * fac_sig_rzw1rxw2_i;

                    sig_zz_rzw1ryw2_r[i] += 3.0 * fac_sig_rzw1ryw2_r;

                    sig_zz_rzw1ryw2_i[i] += 3.0 * fac_sig_rzw1ryw2_i;

                    sig_zz_rzw1rzw2_r[i] += 3.0 * fac_sig_rzw1rzw2_r;

                    sig_zz_rzw1rzw2_i[i] += 3.0 * fac_sig_rzw1rzw2_i;

                                   

                    //  sig_xy

                    rhow1_r = rhoBx_rho_r[i];

                    rxw1_r = rhoBx_grad_x_r[i];

                    ryw1_r = rhoBx_grad_y_r[i];

                    rzw1_r = rhoBx_grad_z_r[i];

                    rhow1_i = rhoBx_rho_i[i];

                    rxw1_i = rhoBx_grad_x_i[i];

                    ryw1_i = rhoBx_grad_y_i[i];

                    rzw1_i = rhoBx_grad_z_i[i];


                    rhow2_r = rhoBy_rho_r[i];

                    rxw2_r = rhoBy_grad_x_r[i];

                    ryw2_r = rhoBy_grad_y_r[i];

                    rzw2_r = rhoBy_grad_z_r[i];

                    rhow2_i = rhoBy_rho_i[i];

                    rxw2_i = rhoBy_grad_x_i[i];

                    ryw2_i = rhoBy_grad_y_i[i];

                    rzw2_i = rhoBy_grad_z_i[i];


                    sig_xy_rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);


                    sig_xy_rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    sig_xy_rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    sig_xy_rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    sig_xy_ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    sig_xy_ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    sig_xy_rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    sig_xy_rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    sig_xy_rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    sig_xy_rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    sig_xy_rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    sig_xy_rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    sig_xy_rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    sig_xy_rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    sig_xy_ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    sig_xy_ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    sig_xy_ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    sig_xy_ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    sig_xy_ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    sig_xy_ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    sig_xy_rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    sig_xy_rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    sig_xy_rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    sig_xy_rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    sig_xy_rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    sig_xy_rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    sig_xy_rhow1rhow2_r[i] *= 6.0 ; 

                    sig_xy_rhow1rhow2_i[i] *= 6.0 ;

                    sig_xy_rxw1rhow2_r[i] *= 6.0 ; 

                    sig_xy_rxw1rhow2_i[i] *= 6.0 ; 

                    sig_xy_ryw1rhow2_r[i] *= 6.0 ; 

                    sig_xy_ryw1rhow2_i[i] *= 6.0 ; 

                    sig_xy_rzw1rhow2_r[i] *= 6.0 ; 

                    sig_xy_rzw1rhow2_i[i] *= 6.0 ; 

                    sig_xy_rxw1rxw2_r[i] *= 6.0 ; 

                    sig_xy_rxw1rxw2_i[i] *= 6.0 ; 

                    sig_xy_rxw1ryw2_r[i] *= 6.0 ;

                    sig_xy_rxw1ryw2_i[i] *= 6.0 ;

                    sig_xy_rxw1rzw2_r[i] *= 6.0 ; 

                    sig_xy_rxw1rzw2_i[i] *= 6.0 ; 

                    sig_xy_ryw1rxw2_r[i] *= 6.0 ; 

                    sig_xy_ryw1rxw2_i[i] *= 6.0 ; 

                    sig_xy_ryw1ryw2_r[i] *= 6.0 ; 

                    sig_xy_ryw1ryw2_i[i] *= 6.0 ; 

                    sig_xy_ryw1rzw2_r[i] *= 6.0 ; 

                    sig_xy_ryw1rzw2_i[i] *= 6.0 ; 

                    sig_xy_rzw1rxw2_r[i] *= 6.0 ; 

                    sig_xy_rzw1rxw2_i[i] *= 6.0 ; 

                    sig_xy_rzw1ryw2_r[i] *= 6.0 ; 

                    sig_xy_rzw1ryw2_i[i] *= 6.0 ; 

                    sig_xy_rzw1rzw2_r[i] *= 6.0 ; 

                    sig_xy_rzw1rzw2_i[i] *= 6.0 ; 



                    //  sig_xz

                    rhow1_r = rhoBx_rho_r[i];

                    rxw1_r = rhoBx_grad_x_r[i];

                    ryw1_r = rhoBx_grad_y_r[i];

                    rzw1_r = rhoBx_grad_z_r[i];

                    rhow1_i = rhoBx_rho_i[i];

                    rxw1_i = rhoBx_grad_x_i[i];

                    ryw1_i = rhoBx_grad_y_i[i];

                    rzw1_i = rhoBx_grad_z_i[i];


                    rhow2_r = rhoBz_rho_r[i];

                    rxw2_r = rhoBz_grad_x_r[i];

                    ryw2_r = rhoBz_grad_y_r[i];

                    rzw2_r = rhoBz_grad_z_r[i];

                    rhow2_i = rhoBz_rho_i[i];

                    rxw2_i = rhoBz_grad_x_i[i];

                    ryw2_i = rhoBz_grad_y_i[i];

                    rzw2_i = rhoBz_grad_z_i[i];




                    sig_xz_rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);


                    sig_xz_rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    sig_xz_rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    sig_xz_rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    sig_xz_ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    sig_xz_ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    sig_xz_rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    sig_xz_rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    sig_xz_rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    sig_xz_rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    sig_xz_rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    sig_xz_rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    sig_xz_rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    sig_xz_rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    sig_xz_ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    sig_xz_ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    sig_xz_ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    sig_xz_ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    sig_xz_ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    sig_xz_ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    sig_xz_rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    sig_xz_rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    sig_xz_rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    sig_xz_rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    sig_xz_rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    sig_xz_rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    sig_xz_rhow1rhow2_r[i] *= 6.0 ; 

                    sig_xz_rhow1rhow2_i[i] *= 6.0 ;

                    sig_xz_rxw1rhow2_r[i] *= 6.0 ; 

                    sig_xz_rxw1rhow2_i[i] *= 6.0 ; 

                    sig_xz_ryw1rhow2_r[i] *= 6.0 ; 

                    sig_xz_ryw1rhow2_i[i] *= 6.0 ; 

                    sig_xz_rzw1rhow2_r[i] *= 6.0 ; 

                    sig_xz_rzw1rhow2_i[i] *= 6.0 ; 

                    sig_xz_rxw1rxw2_r[i] *= 6.0 ; 

                    sig_xz_rxw1rxw2_i[i] *= 6.0 ; 

                    sig_xz_rxw1ryw2_r[i] *= 6.0 ;

                    sig_xz_rxw1ryw2_i[i] *= 6.0 ;

                    sig_xz_rxw1rzw2_r[i] *= 6.0 ; 

                    sig_xz_rxw1rzw2_i[i] *= 6.0 ; 

                    sig_xz_ryw1rxw2_r[i] *= 6.0 ; 

                    sig_xz_ryw1rxw2_i[i] *= 6.0 ; 

                    sig_xz_ryw1ryw2_r[i] *= 6.0 ; 

                    sig_xz_ryw1ryw2_i[i] *= 6.0 ; 

                    sig_xz_ryw1rzw2_r[i] *= 6.0 ; 

                    sig_xz_ryw1rzw2_i[i] *= 6.0 ; 

                    sig_xz_rzw1rxw2_r[i] *= 6.0 ; 

                    sig_xz_rzw1rxw2_i[i] *= 6.0 ; 

                    sig_xz_rzw1ryw2_r[i] *= 6.0 ; 

                    sig_xz_rzw1ryw2_i[i] *= 6.0 ; 

                    sig_xz_rzw1rzw2_r[i] *= 6.0 ; 

                    sig_xz_rzw1rzw2_i[i] *= 6.0 ; 



                    //  sig_yz

                    rhow1_r = rhoBy_rho_r[i];

                    rxw1_r = rhoBy_grad_x_r[i];

                    ryw1_r = rhoBy_grad_y_r[i];

                    rzw1_r = rhoBy_grad_z_r[i];

                    rhow1_i = rhoBy_rho_i[i];

                    rxw1_i = rhoBy_grad_x_i[i];

                    ryw1_i = rhoBy_grad_y_i[i];

                    rzw1_i = rhoBy_grad_z_i[i];


                    rhow2_r = rhoBz_rho_r[i];

                    rxw2_r = rhoBz_grad_x_r[i];

                    ryw2_r = rhoBz_grad_y_r[i];

                    rzw2_r = rhoBz_grad_z_r[i];

                    rhow2_i = rhoBz_rho_i[i];

                    rxw2_i = rhoBz_grad_x_i[i];

                    ryw2_i = rhoBz_grad_y_i[i];

                    rzw2_i = rhoBz_grad_z_i[i];




                    sig_yz_rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);


                    sig_yz_rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    sig_yz_rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    sig_yz_rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    sig_yz_ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    sig_yz_ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    sig_yz_rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    sig_yz_rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    sig_yz_rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    sig_yz_rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    sig_yz_rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    sig_yz_rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    sig_yz_rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    sig_yz_rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    sig_yz_ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    sig_yz_ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    sig_yz_ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    sig_yz_ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    sig_yz_ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    sig_yz_ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    sig_yz_rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    sig_yz_rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    sig_yz_rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    sig_yz_rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    sig_yz_rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    sig_yz_rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    sig_yz_rhow1rhow2_r[i] *= 6.0 ; 

                    sig_yz_rhow1rhow2_i[i] *= 6.0 ;

                    sig_yz_rxw1rhow2_r[i] *= 6.0 ; 

                    sig_yz_rxw1rhow2_i[i] *= 6.0 ; 

                    sig_yz_ryw1rhow2_r[i] *= 6.0 ; 

                    sig_yz_ryw1rhow2_i[i] *= 6.0 ; 

                    sig_yz_rzw1rhow2_r[i] *= 6.0 ; 

                    sig_yz_rzw1rhow2_i[i] *= 6.0 ; 

                    sig_yz_rxw1rxw2_r[i] *= 6.0 ; 

                    sig_yz_rxw1rxw2_i[i] *= 6.0 ; 

                    sig_yz_rxw1ryw2_r[i] *= 6.0 ;

                    sig_yz_rxw1ryw2_i[i] *= 6.0 ;

                    sig_yz_rxw1rzw2_r[i] *= 6.0 ; 

                    sig_yz_rxw1rzw2_i[i] *= 6.0 ; 

                    sig_yz_ryw1rxw2_r[i] *= 6.0 ; 

                    sig_yz_ryw1rxw2_i[i] *= 6.0 ; 

                    sig_yz_ryw1ryw2_r[i] *= 6.0 ; 

                    sig_yz_ryw1ryw2_i[i] *= 6.0 ; 

                    sig_yz_ryw1rzw2_r[i] *= 6.0 ; 

                    sig_yz_ryw1rzw2_i[i] *= 6.0 ; 

                    sig_yz_rzw1rxw2_r[i] *= 6.0 ; 

                    sig_yz_rzw1rxw2_i[i] *= 6.0 ; 

                    sig_yz_rzw1ryw2_r[i] *= 6.0 ; 

                    sig_yz_rzw1ryw2_i[i] *= 6.0 ; 

                    sig_yz_rzw1rzw2_r[i] *= 6.0 ; 

                    sig_yz_rzw1rzw2_i[i] *= 6.0 ; 


                    //  lamtau_xx

                    rhow1_r = rhoBx_rho_r[i];

                    rxw1_r = rhoBx_grad_x_r[i];

                    ryw1_r = rhoBx_grad_y_r[i];

                    rzw1_r = rhoBx_grad_z_r[i];

                    rhow1_i = rhoBx_rho_i[i];

                    rxw1_i = rhoBx_grad_x_i[i];

                    ryw1_i = rhoBx_grad_y_i[i];

                    rzw1_i = rhoBx_grad_z_i[i];


                    rhow2_r = rhoCx_rho_r[i];

                    rxw2_r = rhoCx_grad_x_r[i];

                    ryw2_r = rhoCx_grad_y_r[i];

                    rzw2_r = rhoCx_grad_z_r[i];

                    rhow2_i = rhoCx_rho_i[i];

                    rxw2_i = rhoCx_grad_x_i[i];

                    ryw2_i = rhoCx_grad_y_i[i];

                    rzw2_i = rhoCx_grad_z_i[i];



                    lamtau_xx_rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);


                    lamtau_xx_rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    lamtau_xx_rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    lamtau_xx_rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    lamtau_xx_ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    lamtau_xx_ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    lamtau_xx_rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    lamtau_xx_rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    lamtau_xx_rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    lamtau_xx_rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    lamtau_xx_rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    lamtau_xx_rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    lamtau_xx_rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    lamtau_xx_rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    lamtau_xx_ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    lamtau_xx_ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    lamtau_xx_ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    lamtau_xx_ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    lamtau_xx_ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    lamtau_xx_ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    lamtau_xx_rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    lamtau_xx_rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    lamtau_xx_rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    lamtau_xx_rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    lamtau_xx_rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    lamtau_xx_rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;



                    //  lamtau_yy

                    rhow1_r = rhoBy_rho_r[i];

                    rxw1_r = rhoBy_grad_x_r[i];

                    ryw1_r = rhoBy_grad_y_r[i];

                    rzw1_r = rhoBy_grad_z_r[i];

                    rhow1_i = rhoBy_rho_i[i];

                    rxw1_i = rhoBy_grad_x_i[i];

                    ryw1_i = rhoBy_grad_y_i[i];

                    rzw1_i = rhoBy_grad_z_i[i];


                    rhow2_r = rhoCy_rho_r[i];

                    rxw2_r = rhoCy_grad_x_r[i];

                    ryw2_r = rhoCy_grad_y_r[i];

                    rzw2_r = rhoCy_grad_z_r[i];

                    rhow2_i = rhoCy_rho_i[i];

                    rxw2_i = rhoCy_grad_x_i[i];

                    ryw2_i = rhoCy_grad_y_i[i];

                    rzw2_i = rhoCy_grad_z_i[i];




                    lamtau_yy_rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);


                    lamtau_yy_rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    lamtau_yy_rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    lamtau_yy_rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    lamtau_yy_ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    lamtau_yy_ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    lamtau_yy_rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    lamtau_yy_rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    lamtau_yy_rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    lamtau_yy_rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    lamtau_yy_rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    lamtau_yy_rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    lamtau_yy_rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    lamtau_yy_rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    lamtau_yy_ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    lamtau_yy_ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    lamtau_yy_ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    lamtau_yy_ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    lamtau_yy_ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    lamtau_yy_ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    lamtau_yy_rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    lamtau_yy_rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    lamtau_yy_rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    lamtau_yy_rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    lamtau_yy_rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    lamtau_yy_rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;



                    //  lamtau_zz

                    rhow1_r = rhoBz_rho_r[i];

                    rxw1_r = rhoBz_grad_x_r[i];

                    ryw1_r = rhoBz_grad_y_r[i];

                    rzw1_r = rhoBz_grad_z_r[i];

                    rhow1_i = rhoBz_rho_i[i];

                    rxw1_i = rhoBz_grad_x_i[i];

                    ryw1_i = rhoBz_grad_y_i[i];

                    rzw1_i = rhoBz_grad_z_i[i];


                    rhow2_r = rhoCz_rho_r[i];

                    rxw2_r = rhoCz_grad_x_r[i];

                    ryw2_r = rhoCz_grad_y_r[i];

                    rzw2_r = rhoCz_grad_z_r[i];

                    rhow2_i = rhoCz_rho_i[i];

                    rxw2_i = rhoCz_grad_x_i[i];

                    ryw2_i = rhoCz_grad_y_i[i];

                    rzw2_i = rhoCz_grad_z_i[i];



                    lamtau_zz_rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);


                    lamtau_zz_rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    lamtau_zz_rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    lamtau_zz_rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    lamtau_zz_ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    lamtau_zz_ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    lamtau_zz_rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    lamtau_zz_rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    lamtau_zz_rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    lamtau_zz_rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    lamtau_zz_rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    lamtau_zz_rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    lamtau_zz_rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    lamtau_zz_rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    lamtau_zz_ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    lamtau_zz_ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    lamtau_zz_ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    lamtau_zz_ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    lamtau_zz_ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    lamtau_zz_ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    lamtau_zz_rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    lamtau_zz_rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    lamtau_zz_rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    lamtau_zz_rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    lamtau_zz_rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    lamtau_zz_rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    // Sigma factors

                    fac_lamtau_rhow1rhow2_r  =  (lamtau_xx_rhow1rhow2_r[i] +  lamtau_yy_rhow1rhow2_r[i] +  lamtau_zz_rhow1rhow2_r[i]) ;  
  
  
                    fac_lamtau_rhow1rhow2_i   =  (lamtau_xx_rhow1rhow2_i[i] +  lamtau_yy_rhow1rhow2_i[i] +  lamtau_zz_rhow1rhow2_i[i]) ;  
  
  
                    fac_lamtau_rxw1rhow2_r  = (lamtau_xx_rxw1rhow2_r[i] + lamtau_yy_rxw1rhow2_r[i] + lamtau_zz_rxw1rhow2_r[i] );
  
                                             
                    fac_lamtau_rxw1rhow2_i  =  (lamtau_xx_rxw1rhow2_i[i] + lamtau_yy_rxw1rhow2_i[i] + lamtau_zz_rxw1rhow2_i[i] );
                                           
  
                    fac_lamtau_ryw1rhow2_r  = (lamtau_xx_ryw1rhow2_r[i] + lamtau_yy_ryw1rhow2_r[i] + lamtau_zz_ryw1rhow2_r[i] );
  
  
                    fac_lamtau_ryw1rhow2_i  =  (lamtau_xx_ryw1rhow2_i[i] + lamtau_yy_ryw1rhow2_i[i] + lamtau_zz_ryw1rhow2_i[i] );
   
  
                    fac_lamtau_rzw1rhow2_r  = (lamtau_xx_rzw1rhow2_r[i] + lamtau_yy_rzw1rhow2_r[i] + lamtau_zz_rzw1rhow2_r[i] );
  
  
                    fac_lamtau_rzw1rhow2_i  =  (lamtau_xx_rzw1rhow2_i[i] + lamtau_yy_rzw1rhow2_i[i] + lamtau_zz_rzw1rhow2_i[i] );
                                             
  
                    fac_lamtau_rxw1rxw2_r  = ( lamtau_xx_rxw1rxw2_r[i] + lamtau_yy_rxw1rxw2_r[i] + lamtau_zz_rxw1rxw2_r[i] );
  
    
                    fac_lamtau_rxw1rxw2_i  = ( lamtau_xx_rxw1rxw2_i[i] + lamtau_yy_rxw1rxw2_i[i] + lamtau_zz_rxw1rxw2_i[i] );
                              
  
                    fac_lamtau_rxw1ryw2_r  = ( lamtau_xx_rxw1ryw2_r[i] + lamtau_yy_rxw1ryw2_r[i] + lamtau_zz_rxw1ryw2_r[i] );
  
  
                    fac_lamtau_rxw1ryw2_i  = ( lamtau_xx_rxw1ryw2_i[i] + lamtau_yy_rxw1ryw2_i[i] + lamtau_zz_rxw1ryw2_i[i] );
  
  
                    fac_lamtau_rxw1rzw2_r  = (lamtau_xx_rxw1rzw2_r[i] + lamtau_yy_rxw1rzw2_r[i] + lamtau_zz_rxw1rzw2_r[i] );
  
  
                    fac_lamtau_rxw1rzw2_i  = (lamtau_xx_rxw1rzw2_i[i] + lamtau_yy_rxw1rzw2_i[i] + lamtau_zz_rxw1rzw2_i[i] );
  
                                  
                    fac_lamtau_ryw1rxw2_r  = (lamtau_xx_ryw1rxw2_r[i] + lamtau_yy_ryw1rxw2_r[i]  + lamtau_zz_ryw1rxw2_r[i]);
                                      
  
                    fac_lamtau_ryw1rxw2_i  = (lamtau_xx_ryw1rxw2_i[i] + lamtau_yy_ryw1rxw2_i[i]  + lamtau_zz_ryw1rxw2_i[i]);
  
  
                    fac_lamtau_ryw1ryw2_r  = (lamtau_xx_ryw1ryw2_r[i] + lamtau_yy_ryw1ryw2_r[i]  + lamtau_zz_ryw1ryw2_r[i]);

                                        
                    fac_lamtau_ryw1ryw2_i  = (lamtau_xx_ryw1ryw2_i[i] + lamtau_yy_ryw1ryw2_i[i]  + lamtau_zz_ryw1ryw2_i[i] );
  

                    fac_lamtau_ryw1rzw2_r  = (lamtau_xx_ryw1rzw2_r[i] + lamtau_yy_ryw1rzw2_r[i]  + lamtau_zz_ryw1rzw2_r[i]);
  

                    fac_lamtau_ryw1rzw2_i  = (lamtau_xx_ryw1rzw2_i[i] + lamtau_yy_ryw1rzw2_i[i]  + lamtau_zz_ryw1rzw2_i[i]);
  
  
                    fac_lamtau_rzw1rxw2_r  = (lamtau_xx_rzw1rxw2_r[i] + lamtau_yy_rzw1rxw2_r[i]  + lamtau_zz_rzw1rxw2_r[i]);
  
  
                    fac_lamtau_rzw1rxw2_i  = (lamtau_xx_rzw1rxw2_i[i] + lamtau_yy_rzw1rxw2_i[i]  + lamtau_zz_rzw1rxw2_i[i]);
  
  
                    fac_lamtau_rzw1ryw2_r  = (lamtau_xx_rzw1ryw2_r[i] + lamtau_yy_rzw1ryw2_r[i]  + lamtau_zz_rzw1ryw2_r[i]);
  
  
                    fac_lamtau_rzw1ryw2_i  = (lamtau_xx_rzw1ryw2_i[i] + lamtau_yy_rzw1ryw2_i[i]  + lamtau_zz_rzw1ryw2_i[i]);
  
  
                    fac_lamtau_rzw1rzw2_r = (lamtau_xx_rzw1rzw2_r[i] + lamtau_yy_rzw1rzw2_r[i]  + lamtau_zz_rzw1rzw2_r[i]);


                    fac_lamtau_rzw1rzw2_i = (lamtau_xx_rzw1rzw2_i[i] + lamtau_yy_rzw1rzw2_i[i]  + lamtau_zz_rzw1rzw2_i[i]);


                    lamtau_xx_rhow1rhow2_r[i] *= 12.0 ; 

                    lamtau_xx_rhow1rhow2_i[i] *= 12.0 ;

                    lamtau_xx_rxw1rhow2_r[i] *= 12.0 ; 

                    lamtau_xx_rxw1rhow2_i[i] *= 12.0 ; 

                    lamtau_xx_ryw1rhow2_r[i] *= 12.0 ; 

                    lamtau_xx_ryw1rhow2_i[i] *= 12.0 ; 

                    lamtau_xx_rzw1rhow2_r[i] *= 12.0 ; 

                    lamtau_xx_rzw1rhow2_i[i] *= 12.0 ; 

                    lamtau_xx_rxw1rxw2_r[i] *= 12.0 ; 

                    lamtau_xx_rxw1rxw2_i[i] *= 12.0 ; 

                    lamtau_xx_rxw1ryw2_r[i] *= 12.0 ;

                    lamtau_xx_rxw1ryw2_i[i] *= 12.0 ;

                    lamtau_xx_rxw1rzw2_r[i] *= 12.0 ; 

                    lamtau_xx_rxw1rzw2_i[i] *= 12.0 ; 

                    lamtau_xx_ryw1rxw2_r[i] *= 12.0 ; 

                    lamtau_xx_ryw1rxw2_i[i] *= 12.0 ; 

                    lamtau_xx_ryw1ryw2_r[i] *= 12.0 ; 

                    lamtau_xx_ryw1ryw2_i[i] *= 12.0 ; 

                    lamtau_xx_ryw1rzw2_r[i] *= 12.0 ; 

                    lamtau_xx_ryw1rzw2_i[i] *= 12.0 ; 

                    lamtau_xx_rzw1rxw2_r[i] *= 12.0 ; 

                    lamtau_xx_rzw1rxw2_i[i] *= 12.0 ; 

                    lamtau_xx_rzw1ryw2_r[i] *= 12.0 ; 

                    lamtau_xx_rzw1ryw2_i[i] *= 12.0 ; 

                    lamtau_xx_rzw1rzw2_r[i] *= 12.0 ; 

                    lamtau_xx_rzw1rzw2_i[i] *= 12.0 ; 



                    lamtau_yy_rhow1rhow2_r[i] *= 12.0 ; 

                    lamtau_yy_rhow1rhow2_i[i] *= 12.0 ;

                    lamtau_yy_rxw1rhow2_r[i] *= 12.0 ; 

                    lamtau_yy_rxw1rhow2_i[i] *= 12.0 ; 

                    lamtau_yy_ryw1rhow2_r[i] *= 12.0 ; 

                    lamtau_yy_ryw1rhow2_i[i] *= 12.0 ; 

                    lamtau_yy_rzw1rhow2_r[i] *= 12.0 ; 

                    lamtau_yy_rzw1rhow2_i[i] *= 12.0 ; 

                    lamtau_yy_rxw1rxw2_r[i] *= 12.0 ; 

                    lamtau_yy_rxw1rxw2_i[i] *= 12.0 ; 

                    lamtau_yy_rxw1ryw2_r[i] *= 12.0 ;

                    lamtau_yy_rxw1ryw2_i[i] *= 12.0 ;

                    lamtau_yy_rxw1rzw2_r[i] *= 12.0 ; 

                    lamtau_yy_rxw1rzw2_i[i] *= 12.0 ; 

                    lamtau_yy_ryw1rxw2_r[i] *= 12.0 ; 

                    lamtau_yy_ryw1rxw2_i[i] *= 12.0 ; 

                    lamtau_yy_ryw1ryw2_r[i] *= 12.0 ; 

                    lamtau_yy_ryw1ryw2_i[i] *= 12.0 ; 

                    lamtau_yy_ryw1rzw2_r[i] *= 12.0 ; 

                    lamtau_yy_ryw1rzw2_i[i] *= 12.0 ; 

                    lamtau_yy_rzw1rxw2_r[i] *= 12.0 ; 

                    lamtau_yy_rzw1rxw2_i[i] *= 12.0 ; 

                    lamtau_yy_rzw1ryw2_r[i] *= 12.0 ; 

                    lamtau_yy_rzw1ryw2_i[i] *= 12.0 ; 

                    lamtau_yy_rzw1rzw2_r[i] *= 12.0 ; 

                    lamtau_yy_rzw1rzw2_i[i] *= 12.0 ; 




                    lamtau_zz_rhow1rhow2_r[i] *= 12.0 ; 

                    lamtau_zz_rhow1rhow2_i[i] *= 12.0 ;

                    lamtau_zz_rxw1rhow2_r[i] *= 12.0 ; 

                    lamtau_zz_rxw1rhow2_i[i] *= 12.0 ; 

                    lamtau_zz_ryw1rhow2_r[i] *= 12.0 ; 

                    lamtau_zz_ryw1rhow2_i[i] *= 12.0 ; 

                    lamtau_zz_rzw1rhow2_r[i] *= 12.0 ; 

                    lamtau_zz_rzw1rhow2_i[i] *= 12.0 ; 

                    lamtau_zz_rxw1rxw2_r[i] *= 12.0 ; 

                    lamtau_zz_rxw1rxw2_i[i] *= 12.0 ; 

                    lamtau_zz_rxw1ryw2_r[i] *= 12.0 ;

                    lamtau_zz_rxw1ryw2_i[i] *= 12.0 ;

                    lamtau_zz_rxw1rzw2_r[i] *= 12.0 ; 

                    lamtau_zz_rxw1rzw2_i[i] *= 12.0 ; 

                    lamtau_zz_ryw1rxw2_r[i] *= 12.0 ; 

                    lamtau_zz_ryw1rxw2_i[i] *= 12.0 ; 

                    lamtau_zz_ryw1ryw2_r[i] *= 12.0 ; 

                    lamtau_zz_ryw1ryw2_i[i] *= 12.0 ; 

                    lamtau_zz_ryw1rzw2_r[i] *= 12.0 ; 

                    lamtau_zz_ryw1rzw2_i[i] *= 12.0 ; 

                    lamtau_zz_rzw1rxw2_r[i] *= 12.0 ; 

                    lamtau_zz_rzw1rxw2_i[i] *= 12.0 ; 

                    lamtau_zz_rzw1ryw2_r[i] *= 12.0 ; 

                    lamtau_zz_rzw1ryw2_i[i] *= 12.0 ; 

                    lamtau_zz_rzw1rzw2_r[i] *= 12.0 ; 

                    lamtau_zz_rzw1rzw2_i[i] *= 12.0 ; 





                    // Add factors 


                    lamtau_xx_rhow1rhow2_r[i] += 6.0 * fac_lamtau_rhow1rhow2_r;

                    lamtau_xx_rhow1rhow2_i[i] += 6.0 * fac_lamtau_rhow1rhow2_i;

                    lamtau_xx_rxw1rhow2_r[i] += 6.0 * fac_lamtau_rxw1rhow2_r;

                    lamtau_xx_rxw1rhow2_i[i] += 6.0 * fac_lamtau_rxw1rhow2_i;

                    lamtau_xx_ryw1rhow2_r[i] += 6.0 * fac_lamtau_ryw1rhow2_r;

                    lamtau_xx_ryw1rhow2_i[i] += 6.0 * fac_lamtau_ryw1rhow2_i;

                    lamtau_xx_rzw1rhow2_r[i] += 6.0 * fac_lamtau_rzw1rhow2_r;

                    lamtau_xx_rzw1rhow2_i[i] += 6.0 * fac_lamtau_rzw1rhow2_i;

                    lamtau_xx_rxw1rxw2_r[i] += 6.0 * fac_lamtau_rxw1rxw2_r;

                    lamtau_xx_rxw1rxw2_i[i] += 6.0 * fac_lamtau_rxw1rxw2_i;

                    lamtau_xx_rxw1ryw2_r[i] += 6.0 * fac_lamtau_rxw1ryw2_r;

                    lamtau_xx_rxw1ryw2_i[i] += 6.0 * fac_lamtau_rxw1ryw2_i;

                    lamtau_xx_rxw1rzw2_r[i] += 6.0 * fac_lamtau_rxw1rzw2_r;

                    lamtau_xx_rxw1rzw2_i[i] += 6.0 * fac_lamtau_rxw1rzw2_i;

                    lamtau_xx_ryw1rxw2_r[i] += 6.0 * fac_lamtau_ryw1rxw2_r;

                    lamtau_xx_ryw1rxw2_i[i] += 6.0 * fac_lamtau_ryw1rxw2_i;

                    lamtau_xx_ryw1ryw2_r[i] += 6.0 * fac_lamtau_ryw1ryw2_r;

                    lamtau_xx_ryw1ryw2_i[i] += 6.0 * fac_lamtau_ryw1ryw2_i;

                    lamtau_xx_ryw1rzw2_r[i] += 6.0 * fac_lamtau_ryw1rzw2_r;

                    lamtau_xx_ryw1rzw2_i[i] += 6.0 * fac_lamtau_ryw1rzw2_i;

                    lamtau_xx_rzw1rxw2_r[i] += 6.0 * fac_lamtau_rzw1rxw2_r;

                    lamtau_xx_rzw1rxw2_i[i] += 6.0 * fac_lamtau_rzw1rxw2_i;

                    lamtau_xx_rzw1ryw2_r[i] += 6.0 * fac_lamtau_rzw1ryw2_r;

                    lamtau_xx_rzw1ryw2_i[i] += 6.0 * fac_lamtau_rzw1ryw2_i;

                    lamtau_xx_rzw1rzw2_r[i] += 6.0 * fac_lamtau_rzw1rzw2_r;

                    lamtau_xx_rzw1rzw2_i[i] += 6.0 * fac_lamtau_rzw1rzw2_i;


                    lamtau_yy_rhow1rhow2_r[i] += 6.0 * fac_lamtau_rhow1rhow2_r;

                    lamtau_yy_rhow1rhow2_i[i] += 6.0 * fac_lamtau_rhow1rhow2_i;

                    lamtau_yy_rxw1rhow2_r[i] += 6.0 * fac_lamtau_rxw1rhow2_r;

                    lamtau_yy_rxw1rhow2_i[i] += 6.0 * fac_lamtau_rxw1rhow2_i;

                    lamtau_yy_ryw1rhow2_r[i] += 6.0 * fac_lamtau_ryw1rhow2_r;

                    lamtau_yy_ryw1rhow2_i[i] += 6.0 * fac_lamtau_ryw1rhow2_i;

                    lamtau_yy_rzw1rhow2_r[i] += 6.0 * fac_lamtau_rzw1rhow2_r;

                    lamtau_yy_rzw1rhow2_i[i] += 6.0 * fac_lamtau_rzw1rhow2_i;

                    lamtau_yy_rxw1rxw2_r[i] += 6.0 * fac_lamtau_rxw1rxw2_r;

                    lamtau_yy_rxw1rxw2_i[i] += 6.0 * fac_lamtau_rxw1rxw2_i;

                    lamtau_yy_rxw1ryw2_r[i] += 6.0 * fac_lamtau_rxw1ryw2_r;

                    lamtau_yy_rxw1ryw2_i[i] += 6.0 * fac_lamtau_rxw1ryw2_i;

                    lamtau_yy_rxw1rzw2_r[i] += 6.0 * fac_lamtau_rxw1rzw2_r;

                    lamtau_yy_rxw1rzw2_i[i] += 6.0 * fac_lamtau_rxw1rzw2_i;

                    lamtau_yy_ryw1rxw2_r[i] += 6.0 * fac_lamtau_ryw1rxw2_r;

                    lamtau_yy_ryw1rxw2_i[i] += 6.0 * fac_lamtau_ryw1rxw2_i;

                    lamtau_yy_ryw1ryw2_r[i] += 6.0 * fac_lamtau_ryw1ryw2_r;

                    lamtau_yy_ryw1ryw2_i[i] += 6.0 * fac_lamtau_ryw1ryw2_i;

                    lamtau_yy_ryw1rzw2_r[i] += 6.0 * fac_lamtau_ryw1rzw2_r;

                    lamtau_yy_ryw1rzw2_i[i] += 6.0 * fac_lamtau_ryw1rzw2_i;

                    lamtau_yy_rzw1rxw2_r[i] += 6.0 * fac_lamtau_rzw1rxw2_r;

                    lamtau_yy_rzw1rxw2_i[i] += 6.0 * fac_lamtau_rzw1rxw2_i;

                    lamtau_yy_rzw1ryw2_r[i] += 6.0 * fac_lamtau_rzw1ryw2_r;

                    lamtau_yy_rzw1ryw2_i[i] += 6.0 * fac_lamtau_rzw1ryw2_i;

                    lamtau_yy_rzw1rzw2_r[i] += 6.0 * fac_lamtau_rzw1rzw2_r;

                    lamtau_yy_rzw1rzw2_i[i] += 6.0 * fac_lamtau_rzw1rzw2_i;


                    lamtau_zz_rhow1rhow2_r[i] += 6.0 * fac_lamtau_rhow1rhow2_r;

                    lamtau_zz_rhow1rhow2_i[i] += 6.0 * fac_lamtau_rhow1rhow2_i;

                    lamtau_zz_rxw1rhow2_r[i] += 6.0 * fac_lamtau_rxw1rhow2_r;

                    lamtau_zz_rxw1rhow2_i[i] += 6.0 * fac_lamtau_rxw1rhow2_i;

                    lamtau_zz_ryw1rhow2_r[i] += 6.0 * fac_lamtau_ryw1rhow2_r;

                    lamtau_zz_ryw1rhow2_i[i] += 6.0 * fac_lamtau_ryw1rhow2_i;

                    lamtau_zz_rzw1rhow2_r[i] += 6.0 * fac_lamtau_rzw1rhow2_r;

                    lamtau_zz_rzw1rhow2_i[i] += 6.0 * fac_lamtau_rzw1rhow2_i;

                    lamtau_zz_rxw1rxw2_r[i] += 6.0 * fac_lamtau_rxw1rxw2_r;

                    lamtau_zz_rxw1rxw2_i[i] += 6.0 * fac_lamtau_rxw1rxw2_i;

                    lamtau_zz_rxw1ryw2_r[i] += 6.0 * fac_lamtau_rxw1ryw2_r;

                    lamtau_zz_rxw1ryw2_i[i] += 6.0 * fac_lamtau_rxw1ryw2_i;

                    lamtau_zz_rxw1rzw2_r[i] += 6.0 * fac_lamtau_rxw1rzw2_r;

                    lamtau_zz_rxw1rzw2_i[i] += 6.0 * fac_lamtau_rxw1rzw2_i;

                    lamtau_zz_ryw1rxw2_r[i] += 6.0 * fac_lamtau_ryw1rxw2_r;

                    lamtau_zz_ryw1rxw2_i[i] += 6.0 * fac_lamtau_ryw1rxw2_i;

                    lamtau_zz_ryw1ryw2_r[i] += 6.0 * fac_lamtau_ryw1ryw2_r;

                    lamtau_zz_ryw1ryw2_i[i] += 6.0 * fac_lamtau_ryw1ryw2_i;

                    lamtau_zz_ryw1rzw2_r[i] += 6.0 * fac_lamtau_ryw1rzw2_r;

                    lamtau_zz_ryw1rzw2_i[i] += 6.0 * fac_lamtau_ryw1rzw2_i;

                    lamtau_zz_rzw1rxw2_r[i] += 6.0 * fac_lamtau_rzw1rxw2_r;

                    lamtau_zz_rzw1rxw2_i[i] += 6.0 * fac_lamtau_rzw1rxw2_i;

                    lamtau_zz_rzw1ryw2_r[i] += 6.0 * fac_lamtau_rzw1ryw2_r;

                    lamtau_zz_rzw1ryw2_i[i] += 6.0 * fac_lamtau_rzw1ryw2_i;

                    lamtau_zz_rzw1rzw2_r[i] += 6.0 * fac_lamtau_rzw1rzw2_r;

                    lamtau_zz_rzw1rzw2_i[i] += 6.0 * fac_lamtau_rzw1rzw2_i;


                    //  lamtau_xy (1)

                    rhow1_r = rhoBx_rho_r[i];

                    rxw1_r = rhoBx_grad_x_r[i];

                    ryw1_r = rhoBx_grad_y_r[i];

                    rzw1_r = rhoBx_grad_z_r[i];

                    rhow1_i = rhoBx_rho_i[i];

                    rxw1_i = rhoBx_grad_x_i[i];

                    ryw1_i = rhoBx_grad_y_i[i];

                    rzw1_i = rhoBx_grad_z_i[i];


                    rhow2_r = rhoCy_rho_r[i];

                    rxw2_r = rhoCy_grad_x_r[i];

                    ryw2_r = rhoCy_grad_y_r[i];

                    rzw2_r = rhoCy_grad_z_r[i];

                    rhow2_i = rhoCy_rho_i[i];

                    rxw2_i = rhoCy_grad_x_i[i];

                    ryw2_i = rhoCy_grad_y_i[i];

                    rzw2_i = rhoCy_grad_z_i[i];



                    lamtau_xy_rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);


                    lamtau_xy_rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    lamtau_xy_rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    lamtau_xy_rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    lamtau_xy_ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    lamtau_xy_ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    lamtau_xy_rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    lamtau_xy_rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    lamtau_xy_rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    lamtau_xy_rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    lamtau_xy_rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    lamtau_xy_rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    lamtau_xy_rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    lamtau_xy_rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    lamtau_xy_ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    lamtau_xy_ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    lamtau_xy_ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    lamtau_xy_ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    lamtau_xy_ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    lamtau_xy_ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    lamtau_xy_rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    lamtau_xy_rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    lamtau_xy_rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    lamtau_xy_rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    lamtau_xy_rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    lamtau_xy_rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    //  lamtau_xy (2)

                    rhow1_r = rhoBy_rho_r[i];

                    rxw1_r = rhoBy_grad_x_r[i];

                    ryw1_r = rhoBy_grad_y_r[i];

                    rzw1_r = rhoBy_grad_z_r[i];

                    rhow1_i = rhoBy_rho_i[i];

                    rxw1_i = rhoBy_grad_x_i[i];

                    ryw1_i = rhoBy_grad_y_i[i];

                    rzw1_i = rhoBy_grad_z_i[i];


                    rhow2_r = rhoCx_rho_r[i];

                    rxw2_r = rhoCx_grad_x_r[i];

                    ryw2_r = rhoCx_grad_y_r[i];

                    rzw2_r = rhoCx_grad_z_r[i];

                    rhow2_i = rhoCx_rho_i[i];

                    rxw2_i = rhoCx_grad_x_i[i];

                    ryw2_i = rhoCx_grad_y_i[i];

                    rzw2_i = rhoCx_grad_z_i[i];


                    lamtau_xy_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);


                    lamtau_xy_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    lamtau_xy_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    lamtau_xy_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    lamtau_xy_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    lamtau_xy_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    lamtau_xy_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    lamtau_xy_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    lamtau_xy_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    lamtau_xy_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    lamtau_xy_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    lamtau_xy_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    lamtau_xy_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    lamtau_xy_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    lamtau_xy_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    lamtau_xy_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    lamtau_xy_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    lamtau_xy_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    lamtau_xy_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    lamtau_xy_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    lamtau_xy_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    lamtau_xy_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    lamtau_xy_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    lamtau_xy_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    lamtau_xy_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    lamtau_xy_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;



                    lamtau_xy_rhow1rhow2_r[i] *= 6.0 ; 

                    lamtau_xy_rhow1rhow2_i[i] *= 6.0 ;

                    lamtau_xy_rxw1rhow2_r[i] *= 6.0 ; 

                    lamtau_xy_rxw1rhow2_i[i] *= 6.0 ; 

                    lamtau_xy_ryw1rhow2_r[i] *= 6.0 ; 

                    lamtau_xy_ryw1rhow2_i[i] *= 6.0 ; 

                    lamtau_xy_rzw1rhow2_r[i] *= 6.0 ; 

                    lamtau_xy_rzw1rhow2_i[i] *= 6.0 ; 

                    lamtau_xy_rxw1rxw2_r[i] *= 6.0 ; 

                    lamtau_xy_rxw1rxw2_i[i] *= 6.0 ; 

                    lamtau_xy_rxw1ryw2_r[i] *= 6.0 ;

                    lamtau_xy_rxw1ryw2_i[i] *= 6.0 ;

                    lamtau_xy_rxw1rzw2_r[i] *= 6.0 ; 

                    lamtau_xy_rxw1rzw2_i[i] *= 6.0 ; 

                    lamtau_xy_ryw1rxw2_r[i] *= 6.0 ; 

                    lamtau_xy_ryw1rxw2_i[i] *= 6.0 ; 

                    lamtau_xy_ryw1ryw2_r[i] *= 6.0 ; 

                    lamtau_xy_ryw1ryw2_i[i] *= 6.0 ; 

                    lamtau_xy_ryw1rzw2_r[i] *= 6.0 ; 

                    lamtau_xy_ryw1rzw2_i[i] *= 6.0 ; 

                    lamtau_xy_rzw1rxw2_r[i] *= 6.0 ; 

                    lamtau_xy_rzw1rxw2_i[i] *= 6.0 ; 

                    lamtau_xy_rzw1ryw2_r[i] *= 6.0 ; 

                    lamtau_xy_rzw1ryw2_i[i] *= 6.0 ; 

                    lamtau_xy_rzw1rzw2_r[i] *= 6.0 ; 

                    lamtau_xy_rzw1rzw2_i[i] *= 6.0 ; 


                    //  lamtau_xz (1)

                    rhow1_r = rhoBx_rho_r[i];

                    rxw1_r = rhoBx_grad_x_r[i];

                    ryw1_r = rhoBx_grad_y_r[i];

                    rzw1_r = rhoBx_grad_z_r[i];

                    rhow1_i = rhoBx_rho_i[i];

                    rxw1_i = rhoBx_grad_x_i[i];

                    ryw1_i = rhoBx_grad_y_i[i];

                    rzw1_i = rhoBx_grad_z_i[i];


                    rhow2_r = rhoCz_rho_r[i];

                    rxw2_r = rhoCz_grad_x_r[i];

                    ryw2_r = rhoCz_grad_y_r[i];

                    rzw2_r = rhoCz_grad_z_r[i];

                    rhow2_i = rhoCz_rho_i[i];

                    rxw2_i = rhoCz_grad_x_i[i];

                    ryw2_i = rhoCz_grad_y_i[i];

                    rzw2_i = rhoCz_grad_z_i[i];



                    lamtau_xz_rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);


                    lamtau_xz_rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    lamtau_xz_rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    lamtau_xz_rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    lamtau_xz_ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    lamtau_xz_ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    lamtau_xz_rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    lamtau_xz_rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    lamtau_xz_rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    lamtau_xz_rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    lamtau_xz_rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    lamtau_xz_rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    lamtau_xz_rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    lamtau_xz_rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    lamtau_xz_ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    lamtau_xz_ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    lamtau_xz_ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    lamtau_xz_ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    lamtau_xz_ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    lamtau_xz_ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    lamtau_xz_rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    lamtau_xz_rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    lamtau_xz_rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    lamtau_xz_rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    lamtau_xz_rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    lamtau_xz_rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    //  lamtau_xz (2)

                    rhow1_r = rhoBz_rho_r[i];

                    rxw1_r = rhoBz_grad_x_r[i];

                    ryw1_r = rhoBz_grad_y_r[i];

                    rzw1_r = rhoBz_grad_z_r[i];

                    rhow1_i = rhoBz_rho_i[i];

                    rxw1_i = rhoBz_grad_x_i[i];

                    ryw1_i = rhoBz_grad_y_i[i];

                    rzw1_i = rhoBz_grad_z_i[i];


                    rhow2_r = rhoCx_rho_r[i];

                    rxw2_r = rhoCx_grad_x_r[i];

                    ryw2_r = rhoCx_grad_y_r[i];

                    rzw2_r = rhoCx_grad_z_r[i];

                    rhow2_i = rhoCx_rho_i[i];

                    rxw2_i = rhoCx_grad_x_i[i];

                    ryw2_i = rhoCx_grad_y_i[i];

                    rzw2_i = rhoCx_grad_z_i[i];


                    lamtau_xz_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);


                    lamtau_xz_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    lamtau_xz_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    lamtau_xz_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    lamtau_xz_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    lamtau_xz_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    lamtau_xz_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    lamtau_xz_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    lamtau_xz_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    lamtau_xz_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    lamtau_xz_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    lamtau_xz_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    lamtau_xz_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    lamtau_xz_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    lamtau_xz_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    lamtau_xz_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    lamtau_xz_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    lamtau_xz_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    lamtau_xz_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    lamtau_xz_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    lamtau_xz_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    lamtau_xz_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    lamtau_xz_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    lamtau_xz_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    lamtau_xz_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    lamtau_xz_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    lamtau_xz_rhow1rhow2_r[i] *= 6.0 ; 

                    lamtau_xz_rhow1rhow2_i[i] *= 6.0 ;

                    lamtau_xz_rxw1rhow2_r[i] *= 6.0 ; 

                    lamtau_xz_rxw1rhow2_i[i] *= 6.0 ; 

                    lamtau_xz_ryw1rhow2_r[i] *= 6.0 ; 

                    lamtau_xz_ryw1rhow2_i[i] *= 6.0 ; 

                    lamtau_xz_rzw1rhow2_r[i] *= 6.0 ; 

                    lamtau_xz_rzw1rhow2_i[i] *= 6.0 ; 

                    lamtau_xz_rxw1rxw2_r[i] *= 6.0 ; 

                    lamtau_xz_rxw1rxw2_i[i] *= 6.0 ; 

                    lamtau_xz_rxw1ryw2_r[i] *= 6.0 ;

                    lamtau_xz_rxw1ryw2_i[i] *= 6.0 ;

                    lamtau_xz_rxw1rzw2_r[i] *= 6.0 ; 

                    lamtau_xz_rxw1rzw2_i[i] *= 6.0 ; 

                    lamtau_xz_ryw1rxw2_r[i] *= 6.0 ; 

                    lamtau_xz_ryw1rxw2_i[i] *= 6.0 ; 

                    lamtau_xz_ryw1ryw2_r[i] *= 6.0 ; 

                    lamtau_xz_ryw1ryw2_i[i] *= 6.0 ; 

                    lamtau_xz_ryw1rzw2_r[i] *= 6.0 ; 

                    lamtau_xz_ryw1rzw2_i[i] *= 6.0 ; 

                    lamtau_xz_rzw1rxw2_r[i] *= 6.0 ; 

                    lamtau_xz_rzw1rxw2_i[i] *= 6.0 ; 

                    lamtau_xz_rzw1ryw2_r[i] *= 6.0 ; 

                    lamtau_xz_rzw1ryw2_i[i] *= 6.0 ; 

                    lamtau_xz_rzw1rzw2_r[i] *= 6.0 ; 

                    lamtau_xz_rzw1rzw2_i[i] *= 6.0 ; 


                    //  lamtau_yz (1)

                    rhow1_r = rhoBy_rho_r[i];

                    rxw1_r = rhoBy_grad_x_r[i];

                    ryw1_r = rhoBy_grad_y_r[i];

                    rzw1_r = rhoBy_grad_z_r[i];

                    rhow1_i = rhoBy_rho_i[i];

                    rxw1_i = rhoBy_grad_x_i[i];

                    ryw1_i = rhoBy_grad_y_i[i];

                    rzw1_i = rhoBy_grad_z_i[i];


                    rhow2_r = rhoCz_rho_r[i];

                    rxw2_r = rhoCz_grad_x_r[i];

                    ryw2_r = rhoCz_grad_y_r[i];

                    rzw2_r = rhoCz_grad_z_r[i];

                    rhow2_i = rhoCz_rho_i[i];

                    rxw2_i = rhoCz_grad_x_i[i];

                    ryw2_i = rhoCz_grad_y_i[i];

                    rzw2_i = rhoCz_grad_z_i[i];



                    lamtau_yz_rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);


                    lamtau_yz_rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    lamtau_yz_rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    lamtau_yz_rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    lamtau_yz_ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    lamtau_yz_ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    lamtau_yz_rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    lamtau_yz_rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    lamtau_yz_rxw1rxw2_r[i] = rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    lamtau_yz_rxw1rxw2_i[i] = rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    lamtau_yz_rxw1ryw2_r[i] = rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    lamtau_yz_rxw1ryw2_i[i] = rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    lamtau_yz_rxw1rzw2_r[i] = rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    lamtau_yz_rxw1rzw2_i[i] = rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    lamtau_yz_ryw1rxw2_r[i] = ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    lamtau_yz_ryw1rxw2_i[i] = ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    lamtau_yz_ryw1ryw2_r[i] = ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    lamtau_yz_ryw1ryw2_i[i] = ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    lamtau_yz_ryw1rzw2_r[i] = ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    lamtau_yz_ryw1rzw2_i[i] = ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    lamtau_yz_rzw1rxw2_r[i] = rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    lamtau_yz_rzw1rxw2_i[i] = rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    lamtau_yz_rzw1ryw2_r[i] = rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    lamtau_yz_rzw1ryw2_i[i] = rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    lamtau_yz_rzw1rzw2_r[i] = rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    lamtau_yz_rzw1rzw2_i[i] = rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;



                    //  lamtau_yz (2)

                    rhow1_r = rhoBz_rho_r[i];

                    rxw1_r = rhoBz_grad_x_r[i];

                    ryw1_r = rhoBz_grad_y_r[i];

                    rzw1_r = rhoBz_grad_z_r[i];

                    rhow1_i = rhoBz_rho_i[i];

                    rxw1_i = rhoBz_grad_x_i[i];

                    ryw1_i = rhoBz_grad_y_i[i];

                    rzw1_i = rhoBz_grad_z_i[i];


                    rhow2_r = rhoCy_rho_r[i];

                    rxw2_r = rhoCy_grad_x_r[i];

                    ryw2_r = rhoCy_grad_y_r[i];

                    rzw2_r = rhoCy_grad_z_r[i];

                    rhow2_i = rhoCy_rho_i[i];

                    rxw2_i = rhoCy_grad_x_i[i];

                    ryw2_i = rhoCy_grad_y_i[i];

                    rzw2_i = rhoCy_grad_z_i[i];



                    lamtau_yz_rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);
                    

                    lamtau_yz_rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);


                    lamtau_yz_rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);


                    lamtau_yz_rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);


                    lamtau_yz_ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);


                    lamtau_yz_ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);


                    lamtau_yz_rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);


                    lamtau_yz_rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);


                    lamtau_yz_rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;


                    lamtau_yz_rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;


                    lamtau_yz_rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;


                    lamtau_yz_rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;


                    lamtau_yz_rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    lamtau_yz_rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    lamtau_yz_ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    lamtau_yz_ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    lamtau_yz_ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    lamtau_yz_ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    lamtau_yz_ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    lamtau_yz_ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    lamtau_yz_rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    lamtau_yz_rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    lamtau_yz_rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    lamtau_yz_rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    lamtau_yz_rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    lamtau_yz_rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;


                    lamtau_yz_rhow1rhow2_r[i] *= 6.0 ; 

                    lamtau_yz_rhow1rhow2_i[i] *= 6.0 ;

                    lamtau_yz_rxw1rhow2_r[i] *= 6.0 ; 

                    lamtau_yz_rxw1rhow2_i[i] *= 6.0 ; 

                    lamtau_yz_ryw1rhow2_r[i] *= 6.0 ; 

                    lamtau_yz_ryw1rhow2_i[i] *= 6.0 ; 

                    lamtau_yz_rzw1rhow2_r[i] *= 6.0 ; 

                    lamtau_yz_rzw1rhow2_i[i] *= 6.0 ; 

                    lamtau_yz_rxw1rxw2_r[i] *= 6.0 ; 

                    lamtau_yz_rxw1rxw2_i[i] *= 6.0 ; 

                    lamtau_yz_rxw1ryw2_r[i] *= 6.0 ;

                    lamtau_yz_rxw1ryw2_i[i] *= 6.0 ;

                    lamtau_yz_rxw1rzw2_r[i] *= 6.0 ; 

                    lamtau_yz_rxw1rzw2_i[i] *= 6.0 ; 

                    lamtau_yz_ryw1rxw2_r[i] *= 6.0 ; 

                    lamtau_yz_ryw1rxw2_i[i] *= 6.0 ; 

                    lamtau_yz_ryw1ryw2_r[i] *= 6.0 ; 

                    lamtau_yz_ryw1ryw2_i[i] *= 6.0 ; 

                    lamtau_yz_ryw1rzw2_r[i] *= 6.0 ; 

                    lamtau_yz_ryw1rzw2_i[i] *= 6.0 ; 

                    lamtau_yz_rzw1rxw2_r[i] *= 6.0 ; 

                    lamtau_yz_rzw1rxw2_i[i] *= 6.0 ; 

                    lamtau_yz_rzw1ryw2_r[i] *= 6.0 ; 

                    lamtau_yz_rzw1ryw2_i[i] *= 6.0 ; 

                    lamtau_yz_rzw1rzw2_r[i] *= 6.0 ; 

                    lamtau_yz_rzw1rzw2_i[i] *= 6.0 ; 


                }
            }
        }    
    if (fstr::upcase(quadMode) == "CRF_I")
        {
            // This routine is for computing the Fbc, Fbd, Fcd first-order fock matrices for the general cubic response function

        for (int32_t j = 0; j < numdens / 6; j++)
        {

                // B 
                auto rhowB_r = rwDensityGrid.alphaDensity(6 * j);

                auto gradwB_x_r = rwDensityGrid.alphaDensityGradientX(6 * j);

                auto gradwB_y_r = rwDensityGrid.alphaDensityGradientY(6 * j);

                auto gradwB_z_r = rwDensityGrid.alphaDensityGradientZ(6 * j);


                auto rhowB_i = rwDensityGrid.alphaDensity(6 * j + 1);

                auto gradwB_x_i = rwDensityGrid.alphaDensityGradientX(6 * j + 1);

                auto gradwB_y_i = rwDensityGrid.alphaDensityGradientY(6 * j + 1);

                auto gradwB_z_i = rwDensityGrid.alphaDensityGradientZ(6 * j + 1);

                // C

                auto rhowC_r = rwDensityGrid.alphaDensity(6 * j + 2);

                auto gradwC_x_r = rwDensityGrid.alphaDensityGradientX(6 * j + 2);

                auto gradwC_y_r = rwDensityGrid.alphaDensityGradientY(6 * j + 2);

                auto gradwC_z_r = rwDensityGrid.alphaDensityGradientZ(6 * j + 2);


                auto rhowC_i = rwDensityGrid.alphaDensity(6 * j + 3);

                auto gradwC_x_i = rwDensityGrid.alphaDensityGradientX(6 * j + 3);

                auto gradwC_y_i = rwDensityGrid.alphaDensityGradientY(6 * j + 3);

                auto gradwC_z_i = rwDensityGrid.alphaDensityGradientZ(6 * j + 3);

                // D 

                auto rhowD_r = rwDensityGrid.alphaDensity(6 * j + 4) ;

                auto gradwD_x_r = rwDensityGrid.alphaDensityGradientX(6 * j + 4);

                auto gradwD_y_r = rwDensityGrid.alphaDensityGradientY(6 * j + 4);

                auto gradwD_z_r = rwDensityGrid.alphaDensityGradientZ(6 * j + 4);


                auto rhowD_i = rwDensityGrid.alphaDensity(6 * j + 5);

                auto gradwD_x_i = rwDensityGrid.alphaDensityGradientX(6 * j + 5);

                auto gradwD_y_i = rwDensityGrid.alphaDensityGradientY(6 * j + 5);

                auto gradwD_z_i = rwDensityGrid.alphaDensityGradientZ(6 * j + 5);


                // Perturbed density products to be stored

                //BC

                auto rhowBrhowC_r= rhow1rhow2(6 * j);

                auto rhowBrhowC_i= rhow1rhow2(6 * j + 1);

                // rxwBrhowC

                auto rxwBrhowC_r = rxw1rhow2(6 * j);

                auto rxwBrhowC_i = rxw1rhow2(6 * j + 1);

                // rywBrhowC

                auto rywBrhowC_r = ryw1rhow2(6 * j);

                auto rywBrhowC_i = ryw1rhow2(6 * j + 1);

                // rzwBrhowC

                auto rzwBrhowC_r = rzw1rhow2(6 * j);

                auto rzwBrhowC_i = rzw1rhow2(6 * j + 1);

                // rAwBrBwC

                auto rxwBrxwC_r =rxw1rxw2(6 * j);

                auto rxwBrxwC_i =rxw1rxw2(6 * j + 1);

                auto rxwBrywC_r =rxw1ryw2(6 * j);

                auto rxwBrywC_i =rxw1ryw2(6 * j + 1);

                auto rxwBrzwC_r =rxw1rzw2(6 * j);

                auto rxwBrzwC_i =rxw1rzw2(6 * j + 1);

                auto rywBrxwC_r =ryw1rxw2(6 * j);

                auto rywBrxwC_i =ryw1rxw2(6 * j + 1);

                auto rywBrywC_r =ryw1ryw2(6 * j);

                auto rywBrywC_i =ryw1ryw2(6 * j + 1);

                auto rywBrzwC_r =ryw1rzw2(6 * j);

                auto rywBrzwC_i =ryw1rzw2(6 * j + 1);

                auto rzwBrxwC_r =rzw1rxw2(6 * j);

                auto rzwBrxwC_i =rzw1rxw2(6 * j + 1);

                auto rzwBrywC_r =rzw1ryw2(6 * j);

                auto rzwBrywC_i =rzw1ryw2(6 * j + 1);

                auto rzwBrzwC_r =rzw1rzw2(6 * j);

                auto rzwBrzwC_i =rzw1rzw2(6 * j + 1);


            // BD

                auto rhowBrhowD_r= rhow1rhow2(6 * j  + 2);

                auto rhowBrhowD_i= rhow1rhow2(6 *j  + 3);

                // rxwBrhowD

                auto rxwBrhowD_r = rxw1rhow2(6 * j  + 2);

                auto rxwBrhowD_i = rxw1rhow2(6 *j  + 3);

                // rywBrhowD

                auto rywBrhowD_r = ryw1rhow2(6 * j  + 2);

                auto rywBrhowD_i = ryw1rhow2(6 *j  + 3);

                // rzwBrhowD

                auto rzwBrhowD_r = rzw1rhow2(6 * j  + 2);

                auto rzwBrhowD_i = rzw1rhow2(6 *j  + 3);

                // rAwBrBwD

                auto rxwBrxwD_r =rxw1rxw2(6 * j  + 2);

                auto rxwBrxwD_i =rxw1rxw2(6 *j  + 3);

                auto rxwBrywD_r =rxw1ryw2(6 * j  + 2);

                auto rxwBrywD_i =rxw1ryw2(6 *j  + 3);

                auto rxwBrzwD_r =rxw1rzw2(6 * j  + 2);

                auto rxwBrzwD_i =rxw1rzw2(6 *j  + 3);

                auto rywBrxwD_r =ryw1rxw2(6 * j  + 2);

                auto rywBrxwD_i =ryw1rxw2(6 *j  + 3);

                auto rywBrywD_r =ryw1ryw2(6 * j  + 2);

                auto rywBrywD_i =ryw1ryw2(6 *j  + 3);

                auto rywBrzwD_r =ryw1rzw2(6 * j  + 2);

                auto rywBrzwD_i =ryw1rzw2(6 *j  + 3);

                auto rzwBrxwD_r =rzw1rxw2(6 * j  + 2);

                auto rzwBrxwD_i =rzw1rxw2(6 *j  + 3);

                auto rzwBrywD_r =rzw1ryw2(6 * j  + 2);

                auto rzwBrywD_i =rzw1ryw2(6 *j  + 3);

                auto rzwBrzwD_r =rzw1rzw2(6 * j  + 2);

                auto rzwBrzwD_i =rzw1rzw2(6 *j  + 3);
                

                // CD

                auto rhowCrhowD_r= rhow1rhow2(6 * j  + 4);

                auto rhowCrhowD_i= rhow1rhow2(6 *j  + 5);

                // rxwCrhowD

                auto rxwCrhowD_r = rxw1rhow2(6 * j  + 4);

                auto rxwCrhowD_i = rxw1rhow2(6 *j  + 5);

                // rywCrhowD

                auto rywCrhowD_r = ryw1rhow2(6 * j  + 4);

                auto rywCrhowD_i = ryw1rhow2(6 *j  + 5);

                // rzwCrhowD

                auto rzwCrhowD_r = rzw1rhow2(6 * j  + 4);

                auto rzwCrhowD_i = rzw1rhow2(6 *j  + 5);

                // rAwCrBwD

                auto rxwCrxwD_r =rxw1rxw2(6 * j  + 4);

                auto rxwCrxwD_i =rxw1rxw2(6 *j  + 5);

                auto rxwCrywD_r =rxw1ryw2(6 * j  + 4);

                auto rxwCrywD_i =rxw1ryw2(6 *j  + 5);

                auto rxwCrzwD_r =rxw1rzw2(6 * j  + 4);

                auto rxwCrzwD_i =rxw1rzw2(6 *j  + 5);

                auto rywCrxwD_r =ryw1rxw2(6 * j  + 4);

                auto rywCrxwD_i =ryw1rxw2(6 *j  + 5);

                auto rywCrywD_r =ryw1ryw2(6 * j  + 4);

                auto rywCrywD_i =ryw1ryw2(6 *j  + 5);

                auto rywCrzwD_r =ryw1rzw2(6 * j  + 4);

                auto rywCrzwD_i =ryw1rzw2(6 *j  + 5);

                auto rzwCrxwD_r =rzw1rxw2(6 * j  + 4);

                auto rzwCrxwD_i =rzw1rxw2(6 *j  + 5);

                auto rzwCrywD_r =rzw1ryw2(6 * j  + 4);

                auto rzwCrywD_i =rzw1ryw2(6 *j  + 5);

                auto rzwCrzwD_r =rzw1rzw2(6 * j  + 4);

                auto rzwCrzwD_i =rzw1rzw2(6 *j  + 5);


                for (int32_t i = 0; i < npoints; i++)
                {
                    // RW1 densities


                    // B 

                    double rxwB_r = gradwB_x_r[i];

                    double rywB_r = gradwB_y_r[i];

                    double rzwB_r = gradwB_z_r[i];

                    double rxwB_i = gradwB_x_i[i];

                    double rywB_i = gradwB_y_i[i];

                    double rzwB_i = gradwB_z_i[i];

                    // C

                    double rxwC_r = gradwC_x_r[i];

                    double rywC_r = gradwC_y_r[i];

                    double rzwC_r = gradwC_z_r[i];

                    double rxwC_i = gradwC_x_i[i];

                    double rywC_i = gradwC_y_i[i];

                    double rzwC_i = gradwC_z_i[i];


                    // D 
                    double rxwD_r = gradwD_x_r[i];

                    double rywD_r = gradwD_y_r[i];

                    double rzwD_r = gradwD_z_r[i];

                    double rxwD_i = gradwD_x_i[i];

                    double rywD_i = gradwD_y_i[i];

                    double rzwD_i = gradwD_z_i[i];


                    // BC

                    rhowBrhowC_r[i] = 2.0 * (rhowB_r[i] * rhowC_r[i] - rhowB_i[i] * rhowC_i[i]);

                    rhowBrhowC_i[i] = 2.0 * (rhowB_r[i] * rhowC_i[i] + rhowB_i[i] * rhowC_r[i]);

                    rxwBrhowC_r[i] = 2.0 * (rxwB_r * rhowC_r[i] - rxwB_i * rhowC_i[i]

                                            + rxwC_r * rhowB_r[i] - rxwC_i * rhowB_i[i]);

                    rxwBrhowC_i[i] = 2.0 * (rxwB_r * rhowC_i[i] + rxwB_i * rhowC_r[i]

                                            + rxwC_r * rhowB_i[i] + rxwC_i * rhowB_r[i]);

                    rywBrhowC_r[i] = 2.0 * (rywB_r * rhowC_r[i] - rywB_i * rhowC_i[i]

                                            + rywC_r * rhowB_r[i] - rywC_i * rhowB_i[i]);

                    rywBrhowC_i[i] = 2.0 * (rywB_r * rhowC_i[i] + rywB_i * rhowC_r[i]

                                            + rywC_r * rhowB_i[i] + rywC_i * rhowB_r[i]);

                    rzwBrhowC_r[i] = 2.0 * (rzwB_r * rhowC_r[i] - rzwB_i * rhowC_i[i]

                                            + rzwC_r * rhowB_r[i] - rzwC_i * rhowB_i[i]);

                    rzwBrhowC_i[i] = 2.0 * (rzwB_r * rhowC_i[i] + rzwB_i * rhowC_r[i]

                                            + rzwC_r * rhowB_i[i] + rzwC_i * rhowB_r[i]);

                    rxwBrxwC_r[i] = rxwB_r * rxwC_r - rxwB_i * rxwC_i

                                    + rxwC_r * rxwB_r - rxwC_i * rxwB_i;

                    rxwBrxwC_i[i] = rxwB_r * rxwC_i + rxwB_r * rxwC_i

                                    + rxwC_r * rxwB_i + rxwC_r * rxwB_i;

                    rxwBrywC_r[i] = rxwB_r * rywC_r - rxwB_i * rywC_i

                                    + rxwC_r * rywB_r - rxwC_i * rywB_i;

                    rxwBrywC_i[i] = rxwB_r * rywC_i + rxwB_r * rywC_i

                                    + rxwC_r * rywB_i + rxwC_r * rywB_i;

                    rxwBrzwC_r[i] = rxwB_r * rzwC_r - rxwB_i * rzwC_i

                                    + rxwC_r * rzwB_r - rxwC_i * rzwB_i;

                    rxwBrzwC_i[i] = rxwB_r * rzwC_i + rxwB_r * rzwC_i

                                    + rxwC_r * rzwB_i + rxwC_r * rzwB_i;

                    rywBrxwC_r[i] = rywB_r * rxwC_r - rywB_i * rxwC_i

                                    + rywC_r * rxwB_r - rywC_i * rxwB_i;

                    rywBrxwC_i[i] = rywB_r * rxwC_i + rywB_r * rxwC_i

                                    + rywC_r * rxwB_i + rywC_r * rxwB_i;

                    rywBrywC_r[i] = rywB_r * rywC_r - rywB_i * rywC_i

                                    + rywC_r * rywB_r - rywC_i * rywB_i;

                    rywBrywC_i[i] = rywB_r * rywC_i + rywB_r * rywC_i

                                    + rywC_r * rywB_i + rywC_r * rywB_i;

                    rywBrzwC_r[i] = rywB_r * rzwC_r - rywB_i * rzwC_i

                                    + rywC_r * rzwB_r - rywC_i * rzwB_i;

                    rywBrzwC_i[i] = rywB_r * rzwC_i + rywB_r * rzwC_i

                                    + rywC_r * rzwB_i + rywC_r * rzwB_i;

                    rzwBrxwC_r[i] = rzwB_r * rxwC_r - rzwB_i * rxwC_i

                                    + rzwC_r * rxwB_r - rzwC_i * rxwB_i;

                    rzwBrxwC_i[i] = rzwB_r * rxwC_i + rzwB_r * rxwC_i

                                    + rzwC_r * rxwB_i + rzwC_r * rxwB_i;

                    rzwBrywC_r[i] = rzwB_r * rywC_r - rzwB_i * rywC_i

                                    + rzwC_r * rywB_r - rzwC_i * rywB_i;

                    rzwBrywC_i[i] = rzwB_r * rywC_i + rzwB_r * rywC_i

                                    + rzwC_r * rywB_i + rzwC_r * rywB_i;

                    rzwBrzwC_r[i] = rzwB_r * rzwC_r - rzwB_i * rzwC_i

                                    + rzwC_r * rzwB_r - rzwC_i * rzwB_i;

                    rzwBrzwC_i[i] = rzwB_r * rzwC_i + rzwB_r * rzwC_i

                                    + rzwC_r * rzwB_i + rzwC_r * rzwB_i;


                     // BD

                    rhowBrhowD_r[i] = 2.0 * (rhowB_r[i] * rhowD_r[i] - rhowB_i[i] * rhowD_i[i]);

                    rhowBrhowD_i[i] = 2.0 * (rhowB_r[i] * rhowD_i[i] + rhowB_i[i] * rhowD_r[i]);

                    rxwBrhowD_r[i] = 2.0 * (rxwB_r * rhowD_r[i] - rxwB_i * rhowD_i[i]

                                            + rxwD_r * rhowB_r[i] - rxwD_i * rhowB_i[i]);

                    rxwBrhowD_i[i] = 2.0 * (rxwB_r * rhowD_i[i] + rxwB_i * rhowD_r[i]

                                            + rxwD_r * rhowB_i[i] + rxwD_i * rhowB_r[i]);

                    rywBrhowD_r[i] = 2.0 * (rywB_r * rhowD_r[i] - rywB_i * rhowD_i[i]

                                            + rywD_r * rhowB_r[i] - rywD_i * rhowB_i[i]);

                    rywBrhowD_i[i] = 2.0 * (rywB_r * rhowD_i[i] + rywB_i * rhowD_r[i]

                                            + rywD_r * rhowB_i[i] + rywD_i * rhowB_r[i]);

                    rzwBrhowD_r[i] = 2.0 * (rzwB_r * rhowD_r[i] - rzwB_i * rhowD_i[i]

                                            + rzwD_r * rhowB_r[i] - rzwD_i * rhowB_i[i]);

                    rzwBrhowD_i[i] = 2.0 * (rzwB_r * rhowD_i[i] + rzwB_i * rhowD_r[i]

                                            + rzwD_r * rhowB_i[i] + rzwD_i * rhowB_r[i]);

                    rxwBrxwD_r[i] = rxwB_r * rxwD_r - rxwB_i * rxwD_i

                                    + rxwD_r * rxwB_r - rxwD_i * rxwB_i;

                    rxwBrxwD_i[i] = rxwB_r * rxwD_i + rxwB_r * rxwD_i

                                    + rxwD_r * rxwB_i + rxwD_r * rxwB_i;

                    rxwBrywD_r[i] = rxwB_r * rywD_r - rxwB_i * rywD_i

                                    + rxwD_r * rywB_r - rxwD_i * rywB_i;

                    rxwBrywD_i[i] = rxwB_r * rywD_i + rxwB_r * rywD_i

                                    + rxwD_r * rywB_i + rxwD_r * rywB_i;

                    rxwBrzwD_r[i] = rxwB_r * rzwD_r - rxwB_i * rzwD_i

                                    + rxwD_r * rzwB_r - rxwD_i * rzwB_i;

                    rxwBrzwD_i[i] = rxwB_r * rzwD_i + rxwB_r * rzwD_i

                                    + rxwD_r * rzwB_i + rxwD_r * rzwB_i;

                    rywBrxwD_r[i] = rywB_r * rxwD_r - rywB_i * rxwD_i

                                    + rywD_r * rxwB_r - rywD_i * rxwB_i;

                    rywBrxwD_i[i] = rywB_r * rxwD_i + rywB_r * rxwD_i

                                    + rywD_r * rxwB_i + rywD_r * rxwB_i;

                    rywBrywD_r[i] = rywB_r * rywD_r - rywB_i * rywD_i

                                    + rywD_r * rywB_r - rywD_i * rywB_i;

                    rywBrywD_i[i] = rywB_r * rywD_i + rywB_r * rywD_i

                                    + rywD_r * rywB_i + rywD_r * rywB_i;

                    rywBrzwD_r[i] = rywB_r * rzwD_r - rywB_i * rzwD_i

                                    + rywD_r * rzwB_r - rywD_i * rzwB_i;

                    rywBrzwD_i[i] = rywB_r * rzwD_i + rywB_r * rzwD_i

                                    + rywD_r * rzwB_i + rywD_r * rzwB_i;

                    rzwBrxwD_r[i] = rzwB_r * rxwD_r - rzwB_i * rxwD_i

                                    + rzwD_r * rxwB_r - rzwD_i * rxwB_i;

                    rzwBrxwD_i[i] = rzwB_r * rxwD_i + rzwB_r * rxwD_i

                                    + rzwD_r * rxwB_i + rzwD_r * rxwB_i;

                    rzwBrywD_r[i] = rzwB_r * rywD_r - rzwB_i * rywD_i

                                    + rzwD_r * rywB_r - rzwD_i * rywB_i;

                    rzwBrywD_i[i] = rzwB_r * rywD_i + rzwB_r * rywD_i

                                    + rzwD_r * rywB_i + rzwD_r * rywB_i;

                    rzwBrzwD_r[i] = rzwB_r * rzwD_r - rzwB_i * rzwD_i

                                    + rzwD_r * rzwB_r - rzwD_i * rzwB_i;

                    rzwBrzwD_i[i] = rzwB_r * rzwD_i + rzwB_r * rzwD_i

                                    + rzwD_r * rzwB_i + rzwD_r * rzwB_i;


                     // CD

                    rhowCrhowD_r[i] = 2.0 * (rhowC_r[i] * rhowD_r[i] - rhowC_i[i] * rhowD_i[i]);

                    rhowCrhowD_i[i] = 2.0 * (rhowC_r[i] * rhowD_i[i] + rhowC_i[i] * rhowD_r[i]);

                    rxwCrhowD_r[i] = 2.0 * (rxwC_r * rhowD_r[i] - rxwC_i * rhowD_i[i]

                                            + rxwD_r * rhowC_r[i] - rxwD_i * rhowC_i[i]);

                    rxwCrhowD_i[i] = 2.0 * (rxwC_r * rhowD_i[i] + rxwC_i * rhowD_r[i]

                                            + rxwD_r * rhowC_i[i] + rxwD_i * rhowC_r[i]);

                    rywCrhowD_r[i] = 2.0 * (rywC_r * rhowD_r[i] - rywC_i * rhowD_i[i]

                                            + rywD_r * rhowC_r[i] - rywD_i * rhowC_i[i]);

                    rywCrhowD_i[i] = 2.0 * (rywC_r * rhowD_i[i] + rywC_i * rhowD_r[i]

                                            + rywD_r * rhowC_i[i] + rywD_i * rhowC_r[i]);

                    rzwCrhowD_r[i] = 2.0 * (rzwC_r * rhowD_r[i] - rzwC_i * rhowD_i[i]

                                            + rzwD_r * rhowC_r[i] - rzwD_i * rhowC_i[i]);

                    rzwCrhowD_i[i] = 2.0 * (rzwC_r * rhowD_i[i] + rzwC_i * rhowD_r[i]

                                            + rzwD_r * rhowC_i[i] + rzwD_i * rhowC_r[i]);

                    rxwCrxwD_r[i] = rxwC_r * rxwD_r - rxwC_i * rxwD_i

                                    + rxwD_r * rxwC_r - rxwD_i * rxwC_i;

                    rxwCrxwD_i[i] = rxwC_r * rxwD_i + rxwC_r * rxwD_i

                                    + rxwD_r * rxwC_i + rxwD_r * rxwC_i;

                    rxwCrywD_r[i] = rxwC_r * rywD_r - rxwC_i * rywD_i

                                    + rxwD_r * rywC_r - rxwD_i * rywC_i;

                    rxwCrywD_i[i] = rxwC_r * rywD_i + rxwC_r * rywD_i

                                    + rxwD_r * rywC_i + rxwD_r * rywC_i;

                    rxwCrzwD_r[i] = rxwC_r * rzwD_r - rxwC_i * rzwD_i

                                    + rxwD_r * rzwC_r - rxwD_i * rzwC_i;

                    rxwCrzwD_i[i] = rxwC_r * rzwD_i + rxwC_r * rzwD_i

                                    + rxwD_r * rzwC_i + rxwD_r * rzwC_i;

                    rywCrxwD_r[i] = rywC_r * rxwD_r - rywC_i * rxwD_i

                                    + rywD_r * rxwC_r - rywD_i * rxwC_i;

                    rywCrxwD_i[i] = rywC_r * rxwD_i + rywC_r * rxwD_i

                                    + rywD_r * rxwC_i + rywD_r * rxwC_i;

                    rywCrywD_r[i] = rywC_r * rywD_r - rywC_i * rywD_i

                                    + rywD_r * rywC_r - rywD_i * rywC_i;

                    rywCrywD_i[i] = rywC_r * rywD_i + rywC_r * rywD_i

                                    + rywD_r * rywC_i + rywD_r * rywC_i;

                    rywCrzwD_r[i] = rywC_r * rzwD_r - rywC_i * rzwD_i

                                    + rywD_r * rzwC_r - rywD_i * rzwC_i;

                    rywCrzwD_i[i] = rywC_r * rzwD_i + rywC_r * rzwD_i

                                    + rywD_r * rzwC_i + rywD_r * rzwC_i;

                    rzwCrxwD_r[i] = rzwC_r * rxwD_r - rzwC_i * rxwD_i

                                    + rzwD_r * rxwC_r - rzwD_i * rxwC_i;

                    rzwCrxwD_i[i] = rzwC_r * rxwD_i + rzwC_r * rxwD_i

                                    + rzwD_r * rxwC_i + rzwD_r * rxwC_i;

                    rzwCrywD_r[i] = rzwC_r * rywD_r - rzwC_i * rywD_i

                                    + rzwD_r * rywC_r - rzwD_i * rywC_i;

                    rzwCrywD_i[i] = rzwC_r * rywD_i + rzwC_r * rywD_i

                                    + rzwD_r * rywC_i + rzwD_r * rywC_i;

                    rzwCrzwD_r[i] = rzwC_r * rzwD_r - rzwC_i * rzwD_i

                                    + rzwD_r * rzwC_r - rzwD_i * rzwC_i;

                    rzwCrzwD_i[i] = rzwC_r * rzwD_i + rzwC_r * rzwD_i

                                    + rzwD_r * rzwC_i + rzwD_r * rzwC_i;


                }
            }

        }
        if (fstr::upcase(quadMode) == "CRF_II")
        {
            for (int32_t j = 0; j < numdens / 2; j++)
            {
                // Perturbed densities

                // Perturbed density products to be stored

                auto rhow1rhow2_r = rhow1rhow2(2 * j);

                auto rhow1rhow2_i = rhow1rhow2(2 * j + 1);

                // rxw1rhow2

                auto rxw1rhow2_r = rxw1rhow2(2 * j);

                auto rxw1rhow2_i = rxw1rhow2(2 * j + 1);

                // ryw1rhow2

                auto ryw1rhow2_r = ryw1rhow2(2 * j);

                auto ryw1rhow2_i = ryw1rhow2(2 * j + 1);

                // rzw1rhow2

                auto rzw1rhow2_r = rzw1rhow2(2 * j);

                auto rzw1rhow2_i = rzw1rhow2(2 * j + 1);

                // rAw1rBw2

                auto rxw1rxw2_r = rxw1rxw2(2 * j);

                auto rxw1rxw2_i = rxw1rxw2(2 * j + 1);

                auto rxw1ryw2_r = rxw1ryw2(2 * j);

                auto rxw1ryw2_i = rxw1ryw2(2 * j + 1);

                auto rxw1rzw2_r = rxw1rzw2(2 * j);

                auto rxw1rzw2_i = rxw1rzw2(2 * j + 1);

                auto ryw1rxw2_r = ryw1rxw2(2 * j);

                auto ryw1rxw2_i = ryw1rxw2(2 * j + 1);

                auto ryw1ryw2_r = ryw1ryw2(2 * j);

                auto ryw1ryw2_i = ryw1ryw2(2 * j + 1);

                auto ryw1rzw2_r = ryw1rzw2(2 * j);

                auto ryw1rzw2_i = ryw1rzw2(2 * j + 1);

                auto rzw1rxw2_r = rzw1rxw2(2 * j);

                auto rzw1rxw2_i = rzw1rxw2(2 * j + 1);

                auto rzw1ryw2_r = rzw1ryw2(2 * j);

                auto rzw1ryw2_i = rzw1ryw2(2 * j + 1);

                auto rzw1rzw2_r = rzw1rzw2(2 * j);

                auto rzw1rzw2_i = rzw1rzw2(2 * j + 1);

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

                auto gradBCx_r = rwDensityGrid.alphaDensityGradientX(12 * j + 6);

                auto gradBCy_r = rwDensityGrid.alphaDensityGradientY(12 * j + 6);

                auto gradBCz_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 6);


                auto rhoBC_i = rwDensityGrid.alphaDensity(12 * j + 7);

                auto gradBCx_i = rwDensityGrid.alphaDensityGradientX(12 * j + 7);

                auto gradBCy_i = rwDensityGrid.alphaDensityGradientY(12 * j + 7);

                auto gradBCz_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 7);


                auto rhoBD_r = rwDensityGrid.alphaDensity(12 * j + 8);

                auto gradBDx_r = rwDensityGrid.alphaDensityGradientX(12 * j + 8);

                auto gradBDy_r = rwDensityGrid.alphaDensityGradientY(12 * j + 8);

                auto gradBDz_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 8);


                auto rhoBD_i = rwDensityGrid.alphaDensity(12 * j + 9);

                auto gradBDx_i = rwDensityGrid.alphaDensityGradientX(12 * j + 9);

                auto gradBDy_i = rwDensityGrid.alphaDensityGradientY(12 * j + 9);

                auto gradBDz_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 9);

            
                auto rhoCD_r = rwDensityGrid.alphaDensity(12 * j + 10);

                auto gradCDx_r = rwDensityGrid.alphaDensityGradientX(12 * j + 10);

                auto gradCDy_r = rwDensityGrid.alphaDensityGradientY(12 * j + 10);

                auto gradCDz_r = rwDensityGrid.alphaDensityGradientZ(12 * j + 10);


                auto rhoCD_i = rwDensityGrid.alphaDensity(12 * j + 11);

                auto gradCDx_i = rwDensityGrid.alphaDensityGradientX(12 * j + 11);

                auto gradCDy_i = rwDensityGrid.alphaDensityGradientY(12 * j + 11);

                auto gradCDz_i = rwDensityGrid.alphaDensityGradientZ(12 * j + 11);


                for (int32_t i = 0; i < npoints; i++)
                {
                    // B CD

                    double rhow1_r = rhoB_r[i];

                    double rxw1_r = gradB_x_r[i];

                    double ryw1_r = gradB_y_r[i];

                    double rzw1_r = gradB_z_r[i];

                    double rhow1_i = rhoB_i[i];

                    double rxw1_i = gradB_x_i[i];

                    double ryw1_i = gradB_y_i[i];

                    double rzw1_i = gradB_z_i[i];

                    // RW2 densities

                    double rhow2_r = rhoCD_r[i];

                    double rxw2_r = gradCDx_r[i];

                    double ryw2_r = gradCDy_r[i];

                    double rzw2_r = gradCDz_r[i];

                    double rhow2_i = rhoCD_i[i];

                    double rxw2_i = gradCDx_i[i];

                    double ryw2_i = gradCDy_i[i];

                    double rzw2_i = gradCDz_i[i];

                    // Densities for terms 1-3
                    

                    rhow1rhow2_r[i] = 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    rhow1rhow2_i[i] = 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);

                    rxw1rhow2_r[i] = 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);

                    rxw1rhow2_i[i] = 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);

                    ryw1rhow2_r[i] = 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);

                    ryw1rhow2_i[i] = 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);

                    rzw1rhow2_r[i] = 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);

                    rzw1rhow2_i[i] = 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);

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

                    // C BD
                    
                    rhow1_r = rhoC_r[i];

                    rxw1_r = gradC_x_r[i];

                    ryw1_r = gradC_y_r[i];

                    rzw1_r = gradC_z_r[i];

                    rhow1_i = rhoC_i[i];

                    rxw1_i = gradC_x_i[i];

                    ryw1_i = gradC_y_i[i];

                    rzw1_i = gradC_z_i[i];


                    rhow2_r = rhoBD_r[i];

                    rxw2_r = gradBDx_r[i];

                    ryw2_r = gradBDy_r[i];

                    rzw2_r = gradBDz_r[i];

                    rhow2_i = rhoBD_i[i];

                    rxw2_i = gradBDx_i[i];

                    ryw2_i = gradBDy_i[i];

                    rzw2_i = gradBDz_i[i];

                    // Densities for terms 1-3
                    

                    rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);

                    rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);

                    rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);

                    ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);

                    ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);

                    rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);

                    rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);

                    rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;

                    rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;

                    rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;

                    rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;

                    rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;

                    // D BC
                    
                    rhow1_r = rhoD_r[i];

                    rxw1_r = gradD_x_r[i];

                    ryw1_r = gradD_y_r[i];

                    rzw1_r = gradD_z_r[i];

                    rhow1_i = rhoD_i[i];

                    rxw1_i = gradD_x_i[i];

                    ryw1_i = gradD_y_i[i];

                    rzw1_i = gradD_z_i[i];


                    rhow2_r = rhoBC_r[i];

                    rxw2_r = gradBCx_r[i];

                    ryw2_r = gradBCy_r[i];

                    rzw2_r = gradBCz_r[i];

                    rhow2_i = rhoBC_i[i];

                    rxw2_i = gradBCx_i[i];

                    ryw2_i = gradBCy_i[i];

                    rzw2_i = gradBCz_i[i];

                    // Densities for terms 1-3
                    

                    rhow1rhow2_r[i] += 2.0 * (rhow1_r * rhow2_r - rhow1_i * rhow2_i);

                    rhow1rhow2_i[i] += 2.0 * (rhow1_r * rhow2_i + rhow1_i * rhow2_r);

                    rxw1rhow2_r[i] += 2.0 * (rxw1_r * rhow2_r - rxw1_i * rhow2_i

                                            + rxw2_r * rhow1_r - rxw2_i * rhow1_i);

                    rxw1rhow2_i[i] += 2.0 * (rxw1_r * rhow2_i + rxw1_i * rhow2_r

                                            + rxw2_r * rhow1_i + rxw2_i * rhow1_r);

                    ryw1rhow2_r[i] += 2.0 * (ryw1_r * rhow2_r - ryw1_i * rhow2_i

                                            + ryw2_r * rhow1_r - ryw2_i * rhow1_i);

                    ryw1rhow2_i[i] += 2.0 * (ryw1_r * rhow2_i + ryw1_i * rhow2_r

                                            + ryw2_r * rhow1_i + ryw2_i * rhow1_r);

                    rzw1rhow2_r[i] += 2.0 * (rzw1_r * rhow2_r - rzw1_i * rhow2_i

                                            + rzw2_r * rhow1_r - rzw2_i * rhow1_i);

                    rzw1rhow2_i[i] += 2.0 * (rzw1_r * rhow2_i + rzw1_i * rhow2_r

                                            + rzw2_r * rhow1_i + rzw2_i * rhow1_r);

                    rxw1rxw2_r[i] += rxw1_r * rxw2_r - rxw1_i * rxw2_i

                                    + rxw2_r * rxw1_r - rxw2_i * rxw1_i;

                    rxw1rxw2_i[i] += rxw1_r * rxw2_i + rxw1_r * rxw2_i

                                    + rxw2_r * rxw1_i + rxw2_r * rxw1_i;

                    rxw1ryw2_r[i] += rxw1_r * ryw2_r - rxw1_i * ryw2_i

                                    + rxw2_r * ryw1_r - rxw2_i * ryw1_i;

                    rxw1ryw2_i[i] += rxw1_r * ryw2_i + rxw1_r * ryw2_i

                                    + rxw2_r * ryw1_i + rxw2_r * ryw1_i;

                    rxw1rzw2_r[i] += rxw1_r * rzw2_r - rxw1_i * rzw2_i

                                    + rxw2_r * rzw1_r - rxw2_i * rzw1_i;

                    rxw1rzw2_i[i] += rxw1_r * rzw2_i + rxw1_r * rzw2_i

                                    + rxw2_r * rzw1_i + rxw2_r * rzw1_i;

                    ryw1rxw2_r[i] += ryw1_r * rxw2_r - ryw1_i * rxw2_i

                                    + ryw2_r * rxw1_r - ryw2_i * rxw1_i;

                    ryw1rxw2_i[i] += ryw1_r * rxw2_i + ryw1_r * rxw2_i

                                    + ryw2_r * rxw1_i + ryw2_r * rxw1_i;

                    ryw1ryw2_r[i] += ryw1_r * ryw2_r - ryw1_i * ryw2_i

                                    + ryw2_r * ryw1_r - ryw2_i * ryw1_i;

                    ryw1ryw2_i[i] += ryw1_r * ryw2_i + ryw1_r * ryw2_i

                                    + ryw2_r * ryw1_i + ryw2_r * ryw1_i;

                    ryw1rzw2_r[i] += ryw1_r * rzw2_r - ryw1_i * rzw2_i

                                    + ryw2_r * rzw1_r - ryw2_i * rzw1_i;

                    ryw1rzw2_i[i] += ryw1_r * rzw2_i + ryw1_r * rzw2_i

                                    + ryw2_r * rzw1_i + ryw2_r * rzw1_i;

                    rzw1rxw2_r[i] += rzw1_r * rxw2_r - rzw1_i * rxw2_i

                                    + rzw2_r * rxw1_r - rzw2_i * rxw1_i;

                    rzw1rxw2_i[i] += rzw1_r * rxw2_i + rzw1_r * rxw2_i

                                    + rzw2_r * rxw1_i + rzw2_r * rxw1_i;

                    rzw1ryw2_r[i] += rzw1_r * ryw2_r - rzw1_i * ryw2_i

                                    + rzw2_r * ryw1_r - rzw2_i * ryw1_i;

                    rzw1ryw2_i[i] += rzw1_r * ryw2_i + rzw1_r * ryw2_i

                                    + rzw2_r * ryw1_i + rzw2_r * ryw1_i;

                    rzw1rzw2_r[i] += rzw1_r * rzw2_r - rzw1_i * rzw2_i

                                    + rzw2_r * rzw1_r - rzw2_i * rzw1_i;

                    rzw1rzw2_i[i] += rzw1_r * rzw2_i + rzw1_r * rzw2_i

                                    + rzw2_r * rzw1_i + rzw2_r * rzw1_i;

                }
            }
        }
        if (fstr::upcase(quadMode) == "SHG")
        {
            for (int32_t j = 0; j < numdens / 12; j++)

            {
                // Density products to be stored

                // RhoRho part

                auto rho_sig_x_r = rhow1rhow2(12 * j);

                auto rho_sig_x_i = rhow1rhow2(12 * j + 1);

                auto rho_sig_y_r = rhow1rhow2(12 * j + 2);

                auto rho_sig_y_i = rhow1rhow2(12 * j + 3);

                auto rho_sig_z_r = rhow1rhow2(12 * j + 4);

                auto rho_sig_z_i = rhow1rhow2(12 * j + 5);

                auto rho_lam_xy_r = rhow1rhow2(12 * j + 6);

                auto rho_lam_xy_i = rhow1rhow2(12 * j + 7);

                auto rho_lam_xz_r = rhow1rhow2(12 * j + 8);

                auto rho_lam_xz_i = rhow1rhow2(12 * j + 9);

                auto rho_lam_yz_r = rhow1rhow2(12 * j + 10);

                auto rho_lam_yz_i = rhow1rhow2(12 * j + 11);

                // rxw1rhow part

                auto rxw1rhow2_sig_x_r = rxw1rhow2(12 * j);

                auto rxw1rhow2_sig_x_i = rxw1rhow2(12 * j + 1);

                auto rxw1rhow2_sig_y_r = rxw1rhow2(12 * j + 2);

                auto rxw1rhow2_sig_y_i = rxw1rhow2(12 * j + 3);

                auto rxw1rhow2_sig_z_r = rxw1rhow2(12 * j + 4);

                auto rxw1rhow2_sig_z_i = rxw1rhow2(12 * j + 5);

                auto rxw1rhow2_lam_xy_r = rxw1rhow2(12 * j + 6);

                auto rxw1rhow2_lam_xy_i = rxw1rhow2(12 * j + 7);

                auto rxw1rhow2_lam_xz_r = rxw1rhow2(12 * j + 8);

                auto rxw1rhow2_lam_xz_i = rxw1rhow2(12 * j + 9);

                auto rxw1rhow2_lam_yz_r = rxw1rhow2(12 * j + 10);

                auto rxw1rhow2_lam_yz_i = rxw1rhow2(12 * j + 11);

                // ryw1rhow part

                auto ryw1rhow2_sig_x_r = ryw1rhow2(12 * j);

                auto ryw1rhow2_sig_x_i = ryw1rhow2(12 * j + 1);

                auto ryw1rhow2_sig_y_r = ryw1rhow2(12 * j + 2);

                auto ryw1rhow2_sig_y_i = ryw1rhow2(12 * j + 3);

                auto ryw1rhow2_sig_z_r = ryw1rhow2(12 * j + 4);

                auto ryw1rhow2_sig_z_i = ryw1rhow2(12 * j + 5);

                auto ryw1rhow2_lam_xy_r = ryw1rhow2(12 * j + 6);

                auto ryw1rhow2_lam_xy_i = ryw1rhow2(12 * j + 7);

                auto ryw1rhow2_lam_xz_r = ryw1rhow2(12 * j + 8);

                auto ryw1rhow2_lam_xz_i = ryw1rhow2(12 * j + 9);

                auto ryw1rhow2_lam_yz_r = ryw1rhow2(12 * j + 10);

                auto ryw1rhow2_lam_yz_i = ryw1rhow2(12 * j + 11);

                // rzw1rhow part

                auto rzw1rhow2_sig_x_r = rzw1rhow2(12 * j);

                auto rzw1rhow2_sig_x_i = rzw1rhow2(12 * j + 1);

                auto rzw1rhow2_sig_y_r = rzw1rhow2(12 * j + 2);

                auto rzw1rhow2_sig_y_i = rzw1rhow2(12 * j + 3);

                auto rzw1rhow2_sig_z_r = rzw1rhow2(12 * j + 4);

                auto rzw1rhow2_sig_z_i = rzw1rhow2(12 * j + 5);

                auto rzw1rhow2_lam_xy_r = rzw1rhow2(12 * j + 6);

                auto rzw1rhow2_lam_xy_i = rzw1rhow2(12 * j + 7);

                auto rzw1rhow2_lam_xz_r = rzw1rhow2(12 * j + 8);

                auto rzw1rhow2_lam_xz_i = rzw1rhow2(12 * j + 9);

                auto rzw1rhow2_lam_yz_r = rzw1rhow2(12 * j + 10);

                auto rzw1rhow2_lam_yz_i = rzw1rhow2(12 * j + 11);

                // rxw1rxw2 part

                auto rxw1rxw2_sig_x_r = rxw1rxw2(12 * j);

                auto rxw1rxw2_sig_x_i = rxw1rxw2(12 * j + 1);

                auto rxw1rxw2_sig_y_r = rxw1rxw2(12 * j + 2);

                auto rxw1rxw2_sig_y_i = rxw1rxw2(12 * j + 3);

                auto rxw1rxw2_sig_z_r = rxw1rxw2(12 * j + 4);

                auto rxw1rxw2_sig_z_i = rxw1rxw2(12 * j + 5);

                auto rxw1rxw2_lam_xy_r = rxw1rxw2(12 * j + 6);

                auto rxw1rxw2_lam_xy_i = rxw1rxw2(12 * j + 7);

                auto rxw1rxw2_lam_xz_r = rxw1rxw2(12 * j + 8);

                auto rxw1rxw2_lam_xz_i = rxw1rxw2(12 * j + 9);

                auto rxw1rxw2_lam_yz_r = rxw1rxw2(12 * j + 10);

                auto rxw1rxw2_lam_yz_i = rxw1rxw2(12 * j + 11);

                // rxw1ryw2 part

                auto rxw1ryw2_sig_x_r = rxw1ryw2(12 * j);

                auto rxw1ryw2_sig_x_i = rxw1ryw2(12 * j + 1);

                auto rxw1ryw2_sig_y_r = rxw1ryw2(12 * j + 2);

                auto rxw1ryw2_sig_y_i = rxw1ryw2(12 * j + 3);

                auto rxw1ryw2_sig_z_r = rxw1ryw2(12 * j + 4);

                auto rxw1ryw2_sig_z_i = rxw1ryw2(12 * j + 5);

                auto rxw1ryw2_lam_xy_r = rxw1ryw2(12 * j + 6);

                auto rxw1ryw2_lam_xy_i = rxw1ryw2(12 * j + 7);

                auto rxw1ryw2_lam_xz_r = rxw1ryw2(12 * j + 8);

                auto rxw1ryw2_lam_xz_i = rxw1ryw2(12 * j + 9);

                auto rxw1ryw2_lam_yz_r = rxw1ryw2(12 * j + 10);

                auto rxw1ryw2_lam_yz_i = rxw1ryw2(12 * j + 11);

                // rxw1rzw2 part

                auto rxw1rzw2_sig_x_r = rxw1rzw2(12 * j);

                auto rxw1rzw2_sig_x_i = rxw1rzw2(12 * j + 1);

                auto rxw1rzw2_sig_y_r = rxw1rzw2(12 * j + 2);

                auto rxw1rzw2_sig_y_i = rxw1rzw2(12 * j + 3);

                auto rxw1rzw2_sig_z_r = rxw1rzw2(12 * j + 4);

                auto rxw1rzw2_sig_z_i = rxw1rzw2(12 * j + 5);

                auto rxw1rzw2_lam_xy_r = rxw1rzw2(12 * j + 6);

                auto rxw1rzw2_lam_xy_i = rxw1rzw2(12 * j + 7);

                auto rxw1rzw2_lam_xz_r = rxw1rzw2(12 * j + 8);

                auto rxw1rzw2_lam_xz_i = rxw1rzw2(12 * j + 9);

                auto rxw1rzw2_lam_yz_r = rxw1rzw2(12 * j + 10);

                auto rxw1rzw2_lam_yz_i = rxw1rzw2(12 * j + 11);

                // ryw1rxw2 part

                auto ryw1rxw2_sig_x_r = ryw1rxw2(12 * j);

                auto ryw1rxw2_sig_x_i = ryw1rxw2(12 * j + 1);

                auto ryw1rxw2_sig_y_r = ryw1rxw2(12 * j + 2);

                auto ryw1rxw2_sig_y_i = ryw1rxw2(12 * j + 3);

                auto ryw1rxw2_sig_z_r = ryw1rxw2(12 * j + 4);

                auto ryw1rxw2_sig_z_i = ryw1rxw2(12 * j + 5);

                auto ryw1rxw2_lam_xy_r = ryw1rxw2(12 * j + 6);

                auto ryw1rxw2_lam_xy_i = ryw1rxw2(12 * j + 7);

                auto ryw1rxw2_lam_xz_r = ryw1rxw2(12 * j + 8);

                auto ryw1rxw2_lam_xz_i = ryw1rxw2(12 * j + 9);

                auto ryw1rxw2_lam_yz_r = ryw1rxw2(12 * j + 10);

                auto ryw1rxw2_lam_yz_i = ryw1rxw2(12 * j + 11);

                // ryw1ryw2 part

                auto ryw1ryw2_sig_x_r = ryw1ryw2(12 * j);

                auto ryw1ryw2_sig_x_i = ryw1ryw2(12 * j + 1);

                auto ryw1ryw2_sig_y_r = ryw1ryw2(12 * j + 2);

                auto ryw1ryw2_sig_y_i = ryw1ryw2(12 * j + 3);

                auto ryw1ryw2_sig_z_r = ryw1ryw2(12 * j + 4);

                auto ryw1ryw2_sig_z_i = ryw1ryw2(12 * j + 5);

                auto ryw1ryw2_lam_xy_r = ryw1ryw2(12 * j + 6);

                auto ryw1ryw2_lam_xy_i = ryw1ryw2(12 * j + 7);

                auto ryw1ryw2_lam_xz_r = ryw1ryw2(12 * j + 8);

                auto ryw1ryw2_lam_xz_i = ryw1ryw2(12 * j + 9);

                auto ryw1ryw2_lam_yz_r = ryw1ryw2(12 * j + 10);

                auto ryw1ryw2_lam_yz_i = ryw1ryw2(12 * j + 11);

                // ryw1rzw2 part

                auto ryw1rzw2_sig_x_r = ryw1rzw2(12 * j);

                auto ryw1rzw2_sig_x_i = ryw1rzw2(12 * j + 1);

                auto ryw1rzw2_sig_y_r = ryw1rzw2(12 * j + 2);

                auto ryw1rzw2_sig_y_i = ryw1rzw2(12 * j + 3);

                auto ryw1rzw2_sig_z_r = ryw1rzw2(12 * j + 4);

                auto ryw1rzw2_sig_z_i = ryw1rzw2(12 * j + 5);

                auto ryw1rzw2_lam_xy_r = ryw1rzw2(12 * j + 6);

                auto ryw1rzw2_lam_xy_i = ryw1rzw2(12 * j + 7);

                auto ryw1rzw2_lam_xz_r = ryw1rzw2(12 * j + 8);

                auto ryw1rzw2_lam_xz_i = ryw1rzw2(12 * j + 9);

                auto ryw1rzw2_lam_yz_r = ryw1rzw2(12 * j + 10);

                auto ryw1rzw2_lam_yz_i = ryw1rzw2(12 * j + 11);

                // rzw1rxw2 part

                auto rzw1rxw2_sig_x_r = rzw1rxw2(12 * j);

                auto rzw1rxw2_sig_x_i = rzw1rxw2(12 * j + 1);

                auto rzw1rxw2_sig_y_r = rzw1rxw2(12 * j + 2);

                auto rzw1rxw2_sig_y_i = rzw1rxw2(12 * j + 3);

                auto rzw1rxw2_sig_z_r = rzw1rxw2(12 * j + 4);

                auto rzw1rxw2_sig_z_i = rzw1rxw2(12 * j + 5);

                auto rzw1rxw2_lam_xy_r = rzw1rxw2(12 * j + 6);

                auto rzw1rxw2_lam_xy_i = rzw1rxw2(12 * j + 7);

                auto rzw1rxw2_lam_xz_r = rzw1rxw2(12 * j + 8);

                auto rzw1rxw2_lam_xz_i = rzw1rxw2(12 * j + 9);

                auto rzw1rxw2_lam_yz_r = rzw1rxw2(12 * j + 10);

                auto rzw1rxw2_lam_yz_i = rzw1rxw2(12 * j + 11);

                // rzw1ryw2 part

                auto rzw1ryw2_sig_x_r = rzw1ryw2(12 * j);

                auto rzw1ryw2_sig_x_i = rzw1ryw2(12 * j + 1);

                auto rzw1ryw2_sig_y_r = rzw1ryw2(12 * j + 2);

                auto rzw1ryw2_sig_y_i = rzw1ryw2(12 * j + 3);

                auto rzw1ryw2_sig_z_r = rzw1ryw2(12 * j + 4);

                auto rzw1ryw2_sig_z_i = rzw1ryw2(12 * j + 5);

                auto rzw1ryw2_lam_xy_r = rzw1ryw2(12 * j + 6);

                auto rzw1ryw2_lam_xy_i = rzw1ryw2(12 * j + 7);

                auto rzw1ryw2_lam_xz_r = rzw1ryw2(12 * j + 8);

                auto rzw1ryw2_lam_xz_i = rzw1ryw2(12 * j + 9);

                auto rzw1ryw2_lam_yz_r = rzw1ryw2(12 * j + 10);

                auto rzw1ryw2_lam_yz_i = rzw1ryw2(12 * j + 11);

                // rzw1rzw2 part

                auto rzw1rzw2_sig_x_r = rzw1rzw2(12 * j);

                auto rzw1rzw2_sig_x_i = rzw1rzw2(12 * j + 1);

                auto rzw1rzw2_sig_y_r = rzw1rzw2(12 * j + 2);

                auto rzw1rzw2_sig_y_i = rzw1rzw2(12 * j + 3);

                auto rzw1rzw2_sig_z_r = rzw1rzw2(12 * j + 4);

                auto rzw1rzw2_sig_z_i = rzw1rzw2(12 * j + 5);

                auto rzw1rzw2_lam_xy_r = rzw1rzw2(12 * j + 6);

                auto rzw1rzw2_lam_xy_i = rzw1rzw2(12 * j + 7);

                auto rzw1rzw2_lam_xz_r = rzw1rzw2(12 * j + 8);

                auto rzw1rzw2_lam_xz_i = rzw1rzw2(12 * j + 9);

                auto rzw1rzw2_lam_yz_r = rzw1rzw2(12 * j + 10);

                auto rzw1rzw2_lam_yz_i = rzw1rzw2(12 * j + 11);

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

                auto rho_sig_x_r = rhow1rhow2(6 * j);

                auto rho_sig_y_r = rhow1rhow2(6 * j + 1);

                auto rho_sig_z_r = rhow1rhow2(6 * j + 2);

                auto rho_lam_xy_r = rhow1rhow2(6 * j + 3);

                auto rho_lam_xz_r = rhow1rhow2(6 * j + 4);

                auto rho_lam_yz_r = rhow1rhow2(6 * j + 5);

                // rxw1rhow part

                auto rxw1rhow2_sig_x_r = rxw1rhow2(6 * j);

                auto rxw1rhow2_sig_y_r = rxw1rhow2(6 * j + 1);

                auto rxw1rhow2_sig_z_r = rxw1rhow2(6 * j + 2);

                auto rxw1rhow2_lam_xy_r = rxw1rhow2(6 * j + 3);

                auto rxw1rhow2_lam_xz_r = rxw1rhow2(6 * j + 4);

                auto rxw1rhow2_lam_yz_r = rxw1rhow2(6 * j + 5);

                // ryw1rhow part

                auto ryw1rhow2_sig_x_r = ryw1rhow2(6 * j);

                auto ryw1rhow2_sig_y_r = ryw1rhow2(6 * j + 1);

                auto ryw1rhow2_sig_z_r = ryw1rhow2(6 * j + 2);

                auto ryw1rhow2_lam_xy_r = ryw1rhow2(6 * j + 3);

                auto ryw1rhow2_lam_xz_r = ryw1rhow2(6 * j + 4);

                auto ryw1rhow2_lam_yz_r = ryw1rhow2(6 * j + 5);

                // rzw1rhow part

                auto rzw1rhow2_sig_x_r = rzw1rhow2(6 * j);

                auto rzw1rhow2_sig_y_r = rzw1rhow2(6 * j + 1);

                auto rzw1rhow2_sig_z_r = rzw1rhow2(6 * j + 2);

                auto rzw1rhow2_lam_xy_r = rzw1rhow2(6 * j + 3);

                auto rzw1rhow2_lam_xz_r = rzw1rhow2(6 * j + 4);

                auto rzw1rhow2_lam_yz_r = rzw1rhow2(6 * j + 5);

                // rxw1rxw2 part

                auto rxw1rxw2_sig_x_r = rxw1rxw2(6 * j);

                auto rxw1rxw2_sig_y_r = rxw1rxw2(6 * j + 1);

                auto rxw1rxw2_sig_z_r = rxw1rxw2(6 * j + 2);

                auto rxw1rxw2_lam_xy_r = rxw1rxw2(6 * j + 3);

                auto rxw1rxw2_lam_xz_r = rxw1rxw2(6 * j + 4);

                auto rxw1rxw2_lam_yz_r = rxw1rxw2(6 * j + 5);

                // rxw1ryw2 part

                auto rxw1ryw2_sig_x_r = rxw1ryw2(6 * j);

                auto rxw1ryw2_sig_y_r = rxw1ryw2(6 * j + 1);

                auto rxw1ryw2_sig_z_r = rxw1ryw2(6 * j + 2);

                auto rxw1ryw2_lam_xy_r = rxw1ryw2(6 * j + 3);

                auto rxw1ryw2_lam_xz_r = rxw1ryw2(6 * j + 4);

                auto rxw1ryw2_lam_yz_r = rxw1ryw2(6 * j + 5);

                // rxw1rzw2 part

                auto rxw1rzw2_sig_x_r = rxw1rzw2(6 * j);

                auto rxw1rzw2_sig_y_r = rxw1rzw2(6 * j + 1);

                auto rxw1rzw2_sig_z_r = rxw1rzw2(6 * j + 2);

                auto rxw1rzw2_lam_xy_r = rxw1rzw2(6 * j + 3);

                auto rxw1rzw2_lam_xz_r = rxw1rzw2(6 * j + 4);

                auto rxw1rzw2_lam_yz_r = rxw1rzw2(6 * j + 5);

                // ryw1rxw2 part

                auto ryw1rxw2_sig_x_r = ryw1rxw2(6 * j);

                auto ryw1rxw2_sig_y_r = ryw1rxw2(6 * j + 1);

                auto ryw1rxw2_sig_z_r = ryw1rxw2(6 * j + 2);

                auto ryw1rxw2_lam_xy_r = ryw1rxw2(6 * j + 3);

                auto ryw1rxw2_lam_xz_r = ryw1rxw2(6 * j + 4);

                auto ryw1rxw2_lam_yz_r = ryw1rxw2(6 * j + 5);

                // ryw1ryw2 part

                auto ryw1ryw2_sig_x_r = ryw1ryw2(6 * j);

                auto ryw1ryw2_sig_y_r = ryw1ryw2(6 * j + 1);

                auto ryw1ryw2_sig_z_r = ryw1ryw2(6 * j + 2);

                auto ryw1ryw2_lam_xy_r = ryw1ryw2(6 * j + 3);

                auto ryw1ryw2_lam_xz_r = ryw1ryw2(6 * j + 4);

                auto ryw1ryw2_lam_yz_r = ryw1ryw2(6 * j + 5);

                // ryw1rzw2 part

                auto ryw1rzw2_sig_x_r = ryw1rzw2(6 * j);

                auto ryw1rzw2_sig_y_r = ryw1rzw2(6 * j + 1);

                auto ryw1rzw2_sig_z_r = ryw1rzw2(6 * j + 2);

                auto ryw1rzw2_lam_xy_r = ryw1rzw2(6 * j + 3);

                auto ryw1rzw2_lam_xz_r = ryw1rzw2(6 * j + 4);

                auto ryw1rzw2_lam_yz_r = ryw1rzw2(6 * j + 5);

                // rzw1rxw2 part

                auto rzw1rxw2_sig_x_r = rzw1rxw2(6 * j);

                auto rzw1rxw2_sig_y_r = rzw1rxw2(6 * j + 1);

                auto rzw1rxw2_sig_z_r = rzw1rxw2(6 * j + 2);

                auto rzw1rxw2_lam_xy_r = rzw1rxw2(6 * j + 3);

                auto rzw1rxw2_lam_xz_r = rzw1rxw2(6 * j + 4);

                auto rzw1rxw2_lam_yz_r = rzw1rxw2(6 * j + 5);

                // rzw1ryw2 part

                auto rzw1ryw2_sig_x_r = rzw1ryw2(6 * j);

                auto rzw1ryw2_sig_y_r = rzw1ryw2(6 * j + 1);

                auto rzw1ryw2_sig_z_r = rzw1ryw2(6 * j + 2);

                auto rzw1ryw2_lam_xy_r = rzw1ryw2(6 * j + 3);

                auto rzw1ryw2_lam_xz_r = rzw1ryw2(6 * j + 4);

                auto rzw1ryw2_lam_yz_r = rzw1ryw2(6 * j + 5);

                // rzw1rzw2 part

                auto rzw1rzw2_sig_x_r = rzw1rzw2(6 * j);

                auto rzw1rzw2_sig_y_r = rzw1rzw2(6 * j + 1);

                auto rzw1rzw2_sig_z_r = rzw1rzw2(6 * j + 2);

                auto rzw1rzw2_lam_xy_r = rzw1rzw2(6 * j + 3);

                auto rzw1rzw2_lam_xz_r = rzw1rzw2(6 * j + 4);

                auto rzw1rzw2_lam_yz_r = rzw1rzw2(6 * j + 5);

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

                auto rhow1rhow2_r = rhow1rhow2(2 * j);

                auto rhow1rhow2_i = rhow1rhow2(2 * j + 1);

                // rxw1rhow2

                auto rxw1rhow2_r = rxw1rhow2(2 * j);

                auto rxw1rhow2_i = rxw1rhow2(2 * j + 1);

                // ryw1rhow2

                auto ryw1rhow2_r = ryw1rhow2(2 * j);

                auto ryw1rhow2_i = ryw1rhow2(2 * j + 1);

                // rzw1rhow2

                auto rzw1rhow2_r = rzw1rhow2(2 * j);

                auto rzw1rhow2_i = rzw1rhow2(2 * j + 1);

                // rAw1rBw2

                auto rxw1rxw2_r = rxw1rxw2(2 * j);

                auto rxw1rxw2_i = rxw1rxw2(2 * j + 1);

                auto rxw1ryw2_r = rxw1ryw2(2 * j);

                auto rxw1ryw2_i = rxw1ryw2(2 * j + 1);

                auto rxw1rzw2_r = rxw1rzw2(2 * j);

                auto rxw1rzw2_i = rxw1rzw2(2 * j + 1);

                auto ryw1rxw2_r = ryw1rxw2(2 * j);

                auto ryw1rxw2_i = ryw1rxw2(2 * j + 1);

                auto ryw1ryw2_r = ryw1ryw2(2 * j);

                auto ryw1ryw2_i = ryw1ryw2(2 * j + 1);

                auto ryw1rzw2_r = ryw1rzw2(2 * j);

                auto ryw1rzw2_i = ryw1rzw2(2 * j + 1);

                auto rzw1rxw2_r = rzw1rxw2(2 * j);

                auto rzw1rxw2_i = rzw1rxw2(2 * j + 1);

                auto rzw1ryw2_r = rzw1ryw2(2 * j);

                auto rzw1ryw2_i = rzw1ryw2(2 * j + 1);

                auto rzw1rzw2_r = rzw1rzw2(2 * j);

                auto rzw1rzw2_i = rzw1rzw2(2 * j + 1);

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
