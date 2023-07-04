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

#include "SADGuessDriver.hpp"

#include <algorithm>

#include <mpi.h>

#include "AODensityMatrix.hpp"
#include "DenseDiagonalizer.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DensityMatrixType.hpp"
#include "ErrorHandler.hpp"
#include "MpiFunc.hpp"
#include "OverlapIntegralsDriver.hpp"
#include "PartialCharges.hpp"
#include "StringFormat.hpp"

CSADGuessDriver::CSADGuessDriver(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;
}

CSADGuessDriver::~CSADGuessDriver()
{
}

std::vector<double>
CSADGuessDriver::_getOcc1s(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s
    return std::vector<double>({occ});
}

std::vector<double>
CSADGuessDriver::_getOcc2s(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s
    return std::vector<double>({1.0, occ});
}

std::vector<double>
CSADGuessDriver::_getOcc2s2p(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s   2p-1 2p0  2p+1
    return std::vector<double>({1.0, occ, occ, occ, occ});
}

std::vector<double>
CSADGuessDriver::_getOcc3s(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s   3s   2p-1 2p0  2p+1
    return std::vector<double>({1.0, 1.0, occ, 1.0, 1.0, 1.0});
}

std::vector<double>
CSADGuessDriver::_getOcc3s3p(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s   3s   2p-1 3p-1 2p0  3p0  2p+1 3p+1
    return std::vector<double>({1.0, 1.0, occ, 1.0, occ, 1.0, occ, 1.0, occ});
}

std::vector<double>
CSADGuessDriver::_getOcc4s(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s   3s   4s   2p-1 3p-1 2p0  3p0  2p+1 3p+1
    return std::vector<double>({1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
}

std::vector<double>
CSADGuessDriver::_getOcc3d(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s   3s   4s  2p-1 3p-1 2p0  3p0  2p+1 3p+1 3d-2 3d-1 3d0  3d+1 3d+2
    return std::vector<double>({1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, occ, occ, occ, occ, occ});
}

std::vector<double>
CSADGuessDriver::_getOcc4s4p(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s   3s   4s  2p-1 3p-1 4p-1 2p0  3p0  4p0  2p+1 3p+1 4p+1 3d-2 3d-1 3d0  3d+1 3d+2
    return std::vector<double>({1.0, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0});
}

std::vector<double>
CSADGuessDriver::_getOcc5s(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s   3s   4s   5s  2p-1 3p-1 4p-1 2p0  3p0  4p0  2p+1 3p+1 4p+1 3d-2 3d-1 3d0  3d+1 3d+2
    return std::vector<double>({1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
}

std::vector<double>
CSADGuessDriver::_getOcc4d(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s   3s   4s   5s  2p-1 3p-1 4p-1 2p0  3p0  4p0  2p+1 3p+1 4p+1 3d-2 4d-2 3d-1 4d-1 3d0  4d0  3d+1 4d+1 3d+2 4d+2
    return std::vector<double>({1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, occ, 1.0, occ, 1.0, occ, 1.0, occ, 1.0, occ});
}

std::vector<double>
CSADGuessDriver::_getOcc5s5p(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s   3s   4s   5s  2p-1 3p-1 4p-1 5p-1 2p0  3p0  4p0  5p0  2p+1 3p+1 4p+1 5p+1 3d-2 4d-2 3d-1 4d-1 3d0  4d0  3d+1 4d+1 3d+2 4d+2
    return std::vector<double>({1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
}

std::vector<double>
CSADGuessDriver::_getOcc6s(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s   3s   4s   5s  6s   2p-1 3p-1 4p-1 5p-1 2p0  3p0  4p0  5p0  2p+1 3p+1 4p+1 5p+1 3d-2 4d-2 3d-1 4d-1 3d0  4d0  3d+1 4d+1 3d+2 4d+2
    return std::vector<double>({1.0, 1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
}

std::vector<double>
CSADGuessDriver::_getOcc4f(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s   3s   4s   5s  6s   2p-1 3p-1 4p-1 5p-1 2p0  3p0  4p0  5p0  2p+1 3p+1 4p+1 5p+1 3d-2 4d-2 3d-1 4d-1 3d0  4d0  3d+1 4d+1 3d+2 4d+2 4f-3 4f-2 4f-1 4f0  4f+1 4f+2 4f+3
    return std::vector<double>({1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, occ, occ, occ, occ, occ, occ, occ});
}

std::vector<double>
CSADGuessDriver::_getOcc5d(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s   3s   4s   5s  6s   2p-1 3p-1 4p-1 5p-1 2p0  3p0  4p0  5p0  2p+1 3p+1 4p+1 5p+1 3d-2 4d-2 5d-2 3d-1 4d-1 5d-1 3d0  4d0  5d0  3d+1 4d+1 5d+1 3d+2 4d+2 5d+2 4f-3 4f-2 4f-1 4f0  4f+1 4f+2 4f+3
    return std::vector<double>({1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
}

std::vector<double>
CSADGuessDriver::_getOcc6s6p(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s   3s   4s   5s  6s   2p-1 3p-1 4p-1 5p-1 6s-1 2p0  3p0  4p0  5p0  6s0  2p+1 3p+1 4p+1 5p+1 6s+1 3d-2 4d-2 5d-2 3d-1 4d-1 5d-1 3d0  4d0  5d0  3d+1 4d+1 5d+1 3d+2 4d+2 5d+2 4f-3 4f-2 4f-1 4f0  4f+1 4f+2 4f+3
    return std::vector<double>({1.0, 1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
}

std::vector<double>
CSADGuessDriver::_getOcc7s(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s   3s   4s   5s  6s   7s   2p-1 3p-1 4p-1 5p-1 6s-1 2p0  3p0  4p0  5p0  6s0  2p+1 3p+1 4p+1 5p+1 6s+1 3d-2 4d-2 5d-2 3d-1 4d-1 5d-1 3d0  4d0  5d0  3d+1 4d+1 5d+1 3d+2 4d+2 5d+2 4f-3 4f-2 4f-1 4f0  4f+1 4f+2 4f+3
    return std::vector<double>({1.0, 1.0, 1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
}

std::vector<double>
CSADGuessDriver::_getOcc5f(const double nocc) const
{
    const double occ = std::max(0.0, nocc);

    //                           1s   2s   3s   4s   5s  6s   7s   2p-1 3p-1 4p-1 5p-1 6s-1 2p0  3p0  4p0  5p0  6s0  2p+1 3p+1 4p+1 5p+1 6s+1 3d-2 4d-2 5d-2 3d-1 4d-1 5d-1 3d0  4d0  5d0  3d+1 4d+1 5d+1 3d+2 4d+2 5d+2 4f-3 5f-3 4f-2 5f-2 4f-1 5f-1 4f0  5f0  4f+1 5f+1 4f+2 5f+2 4f+3 5f+3
    return std::vector<double>({1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, occ, 1.0, occ, 1.0, occ, 1.0, occ, 1.0, occ, 1.0, occ, 1.0, occ});
}

std::vector<double>
CSADGuessDriver::getOccupationNumbersForElement(const int32_t elem_id, const double nelec) const
{
    switch (elem_id)
    {
            // dummy atom

        case 0:
            return std::vector<double>();

            // H,He

        case 1:
            return _getOcc1s(0.5 + nelec);

        case 2:
            return _getOcc1s(1.0 + nelec);

            // Li,Be

        case 3:
            return _getOcc2s(0.5 + nelec);

        case 4:
            return _getOcc2s(1.0 + nelec);

            // B,C,N,O,F,Ne

        case 5:
            return _getOcc2s2p(0.375 + nelec / 4.0);

        case 6:
            return _getOcc2s2p(0.500 + nelec / 4.0);

        case 7:
            return _getOcc2s2p(0.625 + nelec / 4.0);

        case 8:
            return _getOcc2s2p(0.750 + nelec / 4.0);

        case 9:
            return _getOcc2s2p(0.875 + nelec / 4.0);

        case 10:
            return _getOcc2s2p(1.000 + nelec / 4.0);

            // Na,Mg

        case 11:
            return _getOcc3s(0.5 + nelec);

        case 12:
            return _getOcc3s(1.0 + nelec);

            // Al,Si,P,S,Cl,Ar

        case 13:
            return _getOcc3s3p(0.375 + nelec / 4.0);

        case 14:
            return _getOcc3s3p(0.500 + nelec / 4.0);

        case 15:
            return _getOcc3s3p(0.625 + nelec / 4.0);

        case 16:
            return _getOcc3s3p(0.750 + nelec / 4.0);

        case 17:
            return _getOcc3s3p(0.875 + nelec / 4.0);

        case 18:
            return _getOcc3s3p(1.000 + nelec / 4.0);

            // K,Ca

        case 19:
            return _getOcc4s(0.5 + nelec);

        case 20:
            return _getOcc4s(1.0 + nelec);

            // Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn

        case 21:
            return _getOcc3d(0.1 + nelec / 5.0);

        case 22:
            return _getOcc3d(0.2 + nelec / 5.0);

        case 23:
            return _getOcc3d(0.3 + nelec / 5.0);

        case 24:
            return _getOcc3d(0.4 + nelec / 5.0);

        case 25:
            return _getOcc3d(0.5 + nelec / 5.0);

        case 26:
            return _getOcc3d(0.6 + nelec / 5.0);

        case 27:
            return _getOcc3d(0.7 + nelec / 5.0);

        case 28:
            return _getOcc3d(0.8 + nelec / 5.0);

        case 29:
            return _getOcc3d(0.9 + nelec / 5.0);

        case 30:
            return _getOcc3d(1.0 + nelec / 5.0);

            // Ga,Ge,As,Se,Br,Kr

        case 31:
            return _getOcc4s4p(0.375 + nelec / 4.0);

        case 32:
            return _getOcc4s4p(0.500 + nelec / 4.0);

        case 33:
            return _getOcc4s4p(0.625 + nelec / 4.0);

        case 34:
            return _getOcc4s4p(0.750 + nelec / 4.0);

        case 35:
            return _getOcc4s4p(0.875 + nelec / 4.0);

        case 36:
            return _getOcc4s4p(1.000 + nelec / 4.0);

            // Rb,Sr

        case 37:
            return _getOcc5s(0.5 + nelec);

        case 38:
            return _getOcc5s(1.0 + nelec);

            // Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd (39-48)

        case 39:
            return _getOcc4d(0.1 + nelec / 5.0);

        case 40:
            return _getOcc4d(0.2 + nelec / 5.0);

        case 41:
            return _getOcc4d(0.3 + nelec / 5.0);

        case 42:
            return _getOcc4d(0.4 + nelec / 5.0);

        case 43:
            return _getOcc4d(0.5 + nelec / 5.0);

        case 44:
            return _getOcc4d(0.6 + nelec / 5.0);

        case 45:
            return _getOcc4d(0.7 + nelec / 5.0);

        case 46:
            return _getOcc4d(0.8 + nelec / 5.0);

        case 47:
            return _getOcc4d(0.9 + nelec / 5.0);

        case 48:
            return _getOcc4d(1.0 + nelec / 5.0);

            // In,Sn,Sb,Te,I,Xe

        case 49:
            return _getOcc5s5p(0.375 + nelec / 4.0);

        case 50:
            return _getOcc5s5p(0.500 + nelec / 4.0);

        case 51:
            return _getOcc5s5p(0.625 + nelec / 4.0);

        case 52:
            return _getOcc5s5p(0.750 + nelec / 4.0);

        case 53:
            return _getOcc5s5p(0.875 + nelec / 4.0);

        case 54:
            return _getOcc5s5p(1.000 + nelec / 4.0);

            // Cs,Ba (55-56)

        case 55:
            return _getOcc6s(0.5 + nelec);

        case 56:
            return _getOcc6s(1.0 + nelec);

            // La,Ce,Pr,Nd,Pm,Sm,Eu,Gb,Tb,Dy,Ho,Er,Tm,Yb

        case 57:
            return _getOcc4f( 1.0 / 14.0 + nelec / 7.0);

        case 58:
            return _getOcc4f( 2.0 / 14.0 + nelec / 7.0);

        case 59:
            return _getOcc4f( 3.0 / 14.0 + nelec / 7.0);

        case 60:
            return _getOcc4f( 4.0 / 14.0 + nelec / 7.0);

        case 61:
            return _getOcc4f( 5.0 / 14.0 + nelec / 7.0);

        case 62:
            return _getOcc4f( 6.0 / 14.0 + nelec / 7.0);

        case 63:
            return _getOcc4f( 7.0 / 14.0 + nelec / 7.0);

        case 64:
            return _getOcc4f( 8.0 / 14.0 + nelec / 7.0);

        case 65:
            return _getOcc4f( 9.0 / 14.0 + nelec / 7.0);

        case 66:
            return _getOcc4f(10.0 / 14.0 + nelec / 7.0);

        case 67:
            return _getOcc4f(11.0 / 14.0 + nelec / 7.0);

        case 68:
            return _getOcc4f(12.0 / 14.0 + nelec / 7.0);

        case 69:
            return _getOcc4f(13.0 / 14.0 + nelec / 7.0);

        case 70:
            return _getOcc4f(14.0 / 14.0 + nelec / 7.0);

            // Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg

        case 71:
            return _getOcc5d(0.1 + nelec / 5.0);

        case 72:
            return _getOcc5d(0.2 + nelec / 5.0);

        case 73:
            return _getOcc5d(0.3 + nelec / 5.0);

        case 74:
            return _getOcc5d(0.4 + nelec / 5.0);

        case 75:
            return _getOcc5d(0.5 + nelec / 5.0);

        case 76:
            return _getOcc5d(0.6 + nelec / 5.0);

        case 77:
            return _getOcc5d(0.7 + nelec / 5.0);

        case 78:
            return _getOcc5d(0.8 + nelec / 5.0);

        case 79:
            return _getOcc5d(0.9 + nelec / 5.0);

        case 80:
            return _getOcc5d(1.0 + nelec / 5.0);

            // Tl,Pb,Bi,Po,At,Rn (81-86)

        case 81:
            return _getOcc6s6p(0.375 + nelec / 4.0);

        case 82:
            return _getOcc6s6p(0.500 + nelec / 4.0);

        case 83:
            return _getOcc6s6p(0.625 + nelec / 4.0);

        case 84:
            return _getOcc6s6p(0.750 + nelec / 4.0);

        case 85:
            return _getOcc6s6p(0.875 + nelec / 4.0);

        case 86:
            return _getOcc6s6p(1.000 + nelec / 4.0);

            // Fr,Ra (87-88)

        case 87:
            return _getOcc7s(0.5 + nelec);

        case 88:
            return _getOcc7s(1.0 + nelec);

            // Ac,Th,Pa,U,Np,Pu,Am,Cm (89-96)

        case 89:
            return _getOcc5f( 1.0 / 14.0 + nelec / 7.0);

        case 90:
            return _getOcc5f( 2.0 / 14.0 + nelec / 7.0);

        case 91:
            return _getOcc5f( 3.0 / 14.0 + nelec / 7.0);

        case 92:
            return _getOcc5f( 4.0 / 14.0 + nelec / 7.0);

        case 93:
            return _getOcc5f( 5.0 / 14.0 + nelec / 7.0);

        case 94:
            return _getOcc5f( 6.0 / 14.0 + nelec / 7.0);

        case 95:
            return _getOcc5f( 7.0 / 14.0 + nelec / 7.0);

        case 96:
            return _getOcc5f( 8.0 / 14.0 + nelec / 7.0);

            // 

        default:
            return std::vector<double>();
    }
}

std::vector<std::vector<double>>
CSADGuessDriver::getAlphaBetaOccupationNumbersForMolecule(const CMolecule& molecule,
                                                          const double     netCharge,
                                                          const double     numberOfUnpairedElectrons) const
{
    std::vector<std::vector<double>> occ_numbers;

    auto natoms = molecule.getNumberOfAtoms();

    auto partialcharges = parchg::getPartialCharges(molecule, netCharge);

    auto idselem = molecule.getIdsElemental();

    int32_t sum_id_elems = 0, sum_unpaired_electrons_on_atoms = 0;

    bool use_hint_for_unpaired_electrons = (static_cast<int32_t>(_numberOfUnpairedElectronsOnAtoms.size()) == natoms);

    for (int32_t i = 0; i < natoms; i++)
    {
        sum_id_elems += idselem[i];

        if (use_hint_for_unpaired_electrons)
        {
            sum_unpaired_electrons_on_atoms += _numberOfUnpairedElectronsOnAtoms[i];
        }
    }

    for (int32_t i = 0; i < natoms; i++)
    {
        auto weight = static_cast<double>(idselem[i]) / sum_id_elems;

        auto alpha_elec = 0.5 * (numberOfUnpairedElectrons - sum_unpaired_electrons_on_atoms) * weight;
        auto beta_elec = -0.5 * (numberOfUnpairedElectrons - sum_unpaired_electrons_on_atoms) * weight;

        if (use_hint_for_unpaired_electrons)
        {
            alpha_elec += 0.5 * _numberOfUnpairedElectronsOnAtoms[i];
            beta_elec -= 0.5 * _numberOfUnpairedElectronsOnAtoms[i];
        }

        auto alpha_occ = getOccupationNumbersForElement(idselem[i], -partialcharges[i] * 0.5 + alpha_elec);
        auto beta_occ = getOccupationNumbersForElement(idselem[i], -partialcharges[i] * 0.5 + beta_elec);

        occ_numbers.push_back(alpha_occ);
        occ_numbers.push_back(beta_occ);
    }

    return occ_numbers;
}

std::vector<std::vector<int32_t>>
CSADGuessDriver::getAOIndicesOfAtoms(const CMolecule&       molecule,
                                     const CMolecularBasis& basis) const
{
    std::vector<std::vector<int32_t>> aoinds_atoms;

    int32_t natoms = molecule.getNumberOfAtoms();

    for (int32_t atomidx = 0; atomidx < natoms; atomidx++)
    {
        aoinds_atoms.push_back(std::vector<int32_t>());
    }

    int32_t max_angl = basis.getMolecularMaxAngularMomentum(molecule);

    for (int32_t aoidx = 0, angl = 0; angl <= max_angl; angl++)
    {
        for (int32_t s = -angl; s <= angl; s++)
        {
            for (int32_t atomidx = 0; atomidx < natoms; atomidx++)
            {
                int32_t idelem = molecule.getIdsElemental()[atomidx];

                int32_t nao = basis.getNumberOfBasisFunctions(idelem, angl);

                for (int32_t i = 0; i < nao; i++, aoidx++)
                {
                    aoinds_atoms[atomidx].push_back(aoidx);
                }
            }
        }
    }

    return aoinds_atoms;
}

void
CSADGuessDriver::setNumberOfUnpairedElectronsOnAtoms(const std::vector<double>& numUnpairedElectrons)
{
    _numberOfUnpairedElectronsOnAtoms.clear();

    for (int32_t i = 0; i < static_cast<int32_t>(numUnpairedElectrons.size()); i++)
    {
        _numberOfUnpairedElectronsOnAtoms.push_back(numUnpairedElectrons[i]);
    }
}

CAODensityMatrix
CSADGuessDriver::compute(const CMolecule&       molecule,
                         const CMolecularBasis& basis_1,
                         const CMolecularBasis& basis_2,
                         const std::string&     densityType) const
{
    CAODensityMatrix dsad;

    if (_locRank == mpi::master())
    {
        // generate SAD guess

        COverlapIntegralsDriver ovldrv(_locComm);

        auto S12 = ovldrv.compute(molecule, basis_1, basis_2);

        auto S22 = ovldrv.compute(molecule, basis_2);

        dsad = _compSADGuess(molecule, basis_1, basis_2, S12, S22, densityType);
    }

    return dsad;
}

CAODensityMatrix
CSADGuessDriver::_compSADGuess(const CMolecule&       molecule,
                               const CMolecularBasis& basis_1,
                               const CMolecularBasis& basis_2,
                               const COverlapMatrix&  S12,
                               const COverlapMatrix&  S22,
                               const std::string&     densityType) const
{
    auto natoms = molecule.getNumberOfAtoms();

    auto nao_1 = S12.getNumberOfRows();

    auto nao_2 = S12.getNumberOfColumns();

    // sanity checks

    std::string err_ovl_size("SADGuessDriver.compute: Mismatch between overlap matrices");

    errors::assertMsgCritical(nao_2 == S22.getNumberOfRows(), err_ovl_size);

    errors::assertMsgCritical(nao_2 == S22.getNumberOfColumns(), err_ovl_size);

    // AO indices for atoms

    auto aoinds_atoms_1 = getAOIndicesOfAtoms(molecule, basis_1);

    auto aoinds_atoms_2 = getAOIndicesOfAtoms(molecule, basis_2);

    // more sanity checks

    int32_t count_ao_1 = 0;

    int32_t count_ao_2 = 0;

    for (int32_t atomidx = 0; atomidx < natoms; atomidx++)
    {
        count_ao_1 += static_cast<int32_t>(aoinds_atoms_1[atomidx].size());

        count_ao_2 += static_cast<int32_t>(aoinds_atoms_2[atomidx].size());
    }

    std::string err_bas_size("SADGuessDriver.compute: Mismatch between basis set and overlap matrix");

    errors::assertMsgCritical(count_ao_1 == nao_1 && count_ao_2 == nao_2, err_bas_size);

    errors::assertMsgCritical(count_ao_2 == S22.getNumberOfRows(), err_bas_size);

    // number of excessive electrons

    double charge = molecule.getCharge();

    double mult_1 = static_cast<double>(molecule.getMultiplicity() - 1);

    // occupation numbers

    auto alpha_beta_occ = getAlphaBetaOccupationNumbersForMolecule(molecule, charge, mult_1);

    std::vector<std::vector<double>> alpha_occ, beta_occ;

    for (int32_t i = 0; i < natoms; i++)
    {
        alpha_occ.push_back(alpha_beta_occ[i * 2 + 0]);
        beta_occ.push_back(alpha_beta_occ[i * 2 + 1]);
    }

    // C_SAD matrix

    CDenseMatrix csad_alpha(nao_2, nao_1);

    CDenseMatrix csad_beta(nao_2, nao_1);

    csad_alpha.zero();

    csad_beta.zero();

    #pragma omp parallel for schedule(dynamic)
    for (int32_t atomidx = 0; atomidx < natoms; atomidx++)
    {
        // AO indices for this atom

        const std::vector<int32_t>& aoinds_1 = aoinds_atoms_1[atomidx];

        const std::vector<int32_t>& aoinds_2 = aoinds_atoms_2[atomidx];

        // set up AO indices dimensions

        auto naodim_1 = static_cast<int32_t>(aoinds_1.size());

        auto naodim_2 = static_cast<int32_t>(aoinds_2.size());

        // size checking

        std::string err_ao_size("SADGuessDriver.compute: Mismatch between basis set and occupation number");

        errors::assertMsgCritical(alpha_occ[atomidx].size() == aoinds_1.size(), err_ao_size);

        errors::assertMsgCritical(beta_occ[atomidx].size() == aoinds_1.size(), err_ao_size);

        // atomic block of AOs

        CDenseMatrix block_12(naodim_1, naodim_2);

        CDenseMatrix block_22(naodim_2, naodim_2);

        block_12.zero();

        block_22.zero();

        for (int32_t i = 0; i < naodim_1; i++)
        {
            for (int32_t j = 0; j < naodim_2; j++)
            {
                block_12.values()[i * naodim_2 + j] =

                    S12.values()[aoinds_1[i] * nao_2 + aoinds_2[j]];
            }
        }

        for (int32_t i = 0; i < naodim_2; i++)
        {
            for (int32_t j = 0; j < naodim_2; j++)
            {
                block_22.values()[i * naodim_2 + j] =

                    S22.values()[aoinds_2[i] * nao_2 + aoinds_2[j]];
            }
        }

        // A = S12' C1(identity)

        CDenseMatrix mat_c1(naodim_1, naodim_1);

        mat_c1.zero();

        for (int32_t i = 0; i < naodim_1; i++)
        {
            mat_c1.values()[i * naodim_1 + i] = 1.0;
        }

        auto mat_a = denblas::multAtB(block_12, mat_c1);

        // S22^-1

        CDenseDiagonalizer diagdrv;

        diagdrv.diagonalize(block_22);

        std::string err_diag("SADGuessDriver.compute: Matrix diagonalization failed");

        errors::assertMsgCritical(diagdrv.getState(), err_diag);

        auto block_22_inv = diagdrv.getInvertedMatrix();

        // M = A' S22^-1 A

        auto prod = denblas::multAB(block_22_inv, mat_a);

        auto mat_m = denblas::multAtB(mat_a, prod);

        // M^-1/2

        diagdrv.diagonalize(mat_m);

        errors::assertMsgCritical(diagdrv.getState(), err_diag);

        auto mat_m_invsqrt = diagdrv.getInvertedSqrtMatrix();

        // C2 = S22^-1 A M^-1/2

        prod = denblas::multAB(mat_a, mat_m_invsqrt);

        auto mat_c2 = denblas::multAB(block_22_inv, prod);

        // update csad_alpha

        for (int32_t j = 0; j < naodim_2; j++)
        {
            for (int32_t i = 0; i < naodim_1; i++)
            {
                csad_alpha.values()[aoinds_2[j] * nao_1 + aoinds_1[i]] =

                    mat_c2.values()[j * naodim_1 + i] * sqrt(alpha_occ[atomidx][i]);
            }
        }

        // update csad_beta

        if (fstr::upcase(densityType) == std::string("UNRESTRICTED"))
        {
            for (int32_t j = 0; j < naodim_2; j++)
            {
                for (int32_t i = 0; i < naodim_1; i++)
                {
                    csad_beta.values()[aoinds_2[j] * nao_1 + aoinds_1[i]] =

                        mat_c2.values()[j * naodim_1 + i] * sqrt(beta_occ[atomidx][i]);
                }
            }
        }
    }

    // D_SAD density matrix

    std::vector<CDenseMatrix> dsad;

    if (fstr::upcase(densityType) == std::string("RESTRICTED"))
    {
        dsad.push_back(denblas::multABt(csad_alpha, csad_alpha));

        return CAODensityMatrix(dsad, denmat::rest);
    }
    else if (fstr::upcase(densityType) == std::string("UNRESTRICTED"))
    {
        dsad.push_back(denblas::multABt(csad_alpha, csad_alpha));

        dsad.push_back(denblas::multABt(csad_beta, csad_beta));

        return CAODensityMatrix(dsad, denmat::unrest);
    }
    else
    {
        return CAODensityMatrix();
    }
}
