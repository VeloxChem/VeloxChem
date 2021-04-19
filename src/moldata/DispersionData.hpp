//
//                           VELOXCHEM 1.0-RC
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
//
//  This file contains derivative work of dftd4 (v2.4.0):
//  Copyright © 2017-2019 Stefan Grimme, Sebastian Ehlert, Eike Caldeweyher

#ifndef DispersionData_hpp
#define DispersionData_hpp

#include <cstdint>
#include <vector>

namespace dispdata {  // dispdata namespace

/**
 Creates the refn data (reference: dftd4 v2.4.0).

 @return a vector of refn data.
 */
std::vector<int32_t> getRefN();

/**
 Creates the Zeff data (reference: dftd4 v2.4.0).

 @return a vector of Zeff data.
 */
std::vector<int32_t> getZeff();

/**
 Creates the chemical hardness data (reference: dftd4 v2.4.0).

 @return a vector of chemical hardness data.
 */
std::vector<double> getChemicalHardness();

/**
 Creates the r4r2 (sqrt_z_r4_over_r2) data (reference: dftd4 v2.4.0).

 @return a vector of r4r2 data.
 */
std::vector<double> getR4R2();

/**
 Creates the refsys data (reference: dftd4 v2.4.0).

 @return a 2d vector of refsys data.
 */
std::vector<std::vector<int32_t>> getRefSys();

/**
 Creates the clsh data (reference: dftd4 v2.4.0).

 @return a 2d vector of clsh data.
 */
std::vector<std::vector<double>> getClsH();

/**
 Creates the clsq data (reference: dftd4 v2.4.0).

 @return a 2d vector of clsq data.
 */
std::vector<std::vector<double>> getClsQ();

/**
 Creates the ascale data (reference: dftd4 v2.4.0).

 @return a 2d vector of ascale data.
 */
std::vector<std::vector<double>> getAscale();

/**
 Creates the refcn data (reference: dftd4 v2.4.0).

 @return a 2d vector of refcn data.
 */
std::vector<std::vector<double>> getRefCN();

/**
 Creates the refcovcn data (reference: dftd4 v2.4.0).

 @return a 2d vector of refcovcn data.
 */
std::vector<std::vector<double>> getRefCovCN();

/**
 Creates the hcount data (reference: dftd4 v2.4.0).

 @return a 2d vector of hcount data.
 */
std::vector<std::vector<double>> getHcount();

/**
 Creates the alphaiw data (reference: dftd4 v2.4.0).

 @return a 3d vector of alphaiw data.
 */
std::vector<std::vector<std::vector<double>>> getAlphaiw();

/**
 Gets the sscale value (reference: dftd4 v2.4.0).

 @param index the index.
 @return the sscale value.
 */
double getSscale(int32_t index);

/**
 Gets the secaiw data (reference: dftd4 v2.4.0).

 @param index the index.
 @return a vector of secaiw data.
 */
std::vector<double> getSecaiw(int32_t index);

}  // namespace dispdata

#endif /* DispersionData_hpp */
