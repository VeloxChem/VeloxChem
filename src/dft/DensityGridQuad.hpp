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

#ifndef DensityGridQuad_hpp
#define DensityGridQuad_hpp

#include <cstdint>
#include <ostream>
#include <vector>

#include "DensityGrid.hpp"
#include "DensityGridType.hpp"
#include "MemBlock2D.hpp"
#include "MolecularGrid.hpp"
#include "XCFuncType.hpp"

/**
 Class CDensityGridQuad class implements density grid.

 @author K. Ahmadzadeh
 */
class CDensityGridQuad
{
    /**
     The type of density grid.
     */
    dengrid _gridType;

    /**
     The number of associated density matrices.
     */
    int32_t _nDensityMatrices;

    /**
     The density variables values at grid points.
     */
    CMemBlock2D<double> _densityValues;

   public:
    /**
     Creates an empty density grid object.
     */
    CDensityGridQuad();

    /**
     Creates a density grid object.

     @param densityValues the 2D memory block object with density values data.
     @param gridType the type of density grid.
     */
    CDensityGridQuad(const CMemBlock2D<double>& densityValues, const dengrid gridType);

    /**
     Creates a density grid object.

     @param nGridPoints the number of grid points with density values.
     @param nDensityMatrices the number of AO density matrices.
     @param xcFuncType the type of exchange-correlation functional.
     @param gridType the type of density grid.
     */
    CDensityGridQuad(const int32_t nGridPoints, const int32_t nDensityMatrices, const xcfun xcFuncType, const dengrid gridType);

    /**
     Creates a density grid object by copying other density grid object.

     @param source the density grid object.
     */
    CDensityGridQuad(const CDensityGridQuad& source);

    /**
     Creates a density grid object by by moving other density grid object.

     @param source the density grid object.
     */
    CDensityGridQuad(CDensityGridQuad&& source) noexcept;

    /**
     Destroys a density grid object.
     */
    ~CDensityGridQuad();

    /**
     Assigns a density grid object by copying other density grid object.

     @param source the density grid object.
     */
    CDensityGridQuad& operator=(const CDensityGridQuad& source);

    /**
     Assigns a density grid object by moving other density grid object.

     @param source the density grid object.
     */
    CDensityGridQuad& operator=(CDensityGridQuad&& source) noexcept;

    /**
     Compares density grid object with other density grid object.

     @param other the density grid object.
     @return true if density grid objects are equal, false otherwise.
     */
    bool operator==(const CDensityGridQuad& other) const;

    /**
     Compares density grid object with other density grid object.

     @param other the density grid object.
     @return true if density grid objects are not equal, false otherwise.
     */
    bool operator!=(const CDensityGridQuad& other) const;

    /**
     Initialize density values at grid point to zero.
     */
    void zero();

    /**
     Gets number of grid points in density grid object.

     @return the number of grid points.
     */
    int32_t getNumberOfGridPoints() const;

    /**
     Gets number of densitu matrices associated with density grid.

     @return the number of density matrices.
     */
    int32_t getNumberOfDensityMatrices() const;

    /**
     Gets type of density grid type object.

     @return the type of density grid.
     */
    dengrid getDensityGridType() const;

    const double* gam(const int32_t iDensityMatrix) const;
    double*       gam(const int32_t iDensityMatrix);

    const double* gamX(const int32_t iDensityMatrix) const;
    double*       gamX(const int32_t iDensityMatrix);
    const double* gamY(const int32_t iDensityMatrix) const;
    double*       gamY(const int32_t iDensityMatrix);
    const double* gamZ(const int32_t iDensityMatrix) const;
    double*       gamZ(const int32_t iDensityMatrix);

    const double* gamXX(const int32_t iDensityMatrix) const;
    double*       gamXX(const int32_t iDensityMatrix);
    const double* gamXY(const int32_t iDensityMatrix) const;
    double*       gamXY(const int32_t iDensityMatrix);
    const double* gamXZ(const int32_t iDensityMatrix) const;
    double*       gamXZ(const int32_t iDensityMatrix);

    const double* gamYX(const int32_t iDensityMatrix) const;
    double*       gamYX(const int32_t iDensityMatrix);
    const double* gamYY(const int32_t iDensityMatrix) const;
    double*       gamYY(const int32_t iDensityMatrix);
    const double* gamYZ(const int32_t iDensityMatrix) const;
    double*       gamYZ(const int32_t iDensityMatrix);

    const double* gamZX(const int32_t iDensityMatrix) const;
    double*       gamZX(const int32_t iDensityMatrix);
    const double* gamZY(const int32_t iDensityMatrix) const;
    double*       gamZY(const int32_t iDensityMatrix);
    const double* gamZZ(const int32_t iDensityMatrix) const;
    double*       gamZZ(const int32_t iDensityMatrix);

    const double* rt_gam(const int32_t iDensityMatrix) const;
    double*       rt_gam(const int32_t iDensityMatrix);
    const double* rl_gam(const int32_t iDensityMatrix) const;
    double*       rl_gam(const int32_t iDensityMatrix);
    const double* tt_gam(const int32_t iDensityMatrix) const;
    double*       tt_gam(const int32_t iDensityMatrix);
    const double* tl_gam(const int32_t iDensityMatrix) const;
    double*       tl_gam(const int32_t iDensityMatrix);
    const double* ll_gam(const int32_t iDensityMatrix) const;
    double*       ll_gam(const int32_t iDensityMatrix);

    const double* st_gamX(const int32_t iDensityMatrix) const;
    double*       st_gamX(const int32_t iDensityMatrix);
    const double* st_gamY(const int32_t iDensityMatrix) const;
    double*       st_gamY(const int32_t iDensityMatrix);
    const double* st_gamZ(const int32_t iDensityMatrix) const;
    double*       st_gamZ(const int32_t iDensityMatrix);

    const double* sl_gamX(const int32_t iDensityMatrix) const;
    double*       sl_gamX(const int32_t iDensityMatrix);
    const double* sl_gamY(const int32_t iDensityMatrix) const;
    double*       sl_gamY(const int32_t iDensityMatrix);
    const double* sl_gamZ(const int32_t iDensityMatrix) const;
    double*       sl_gamZ(const int32_t iDensityMatrix);

    inline double
    prod2_r(double B_r, double B_i, double C_r, double C_i)
    {
        return (B_r * C_r - B_i * C_i);
    }

    inline double
    prod2_i(double B_r, double B_i, double C_r, double C_i)
    {
        return (B_i * C_r + B_r * C_i);
    }

    /**
     Generates products of one-time transformed densities to be used for quadratic response.

     @param rwDensityGrid one-time transformed densities evaluated on the grid
     @param xcFuncType the type of exchange-correlation functional.
     @param numDens the number of densities.
     @param quadMode a string that specifies which densities should be combined.
     */
    void DensityProd(const CDensityGrid& rwDensityGrid, const xcfun xcFuncType, const int32_t numDens, const std::string& quadMode);

    /**
     Converts density grid object to text and insert it into output text
     stream.

     @param output the output text stream.
     @param source the density grid.
     */
    friend std::ostream& operator<<(std::ostream& output, const CDensityGridQuad& source);
};

#endif /* DensityGridQuad_hpp */
