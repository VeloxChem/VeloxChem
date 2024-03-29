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

#ifndef DensityGridCubic_hpp
#define DensityGridCubic_hpp

#include <cstdint>
#include <ostream>
#include <vector>

#include "DensityGrid.hpp"
#include "DensityGridType.hpp"
#include "MemBlock2D.hpp"
#include "MolecularGrid.hpp"
#include "XCFuncType.hpp"

/**
 Class CDensityGridCubic class implements density grid.

 @author K. Ahmadzadeh
 */
class CDensityGridCubic
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
    CDensityGridCubic();

    /**
     Creates a density grid object.

     @param densityValues the 2D memory block object with density values data.
     @param gridType the type of density grid.
     */
    CDensityGridCubic(const CMemBlock2D<double>& densityValues, const dengrid gridType);

    /**
     Creates a density grid object.

     @param nGridPoints the number of grid points with density values.
     @param nDensityMatrices the number of AO density matrices.
     @param xcFuncType the type of exchange-correlation functional.
     @param gridType the type of density grid.
     */
    CDensityGridCubic(const int32_t nGridPoints, const int32_t nDensityMatrices, const xcfun xcFuncType, const dengrid gridType);

    /**
     Creates a density grid object by copying other density grid object.

     @param source the density grid object.
     */
    CDensityGridCubic(const CDensityGridCubic& source);

    /**
     Creates a density grid object by by moving other density grid object.

     @param source the density grid object.
     */
    CDensityGridCubic(CDensityGridCubic&& source) noexcept;

    /**
     Destroys a density grid object.
     */
    ~CDensityGridCubic();

    /**
     Assigns a density grid object by copying other density grid object.

     @param source the density grid object.
     */
    CDensityGridCubic& operator=(const CDensityGridCubic& source);

    /**
     Assigns a density grid object by moving other density grid object.

     @param source the density grid object.
     */
    CDensityGridCubic& operator=(CDensityGridCubic&& source) noexcept;

    /**
     Compares density grid object with other density grid object.

     @param other the density grid object.
     @return true if density grid objects are equal, false otherwise.
     */
    bool operator==(const CDensityGridCubic& other) const;

    /**
     Compares density grid object with other density grid object.

     @param other the density grid object.
     @return true if density grid objects are not equal, false otherwise.
     */
    bool operator!=(const CDensityGridCubic& other) const;

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

    const double* pi(const int32_t iDensityMatrix) const;
    double*       pi(const int32_t iDensityMatrix);

    const double* piX(const int32_t iDensityMatrix) const;
    double*       piX(const int32_t iDensityMatrix);
    const double* piY(const int32_t iDensityMatrix) const;
    double*       piY(const int32_t iDensityMatrix);
    const double* piZ(const int32_t iDensityMatrix) const;
    double*       piZ(const int32_t iDensityMatrix);

    const double* gam2(const int32_t iDensityMatrix) const;
    double*       gam2(const int32_t iDensityMatrix);

    const double* gam2X(const int32_t iDensityMatrix) const;
    double*       gam2X(const int32_t iDensityMatrix);
    const double* gam2Y(const int32_t iDensityMatrix) const;
    double*       gam2Y(const int32_t iDensityMatrix);
    const double* gam2Z(const int32_t iDensityMatrix) const;
    double*       gam2Z(const int32_t iDensityMatrix);

    const double* gam2XX(const int32_t iDensityMatrix) const;
    double*       gam2XX(const int32_t iDensityMatrix);
    const double* gam2XY(const int32_t iDensityMatrix) const;
    double*       gam2XY(const int32_t iDensityMatrix);
    const double* gam2XZ(const int32_t iDensityMatrix) const;
    double*       gam2XZ(const int32_t iDensityMatrix);

    const double* gam2YX(const int32_t iDensityMatrix) const;
    double*       gam2YX(const int32_t iDensityMatrix);
    const double* gam2YY(const int32_t iDensityMatrix) const;
    double*       gam2YY(const int32_t iDensityMatrix);
    const double* gam2YZ(const int32_t iDensityMatrix) const;
    double*       gam2YZ(const int32_t iDensityMatrix);

    const double* gam2ZX(const int32_t iDensityMatrix) const;
    double*       gam2ZX(const int32_t iDensityMatrix);
    const double* gam2ZY(const int32_t iDensityMatrix) const;
    double*       gam2ZY(const int32_t iDensityMatrix);
    const double* gam2ZZ(const int32_t iDensityMatrix) const;
    double*       gam2ZZ(const int32_t iDensityMatrix);

    const double* rt_gam2(const int32_t iDensityMatrix) const;
    double*       rt_gam2(const int32_t iDensityMatrix);
    const double* rl_gam2(const int32_t iDensityMatrix) const;
    double*       rl_gam2(const int32_t iDensityMatrix);
    const double* tt_gam2(const int32_t iDensityMatrix) const;
    double*       tt_gam2(const int32_t iDensityMatrix);
    const double* tl_gam2(const int32_t iDensityMatrix) const;
    double*       tl_gam2(const int32_t iDensityMatrix);
    const double* ll_gam2(const int32_t iDensityMatrix) const;
    double*       ll_gam2(const int32_t iDensityMatrix);

    const double* st_gam2X(const int32_t iDensityMatrix) const;
    double*       st_gam2X(const int32_t iDensityMatrix);
    const double* st_gam2Y(const int32_t iDensityMatrix) const;
    double*       st_gam2Y(const int32_t iDensityMatrix);
    const double* st_gam2Z(const int32_t iDensityMatrix) const;
    double*       st_gam2Z(const int32_t iDensityMatrix);

    const double* sl_gam2X(const int32_t iDensityMatrix) const;
    double*       sl_gam2X(const int32_t iDensityMatrix);
    const double* sl_gam2Y(const int32_t iDensityMatrix) const;
    double*       sl_gam2Y(const int32_t iDensityMatrix);
    const double* sl_gam2Z(const int32_t iDensityMatrix) const;
    double*       sl_gam2Z(const int32_t iDensityMatrix);

    const double* rrt_pi(const int32_t iDensityMatrix) const;
    double*       rrt_pi(const int32_t iDensityMatrix);
    const double* rrl_pi(const int32_t iDensityMatrix) const;
    double*       rrl_pi(const int32_t iDensityMatrix);
    const double* rtt_pi(const int32_t iDensityMatrix) const;
    double*       rtt_pi(const int32_t iDensityMatrix);
    const double* rtl_pi(const int32_t iDensityMatrix) const;
    double*       rtl_pi(const int32_t iDensityMatrix);
    const double* rll_pi(const int32_t iDensityMatrix) const;
    double*       rll_pi(const int32_t iDensityMatrix);
    const double* ttt_pi(const int32_t iDensityMatrix) const;
    double*       ttt_pi(const int32_t iDensityMatrix);
    const double* ttl_pi(const int32_t iDensityMatrix) const;
    double*       ttl_pi(const int32_t iDensityMatrix);
    const double* tll_pi(const int32_t iDensityMatrix) const;
    double*       tll_pi(const int32_t iDensityMatrix);
    const double* lll_pi(const int32_t iDensityMatrix) const;
    double*       lll_pi(const int32_t iDensityMatrix);

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

    const double* piXX(const int32_t iDensityMatrix) const;
    double*       piXX(const int32_t iDensityMatrix);
    const double* piXY(const int32_t iDensityMatrix) const;
    double*       piXY(const int32_t iDensityMatrix);
    const double* piXZ(const int32_t iDensityMatrix) const;
    double*       piXZ(const int32_t iDensityMatrix);

    const double* piYX(const int32_t iDensityMatrix) const;
    double*       piYX(const int32_t iDensityMatrix);
    const double* piYY(const int32_t iDensityMatrix) const;
    double*       piYY(const int32_t iDensityMatrix);
    const double* piYZ(const int32_t iDensityMatrix) const;
    double*       piYZ(const int32_t iDensityMatrix);

    const double* piZX(const int32_t iDensityMatrix) const;
    double*       piZX(const int32_t iDensityMatrix);
    const double* piZY(const int32_t iDensityMatrix) const;
    double*       piZY(const int32_t iDensityMatrix);
    const double* piZZ(const int32_t iDensityMatrix) const;
    double*       piZZ(const int32_t iDensityMatrix);

    const double* piXXX(const int32_t iDensityMatrix) const;
    double*       piXXX(const int32_t iDensityMatrix);
    const double* piXXY(const int32_t iDensityMatrix) const;
    double*       piXXY(const int32_t iDensityMatrix);
    const double* piXXZ(const int32_t iDensityMatrix) const;
    double*       piXXZ(const int32_t iDensityMatrix);

    const double* piXYX(const int32_t iDensityMatrix) const;
    double*       piXYX(const int32_t iDensityMatrix);
    const double* piXYY(const int32_t iDensityMatrix) const;
    double*       piXYY(const int32_t iDensityMatrix);
    const double* piXYZ(const int32_t iDensityMatrix) const;
    double*       piXYZ(const int32_t iDensityMatrix);

    const double* piXZX(const int32_t iDensityMatrix) const;
    double*       piXZX(const int32_t iDensityMatrix);
    const double* piXZY(const int32_t iDensityMatrix) const;
    double*       piXZY(const int32_t iDensityMatrix);
    const double* piXZZ(const int32_t iDensityMatrix) const;
    double*       piXZZ(const int32_t iDensityMatrix);

    const double* piYXX(const int32_t iDensityMatrix) const;
    double*       piYXX(const int32_t iDensityMatrix);
    const double* piYXY(const int32_t iDensityMatrix) const;
    double*       piYXY(const int32_t iDensityMatrix);
    const double* piYXZ(const int32_t iDensityMatrix) const;
    double*       piYXZ(const int32_t iDensityMatrix);

    const double* piYYX(const int32_t iDensityMatrix) const;
    double*       piYYX(const int32_t iDensityMatrix);
    const double* piYYY(const int32_t iDensityMatrix) const;
    double*       piYYY(const int32_t iDensityMatrix);
    const double* piYYZ(const int32_t iDensityMatrix) const;
    double*       piYYZ(const int32_t iDensityMatrix);

    const double* piYZX(const int32_t iDensityMatrix) const;
    double*       piYZX(const int32_t iDensityMatrix);
    const double* piYZY(const int32_t iDensityMatrix) const;
    double*       piYZY(const int32_t iDensityMatrix);
    const double* piYZZ(const int32_t iDensityMatrix) const;
    double*       piYZZ(const int32_t iDensityMatrix);

    const double* piZXX(const int32_t iDensityMatrix) const;
    double*       piZXX(const int32_t iDensityMatrix);
    const double* piZXY(const int32_t iDensityMatrix) const;
    double*       piZXY(const int32_t iDensityMatrix);
    const double* piZXZ(const int32_t iDensityMatrix) const;
    double*       piZXZ(const int32_t iDensityMatrix);

    const double* piZYX(const int32_t iDensityMatrix) const;
    double*       piZYX(const int32_t iDensityMatrix);
    const double* piZYY(const int32_t iDensityMatrix) const;
    double*       piZYY(const int32_t iDensityMatrix);
    const double* piZYZ(const int32_t iDensityMatrix) const;
    double*       piZYZ(const int32_t iDensityMatrix);

    const double* piZZX(const int32_t iDensityMatrix) const;
    double*       piZZX(const int32_t iDensityMatrix);
    const double* piZZY(const int32_t iDensityMatrix) const;
    double*       piZZY(const int32_t iDensityMatrix);
    const double* piZZZ(const int32_t iDensityMatrix) const;
    double*       piZZZ(const int32_t iDensityMatrix);

    const double* gam(const int32_t iDensityMatrix) const;
    double*       gam(const int32_t iDensityMatrix);

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

    const double* gamX(const int32_t iDensityMatrix) const;
    double*       gamX(const int32_t iDensityMatrix);
    const double* gamY(const int32_t iDensityMatrix) const;
    double*       gamY(const int32_t iDensityMatrix);
    const double* gamZ(const int32_t iDensityMatrix) const;
    double*       gamZ(const int32_t iDensityMatrix);

    const double* rt_piX(const int32_t iDensityMatrix) const;
    double*       rt_piX(const int32_t iDensityMatrix);
    const double* rt_piY(const int32_t iDensityMatrix) const;
    double*       rt_piY(const int32_t iDensityMatrix);
    const double* rt_piZ(const int32_t iDensityMatrix) const;
    double*       rt_piZ(const int32_t iDensityMatrix);

    const double* rl_piX(const int32_t iDensityMatrix) const;
    double*       rl_piX(const int32_t iDensityMatrix);
    const double* rl_piY(const int32_t iDensityMatrix) const;
    double*       rl_piY(const int32_t iDensityMatrix);
    const double* rl_piZ(const int32_t iDensityMatrix) const;
    double*       rl_piZ(const int32_t iDensityMatrix);

    const double* ll_piX(const int32_t iDensityMatrix) const;
    double*       ll_piX(const int32_t iDensityMatrix);
    const double* ll_piY(const int32_t iDensityMatrix) const;
    double*       ll_piY(const int32_t iDensityMatrix);
    const double* ll_piZ(const int32_t iDensityMatrix) const;
    double*       ll_piZ(const int32_t iDensityMatrix);

    const double* tt_piX(const int32_t iDensityMatrix) const;
    double*       tt_piX(const int32_t iDensityMatrix);
    const double* tt_piY(const int32_t iDensityMatrix) const;
    double*       tt_piY(const int32_t iDensityMatrix);
    const double* tt_piZ(const int32_t iDensityMatrix) const;
    double*       tt_piZ(const int32_t iDensityMatrix);

    const double* tl_piX(const int32_t iDensityMatrix) const;
    double*       tl_piX(const int32_t iDensityMatrix);
    const double* tl_piY(const int32_t iDensityMatrix) const;
    double*       tl_piY(const int32_t iDensityMatrix);
    const double* tl_piZ(const int32_t iDensityMatrix) const;
    double*       tl_piZ(const int32_t iDensityMatrix);

    const double* l_piXX(const int32_t iDensityMatrix) const;
    double*       l_piXX(const int32_t iDensityMatrix);
    const double* l_piXY(const int32_t iDensityMatrix) const;
    double*       l_piXY(const int32_t iDensityMatrix);
    const double* l_piXZ(const int32_t iDensityMatrix) const;
    double*       l_piXZ(const int32_t iDensityMatrix);

    const double* l_piYX(const int32_t iDensityMatrix) const;
    double*       l_piYX(const int32_t iDensityMatrix);
    const double* l_piYY(const int32_t iDensityMatrix) const;
    double*       l_piYY(const int32_t iDensityMatrix);
    const double* l_piYZ(const int32_t iDensityMatrix) const;
    double*       l_piYZ(const int32_t iDensityMatrix);

    const double* l_piZX(const int32_t iDensityMatrix) const;
    double*       l_piZX(const int32_t iDensityMatrix);
    const double* l_piZY(const int32_t iDensityMatrix) const;
    double*       l_piZY(const int32_t iDensityMatrix);
    const double* l_piZZ(const int32_t iDensityMatrix) const;
    double*       l_piZZ(const int32_t iDensityMatrix);

    const double* t_piXX(const int32_t iDensityMatrix) const;
    double*       t_piXX(const int32_t iDensityMatrix);
    const double* t_piXY(const int32_t iDensityMatrix) const;
    double*       t_piXY(const int32_t iDensityMatrix);
    const double* t_piXZ(const int32_t iDensityMatrix) const;
    double*       t_piXZ(const int32_t iDensityMatrix);

    const double* t_piYX(const int32_t iDensityMatrix) const;
    double*       t_piYX(const int32_t iDensityMatrix);
    const double* t_piYY(const int32_t iDensityMatrix) const;
    double*       t_piYY(const int32_t iDensityMatrix);
    const double* t_piYZ(const int32_t iDensityMatrix) const;
    double*       t_piYZ(const int32_t iDensityMatrix);

    const double* t_piZX(const int32_t iDensityMatrix) const;
    double*       t_piZX(const int32_t iDensityMatrix);
    const double* t_piZY(const int32_t iDensityMatrix) const;
    double*       t_piZY(const int32_t iDensityMatrix);
    const double* t_piZZ(const int32_t iDensityMatrix) const;
    double*       t_piZZ(const int32_t iDensityMatrix);

    inline double
    prod3_r(double B_r, double B_i, double C_r, double C_i, double D_r, double D_i)
    {
        return (B_r * C_r * D_r - B_i * D_i * C_r - C_i * D_i * B_r - B_i * C_i * D_r);
    }

    inline double
    prod3_i(double B_r, double B_i, double C_r, double C_i, double D_r, double D_i)
    {
        return (-B_i * C_i * D_i + B_i * C_r * D_r + C_i * B_r * D_r + D_i * B_r * C_r);
    }

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

     @param densityGridab the screened density grid.
     @param molecularGridab the screened molecular grid.
     @param rwDensityGrid  one-time transformed densities evaluated on the grid
     @param densityThreshold the screening threshold for density values.
     @param numdens the number of densities.
     @param quadMode a string that specifies which densities should be combined.
     */
    void DensityProd(const CDensityGrid& rwDensityGrid,
                     const CDensityGrid& rw2DensityGrid,
                     const xcfun         xcFuncType,
                     const int32_t       numdens,
                     const std::string&  quadMode);

    /**
     Converts density grid object to text and insert it into output text
     stream.

     @param output the output text stream.
     @param source the density grid.
     */
    friend std::ostream& operator<<(std::ostream& output, const CDensityGridCubic& source);
};

#endif /* DensityGridQuad_hpp */
