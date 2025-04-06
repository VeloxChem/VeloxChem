//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef DensityGridCubic_hpp
#define DensityGridCubic_hpp

#include "DensityGrid.hpp"
#include "DensityGridData2D.hpp"
#include "XCFunctionalType.hpp"

/**
 Class CDensityGridCubic class implements density grid.
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
    int _nDensityMatrices;

    /**
     The density variables values at grid points.
     */
    CDensityGridData2D _densityValues;

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
    CDensityGridCubic(const CDensityGridData2D& densityValues, const dengrid gridType);

    /**
     Creates a density grid object.

     @param nGridPoints the number of grid points with density values.
     @param nDensityMatrices the number of AO density matrices.
     @param xcFuncType the type of exchange-correlation functional.
     @param gridType the type of density grid.
     */
    CDensityGridCubic(const int nGridPoints, const int nDensityMatrices, const xcfun xcFuncType, const dengrid gridType);

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
    int getNumberOfGridPoints() const;

    /**
     Gets number of densitu matrices associated with density grid.

     @return the number of density matrices.
     */
    int getNumberOfDensityMatrices() const;

    /**
     Gets type of density grid type object.

     @return the type of density grid.
     */
    dengrid getDensityGridType() const;

    const double* pi(const int iDensityMatrix) const;
    double*       pi(const int iDensityMatrix);

    const double* piX(const int iDensityMatrix) const;
    double*       piX(const int iDensityMatrix);
    const double* piY(const int iDensityMatrix) const;
    double*       piY(const int iDensityMatrix);
    const double* piZ(const int iDensityMatrix) const;
    double*       piZ(const int iDensityMatrix);

    const double* gam2(const int iDensityMatrix) const;
    double*       gam2(const int iDensityMatrix);

    const double* gam2X(const int iDensityMatrix) const;
    double*       gam2X(const int iDensityMatrix);
    const double* gam2Y(const int iDensityMatrix) const;
    double*       gam2Y(const int iDensityMatrix);
    const double* gam2Z(const int iDensityMatrix) const;
    double*       gam2Z(const int iDensityMatrix);

    const double* gam2XX(const int iDensityMatrix) const;
    double*       gam2XX(const int iDensityMatrix);
    const double* gam2XY(const int iDensityMatrix) const;
    double*       gam2XY(const int iDensityMatrix);
    const double* gam2XZ(const int iDensityMatrix) const;
    double*       gam2XZ(const int iDensityMatrix);

    const double* gam2YX(const int iDensityMatrix) const;
    double*       gam2YX(const int iDensityMatrix);
    const double* gam2YY(const int iDensityMatrix) const;
    double*       gam2YY(const int iDensityMatrix);
    const double* gam2YZ(const int iDensityMatrix) const;
    double*       gam2YZ(const int iDensityMatrix);

    const double* gam2ZX(const int iDensityMatrix) const;
    double*       gam2ZX(const int iDensityMatrix);
    const double* gam2ZY(const int iDensityMatrix) const;
    double*       gam2ZY(const int iDensityMatrix);
    const double* gam2ZZ(const int iDensityMatrix) const;
    double*       gam2ZZ(const int iDensityMatrix);

    const double* rt_gam2(const int iDensityMatrix) const;
    double*       rt_gam2(const int iDensityMatrix);
    const double* rl_gam2(const int iDensityMatrix) const;
    double*       rl_gam2(const int iDensityMatrix);
    const double* tt_gam2(const int iDensityMatrix) const;
    double*       tt_gam2(const int iDensityMatrix);
    const double* tl_gam2(const int iDensityMatrix) const;
    double*       tl_gam2(const int iDensityMatrix);
    const double* ll_gam2(const int iDensityMatrix) const;
    double*       ll_gam2(const int iDensityMatrix);

    const double* st_gam2X(const int iDensityMatrix) const;
    double*       st_gam2X(const int iDensityMatrix);
    const double* st_gam2Y(const int iDensityMatrix) const;
    double*       st_gam2Y(const int iDensityMatrix);
    const double* st_gam2Z(const int iDensityMatrix) const;
    double*       st_gam2Z(const int iDensityMatrix);

    const double* sl_gam2X(const int iDensityMatrix) const;
    double*       sl_gam2X(const int iDensityMatrix);
    const double* sl_gam2Y(const int iDensityMatrix) const;
    double*       sl_gam2Y(const int iDensityMatrix);
    const double* sl_gam2Z(const int iDensityMatrix) const;
    double*       sl_gam2Z(const int iDensityMatrix);

    const double* rrt_pi(const int iDensityMatrix) const;
    double*       rrt_pi(const int iDensityMatrix);
    const double* rrl_pi(const int iDensityMatrix) const;
    double*       rrl_pi(const int iDensityMatrix);
    const double* rtt_pi(const int iDensityMatrix) const;
    double*       rtt_pi(const int iDensityMatrix);
    const double* rtl_pi(const int iDensityMatrix) const;
    double*       rtl_pi(const int iDensityMatrix);
    const double* rll_pi(const int iDensityMatrix) const;
    double*       rll_pi(const int iDensityMatrix);
    const double* ttt_pi(const int iDensityMatrix) const;
    double*       ttt_pi(const int iDensityMatrix);
    const double* ttl_pi(const int iDensityMatrix) const;
    double*       ttl_pi(const int iDensityMatrix);
    const double* tll_pi(const int iDensityMatrix) const;
    double*       tll_pi(const int iDensityMatrix);
    const double* lll_pi(const int iDensityMatrix) const;
    double*       lll_pi(const int iDensityMatrix);

    const double* rt_gam(const int iDensityMatrix) const;
    double*       rt_gam(const int iDensityMatrix);
    const double* rl_gam(const int iDensityMatrix) const;
    double*       rl_gam(const int iDensityMatrix);
    const double* tt_gam(const int iDensityMatrix) const;
    double*       tt_gam(const int iDensityMatrix);
    const double* tl_gam(const int iDensityMatrix) const;
    double*       tl_gam(const int iDensityMatrix);
    const double* ll_gam(const int iDensityMatrix) const;
    double*       ll_gam(const int iDensityMatrix);

    const double* piXX(const int iDensityMatrix) const;
    double*       piXX(const int iDensityMatrix);
    const double* piXY(const int iDensityMatrix) const;
    double*       piXY(const int iDensityMatrix);
    const double* piXZ(const int iDensityMatrix) const;
    double*       piXZ(const int iDensityMatrix);

    const double* piYX(const int iDensityMatrix) const;
    double*       piYX(const int iDensityMatrix);
    const double* piYY(const int iDensityMatrix) const;
    double*       piYY(const int iDensityMatrix);
    const double* piYZ(const int iDensityMatrix) const;
    double*       piYZ(const int iDensityMatrix);

    const double* piZX(const int iDensityMatrix) const;
    double*       piZX(const int iDensityMatrix);
    const double* piZY(const int iDensityMatrix) const;
    double*       piZY(const int iDensityMatrix);
    const double* piZZ(const int iDensityMatrix) const;
    double*       piZZ(const int iDensityMatrix);

    const double* piXXX(const int iDensityMatrix) const;
    double*       piXXX(const int iDensityMatrix);
    const double* piXXY(const int iDensityMatrix) const;
    double*       piXXY(const int iDensityMatrix);
    const double* piXXZ(const int iDensityMatrix) const;
    double*       piXXZ(const int iDensityMatrix);

    const double* piXYX(const int iDensityMatrix) const;
    double*       piXYX(const int iDensityMatrix);
    const double* piXYY(const int iDensityMatrix) const;
    double*       piXYY(const int iDensityMatrix);
    const double* piXYZ(const int iDensityMatrix) const;
    double*       piXYZ(const int iDensityMatrix);

    const double* piXZX(const int iDensityMatrix) const;
    double*       piXZX(const int iDensityMatrix);
    const double* piXZY(const int iDensityMatrix) const;
    double*       piXZY(const int iDensityMatrix);
    const double* piXZZ(const int iDensityMatrix) const;
    double*       piXZZ(const int iDensityMatrix);

    const double* piYXX(const int iDensityMatrix) const;
    double*       piYXX(const int iDensityMatrix);
    const double* piYXY(const int iDensityMatrix) const;
    double*       piYXY(const int iDensityMatrix);
    const double* piYXZ(const int iDensityMatrix) const;
    double*       piYXZ(const int iDensityMatrix);

    const double* piYYX(const int iDensityMatrix) const;
    double*       piYYX(const int iDensityMatrix);
    const double* piYYY(const int iDensityMatrix) const;
    double*       piYYY(const int iDensityMatrix);
    const double* piYYZ(const int iDensityMatrix) const;
    double*       piYYZ(const int iDensityMatrix);

    const double* piYZX(const int iDensityMatrix) const;
    double*       piYZX(const int iDensityMatrix);
    const double* piYZY(const int iDensityMatrix) const;
    double*       piYZY(const int iDensityMatrix);
    const double* piYZZ(const int iDensityMatrix) const;
    double*       piYZZ(const int iDensityMatrix);

    const double* piZXX(const int iDensityMatrix) const;
    double*       piZXX(const int iDensityMatrix);
    const double* piZXY(const int iDensityMatrix) const;
    double*       piZXY(const int iDensityMatrix);
    const double* piZXZ(const int iDensityMatrix) const;
    double*       piZXZ(const int iDensityMatrix);

    const double* piZYX(const int iDensityMatrix) const;
    double*       piZYX(const int iDensityMatrix);
    const double* piZYY(const int iDensityMatrix) const;
    double*       piZYY(const int iDensityMatrix);
    const double* piZYZ(const int iDensityMatrix) const;
    double*       piZYZ(const int iDensityMatrix);

    const double* piZZX(const int iDensityMatrix) const;
    double*       piZZX(const int iDensityMatrix);
    const double* piZZY(const int iDensityMatrix) const;
    double*       piZZY(const int iDensityMatrix);
    const double* piZZZ(const int iDensityMatrix) const;
    double*       piZZZ(const int iDensityMatrix);

    const double* gam(const int iDensityMatrix) const;
    double*       gam(const int iDensityMatrix);

    const double* st_gamX(const int iDensityMatrix) const;
    double*       st_gamX(const int iDensityMatrix);
    const double* st_gamY(const int iDensityMatrix) const;
    double*       st_gamY(const int iDensityMatrix);
    const double* st_gamZ(const int iDensityMatrix) const;
    double*       st_gamZ(const int iDensityMatrix);

    const double* sl_gamX(const int iDensityMatrix) const;
    double*       sl_gamX(const int iDensityMatrix);
    const double* sl_gamY(const int iDensityMatrix) const;
    double*       sl_gamY(const int iDensityMatrix);
    const double* sl_gamZ(const int iDensityMatrix) const;
    double*       sl_gamZ(const int iDensityMatrix);

    const double* gamXX(const int iDensityMatrix) const;
    double*       gamXX(const int iDensityMatrix);
    const double* gamXY(const int iDensityMatrix) const;
    double*       gamXY(const int iDensityMatrix);
    const double* gamXZ(const int iDensityMatrix) const;
    double*       gamXZ(const int iDensityMatrix);

    const double* gamYX(const int iDensityMatrix) const;
    double*       gamYX(const int iDensityMatrix);
    const double* gamYY(const int iDensityMatrix) const;
    double*       gamYY(const int iDensityMatrix);
    const double* gamYZ(const int iDensityMatrix) const;
    double*       gamYZ(const int iDensityMatrix);

    const double* gamZX(const int iDensityMatrix) const;
    double*       gamZX(const int iDensityMatrix);
    const double* gamZY(const int iDensityMatrix) const;
    double*       gamZY(const int iDensityMatrix);
    const double* gamZZ(const int iDensityMatrix) const;
    double*       gamZZ(const int iDensityMatrix);

    const double* gamX(const int iDensityMatrix) const;
    double*       gamX(const int iDensityMatrix);
    const double* gamY(const int iDensityMatrix) const;
    double*       gamY(const int iDensityMatrix);
    const double* gamZ(const int iDensityMatrix) const;
    double*       gamZ(const int iDensityMatrix);

    const double* rt_piX(const int iDensityMatrix) const;
    double*       rt_piX(const int iDensityMatrix);
    const double* rt_piY(const int iDensityMatrix) const;
    double*       rt_piY(const int iDensityMatrix);
    const double* rt_piZ(const int iDensityMatrix) const;
    double*       rt_piZ(const int iDensityMatrix);

    const double* rl_piX(const int iDensityMatrix) const;
    double*       rl_piX(const int iDensityMatrix);
    const double* rl_piY(const int iDensityMatrix) const;
    double*       rl_piY(const int iDensityMatrix);
    const double* rl_piZ(const int iDensityMatrix) const;
    double*       rl_piZ(const int iDensityMatrix);

    const double* ll_piX(const int iDensityMatrix) const;
    double*       ll_piX(const int iDensityMatrix);
    const double* ll_piY(const int iDensityMatrix) const;
    double*       ll_piY(const int iDensityMatrix);
    const double* ll_piZ(const int iDensityMatrix) const;
    double*       ll_piZ(const int iDensityMatrix);

    const double* tt_piX(const int iDensityMatrix) const;
    double*       tt_piX(const int iDensityMatrix);
    const double* tt_piY(const int iDensityMatrix) const;
    double*       tt_piY(const int iDensityMatrix);
    const double* tt_piZ(const int iDensityMatrix) const;
    double*       tt_piZ(const int iDensityMatrix);

    const double* tl_piX(const int iDensityMatrix) const;
    double*       tl_piX(const int iDensityMatrix);
    const double* tl_piY(const int iDensityMatrix) const;
    double*       tl_piY(const int iDensityMatrix);
    const double* tl_piZ(const int iDensityMatrix) const;
    double*       tl_piZ(const int iDensityMatrix);

    const double* l_piXX(const int iDensityMatrix) const;
    double*       l_piXX(const int iDensityMatrix);
    const double* l_piXY(const int iDensityMatrix) const;
    double*       l_piXY(const int iDensityMatrix);
    const double* l_piXZ(const int iDensityMatrix) const;
    double*       l_piXZ(const int iDensityMatrix);

    const double* l_piYX(const int iDensityMatrix) const;
    double*       l_piYX(const int iDensityMatrix);
    const double* l_piYY(const int iDensityMatrix) const;
    double*       l_piYY(const int iDensityMatrix);
    const double* l_piYZ(const int iDensityMatrix) const;
    double*       l_piYZ(const int iDensityMatrix);

    const double* l_piZX(const int iDensityMatrix) const;
    double*       l_piZX(const int iDensityMatrix);
    const double* l_piZY(const int iDensityMatrix) const;
    double*       l_piZY(const int iDensityMatrix);
    const double* l_piZZ(const int iDensityMatrix) const;
    double*       l_piZZ(const int iDensityMatrix);

    const double* t_piXX(const int iDensityMatrix) const;
    double*       t_piXX(const int iDensityMatrix);
    const double* t_piXY(const int iDensityMatrix) const;
    double*       t_piXY(const int iDensityMatrix);
    const double* t_piXZ(const int iDensityMatrix) const;
    double*       t_piXZ(const int iDensityMatrix);

    const double* t_piYX(const int iDensityMatrix) const;
    double*       t_piYX(const int iDensityMatrix);
    const double* t_piYY(const int iDensityMatrix) const;
    double*       t_piYY(const int iDensityMatrix);
    const double* t_piYZ(const int iDensityMatrix) const;
    double*       t_piYZ(const int iDensityMatrix);

    const double* t_piZX(const int iDensityMatrix) const;
    double*       t_piZX(const int iDensityMatrix);
    const double* t_piZY(const int iDensityMatrix) const;
    double*       t_piZY(const int iDensityMatrix);
    const double* t_piZZ(const int iDensityMatrix) const;
    double*       t_piZZ(const int iDensityMatrix);

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
     @param cubeMode a string that specifies which densities should be combined.
     */
    void DensityProd(const CDensityGrid& rwDensityGrid,
                     const CDensityGrid& rw2DensityGrid,
                     const xcfun         xcFuncType,
                     const int           numdens,
                     const std::string&  cubeMode);

    void DensityProdForLDA_TPA(const CDensityGrid& rwDensityGrid,
                               const CDensityGrid& rw2DensityGrid,
                               const int           numdens);

    void DensityProdForLDA_CRF(const CDensityGrid& rwDensityGrid,
                               const CDensityGrid& rw2DensityGrid,
                               const int           numdens);

    void DensityProdForLDA_3PA(const CDensityGrid& rwDensityGrid,
                               const CDensityGrid& rw2DensityGrid,
                               const int           numdens);

    void DensityProdForGGA_TPA(const CDensityGrid& rwDensityGrid,
                               const CDensityGrid& rw2DensityGrid,
                               const int           numdens);

    void DensityProdForGGA_CRF(const CDensityGrid& rwDensityGrid,
                               const CDensityGrid& rw2DensityGrid,
                               const int           numdens);

    void DensityProdForGGA_3PA(const CDensityGrid& rwDensityGrid,
                               const CDensityGrid& rw2DensityGrid,
                               const int           numdens);

    void DensityProdForMGGA_TPA(const CDensityGrid& rwDensityGrid,
                                const CDensityGrid& rw2DensityGrid,
                                const int           numdens);

    void DensityProdForMGGA_CRF(const CDensityGrid& rwDensityGrid,
                                const CDensityGrid& rw2DensityGrid,
                                const int           numdens);

    void DensityProdForMGGA_3PA(const CDensityGrid& rwDensityGrid,
                                const CDensityGrid& rw2DensityGrid,
                                const int           numdens);
};

#endif /* DensityGridQuad_hpp */
