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

#ifndef DensityGridQuad_hpp
#define DensityGridQuad_hpp

#include <cstdint>
#include <vector>

#include "DensityGrid.hpp"
#include "DensityGridData2D.hpp"
#include "XCFunctionalType.hpp"

/**
 Class CDensityGridQuad class implements density grid.
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
    int _nDensityMatrices;

    /**
     The density variables values at grid points.
     */
    CDensityGridData2D _densityValues;

   public:
    /**
     Creates an empty density grid object.
     */
    CDensityGridQuad();

    /**
     Creates a density grid object.

     @param nGridPoints the number of grid points with density values.
     @param nDensityMatrices the number of AO density matrices.
     @param xcFuncType the type of exchange-correlation functional.
     @param gridType the type of density grid.
     */
    CDensityGridQuad(const int nGridPoints, const int nDensityMatrices, const xcfun xcFuncType, const dengrid gridType);

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

    const double* gam(const int iDensityMatrix) const;
    double*       gam(const int iDensityMatrix);

    const double* gamX(const int iDensityMatrix) const;
    double*       gamX(const int iDensityMatrix);
    const double* gamY(const int iDensityMatrix) const;
    double*       gamY(const int iDensityMatrix);
    const double* gamZ(const int iDensityMatrix) const;
    double*       gamZ(const int iDensityMatrix);

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
    void DensityProd(const CDensityGrid& rwDensityGrid, const xcfun xcFuncType, const int numDens, const std::string& quadMode);

    void DensityProdForLDA(const CDensityGrid& rwDensityGrid, const int numDens, const std::string& quadMode);
    void DensityProdForGGA(const CDensityGrid& rwDensityGrid, const int numDens, const std::string& quadMode);
    void DensityProdForMGGA(const CDensityGrid& rwDensityGrid, const int numDens, const std::string& quadMode);
};

#endif /* DensityGridQuad_hpp */
