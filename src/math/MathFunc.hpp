//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MathFunc_hpp
#define MathFunc_hpp

#include <cmath>
#include <cstdint>

namespace mathfunc {  // mathfunc namespace

/**
 Sets all elements of real numbers vector to zero.

 @param vector the vector of real numbers.
 @param nElements the number of elements in vector.
 */
void zero(double* vector, const int32_t nElements);

/**
 Sets all elements of integer numbers vector to zero.

 @param vector the vector of integer numbers.
 @param nElements the number of elements in vector.
 */
void zero(int32_t* vector, const int32_t nElements);

/**
 Sets all elements of real numbers vector to specific value.

 @param vector the vector of real numbers.
 @param value the value of element.
 @param nElements the number of elements in vector.
 */
void set_to(double* vector, const double value, const int32_t nElements);

/**
 Sets all elements of integer numbers vector to specific value.

 @param vector the vector of integer numbers.
 @param value the value of element.
 @param nElements the number of elements in vector.
 */
void set_to(int32_t* vector, const int32_t value, const int32_t nElements);

/**
 Computes sum of all elements in real numbers vector.

 @param vector the vector of real numbers.
 @param nElements the number of elements in vector.
 @return sum of all elements in vector.
 */
double sum(const double* vector, const int32_t nElements);

/**
 Computes sum of all elements in integer numbers vector.

 @param vector the vector of integer numbers.
 @param nElements the number of elements in vector.
 @return sum of all elements in vector.
 */
int32_t sum(const int32_t* vector, const int32_t nElements);

/**
 Scales all elements of real numbers vector by specific factor.

 @param vector the vector of real numbers.
 @param factor the scaling factor.
 @param nElements the number of elements in vector.
 */
void scale(double* vector, const double factor, const int32_t nElements);

/**
 Adds vector scaled by factor to other vector i.e. va = va + f * vb.

 @param aVector the vector of real numbers.
 @param bVector the vector of real numbers.
 @param factor the scaling factor.
 @param nElements the number of elements in vector.
 */
void add_scaled(double* aVector, const double* bVector, const double factor, const int32_t nElements);

/**
 Determines largest element in real numbers vector.

 @param vector the vector of real numbers.
 @param nElements the number of elements in vector.
 @return the largest element in real numbers vector.
 */
double max(const double* vector, const int32_t nElements);

/**
 Determines largest element in integer numbers vector.

 @param vector the vector of integer numbers.
 @param nElements the number of elements in vector.
 @return the largest element in integer numbers vector.
 */
int32_t max(const int32_t* vector, const int32_t nElements);

/**
 Normalizes vector of real numbers.

 @param vector the vector of real numbers.
 @param nElements the number of elements in vector.
 */
void normalize(double* vector, const int32_t nElements);

/**
 Sets indexes vector using size vector.

 @param aVector the indexes vector.
 @param bVector the sizes vector.
 @param nElements the number of elements in vectors.
 */
void indexes(int32_t* aVector, const int32_t* bVector, const int32_t nElements);

/**
 Sets indexes vector using size vector and offset.

 @param aVector the indexes vector.
 @param bVector the sizes vector.
 @param offset the offset of first index.
 @param nElements the number of elements in vectors.
 */
void indexes(int32_t* aVector, const int32_t* bVector, const int32_t offset, const int32_t nElements);

/**
 Sets ordering vector for given vector of binary values (0 or 1) by storing
 all indexes of binary vector elements equal to 1.

 @param aVector the ordering vector.
 @param bVector the binary vector.
 @param nElements the number of elements in vectors.
 */
void ordering(int32_t* aVector, const int32_t* bVector, const int32_t nElements);

/**
 Computes distance between two 3D vectors.

 @param aCoordX the Cartesian X coordinate of first vector.
 @param aCoordY the Cartesian Y coordinate of first vector.
 @param aCoordZ the Cartesian Z coordinate of first vector.
 @param bCoordX the Cartesian X coordinate of second vector.
 @param bCoordY the Cartesian Y coordinate of second vector.
 @param bCoordZ the Cartesian Z coordinate of second vector.
 @return the distance between vectors.
 */
inline double
distance(const double aCoordX, const double aCoordY, const double aCoordZ, const double bCoordX, const double bCoordY, const double bCoordZ)
{
    auto rx = aCoordX - bCoordX;

    auto ry = aCoordY - bCoordY;

    auto rz = aCoordZ - bCoordZ;

    return std::sqrt(rx * rx + ry * ry + rz * rz);
}

/**
 Computes distances between reference point A and vector of B points.

 @param abDistancesX the vector of distances R(AB)_x = A_x - B_x.
 @param abDistancesY the vector of distances R(AB)_y = A_y - B_y.
 @param abDistancesZ the vector of distances R(AB)_z = A_z - B_z.
 @param aCoordX the Cartesian X coordinate of point A.
 @param aCoordY the Cartesian Y coordinate of point A.
 @param aCoordZ the Cartesian Z coordinate of point A.
 @param bCoordsX the vector of Cartesian X coordinates of points B.
 @param bCoordsY the vector of Cartesian Y coordinates of points B.
 @param bCoordsZ the vector of Cartesian Z coordinates of points B.
 @param nElements the number of points B.
 */
void distances(double*       abDistancesX,
               double*       abDistancesY,
               double*       abDistancesZ,
               const double  aCoordX,
               const double  aCoordY,
               const double  aCoordZ,
               const double* bCoordsX,
               const double* bCoordsY,
               const double* bCoordsZ,
               const int32_t nElements);

/**
 Computes Chebtshev quadrature of second kind in [-1,1] interval.

 @param coordinates the vector of quadature coordinates.
 @param weights the vector of quadrature weights.
 @param nPoints the number of points in quadrature.
 */
void quadChebyshevOfKindTwo(double* coordinates, double* weights, const int32_t nPoints);

/**
 Copies integer numbers from one vector to another vector.

 @param aVector the destination vector.
 @param aPosition the position of first copied element in destination vector.
 @param bVector the source vector.
 @param bPosition the position of first copied element in source vector.
 @param nElements the number of elements.
 */
void copy(int32_t* aVector, const int32_t aPosition, const int32_t* bVector, const int32_t bPosition, const int32_t nElements);

/**
 Copies real numbers from one vector to another vector.

 @param aVector the destination vector.
 @param aPosition the position of first copied element in destination vector.
 @param bVector the source vector.
 @param bPosition the position of first copied element in source vector.
 @param nElements the number of elements.
 */
void copy(double* aVector, const int32_t aPosition, const double* bVector, const int32_t bPosition, const int32_t nElements);

/**
 Determines maximum number of components for tensor of given order.

 @param order the order of tensor.
 @return the number of components.
 */
int32_t maxTensorComponents(const int32_t order);

}  // namespace mathfunc

#endif /* MathFunc_hpp */
