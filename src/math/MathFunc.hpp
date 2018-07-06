//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef MathFunc_hpp
#define MathFunc_hpp

#include <cstdint>
#include <cmath>

namespace mathfunc { // mathfunc namespace

    /**
     Sets all elements of real numbers vector to zero.

     @param vector the vector of real numbers.
     @param nElements the number of elements in vector.
     */
    void zero(      double* vector,
              const int32_t nElements);
    
    /**
     Sets all elements of integer numbers vector to zero.
     
     @param vector the vector of integer numbers.
     @param nElements the number of elements in vector.
     */
    void zero(      int32_t* vector,
              const int32_t nElements);

    /**
     Sets all elements of real numbers vector to specific value.

     @param vector the vector of real numbers.
     @param value the value of element.
     @param nElements the number of elements in vector.
     */
    void set_to(      double* vector,
                const double  value,
                const int32_t nElements);
    
    /**
     Sets all elements of integer numbers vector to specific value.
     
     @param vector the vector of integer numbers.
     @param value the value of element.
     @param nElements the number of elements in vector.
     */
    void set_to(      int32_t* vector,
                const int32_t  value,
                const int32_t  nElements);

    /**
     Computes sum of all elements in real numbers vector.

     @param vector the vector of real numbers.
     @param nElements the number of elements in vector.
     @return sum of all elements in vector.
     */
    double sum(const double* vector,
               const int32_t nElements);

    /**
     Computes sum of all elements in integer numbers vector.

     @param vector the vector of integer numbers.
     @param nElements the number of elements in vector.
     @return sum of all elements in vector.
     */
    int32_t sum(const int32_t* vector,
                const int32_t  nElements);
    
    /**
     Determines largest element in real numbers vector.

     @param vector the vector of real numbers.
     @param nElements the number of elements in vector.
     @return the largest element in real numbers vector.
     */
    double max(const double* vector,
               const int32_t nElements);
    
    /**
     Determines largest element in integer numbers vector.
     
     @param vector the vector of integer numbers.
     @param nElements the number of elements in vector.
     @return the largest element in integer numbers vector.
     */
    int32_t max(const int32_t* vector,
                const int32_t  nElements);

    /**
     Normalizes vector of real numbers.

     @param vector the vector of real numbers.
     @param nElements the number of elements in vector.
     */
    void normalize(      double* vector,
                   const int32_t nElements);
    
    /**
     Sets indexes vector using size vector.

     @param aVector the indexes vector.
     @param bVector the sizes vector.
     @param nElements the number of elements in vectors.
     */
    void indexes(      int32_t* aVector,
                 const int32_t* bVector,
                 const int32_t  nElements);

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
    inline double distance(const double aCoordX,
                           const double aCoordY,
                           const double aCoordZ,
                           const double bCoordX,
                           const double bCoordY,
                           const double bCoordZ)
    {
        auto rx = aCoordX - bCoordX;

        auto ry = aCoordY - bCoordY;

        auto rz = aCoordZ - bCoordZ;

        return std::sqrt(rx * rx  + ry * ry + rz * rz);
    }
    
    /**
     Computes Chebtshev quadrature of second kind in [-1,1] interval.

     @param coordinates the vector of quadature coordinates.
     @param weights the vector of quadrature weights.
     @param nPoints the number of points in quadrature.
     */
    void quadChebyshevOfKindTwo(      double* coordinates,
                                      double* weights,
                                const int32_t nPoints);

} // mathfunc namespace

#endif /* MathFunc_hpp */
