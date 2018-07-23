//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef OneIntsFunc_hpp
#define OneIntsFunc_hpp

#include <cstdint>
#include <cmath>

#include "MemBlock2D.hpp"
#include "GtoBlock.hpp"

namespace intsfunc { // intsfunc namespace
    
    /**
     Computes distances between specific contracted GTO on bra side with all
     contracted GTOs on ket side.

     @param abDistances the vector of Cartesian R(AB) = A - B distances.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesAB(      CMemBlock2D<double>& abDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto);
    
    /**
     Computes Obara-Saika factors for overlap integrals.

     @param osFactors the vector of Obara-Saika factors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compFactorsForOverlap(      CMemBlock2D<double>& osFactors,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto);
    
    /**
     Computes vector of distances between center P of combined primitive GTOs
     and center A of primitive GTO on bra side.
     
     @param paDistances the vector of Cartesian R(PA) = P - A distances.
     @param abDistances the vector of Cartesian R(AB) = A - B distances.
     @param osFactors the vector of Obara-Saika factors.
     @param nFactors the fundamental dimensions of Obara-Saika factors vectors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesPA(      CMemBlock2D<double>& paDistances,
                         const CMemBlock2D<double>& abDistances,
                         const CMemBlock2D<double>& osFactors,
                         const int32_t              nFactors,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto);
    
    /**
     Computes vector of distances between center P of combined primitive GTOs
     and center B of primitive GTO on ket side.
     
     @param pbDistances the vector of Cartesian R(PB) = P - B distances.
     @param abDistances the vector of Cartesian R(AB) = A - B distances.
     @param osFactors the vector of Obara-Saika factors.
     @param nFactors the fundamental dimensions of Obara-Saika factors vectors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesPB(      CMemBlock2D<double>& pbDistances,
                         const CMemBlock2D<double>& abDistances,
                         const CMemBlock2D<double>& osFactors,
                         const int32_t              nFactors,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto);
    
    /**
     Computes vectors of Xi and Zeta factors for combination of specific
     primitive Gaussian exponent on bra side and vector of primitive Gaussian
     exponents on ket side.
     
     Xi = 1 / (e_bra + e_ket) and Zeta = e_bra * e_ket / (e_bra + e_ket)
     
     @param factorsXi the vector of Xi factors.
     @param factorsZeta the vector of Zeta factors.
     @param braExponent the primitive Gaussian exponent on bra side.
     @param ketExponents the vector of primitive Gaussian exponents on ket side.
     @param nElements the number of elements in vector of primitive Gaussian
            exponents on ket side.
     */
    void compXiAndZeta(      double* factorsXi,
                             double* factorsZeta,
                       const double  braExponent,
                       const double* ketExponents,
                       const int32_t nElements);
    
    /**
     Computes vector of distances between center P of combined primitive
     Gaussian and center A of primitive Gaussian on bra side.

     @param xDistancesPA the vector of distances RPA_x = P_x - A_x.
     @param yDistancesPA the vector of distances RPA_y = P_y - A_y.
     @param zDistancesPA the vector of distances RPA_z = P_z - A_z.
     @param ketExponents the vector of primitive Gaussian exponents on ket side.
     @param factorsXi the vector of Xi factors.
     @param xDistancesAB the vector of distances RAB_x = A_x - B_x.
     @param yDistancesAB the vector of distances RAB_y = A_y - B_y.
     @param zDistancesAB the vector of distances RAB_z = A_z - B_z.
     @param nElements the number of elements in vector of primitive Gaussian
            exponents on ket side.
     */
    void compDistancesPA(      double* xDistancesPA,
                               double* yDistancesPA,
                               double* zDistancesPA,
                         const double* ketExponents,
                         const double* factorsXi,
                         const double* xDistancesAB,
                         const double* yDistancesAB,
                         const double* zDistancesAB,
                         const int32_t nElements);
    
    /**
     Computes vector of distances between center P of combined primitive
     Gaussian and center B of primitive Gaussian on ket side.
     
     @param xDistancesPB the vector of distances RPB_x = P_x - B_x.
     @param yDistancesPB the vector of distances RPB_y = P_y - B_y.
     @param zDistancesPB the vector of distances RPB_z = P_z - B_z.
     @prama braExponent the primitive Gaussian exponent on bra side.
     @param factorsXi the vector of Xi factors.
     @param xDistancesAB the vector of distances RAB_x = A_x - B_x.
     @param yDistancesAB the vector of distances RAB_y = A_y - B_y.
     @param zDistancesAB the vector of distances RAB_z = A_z - B_z.
     @param nElements the number of elements in vector of primitive Gaussian
     exponents on ket side.
     */
    void compDistancesPB(      double* xDistancesPB,
                               double* yDistancesPB,
                               double* zDistancesPB,
                         const double  braExponent,
                         const double* factorsXi,
                         const double* xDistancesAB,
                         const double* yDistancesAB,
                         const double* zDistancesAB,
                         const int32_t nElements);
    
} // intsfunc namespace

#endif /* OneIntsFunc_hpp */
