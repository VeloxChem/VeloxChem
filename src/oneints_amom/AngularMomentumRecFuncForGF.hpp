//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include <cstdint>

#include "GtoBlock.hpp"
#include "MemBlock2D.hpp"
#include "RecursionMap.hpp"

namespace amomrecfunc {  // amomrecfunc namespace

/**
Computes batch of primitive (G|L|F) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGF(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto);

/**
Computes sub-batch (0,50) of primitive (G|L|F) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGF_0_50(CMemBlock2D<double>&       primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

/**
Computes sub-batch (50,100) of primitive (G|L|F) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGF_50_100(CMemBlock2D<double>&       primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto);

/**
Computes sub-batch (100,150) of primitive (G|L|F) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGF_100_150(CMemBlock2D<double>&       primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

/**
Computes sub-batch (150,200) of primitive (G|L|F) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGF_150_200(CMemBlock2D<double>&       primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

/**
Computes sub-batch (200,250) of primitive (G|L|F) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGF_200_250(CMemBlock2D<double>&       primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

/**
Computes sub-batch (250,300) of primitive (G|L|F) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGF_250_300(CMemBlock2D<double>&       primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

/**
Computes sub-batch (300,350) of primitive (G|L|F) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGF_300_350(CMemBlock2D<double>&       primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

/**
Computes sub-batch (350,400) of primitive (G|L|F) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGF_350_400(CMemBlock2D<double>&       primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

/**
Computes sub-batch (400,450) of primitive (G|L|F) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGF_400_450(CMemBlock2D<double>&       primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

}  // namespace amomrecfunc
