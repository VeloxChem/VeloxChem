//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include <cstdint>

#include "GtoBlock.hpp"
#include "MemBlock2D.hpp"
#include "RecursionMap.hpp"

namespace amomrecfunc {  // amomrecfunc namespace

/**
Computes batch of primitive (D|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDD(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto);

/**
Computes sub-batch (0,36) of primitive (D|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDD_0_36(CMemBlock2D<double>&       primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

/**
Computes sub-batch (36,72) of primitive (D|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDD_36_72(CMemBlock2D<double>&       primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

/**
Computes sub-batch (72,108) of primitive (D|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDD_72_108(CMemBlock2D<double>&       primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto);

/**
Computes batch of primitive (D|L|F) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDF(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto);

/**
Computes sub-batch (0,45) of primitive (D|L|F) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDF_0_45(CMemBlock2D<double>&       primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

/**
Computes sub-batch (45,90) of primitive (D|L|F) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDF_45_90(CMemBlock2D<double>&       primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

/**
Computes sub-batch (90,135) of primitive (D|L|F) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDF_90_135(CMemBlock2D<double>&       primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto);

/**
Computes sub-batch (135,180) of primitive (D|L|F) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDF_135_180(CMemBlock2D<double>&       primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

/**
Computes batch of primitive (F|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForFD(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto);

/**
Computes sub-batch (0,45) of primitive (F|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForFD_0_45(CMemBlock2D<double>&       primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

/**
Computes sub-batch (45,90) of primitive (F|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForFD_45_90(CMemBlock2D<double>&       primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

/**
Computes sub-batch (90,135) of primitive (F|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForFD_90_135(CMemBlock2D<double>&       primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto);

/**
Computes sub-batch (135,180) of primitive (F|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForFD_135_180(CMemBlock2D<double>&       primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

/**
Computes batch of primitive (D|L|G) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDG(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto);

/**
Computes sub-batch (0,45) of primitive (D|L|G) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDG_0_45(CMemBlock2D<double>&       primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

/**
Computes sub-batch (45,90) of primitive (D|L|G) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDG_45_90(CMemBlock2D<double>&       primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

/**
Computes sub-batch (90,135) of primitive (D|L|G) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDG_90_135(CMemBlock2D<double>&       primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto);

/**
Computes sub-batch (135,180) of primitive (D|L|G) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDG_135_180(CMemBlock2D<double>&       primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

/**
Computes sub-batch (180,225) of primitive (D|L|G) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDG_180_225(CMemBlock2D<double>&       primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

/**
Computes sub-batch (225,270) of primitive (D|L|G) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForDG_225_270(CMemBlock2D<double>&       primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

/**
Computes batch of primitive (G|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGD(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto);

/**
Computes sub-batch (0,45) of primitive (G|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGD_0_45(CMemBlock2D<double>&       primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

/**
Computes sub-batch (45,90) of primitive (G|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGD_45_90(CMemBlock2D<double>&       primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

/**
Computes sub-batch (90,135) of primitive (G|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGD_90_135(CMemBlock2D<double>&       primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto);

/**
Computes sub-batch (135,180) of primitive (G|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGD_135_180(CMemBlock2D<double>&       primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

/**
Computes sub-batch (180,225) of primitive (G|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGD_180_225(CMemBlock2D<double>&       primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

/**
Computes sub-batch (225,270) of primitive (G|L|D) angular momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compAngularMomentumForGD_225_270(CMemBlock2D<double>&       primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

}  // namespace amomrecfunc
