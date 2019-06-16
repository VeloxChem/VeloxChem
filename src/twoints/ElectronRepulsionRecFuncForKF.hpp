//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include<cstdint>

#include "MemBlock2D.hpp"
#include "GtoPairsBlock.hpp"
#include "RecursionMap.hpp"

namespace erirecfunc { // erirecfunc namespace

    /**
    Computes batch of primitive (SK|G|SF) electron repulsion integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param wpDistances the vector of distances R(WP) = W - P.
    @param braGtoPairsBlock the GTOs psirs block on bra side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
    @param iContrPair the index of contracted GTOs apir on bra side.
    */
    void compElectronRepulsionForSKSF(      CMemBlock2D<double>& primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair);

    /**
    Computes sub-batch (0,45) of primitive (SK|G|SF) electron repulsion integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param wpDistances the vector of distances R(WP) = W - P.
    @param braGtoPairsBlock the GTOs psirs block on bra side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
    @param iContrPair the index of contracted GTOs apir on bra side.
    */
    void compElectronRepulsionForSKSF_0_45(      CMemBlock2D<double>& primBuffer,
                                           const CRecursionMap&       recursionMap,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& wpDistances,
                                           const CGtoPairsBlock&      braGtoPairsBlock,
                                           const CGtoPairsBlock&      ketGtoPairsBlock,
                                           const int32_t              nKetPrimPairs,
                                           const int32_t              iContrPair);

    /**
    Computes sub-batch (45,90) of primitive (SK|G|SF) electron repulsion integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param wpDistances the vector of distances R(WP) = W - P.
    @param braGtoPairsBlock the GTOs psirs block on bra side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
    @param iContrPair the index of contracted GTOs apir on bra side.
    */
    void compElectronRepulsionForSKSF_45_90(      CMemBlock2D<double>& primBuffer,
                                            const CRecursionMap&       recursionMap,
                                            const CMemBlock2D<double>& osFactors,
                                            const CMemBlock2D<double>& wpDistances,
                                            const CGtoPairsBlock&      braGtoPairsBlock,
                                            const CGtoPairsBlock&      ketGtoPairsBlock,
                                            const int32_t              nKetPrimPairs,
                                            const int32_t              iContrPair);

    /**
    Computes sub-batch (90,135) of primitive (SK|G|SF) electron repulsion integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param wpDistances the vector of distances R(WP) = W - P.
    @param braGtoPairsBlock the GTOs psirs block on bra side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
    @param iContrPair the index of contracted GTOs apir on bra side.
    */
    void compElectronRepulsionForSKSF_90_135(      CMemBlock2D<double>& primBuffer,
                                             const CRecursionMap&       recursionMap,
                                             const CMemBlock2D<double>& osFactors,
                                             const CMemBlock2D<double>& wpDistances,
                                             const CGtoPairsBlock&      braGtoPairsBlock,
                                             const CGtoPairsBlock&      ketGtoPairsBlock,
                                             const int32_t              nKetPrimPairs,
                                             const int32_t              iContrPair);

    /**
    Computes sub-batch (135,180) of primitive (SK|G|SF) electron repulsion integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param wpDistances the vector of distances R(WP) = W - P.
    @param braGtoPairsBlock the GTOs psirs block on bra side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
    @param iContrPair the index of contracted GTOs apir on bra side.
    */
    void compElectronRepulsionForSKSF_135_180(      CMemBlock2D<double>& primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (180,225) of primitive (SK|G|SF) electron repulsion integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param wpDistances the vector of distances R(WP) = W - P.
    @param braGtoPairsBlock the GTOs psirs block on bra side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
    @param iContrPair the index of contracted GTOs apir on bra side.
    */
    void compElectronRepulsionForSKSF_180_225(      CMemBlock2D<double>& primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (225,270) of primitive (SK|G|SF) electron repulsion integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param wpDistances the vector of distances R(WP) = W - P.
    @param braGtoPairsBlock the GTOs psirs block on bra side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
    @param iContrPair the index of contracted GTOs apir on bra side.
    */
    void compElectronRepulsionForSKSF_225_270(      CMemBlock2D<double>& primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (270,315) of primitive (SK|G|SF) electron repulsion integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param wpDistances the vector of distances R(WP) = W - P.
    @param braGtoPairsBlock the GTOs psirs block on bra side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
    @param iContrPair the index of contracted GTOs apir on bra side.
    */
    void compElectronRepulsionForSKSF_270_315(      CMemBlock2D<double>& primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (315,360) of primitive (SK|G|SF) electron repulsion integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param wpDistances the vector of distances R(WP) = W - P.
    @param braGtoPairsBlock the GTOs psirs block on bra side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
    @param iContrPair the index of contracted GTOs apir on bra side.
    */
    void compElectronRepulsionForSKSF_315_360(      CMemBlock2D<double>& primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);


} // erirecfunc namespace

