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
    Computes batch of primitive (SL|G|SH) electron repulsion integrals and stores
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
    void compElectronRepulsionForSLSH(      CMemBlock2D<double>* primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair);

    /**
    Computes sub-batch (0,95) of primitive (SL|G|SH) electron repulsion integrals and stores
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
    void compElectronRepulsionForSLSH_0_95(      CMemBlock2D<double>* primBuffer,
                                           const CRecursionMap&       recursionMap,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& wpDistances,
                                           const CGtoPairsBlock&      braGtoPairsBlock,
                                           const CGtoPairsBlock&      ketGtoPairsBlock,
                                           const int32_t              nKetPrimPairs,
                                           const int32_t              iContrPair);

    /**
    Computes sub-batch (95,190) of primitive (SL|G|SH) electron repulsion integrals and stores
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
    void compElectronRepulsionForSLSH_95_190(      CMemBlock2D<double>* primBuffer,
                                             const CRecursionMap&       recursionMap,
                                             const CMemBlock2D<double>& osFactors,
                                             const CMemBlock2D<double>& wpDistances,
                                             const CGtoPairsBlock&      braGtoPairsBlock,
                                             const CGtoPairsBlock&      ketGtoPairsBlock,
                                             const int32_t              nKetPrimPairs,
                                             const int32_t              iContrPair);

    /**
    Computes sub-batch (190,285) of primitive (SL|G|SH) electron repulsion integrals and stores
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
    void compElectronRepulsionForSLSH_190_285(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (285,380) of primitive (SL|G|SH) electron repulsion integrals and stores
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
    void compElectronRepulsionForSLSH_285_380(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (380,475) of primitive (SL|G|SH) electron repulsion integrals and stores
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
    void compElectronRepulsionForSLSH_380_475(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (475,569) of primitive (SL|G|SH) electron repulsion integrals and stores
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
    void compElectronRepulsionForSLSH_475_569(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (569,663) of primitive (SL|G|SH) electron repulsion integrals and stores
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
    void compElectronRepulsionForSLSH_569_663(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (663,757) of primitive (SL|G|SH) electron repulsion integrals and stores
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
    void compElectronRepulsionForSLSH_663_757(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (757,851) of primitive (SL|G|SH) electron repulsion integrals and stores
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
    void compElectronRepulsionForSLSH_757_851(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (851,945) of primitive (SL|G|SH) electron repulsion integrals and stores
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
    void compElectronRepulsionForSLSH_851_945(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);


} // erirecfunc namespace

