//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include<cstdint>

#include "MemBlock2D.hpp"
#include "GtoPairsBlock.hpp"
#include "RecursionMap.hpp"

namespace erirecfunc { // erirecfunc namespace

    /**
    Computes batch of primitive (SH|G|SI) electron repulsion integrals and stores
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
    void compElectronRepulsionForSHSI(      CMemBlock2D<double>* primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair);

    /**
    Computes sub-batch (0,98) of primitive (SH|G|SI) electron repulsion integrals and stores
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
    void compElectronRepulsionForSHSI_0_98(      CMemBlock2D<double>* primBuffer,
                                           const CRecursionMap&       recursionMap,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& wpDistances,
                                           const CGtoPairsBlock&      braGtoPairsBlock,
                                           const CGtoPairsBlock&      ketGtoPairsBlock,
                                           const int32_t              nKetPrimPairs,
                                           const int32_t              iContrPair);

    /**
    Computes sub-batch (98,196) of primitive (SH|G|SI) electron repulsion integrals and stores
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
    void compElectronRepulsionForSHSI_98_196(      CMemBlock2D<double>* primBuffer,
                                             const CRecursionMap&       recursionMap,
                                             const CMemBlock2D<double>& osFactors,
                                             const CMemBlock2D<double>& wpDistances,
                                             const CGtoPairsBlock&      braGtoPairsBlock,
                                             const CGtoPairsBlock&      ketGtoPairsBlock,
                                             const int32_t              nKetPrimPairs,
                                             const int32_t              iContrPair);

    /**
    Computes sub-batch (196,294) of primitive (SH|G|SI) electron repulsion integrals and stores
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
    void compElectronRepulsionForSHSI_196_294(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (294,392) of primitive (SH|G|SI) electron repulsion integrals and stores
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
    void compElectronRepulsionForSHSI_294_392(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (392,490) of primitive (SH|G|SI) electron repulsion integrals and stores
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
    void compElectronRepulsionForSHSI_392_490(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (490,588) of primitive (SH|G|SI) electron repulsion integrals and stores
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
    void compElectronRepulsionForSHSI_490_588(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);


} // erirecfunc namespace

