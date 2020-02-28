//
//                           VELOXCHEM 1.0-RC
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
    Computes batch of primitive (SI|G|SL) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISL(      CMemBlock2D<double>* primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair);

    /**
    Computes sub-batch (0,97) of primitive (SI|G|SL) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISL_0_97(      CMemBlock2D<double>* primBuffer,
                                           const CRecursionMap&       recursionMap,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& wpDistances,
                                           const CGtoPairsBlock&      braGtoPairsBlock,
                                           const CGtoPairsBlock&      ketGtoPairsBlock,
                                           const int32_t              nKetPrimPairs,
                                           const int32_t              iContrPair);

    /**
    Computes sub-batch (97,194) of primitive (SI|G|SL) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISL_97_194(      CMemBlock2D<double>* primBuffer,
                                             const CRecursionMap&       recursionMap,
                                             const CMemBlock2D<double>& osFactors,
                                             const CMemBlock2D<double>& wpDistances,
                                             const CGtoPairsBlock&      braGtoPairsBlock,
                                             const CGtoPairsBlock&      ketGtoPairsBlock,
                                             const int32_t              nKetPrimPairs,
                                             const int32_t              iContrPair);

    /**
    Computes sub-batch (194,291) of primitive (SI|G|SL) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISL_194_291(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (291,388) of primitive (SI|G|SL) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISL_291_388(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (388,485) of primitive (SI|G|SL) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISL_388_485(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (485,582) of primitive (SI|G|SL) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISL_485_582(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (582,679) of primitive (SI|G|SL) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISL_582_679(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (679,776) of primitive (SI|G|SL) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISL_679_776(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (776,873) of primitive (SI|G|SL) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISL_776_873(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (873,970) of primitive (SI|G|SL) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISL_873_970(      CMemBlock2D<double>* primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (970,1067) of primitive (SI|G|SL) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISL_970_1067(      CMemBlock2D<double>* primBuffer,
                                               const CRecursionMap&       recursionMap,
                                               const CMemBlock2D<double>& osFactors,
                                               const CMemBlock2D<double>& wpDistances,
                                               const CGtoPairsBlock&      braGtoPairsBlock,
                                               const CGtoPairsBlock&      ketGtoPairsBlock,
                                               const int32_t              nKetPrimPairs,
                                               const int32_t              iContrPair);

    /**
    Computes sub-batch (1067,1164) of primitive (SI|G|SL) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISL_1067_1164(      CMemBlock2D<double>* primBuffer,
                                                const CRecursionMap&       recursionMap,
                                                const CMemBlock2D<double>& osFactors,
                                                const CMemBlock2D<double>& wpDistances,
                                                const CGtoPairsBlock&      braGtoPairsBlock,
                                                const CGtoPairsBlock&      ketGtoPairsBlock,
                                                const int32_t              nKetPrimPairs,
                                                const int32_t              iContrPair);

    /**
    Computes sub-batch (1164,1260) of primitive (SI|G|SL) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISL_1164_1260(      CMemBlock2D<double>* primBuffer,
                                                const CRecursionMap&       recursionMap,
                                                const CMemBlock2D<double>& osFactors,
                                                const CMemBlock2D<double>& wpDistances,
                                                const CGtoPairsBlock&      braGtoPairsBlock,
                                                const CGtoPairsBlock&      ketGtoPairsBlock,
                                                const int32_t              nKetPrimPairs,
                                                const int32_t              iContrPair);


} // erirecfunc namespace

