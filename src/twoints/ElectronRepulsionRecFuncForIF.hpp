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
    Computes batch of primitive (SI|G|SF) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISF(      CMemBlock2D<double>& primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair);

    /**
    Computes sub-batch (0,47) of primitive (SI|G|SF) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISF_0_47(      CMemBlock2D<double>& primBuffer,
                                           const CRecursionMap&       recursionMap,
                                           const CMemBlock2D<double>& osFactors,
                                           const CMemBlock2D<double>& wpDistances,
                                           const CGtoPairsBlock&      braGtoPairsBlock,
                                           const CGtoPairsBlock&      ketGtoPairsBlock,
                                           const int32_t              nKetPrimPairs,
                                           const int32_t              iContrPair);

    /**
    Computes sub-batch (47,94) of primitive (SI|G|SF) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISF_47_94(      CMemBlock2D<double>& primBuffer,
                                            const CRecursionMap&       recursionMap,
                                            const CMemBlock2D<double>& osFactors,
                                            const CMemBlock2D<double>& wpDistances,
                                            const CGtoPairsBlock&      braGtoPairsBlock,
                                            const CGtoPairsBlock&      ketGtoPairsBlock,
                                            const int32_t              nKetPrimPairs,
                                            const int32_t              iContrPair);

    /**
    Computes sub-batch (94,141) of primitive (SI|G|SF) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISF_94_141(      CMemBlock2D<double>& primBuffer,
                                             const CRecursionMap&       recursionMap,
                                             const CMemBlock2D<double>& osFactors,
                                             const CMemBlock2D<double>& wpDistances,
                                             const CGtoPairsBlock&      braGtoPairsBlock,
                                             const CGtoPairsBlock&      ketGtoPairsBlock,
                                             const int32_t              nKetPrimPairs,
                                             const int32_t              iContrPair);

    /**
    Computes sub-batch (141,188) of primitive (SI|G|SF) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISF_141_188(      CMemBlock2D<double>& primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (188,234) of primitive (SI|G|SF) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISF_188_234(      CMemBlock2D<double>& primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);

    /**
    Computes sub-batch (234,280) of primitive (SI|G|SF) electron repulsion integrals and stores
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
    void compElectronRepulsionForSISF_234_280(      CMemBlock2D<double>& primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& wpDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetPrimPairs,
                                              const int32_t              iContrPair);


} // erirecfunc namespace

