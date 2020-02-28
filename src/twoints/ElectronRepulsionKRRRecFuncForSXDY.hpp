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

namespace erikrrfunc { // erikrrfunc namespace

    /**
    Computes batch of contracted (SX|G|DD) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXDD(      CMemBlock2D<double>& ketBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& cdDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair);

    /**
    Computes batch of contracted (SX|G|DF) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXDF(      CMemBlock2D<double>& ketBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& cdDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair);

    /**
    Computes batch of contracted (SX|G|DG) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXDG(      CMemBlock2D<double>& ketBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& cdDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair);

    /**
    Computes batch of contracted (SX|G|DH) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXDH(      CMemBlock2D<double>& ketBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& cdDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair);

    /**
    Computes sub-batch (0,63) of contracted (SX|G|DH) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXDH_0_63(      CMemBlock2D<double>& ketBuffer,
                                           const CRecursionMap&       recursionMap,
                                           const CMemBlock2D<double>& cdDistances,
                                           const CGtoPairsBlock&      braGtoPairsBlock,
                                           const CGtoPairsBlock&      ketGtoPairsBlock,
                                           const int32_t              nKetContrPairs,
                                           const int32_t              iContrPair);

    /**
    Computes sub-batch (63,126) of contracted (SX|G|DH) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXDH_63_126(      CMemBlock2D<double>& ketBuffer,
                                             const CRecursionMap&       recursionMap,
                                             const CMemBlock2D<double>& cdDistances,
                                             const CGtoPairsBlock&      braGtoPairsBlock,
                                             const CGtoPairsBlock&      ketGtoPairsBlock,
                                             const int32_t              nKetContrPairs,
                                             const int32_t              iContrPair);

    /**
    Computes batch of contracted (SX|G|DI) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXDI(      CMemBlock2D<double>& ketBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& cdDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair);

    /**
    Computes sub-batch (0,84) of contracted (SX|G|DI) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXDI_0_84(      CMemBlock2D<double>& ketBuffer,
                                           const CRecursionMap&       recursionMap,
                                           const CMemBlock2D<double>& cdDistances,
                                           const CGtoPairsBlock&      braGtoPairsBlock,
                                           const CGtoPairsBlock&      ketGtoPairsBlock,
                                           const int32_t              nKetContrPairs,
                                           const int32_t              iContrPair);

    /**
    Computes sub-batch (84,168) of contracted (SX|G|DI) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXDI_84_168(      CMemBlock2D<double>& ketBuffer,
                                             const CRecursionMap&       recursionMap,
                                             const CMemBlock2D<double>& cdDistances,
                                             const CGtoPairsBlock&      braGtoPairsBlock,
                                             const CGtoPairsBlock&      ketGtoPairsBlock,
                                             const int32_t              nKetContrPairs,
                                             const int32_t              iContrPair);


} // erikrrfunc namespace

