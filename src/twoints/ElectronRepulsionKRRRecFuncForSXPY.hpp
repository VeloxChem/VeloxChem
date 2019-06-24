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

namespace erikrrfunc { // erikrrfunc namespace

    /**
    Computes batch of contracted (SX|G|PP) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXPP(      CMemBlock2D<double>& ketBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& cdDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair);

    /**
    Computes batch of contracted (SX|G|PD) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXPD(      CMemBlock2D<double>& ketBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& cdDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair);

    /**
    Computes batch of contracted (SX|G|PF) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXPF(      CMemBlock2D<double>& ketBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& cdDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair);

    /**
    Computes batch of contracted (SX|G|PG) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXPG(      CMemBlock2D<double>& ketBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& cdDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair);

    /**
    Computes batch of contracted (SX|G|PH) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXPH(      CMemBlock2D<double>& ketBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& cdDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair);

    /**
    Computes batch of contracted (SX|G|PI) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXPI(      CMemBlock2D<double>& ketBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& cdDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair);

    /**
    Computes batch of contracted (SX|G|PK) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXPK(      CMemBlock2D<double>& ketBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& cdDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair);

    /**
    Computes sub-batch (0,54) of contracted (SX|G|PK) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXPK_0_54(      CMemBlock2D<double>& ketBuffer,
                                           const CRecursionMap&       recursionMap,
                                           const CMemBlock2D<double>& cdDistances,
                                           const CGtoPairsBlock&      braGtoPairsBlock,
                                           const CGtoPairsBlock&      ketGtoPairsBlock,
                                           const int32_t              nKetContrPairs,
                                           const int32_t              iContrPair);

    /**
    Computes sub-batch (54,108) of contracted (SX|G|PK) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXPK_54_108(      CMemBlock2D<double>& ketBuffer,
                                             const CRecursionMap&       recursionMap,
                                             const CMemBlock2D<double>& cdDistances,
                                             const CGtoPairsBlock&      braGtoPairsBlock,
                                             const CGtoPairsBlock&      ketGtoPairsBlock,
                                             const int32_t              nKetContrPairs,
                                             const int32_t              iContrPair);


} // erikrrfunc namespace

