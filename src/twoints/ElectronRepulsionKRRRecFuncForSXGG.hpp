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

namespace erikrrfunc { // erikrrfunc namespace

    /**
    Computes batch of contracted (SX|G|GG) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXGG(      CMemBlock2D<double>& ketBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& cdDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair);

    /**
    Computes sub-batch (0,75) of contracted (SX|G|GG) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXGG_0_75(      CMemBlock2D<double>& ketBuffer,
                                           const CRecursionMap&       recursionMap,
                                           const CMemBlock2D<double>& cdDistances,
                                           const CGtoPairsBlock&      braGtoPairsBlock,
                                           const CGtoPairsBlock&      ketGtoPairsBlock,
                                           const int32_t              nKetContrPairs,
                                           const int32_t              iContrPair);

    /**
    Computes sub-batch (75,150) of contracted (SX|G|GG) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXGG_75_150(      CMemBlock2D<double>& ketBuffer,
                                             const CRecursionMap&       recursionMap,
                                             const CMemBlock2D<double>& cdDistances,
                                             const CGtoPairsBlock&      braGtoPairsBlock,
                                             const CGtoPairsBlock&      ketGtoPairsBlock,
                                             const int32_t              nKetContrPairs,
                                             const int32_t              iContrPair);

    /**
    Computes sub-batch (150,225) of contracted (SX|G|GG) electron repulsion integrals and stores
    results in integrals buffer.

    @param ketBuffer the horizontal recursion buffer for ket side.
    @param recursionMap the recursion map object.
    @param cdDistances the vector of distances R(CD) = C - D.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForSXGG_150_225(      CMemBlock2D<double>& ketBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& cdDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetContrPairs,
                                              const int32_t              iContrPair);


} // erikrrfunc namespace

