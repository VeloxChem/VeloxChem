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

namespace eribrrfunc { // eribrrfunc namespace

    /**
    Computes batch of contracted (GG|G|XY electron repulsion integrals and stores
    results in integrals buffer.

    @param braBuffer the horizontal recursion buffer for bra side.
    @param recursionMap the recursion map object.
    @param abDistances the vector of distances R(AB) = A - B.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForGGXY(      CMemBlock2D<double>& braBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& abDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair);

    /**
    Computes sub-batch (0,75) of contracted (GG|G|XY electron repulsion integrals and stores
    results in integrals buffer.

    @param braBuffer the horizontal recursion buffer for bra side.
    @param recursionMap the recursion map object.
    @param abDistances the vector of distances R(AB) = A - B.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForGGXY_0_75(      CMemBlock2D<double>& braBuffer,
                                           const CRecursionMap&       recursionMap,
                                           const CMemBlock2D<double>& abDistances,
                                           const CGtoPairsBlock&      braGtoPairsBlock,
                                           const CGtoPairsBlock&      ketGtoPairsBlock,
                                           const int32_t              nKetContrPairs,
                                           const int32_t              iContrPair);

    /**
    Computes sub-batch (75,150) of contracted (GG|G|XY electron repulsion integrals and stores
    results in integrals buffer.

    @param braBuffer the horizontal recursion buffer for bra side.
    @param recursionMap the recursion map object.
    @param abDistances the vector of distances R(AB) = A - B.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForGGXY_75_150(      CMemBlock2D<double>& braBuffer,
                                             const CRecursionMap&       recursionMap,
                                             const CMemBlock2D<double>& abDistances,
                                             const CGtoPairsBlock&      braGtoPairsBlock,
                                             const CGtoPairsBlock&      ketGtoPairsBlock,
                                             const int32_t              nKetContrPairs,
                                             const int32_t              iContrPair);

    /**
    Computes sub-batch (150,225) of contracted (GG|G|XY electron repulsion integrals and stores
    results in integrals buffer.

    @param braBuffer the horizontal recursion buffer for bra side.
    @param recursionMap the recursion map object.
    @param abDistances the vector of distances R(AB) = A - B.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForGGXY_150_225(      CMemBlock2D<double>& braBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& abDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetContrPairs,
                                              const int32_t              iContrPair);


} // eribrrfunc namespace

