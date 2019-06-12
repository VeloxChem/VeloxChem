//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include<cstdint>

#include "MemBlock2D.hpp"
#include "GtoBlock.hpp"
#include "RecursionMap.hpp"

namespace lmomrecfunc { // lmomrecfunc namespace

    /**
    Computes batch of primitive (G|P|G) linear momentum integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compLinearMomentumForGG(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto);

    /**
    Computes sub-batch (0,49) of primitive (G|P|G) linear momentum integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compLinearMomentumForGG_0_49(      CMemBlock2D<double>& primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes sub-batch (49,98) of primitive (G|P|G) linear momentum integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compLinearMomentumForGG_49_98(      CMemBlock2D<double>& primBuffer,
                                       const CRecursionMap&       recursionMap,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (98,147) of primitive (G|P|G) linear momentum integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compLinearMomentumForGG_98_147(      CMemBlock2D<double>& primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& paDistances,
                                        const CGtoBlock&           braGtoBlock,
                                        const CGtoBlock&           ketGtoBlock,
                                        const int32_t              iContrGto);

    /**
    Computes sub-batch (147,195) of primitive (G|P|G) linear momentum integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compLinearMomentumForGG_147_195(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (195,243) of primitive (G|P|G) linear momentum integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compLinearMomentumForGG_195_243(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (243,291) of primitive (G|P|G) linear momentum integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compLinearMomentumForGG_243_291(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (291,339) of primitive (G|P|G) linear momentum integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compLinearMomentumForGG_291_339(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (339,387) of primitive (G|P|G) linear momentum integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compLinearMomentumForGG_339_387(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (387,435) of primitive (G|P|G) linear momentum integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compLinearMomentumForGG_387_435(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (435,483) of primitive (G|P|G) linear momentum integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compLinearMomentumForGG_435_483(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (483,531) of primitive (G|P|G) linear momentum integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compLinearMomentumForGG_483_531(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (531,579) of primitive (G|P|G) linear momentum integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compLinearMomentumForGG_531_579(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (579,627) of primitive (G|P|G) linear momentum integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compLinearMomentumForGG_579_627(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (627,675) of primitive (G|P|G) linear momentum integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param recursionMap the recursion map object.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compLinearMomentumForGG_627_675(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);


} // lmomrecfunc namespace

