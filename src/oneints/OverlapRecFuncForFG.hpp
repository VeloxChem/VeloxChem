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

namespace ovlrecfunc { // ovlrecfunc namespace

    /**
    Computes batch of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CMemBlock2D<double>& pa2pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto);

    /**
    Computes sub-batch (0,5) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_0_5(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pbDistances,
                              const CMemBlock2D<double>& pa2pbDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto);

    /**
    Computes sub-batch (5,10) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_5_10(      CMemBlock2D<double>& primBuffer,
                               const CMemBlock2D<double>& auxBuffer,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pbDistances,
                               const CMemBlock2D<double>& pa2pbDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto);

    /**
    Computes sub-batch (10,15) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_10_15(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (15,20) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_15_20(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (20,25) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_20_25(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (25,30) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_25_30(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (30,35) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_30_35(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (35,40) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_35_40(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (40,45) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_40_45(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (45,50) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_45_50(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (50,55) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_50_55(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (55,60) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_55_60(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (60,65) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_60_65(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (65,70) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_65_70(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (70,75) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_70_75(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (75,80) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_75_80(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (80,85) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_80_85(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (85,90) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_85_90(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (90,95) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_90_95(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (95,100) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_95_100(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto);

    /**
    Computes sub-batch (100,105) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_100_105(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pbDistances,
                                  const CMemBlock2D<double>& pa2pbDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);

    /**
    Computes sub-batch (105,110) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_105_110(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pbDistances,
                                  const CMemBlock2D<double>& pa2pbDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);

    /**
    Computes sub-batch (110,115) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_110_115(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pbDistances,
                                  const CMemBlock2D<double>& pa2pbDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);

    /**
    Computes sub-batch (115,120) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_115_120(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pbDistances,
                                  const CMemBlock2D<double>& pa2pbDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);

    /**
    Computes sub-batch (120,125) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_120_125(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pbDistances,
                                  const CMemBlock2D<double>& pa2pbDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);

    /**
    Computes sub-batch (125,130) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_125_130(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pbDistances,
                                  const CMemBlock2D<double>& pa2pbDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);

    /**
    Computes sub-batch (130,135) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_130_135(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pbDistances,
                                  const CMemBlock2D<double>& pa2pbDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);

    /**
    Computes sub-batch (135,140) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_135_140(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pbDistances,
                                  const CMemBlock2D<double>& pa2pbDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);

    /**
    Computes sub-batch (140,145) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_140_145(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pbDistances,
                                  const CMemBlock2D<double>& pa2pbDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);

    /**
    Computes sub-batch (145,150) of primitive (F||G) overlap integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compOverlapForFG_145_150(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pbDistances,
                                  const CMemBlock2D<double>& pa2pbDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);


} // ovlrecfunc namespace

