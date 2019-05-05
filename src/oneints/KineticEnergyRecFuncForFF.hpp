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

namespace kinrecfunc { // kinrecfunc namespace

    /**
    Computes batch of primitive (F|T|F) kinetic energy integrals and stores
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
    void compKineticEnergyForFF(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (0,10) of primitive (F|T|F) kinetic energy integrals and stores
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
    void compKineticEnergyForFF_0_10(      CMemBlock2D<double>& primBuffer,
                                     const CMemBlock2D<double>& auxBuffer,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
                                     const CMemBlock2D<double>& pbDistances,
                                     const CMemBlock2D<double>& pa2pbDistances,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto);

    /**
    Computes sub-batch (10,20) of primitive (F|T|F) kinetic energy integrals and stores
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
    void compKineticEnergyForFF_10_20(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pbDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes sub-batch (20,30) of primitive (F|T|F) kinetic energy integrals and stores
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
    void compKineticEnergyForFF_20_30(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pbDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes sub-batch (30,40) of primitive (F|T|F) kinetic energy integrals and stores
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
    void compKineticEnergyForFF_30_40(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pbDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes sub-batch (40,50) of primitive (F|T|F) kinetic energy integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the set of distance tensors R(R(PA) = P - A.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compKineticEnergyForFF_40_50(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes sub-batch (50,60) of primitive (F|T|F) kinetic energy integrals and stores
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
    void compKineticEnergyForFF_50_60(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pbDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes sub-batch (60,70) of primitive (F|T|F) kinetic energy integrals and stores
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
    void compKineticEnergyForFF_60_70(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pbDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes sub-batch (70,80) of primitive (F|T|F) kinetic energy integrals and stores
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
    void compKineticEnergyForFF_70_80(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pbDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes sub-batch (80,90) of primitive (F|T|F) kinetic energy integrals and stores
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
    void compKineticEnergyForFF_80_90(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pbDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes sub-batch (90,100) of primitive (F|T|F) kinetic energy integrals and stores
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
    void compKineticEnergyForFF_90_100(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pa2pbDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);


} // kinrecfunc namespace

