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
    Computes batch of primitive (P|T|P) kinetic energy integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compKineticEnergyForPP(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes batch of primitive (P|T|D) kinetic energy integrals and stores
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
    void compKineticEnergyForPD(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (0,9) of primitive (P|T|D) kinetic energy integrals and stores
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
    void compKineticEnergyForPD_0_9(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pa2pbDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

    /**
    Computes sub-batch (9,18) of primitive (P|T|D) kinetic energy integrals and stores
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
    void compKineticEnergyForPD_9_18(      CMemBlock2D<double>& primBuffer,
                                     const CMemBlock2D<double>& auxBuffer,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
                                     const CMemBlock2D<double>& pbDistances,
                                     const CMemBlock2D<double>& pa2pbDistances,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto);

    /**
    Computes batch of primitive (D|T|P) kinetic energy integrals and stores
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
    void compKineticEnergyForDP(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (0,9) of primitive (D|T|P) kinetic energy integrals and stores
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
    void compKineticEnergyForDP_0_9(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pa2pbDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

    /**
    Computes sub-batch (9,18) of primitive (D|T|P) kinetic energy integrals and stores
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
    void compKineticEnergyForDP_9_18(      CMemBlock2D<double>& primBuffer,
                                     const CMemBlock2D<double>& auxBuffer,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
                                     const CMemBlock2D<double>& pbDistances,
                                     const CMemBlock2D<double>& pa2pbDistances,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto);

    /**
    Computes batch of primitive (P|T|F) kinetic energy integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compKineticEnergyForPF(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (0,15) of primitive (P|T|F) kinetic energy integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compKineticEnergyForPF_0_15(      CMemBlock2D<double>& primBuffer,
                                     const CMemBlock2D<double>& auxBuffer,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& pbDistances,
                                     const CMemBlock2D<double>& pa2pbDistances,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto);

    /**
    Computes sub-batch (15,30) of primitive (P|T|F) kinetic energy integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param pbDistances the set of distance tensors R(R(PB) = P - B.
    @param pa2pbDistances the set of products of distance tensors R(PA)xR(PB).
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compKineticEnergyForPF_15_30(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& pbDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes batch of primitive (F|T|P) kinetic energy integrals and stores
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
    void compKineticEnergyForFP(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (0,15) of primitive (F|T|P) kinetic energy integrals and stores
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
    void compKineticEnergyForFP_0_15(      CMemBlock2D<double>& primBuffer,
                                     const CMemBlock2D<double>& auxBuffer,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
                                     const CMemBlock2D<double>& pa2pbDistances,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto);

    /**
    Computes sub-batch (15,30) of primitive (F|T|P) kinetic energy integrals and stores
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
    void compKineticEnergyForFP_15_30(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes batch of primitive (P|T|G) kinetic energy integrals and stores
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
    void compKineticEnergyForPG(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (0,9) of primitive (P|T|G) kinetic energy integrals and stores
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
    void compKineticEnergyForPG_0_9(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pa2pbDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

    /**
    Computes sub-batch (9,18) of primitive (P|T|G) kinetic energy integrals and stores
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
    void compKineticEnergyForPG_9_18(      CMemBlock2D<double>& primBuffer,
                                     const CMemBlock2D<double>& auxBuffer,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
                                     const CMemBlock2D<double>& pbDistances,
                                     const CMemBlock2D<double>& pa2pbDistances,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto);

    /**
    Computes sub-batch (18,27) of primitive (P|T|G) kinetic energy integrals and stores
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
    void compKineticEnergyForPG_18_27(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pbDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes sub-batch (27,36) of primitive (P|T|G) kinetic energy integrals and stores
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
    void compKineticEnergyForPG_27_36(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pbDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes sub-batch (36,45) of primitive (P|T|G) kinetic energy integrals and stores
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
    void compKineticEnergyForPG_36_45(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pbDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes batch of primitive (G|T|P) kinetic energy integrals and stores
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
    void compKineticEnergyForGP(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto);

    /**
    Computes sub-batch (0,9) of primitive (G|T|P) kinetic energy integrals and stores
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
    void compKineticEnergyForGP_0_9(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pa2pbDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

    /**
    Computes sub-batch (9,18) of primitive (G|T|P) kinetic energy integrals and stores
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
    void compKineticEnergyForGP_9_18(      CMemBlock2D<double>& primBuffer,
                                     const CMemBlock2D<double>& auxBuffer,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
                                     const CMemBlock2D<double>& pbDistances,
                                     const CMemBlock2D<double>& pa2pbDistances,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto);

    /**
    Computes sub-batch (18,27) of primitive (G|T|P) kinetic energy integrals and stores
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
    void compKineticEnergyForGP_18_27(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes sub-batch (27,36) of primitive (G|T|P) kinetic energy integrals and stores
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
    void compKineticEnergyForGP_27_36(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pbDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

    /**
    Computes sub-batch (36,45) of primitive (G|T|P) kinetic energy integrals and stores
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
    void compKineticEnergyForGP_36_45(      CMemBlock2D<double>& primBuffer,
                                      const CMemBlock2D<double>& auxBuffer,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pbDistances,
                                      const CMemBlock2D<double>& pa2pbDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);


} // kinrecfunc namespace

