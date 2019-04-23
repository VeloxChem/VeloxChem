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

namespace npotrecfunc { // npotrecfunc namespace

    /**
    Computes batch of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

    /**
    Computes sub-batch (0,5) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_0_5(      CMemBlock2D<double>& primBuffer,
                                       const CMemBlock2D<double>& auxBuffer,
                                       const CMemBlock2D<double>& osFactors,
                                       const CMemBlock2D<double>& paDistances,
                                       const CMemBlock2D<double>& pbDistances,
                                       const CMemBlock2D<double>& pcDistances,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock,
                                       const int32_t              iContrGto);

    /**
    Computes sub-batch (5,10) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_5_10(      CMemBlock2D<double>& primBuffer,
                                        const CMemBlock2D<double>& auxBuffer,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& paDistances,
                                        const CMemBlock2D<double>& pbDistances,
                                        const CMemBlock2D<double>& pcDistances,
                                        const CGtoBlock&           braGtoBlock,
                                        const CGtoBlock&           ketGtoBlock,
                                        const int32_t              iContrGto);

    /**
    Computes sub-batch (10,15) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_10_15(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (15,20) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_15_20(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (20,25) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_20_25(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (25,30) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_25_30(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (30,35) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_30_35(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (35,40) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_35_40(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (40,45) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_40_45(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (45,50) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_45_50(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (50,55) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_50_55(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (55,60) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_55_60(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (60,65) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_60_65(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (65,70) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_65_70(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (70,75) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_70_75(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (75,80) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_75_80(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (80,85) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_80_85(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (85,90) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_85_90(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (90,95) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_90_95(      CMemBlock2D<double>& primBuffer,
                                         const CMemBlock2D<double>& auxBuffer,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& paDistances,
                                         const CMemBlock2D<double>& pbDistances,
                                         const CMemBlock2D<double>& pcDistances,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const int32_t              iContrGto);

    /**
    Computes sub-batch (95,100) of primitive (F|A|F) nuclear potential integrals and stores
    results in primitives buffer.

    @param primBuffer the primitives buffer.
    @param auxBuffer the auxilaries buffer.
    @param osFactors the Obara-Saika recursion factors.
    @param paDistances the vector of distances R(PA) = P - A.
    @param pbDistances the vector of distances R(PB) = P - B.
    @param pcDistances the vector of distances R(PC) = P - C.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compNuclearPotentialForFF_95_100(      CMemBlock2D<double>& primBuffer,
                                          const CMemBlock2D<double>& auxBuffer,
                                          const CMemBlock2D<double>& osFactors,
                                          const CMemBlock2D<double>& paDistances,
                                          const CMemBlock2D<double>& pbDistances,
                                          const CMemBlock2D<double>& pcDistances,
                                          const CGtoBlock&           braGtoBlock,
                                          const CGtoBlock&           ketGtoBlock,
                                          const int32_t              iContrGto);


} // npotrecfunc namespace

