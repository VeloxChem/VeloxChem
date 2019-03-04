//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef KineticEnergyRecFunc_hpp
#define KineticEnergyRecFunc_hpp

#include "MemBlock2D.hpp"
#include "GtoBlock.hpp"
#include "VecIndexes.hpp"

namespace kinrecfunc { // kinrecfunc namespace
    
    /**
     Computes batch of primitive (S|T|S) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param abDistances the vector of distances R(AB) = A - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForSS(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  abDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|T|P) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForSP(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  pbDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|T|S) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForPS(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|T|P) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForPP(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|T|D) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForSD(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  pbDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|T|S) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForDS(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|T|D) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForPD(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|T|P) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForDP(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|T|D) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForDD(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|T|F) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForSF(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  pbDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|T|S) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForFS(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|T|F) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForPF(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|T|P) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForFP(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|T|F) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForDF(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|T|D) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForFD(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|T|F) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForFF(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|T|G) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForSG(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  pbDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|T|S) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForGS(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|T|G) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForPG(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|T|P) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForGP(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|T|G) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForDG(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|T|D) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForGD(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|T|G) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForFG(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|T|F) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForGF(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|T|G) kinetic energy integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compKineticEnergyForGG(      CMemBlock2D<double>&  primBuffer,
                                const CVecThreeIndexes&     recPattern,
                                const std::vector<int32_t>& recIndexes,
                                const CMemBlock2D<double>&  osFactors,
                                const CMemBlock2D<double>&  paDistances,
                                const CGtoBlock&            braGtoBlock,
                                const CGtoBlock&            ketGtoBlock,
                                const int32_t               iContrGto);
    
} // kinrecfunc namespace

#endif /* KineticEnergyRecFunc_hpp */
