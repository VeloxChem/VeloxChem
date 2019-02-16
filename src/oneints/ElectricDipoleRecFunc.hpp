//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef ElectricDipoleRecFunc_hpp
#define ElectricDipoleRecFunc_hpp

#include "MemBlock2D.hpp"
#include "GtoBlock.hpp"
#include "VecIndexes.hpp"

namespace ediprecfunc { // ediprecfunc namespace
    
    /**
     Computes batch of primitive (S|M|S) electric dipole integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param abDistances the vector of distances R(AB) = A - B.
     @param pcDistances the vector of distances R(PC) = P - C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectricDipoleForSS(      CMemBlock2D<double>&  primBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  osFactors,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CMemBlock2D<double>&  pcDistances,
                                 const CGtoBlock&            braGtoBlock,
                                 const CGtoBlock&            ketGtoBlock,
                                 const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|M|P) electric dipole integrals and stores
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
    void
    compElectricDipoleForSP(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  pbDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|M|S) electric dipole integrals and stores
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
    void
    compElectricDipoleForPS(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|M|P) electric dipole integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void
    compElectricDipoleForPP(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CMemBlock2D<double>&  pbDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|M|D) electric dipole integrals and stores
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
    void
    compElectricDipoleForSD(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  pbDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|M|S) electric dipole integrals and stores
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
    void
    compElectricDipoleForDS(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|M|D) electric dipole integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void
    compElectricDipoleForPD(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CMemBlock2D<double>&  pbDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|M|P) electric dipole integrals and stores
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
    void
    compElectricDipoleForDP(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|M|F) electric dipole integrals and stores
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
    void
    compElectricDipoleForSF(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  pbDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|M|S) electric dipole integrals and stores
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
    void
    compElectricDipoleForFS(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|M|D) electric dipole integrals and stores
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
    void
    compElectricDipoleForDD(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|M|F) electric dipole integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void
    compElectricDipoleForPF(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CMemBlock2D<double>&  pbDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|M|P) electric dipole integrals and stores
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
    void
    compElectricDipoleForFP(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|M|G) electric dipole integrals and stores
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
    void
    compElectricDipoleForSG(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  pbDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|M|S) electric dipole integrals and stores
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
    void
    compElectricDipoleForGS(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|M|F) electric dipole integrals and stores
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
    void
    compElectricDipoleForDF(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|M|D) electric dipole integrals and stores
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
    void
    compElectricDipoleForFD(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|M|F) electric dipole integrals and stores
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
    void
    compElectricDipoleForFF(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|M|G) electric dipole integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void
    compElectricDipoleForPG(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CMemBlock2D<double>&  pbDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|M|P) electric dipole integrals and stores
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
    void
    compElectricDipoleForGP(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|M|G) electric dipole integrals and stores
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
    void
    compElectricDipoleForDG(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|M|D) electric dipole integrals and stores
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
    void
    compElectricDipoleForGD(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|M|G) electric dipole integrals and stores
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
    void
    compElectricDipoleForFG(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|M|F) electric dipole integrals and stores
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
    void
    compElectricDipoleForGF(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|M|G) electric dipole integrals and stores
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
    void
    compElectricDipoleForGG(      CMemBlock2D<double>&  primBuffer,
                            const CVecThreeIndexes&     recPattern,
                            const std::vector<int32_t>& recIndexes,
                            const CMemBlock2D<double>&  osFactors,
                            const CMemBlock2D<double>&  paDistances,
                            const CGtoBlock&            braGtoBlock,
                            const CGtoBlock&            ketGtoBlock,
                            const int32_t               iContrGto);
    
} // ediprecfunc namespace

#endif /* ElectricDipoleRecFunc_hpp */
