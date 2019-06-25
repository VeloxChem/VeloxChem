//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef ElectronRepulsionIntegralsDriver_hpp
#define ElectronRepulsionIntegralsDriver_hpp

#include <cstdint>

#include "mpi.h"

#ifdef MAC_OS_OMP
#include "/opt/intel/compilers_and_libraries/mac/include/omp.h"
#else
#include "omp.h"
#endif

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "GtoPairsBlock.hpp"
#include "VecIndexes.hpp"
#include "BoysFunction.hpp"
#include "GtoPairsContainer.hpp"
#include "TwoIntsDistributor.hpp"
#include "ScreeningContainer.hpp"
#include "VecMemBlocks.hpp"
#include "AODensityMatrix.hpp"
#include "AOFockMatrix.hpp"
#include "RecursionMap.hpp"

/**
 Class CElectronicRepulsionIntegralsDriver computes electron repulsion
 <f(r)g(r)| 1/|r-r'||h(r') i(r')> integrals.
 
 @author Z. Rinkevicius
 */
class CElectronRepulsionIntegralsDriver
{
    /**
     The rank of associated local MPI process.
     */
    int32_t _locRank;
    
    /**
     The total number of local MPI processes.
     */
    int32_t _locNodes;
    
    /**
     The MPI communicator.
     */
    MPI_Comm _locComm;
    
    /**
     Computes electronic repulsion integrals for combination of GTOs pairs
     blocks.
     
     @param distPattern the ponter to integrals distribution pattern.
     @param intsScreener the integrals screener object.
     @param braGtoPairsBlock the GTOs pairsblock on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void _compElectronRepulsionForGtoPairsBlocks(      CTwoIntsDistribution&   distPattern,
                                                 const CCauchySchwarzScreener& intsScreener,
                                                 const CGtoPairsBlock&         braGtoPairsBlock,
                                                 const CGtoPairsBlock&         ketGtoPairsBlock) const;
    
    /**
     Sets horrizontal recursion map for bra side.
     
     @param angularMomentumA the angular momentum of center A.
     @param angularMomentumB the angular momentum of center B.
     @param angularMomentumC the angular momentum of center C.
     @param angularMomentumD the angular momentum od center D-
     @return the recursion map.
     */
    CRecursionMap _setBraHorizontalRecursionPattern(const int32_t angularMomentumA,
                                                    const int32_t angularMomentumB,
                                                    const int32_t angularMomentumC,
                                                    const int32_t angularMomentumD) const;
    
    /**
     Sets horrizontal recursion map for ket side.
     
     @param angularMomentumA the angular momentum of center A.
     @param angularMomentumB the angular momentum of center B.
     @param angularMomentumC the angular momentum of center C.
     @param angularMomentumD the angular momentum od center D.
     @return the recursion map.
     */
    CRecursionMap _setKetHorizontalRecursionPattern(const int32_t angularMomentumA,
                                                    const int32_t angularMomentumB,
                                                    const int32_t angularMomentumC,
                                                    const int32_t angularMomentumD) const;
    
    /**
     Sets Obara-Saika vertical recursion pattern for given set of leading terms.
     
     @param leadTerms the recursion map from horizontal recursion on bra and ket sides.
     @param maxNumberOfPrimPairs the maximum number of primitive pairs.
     @return the vector of three indexes object with recursion pattern.
     */
    CRecursionMap _setVerticalRecursionMap(const CRecursionMap& leadTerms,
                                           const int32_t        maxNumberOfPrimPairs) const;
    
    /**
     Computes batch of primitive electron repulsion integrals using Obara-Saika
     recursion and stores results in primitives buffer.
     Reference: S. Obara, A. Saika, J. Chem. Phys. 84, 3963 (1986).
     
     Batch size: (one contracted GTOs pair on bra side) x (all contracted GTOs
                  pairs on ket side).
     
     @param primBuffer the pointer to array of primitive buffers.
     @param recursionMap the vertical recursion map for bra and ket sides.
     @param bfTable the Boys function evaluator.
     @param bfArguments the vector of Boys function arguments.
     @param bfValues the vector of Boys function values.
     @param bfOrder the order of Boys function.
     @param osFactors the Obara-Saika recursion factors.
     @param pqDistances the vector of distances R(PQ) = P - Q.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param wqDistances the vector of distances R(WQ) = W - Q.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void _compPrimElectronRepulsionInts(      CMemBlock2D<double>*  primBuffer,
                                        const CRecursionMap&        recursionMap, 
                                        const CBoysFunction&        bfTable,
                                              CMemBlock<double>&    bfArguments,
                                              CMemBlock2D<double>&  bfValues,
                                        const int32_t               bfOrder,
                                        const CMemBlock2D<double>&  osFactors,
                                        const CMemBlock2D<double>&  pqDistances,
                                        const CMemBlock2D<double>&  wpDistances,
                                        const CMemBlock2D<double>&  wqDistances,
                                        const CGtoPairsBlock&       braGtoPairsBlock,
                                        const CGtoPairsBlock&       ketGtoPairsBlock,
                                        const int32_t               nKetPrimPairs,
                                        const int32_t               iContrPair) const;
    
    /**
     Applies horizontal recursion on ket side of contracted integrals buffer.
     
     @param ketBuffer the horizontal recursion buffer.
     @param recursionMap the horizontal recursion map for ket side.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void _applyHRRonKet(      CMemBlock2D<double>&  ketBuffer,
                        const CRecursionMap&        recursionMap,
                        const CMemBlock2D<double>&  cdDistances,
                        const CGtoPairsBlock&       braGtoPairsBlock,
                        const CGtoPairsBlock&       ketGtoPairsBlock,
                        const int32_t               nKetContrPairs,
                        const int32_t               iContrPair) const;
    
    /**
     Applies horizontal recursion on bra side of contracted integrals buffer.
     
     @param braBuffer the horizontal recursion buffer.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param recursionMap the horizontal recursion map for ket side.
     @param abDistances the vector of distances R(AB) = A - B.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetContrPairs the number of contractes GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void _applyHRRonBra(      CMemBlock2D<double>&  braBuffer,
                        const CRecursionMap&        recursionMap,
                        const CMemBlock2D<double>&  abDistances,
                        const CGtoPairsBlock&       braGtoPairsBlock,
                        const CGtoPairsBlock&       ketGtoPairsBlock,
                        const int32_t               nKetContrPairs,
                        const int32_t               iContrPair) const;
    
    /**
     Comutes electron repulsion integrals and stores them into AO Fock matrix
     for two GTOs pairs containers.
     
     @param aoFockMatrix the AO Fock matrix.
     @param aoDensityMatrix the AO density matrix.
     @param braGtoPairsContainer the GTOs pairs container on bra side.
     @param ketGtoPairsContainer the GTOs pairs container on ket side.
     @param screeningContainer the screening container object. 
     */
    void _compElectronRepulsionIntegrals(      CAOFockMatrix&       aoFockMatrix,
                                         const CAODensityMatrix&    aoDensityMatrix,
                                         const CGtoPairsContainer*  braGtoPairsContainer,
                                         const CGtoPairsContainer*  ketGtoPairsContainer,
                                         const CScreeningContainer* screeningContainer) const;

    /**
     Comutes maximum Q values of electron repulsion integrals for GTOs pairs
     block.

     @param qValuesBuffer the Q values buffer.
     @param gtoPairsBlock the GTOs pairs block.
     */
    void _compMaxQValuesForGtoPairsBlock(      double*         qValuesBuffer,
                                         const CGtoPairsBlock& gtoPairsBlock) const;
    
    /**
     Allocates Q values buffer for GTOs pair blocks container.

     @param gtoPairsContainer the GTOs pair blocks container.
     @return the Q values buffer as vector of mememory block objects.
     */
     CVecMemBlock<double> _getQValuesBuffer(const CGtoPairsContainer& gtoPairsContainer) const;
    
    /**
     Creates tasks execution grid for computation of AO Fock matrix on
     specific MPI process.

     @param nBraGtoPairsBlocks the number of GTOs pairs blocks on bra side.
     @param nKetGtoPairsBlocks the number of GTOs pairs blocks on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
             blocks.
     @return the vector of tasks with elements 1 or 0: 1  execute task on
             this MPI process, 0 skip task on this MPI process.
     */
    CMemBlock<int32_t> _setTasksGrid(const int32_t nBraGtoPairsBlocks,
                                     const int32_t nKetGtoPairsBlocks,
                                     const bool    isBraEqualKet) const;
        
public:
    
    /**
     Creates an electron repulsion integrals driver object using MPI info.
     
     @param comm the MPI communicator.
     */
    CElectronRepulsionIntegralsDriver(MPI_Comm comm);
    
    /**
     Destroys an electron repulsion integrals driver object.
     */
    ~CElectronRepulsionIntegralsDriver();
    
    /**
     Computes electron repulsion integrals and stores them in AO Fock matrix
     for molecule with specific AO basis set.
     
     @param aoFockMatrix the AO Fock matrix.
     @param aoDensityMatrix the AO density matrix.
     @param molecule the molecule.
     @param aoBasis the molecular AO basis.
     @param screeningContainer the screening container object.
     */
    void compute(      CAOFockMatrix&       aoFockMatrix,
                 const CAODensityMatrix&    aoDensityMatrix,
                 const CMolecule&           molecule,
                 const CMolecularBasis&     aoBasis,
                 const CScreeningContainer& screeningContainer) const;
    
    /**
     Computes Q values for electron repulsion integrals for molecule with
     specific AO basis set and stores results in screening container object.
     
     @param screeningScheme the screening scheme for screening container object.
     @param threshold the screening threshold for screening container object.
     @param molecule the molecule.
     @param aoBasis the molecular AO basis.
     @return the screening container with Q values.
     */
    CScreeningContainer compute(const ericut           screeningScheme,
                                const double           threshold,
                                const CMolecule&       molecule,
                                const CMolecularBasis& aoBasis) const;
    
    /**
     Computes electron repulsion integrals blocks for pair of GTOs pairs blocks
     and stores them into integrals batch.
     
     @param intsBatch the pointer to integrals batch buffer.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compute(      double*         intsBatch,
                 const CGtoPairsBlock& braGtoPairsBlock,
                 const CGtoPairsBlock& ketGtoPairsBlock) const;
    
    /**
     Comutes maximum Q values of electron repulsion integrals for two GTOs
     pairs containers and store results in bra and ket Q values buffers.
     
     @param braQValuesBuffer the Q values buffer on bra side.
     @param ketQValuesBuffer the Q values buffer on ket side.
     @param braGtoPairsContainer the GTOs pairs container on bra side.
     @param ketGtoPairsContainer the GTOs pairs container on ket side.
     */
    void computeMaxQValues(      CVecMemBlock<double>* braQValuesBuffer,
                                 CVecMemBlock<double>* ketQValuesBuffer,
                           const CGtoPairsContainer*   braGtoPairsContainer,
                           const CGtoPairsContainer*   ketGtoPairsContainer) const;
};


#endif /* ElectronRepulsionIntegralsDriver_hpp */
