//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef GenFunc_hpp
#define GenFunc_hpp

#include <cstdint>
#include <vector>

#include "MemBlock.hpp"
#include "MemBlock2D.hpp"
#include "SphericalMomentum.hpp"
#include "VecIndexes.hpp"
#include "GtoBlock.hpp"
#include "GtoContainer.hpp"
#include "SparseMatrix.hpp"
#include "GtoPairsBlock.hpp"


namespace genfunc { // genfunc namespace
    
    /**
     Contracts primitive data vectors to contracted data vectors.

     @param contrData the contracted data vectors.
     @param primData the primitive data vectors.
     @param contrIndex the index of first contracted data vector.
     @param primIndex the index of first primitive data vector.
     @param startPositions the vector of start positions in contraction pattern.
     @param endPositions the vector of end positions in contractrion pattern.
     @param nElements the number of elements in individual contracted data
            vector.
     @param nBlocks the number of contracted vectors.
     */
    void contract(      CMemBlock2D<double>& contrData,
                  const CMemBlock2D<double>& primData,
                  const int32_t              contrIndex,
                  const int32_t              primIndex,
                  const int32_t*             startPositions,
                  const int32_t*             endPositions,
                  const int32_t              nElements,
                  const int32_t              nBlocks);
    
    
    /**
     Contracts primitive data vectors to contracted data vectors using two step
     procedure. NOTE: Primitive data is destroyed during contraction process.

     @param contrData the contracted data vectors.
     @param primData the primitive data vectors.
     @param primIndex the index of first primitive data vector.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void contract(      CMemBlock2D<double>& contrData,
                        CMemBlock2D<double>& primData,
                  const int32_t              primIndex,
                  const CGtoBlock&           braGtoBlock,
                  const CGtoBlock&           ketGtoBlock,
                  const int32_t              iContrGto);
    
    /**
     Contracts set of primitive data vectoes to contracted data vectors using
     two step procedure. NOTE: Primitive data is destroyed during contraction
     process.

     @param contrData the contracted data vectors.
     @param primData the primitive data vectors.
     @param contrPattern the contracted data vectors distribution pattern.
     @param contrIndexes the contracted data vectors indexing pattern.
     @param primPattern the primitive data vectors distribution pattern.
     @param primIndexes the primitive data vectors indexing pattern,
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void contract(      CMemBlock2D<double>&  contrData,
                        CMemBlock2D<double>&  primData,
                  const CVecThreeIndexes&     contrPattern,
                  const std::vector<int32_t>& contrIndexes,
                  const CVecThreeIndexes&     primPattern,
                  const std::vector<int32_t>& primIndexes,
                  const CGtoBlock&            braGtoBlock,
                  const CGtoPairsBlock&       ketGtoPairsBlock,
                  const int32_t               iContrGto);
    
    /**
     Transforms Cartesian data vectors to spherical data vectors.

     @param spherData the spherical data vectors.
     @param cartData the Cartesian data vectors.
     @param spherMomentum the spherical momentum object.
     @param spherIndex the position of spherical data vectors.
     @param cartIndex the position of Cartesian data vectors.
     @param nElements the number elements in individual data vector.
     @param nBlocks the number of data vectors per spherical momentum component.
     */
    void transform(      CMemBlock2D<double>& spherData,
                   const CMemBlock2D<double>& cartData,
                   const CSphericalMomentum&  spherMomentum,
                   const int32_t              spherIndex,
                   const int32_t              cartIndex,
                   const int32_t              nElements,
                   const int32_t              nBlocks);
    
    /**
     Transforms Cartesian integrals to spherical integrals.
     
     @param spherData the spherical data vectors.
     @param cartData the Cartesian data vectors.
     @param braMomentum the spherical momentum object for bra side.
     @param ketMomentum the spherical momentum object for ket side.
     @param spherIndex the position of spherical data vectors.
     @param cartIndex the position of Cartesian data vectors.
     @param nElements the number elements in individual data vector.
     */
    void transform(      CMemBlock2D<double>& spherData,
                   const CMemBlock2D<double>& cartData,
                   const CSphericalMomentum&  braMomentum,
                   const CSphericalMomentum&  ketMomentum,
                   const int32_t              spherIndex,
                   const int32_t              cartIndex,
                   const int32_t              nElements);
    
    /**
     Transforms set Cartesian data vectors to set spherical data vectors by
     applying bra side transformation: <cart|f(x)|cart> to <spher|f(x)|cart>.

     @param spherData the spherical data vectors.
     @param cartData the Cartesian data
     @param spherMomentum the spherical momentum object for bra side.
     @param spherPattern the distribution pattern of spherical data vectors.
     @param spherIndexes the indexing vector of spherical data vectors.
     @param cartPattern the distribution pattern of Cartesian data vectors.
     @param cartIndexes the indexing vector of Cartesian data vectors.
     @param nElements the number of elements in data vectors.
     */
    void transform(      CMemBlock2D<double>&  spherData,
                   const CMemBlock2D<double>&  cartData,
                   const CSphericalMomentum&   spherMomentum,
                   const CVecThreeIndexes&     spherPattern,
                   const std::vector<int32_t>& spherIndexes,
                   const CVecThreeIndexes&     cartPattern,
                   const std::vector<int32_t>& cartIndexes,
                   const int32_t               nElements);
    
    /**
     Transforms Cartesian half-transformed integrals to spherical integrals by
     applying ket side transformation: <cart|f(x)|spher> to <spher|f(x)|spher>.
     
     @param spherData the spherical data vectors.
     @param cartData the Cartesian data vectors.
     @param ketMomentumC the spherical momentum object for C center on ket side.
     @param ketMomentumD the spherical momentum object for D center on ket side.
     @param spherIndex the position of spherical data vectors.
     @param cartIndex the position of Cartesian data vectors.
     @param nElements the number elements in individual data vector.
     @param nBlocks the number of data vectors on bra side.
     */
    void transform(      CMemBlock2D<double>& spherData,
                   const CMemBlock2D<double>& cartData,
                   const CSphericalMomentum&  ketMomentumC,
                   const CSphericalMomentum&  ketMomentumD,
                   const int32_t              spherIndex,
                   const int32_t              cartIndex,
                   const int32_t              nElements,
                   const int32_t              nBlocks);
    
    /**
     Reshapes batch of integrals batch (1 contracted GTO on bra side x all
     contracted GTOs on ket side) into compressed data and stores it in sparse
     matrix.

     @param sparseMatrix the sparse matrix.
     @param rowValues the temporary buffer for values of integrals.
     @param colIndexes the temporary buffer for indexes of integrals.
     @param spherData the batch of integrals.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compress(      CSparseMatrix&       sparseMatrix,
                        CMemBlock<double>&   rowValues,
                        CMemBlock<int32_t>&  colIndexes,
                  const CMemBlock2D<double>& spherData,
                  const CGtoBlock&           braGtoBlock,
                  const CGtoBlock&           ketGtoBlock,
                  const int32_t              iContrGto);
    
    /**
     Creates sparse matrix by distributing list of sparse matrices according to
     pattern given by GTOs containers from bra and ket sides.

     @param listOfMatrices the list of sparse matrices.
     @param braGtoContainer the GTOs container for bra side.
     @param ketGtoContainer the GTOs container for ket side.
     @return the sparse matrix.
     */
    CSparseMatrix distribute(const CSparseMatrix* listOfMatrices,
                             const CGtoContainer* braGtoContainer,
                             const CGtoContainer* ketGtoContainer);
    
    /**
     Distributes batch of integrals (contracted GTO on bra side x all contracted
     GTOs on ket side) in integrals matrix.

     @param values the integrals matrix.
     @param spherData the batch of integrals.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketAngularMomentum the angular momentum of ket side.
     @param nBraContrGtos the number of contracted GTOs on bra side.
     @param nKetContrGtos the number of contracted GTOs on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void distribute(      double*              values,
                    const CMemBlock2D<double>& spherData,
                    const int32_t              braAngularMomentum,
                    const int32_t              ketAngularMomentum,
                    const int32_t              nBraContrGtos,
                    const int32_t              nKetContrGtos,
                    const int32_t              iContrGto);
    
    /**
     Distributes batch of integrals (contracted GTO on bra side x all contracted
     GTOs on ket side) in integrals matrix.

     @param values the integrals matrix.
     @param spherData the batch of integrals.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param isBraEqualKet the flag indicating equality of bra and ket sides.
     @param nColumns the number of columns in integrals matrix.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void distribute(      double*              values,
                    const CMemBlock2D<double>& spherData,
                    const CGtoBlock&           braGtoBlock,
                    const CGtoBlock&           ketGtoBlock,
                    const bool                 isBraEqualKet,
                    const int32_t              nColumns,
                    const int32_t              iContrGto);
    
    /**
     Checks if two indexes object is inside vector of two indexes objects.

     @param vector the vector of two indexes objects.
     @param pair the two indexes object.
     @return true if two indexes object is found in vector of two indexes
             objects, false - otherwise.
     */
    bool isInVector(const CVecTwoIndexes& vector,
                    const CTwoIndexes&    pair);
    
    /**
     Checks if three indexes object is inside vector of three indexes objects.
     
     @param vector the vector of three indexes objects.
     @param triple the three indexes object.
     @return true if three indexes object is found in vector of three indexes
             objects, false - otherwise.
     */
    bool isInVector(const CVecThreeIndexes& vector,
                    const CThreeIndexes&    triple);
    
    /**
     Checks if four indexes object is inside vector of four indexes objects.
     
     @param vector the vector of four indexes objects.
     @param quadruple the four indexes object.
     @return true if four indexes object is found in vector of four indexes
     objects, false - otherwise.
     */
    bool isInVector(const CVecFourIndexes& vector,
                    const CFourIndexes&    quadruple);
    
    /**
     Conditionally adds two indexes object to vector of two indexes objects.
     Addition is skipped if two indexes object is not valid indexing pair or
     already resides in vector of two indexes objects.

     @param vector the vector of two indexes objects.
     @param pair the two indexes object.
     @return true if two indexes object is added to vector of two indexes
             objects, false otherwise.
     */
    bool addValidAndUniquePair(      CVecTwoIndexes& vector,
                               const CTwoIndexes&    pair);
    
    /**
    Conditionally adds three indexes object to vector of three indexes objects.
    Addition is skipped if three indexes object is not valid indexing triple or
    already resides in vector of three indexes objects.
        
    @param vector the vector of three indexes objects.
    @param triple the three indexes object.
    @return true if three indexes object is added to vector of three indexes
            objects, false otherwise.
    */
    bool addValidAndUniqueTriple(      CVecThreeIndexes& vector,
                                 const CThreeIndexes&    triple);
    
    /**
     Conditionally adds four indexes object to vector of four indexes objects.
     Addition is skipped if four indexes object is not valid indexing quadruple
     or already resides in vector of four indexes objects.
     
     @param vector the vector of three indexes objects.
     @param quadruple the four indexes object.
     @return true if four indexes object is added to vector of four indexes
     objects, false otherwise.
     */
    bool addValidAndUniqueQuadruple(      CVecFourIndexes& vector,
                                    const CFourIndexes&    quadruple);
    
    /**
     Finds index from vector of indexes associated with two indexes object in
     vector of two indexes objects.

     @param indexes the vector of indexes.
     @param vector the vector of two indexes objects.
     @param pair the two indexes object.
     @return the index assocated with two indexes object.
     */
    int32_t findPairIndex(const std::vector<int32_t>& indexes,
                          const CVecTwoIndexes&       vector,
                          const CTwoIndexes&          pair);
    
    /**
     Finds index from vector of indexes associated with three indexes object in
     vector of three indexes objects.
     
     @param indexes the vector of indexes.
     @param vector the vector of three indexes objects.
     @param triple the three indexes object.
     @return the index assocated with three indexes object.
     */
    int32_t findTripleIndex(const std::vector<int32_t>& indexes,
                            const CVecThreeIndexes&     vector,
                            const CThreeIndexes&        triple);
    
    /**
     Gets maximum order of indexes pair in vector of triple indexes objects,
     where third index in triple object is used to store order information.

     @param vector the vector of three indexes objects.
     @param firstIndex the first index of indexes pair.
     @param secondIndex the second index of indexes pair.
     @return the maximum order.
     */
    int32_t maxOrderOfPair(const CVecThreeIndexes& vector,
                           const int32_t           firstIndex,
                           const int32_t           secondIndex);
    
    /**
     Creates vector of three indexes objects of (x,y,0) type coresponding to
     all three index objects of (x,0,y) type found in vector of three indexes
     objects.
     
     @param vector the vector of three indexes objects.
     @return the vector of three indexes objects of (x,y,0) type.
     */
    CVecThreeIndexes getPairsFromTripleIndexes(const CVecThreeIndexes& vector);
    
    /**
     Creates vector of three indexes objects of (x,y,z) type coresponding to
     all four index objects of (0,x,y,z) type found in vector of four indexes
     objects.
     
     @param vector the vector of four indexes objects.
     @return the vector of three indexes objects of (x,y,z) type.
     */
    CVecThreeIndexes getTriplesFromQuadrupleIndexes(const CVecFourIndexes& vector);

} // genfunc namespace

#endif /* GenFunc_hpp */
