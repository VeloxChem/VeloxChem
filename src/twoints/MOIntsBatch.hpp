//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MOIntsBatch_hpp
#define MOIntsBatch_hpp

#include <cstdint>
#include <vector>

#include "DenseMatrix.hpp"
#include "TwoIndexes.hpp"
#include "MOIntsType.hpp"
#include "AOFockMatrix.hpp"

/**
 Class CMOIntsBatch stores batch of MO integrals and provides set of methods
 for handling of MO integrals data.
 
 @author Z. Rinkevicius
 */
class CMOIntsBatch
{
    /**
     The set of MO integrals for specific generator pairs.
     */
    std::vector<CDenseMatrix> _moIntegrals;
    
    /**
     The set of generator pairs.
     */
    std::vector<CTwoIndexes> _generatorPairs;
    
    /**
     The external indexes of MO integrals.
     */
    CTwoIndexes _externalIndexes;
    
    /**
     The type of integrals batch.
     */
    
    moints _batchType;
    
public:
    
    /**
     Creates an empty MO integrals batch object.
     */
    CMOIntsBatch();
    
    /**
     Creates a MO integrals batch object.
     
     @param moIntegrals the set of MO integrals <ij|xy> (x,y - external indexes).
     @param generatorPairs the set of generator pairs (i,j).
     @param externalIndexes the positions of external indexes in <ij|xy>
             integrals.
     @param batchType the type of MO integrals batch.
     */
    CMOIntsBatch(const std::vector<CDenseMatrix>& moIntegrals,
                 const std::vector<CTwoIndexes>&  generatorPairs,
                 const CTwoIndexes&               externalIndexes,
                 const moints                     batchType);
    
    /**
     Creates a MO integrals batch object by copying other MO integrals batch
     object.
     
     @param source the MO integrals batch object.
     */
    CMOIntsBatch(const CMOIntsBatch& source);
    
    /**
     Creates a MO integrals batch object by moving other MO integrals batch
     object.
     
     @param source the MO integrals batch object.
     */
    CMOIntsBatch(CMOIntsBatch&& source) noexcept;
    
    /**
     Destroys a MO integrals batch object.
     */
    ~CMOIntsBatch();
    
    /**
     Assigns a MO integrals batch object by copying other MO integrals batch '
     object.
     
     @param source the MO integrals batch object.
     */
    CMOIntsBatch& operator=(const CMOIntsBatch& source);
    
    /**
     Assigns a MO integrals batch object by moving other MO integrals batch
     object.
     
     @param source the MO integrals batch object.
     */
    CMOIntsBatch& operator=(CMOIntsBatch&& source) noexcept;
    
    /**
     Compares MO integrals batch object with other MO integrals batch object.
     
     @param other the MO integrals batch object.
     @return true if MO integrals batch objects are equal, false otherwise.
     */
    bool operator==(const CMOIntsBatch& other) const;
    
    /**
     Compares MO integrals batch object with other MO integrals batch object.
     
     @param other the MO integrals batch object.
     @return true if MO integrals batch objects are not equal, false otherwise.
     */
    bool operator!=(const CMOIntsBatch& other) const;
    
    /**
     Adds set of mo integrals to batch by transforming AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param braVector the transformation matrix for bra side.
     @param ketVector the transformation matrix for ket side.
     @param braIndexes the indexes of C_i orbitals.
     @param ketIndexes the indexes of C_j orbitals.
     */
    void append(const CAOFockMatrix&        aoFockMatrix,
                const CDenseMatrix&         braVector,
                const CDenseMatrix&         ketVector,
                const std::vector<int32_t>& braIndexes,
                const std::vector<int32_t>& ketIndexes);
    
    /**
     Sets type of MO integrals batch.
     
     @param batchType the type of MO integrals batch.
     */
    void setBatchType(const moints batchType);
    
    /**
     Sets posiitons of external indexes in in <ij|xy> integrals.
     
     @param externalIndexes the positions of external indexes.
     */
    void setExternalIndexes(const CTwoIndexes& externalIndexes);
    
    /**
     Converts MO integrals batch object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the MO integrals batch object.
     */
    friend std::ostream& operator<<(      std::ostream& output,
                                    const CMOIntsBatch& source);
};


#endif /* MOIntsBatch_hpp */
