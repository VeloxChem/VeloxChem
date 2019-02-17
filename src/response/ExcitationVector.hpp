//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef ExcitationVector_hpp
#define ExcitationVector_hpp

#include <cstdint>
#include <ostream>
#include <vector>

#include "SpinBlock.hpp"
#include "MemBlock.hpp"

/**
Class CExcitationVector class stores information on one particle excitation and
provides set of methods for manipulating their data.

@author Z. Rinkevicius
*/
class CExcitationVector
{
    /**
     The type of generating excitations.
     */
    szblock _excitationType;
    
    /**
     The vector of molecular orbital indexes for bra side.
     */
    CMemBlock<int32_t> _braIndexes;
    
    /**
     The vector of molecular orbital indexes for ket side.
     */
    CMemBlock<int32_t> _ketIndexes;
    
    /**
     The vector of coefficients associated with one particle excitations.
     */
    CMemBlock<double> _coefficents;
    
public:
    
    /**
     Creates an empty excitation vector object.
     */
    CExcitationVector();
    
    /**
     Creates a excitation vector object.
     
     @param excitationType the single particle excitation type.
     @param braIndexes the vector of molecular orbital indexes on bra side.
     @param ketIndexes the vector of molecular orbital indexes on ket side.
     @param coefficients the vector of coefficients associates with single
            particle excitations.
     */
    CExcitationVector(const szblock               excitationType,
                      const std::vector<int32_t>& braIndexes,
                      const std::vector<int32_t>& ketIndexes,
                      const std::vector<double>&  coefficients);
    
    /**
     Creates a excitation vector object by copying other excitation vector object.
     
     @param source the excitation vector object.
     */
    CExcitationVector(const CExcitationVector& source);
    
    /**
     Creates a excitation vector object by moving other excitation vector object.
     
     @param source the excitation vector object.
     */
    CExcitationVector(CExcitationVector&& source) noexcept;
    
    /**
     Destroys a excitation vector object.
     */
    ~CExcitationVector();
    
    /**
     Assigns a excitation vector object by copying other excitation vector object.
     
     @param source the excitation vector object.
     */
    CExcitationVector& operator=(const CExcitationVector& source);
    
    /**
     Assigns a excitation vector object by moving other excitation vector object.
     
     @param source the excitation vector object.
     */
    CExcitationVector& operator=(CExcitationVector&& source) noexcept;
    
    /**
     Compares excitation vector object with other excitation vector object.
     
     @param other the excitation vector object.
     @return true if excitation vector objects are equal, false otherwise.
     */
    bool operator==(const CExcitationVector& other) const;
    
    /**
     Compares excitation vector object with other excitation vector object.
     
     @param other the excitation vector object.
     @return true if excitation vector objects are not equal, false otherwise.
     */
    bool operator!=(const CExcitationVector& other) const;
    
    /**
     Converts excitation vector object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the excitation vector object.
     */
    friend std::ostream& operator<<(      std::ostream&      output,
                                    const CExcitationVector& source);
};

#endif /* ExcitationVector_hpp */
