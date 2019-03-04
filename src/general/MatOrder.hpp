//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MatOrder_hpp
#define MatOrder_hpp

/**
 Enumerate class matorder:
 
 Defines all allowed key values for matrix data ordering in memory:
 matorder::col_major - the column major ordering of data
 matorder::row_major - the row major ordering of data
 */

enum class matorder : int32_t
{
    col_major,
    row_major
};


#endif /* MatOrder_hpp */
