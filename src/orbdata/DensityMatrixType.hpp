//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef DensityMatrixType_hpp
#define DensityMatrixType_hpp

#include <string>

/**
 Enumerate class denmat:
 
 Defines supported density matrix types:
 denmat::rest   - the restricted density matrix
 denmat::unrest - the unrestricted density matrix
 denmat::rmoij  - the restricted C_i C_j^T matrix
 denmat::umoij  - the unrestricted C_i C_j^T matrix
 */
enum class denmat
{
    rest,
    unrest,
    rmoij,
    umoij
};

/**
 Converts key value of density matrix type to integer number.
 
 @param denMatrix the enumerate class value.
 @return the integer number.
 */
inline int32_t to_int(const denmat denMatrix)
{
    return static_cast<int32_t>(denMatrix);
}

/**
 Converts integer key value to density matrix type.
 
 @param keyValue the integer key value.
 @return the density matrix type.
 */
inline denmat to_denmat(const int32_t keyValue)
{
    if (keyValue == to_int(denmat::rest)) return denmat::rest;
    
    if (keyValue == to_int(denmat::unrest)) return denmat::unrest;
    
    if (keyValue == to_int(denmat::rmoij)) return denmat::rmoij;
    
    if (keyValue == to_int(denmat::umoij)) return denmat::umoij;
    
    return denmat::rest;
}

/**
 Converts enumerate class value to it's string label.
 
 @param denMatrix the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string to_string(const denmat denMatrix)
{
    if (denMatrix == denmat::rest)
    {
        return std::string("Restricted Density Matrix");
    }
    
    if (denMatrix == denmat::unrest)
    {
        return std::string("Unrestricted Density Matrix");
    }
    
    if (denMatrix == denmat::rmoij)
    {
        return std::string("Restricted C_iC_j^T Density Matrix");
    }
    
    if (denMatrix == denmat::umoij)
    {
        return std::string("Unrestricted C_iC_j^T Density Matrix");
    }
   
    return std::string("UNKNOWN");
}

#endif /* DensityMatrixType_hpp */
