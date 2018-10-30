//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef FockMatrixType_hpp
#define FockMatrixType_hpp

#include <string>

/**
 Enumerate class fockmat:
 
 Defines supported Fock matrix types:
 fockmat::restjk  - the restricted Fock matrix (Coulomb + exchange)
 fockmat::restjkx - the restricted scaled Fock matrix (Coulomb + scaled exchange)
 fockmat::restj   - the restricted Coulomb matrix
 fockmat::restk   - the restricted exchange matrix
 fockmat::restkx  - the restricted scaled exchange matrix
 */
enum class fockmat
{
    restjk,
    restjkx,
    restj,
    restk,
    restkx
};

/**
 Converts enumerate class value to it's string label.
 
 @param fockMatrix the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string to_string(const fockmat fockMatrix)
{
    if (fockMatrix == fockmat::restjk)
    {
        return std::string("Restricted J + K Matrix");
    }
   
    if (fockMatrix == fockmat::restjkx)
    {
        return std::string("Restricted J + xK Matrix");
    }
    
    if (fockMatrix == fockmat::restj)
    {
        return std::string("Restricted J Matrix");
    }
    
    if (fockMatrix == fockmat::restk)
    {
        return std::string("Restricted K Matrix");
    }
    
    if (fockMatrix == fockmat::restkx)
    {
        return std::string("Restricted xK Matrix");
    }
    
    return std::string("UNKNOWN");
}

#endif /* FockMatrixType_hpp */
