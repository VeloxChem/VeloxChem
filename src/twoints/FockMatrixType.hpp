//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

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
 fockmat::rgenjk  - the restricted general Fock matrix (Coulomb + exchange)
 fockmat::rgenjkx - the restricted scaled general Fock matrix (Coulomb + scaled exchange)
 fockmat::rgenj   - the restricted general Coulomb matrix
 fockmat::rgenk   - the restricted general exchange matrix
 fockmat::rgenkx  - the restricted scaled general exchange matrix
 fockmat::unrestjk  - the unrestricted Fock matrix (Coulomb + exchange)
 fockmat::unrestjkx - the unrestricted scaled Fock matrix (Coulomb + scaled exchange)
 fockmat::unrestj   - the unrestricted Coulomb matrix
 */
enum class fockmat
{
    restjk,
    restjkx,
    restj,
    restk,
    restkx,
    rgenjk,
    rgenjkx,
    rgenj,
    rgenk,
    rgenkx,
    unrestjk,
    unrestj,
    unrestjkx
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
        return std::string("Restricted 2J + K Matrix");
    }
   
    if (fockMatrix == fockmat::restjkx)
    {
        return std::string("Restricted 2J + xK Matrix");
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
    
    if (fockMatrix == fockmat::rgenjk)
    {
        return std::string("Restricted general 2J + K Matrix");
    }
    
    if (fockMatrix == fockmat::rgenjkx)
    {
        return std::string("Restricted general 2J + xK Matrix");
    }
    
    if (fockMatrix == fockmat::rgenj)
    {
        return std::string("Restricted general J Matrix");
    }
    
    if (fockMatrix == fockmat::rgenk)
    {
        return std::string("Restricted general K Matrix");
    }
    
    if (fockMatrix == fockmat::rgenkx)
    {
        return std::string("Restricted general xK Matrix");
    }

    if (fockMatrix == fockmat::unrestjk)
    {
        return std::string("Unrestricted 2J + K Matrix");
    }

    if (fockMatrix == fockmat::unrestj)
    {
        return std::string("Unrestricted 2J Matrix");
    }
    
    if (fockMatrix == fockmat::unrestjkx)
    {
        return std::string("Unrestricted 2J + xK Matrix");
    }
    
    return std::string("UNKNOWN");
}

#endif /* FockMatrixType_hpp */
