#ifndef FockType_hpp
#define FockType_hpp

#include <string>

/**
 Enumerate class fock_t:
 
 Defines supported Fock matrix types:
 fock_t::restjk  - the restricted Fock matrix (Coulomb + exchange)
 fock_t::restjkx - the restricted scaled Fock matrix (Coulomb + scaled exchange)
 fock_t::restj   - the restricted Coulomb matrix
 fock_t::restk   - the restricted exchange matrix
 fock_t::restkx  - the restricted scaled exchange matrix
 fock_t::rgenjk  - the restricted general Fock matrix (Coulomb + exchange)
 fock_t::rgenjkx - the restricted scaled general Fock matrix (Coulomb + scaled exchange)
 fock_t::rgenj   - the restricted general Coulomb matrix
 fock_t::rgenk   - the restricted general exchange matrix
 fock_t::rgenkx  - the restricted scaled general exchange matrix
 fock_t::unrestjk  - the unrestricted Fock matrix (Coulomb + exchange)
 fock_t::unrestjkx - the unrestricted scaled Fock matrix (Coulomb + scaled exchange)
 fock_t::unrestj  - the unrestricted Coulomb matrix
 */
enum class fock_t
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
 
 @param matrix the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string to_string(const fock_t matrix)
{
    if (matrix == fock_t::restjk)
    {
        return std::string("Restricted 2J + K Matrix");
    }
   
    if (matrix == fock_t::restjkx)
    {
        return std::string("Restricted 2J + xK Matrix");
    }
    
    if (matrix == fock_t::restj)
    {
        return std::string("Restricted J Matrix");
    }
    
    if (matrix == fock_t::restk)
    {
        return std::string("Restricted K Matrix");
    }
    
    if (matrix == fock_t::restkx)
    {
        return std::string("Restricted xK Matrix");
    }
    
    if (matrix == fock_t::rgenjk)
    {
        return std::string("Restricted general 2J + K Matrix");
    }
    
    if (matrix == fock_t::rgenjkx)
    {
        return std::string("Restricted general 2J + xK Matrix");
    }
    
    if (matrix == fock_t::rgenj)
    {
        return std::string("Restricted general J Matrix");
    }
    
    if (matrix == fock_t::rgenk)
    {
        return std::string("Restricted general K Matrix");
    }
    
    if (matrix == fock_t::rgenkx)
    {
        return std::string("Restricted general xK Matrix");
    }

    if (matrix == fock_t::unrestjk)
    {
        return std::string("Unrestricted 2J + K Matrix");
    }

    if (matrix == fock_t::unrestj)
    {
        return std::string("Unrestricted 2J Matrix");
    }
    
    if (matrix == fock_t::unrestjkx)
    {
        return std::string("Unrestricted 2J + xK Matrix");
    }
    
    return std::string("UNKNOWN");
}

#endif /* matrixType_hpp */
