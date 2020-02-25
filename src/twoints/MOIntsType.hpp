//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef MOIntsType_hpp
#define MOIntsType_hpp

#include <string>

/**
 Enumerate class moints:
 
 Defines supported types of molecular integrals:
 moints::oooo - the (oo|oo) integrals
 moints::ooov - the (oo|ov) integrals
 moints::oovv - the (oo|vv) integrals
 moints::ovov - the (ov|ov) integrals
 moints::ovvv - the (ov|vv) integrals
 moints::vvvv - the (vv|vv) integrals
 moints::asym_oooo - the <oo||oo> integrals
 moints::asym_ooov - the <oo||ov> integrals
 moints::asym_oovv - the <oo||vv> integrals
 moints::asym_ovov - the <ov||ov> integrals
 moints::asym_ovvv - the <ov||vv> integrals
 moints::asym_vvvv - the <vv||vv> integrals
 */
enum class moints
{
    oooo,
    ooov,
    oovv,
    ovov,
    ovvv,
    vvvv,
    asym_oooo,
    asym_ooov,
    asym_oovv,
    asym_ovov,
    asym_ovvv,
    asym_vvvv
};

/**
 Converts enumerate class value to it's string label.
 
 @param moIntsType the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string to_string(const moints moIntsType)
{
    if (moIntsType == moints::oooo)
    {
        return std::string("(oo|oo)");
    }

    if (moIntsType == moints::ooov)
    {
        return std::string("(oo|ov)");
    }
    
    if (moIntsType == moints::oovv)
    {
        return std::string("(oo|vv)");
    }
    
    if (moIntsType == moints::ovov)
    {
        return std::string("(ov|ov)");
    }
    
    if (moIntsType == moints::ovvv)
    {
        return std::string("(ov|vv)");
    }
    
    if (moIntsType == moints::vvvv)
    {
        return std::string("(vv|vv)");
    }
    
    if (moIntsType == moints::asym_oooo)
    {
        return std::string("<oo||oo>");
    }
    
    if (moIntsType == moints::asym_ooov)
    {
        return std::string("<oo||ov>");
    }
    
    if (moIntsType == moints::asym_oovv)
    {
        return std::string("<oo||vv>");
    }
    
    if (moIntsType == moints::asym_ovov)
    {
        return std::string("<ov||ov>");
    }
    
    if (moIntsType == moints::asym_ovvv)
    {
        return std::string("<ov||vv>");
    }
    
    if (moIntsType == moints::asym_vvvv)
    {
        return std::string("<vv||vv>");
    }
    
    return std::string("UNKNOWN");
}

/**
 Checks if enumerate class value corresponds to antisymmetrized integrals.

 @param moIntsType  the enumerate class value.
 @return true if enumerate class value corresponds to antisymmetrized integrals,
         false otherwise.
 */
inline bool isAntisymmetrizedIntegrals(const moints moIntsType)
{
    if (moIntsType == moints::asym_oooo) return true;
    
    if (moIntsType == moints::asym_ooov) return true;
    
    if (moIntsType == moints::asym_oovv) return true;
    
    if (moIntsType == moints::asym_ovov) return true;
    
    if (moIntsType == moints::asym_ovvv) return true;
    
    if (moIntsType == moints::asym_vvvv) return true;
    
    return false;
}


#endif /* MOIntsType_hpp */
