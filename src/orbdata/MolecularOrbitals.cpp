//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "MolecularOrbitals.hpp"

#include "ErrorHandler.hpp"

CMolecularOrbitals::CMolecularOrbitals()

    : _orbitalsType(morb::rest)

    , _orbitals(std::vector<CDenseMatrix>())
{
    
}

CMolecularOrbitals::CMolecularOrbitals(const std::vector<CDenseMatrix>& orbitals,
                                       const morb                       orbitalsType)

    : _orbitalsType(orbitalsType)

    , _orbitals(orbitals)
{
    if (orbitalsType == morb::unrest)
    {
        errors::assertMsgCritical(
                                  orbitals.size() % 2 == 0,
                                  "MolecularOrbitals - Odd number of matrices for unrestricted molecular orbitals"
                                  );
    }
}

CMolecularOrbitals::CMolecularOrbitals(const CMolecularOrbitals& source)

    : _orbitalsType(source._orbitalsType)

    , _orbitals(source._orbitals)
{
    
}

CMolecularOrbitals::CMolecularOrbitals(CMolecularOrbitals&& source) noexcept

    : _orbitalsType(std::move(source._orbitalsType))

    , _orbitals(std::move(source._orbitals))
{
    
}

CMolecularOrbitals::~CMolecularOrbitals()
{
    
}

CMolecularOrbitals&
CMolecularOrbitals::operator=(const CMolecularOrbitals& source)
{
    if (this == &source) return *this;
    
    _orbitalsType = source._orbitalsType;
    
    _orbitals = source._orbitals;
    
    return *this;
}

CMolecularOrbitals&
CMolecularOrbitals::operator=(CMolecularOrbitals&& source) noexcept
{
    if (this == &source) return *this;
    
    _orbitalsType = std::move(source._orbitalsType);
    
    _orbitals = std::move(source._orbitals);
    
    return *this;
}

bool
CMolecularOrbitals::operator==(const CMolecularOrbitals& other) const
{
    if (_orbitalsType != other._orbitalsType) return false;
    
    if (_orbitals.size() != other._orbitals.size()) return false;
    
    for (size_t i = 0; i < _orbitals.size(); i++)
    {
        if (_orbitals[i] != other._orbitals[i]) return false;
    }
    
    return true;
}

bool
CMolecularOrbitals::operator!=(const CMolecularOrbitals& other) const
{
    return !(*this == other);
}

std::ostream&
operator<<(      std::ostream&       output,
           const CMolecularOrbitals& source)
{
    output << std::endl;
    
    output << "[CMolecularOrbitals (Object):" << &source << "]" << std::endl;
    
    output << "_orbitalsType: " << to_string(source._orbitalsType) << std::endl;
    
    output << "_orbitals: " << std::endl;
    
    for (size_t i = 0; i < source._orbitals.size(); i++)
    {
        output << "_orbitals[" << i << "]: "<< std::endl;
        
        output << source._orbitals[i] << std::endl;
    }
    
    return output;
}

