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

    : _orbitalsType(molorb::rest)

    , _orbitals(std::vector<CDenseMatrix>())
{
    
}

CMolecularOrbitals::CMolecularOrbitals(const std::vector<CDenseMatrix>& orbitals,
                                       const molorb                     orbitalsType)

    : _orbitalsType(orbitalsType)

    , _orbitals(orbitals)
{
    
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

int32_t
CMolecularOrbitals::getNumberOfOrbitalsMatrices() const
{
    // restricted molecular orbital matrix
    
    if (_orbitalsType == molorb::rest)
    {
        return static_cast<int32_t>(_orbitals.size());
    }
    
    // unrestricted molecular orbital matrix
    
    if (_orbitalsType == molorb::unrest)
    {
        return static_cast<int32_t>(_orbitals.size()) / 2;
    }
    
    return 0;
}

molorb
CMolecularOrbitals::getOrbitalsType() const
{
    return _orbitalsType;
}

int32_t
CMolecularOrbitals::getNumberOfRows(const int32_t iOrbitalsMatrix) const
{
    // restricted molecular orbital matrix
    
    if (_orbitalsType == molorb::rest && iOrbitalsMatrix < getNumberOfOrbitalsMatrices())
    {
        return _orbitals[iOrbitalsMatrix].getNumberOfRows();
    }
    
    // unrestricted molecular orbital matrix
    
    if (_orbitalsType == molorb::unrest && iOrbitalsMatrix < getNumberOfOrbitalsMatrices())
    {
        return _orbitals[2 * iOrbitalsMatrix].getNumberOfRows();
    }
    
    return 0;
}

int32_t
CMolecularOrbitals::getNumberOfColumns(const int32_t iOrbitalsMatrix) const
{
    // restricted molecular orbital matrix
    
    if (_orbitalsType == molorb::rest && iOrbitalsMatrix < getNumberOfOrbitalsMatrices())
    {
        return _orbitals[iOrbitalsMatrix].getNumberOfColumns();
    }
    
    // unrestricted molecular orbital matrix
    
    if (_orbitalsType == molorb::unrest && iOrbitalsMatrix < getNumberOfOrbitalsMatrices())
    {
        return _orbitals[2 * iOrbitalsMatrix].getNumberOfColumns();
    }
    
    return 0;
}

const double*
CMolecularOrbitals::totalOrbitals(const int32_t iOrbitalsMatrix) const
{
    if (_orbitalsType == molorb::rest && iOrbitalsMatrix < getNumberOfOrbitalsMatrices())
    {
        return _orbitals[iOrbitalsMatrix].values();
    }
    
    return nullptr;
}

const double*
CMolecularOrbitals::alphaOrbitals(const int32_t iOrbitalsMatrix) const
{
    if (_orbitalsType == molorb::unrest && iOrbitalsMatrix < getNumberOfOrbitalsMatrices())
    {
        return _orbitals[2 * iOrbitalsMatrix].values();
    }
    
    return nullptr;
}

const double*
CMolecularOrbitals::betaOrbitals(const int32_t iOrbitalsMatrix) const
{
    if (_orbitalsType == molorb::unrest && iOrbitalsMatrix < getNumberOfOrbitalsMatrices())
    {
        return _orbitals[2 * iOrbitalsMatrix + 1].values();
    }
    
    return nullptr;
}

std::string
CMolecularOrbitals::getString() const
{
    std::string orb_str;

    orb_str += "Orbitals Type: " + to_string(_orbitalsType) + "\n";

    for (size_t i = 0; i < _orbitals.size(); i++)
    {
        orb_str += _orbitals[i].getString();
    }

    return orb_str;
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

