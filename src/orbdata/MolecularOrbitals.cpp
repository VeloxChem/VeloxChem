//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "MolecularOrbitals.hpp"

#include "ErrorHandler.hpp"
#include "DenseLinearAlgebra.hpp"

CMolecularOrbitals::CMolecularOrbitals()

    : _orbitalsType(molorb::rest)

    , _orbitals(std::vector<CDenseMatrix>())

    , _energies(std::vector<CMemBlock<double>>())
{
    
}

CMolecularOrbitals::CMolecularOrbitals(const std::vector<CDenseMatrix>&      orbitals,
                                       const std::vector<CMemBlock<double>>& energies,
                                       const molorb                          orbitalsType)

    : _orbitalsType(orbitalsType)

    , _orbitals(orbitals)

    , _energies(energies)
{

}

CMolecularOrbitals::CMolecularOrbitals(const CMolecularOrbitals& source)

    : _orbitalsType(source._orbitalsType)

    , _orbitals(source._orbitals)

    , _energies(source._energies)
{
    
}

CMolecularOrbitals::CMolecularOrbitals(CMolecularOrbitals&& source) noexcept

    : _orbitalsType(std::move(source._orbitalsType))

    , _orbitals(std::move(source._orbitals))

    , _energies(std::move(source._energies))
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
    
    _energies = source._energies;
    
    return *this;
}

CMolecularOrbitals&
CMolecularOrbitals::operator=(CMolecularOrbitals&& source) noexcept
{
    if (this == &source) return *this;
    
    _orbitalsType = std::move(source._orbitalsType);
    
    _orbitals = std::move(source._orbitals);
    
    _energies = std::move(source._energies);
    
    return *this;
}

molorb
CMolecularOrbitals::getOrbitalsType() const
{
    return _orbitalsType;
}

int32_t
CMolecularOrbitals::getNumberOfRows() const
{
    return _orbitals[0].getNumberOfRows();
}

int32_t
CMolecularOrbitals::getNumberOfColumns() const
{
    return _orbitals[0].getNumberOfColumns();
}

const double*
CMolecularOrbitals::alphaOrbitals() const
{
    return _orbitals[0].values();
}

const double*
CMolecularOrbitals::betaOrbitals() const
{
    if (_orbitalsType == molorb::unrest)
    {
        return _orbitals[1].values();
    }
    
    return nullptr;
}

const double*
CMolecularOrbitals::alphaEnergies() const
{
    return _energies[0].data();
}

const double*
CMolecularOrbitals::betaEnergies() const
{
    if (_orbitalsType == molorb::unrest)
    {
        return _energies[1].data();
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
    
    if (_energies.size() != other._energies.size()) return false;

    for (size_t i = 0; i < _energies.size(); i++)
    {
        if (_energies[i] != other._energies[i]) return false;
    }

    return true;
}

bool
CMolecularOrbitals::operator!=(const CMolecularOrbitals& other) const
{
    return !(*this == other);
}

CAODensityMatrix
CMolecularOrbitals::getAODensity(const int32_t nElectrons) const
{
    if ((nElectrons % 2) == 0)
    {
        auto ndim = nElectrons / 2;
        
        auto nrow = _orbitals[0].getNumberOfRows();
        
        auto ncol = _orbitals[0].getNumberOfColumns();
        
        if (ndim <= ncol)
        {
            auto cmo = _orbitals[0].slice(0, 0, nrow, ndim);
            
            auto den = denblas::multABt(cmo, cmo);
            
            return CAODensityMatrix({den}, denmat::rest);
        }
    }
    
    return CAODensityMatrix();
}

CAODensityMatrix
CMolecularOrbitals::getAODensity(const int32_t nAlphaElectrons,
                                 const int32_t nBetaElectrons) const
{
    auto ndima = nAlphaElectrons;
    
    auto ndimb = nBetaElectrons;
    
    auto nrowa = _orbitals[0].getNumberOfRows();
    
    auto ncola = _orbitals[0].getNumberOfColumns();
    
    auto nrowb = _orbitals[1].getNumberOfRows();
    
    auto ncolb = _orbitals[1].getNumberOfColumns();
    
    if ((ndima <= ncola) && (ndimb <= ncolb))
    {
        auto cmoa = _orbitals[0].slice(0, 0, nrowa, ndima);
        
        auto dena = denblas::multABt(cmoa, cmoa);
        
        auto cmob = _orbitals[1].slice(0, 0, nrowb, ndimb);
        
        auto denb = denblas::multABt(cmob, cmob);
        
        return CAODensityMatrix({dena, denb}, denmat::unrest);
    }
    
    return CAODensityMatrix();
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
    
    output << "_energies: " << std::endl;
    
    for (size_t i = 0; i < source._energies.size(); i++)
    {
        output << "_energies[" << i << "]: "<< std::endl;
        
        output << source._energies[i] << std::endl;
    }
    
    return output;
}

