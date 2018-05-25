//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "CMMParameters.hpp"

#include <utility>
#include <cmath>
#include <limits>

#include "MathConst.hpp"
#include "MpiFunc.hpp"

CCMMParameters::CCMMParameters()

    : _paramModel(cmmtyp::original)

    , _idElemental(-1)

    , _polarizability(1.0)

    , _exponent(1.0)

    , _refCoordNumber(0.0)

    , _factCoordNumber(0.0)
{
    
}

CCMMParameters::CCMMParameters(const cmmtyp               paramModel,
                               const int32_t              idElemental,
                               const double               polarizability,
                               const double               exponent,
                               const double               refCoordNumber,
                               const double               factCoordNumber,
                               const std::vector<double>& refFrequencies,
                               const std::vector<double>& refGammas,
                               const std::vector<double>& refScaleFactors)

    : _paramModel(paramModel)

    , _idElemental(idElemental)

    , _polarizability(polarizability)

    , _exponent(exponent)

    , _refCoordNumber(refCoordNumber)

    , _factCoordNumber(factCoordNumber)

    , _refFrequencies(refFrequencies)

    , _refGammas(refGammas)

    , _refScaleFactors(refScaleFactors)
{

}

CCMMParameters::CCMMParameters(const CCMMParameters& source)

    : _paramModel(source._paramModel)

    , _idElemental(source._idElemental)

    , _polarizability(source._polarizability)

    , _exponent(source._exponent)

    , _refCoordNumber(source._refCoordNumber)

    , _factCoordNumber(source._factCoordNumber)

    , _refFrequencies(source._refFrequencies)

    , _refGammas(source._refGammas)

    , _refScaleFactors(source._refScaleFactors)
{

}

CCMMParameters::CCMMParameters(CCMMParameters&& source) noexcept

    : _paramModel(std::move(source._paramModel))

    , _idElemental(std::move(source._idElemental))

    , _polarizability(std::move(source._polarizability))

    , _exponent(std::move(source._exponent))

    , _refCoordNumber(std::move(source._refCoordNumber))

    , _factCoordNumber(std::move(source._factCoordNumber))

    , _refFrequencies(std::move(source._refFrequencies))

    , _refGammas(std::move(source._refGammas))

    , _refScaleFactors(std::move(source._refScaleFactors))
{

}

CCMMParameters::~CCMMParameters()
{

}

CCMMParameters&
CCMMParameters::operator=(const CCMMParameters& source)
{
    if (this == &source) return *this;

    _paramModel = source._paramModel;
    
    _idElemental = source._idElemental;
    
    _polarizability = source._polarizability;
    
    _exponent = source._exponent;
    
    _refCoordNumber = source._refCoordNumber;
    
    _factCoordNumber = source._factCoordNumber;
    
    _refFrequencies = source._refFrequencies;
    
    _refGammas = source._refGammas;
    
    _refScaleFactors = source._refScaleFactors;

    return *this;
}

CCMMParameters&
CCMMParameters::operator=(CCMMParameters&& source) noexcept
{
    if (this == &source) return *this;

    _paramModel = std::move(source._paramModel);
    
    _idElemental = std::move(source._idElemental);
    
    _polarizability = std::move(source._polarizability);
    
    _exponent = std::move(source._exponent);
    
    _refCoordNumber = std::move(source._refCoordNumber);
    
    _factCoordNumber = std::move(source._factCoordNumber);
    
    _refFrequencies = std::move(source._refFrequencies);
    
    _refGammas = std::move(source._refGammas);
    
    _refScaleFactors = std::move(source._refScaleFactors);

    return *this;
}

bool
CCMMParameters::operator==(const CCMMParameters& other) const
{
    if (_paramModel != other._paramModel) return false;
    
    if (_idElemental != other._idElemental) return false;
    
    if (std::fabs(_polarizability - other._polarizability) > 1.0e-13)
    {
        return false;
    }
    
    if (std::fabs(_exponent - other._exponent) > 1.0e-13) return false;
    
    if (std::fabs(_refCoordNumber - other._refCoordNumber) > 1.0e-13)
    {
        return false;
    }
    
    if (std::fabs(_factCoordNumber- other._factCoordNumber) > 1.0e-13)
    {
        return false;
    }

    if (_refFrequencies.size() != other._refFrequencies.size()) return false;

    for (size_t i = 0; i < _refFrequencies.size(); i++)
    {
        if (std::fabs(_refFrequencies[i] - other._refFrequencies[i]) > 1.0e-13)
        {
            return false;
        }
    }
    
    if (_refGammas.size() != other._refGammas.size()) return false;
    
    for (size_t i = 0; i < _refGammas.size(); i++)
    {
        if (std::fabs(_refGammas[i] - other._refGammas[i]) > 1.0e-13)
        {
            return false;
        }
    }
    
    if (_refScaleFactors.size() != other._refScaleFactors.size()) return false;
    
    for (size_t i = 0; i < _refScaleFactors.size(); i++)
    {
        if (std::fabs(_refScaleFactors[i] - other._refScaleFactors[i]) > 1.0e-13)
        {
            return false;
        }
    }

    return true;
}

bool
CCMMParameters::operator!=(const CCMMParameters& other) const
{
    return !(*this == other);
}

void
CCMMParameters::setForceFieldModel(const cmmtyp paramModel)
{
    _paramModel = paramModel; 
}

void
CCMMParameters::setIdElemental(const int32_t idElemental)
{
    _idElemental = idElemental;
}

void
CCMMParameters::setPolarizability(const double polarizability)
{
    _polarizability = polarizability;
}

void
CCMMParameters::setExponent(const double exponent)
{
    _exponent = exponent;
}

void
CCMMParameters::setCoordinationNumber(const double refCoordNumber)
{
    _refCoordNumber = refCoordNumber;
}

void
CCMMParameters::setCoordinationFactor(const double factCoordNumber)
{
    _factCoordNumber = factCoordNumber;
}

void
CCMMParameters::addFrequency(const double refFrequency,
                             const double refGamma,
                             const double refScaleFactor)
{
    _refFrequencies.push_back(refFrequency);
    
    _refGammas.push_back(refGamma);
    
    _refScaleFactors.push_back(refScaleFactor); 
}

cmmtyp
CCMMParameters::getForceFieldModel() const
{
    return _paramModel; 
}

int32_t
CCMMParameters::getIdElemental() const
{
    return _idElemental;
}

double
CCMMParameters::getPolarizability() const
{
    return _polarizability;
}

double
CCMMParameters::getExponent() const
{
    return _exponent;
}

double
CCMMParameters::getCoordinationNumber() const
{
    return _refCoordNumber;
}

double
CCMMParameters::getCoordinationFactor() const
{
    return _factCoordNumber;
}

std::ostream&
operator<<(      std::ostream&   output, 
           const CCMMParameters& source)
{
    output << std::endl;

    output << "[CCMMParameters (Object):" << &source << "]" << std::endl;
    
    output << "_paramModel: " << to_string(source._paramModel) << std::endl;
    
    output << "_idElemental: " << source._idElemental << std::endl;
    
    output << "_polarizability: " << source._polarizability << std::endl;
    
    output << "_exponent: " << source._exponent << std::endl;
    
    output << "_refCoordNumber: " << source._refCoordNumber << std::endl;
    
    output << "_factCoordNumber: " << source._factCoordNumber << std::endl;
    
    output << "_refFrequencies:" << std::endl;

    for (size_t i = 0; i < source._refFrequencies.size(); i++)
    {
        output << "_refFrequencies[" << i << "]: " << std::endl;

        output << source._refFrequencies[i] << std::endl;
    }

    output << "_refGammas: " << std::endl;

    for (size_t i = 0; i < source._refGammas.size(); i++)
    {
        output << "_refGammas[" << i  << "]: "<< std::endl;

        output << source._refGammas[i] << std::endl;
    }
    
    output << "_refScaleFactors: " << std::endl;
    
    for (size_t i = 0; i < source._refScaleFactors.size(); i++)
    {
        output << "_refScaleFactors[" << i  << "]: "<< std::endl;
        
        output << source._refScaleFactors[i] << std::endl;
    }

    return output;
}
