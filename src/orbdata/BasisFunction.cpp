//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "BasisFunction.hpp"

#include <utility>
#include <cmath>
#include <limits>

#include "MathConst.hpp"
#include "MpiFunc.hpp"

CBasisFunction::CBasisFunction()

    : _nContrVectors(0)

    , _angularMomentum(-1)
{

}

CBasisFunction::CBasisFunction(const std::vector<double>& exponents,
                               const std::vector<double>& normFactors,
                               const int32_t              nContrVectors,
                               const int32_t              angularMomentum)

    : _exponents(exponents)

    , _normFactors(normFactors)

    , _nContrVectors(nContrVectors)

    , _angularMomentum(angularMomentum)

{

}

CBasisFunction::CBasisFunction(const std::vector<double>& exponents,
                               const std::vector<double>& normFactors,
                               const int32_t              angularMomentum)

    : CBasisFunction(exponents, normFactors, 1, angularMomentum)

{
    
}

CBasisFunction::CBasisFunction(const CBasisFunction& source)

    : _exponents(source._exponents)

    , _normFactors(source._normFactors)

    , _nContrVectors(source._nContrVectors)

    , _angularMomentum(source._angularMomentum)
{

}

CBasisFunction::CBasisFunction(CBasisFunction&& source) noexcept

    : _exponents(std::move(source._exponents))

    , _normFactors(std::move(source._normFactors))

    , _nContrVectors(std::move(source._nContrVectors))

    , _angularMomentum(std::move(source._angularMomentum))
{

}

CBasisFunction::~CBasisFunction()
{

}

CBasisFunction&
CBasisFunction::operator=(const CBasisFunction& source)
{
    if (this == &source) return *this;

    _exponents = source._exponents;

    _normFactors = source._normFactors;
    
    _nContrVectors = source._nContrVectors;

    _angularMomentum = source._angularMomentum;

    return *this;
}

CBasisFunction&
CBasisFunction::operator=(CBasisFunction&& source) noexcept
{
    if (this == &source) return *this;

    _exponents = std::move(source._exponents);

    _normFactors = std::move(source._normFactors);
    
    _nContrVectors = std::move(source._nContrVectors);

    _angularMomentum = std::move(source._angularMomentum);

    return *this;
}

bool
CBasisFunction::operator==(const CBasisFunction& other) const
{
    if (_exponents.size() != other._exponents.size()) return false;

    for (size_t i = 0; i < _exponents.size(); i++)
    {
        if (std::fabs(_exponents[i] - other._exponents[i]) > 1.0e-13)
        {
            return false;
        }
    }

    if (_normFactors.size() != other._normFactors.size()) return false;

    for (size_t i = 0; i < _normFactors.size(); i++)
    {
        if (std::fabs(_normFactors[i] - other._normFactors[i]) > 1.0e-13)
        {
            return false;
        }
    }

    if (_nContrVectors != other._nContrVectors) return false;
    
    if (_angularMomentum != other._angularMomentum) return false;

    return true;
}

bool
CBasisFunction::operator!=(const CBasisFunction& other) const
{
    return !(*this == other);
}

void
CBasisFunction::setExponents(const std::vector<double>& exponents)
{
    _exponents = exponents;
}

void
CBasisFunction::setNormalizationFactors(const std::vector<double>& normFactors)
{
    _normFactors = normFactors;
}

void
CBasisFunction::setAngularMomentum(const int32_t angularMomentum)
{
    _angularMomentum = angularMomentum;
}

void
CBasisFunction::add(const double exponent,
                    const double normFactor)
{
    if (_nContrVectors == 0) _nContrVectors = 1;
    
    _exponents.push_back(exponent);

    _normFactors.push_back(normFactor);
}

void
CBasisFunction::normalize()
{
    // NOTE: Currently implemented for l = 0..6

    if (_angularMomentum > 6) return;

    // uncontracted basis, set expansion coeficient to 1.0

    if (_normFactors.size() == 1) _normFactors[0] = 1.0;

    // normalize primitive GBFs

    _rescale();
    
    // renormalize contracted GTOs
    
    auto edim = static_cast<int32_t>(_exponents.size());

    for (int32_t i = 0; i < _nContrVectors; i++)
    {
        // compute overlap
        
        double ovl = 0.0;
        
        for (int32_t j = 0; j < edim; j++)
        {
            ovl += _overlap(j, j, i);
            
            for (int32_t k = j + 1; k < edim; k++)
            {
                ovl += 2.0 * _overlap(j, k, i);
            }
        }
        
        ovl = 1.0 / std::sqrt(ovl);
        
        // renormaliza primitive BFs
        
        for (int32_t j = 0; j < edim; j++)
        {
            _normFactors[i * edim + j] *= ovl;
        }
        
    }
}

std::vector<double>
CBasisFunction::getExponents() const
{
    return _exponents;
}

std::vector<double>
CBasisFunction::getNormalizationFactors() const
{
    return _normFactors;
}

std::vector<double>
CBasisFunction::getNormalizationFactors(const int32_t iContrVector) const
{
    std::vector<double> normfacts;
    
    if (iContrVector < _nContrVectors)
    {
        auto edim = static_cast<int32_t>(_exponents.size());
        
        auto coff = edim * iContrVector;
        
        for (int32_t i = 0; i < edim; i++)
        {
            normfacts.push_back(_normFactors[coff + i]);
        }
    }
    
    return normfacts;
}

int32_t
CBasisFunction::getAngularMomentum() const
{
    return _angularMomentum;
}

int32_t
CBasisFunction::getNumberOfPrimitiveFunctions() const
{
    return static_cast<int32_t>(_exponents.size());
}

int32_t
CBasisFunction::getNumberOfContractedFunctions() const
{
    return _nContrVectors; 
}

int32_t
CBasisFunction::getNumberOfNormalizationFactors() const
{
    return static_cast<int32_t>(_normFactors.size()); 
}

void
CBasisFunction::_rescale()
{
    auto fpi = 2.0 / mathconst::getPiValue();
    
    auto edim = static_cast<int32_t>(_exponents.size());

    for (int32_t i = 0; i < edim; i++)
    {
        auto fact = std::pow(_exponents[i] * fpi, 0.75);
        
        for (int32_t j = 0; j < _nContrVectors; j++)
        {
            _normFactors[j * edim + i] *= fact;
        }
    }

    if (_angularMomentum == 1)
    {
        for (int32_t i = 0; i < edim; i++)
        {
            auto fact = 2.0 * std::sqrt(_exponents[i]);
            
            for (int32_t j = 0; j < _nContrVectors; j++)
            {
                _normFactors[j * edim + i] *= fact;
            }
        }
        
        return;
    }

    if (_angularMomentum == 2)
    {
        double f = 2.0 / std::sqrt(3.0);
        
        for (int32_t i = 0; i < edim; i++)
        {
            auto fact = f * _exponents[i];
            
            for (int32_t j = 0; j < _nContrVectors; j++)
            {
                _normFactors[j * edim + i] *= fact;
            }
        }

        return;
    }

    if (_angularMomentum == 3)
    {
        double f = 4.0 / std::sqrt(15.0);
        
        for (int32_t i = 0; i < edim; i++)
        {
            auto fact = f * _exponents[i] * std::sqrt(_exponents[i]);
            
            for (int32_t j = 0; j < _nContrVectors; j++)
            {
                _normFactors[j * edim + i] *= fact;
            }
        }

        return;
    }

    if (_angularMomentum == 4)
    {
        double f = 2.0 / std::sqrt(105.0);
        
        for (int32_t i = 0; i < edim; i++)
        {
            auto fact = f * _exponents[i] * _exponents[i];
            
            for (int32_t j = 0; j < _nContrVectors; j++)
            {
                _normFactors[j * edim + i] *= fact;
            }
        }

        return;
    }

    if (_angularMomentum == 5)
    {
        double f = 4.0 / std::sqrt(945.0);
        
        for (int32_t i = 0; i < edim; i++)
        {
            auto fact = f * _exponents[i] * _exponents[i]
            
                      * std::sqrt(_exponents[i]);
            
            for (int32_t j = 0; j < _nContrVectors; j++)
            {
                _normFactors[j * edim + i] *= fact;
            }
        }

        return;
    }

    if (_angularMomentum == 6)
    {
        double f = 4.0 / std::sqrt(10395.0);
        
        for (int32_t i = 0; i < edim; i++)
        {
            auto fact = f * _exponents[i] * _exponents[i]
            
                      * _exponents[i];
            
            for (int32_t j = 0; j < _nContrVectors; j++)
            {
                _normFactors[j * edim + i] *= fact;
            }
        }

        return;
    }
}

double
CBasisFunction::_overlap(const int32_t iComponent,
                         const int32_t jComponent,
                         const int32_t iContrVector) const
{
    auto fab = 1.0 / (_exponents[iComponent] + _exponents[jComponent]);

    auto coff = static_cast<int32_t>(_exponents.size()) * iContrVector;
    
    auto ovl = _normFactors[coff + iComponent] * _normFactors[coff + jComponent]

             * std::pow(mathconst::getPiValue() * fab, 1.5);

    if (_angularMomentum == 0) return ovl;

    if (_angularMomentum == 1) return  0.5 * fab * ovl;
    
    auto fab2 = fab * fab;

    if (_angularMomentum == 2) return 3.0 * fab2 * ovl;
    
    if (_angularMomentum == 3) return 7.5 * fab2 * fab * ovl;

    if (_angularMomentum == 4) return 420.0 * fab2 * fab2 * ovl;
    
    if (_angularMomentum == 5) return 1890.0 * fab2 * fab2 * fab * ovl;
    
    if (_angularMomentum == 6) return 41580.0 * fab2 * fab2 * fab2 * ovl;

    return std::numeric_limits<double>::quiet_NaN();
}

void
CBasisFunction::broadcast(int32_t  rank,
                          MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        mpi::bcast(_angularMomentum, comm);
        
        mpi::bcast(_nContrVectors, comm);

        mpi::bcast(_exponents, rank, comm);

        mpi::bcast(_normFactors, rank, comm);
    }
}

std::ostream&
operator<<(      std::ostream&   output, 
           const CBasisFunction& source)
{
    output << std::endl;

    output << "[CBasisFunction (Object):" << &source << "]" << std::endl;

    output << "_angularMomentum: " << source._angularMomentum << std::endl;
    
    output << "_nContrVectors: " << source._nContrVectors << std::endl;

    output << "_exponents: " << std::endl;

    for (size_t i = 0; i < source._exponents.size(); i++)
    {
        output << "_exponents[" << i << "]: " << std::endl;

        output << source._exponents[i] << std::endl;
    }

    output << "_normFactors: " << std::endl;

    for (size_t i = 0; i < source._normFactors.size(); i++)
    {
        output << "_normFactors[" << i  << "]: "<< std::endl;

        output << source._normFactors[i] << std::endl;
    }

    return output;
}
