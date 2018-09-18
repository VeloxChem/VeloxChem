//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "CauchySchwarzScreener.hpp"

#include <utility>
#include <cmath>

CCauchySchwarzScreener::CCauchySchwarzScreener()

    : _screeningScheme(ericut::qq)

    , _threshold(1.0e-13)
{
    
}

CCauchySchwarzScreener::CCauchySchwarzScreener(const CGtoPairsBlock& braGtoPairsBlock,
                                               const CGtoPairsBlock& ketGtoPairsBlock,
                                               const ericut          screeningScheme,
                                               const double          threshold)
    : _screeningScheme(screeningScheme)

    , _threshold(threshold)

    , _braQValues(CMemBlock<double>(braGtoPairsBlock.getNumberOfScreenedContrPairs()))

    , _ketQValues(CMemBlock<double>(ketGtoPairsBlock.getNumberOfScreenedContrPairs()))
{
    if (_screeningScheme == ericut::qqr)
    {
        _braPairExtends = CMemBlock<double>(braGtoPairsBlock.getNumberOfScreenedContrPairs());
        
        _ketPairExtends = CMemBlock<double>(ketGtoPairsBlock.getNumberOfScreenedContrPairs());
        
        // TODO: add computation of extends....
    }
}

CCauchySchwarzScreener::CCauchySchwarzScreener(const CCauchySchwarzScreener& source)

    : _screeningScheme(source._screeningScheme)

    , _braQValues(source._braQValues)

    , _ketQValues(source._ketQValues)

    , _braPairExtends(source._braPairExtends)

    , _ketPairExtends(source._ketPairExtends)

    , _threshold(source._threshold)
{
    
}

CCauchySchwarzScreener::CCauchySchwarzScreener(CCauchySchwarzScreener&& source) noexcept

    : _screeningScheme(std::move(source._screeningScheme))

    , _braQValues(std::move(source._braQValues))

    , _ketQValues(std::move(source._ketQValues))

    , _braPairExtends(std::move(source._braPairExtends))

    , _ketPairExtends(std::move(source._ketPairExtends))

    , _threshold(std::move(source._threshold))
{
    
}

CCauchySchwarzScreener::~CCauchySchwarzScreener()
{
    
}

CCauchySchwarzScreener&
CCauchySchwarzScreener::operator=(const CCauchySchwarzScreener& source)
{
    if (this == &source) return *this;
    
    _screeningScheme = source._screeningScheme;
    
    _braQValues = source._braQValues;
    
    _ketQValues = source._ketQValues;
    
    _braPairExtends = source._braPairExtends;
    
    _ketPairExtends = source._ketPairExtends;
    
    _threshold = source._threshold;
    
    return *this;
}

CCauchySchwarzScreener&
CCauchySchwarzScreener::operator=(CCauchySchwarzScreener&& source) noexcept
{
    if (this == &source) return *this;
    
    _screeningScheme = std::move(source._screeningScheme);
    
    _braQValues = std::move(source._braQValues);
    
    _ketQValues = std::move(source._ketQValues);
    
    _braPairExtends = std::move(source._braPairExtends);
    
    _ketPairExtends = std::move(source._ketPairExtends);
    
    _threshold = std::move(source._threshold);
    
    return *this;
}

bool
CCauchySchwarzScreener::operator==(const CCauchySchwarzScreener& other) const
{
    if (_screeningScheme != other._screeningScheme) return false;
    
    if (_braQValues != other._braQValues) return false;
    
    if (_ketQValues != other._ketQValues) return false;
    
    if (_braPairExtends != other._braPairExtends) return false;
    
    if (_ketPairExtends != other._ketPairExtends) return false;
    
    if (std::fabs(_threshold - other._threshold) > 1.0e-13) return false;
    
    return true;
}

bool
CCauchySchwarzScreener::operator!=(const CCauchySchwarzScreener& other) const
{
    return !(*this == other);
}

void
CCauchySchwarzScreener::setScreeningScheme(const ericut screeningScheme)
{
    _screeningScheme = screeningScheme;
}

void
CCauchySchwarzScreener::setThreshold(const double threshold)
{
    _threshold = threshold; 
}

ericut
CCauchySchwarzScreener::getScreeningScheme() const
{
    return _screeningScheme;
}

double
CCauchySchwarzScreener::getThreshold() const
{
    return _threshold; 
}

double*
CCauchySchwarzScreener::getBraQValues()
{
    return _braQValues.data();
}

const double*
CCauchySchwarzScreener::getBraQValues() const
{
    return _braQValues.data();
}

double*
CCauchySchwarzScreener::getKetQValues()
{
    return _ketQValues.data();
}

const double*
CCauchySchwarzScreener::getKetQValues() const
{
    return _ketQValues.data();
}

std::ostream&
operator<<(      std::ostream&           output,
           const CCauchySchwarzScreener& source)
{
    output << std::endl;
    
    output << "[CCauchySchwarzScreener (Object):" << &source << "]" << std::endl;
    
    output << "_screeningScheme: " << to_string(source._screeningScheme) << std::endl;
    
    output << "_braQValues: " << source._braQValues << std::endl;
    
    output << "_ketQValues: " << source._ketQValues << std::endl;
    
    output << "_braPairExtends: " << source._braPairExtends << std::endl;
    
    output << "_ketPairExtends: " << source._ketPairExtends << std::endl;
    
    output << "_threshold: " << source._threshold << std::endl;
    
    return output;
}

