//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "CauchySchwarzScreener.hpp"

#include <utility>
#include <cmath>

CCauchySchwarzScreener::CCauchySchwarzScreener()

    : _screeningScheme(ericut::qq)

    , _braQValues(CMemBlock<double>())

    , _ketQValues(CMemBlock<double>())

    , _braPairExtends(CMemBlock<double>())

    , _ketPairExtends(CMemBlock<double>())

    , _threshold(1.0e-15)
{
    
}

CCauchySchwarzScreener::CCauchySchwarzScreener(const CMemBlock<double>& braQValues,
                                               const CMemBlock<double>& ketQValues,
                                               const CGtoPairsBlock&    braGtoPairsBlock,
                                               const CGtoPairsBlock&    ketGtoPairsBlock,
                                               const ericut             screeningScheme,
                                               const double             threshold)
    : _screeningScheme(screeningScheme)

    , _braQValues(braQValues)

    , _ketQValues(ketQValues)

    , _braPairExtends(CMemBlock<double>())

    , _ketPairExtends(CMemBlock<double>())

    , _threshold(threshold)
{
    // compute GTOs pairs entents for QQR screening scheme
    
    if (_screeningScheme == ericut::qqr)
    {
        _braPairExtends = CMemBlock<double>(braGtoPairsBlock.getNumberOfScreenedContrPairs());
        
        _setPairsExtents(_braPairExtends, braGtoPairsBlock);
        
        _ketPairExtends = CMemBlock<double>(ketGtoPairsBlock.getNumberOfScreenedContrPairs());
        
        _setPairsExtents(_ketPairExtends, ketGtoPairsBlock);
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

const double*
CCauchySchwarzScreener::getBraQValues() const
{
    return _braQValues.data();
}

const double*
CCauchySchwarzScreener::getKetQValues() const
{
    return _ketQValues.data();
}

bool
CCauchySchwarzScreener::isEmpty() const
{
    return ((_braQValues.size() == 0) || (_ketQValues.size() == 0)); 
}

void
CCauchySchwarzScreener::setScreeningVector(      CMemBlock<int32_t>& qqVector,
                                           const CMemBlock<double>&  pqDistances,
                                           const bool                isBraEqualKet,
                                           const int32_t             iContrPair) const
{
    // initialize vector
    
    mathfunc::set_to(qqVector.data(), 1, qqVector.size());
    
    // set up pointer to screening vector
    
    auto qqvec = qqVector.data();
    
    // set up pointer to Q values on ket side
    
    auto kqvals = _ketQValues.data();
    
    // set up dimensions on ket side
    
    auto kdim = _ketQValues.size();
    
    if (isBraEqualKet) kdim = iContrPair + 1;
    
    // original Cauchy-Schwarz screening scheme
    
    if ((_screeningScheme == ericut::qq) ||
        (_screeningScheme == ericut::qqden))
    {
        auto fbq = _braQValues.at(iContrPair);
        
        for (int32_t i = 0; i < kdim; i++)
        {
            if ((fbq * kqvals[i]) < _threshold) qqvec[i] = 0;
        }
    }
    
    // distance dependent Cauchy-Schwarz screening scheme
    
    if ((_screeningScheme == ericut::qqr) ||
        (_screeningScheme == ericut::qqrden))
    {
        // set up pointer to effective PQ distances
        
        auto rpq = pqDistances.data();
        
        // set up pointer to GTOs pairs extents on ket side
        
        auto kext = _braPairExtends.data();
        
        // data for bra side
        
        auto bext = _braPairExtends.at(iContrPair);
        
        auto fbq  = _braQValues.at(iContrPair);
    
        for (int32_t i = 0; i < kdim; i++)
        {
            auto r = rpq[i] - bext - kext[i];
            
            double fact = (r > 1.0) ? 1.0 / r : 1.0;
            
            if ((fbq * kqvals[i] * fact) < _threshold) qqvec[i] = 0;
        }
    }
}

void
CCauchySchwarzScreener::setScreeningVector(      CMemBlock<int32_t>& qqVector,
                                           const CMemBlock<int32_t>& qqIndexes,
                                           const CMemBlock<double>&  maxDensityElements,
                                           const CMemBlock<double>&  pqDistances,
                                           const int32_t             nContrPairs,
                                           const int32_t             iContrPair) const
{
    // initialize vector
    
    qqVector.zero();
    
    mathfunc::set_to(qqVector.data(), 1, nContrPairs);
    
    // set up pointer to screening vector
    
    auto qqvec = qqVector.data();
    
    // set up pointer to indexes
    
    auto qqidx = qqIndexes.data();
    
    // set up pointer to Q values on ket side
    
    auto kqvals = _ketQValues.data();
    
    // set up pointer to max. density elements
    
    auto mden = maxDensityElements.data();
    
    // original Cauchy-Schwarz + AO density screening scheme
    
    if (_screeningScheme == ericut::qqden)
    {
        auto fbq = _braQValues.at(iContrPair);
        
        for (int32_t i = 0; i < nContrPairs; i++)
        {
            if ((fbq * kqvals[qqidx[i]] * mden[i]) < _threshold) qqvec[i] = 0;
        }
    }
    
    // distance dependent Cauchy-Schwarz + AO density screening scheme
    
    if (_screeningScheme == ericut::qqrden)
    {
        // set up pointer to effective PQ distances
        
        auto rpq = pqDistances.data();
        
        // set up pointer to GTOs pairs extents on ket side
        
        auto kext = _braPairExtends.data();
        
        // data for bra side
        
        auto bext = _braPairExtends.at(iContrPair);
        
        auto fbq  = _braQValues.at(iContrPair);
        
        for (int32_t i = 0; i < nContrPairs; i++)
        {
            auto idx = qqidx[i];
            
            auto r = rpq[idx] - bext - kext[idx];
            
            double fact = (r > 1.0) ? 1.0 / r : 1.0;
            
            if ((fbq * kqvals[idx] * fact * mden[i]) < _threshold) qqvec[i] = 0;
        }
    }
}

void
CCauchySchwarzScreener::_setPairsExtents(      CMemBlock<double>& gtoPairExtents,
                                         const CGtoPairsBlock&    gtoPairsBlock)
{
    // threshold scaling factor
    
    auto ferf = 1.0 / std::erfc(_threshold);
    
    // set up dimensions of GTOs pairs vector
    
    auto ndim = gtoPairsBlock.getNumberOfScreenedContrPairs();
    
    // set up pointers to 1 / (e_a + e_b) factors
    
    auto foxi = gtoPairsBlock.getFactorsOneOverXi();
    
    // set up pointers to coordinates of primitive P center
    
    auto rpx = gtoPairsBlock.getCoordinatesPX();
    
    auto rpy = gtoPairsBlock.getCoordinatesPY();
    
    auto rpz = gtoPairsBlock.getCoordinatesPZ();
    
    // set up pointers to coordinates of effective centers
    
    auto repx = gtoPairsBlock.getEffectiveCoordinatesPX();
    
    auto repy = gtoPairsBlock.getEffectiveCoordinatesPY();
    
    auto repz = gtoPairsBlock.getEffectiveCoordinatesPZ();
    
    // set pointers to start and end position in primitive GTOs pairs vector
    
    auto spos = gtoPairsBlock.getStartPositions();
    
    auto epos = gtoPairsBlock.getEndPositions();
    
    // loop over GTOs pairs
    
    for (int32_t i = 0; i < ndim; i++)
    {
        double rext = 0.0;

        for (int32_t j = spos[i]; j < epos[i]; j++)
        {
            // distances R(P_eff - P) = P_eff - P
            
            auto dpx = repx[i] - rpx[j];
            
            auto dpy = repy[i] - rpy[j];
            
            auto dpz = repz[i] - rpz[j];
            
            // primitive extent
            
            auto pext = std::sqrt(dpx * dpx +  dpy * dpy + dpz * dpz)
            
                      + std::sqrt(2.0 * foxi[j]) * ferf;
            
            if (pext > rext) rext = pext;
        }
        
        gtoPairExtents.at(i) = rext;
    }
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

