//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "TwoIntsDistributor.hpp"

#include "StringFormat.hpp"

CTwoIntsDistribution::CTwoIntsDistribution()

    : _distPattern(dist2e::batch)

    , _needSyncLock(false)

    , _nRows(0)

    , _nColumns(0)

    , _intsData(nullptr)
{
    
}

CTwoIntsDistribution::CTwoIntsDistribution(      double* intsData,
                                           const int32_t nRows,
                                           const int32_t nColumns,
                                           const dist2e  distPattern)

    : _distPattern(distPattern)

    , _needSyncLock(false)

    , _nRows(nRows)

    , _nColumns(nColumns)

    , _intsData(intsData)
{
    
}

CTwoIntsDistribution::CTwoIntsDistribution(const CTwoIntsDistribution& source)

    : _distPattern(source._distPattern)

    , _needSyncLock(source._needSyncLock)

    , _nRows(source._nRows)

    , _nColumns(source._nColumns)
{
    _intsData = source._intsData;
}

CTwoIntsDistribution::~CTwoIntsDistribution()
{
    
}

CTwoIntsDistribution&
CTwoIntsDistribution::operator=(const CTwoIntsDistribution& source)
{
    if (this == &source) return *this;
    
    _distPattern = source._distPattern;
    
    _needSyncLock = source._needSyncLock;
    
    _nRows = source._nRows;
    
    _nColumns = source._nColumns;
    
    _intsData = source._intsData;
    
    return *this;
}

bool
CTwoIntsDistribution::operator==(const CTwoIntsDistribution& other) const
{
    if (_distPattern != other._distPattern) return false;
    
    if (_needSyncLock != other._needSyncLock) return false;
    
    if (_nRows != other._nRows) return false;
    
    if (_nColumns != other._nColumns) return false;
    
    if (_intsData != other._intsData) return false;
    
    return true;
}

bool
CTwoIntsDistribution::operator!=(const CTwoIntsDistribution& other) const
{
    return !(*this == other);
}

std::ostream&
operator<<(      std::ostream&         output,
           const CTwoIntsDistribution& source)
{
    output << std::endl;
    
    output << "[CTwoIntsDistribution (Object):" << &source << "]" << std::endl;
    
    output << "_distPattern: " << to_string(source._distPattern) << std::endl;
    
    output << "_needSyncLock: " << fstr::to_string(source._needSyncLock) << std::endl;
    
    output << "_nRows: " << source._nRows << std::endl;
    
    output << "_nColumns: " << source._nColumns << std::endl;
    
    output << "_intsData: " << source._intsData << std::endl;
    
    return output;
}
