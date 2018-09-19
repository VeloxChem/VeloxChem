//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "TwoIntsDistributor.hpp"

#include <cmath>

#include "StringFormat.hpp"
#include "AngularMomentum.hpp"

CTwoIntsDistribution::CTwoIntsDistribution()

    : _distPattern(dist2e::batch)

    , _needSyncLock(false)

    , _nRows(0)

    , _nColumns(0)

    , _idGtoPair(-1)

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

    , _idGtoPair(-1)

    , _intsData(intsData)
{
    
}

CTwoIntsDistribution::CTwoIntsDistribution(      double* intsData,
                                           const int32_t nRows,
                                           const int32_t nColumns,
                                           const int32_t idGtoPair,
                                           const dist2e  distPattern)
    : _distPattern(distPattern)

    , _needSyncLock(false)

    , _nRows(nRows)

    , _nColumns(nColumns)

    , _idGtoPair(idGtoPair)

    , _intsData(intsData)
{
    // override identifier of GTOs pair if not used in distributing integrals
    
    if (_distPattern != dist2e::qvalues) _idGtoPair = -1;
}


CTwoIntsDistribution::CTwoIntsDistribution(const CTwoIntsDistribution& source)

    : _distPattern(source._distPattern)

    , _needSyncLock(source._needSyncLock)

    , _nRows(source._nRows)

    , _nColumns(source._nColumns)

    , _idGtoPair(source._idGtoPair)
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
    
    _idGtoPair = source._idGtoPair;
    
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
    
    if (_idGtoPair != other._idGtoPair) return false;
    
    if (_intsData != other._intsData) return false;
    
    return true;
}

bool
CTwoIntsDistribution::operator!=(const CTwoIntsDistribution& other) const
{
    return !(*this == other);
}

bool
CTwoIntsDistribution::needSyncLock() const
{
    return _needSyncLock;
}

void
CTwoIntsDistribution::distribute(const CMemBlock2D<double>& spherInts,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const bool                 isBraEqualKet,
                                 const int32_t              iContrPair)
{
    // distribute two electron integrals into data batch
    
    if (_distPattern == dist2e::batch)
    {
        _distSpherIntsIntoBatch(spherInts, braGtoPairsBlock, ketGtoPairsBlock,
                                isBraEqualKet, iContrPair);
        
        return; 
    }
    
    // distribute two electron integrals into qvalues vector
    
    if (_distPattern == dist2e::qvalues)
    {
        _distSpherIntsIntoQValues(spherInts, braGtoPairsBlock, ketGtoPairsBlock,
                                  isBraEqualKet, iContrPair);
        
        return;
    }
}

void
CTwoIntsDistribution::_distSpherIntsIntoBatch(const CMemBlock2D<double>& spherInts,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const bool                 isBraEqualKet,
                                              const int32_t              iContrPair)
{
    // determine number of angular components in shhell
    
    auto ncomp = angmom::to_SphericalComponents(braGtoPairsBlock.getBraAngularMomentum(),
                                                braGtoPairsBlock.getKetAngularMomentum())
    
               * angmom::to_SphericalComponents(ketGtoPairsBlock.getBraAngularMomentum(),
                                                ketGtoPairsBlock.getKetAngularMomentum());
    
    // determine starting poisition in integrals batch
    
    auto spos = _getStartIndexForBatch(ncomp, iContrPair, isBraEqualKet);
    
    // determine dimensions of GTOs pairs vector on ket side
    
    auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

    if (isBraEqualKet) kdim = iContrPair + 1;
    
    // store integrals into batch of integrals
    
    for (int32_t i = 0; i < kdim; i++)
    {
        for (int32_t j = 0; j < ncomp; j++)
        {
            _intsData[spos + j] = (spherInts.data(j))[i];
        }
        
        spos += ncomp;
    }
}

void
CTwoIntsDistribution::_distSpherIntsIntoQValues(const CMemBlock2D<double>& spherInts,
                                                const CGtoPairsBlock&      braGtoPairsBlock,
                                                const CGtoPairsBlock&      ketGtoPairsBlock,
                                                const bool                 isBraEqualKet,
                                                const int32_t              iContrPair)
{
    // check if dimensions on bra and ket sides are correct
    
    if ((braGtoPairsBlock.getNumberOfScreenedContrPairs() != 1) ||
        (ketGtoPairsBlock.getNumberOfScreenedContrPairs() != 1))
    {
        // error handling code...
        
        return;
    }
    
    // determine number of angular components in shhell
    
    auto ncomp = angmom::to_SphericalComponents(braGtoPairsBlock.getBraAngularMomentum(),
                                                braGtoPairsBlock.getKetAngularMomentum())
    
               * angmom::to_SphericalComponents(ketGtoPairsBlock.getBraAngularMomentum(),
                                                ketGtoPairsBlock.getKetAngularMomentum());
    
    //  determine maximum integral value in shell
    
    auto mval = (spherInts.data(0))[0];
    
    for (int32_t j = 1; j < ncomp; j++)
    {
        auto fval = (spherInts.data(j))[0];
        
        if (fval > mval) mval = fval;
    }
    
    // compute Q value for contracted pair
    
    _intsData[_idGtoPair] = std::sqrt(mval);
}

int32_t
CTwoIntsDistribution::_getStartIndexForBatch(const int32_t nShellComponents,
                                             const int32_t iContrPair,
                                             const bool    isBraEqualKet) const
{
    int32_t idx = 0;
    
    if (isBraEqualKet)
    {
        for (int32_t i = 0; i < iContrPair; i++) idx += i + 1; 
    }
    else
    {
        idx = iContrPair * _nColumns;
    }
    
    return idx * nShellComponents;
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
    
    output << "_idGtoPair: " << source._idGtoPair << std::endl;
    
    output << "_intsData: " << source._intsData << std::endl;
    
    return output;
}
