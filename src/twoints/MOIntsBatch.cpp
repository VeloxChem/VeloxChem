//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "MOIntsBatch.hpp"

#include "DenseLinearAlgebra.hpp"

CMOIntsBatch::CMOIntsBatch()

    : _batchType(moints::oooo)
{
    
}

CMOIntsBatch::CMOIntsBatch(const std::vector<CDenseMatrix>& moIntegrals,
                           const std::vector<CTwoIndexes>&  generatorPairs,
                           const CTwoIndexes&               externalIndexes,
                           const moints                     batchType)

    : _moIntegrals(moIntegrals)

    , _generatorPairs(generatorPairs)

    , _externalIndexes(externalIndexes)

    , _batchType(batchType)
{
    
}

CMOIntsBatch::CMOIntsBatch(const CMOIntsBatch& source)

    : _moIntegrals(source._moIntegrals)

    , _generatorPairs(source._generatorPairs)

    , _externalIndexes(source._externalIndexes)

    , _batchType(source._batchType)
{
    
}

CMOIntsBatch::CMOIntsBatch(CMOIntsBatch&& source) noexcept

    : _moIntegrals(std::move(source._moIntegrals))

    , _generatorPairs(std::move(source._generatorPairs))

    , _externalIndexes(std::move(source._externalIndexes))

    , _batchType(std::move(source._batchType))
{
    
}

CMOIntsBatch::~CMOIntsBatch()
{
    
}

CMOIntsBatch&
CMOIntsBatch::operator=(const CMOIntsBatch& source)
{
    if (this == &source) return *this;
    
    _moIntegrals = source._moIntegrals;
    
    _generatorPairs = source._generatorPairs;
    
    _externalIndexes = source._externalIndexes;
    
    _batchType = source._batchType;
    
    return *this;
}

CMOIntsBatch&
CMOIntsBatch::operator=(CMOIntsBatch&& source) noexcept
{
    if (this == &source) return *this;
    
    _moIntegrals = std::move(source._moIntegrals);
    
    _generatorPairs = std::move(source._generatorPairs);
    
    _externalIndexes = std::move(source._externalIndexes);
    
    _batchType = std::move(source._batchType);
    
    return *this;
}

bool
CMOIntsBatch::operator==(const CMOIntsBatch& other) const
{
    if (_moIntegrals.size() != other._moIntegrals.size()) return false;
    
    for (size_t i = 0; i < _moIntegrals.size(); i++)
    {
        if (_moIntegrals[i] != other._moIntegrals[i]) return false;
    }
    
    if (_generatorPairs.size() != other._generatorPairs.size()) return false;
    
    for (size_t i = 0; i < _generatorPairs.size(); i++)
    {
        if (_generatorPairs[i] != other._generatorPairs[i]) return false;
    }
    
    if (_externalIndexes != other._externalIndexes) return false;
    
    if (_batchType != other._batchType) return false;
    
    return true;
}

bool
CMOIntsBatch::operator!=(const CMOIntsBatch& other) const
{
    return !(*this == other);
}

void
CMOIntsBatch::append(const CAOFockMatrix&        aoFockMatrix,
                     const CDenseMatrix&         braVector,
                     const CDenseMatrix&         ketVector,
                     const std::vector<int32_t>& braIndexes,
                     const std::vector<int32_t>& ketIndexes)
{
    for (int32_t i = 0; i < aoFockMatrix.getNumberOfFockMatrices(); i++)
    {
        auto tmat = denblas::multAB(aoFockMatrix.getReferenceToFock(i),
                                    ketVector);
        
        _moIntegrals.push_back(denblas::multAtB(braVector, tmat));
        
        _generatorPairs.push_back(CTwoIndexes(braIndexes[i], ketIndexes[i]));
    }
}

void
CMOIntsBatch::setBatchType(const moints batchType)
{
    _batchType = batchType; 
}

void
CMOIntsBatch::setExternalIndexes(const CTwoIndexes& externalIndexes)
{
    _externalIndexes = externalIndexes; 
}

const double*
CMOIntsBatch::getBatch(const int32_t iBatch)
{
    return _moIntegrals[iBatch].values();
}

std::ostream&
operator<<(      std::ostream& output,
           const CMOIntsBatch& source)
{
    output << std::endl;
    
    output << "[CMOIntsBatch (Object):" << &source << "]" << std::endl;
    
    output << "_moIntegrals: " << std::endl;
    
    for (size_t i = 0; i < source._moIntegrals.size(); i++)
    {
        output << "_moIntegrals[" << i << "]: " << std::endl;
        
        output << source._moIntegrals[i] << std::endl;
    }
    
    output << "_generatorPairs: " << std::endl;
    
    for (size_t i = 0; i < source._generatorPairs.size(); i++)
    {
        output << "_generatorPairs[" << i << "]: ";
        
        output << source._generatorPairs[i] << std::endl;
    }
    
    output << "_externalIndexes: " << source._externalIndexes << std::endl;
    
    output << "_batchType: "  << to_string(source._batchType) << std::endl;
    
    return output;
}
