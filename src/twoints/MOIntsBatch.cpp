//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

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
        
        if (isAntisymmetrizedIntegrals(_batchType))
        {
            tmat = denblas::multAB(aoFockMatrix.getReferenceToFock(i),
                                   braVector);
            
            _moIntegrals.push_back(denblas::multAtB(ketVector, tmat));
        }
        
        _generatorPairs.push_back(CTwoIndexes(braIndexes[i], ketIndexes[i]));
    }
}

void
CMOIntsBatch::appendMOInts(const CDenseMatrix& moIntegralMatrix,
                           const CTwoIndexes&  twoIndexes)
{
    _moIntegrals.push_back(moIntegralMatrix);

    _generatorPairs.push_back(twoIndexes);
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
CMOIntsBatch::getBatch(const int32_t iBatch) const
{
    return _moIntegrals[iBatch].values();
}

const double*
CMOIntsBatch::getBatchXY(const int32_t iBatch) const
{
    if (isAntisymmetrizedIntegrals(_batchType))
    {
        return _moIntegrals[2 * iBatch].values();
    }
    
    return nullptr;
}

const double*
CMOIntsBatch::getBatchYX(const int32_t iBatch) const
{
    if (isAntisymmetrizedIntegrals(_batchType))
    {
        return _moIntegrals[2 * iBatch + 1].values();
    }
    
    return nullptr;
}

const double*
CMOIntsBatch::getBatch(const CTwoIndexes& iGeneratorPair) const
{
    for (size_t i = 0; i < _moIntegrals.size(); i++)
    {
        if (iGeneratorPair == _generatorPairs[i])
            return _moIntegrals[i].values(); 
    }
    
    return nullptr;
}

moints
CMOIntsBatch::getBatchType() const
{
    return _batchType;
}

CTwoIndexes
CMOIntsBatch::getExternalIndexes() const
{
    return _externalIndexes;
}

std::vector<CTwoIndexes>
CMOIntsBatch::getGeneratorPairs() const
{
    return _generatorPairs;
}

int32_t
CMOIntsBatch::getNumberOfBatches() const
{
    return static_cast<int32_t>(_moIntegrals.size());
}

int32_t
CMOIntsBatch::getNumberOfRows() const
{
    if (_moIntegrals.size() > 0)
    {
        return _moIntegrals[0].getNumberOfRows();
    }
    
    return 0;
}

int32_t
CMOIntsBatch::getNumberOfColumns() const
{
    if (_moIntegrals.size() > 0)
    {
        return _moIntegrals[0].getNumberOfColumns();
    }
    
    return 0;
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
