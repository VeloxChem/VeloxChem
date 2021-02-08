//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "FockSubMatrix.hpp"

#include "AngularMomentum.hpp"

CFockSubMatrix::CFockSubMatrix()

    : _subFockMatrices({})

    , _startPositionsA(CMemBlock<int32_t>())

    , _startPositionsB(CMemBlock<int32_t>())

    , _startPositionsC(CMemBlock<int32_t>())

    , _startPositionsD(CMemBlock<int32_t>())

    , _dimSubMatrixA(0)

    , _dimSubMatrixB(0)

    , _dimSubMatrixC(0)

    , _dimSubMatrixD(0)
{
    
}

CFockSubMatrix::CFockSubMatrix(const std::vector<CMemBlock2D<double>>& subFockMatrices,
                               const CMemBlock<int32_t>&               startPostionsA,
                               const CMemBlock<int32_t>&               startPostionsB,
                               const CMemBlock<int32_t>&               startPostionsC,
                               const CMemBlock<int32_t>&               startPostionsD,
                               const int32_t                           dimSubMatrixA,
                               const int32_t                           dimSubMatrixB,
                               const int32_t                           dimSubMatrixC,
                               const int32_t                           dimSubMatrixD)

    : _subFockMatrices(subFockMatrices)

    , _startPositionsA(startPostionsA)

    , _startPositionsB(startPostionsB)

    , _startPositionsC(startPostionsC)

    , _startPositionsD(startPostionsD)

    , _dimSubMatrixA(dimSubMatrixA)

    , _dimSubMatrixB(dimSubMatrixB)

    , _dimSubMatrixC(dimSubMatrixC)

    , _dimSubMatrixD(dimSubMatrixD)
{
    
}

CFockSubMatrix::CFockSubMatrix(const CGtoPairsBlock& braGtoPairsBlock,
                               const CGtoPairsBlock& ketGtoPairsBlock,
                               const fockmat         fockType)
{
    // set up angular momentum components for bra side
    
    auto aang = braGtoPairsBlock.getBraAngularMomentum();
    
    auto bang = braGtoPairsBlock.getKetAngularMomentum();
    
    auto acomp = angmom::to_SphericalComponents(aang);
    
    auto bcomp = angmom::to_SphericalComponents(bang);
    
    // set up angular momentum components for ket side
    
    auto cang = ketGtoPairsBlock.getBraAngularMomentum();
    
    auto dang = ketGtoPairsBlock.getKetAngularMomentum();
    
    auto ccomp = angmom::to_SphericalComponents(cang);
    
    auto dcomp = angmom::to_SphericalComponents(dang);
    
    // initialize submatrices dimensions
    
    _dimSubMatrixA = braGtoPairsBlock.getNumberOfRowsInBraMatrix();
    
    _dimSubMatrixB = braGtoPairsBlock.getNumberOfRowsInKetMatrix();
    
    _dimSubMatrixC = ketGtoPairsBlock.getNumberOfRowsInBraMatrix();
    
    _dimSubMatrixD = ketGtoPairsBlock.getNumberOfRowsInKetMatrix();
    
    // initialize starting positions
    
    _startPositionsA = CMemBlock<int32_t>(acomp);
    
    _startPositionsB = CMemBlock<int32_t>(bcomp);
    
    _startPositionsC = CMemBlock<int32_t>(ccomp);
    
    _startPositionsD = CMemBlock<int32_t>(dcomp);
    
    for (int32_t i = 0; i < acomp; i++)
    {
        _startPositionsA.at(i) = braGtoPairsBlock.getBraMatrixPosition(i);
    }
    
    for (int32_t i = 0; i < bcomp; i++)
    {
        _startPositionsB.at(i) = braGtoPairsBlock.getKetMatrixPosition(i);
    }
    
    for (int32_t i = 0; i < ccomp; i++)
    {
        _startPositionsC.at(i) = ketGtoPairsBlock.getBraMatrixPosition(i);
    }
    
    for (int32_t i = 0; i < dcomp; i++)
    {
        _startPositionsD.at(i) = ketGtoPairsBlock.getKetMatrixPosition(i);
    }
    
    // allocate submatrices data
    
    _allocSubMatrices(fockType, acomp, bcomp, ccomp, dcomp);
}

CFockSubMatrix::CFockSubMatrix(const CFockSubMatrix& source)

    : _subFockMatrices(source._subFockMatrices)

    , _startPositionsA(source._startPositionsA)

    , _startPositionsB(source._startPositionsB)

    , _startPositionsC(source._startPositionsC)

    , _startPositionsD(source._startPositionsD)

    , _dimSubMatrixA(source._dimSubMatrixA)

    , _dimSubMatrixB(source._dimSubMatrixB)

    , _dimSubMatrixC(source._dimSubMatrixC)

    , _dimSubMatrixD(source._dimSubMatrixD)
{
    
}

CFockSubMatrix::CFockSubMatrix(CFockSubMatrix&& source) noexcept

    : _subFockMatrices(std::move(source._subFockMatrices))

    , _startPositionsA(std::move(source._startPositionsA))

    , _startPositionsB(std::move(source._startPositionsB))

    , _startPositionsC(std::move(source._startPositionsC))

    , _startPositionsD(std::move(source._startPositionsD))

    , _dimSubMatrixA(std::move(source._dimSubMatrixA))

    , _dimSubMatrixB(std::move(source._dimSubMatrixB))

    , _dimSubMatrixC(std::move(source._dimSubMatrixC))

    , _dimSubMatrixD(std::move(source._dimSubMatrixD))
{
    
}

CFockSubMatrix::~CFockSubMatrix()
{
    
}

CFockSubMatrix&
CFockSubMatrix::operator=(const CFockSubMatrix& source)
{
    if (this == &source) return *this;
    
    _subFockMatrices = source._subFockMatrices;
    
    _startPositionsA = source._startPositionsA;
    
    _startPositionsB = source._startPositionsB;
    
    _startPositionsC = source._startPositionsC;
    
    _startPositionsD = source._startPositionsD;
    
    _dimSubMatrixA = source._dimSubMatrixA;
    
    _dimSubMatrixB = source._dimSubMatrixB;
    
    _dimSubMatrixC = source._dimSubMatrixC;
    
    _dimSubMatrixD = source._dimSubMatrixD;
    
    return *this;
}

CFockSubMatrix&
CFockSubMatrix::operator=(CFockSubMatrix&& source) noexcept
{
    if (this == &source) return *this;
    
    _subFockMatrices = std::move(source._subFockMatrices);
    
    _startPositionsA = std::move(source._startPositionsA);
    
    _startPositionsB = std::move(source._startPositionsB);
    
    _startPositionsC = std::move(source._startPositionsC);
    
    _startPositionsD = std::move(source._startPositionsD);
    
    _dimSubMatrixA = std::move(source._dimSubMatrixA);
    
    _dimSubMatrixB = std::move(source._dimSubMatrixB);
    
    _dimSubMatrixC = std::move(source._dimSubMatrixC);
    
    _dimSubMatrixD = std::move(source._dimSubMatrixD);
    
    return *this;
}

bool
CFockSubMatrix::operator==(const CFockSubMatrix& other) const
{
    if (_subFockMatrices.size() != other._subFockMatrices.size())
    {
        return false;
    }
    
    for (size_t i = 0; i < _subFockMatrices.size(); i++)
    {
        if (_subFockMatrices[i] != other._subFockMatrices[i]) return false;
    }
    
    if (_startPositionsA != other._startPositionsA) return false;
    
    if (_startPositionsB != other._startPositionsB) return false;
    
    if (_startPositionsC != other._startPositionsC) return false;
    
    if (_startPositionsD != other._startPositionsD) return false;
    
    if (_dimSubMatrixA != other._dimSubMatrixA) return false;
    
    if (_dimSubMatrixB != other._dimSubMatrixB) return false;
    
    if (_dimSubMatrixC != other._dimSubMatrixC) return false;
    
    if (_dimSubMatrixD != other._dimSubMatrixD) return false;
    
    return true;
}

bool
CFockSubMatrix::operator!=(const CFockSubMatrix& other) const
{
    return !(*this == other);
}

void
CFockSubMatrix::accumulate(      double* aoFockMatrix,
                           const int32_t nColumns,
                           const fockmat fockType) const
{
    // restricted/unrestricted Coulomb + exchange matrix
    
    if ((fockType == fockmat::restjk) || (fockType == fockmat::restjkx) ||
        (fockType == fockmat::unrestjk) || (fockType == fockmat::unrestjkx))
    {
        _addContribution(aoFockMatrix, nColumns, 0,
                         _startPositionsA, _startPositionsB,
                         _dimSubMatrixA, _dimSubMatrixB);
        
        _addContribution(aoFockMatrix, nColumns, 1,
                         _startPositionsC, _startPositionsD,
                         _dimSubMatrixC, _dimSubMatrixD);
        
        _addContribution(aoFockMatrix, nColumns, 2,
                         _startPositionsA, _startPositionsC,
                         _dimSubMatrixA, _dimSubMatrixC);
        
        _addContribution(aoFockMatrix, nColumns, 3,
                         _startPositionsA, _startPositionsD,
                         _dimSubMatrixA, _dimSubMatrixD);
        
        _addContribution(aoFockMatrix, nColumns, 4,
                         _startPositionsB, _startPositionsC,
                         _dimSubMatrixB, _dimSubMatrixC);
        
        _addContribution(aoFockMatrix, nColumns, 5,
                         _startPositionsB, _startPositionsD,
                         _dimSubMatrixB, _dimSubMatrixD);
    }
    
    // restricted Coulomb matrix
    
    if ((fockType == fockmat::restj) || (fockType == fockmat::unrestj))
    {
        _addContribution(aoFockMatrix, nColumns, 0,
                         _startPositionsA, _startPositionsB,
                         _dimSubMatrixA, _dimSubMatrixB);
        
        _addContribution(aoFockMatrix, nColumns, 1,
                         _startPositionsC, _startPositionsD,
                         _dimSubMatrixC, _dimSubMatrixD);
    }
    
    // restricted exchange matrix
    
    if ((fockType == fockmat::restk) || (fockType == fockmat::restkx))
    {
        _addContribution(aoFockMatrix, nColumns, 0,
                         _startPositionsA, _startPositionsC,
                         _dimSubMatrixA, _dimSubMatrixC);
        
        _addContribution(aoFockMatrix, nColumns, 1,
                         _startPositionsA, _startPositionsD,
                         _dimSubMatrixA, _dimSubMatrixD);
        
        _addContribution(aoFockMatrix, nColumns, 2,
                         _startPositionsB, _startPositionsC,
                         _dimSubMatrixB, _dimSubMatrixC);
        
        _addContribution(aoFockMatrix, nColumns, 3,
                         _startPositionsB, _startPositionsD,
                         _dimSubMatrixB, _dimSubMatrixD);
    }
    
    // restricted general Coulomb matrix
    
    if (fockType == fockmat::rgenj)
    {
        _addContribution(aoFockMatrix, nColumns, 0,
                         _startPositionsA, _startPositionsB,
                         _dimSubMatrixA, _dimSubMatrixB);
        
        _addContribution(aoFockMatrix, nColumns, 1,
                         _startPositionsB, _startPositionsA,
                         _dimSubMatrixB, _dimSubMatrixA);
        
        _addContribution(aoFockMatrix, nColumns, 2,
                         _startPositionsC, _startPositionsD,
                         _dimSubMatrixC, _dimSubMatrixD);
        
        _addContribution(aoFockMatrix, nColumns, 3,
                         _startPositionsD, _startPositionsC,
                         _dimSubMatrixD, _dimSubMatrixC);
    }
    
    // restricted general exchange matrix
    
    if ((fockType == fockmat::rgenk) || (fockType == fockmat::rgenkx))
    {
        _addContribution(aoFockMatrix, nColumns, 0,
                         _startPositionsA, _startPositionsC,
                         _dimSubMatrixA, _dimSubMatrixC);
        
        _addContribution(aoFockMatrix, nColumns, 1,
                         _startPositionsC, _startPositionsA,
                         _dimSubMatrixC, _dimSubMatrixA);
        
        _addContribution(aoFockMatrix, nColumns, 2,
                         _startPositionsA, _startPositionsD,
                         _dimSubMatrixA, _dimSubMatrixD);
        
        _addContribution(aoFockMatrix, nColumns, 3,
                         _startPositionsD, _startPositionsA,
                         _dimSubMatrixD, _dimSubMatrixA);
        
        _addContribution(aoFockMatrix, nColumns, 4,
                         _startPositionsB, _startPositionsC,
                         _dimSubMatrixB, _dimSubMatrixC);
        
        _addContribution(aoFockMatrix, nColumns, 5,
                         _startPositionsC, _startPositionsB,
                         _dimSubMatrixC, _dimSubMatrixB);
        
        _addContribution(aoFockMatrix, nColumns, 6,
                         _startPositionsB, _startPositionsD,
                         _dimSubMatrixB, _dimSubMatrixD);
        
        _addContribution(aoFockMatrix, nColumns, 7,
                         _startPositionsD, _startPositionsB,
                         _dimSubMatrixD, _dimSubMatrixB);
    }
    
    // restricted general Coulomb + exchange matrix
    
    if ((fockType == fockmat::rgenjk) || (fockType == fockmat::rgenjkx))
    {
        _addContribution(aoFockMatrix, nColumns, 0,
                         _startPositionsA, _startPositionsB,
                         _dimSubMatrixA, _dimSubMatrixB);
        
        _addContribution(aoFockMatrix, nColumns, 1,
                         _startPositionsB, _startPositionsA,
                         _dimSubMatrixB, _dimSubMatrixA);
        
        _addContribution(aoFockMatrix, nColumns, 2,
                         _startPositionsC, _startPositionsD,
                         _dimSubMatrixC, _dimSubMatrixD);
        
        _addContribution(aoFockMatrix, nColumns, 3,
                         _startPositionsD, _startPositionsC,
                         _dimSubMatrixD, _dimSubMatrixC);
        
        _addContribution(aoFockMatrix, nColumns, 4,
                         _startPositionsA, _startPositionsC,
                         _dimSubMatrixA, _dimSubMatrixC);
        
        _addContribution(aoFockMatrix, nColumns, 5,
                         _startPositionsC, _startPositionsA,
                         _dimSubMatrixC, _dimSubMatrixA);
        
        _addContribution(aoFockMatrix, nColumns, 6,
                         _startPositionsA, _startPositionsD,
                         _dimSubMatrixA, _dimSubMatrixD);
        
        _addContribution(aoFockMatrix, nColumns, 7,
                         _startPositionsD, _startPositionsA,
                         _dimSubMatrixD, _dimSubMatrixA);
        
        _addContribution(aoFockMatrix, nColumns, 8,
                         _startPositionsB, _startPositionsC,
                         _dimSubMatrixB, _dimSubMatrixC);
        
        _addContribution(aoFockMatrix, nColumns, 9,
                         _startPositionsC, _startPositionsB,
                         _dimSubMatrixC, _dimSubMatrixB);
        
        _addContribution(aoFockMatrix, nColumns, 10,
                         _startPositionsB, _startPositionsD,
                         _dimSubMatrixB, _dimSubMatrixD);
        
        _addContribution(aoFockMatrix, nColumns, 11,
                         _startPositionsD, _startPositionsB,
                         _dimSubMatrixD, _dimSubMatrixB);
    }
}

double*
CFockSubMatrix::getSubMatrixData(const int32_t iMatrix,
                                 const int32_t iComponent)
{
    return _subFockMatrices[iMatrix].data(iComponent);
}

const int32_t* 
CFockSubMatrix::getStartPositionsA() const
{
    return _startPositionsA.data();
}

const int32_t*
CFockSubMatrix::getStartPositionsB() const
{
    return _startPositionsB.data();
}

const int32_t*
CFockSubMatrix::getStartPositionsC() const
{
    return _startPositionsC.data();
}

const int32_t*
CFockSubMatrix::getStartPositionsD() const
{
    return _startPositionsD.data();
}

int32_t
CFockSubMatrix::getDimensionsA() const
{
    return _dimSubMatrixA;
}

int32_t
CFockSubMatrix::getDimensionsB() const
{
    return _dimSubMatrixB;
}

int32_t
CFockSubMatrix::getDimensionsC() const
{
    return _dimSubMatrixC;
}

int32_t
CFockSubMatrix::getDimensionsD() const
{
    return _dimSubMatrixD;
}

void
CFockSubMatrix::_allocSubMatrices(const fockmat fockType,
                                  const int32_t nComponentsA,
                                  const int32_t nComponentsB,
                                  const int32_t nComponentsC,
                                  const int32_t nComponentsD)
{
    // restricted Coulomb matrix
    
    if ((fockType == fockmat::restj) || (fockType == fockmat::unrestj))
    {
        auto matdim = _dimSubMatrixA * _dimSubMatrixB;
        
        auto ncomp  = nComponentsA * nComponentsB;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixC * _dimSubMatrixD;
        
        ncomp  = nComponentsC * nComponentsD;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
    }
    
    // restricted exchange matrix
    
    if ((fockType == fockmat::restk) || (fockType == fockmat::restkx))
    {
        auto matdim = _dimSubMatrixA * _dimSubMatrixC;
        
        auto ncomp  = nComponentsA * nComponentsC;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixA * _dimSubMatrixD;
        
        ncomp  = nComponentsA * nComponentsD;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixB * _dimSubMatrixC;
        
        ncomp  = nComponentsB * nComponentsC;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixB * _dimSubMatrixD;
        
        ncomp  = nComponentsB * nComponentsD;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
    }
    
    // restricted/unrestricted Coulomb + exchange matrix
    
    if ((fockType == fockmat::restjk) || (fockType == fockmat::restjkx) ||
        (fockType == fockmat::unrestjk) || (fockType == fockmat::unrestjkx))
    {
        auto matdim = _dimSubMatrixA * _dimSubMatrixB;
        
        auto ncomp  = nComponentsA * nComponentsB;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixC * _dimSubMatrixD;
        
        ncomp  = nComponentsC * nComponentsD;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixA * _dimSubMatrixC;
        
        ncomp  = nComponentsA * nComponentsC;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixA * _dimSubMatrixD;
        
        ncomp  = nComponentsA * nComponentsD;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixB * _dimSubMatrixC;
        
        ncomp  = nComponentsB * nComponentsC;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixB * _dimSubMatrixD;
        
        ncomp  = nComponentsB * nComponentsD;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
    }
    
    // restricted general Coulomb matrix
    
    if (fockType == fockmat::rgenj)
    {
        auto matdim = _dimSubMatrixA * _dimSubMatrixB;
        
        auto ncomp  = nComponentsA * nComponentsB;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixC * _dimSubMatrixD;
        
        ncomp  = nComponentsC * nComponentsD;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
    }
    
    // restricted general exchange matrix
    
    if ((fockType == fockmat::rgenk) || (fockType == fockmat::rgenkx))
    {
        auto matdim = _dimSubMatrixA * _dimSubMatrixC;
        
        auto ncomp  = nComponentsA * nComponentsC;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixA * _dimSubMatrixD;
        
        ncomp  = nComponentsA * nComponentsD;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixB * _dimSubMatrixC;
        
        ncomp  = nComponentsB * nComponentsC;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixB * _dimSubMatrixD;
        
        ncomp  = nComponentsB * nComponentsD;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
    }
    
    // restricted general Coulomb + exchange matrix
    
    if ((fockType == fockmat::rgenjk) || (fockType == fockmat::rgenjkx))
    {
        auto matdim = _dimSubMatrixA * _dimSubMatrixB;
        
        auto ncomp  = nComponentsA * nComponentsB;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixC * _dimSubMatrixD;
        
        ncomp  = nComponentsC * nComponentsD;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixA * _dimSubMatrixC;
        
        ncomp  = nComponentsA * nComponentsC;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixA * _dimSubMatrixD;
        
        ncomp  = nComponentsA * nComponentsD;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixB * _dimSubMatrixC;
        
        ncomp  = nComponentsB * nComponentsC;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        matdim = _dimSubMatrixB * _dimSubMatrixD;
        
        ncomp  = nComponentsB * nComponentsD;
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
        
        _subFockMatrices.push_back(CMemBlock2D<double>(matdim, ncomp));
    }

    // initialize all submatrices to zero
    
    for (size_t i = 0; i < _subFockMatrices.size(); i++)
    {
        _subFockMatrices[i].zero(); 
    }
}

void
CFockSubMatrix::_addContribution(      double*             aoFockMatrix,
                                 const int32_t             nColumns,
                                 const int32_t             iMatrix,
                                 const CMemBlock<int32_t>& braComponents,
                                 const CMemBlock<int32_t>& ketComponents,
                                 const int32_t             braDimensions,
                                 const int32_t             ketDimensions) const
{
    auto bcomp = braComponents.size();
    
    auto kcomp = ketComponents.size();
 
    for (int32_t i = 0; i < bcomp; i++)
    {
        auto istart = braComponents.at(i);
        
        for(int32_t j = 0; j < kcomp; j++)
        {
            auto jstart = ketComponents.at(j);
            
            auto pdat = _subFockMatrices[iMatrix].data(i * kcomp + j);
            
            for (int32_t k = 0; k < braDimensions; k++)
            {
                auto koff = (istart + k) * nColumns;
                
                for (int32_t l = 0; l < ketDimensions; l++)
                {
                    auto loff = koff + jstart + l;
                    
                    aoFockMatrix[loff] += pdat[k * ketDimensions + l];
                }
            }
        }
    }
}

std::ostream&
operator<<(      std::ostream&   output,
           const CFockSubMatrix& source)
{
    output << std::endl;
    
    output << "[CFockSubMatrix (Object):" << &source << "]" << std::endl;
    
    output << "_subFockMatrices: " << std::endl;
    
    for (size_t i = 0; i < source._subFockMatrices.size(); i++)
    {
        output << "_subFockMatrices[" << i << "]: " << std::endl;
        
        output << source._subFockMatrices[i] << std::endl;
    }
    
    output << "_startPositionsA: " << source._startPositionsA << std::endl;
    
    output << "_startPositionsB: " << source._startPositionsB << std::endl;
    
    output << "_startPositionsC: " << source._startPositionsC << std::endl;
    
    output << "_startPositionsD: " << source._startPositionsD << std::endl;
    
    output << "_dimSubMatrixA: " << source._dimSubMatrixA << std::endl;
    
    output << "_dimSubMatrixB: " << source._dimSubMatrixB << std::endl;
    
    output << "_dimSubMatrixC: " << source._dimSubMatrixC << std::endl;
    
    output << "_dimSubMatrixD: " << source._dimSubMatrixD << std::endl;
    
    return output;
}
