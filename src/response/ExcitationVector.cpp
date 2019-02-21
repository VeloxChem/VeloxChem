//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ExcitationVector.hpp"

#include <set>

#include "DenseLinearAlgebra.hpp"

CExcitationVector::CExcitationVector()

    : _excitationType(szblock::aa)
{
    
}

CExcitationVector::CExcitationVector(const szblock               excitationType,
                                     const std::vector<int32_t>& braIndexes,
                                     const std::vector<int32_t>& ketIndexes,
                                     const std::vector<double>&  zCoefficients,
                                     const std::vector<double>&  yCoefficients)

    : _excitationType(excitationType)

    , _braIndexes(CMemBlock<int32_t>(braIndexes))

    , _ketIndexes(CMemBlock<int32_t>(ketIndexes))

    , _zCoefficents(CMemBlock<double>(zCoefficients))

    , _yCoefficents(CMemBlock<double>(yCoefficients))
{
    
}

CExcitationVector::CExcitationVector(const szblock excitationType,
                                     const int32_t braStartPosition,
                                     const int32_t braEndPosition,
                                     const int32_t ketStartPosition,
                                     const int32_t ketEndPosition)

    : _excitationType(excitationType)
{
    auto nbra = braEndPosition - braStartPosition;
    
    auto nket = ketEndPosition - ketStartPosition;
    
    _braIndexes = CMemBlock<int32_t>(nket * nbra);
    
    _ketIndexes = CMemBlock<int32_t>(nket * nbra);
    
    _zCoefficents = CMemBlock<double>(nket * nbra);
    
    _yCoefficents = CMemBlock<double>(nket * nbra);
    
    _zCoefficents.zero();
    
    _yCoefficents.zero();
    
    int32_t idx = 0;
    
    for (int32_t i = braStartPosition; i < braEndPosition; i++)
    {
        for (int32_t j = ketStartPosition; j < ketEndPosition; j++)
        {
            _braIndexes.at(idx) = i;
            
            _ketIndexes.at(idx) = j; 
            
            idx++;
        }
    }
}

CExcitationVector::CExcitationVector(const CExcitationVector& source)

    : _excitationType(source._excitationType)

    , _braIndexes(source._braIndexes)

    , _ketIndexes(source._ketIndexes)

    , _zCoefficents(source._zCoefficents)

    , _yCoefficents(source._yCoefficents)
{
    
}

CExcitationVector::CExcitationVector(CExcitationVector&& source) noexcept

    : _excitationType(std::move(source._excitationType))

    , _braIndexes(std::move(source._braIndexes))

    , _ketIndexes(std::move(source._ketIndexes))

    , _zCoefficents(std::move(source._zCoefficents))

    , _yCoefficents(std::move(source._yCoefficents))
{
    
}

CExcitationVector::~CExcitationVector()
{
    
}

CExcitationVector&
CExcitationVector::operator=(const CExcitationVector& source)
{
    if (this == &source) return *this;
    
    _excitationType = source._excitationType;
    
    _braIndexes = source._braIndexes;
    
    _ketIndexes = source._ketIndexes;
    
    _zCoefficents = source._zCoefficents;
    
    _yCoefficents = source._yCoefficents;
    
    return *this;
}

CExcitationVector&
CExcitationVector::operator=(CExcitationVector&& source) noexcept
{
    if (this == &source) return *this;
    
    _excitationType = std::move(source._excitationType);
    
    _braIndexes = std::move(source._braIndexes);
    
    _ketIndexes = std::move(source._ketIndexes);
    
    _zCoefficents = std::move(source._zCoefficents);
    
    _yCoefficents = std::move(source._yCoefficents);
    
    return *this;
}

bool
CExcitationVector::operator==(const CExcitationVector& other) const
{
    if (_excitationType != other._excitationType) return false;
    
    if (_braIndexes != other._braIndexes) return false;
    
    if (_ketIndexes != other._ketIndexes) return false;
    
    if (_zCoefficents != other._zCoefficents) return false;
    
    if (_yCoefficents != other._yCoefficents) return false;
    
    return true;
}

bool
CExcitationVector::operator!=(const CExcitationVector& other) const
{
    return !(*this == other);
}

void
CExcitationVector::setCoefficientsZY(const CMemBlock<double>& zCoefficients,
                                     const CMemBlock<double>& yCoefficients)
{
    _zCoefficents = zCoefficients;
    
    _yCoefficents = yCoefficients; 
}

double*
CExcitationVector::getCoefficientsZ()
{
    return _zCoefficents.data();
}

const double*
CExcitationVector::getCoefficientsZ() const
{
    return _zCoefficents.data();
}

double*
CExcitationVector::getCoefficientsY()
{
    return _yCoefficents.data();
}

const double*
CExcitationVector::getCoefficientsY() const
{
    return _yCoefficents.data();
}

int32_t
CExcitationVector::getNumberOfExcitations() const
{
    return _zCoefficents.size();
}

std::vector<int32_t>
CExcitationVector::getBraUniqueIndexes() const
{
    std::set<int32_t> exlist;
    
    for (int32_t i = 0; i < _braIndexes.size(); i++)
    {
        exlist.insert(_braIndexes.at(i));
    }
    
    return std::vector<int32_t>(exlist.cbegin(), exlist.cend());
}

std::vector<int32_t>
CExcitationVector::getKetUniqueIndexes() const
{
    std::set<int32_t> exlist;
    
    for (int32_t i = 0; i < _ketIndexes.size(); i++)
    {
        exlist.insert(_ketIndexes.at(i));
    }
    
    return std::vector<int32_t>(exlist.cbegin(), exlist.cend());
}

CDenseMatrix
CExcitationVector::getMatrixZ() const
{
    // get unique indexes vectors
    
    auto bidx = getBraUniqueIndexes();
    
    auto kidx = getKetUniqueIndexes();
    
    //  determine dimension
    
    auto nrow = static_cast<int32_t>(bidx.size());
    
    auto ncol = static_cast<int32_t>(kidx.size());
    
    // allocate and initialize matrix
    
    CDenseMatrix zmat(nrow, ncol);
    
    zmat.zero();
    
    auto zvals = zmat.values();
    
    for (int32_t i = 0; i < getNumberOfExcitations(); i++)
    {
        auto irow = _convertIdentifier(bidx, _braIndexes.at(i));
        
        auto icol = _convertIdentifier(kidx, _ketIndexes.at(i));
        
        zvals[irow * ncol + icol] = _zCoefficents.at(i); 
    }
    
    return zmat;
}

CDenseMatrix
CExcitationVector::getMatrixY() const
{
    // get unique indexes vectors
    
    auto bidx = getBraUniqueIndexes();
    
    auto kidx = getKetUniqueIndexes();
    
    //  determine dimension
    
    auto nrow = static_cast<int32_t>(kidx.size());
    
    auto ncol = static_cast<int32_t>(bidx.size());
    
    // allocate and initialize matrix
    
    CDenseMatrix zmat(nrow, ncol);
    
    zmat.zero();
    
    auto zvals = zmat.values();
    
    for (int32_t i = 0; i < getNumberOfExcitations(); i++)
    {
        auto irow = _convertIdentifier(kidx, _ketIndexes.at(i));
        
        auto icol = _convertIdentifier(bidx, _braIndexes.at(i));
        
        zvals[irow * ncol + icol] = _yCoefficents.at(i);
    }
    
    return zmat;
}

CAODensityMatrix
CExcitationVector::getDensityZ(const CMolecularOrbitals& molecularOrbitals) const
{
    // convert Z vector to matrix
    
    auto zmat = getMatrixZ();
    
    // set up appropiate molecular orbitals blocks
    
    auto bidx = getBraUniqueIndexes();
    
    std::cout << "bra:";
    for (size_t i = 0; i < bidx.size(); i++)
        std::cout << " " << bidx[i];
    std::cout << std::endl;
    
    auto tmo = _getBraOrbitals(molecularOrbitals, getBraUniqueIndexes());
    
    auto tmv = _getKetOrbitals(molecularOrbitals, getKetUniqueIndexes());
    
    std::cout << "tmo:" << tmo.getString();
    
    std::cout << "tmv:" << tmv.getString();
    
    // compute transformed density
    
    auto tden = denblas::multAB(tmo, denblas::multABt(zmat, tmv));

    return CAODensityMatrix({tden}, denmat::rgen);
}

CAODensityMatrix
CExcitationVector::getDensityY(const CMolecularOrbitals& molecularOrbitals) const
{
    // convert Y vector to matrix
    
    auto ymat = getMatrixY();
    
    // set up appropiate molecular orbitals blocks
    
    auto tmo = _getBraOrbitals(molecularOrbitals, getBraUniqueIndexes());
    
    auto tmv = _getBraOrbitals(molecularOrbitals, getKetUniqueIndexes());
    
    // compute transformed density
    
    auto tden = denblas::multAB(tmv, denblas::multABt(ymat, tmo));
    
    return CAODensityMatrix({tden}, denmat::rgen);
}

int32_t
CExcitationVector::_convertIdentifier(const std::vector<int32_t>& identifiers,
                                      const int32_t               index) const
{
    auto dim = static_cast<int32_t>(identifiers.size());
    
    for (int32_t i = 0; i < dim; i++)
    {
        if (identifiers[i] == index) return i;
    }
    
    return -1;
}

CDenseMatrix
CExcitationVector::_getBraOrbitals(const CMolecularOrbitals&   molecularOrbitals,
                                   const std::vector<int32_t>& identifiers) const
{
    if ((_excitationType == szblock::aa) || (_excitationType == szblock::ab))
    {
        return molecularOrbitals.alphaOrbitals(identifiers);
    }
    
    if ((_excitationType == szblock::bb) || (_excitationType == szblock::ba))
    {
        return molecularOrbitals.betaOrbitals(identifiers);
    }
    
    return CDenseMatrix();
}

CDenseMatrix
CExcitationVector::_getKetOrbitals(const CMolecularOrbitals&   molecularOrbitals,
                                   const std::vector<int32_t>& identifiers) const
{
    if ((_excitationType == szblock::aa) || (_excitationType == szblock::ba))
    {
        return molecularOrbitals.alphaOrbitals(identifiers);
    }
    
    if ((_excitationType == szblock::bb) || (_excitationType == szblock::ab))
    {
        return molecularOrbitals.betaOrbitals(identifiers);
    }
    
    return CDenseMatrix();
}

std::ostream&
operator<<(      std::ostream&  output,
           const CExcitationVector& source)
{
    output << std::endl;
    
    output << "[CExcitationVector (Object):" << &source << "]" << std::endl;
    
    output << "_excitationType: " << to_string(source._excitationType) << std::endl;
    
    output << "_braIndexes: " << source._braIndexes << std::endl;
    
    output << "_ketIndexes: " << source._ketIndexes <<  std::endl;
    
    output << "_zCoefficents: " << source._zCoefficents <<  std::endl;
    
    output << "_yCoefficents: " << source._yCoefficents <<  std::endl;
    
    return output;
}
