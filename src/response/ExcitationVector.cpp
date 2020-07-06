//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ExcitationVector.hpp"

#include <set>
#include <sstream>

#include "DenseLinearAlgebra.hpp"
#include "StringFormat.hpp"

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
                                     const int32_t ketEndPosition,
                                     const bool    isTammDancoff)

    : _excitationType(excitationType)
{
    auto nbra = braEndPosition - braStartPosition;

    auto nket = ketEndPosition - ketStartPosition;

    _braIndexes = CMemBlock<int32_t>(nket * nbra);

    _ketIndexes = CMemBlock<int32_t>(nket * nbra);

    _zCoefficents = CMemBlock<double>(nket * nbra);

    _zCoefficents.zero();

    if (!isTammDancoff)
    {
        _yCoefficents = CMemBlock<double>(nket * nbra);

        _yCoefficents.zero();
    }

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

CExcitationVector::CExcitationVector(const szblock excitationType,
                                     const int32_t braStartPosition,
                                     const int32_t braEndPosition,
                                     const int32_t ketStartPosition,
                                     const int32_t ketEndPosition)

    : CExcitationVector(excitationType, braStartPosition, braEndPosition, ketStartPosition, ketEndPosition, false)
{
}

CExcitationVector::CExcitationVector(const std::vector<double>& factors, const std::vector<CExcitationVector>& sources)
{
    auto nvecs = static_cast<int32_t>(factors.size());

    if (nvecs > 0)
    {
        // copy first excitation vector

        _excitationType = sources[0]._excitationType;

        _braIndexes = sources[0]._braIndexes;

        _ketIndexes = sources[0]._ketIndexes;

        _zCoefficents = sources[0]._zCoefficents;

        _yCoefficents = sources[0]._yCoefficents;

        // scale first excitation vector

        mathfunc::scale(_zCoefficents.data(), factors[0], _zCoefficents.size());

        mathfunc::scale(_yCoefficents.data(), factors[0], _yCoefficents.size());

        // add remaining scaled excitation vectors

        for (int32_t i = 1; i < nvecs; i++)
        {
            mathfunc::add_scaled(_zCoefficents.data(), sources[i]._zCoefficents.data(), factors[i], _zCoefficents.size());

            mathfunc::add_scaled(_yCoefficents.data(), sources[i]._yCoefficents.data(), factors[i], _yCoefficents.size());
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
CExcitationVector::setCoefficientsZY(const CMemBlock<double>& zCoefficients, const CMemBlock<double>& yCoefficients)
{
    _zCoefficents = zCoefficients;

    _yCoefficents = yCoefficients;
}

void
CExcitationVector::setCoefficientZ(const double zValue, const int32_t iCoefficient)
{
    _zCoefficents.at(iCoefficient) = zValue;
}

void
CExcitationVector::setCoefficientY(const double yValue, const int32_t iCoefficient)
{
    _yCoefficents.at(iCoefficient) = yValue;
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

double
CExcitationVector::dotCoefficientsZ(const CExcitationVector& other) const
{
    return denblas::dot(_zCoefficents, other._zCoefficents);
}

double
CExcitationVector::dotCoefficientsY(const CExcitationVector& other) const
{
    return denblas::dot(_yCoefficents, other._yCoefficents);
}

int32_t
CExcitationVector::getNumberOfExcitations() const
{
    return _zCoefficents.size();
}

const int32_t*
CExcitationVector::getBraIndexes() const
{
    return _braIndexes.data();
}

const int32_t*
CExcitationVector::getKetIndexes() const
{
    return _ketIndexes.data();
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
    // check empty vector

    if (_zCoefficents.size() == 0)
    {
        return CDenseMatrix();
    }

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
    // check empty vector

    if (_yCoefficents.size() == 0)
    {
        return CDenseMatrix();
    }

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

    auto tmo = _getBraOrbitals(molecularOrbitals, getBraUniqueIndexes());

    auto tmv = _getKetOrbitals(molecularOrbitals, getKetUniqueIndexes());

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

const double*
CExcitationVector::getBraEnergies(const CMolecularOrbitals& molecularOrbitals) const
{
    if ((_excitationType == szblock::aa) || (_excitationType == szblock::ab))
    {
        return molecularOrbitals.alphaEnergies();
    }

    if ((_excitationType == szblock::bb) || (_excitationType == szblock::ba))
    {
        return molecularOrbitals.betaEnergies();
    }

    return nullptr;
}

const double*
CExcitationVector::getKetEnergies(const CMolecularOrbitals& molecularOrbitals) const
{
    if ((_excitationType == szblock::aa) || (_excitationType == szblock::ba))
    {
        return molecularOrbitals.alphaEnergies();
    }

    if ((_excitationType == szblock::bb) || (_excitationType == szblock::ab))
    {
        return molecularOrbitals.betaEnergies();
    }

    return nullptr;
}

std::vector<int32_t>
CExcitationVector::getSmallEnergyIdentifiers(const CMolecularOrbitals& molecularOrbitals, const int32_t nExcitations) const
{
    auto ntot = getNumberOfExcitations();

    std::vector<int32_t> idxvec;

    if (ntot > nExcitations)
    {
        auto beigs = getBraEnergies(molecularOrbitals);

        auto keigs = getKetEnergies(molecularOrbitals);

        std::vector<double> exvals;

        for (int32_t i = 0; i < nExcitations; i++)
        {
            idxvec.push_back(i);

            exvals.push_back(keigs[_ketIndexes.at(i)] - beigs[_braIndexes.at(i)]);
        }

        auto tmax = _getMaxElement(exvals);

        for (int32_t i = nExcitations; i < ntot; i++)
        {
            auto ceval = keigs[_ketIndexes.at(i)] - beigs[_braIndexes.at(i)];

            if (std::get<1>(tmax) > ceval)
            {
                auto midx = std::get<0>(tmax);

                exvals[midx] = ceval;

                idxvec[midx] = i;

                tmax = _getMaxElement(exvals);
            }
        }
    }
    else
    {
        for (int32_t i = 0; i < ntot; i++)
        {
            idxvec.push_back(i);
        }
    }

    return idxvec;
}

CMemBlock<double>
CExcitationVector::getApproximateDiagonal(const CMolecularOrbitals& molecularOrbitals) const
{
    auto ndim = getNumberOfExcitations();

    CMemBlock<double> diagmat(ndim);

    auto beigs = getBraEnergies(molecularOrbitals);

    auto keigs = getKetEnergies(molecularOrbitals);

    auto bidx = _braIndexes.data();

    auto kidx = _ketIndexes.data();

    for (int32_t i = 0; i < ndim; i++)
    {
        diagmat.at(i) = keigs[kidx[i]] - beigs[bidx[i]];
    }

    return diagmat;
}

std::string
CExcitationVector::getString() const
{
    std::stringstream sst("");

    sst << "Z Vector [Dimension " << _zCoefficents.size() << "]\n";

    for (int32_t i = 0; i < _zCoefficents.size(); i++)
    {
        sst << fstr::to_string(_braIndexes.at(i), 6, fmt::right);

        sst << " -> ";

        sst << fstr::to_string(_ketIndexes.at(i), 6, fmt::right);

        sst << " : ";

        sst << fstr::to_string(_zCoefficents.at(i), 8, 15, fmt::right);

        sst << "\n";
    }

    sst << "Y Vector [Dimension " << _yCoefficents.size() << "]\n";

    for (int32_t i = 0; i < _yCoefficents.size(); i++)
    {
        sst << fstr::to_string(_ketIndexes.at(i), 6, fmt::right);

        sst << " -> ";

        sst << fstr::to_string(_braIndexes.at(i), 6, fmt::right);

        sst << " : ";

        sst << fstr::to_string(_yCoefficents.at(i), 8, 15, fmt::right);

        sst << "\n";
    }

    return sst.str();
}

int32_t
CExcitationVector::_convertIdentifier(const std::vector<int32_t>& identifiers, const int32_t index) const
{
    auto dim = static_cast<int32_t>(identifiers.size());

    for (int32_t i = 0; i < dim; i++)
    {
        if (identifiers[i] == index) return i;
    }

    return -1;
}

CDenseMatrix
CExcitationVector::_getBraOrbitals(const CMolecularOrbitals& molecularOrbitals, const std::vector<int32_t>& identifiers) const
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
CExcitationVector::_getKetOrbitals(const CMolecularOrbitals& molecularOrbitals, const std::vector<int32_t>& identifiers) const
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

std::tuple<int32_t, double>
CExcitationVector::_getMaxElement(const std::vector<double>& vectorA) const
{
    int32_t midx = 0;

    double mval = vectorA[0];

    for (int32_t i = 1; i < static_cast<int32_t>(vectorA.size()); i++)
    {
        auto cval = vectorA[i];

        if (mval < cval)
        {
            mval = cval;

            midx = i;
        }
    }

    return std::make_tuple(midx, mval);
}

std::ostream&
operator<<(std::ostream& output, const CExcitationVector& source)
{
    output << std::endl;

    output << "[CExcitationVector (Object):" << &source << "]" << std::endl;

    output << "_excitationType: " << to_string(source._excitationType) << std::endl;

    output << "_braIndexes: " << source._braIndexes << std::endl;

    output << "_ketIndexes: " << source._ketIndexes << std::endl;

    output << "_zCoefficents: " << source._zCoefficents << std::endl;

    output << "_yCoefficents: " << source._yCoefficents << std::endl;

    return output;
}
