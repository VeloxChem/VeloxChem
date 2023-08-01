#include "FockMatrix.hpp"

CFockMatrix::CFockMatrix()

    : _matrix(CMatrix())

    , _exc_scale(0.0)

    , _ftype(fock_t::restjk)
{
    
}

CFockMatrix::CFockMatrix(const CMatrix& matrix,
                         const double   exc_scale,
                         const fock_t   ftype)

    : _matrix(matrix)

    , _exc_scale(exc_scale)

    , _ftype(ftype)
{
    
}

CFockMatrix::CFockMatrix(const CFockMatrix& other)
{
    _matrix = other._matrix;
    
    _exc_scale = other._exc_scale;
    
    _ftype = other._ftype;
}

auto
CFockMatrix::setExchangeScale(const double exc_scale) -> void
{
    _exc_scale = exc_scale;
}

auto
CFockMatrix::setFockType(const fock_t ftype) -> void
{
    _ftype = ftype;
}

auto
CFockMatrix::zero() -> void
{
    _matrix.zero();
}

auto
CFockMatrix::getExchangeScale() const -> double
{
    return _exc_scale;
}

auto
CFockMatrix::getFockType() const -> fock_t
{
    return _ftype;
}

auto
CFockMatrix::getAngularPairs() const -> std::vector<T2Pair>
{
    return _matrix.getAngularPairs();
}

auto
CFockMatrix::getStorageType() const -> mat_t
{
    return _matrix.getType();
}

auto
CFockMatrix::getSubMatrix(const T2Pair& angpair) -> CSubMatrix*
{
    return _matrix.getSubMatrix(angpair);
}

auto
CFockMatrix::getSubMatrix(const T2Pair& angpair) const -> const CSubMatrix*
{
    return _matrix.getSubMatrix(angpair);
}

auto
CFockMatrix::isAngularOrder(const T2Pair& angpair) const -> bool
{
    return _matrix.isAngularOrder(angpair);
}

auto
CFockMatrix::getNumberOfRows() const -> int64_t
{
    return _matrix.getNumberOfRows();
}

auto
CFockMatrix::getNumberOfColumns() const -> int64_t
{
    return _matrix.getNumberOfColumns();
}

auto
CFockMatrix::getFullMatrix() const -> CSubMatrix
{
    return _matrix.getFullMatrix();
}
