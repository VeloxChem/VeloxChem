#include "Matrix.hpp"

#include <algorithm>
#include <cmath>
#include <ranges>

#include "CustomViews.hpp"

CMatrix::CMatrix()

    : _sub_matrices(std::map<std::pair<int, int>, CSubMatrix *>())

    , _mat_type(mat_t::general)
{
}

CMatrix::CMatrix(const std::map<std::pair<int, int>, CSubMatrix> &sub_matrices, const mat_t mat_type)

    : _sub_matrices(std::map<std::pair<int, int>, CSubMatrix *>())

    , _mat_type(mat_type)
{
    std::ranges::for_each(sub_matrices, [&](const auto &mvalue) {
        _sub_matrices.insert({mvalue.first, new CSubMatrix(mvalue.second)});
    });
}

CMatrix::CMatrix(const std::vector<std::pair<int, int>> &ang_pairs,
                 const std::vector<CSubMatrix>          &sub_matrices,
                 const mat_t                             mat_type)

    : _sub_matrices(std::map<std::pair<int, int>, CSubMatrix *>())

    , _mat_type(mat_type)
{
    if (const auto ndims = ang_pairs.size(); ndims == sub_matrices.size())
    {
        std::ranges::for_each(std::views::iota(size_t{0}, ndims), [&](const size_t i) {
            _sub_matrices.insert({ang_pairs[i], new CSubMatrix(sub_matrices[i])});
        });
    }
}

CMatrix::CMatrix(const CMatrix &other)

    : _sub_matrices(std::map<std::pair<int, int>, CSubMatrix *>())

    , _mat_type(other._mat_type)
{
    std::ranges::for_each(other._sub_matrices, [&](const auto &mvalue) {
        _sub_matrices.insert({mvalue.first, new CSubMatrix(*mvalue.second)});
    });
}

CMatrix::CMatrix(CMatrix &&other) noexcept

    : _sub_matrices(std::map<std::pair<int, int>, CSubMatrix *>())

    , _mat_type(mat_t::general)
{
    std::swap(_sub_matrices, other._sub_matrices);

    std::swap(_mat_type, other._mat_type);
}

CMatrix::~CMatrix()
{
    _deallocate();
}

auto
CMatrix::operator=(const CMatrix &other) -> CMatrix &
{
    _deallocate();

    _mat_type = other._mat_type;

    std::ranges::for_each(other._sub_matrices, [&](const auto &mvalue) {
        _sub_matrices.insert({mvalue.first, new CSubMatrix(*mvalue.second)});
    });

    return *this;
}

auto
CMatrix::operator=(CMatrix &&other) noexcept -> CMatrix &
{
    std::swap(_sub_matrices, other._sub_matrices);

    std::swap(_mat_type, other._mat_type);

    return *this;
}

auto
CMatrix::operator==(const CMatrix &other) const -> bool
{
    if (_mat_type != other._mat_type)
    {
        return false;
    }
    else if (angular_pairs() != other.angular_pairs())
    {
        return false;
    }
    else
    {
        auto keys = angular_pairs();

        auto pos = std::ranges::find_if(
            keys, [&](const auto &key) { return *_sub_matrices.at(key) != *other._sub_matrices.at(key); });

        return keys.end() == pos;
    }
}

auto
CMatrix::operator!=(const CMatrix &other) const -> bool
{
    return !(*this == other);
}

auto
CMatrix::operator+(const CMatrix &other) const -> CMatrix
{
    CMatrix mat;

    mat.set_type(_mat_type);

    if (const auto keys = angular_pairs(); keys == other.angular_pairs())
    {
        std::ranges::for_each(keys, [&](const auto &key) {
            mat._sub_matrices.insert({key, new CSubMatrix(*_sub_matrices.at(key) + *other._sub_matrices.at(key))});
        });
    }

    return mat;
}

auto
CMatrix::add(const CSubMatrix &sub_matrix, const std::pair<int, int> &angpair) -> void
{
    _sub_matrices.insert({angpair, new CSubMatrix(sub_matrix)});
}

auto
CMatrix::add(const std::array<size_t, 4> &dimensions, const std::pair<int, int> &angpair) -> void
{
    _sub_matrices.insert({angpair, new CSubMatrix(dimensions)});
}

auto
CMatrix::set_type(const mat_t mat_type) -> void
{
    _mat_type = mat_type;
}

auto
CMatrix::zero() -> void
{
    for (auto &mvalue : _sub_matrices)
        mvalue.second->zero();
}

auto
CMatrix::scale(const double factor) -> void
{
    for (auto &mvalue : _sub_matrices)
        mvalue.second->scale(factor);
}

auto
CMatrix::symmetrize() -> void
{
    if (_mat_type == mat_t::symmetric)
    {
        std::ranges::for_each(_sub_matrices, [&](auto &mvalue) {
            if (mvalue.first.first == mvalue.first.second) mvalue.second->symmetrize();
        });
    }
}

auto
CMatrix::angular_pairs() const -> std::vector<std::pair<int, int>>
{
    std::vector<std::pair<int, int>> angpairs;

    angpairs.reserve(_sub_matrices.size());

    std::ranges::for_each(_sub_matrices, [&](auto &mvalue) { angpairs.push_back(mvalue.first); });

    return angpairs;
}

auto
CMatrix::get_type() const -> mat_t
{
    return _mat_type;
}

auto
CMatrix::sub_matrices() const -> std::vector<CSubMatrix>
{
    std::vector<CSubMatrix> sub_mats;

    sub_mats.reserve(_sub_matrices.size());

    std::ranges::for_each(_sub_matrices, [&](auto &mvalue) { sub_mats.push_back(*mvalue.second); });

    return sub_mats;
}

auto
CMatrix::sub_matrix(const std::pair<int, int> &angpair) -> CSubMatrix *
{
    auto pos = std::ranges::find_if(_sub_matrices, [&](auto &mvalue) { return mvalue.first == angpair; });

    if ((pos == _sub_matrices.end()) && ((_mat_type == mat_t::symmetric) || (_mat_type == mat_t::antisymmetric)))
    {
        const auto rev_angpair = std::pair<int, int>({angpair.second, angpair.first});

        pos = std::ranges::find_if(_sub_matrices, [&](auto &mvalue) { return mvalue.first == rev_angpair; });
    }

    return (pos == _sub_matrices.end()) ? nullptr : pos->second;
}

auto
CMatrix::sub_matrix(const std::pair<int, int> &angpair) const -> const CSubMatrix *
{
    auto pos = std::ranges::find_if(_sub_matrices, [&](auto &mvalue) { return mvalue.first == angpair; });

    if ((pos == _sub_matrices.end()) && ((_mat_type == mat_t::symmetric) || (_mat_type == mat_t::antisymmetric)))
    {
        const auto rev_angpair = std::pair<int, int>({angpair.second, angpair.first});

        pos = std::ranges::find_if(_sub_matrices, [&](auto &mvalue) { return mvalue.first == rev_angpair; });
    }

    return (pos == _sub_matrices.end()) ? nullptr : pos->second;
}

auto
CMatrix::is_angular_order(const std::pair<int, int> &angpair) const -> bool
{
    auto pos = std::ranges::find_if(_sub_matrices, [&](auto &mvalue) { return mvalue.first == angpair; });

    return pos != _sub_matrices.end();
}

auto
CMatrix::number_of_rows() const -> size_t
{
    const auto row_keys = _row_angular_keys();

    const auto col_keys = _column_angular_keys();

    if (row_keys.empty() || col_keys.empty())
    {
        return 0;
    }
    else
    {
        const auto col_ang = *(col_keys.cbegin());

        size_t nrows = 0;

        std::ranges::for_each(row_keys, [&](const auto key) {
            if (const auto submat = sub_matrix({key, col_ang}); submat != nullptr)
            {
                if (is_angular_order({key, col_ang}))
                {
                    nrows += submat->number_of_rows();
                }
                else
                {
                    nrows += submat->number_of_columns();
                }
            }
        });

        return nrows;
    }
}

auto
CMatrix::number_of_columns() const -> size_t
{
    const auto row_keys = _row_angular_keys();

    const auto col_keys = _column_angular_keys();

    if (row_keys.empty() || col_keys.empty())
    {
        return 0;
    }
    else
    {
        const auto row_ang = *(row_keys.cbegin());

        size_t ncols = 0;

        std::ranges::for_each(col_keys, [&](const auto key) {
            if (const auto submat = sub_matrix({row_ang, key}); submat != nullptr)
            {
                if (is_angular_order({row_ang, key}))
                {
                    ncols += submat->number_of_columns();
                }
                else
                {
                    ncols += submat->number_of_rows();
                }
            }
        });

        return ncols;
    }
}

auto
CMatrix::full_matrix() const -> CSubMatrix
{
    auto matrix = CSubMatrix({0, 0, number_of_rows(), number_of_columns()});

    std::ranges::for_each(_sub_matrices, [&](auto &mvalue) {
        const auto submat = mvalue.second;
        const auto dims   = submat->get_dimensions();
        std::ranges::for_each(views::rectangular(dims[2], dims[3]), [&](const auto &index) {
            matrix.at({index.first + dims[0], index.second + dims[1]}) = submat->at(index);
            if (mvalue.first.first != mvalue.first.second)
            {
                if (_mat_type == mat_t::symmetric)
                {
                    matrix.at({index.second + dims[1], index.first + dims[0]}) = submat->at(index);
                }

                if (_mat_type == mat_t::antisymmetric)
                {
                    matrix.at({index.second + dims[1], index.first + dims[0]}) = -submat->at(index);
                }
            }
        });
    });

    return matrix;
}

auto
CMatrix::pointer() const -> const CMatrix *
{
    return this;
}

auto
CMatrix::pointer() -> CMatrix *
{
    return this;
}

auto
CMatrix::_row_angular_keys() const -> std::set<int>
{
    std::set<int> row_keys;

    std::ranges::for_each(_sub_matrices, [&](auto &mvalue) { row_keys.insert(mvalue.first.first); });

    return row_keys;
}

auto
CMatrix::_column_angular_keys() const -> std::set<int>
{
    std::set<int> col_keys;

    std::ranges::for_each(_sub_matrices, [&](auto &mvalue) { col_keys.insert(mvalue.first.second); });

    return col_keys;
}

auto
CMatrix::_deallocate() -> void
{
    std::ranges::for_each(_sub_matrices, [&](auto &mvalue) { delete mvalue.second; });

    _sub_matrices.clear();
}
