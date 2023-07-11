#include "Matrices.hpp"

#include "StringFormat.hpp"
#include "CantorFunc.hpp"

CMatrices::CMatrices()

    : _matrices(std::map<int64_t, CMatrix*>())
{
}

CMatrices::CMatrices(const std::map<int64_t, CMatrix>& matrices)

    : _matrices(std::map<int64_t, CMatrix*>())
{
    for (const auto& mvalue : matrices)
    {
        _matrices.insert({mvalue.first, new CMatrix(mvalue.second)});
    }
}

CMatrices::CMatrices(const CMatrix& matrix, const std::vector<int64_t>& keys)

    : _matrices(std::map<int64_t, CMatrix*>())
{
    for (const auto key : keys)
    {
        _matrices.insert({key, new CMatrix(matrix)});
    }
}

CMatrices::CMatrices(const CMatrices& other)

    : _matrices(std::map<int64_t, CMatrix*>())
{
    for (const auto& mvalue : other._matrices)
    {
        _matrices.insert({mvalue.first, new CMatrix(*mvalue.second)});
    }
}

CMatrices::~CMatrices()
{
    for (auto& mvalue : _matrices)
    {
        delete mvalue.second;
    }
}

auto
CMatrices::add(const CMatrix& matrix, const int64_t key) -> void
{
    _matrices.insert({key, new CMatrix(matrix)});
}

auto
CMatrices::add(const CMatrix& matrix, const std::string& label) -> void
{
    _matrices.insert({_to_key(label), new CMatrix(matrix)});
}

auto
CMatrices::add(const CMatrix& matrix, const int64_t atom, const std::string& label) -> void
{
    const auto key = mathfunc::getCantorIndex({atom, _to_key(label)});
    
    _matrices.insert({key, new CMatrix(matrix)});
}

auto
CMatrices::add(const CMatrix& matrix, const int64_t atom_a, const std::string& label_a, const int64_t atom_b, const std::string& label_b) -> void
{
    const auto key_a = mathfunc::getCantorIndex({atom_a, _to_key(label_a)});
    
    const auto key_b = mathfunc::getCantorIndex({atom_b, _to_key(label_b)});
    
    _matrices.insert({mathfunc::getCantorIndex({key_a, key_b}), new CMatrix(matrix)});
}

auto
CMatrices::add(const CMatrix& matrix, const int64_t atom_a, const std::string& label_a, const int64_t atom_b, const std::string& label_b,
               const int64_t atom_c, const std::string& label_c) -> void
{
    const auto key_a = mathfunc::getCantorIndex({atom_a, _to_key(label_a)});

    const auto key_b = mathfunc::getCantorIndex({atom_b, _to_key(label_b)});
    
    const auto key_c = mathfunc::getCantorIndex({atom_c, _to_key(label_c)});

    const auto key_ab = mathfunc::getCantorIndex({key_a, key_b});
    
    _matrices.insert({mathfunc::getCantorIndex({key_ab, key_c}), new CMatrix(matrix)});
}

auto
CMatrices::add(const CMatrix& matrix, const int64_t atom_a, const std::string& label_a, const int64_t atom_b, const std::string& label_b,
               const int64_t atom_c, const std::string& label_c, const int64_t atom_d, const std::string& label_d) -> void
{
    const auto key_a = mathfunc::getCantorIndex({atom_a, _to_key(label_a)});

    const auto key_b = mathfunc::getCantorIndex({atom_b, _to_key(label_b)});
    
    const auto key_c = mathfunc::getCantorIndex({atom_c, _to_key(label_c)});
    
    const auto key_d = mathfunc::getCantorIndex({atom_d, _to_key(label_d)});

    const auto key_ab = mathfunc::getCantorIndex({key_a, key_b});
    
    const auto key_cd = mathfunc::getCantorIndex({key_c, key_d});
    
    _matrices.insert({mathfunc::getCantorIndex({key_ab, key_cd}), new CMatrix(matrix)});
}

auto
CMatrices::zero() -> void
{
    for (auto& mvalue : _matrices)
    {
        mvalue.second->zero();
    }
}

auto
CMatrices::getKeys() const -> std::vector<int64_t>
{
    if (!_matrices.empty())
    {
        std::vector<int64_t> keys;

        for (const auto& mvalue : _matrices)
        {
            keys.push_back(mvalue.first);
        }

        return keys;
    }
    else
    {
        return std::vector<int64_t>();
    }
}

auto
CMatrices::getMatrix(const int64_t key) -> CMatrix*
{
    for (const auto& mvalue : _matrices)
    {
        if (mvalue.first == key)
        {
            return mvalue.second;
        }
    }

    return nullptr;
}

auto
CMatrices::getMatrix(const int64_t key) const -> const CMatrix*
{
    for (const auto& mvalue : _matrices)
    {
        if (mvalue.first == key)
        {
            return mvalue.second;
        }
    }

    return nullptr;
}

auto
CMatrices::getMatrix(const std::string& label) -> CMatrix*
{
    return getMatrix(_to_key(label));
}

auto
CMatrices::getMatrix(const std::string& label) const -> const CMatrix*
{
    return getMatrix(_to_key(label));
}

auto
CMatrices::getMatrix(const int64_t atom, const std::string& label) -> CMatrix*
{
    const auto key = mathfunc::getCantorIndex({atom, _to_key(label)});
    
    return getMatrix(key);
}

auto
CMatrices::getMatrix(const int64_t atom, const std::string& label) const -> const CMatrix*
{
    const auto key = mathfunc::getCantorIndex({atom, _to_key(label)});
    
    return getMatrix(key);
}

auto
CMatrices::getMatrix(const int64_t atom_a, const std::string& label_a, const int64_t atom_b, const std::string& label_b) -> CMatrix*
{
    const auto key_a = mathfunc::getCantorIndex({atom_a, _to_key(label_a)});

    const auto key_b = mathfunc::getCantorIndex({atom_b, _to_key(label_b)});
    
    return getMatrix(mathfunc::getCantorIndex({key_a, key_b}));
}

auto
CMatrices::getMatrix(const int64_t atom_a, const std::string& label_a, const int64_t atom_b, const std::string& label_b) const -> const CMatrix*
{
    const auto key_a = mathfunc::getCantorIndex({atom_a, _to_key(label_a)});

    const auto key_b = mathfunc::getCantorIndex({atom_b, _to_key(label_b)});
    
    return getMatrix(mathfunc::getCantorIndex({key_a, key_b}));
}

auto
CMatrices::getMatrix(const int64_t atom_a, const std::string& label_a, const int64_t atom_b, const std::string& label_b,
               const int64_t atom_c, const std::string& label_c) -> CMatrix*
{
    const auto key_a = mathfunc::getCantorIndex({atom_a, _to_key(label_a)});

    const auto key_b = mathfunc::getCantorIndex({atom_b, _to_key(label_b)});
    
    const auto key_c = mathfunc::getCantorIndex({atom_c, _to_key(label_c)});

    const auto key_ab = mathfunc::getCantorIndex({key_a, key_b});
    
    return getMatrix(mathfunc::getCantorIndex({key_ab, key_c}));
}

auto
CMatrices::getMatrix(const int64_t atom_a, const std::string& label_a, const int64_t atom_b, const std::string& label_b,
                     const int64_t atom_c, const std::string& label_c) const -> const CMatrix*
{
    const auto key_a = mathfunc::getCantorIndex({atom_a, _to_key(label_a)});

    const auto key_b = mathfunc::getCantorIndex({atom_b, _to_key(label_b)});
    
    const auto key_c = mathfunc::getCantorIndex({atom_c, _to_key(label_c)});

    const auto key_ab = mathfunc::getCantorIndex({key_a, key_b});
    
    return getMatrix(mathfunc::getCantorIndex({key_ab, key_c}));
}

auto
CMatrices::getMatrix(const int64_t atom_a, const std::string& label_a, const int64_t atom_b, const std::string& label_b,
                     const int64_t atom_c, const std::string& label_c, const int64_t atom_d, const std::string& label_d) -> CMatrix*
{
    const auto key_a = mathfunc::getCantorIndex({atom_a, _to_key(label_a)});

    const auto key_b = mathfunc::getCantorIndex({atom_b, _to_key(label_b)});
    
    const auto key_c = mathfunc::getCantorIndex({atom_c, _to_key(label_c)});
    
    const auto key_d = mathfunc::getCantorIndex({atom_d, _to_key(label_d)});

    const auto key_ab = mathfunc::getCantorIndex({key_a, key_b});
    
    const auto key_cd = mathfunc::getCantorIndex({key_c, key_d});
    
    return getMatrix(mathfunc::getCantorIndex({key_ab, key_cd}));
}


auto
CMatrices::getMatrix(const int64_t atom_a, const std::string& label_a, const int64_t atom_b, const std::string& label_b,
                     const int64_t atom_c, const std::string& label_c, const int64_t atom_d, const std::string& label_d) const -> const CMatrix*
{
    const auto key_a = mathfunc::getCantorIndex({atom_a, _to_key(label_a)});

    const auto key_b = mathfunc::getCantorIndex({atom_b, _to_key(label_b)});
    
    const auto key_c = mathfunc::getCantorIndex({atom_c, _to_key(label_c)});
    
    const auto key_d = mathfunc::getCantorIndex({atom_d, _to_key(label_d)});

    const auto key_ab = mathfunc::getCantorIndex({key_a, key_b});
    
    const auto key_cd = mathfunc::getCantorIndex({key_c, key_d});
    
    return getMatrix(mathfunc::getCantorIndex({key_ab, key_cd}));
}

auto
CMatrices::_to_key(const std::string& label) const -> int64_t
{
    const auto index = fstr::upcase(label);

    if (index.size() == 1)
    {
        if (index == "X") return 0;

        if (index == "Y") return 1;

        if (index == "Z") return 2;
    }

    if (index.size() == 2)
    {
        if (index == "XX") return 0;

        if (index == "XY") return 1;

        if (index == "XZ") return 2;

        if (index == "YY") return 3;

        if (index == "YZ") return 4;

        if (index == "ZZ") return 5;
    }

    if (index.size() == 3)
    {
        if (index == "XXX") return 0;

        if (index == "XXY") return 1;

        if (index == "XXZ") return 2;

        if (index == "XYY") return 3;

        if (index == "XYZ") return 4;

        if (index == "XZZ") return 5;

        if (index == "YYY") return 6;

        if (index == "YYZ") return 7;

        if (index == "YZZ") return 8;

        if (index == "ZZZ") return 9;
    }
    
    if (index.size() == 4)
    {
        if (index == "XXXX") return 0;

        if (index == "XXXY") return 1;

        if (index == "XXXZ") return 2;

        if (index == "XXYY") return 3;

        if (index == "XXYZ") return 4;

        if (index == "XXZZ") return 5;

        if (index == "XYYY") return 6;

        if (index == "XYYZ") return 7;

        if (index == "XYZZ") return 8;

        if (index == "XZZZ") return 9;
        
        if (index == "YYYY") return 10;

        if (index == "YYYZ") return 11;

        if (index == "YYZZ") return 12;

        if (index == "YZZZ") return 13;
        
        if (index == "ZZZZ") return 14;
    }
    
    return -1;
}
