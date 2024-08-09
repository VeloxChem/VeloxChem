#include "GtoBlock.hpp"

#include <algorithm>
#include <ranges>

#include "MathFunc.hpp"

CGtoBlock::CGtoBlock()

    : _coordinates{}

    , _exponents{}

    , _norms{}

    , _orb_indices{}

    , _atm_indices{}

    , _angular_momentum{-1}

    , _npgtos{0}
{
}

CGtoBlock::CGtoBlock(const std::vector<TPoint<double>> &coordinates,
                     const std::vector<double>         &exponents,
                     const std::vector<double>         &norms,
                     const std::vector<size_t>         &orb_indices,
                     const std::vector<int>            &atm_indices,
                     const int                          angular_momentum,
                     const int                          npgtos)

    : _coordinates(coordinates)

    , _exponents(exponents)

    , _norms(norms)

    , _orb_indices(orb_indices)

    , _atm_indices(atm_indices)

    , _angular_momentum(angular_momentum)

    , _npgtos(npgtos)
{
}

CGtoBlock::CGtoBlock(const CMolecularBasis &basis,
                     const CMolecule       &molecule,
                     const int              angular_momentum,
                     const int              npgtos)

    : _coordinates{}

    , _exponents{}

    , _norms{}

    , _orb_indices{}

    , _atm_indices{}

    , _angular_momentum{-1}

    , _npgtos{0}
{
    if (const auto gtos = basis.basis_functions(angular_momentum, npgtos); !gtos.empty())
    {
        _angular_momentum = angular_momentum;

        _npgtos = npgtos;

        _orb_indices = basis.index_map(angular_momentum, npgtos);

        _atm_indices = basis.atomic_indices(angular_momentum, npgtos);

        _coordinates.reserve(_atm_indices.size());

        std::ranges::for_each(_atm_indices, [&](const int i) { _coordinates.push_back(molecule.atom_coordinates(i)); });

        const auto ncgtos = static_cast<int>(gtos.size());

        _exponents = std::vector<double>(ncgtos * npgtos, 0.0);

        _norms = std::vector<double>(ncgtos * npgtos, 0.0);

        std::ranges::for_each(std::views::iota(0, ncgtos), [&](const int i) {
            const auto fexps  = gtos[i].get_exponents();
            const auto fnorms = gtos[i].get_normalization_factors();
            std::ranges::for_each(std::views::iota(0, npgtos), [&](const int j) {
                _exponents[j * ncgtos + i] = fexps[j];
                _norms[j * ncgtos + i]     = fnorms[j];
            });
        });
    }
}

CGtoBlock::CGtoBlock(const CMolecularBasis  &basis,
                     const CMolecule        &molecule,
                     const std::vector<int> &atoms,
                     const int               angular_momentum,
                     const int               npgtos)
    : _coordinates{}

    , _exponents{}

    , _norms{}

    , _orb_indices{}

    , _atm_indices{}

    , _angular_momentum{-1}

    , _npgtos{0}
{
    if (const auto gtos = basis.basis_functions(atoms, angular_momentum, npgtos); !gtos.empty())
    {
        _angular_momentum = angular_momentum;

        _npgtos = npgtos;

        _orb_indices = basis.index_map(atoms, angular_momentum, npgtos);

        _atm_indices = basis.atomic_indices(atoms, angular_momentum, npgtos);

        _coordinates.reserve(_atm_indices.size());

        std::ranges::for_each(_atm_indices, [&](const int i) { _coordinates.push_back(molecule.atom_coordinates(i)); });

        const auto ncgtos = static_cast<int>(gtos.size());

        _exponents = std::vector<double>(ncgtos * npgtos, 0.0);

        _norms = std::vector<double>(ncgtos * npgtos, 0.0);

        std::ranges::for_each(std::views::iota(0, ncgtos), [&](const int i) {
            const auto fexps  = gtos[i].get_exponents();
            const auto fnorms = gtos[i].get_normalization_factors();
            std::ranges::for_each(std::views::iota(0, npgtos), [&](const int j) {
                _exponents[j * ncgtos + i] = fexps[j];
                _norms[j * ncgtos + i]     = fnorms[j];
            });
        });
    }
}

CGtoBlock::CGtoBlock(const CGtoBlock &other)

    : _coordinates(other._coordinates)

    , _exponents(other._exponents)

    , _norms(other._norms)

    , _orb_indices(other._orb_indices)

    , _atm_indices(other._atm_indices)

    , _angular_momentum(other._angular_momentum)

    , _npgtos(other._npgtos)
{
}

CGtoBlock::CGtoBlock(CGtoBlock &&other) noexcept

    : _coordinates{}

    , _exponents{}

    , _norms{}

    , _orb_indices{}

    , _atm_indices{}

    , _angular_momentum{-1}

    , _npgtos{0}
{
    std::swap(_coordinates, other._coordinates);

    std::swap(_exponents, other._exponents);

    std::swap(_norms, other._norms);

    std::swap(_orb_indices, other._orb_indices);

    std::swap(_atm_indices, other._atm_indices);

    std::swap(_angular_momentum, other._angular_momentum);

    std::swap(_npgtos, other._npgtos);
}

auto
CGtoBlock::operator=(const CGtoBlock &other) -> CGtoBlock &
{
    _coordinates = other._coordinates;

    _exponents = other._exponents;

    _norms = other._norms;

    _orb_indices = other._orb_indices;

    _atm_indices = other._atm_indices;

    _angular_momentum = other._angular_momentum;

    _npgtos = other._npgtos;

    return *this;
}

auto
CGtoBlock::operator=(CGtoBlock &&other) noexcept -> CGtoBlock &
{
    std::swap(_coordinates, other._coordinates);

    std::swap(_exponents, other._exponents);

    std::swap(_norms, other._norms);

    std::swap(_orb_indices, other._orb_indices);

    std::swap(_atm_indices, other._atm_indices);

    std::swap(_angular_momentum, other._angular_momentum);

    std::swap(_npgtos, other._npgtos);

    return *this;
}

auto
CGtoBlock::operator==(const CGtoBlock &other) const -> bool
{
    if (_angular_momentum != other._angular_momentum)
    {
        return false;
    }
    else if (_npgtos != other._npgtos)
    {
        return false;
    }
    else if (_atm_indices != other._atm_indices)
    {
        return false;
    }
    else if (_orb_indices != other._orb_indices)
    {
        return false;
    }
    else if (!std::ranges::equal(_norms, other._norms, [](auto lhs, auto rhs) -> bool {
                 return mathfunc::equal(lhs, rhs, 1.0e-12, 1.0e-12);
             }))
    {
        return false;
    }
    else if (!std::ranges::equal(_exponents, other._exponents, [](auto lhs, auto rhs) -> bool {
                 return mathfunc::equal(lhs, rhs, 1.0e-12, 1.0e-12);
             }))
    {
        return false;
    }
    else
    {
        return _coordinates == other._coordinates;
    }
}

auto
CGtoBlock::operator!=(const CGtoBlock &other) const -> bool
{
    return !(*this == other);
}

auto
CGtoBlock::coordinates() const -> std::vector<TPoint<double>>
{
    return _coordinates;
}

auto
CGtoBlock::exponents() const -> std::vector<double>
{
    return _exponents;
}

auto
CGtoBlock::normalization_factors() const -> std::vector<double>
{
    return _norms;
}

auto
CGtoBlock::orbital_indices() const -> std::vector<size_t>
{
    return _orb_indices;
}

auto
CGtoBlock::atomic_indices() const -> std::vector<int>
{
    return _atm_indices;
}

auto
CGtoBlock::angular_momentum() const -> int
{
    return _angular_momentum;
}

auto
CGtoBlock::number_of_primitives() const -> int
{
    return _npgtos;
}

auto
CGtoBlock::number_of_basis_functions() const -> int
{
    return static_cast<int>(_coordinates.size());
}
