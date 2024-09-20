#include "GtoBlock.hpp"

#include <algorithm>
#include <ranges>

#include "ErrorHandler.hpp"
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

CGtoBlock::CGtoBlock(const CMolecularBasis &basis, const CMolecule &molecule, const int angular_momentum, const int npgtos)

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
    else if (!std::ranges::equal(_norms, other._norms, [](auto lhs, auto rhs) -> bool { return mathfunc::equal(lhs, rhs, 1.0e-12, 1.0e-12); }))
    {
        return false;
    }
    else if (!std::ranges::equal(
                 _exponents, other._exponents, [](auto lhs, auto rhs) -> bool { return mathfunc::equal(lhs, rhs, 1.0e-12, 1.0e-12); }))
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
CGtoBlock::getAtomicOrbitalsIndexes() const -> std::vector<int>
{
    std::vector<int> ao_inds;

    // go through spherical harmonics components

    for (int comp = 0; comp < _angular_momentum * 2 + 1; comp++)
    {
        // go through CGTOs in this block
        // note that ind starts from 1
        // because _orb_indices[0] is the total number of CGTOs of _angular_momentum
        // which could be larger than the number of CGTOs in this block

        for (size_t ind = 1; ind < _orb_indices.size(); ind++)
        {
            ao_inds.push_back(comp * _orb_indices[0] + _orb_indices[ind]);
        }
    }

    return ao_inds;
}

auto
CGtoBlock::getAtomicOrbitalsIndexesForCartesian(const int ncgtos_d) const -> std::vector<int>
{
    errors::assertMsgCritical(_angular_momentum <= 3, std::string("GtoBlock: getAtomicOrbitalsIndexesForCartesian only supports up to f-orbitals"));

    if (_angular_momentum == 3)
    {
        errors::assertMsgCritical(ncgtos_d > 0, std::string("GtoBlock: getAtomicOrbitalsIndexesForCartesian needs to know ncgtos of d-orbitals"));
    }

    std::vector<int> ao_inds;

    int orb_ind_shift = 0;

    // take into account the shifting of cart_ind due to 6 d Cartesian components
    if (_angular_momentum == 3) orb_ind_shift = ncgtos_d;

    // go through Cartesian components

    for (int comp = 0; comp < (_angular_momentum + 1) * (_angular_momentum + 2) / 2; comp++)
    {
        // go through CGTOs in this block
        // note that ind starts from 1
        // because _orb_indices[0] is the total number of CGTOs of _angular_momentum
        // which could be larger than the number of CGTOs in this block

        for (size_t ind = 1; ind < _orb_indices.size(); ind++)
        {
            ao_inds.push_back(comp * _orb_indices[0] + _orb_indices[ind] + orb_ind_shift);
        }
    }

    return ao_inds;
}

auto
CGtoBlock::getCartesianToSphericalMappingForP() const -> std::unordered_map<int, std::vector<std::pair<int, double>>>
{
    errors::assertMsgCritical(_angular_momentum == 1, std::string("GtoBlock: getCartesianToSphericalMappingForP only works for p-orbitals"));

    std::unordered_map<int, std::vector<std::pair<int, double>>> cart_sph_p;

    // p-1 (0) <- py (1)
    // p_0 (1) <- pz (2)
    // p+1 (2) <- px (0)

    std::unordered_map<int, std::vector<std::pair<int, double>>> cart_sph_comp_map;

    cart_sph_comp_map[0] = std::vector<std::pair<int, double>>({{2, 1.0}});
    cart_sph_comp_map[1] = std::vector<std::pair<int, double>>({{0, 1.0}});
    cart_sph_comp_map[2] = std::vector<std::pair<int, double>>({{1, 1.0}});

    for (const auto& [cart_comp, sph_comp_coef_vec] : cart_sph_comp_map)
    {
        // go through CGTOs in this block
        // note that ind starts from 1
        // because _orb_indices[0] is the total number of CGTOs of _angular_momentum
        // which could be larger than the number of CGTOs in this block

        for (size_t ind = 1; ind < _orb_indices.size(); ind++)
        {
            auto cart_ind = cart_comp * _orb_indices[0] + _orb_indices[ind];

            cart_sph_p[cart_ind] = std::vector<std::pair<int, double>>();

            for (const auto& sph_comp_coef : sph_comp_coef_vec)
            {
                auto sph_comp = sph_comp_coef.first;
                auto sph_coef = sph_comp_coef.second;

                auto sph_ind = sph_comp * _orb_indices[0] + _orb_indices[ind];

                cart_sph_p[cart_ind].push_back(std::pair<int, double>({sph_ind, sph_coef}));
            }
        }
    }

    return cart_sph_p;
}

auto
CGtoBlock::getCartesianToSphericalMappingForD() const -> std::unordered_map<int, std::vector<std::pair<int, double>>>
{
    errors::assertMsgCritical(_angular_momentum == 2, std::string("GtoBlock: getCartesianToSphericalMappingForD only works for d-orbitals"));

    std::unordered_map<int, std::vector<std::pair<int, double>>> cart_sph_d;

    const double f2_3 = 2.0 * std::sqrt(3.0);

    // d-2 (0) <- dxy(1) * f2_3
    // d-1 (1) <- dyz(4) * f2_3
    // d_0 (2) <- dzz(5) * 2.0 - dxx(0) - dyy(3)
    // d+1 (3) <- dxz(2) * f2_3
    // d+2 (4) <- (dxx(0) - dyy(3)) * 0.5 * f2_3

    std::unordered_map<int, std::vector<std::pair<int, double>>> cart_sph_comp_map;

    cart_sph_comp_map[0] = std::vector<std::pair<int, double>>({{2, -1.0}, {4, 0.5 * f2_3}});
    cart_sph_comp_map[1] = std::vector<std::pair<int, double>>({{0, f2_3}});
    cart_sph_comp_map[2] = std::vector<std::pair<int, double>>({{3, f2_3}});
    cart_sph_comp_map[3] = std::vector<std::pair<int, double>>({{2, -1.0}, {4, -0.5 * f2_3}});
    cart_sph_comp_map[4] = std::vector<std::pair<int, double>>({{1, f2_3}});
    cart_sph_comp_map[5] = std::vector<std::pair<int, double>>({{2, 2.0}});

    for (const auto& [cart_comp, sph_comp_coef_vec] : cart_sph_comp_map)
    {
        // go through CGTOs in this block
        // note that ind starts from 1
        // because _orb_indices[0] is the total number of CGTOs of _angular_momentum
        // which could be larger than the number of CGTOs in this block

        for (size_t ind = 1; ind < _orb_indices.size(); ind++)
        {
            auto cart_ind = cart_comp * _orb_indices[0] + _orb_indices[ind];

            cart_sph_d[cart_ind] = std::vector<std::pair<int, double>>();

            for (const auto& sph_comp_coef : sph_comp_coef_vec)
            {
                auto sph_comp = sph_comp_coef.first;
                auto sph_coef = sph_comp_coef.second;

                auto sph_ind = sph_comp * _orb_indices[0] + _orb_indices[ind];

                cart_sph_d[cart_ind].push_back(std::pair<int, double>({sph_ind, sph_coef}));
            }
        }
    }

    return cart_sph_d;
}

auto
CGtoBlock::getCartesianToSphericalMappingForF(const int ncgtos_d) const -> std::unordered_map<int, std::vector<std::pair<int, double>>>
{
    errors::assertMsgCritical(_angular_momentum == 3, std::string("GtoBlock: getCartesianToSphericalMappingForF only works for f-orbitals"));

    std::unordered_map<int, std::vector<std::pair<int, double>>> cart_sph_f;

    const double f3_5 = std::sqrt(2.5);
    const double f3_15 = 2.0 * std::sqrt(15.0);
    const double f3_3 = std::sqrt(1.5);

    // xxx : 0
    // xxy : 1
    // xxz : 2
    // xyy : 3
    // xyz : 4
    // xzz : 5
    // yyy : 6
    // yyz : 7
    // yzz : 8
    // zzz : 9

    // f-3 (0) <- fxxy(1) * 3.0 * f3_5  + fyyy(6) * (-f3_5)
    // f-2 (1) <- fxyz(4) * f3_15
    // f-1 (2) <- fxxy(1) * (-f3_3)     + fyyy(6) * (-f3_3)        + fyzz(8) * 4.0 * f3_3
    // f_0 (3) <- fxxz(2) * (-3.0)      + fyyz(7) * (-3.0)         + fzzz(9) * 2.0
    // f+1 (4) <- fxxx(0) * (-f3_3)     + fxyy(3) * (-f3_3)        + fxzz(5) * 4.0 * f3_3
    // f+2 (5) <- fxxz(2) * 0.5 * f3_15 + fyyz(7) * (-0.5) * f3_15
    // f+3 (6) <- fxxx(0) * f3_5        + fxyy(3) * (-3.0) * f3_5

    std::unordered_map<int, std::vector<std::pair<int, double>>> cart_sph_comp_map;

    cart_sph_comp_map[0] = std::vector<std::pair<int, double>>({{4, -f3_3}, {6, f3_5}});
    cart_sph_comp_map[1] = std::vector<std::pair<int, double>>({{0, 3.0 * f3_5}, {2, -f3_3}});
    cart_sph_comp_map[2] = std::vector<std::pair<int, double>>({{3, -3.0}, {5, 0.5 * f3_15}});
    cart_sph_comp_map[3] = std::vector<std::pair<int, double>>({{4, -f3_3}, {6, -3.0 * f3_5}});
    cart_sph_comp_map[4] = std::vector<std::pair<int, double>>({{1, f3_15}});
    cart_sph_comp_map[5] = std::vector<std::pair<int, double>>({{4, 4.0 * f3_3}});
    cart_sph_comp_map[6] = std::vector<std::pair<int, double>>({{0, -f3_5}, {2, -f3_3}});
    cart_sph_comp_map[7] = std::vector<std::pair<int, double>>({{3, -3.0}, {5, -0.5 * f3_15}});
    cart_sph_comp_map[8] = std::vector<std::pair<int, double>>({{2, 4.0 * f3_3}});
    cart_sph_comp_map[9] = std::vector<std::pair<int, double>>({{3, 2.0}});

    for (const auto& [cart_comp, sph_comp_coef_vec] : cart_sph_comp_map)
    {
        // go through CGTOs in this block
        // note that ind starts from 1
        // because _orb_indices[0] is the total number of CGTOs of _angular_momentum
        // which could be larger than the number of CGTOs in this block

        for (size_t ind = 1; ind < _orb_indices.size(); ind++)
        {
            auto cart_ind = cart_comp * _orb_indices[0] + _orb_indices[ind];

            // take into account the shifting of cart_ind due to 6 d Cartesian components
            cart_ind += ncgtos_d;

            cart_sph_f[cart_ind] = std::vector<std::pair<int, double>>();

            for (const auto& sph_comp_coef : sph_comp_coef_vec)
            {
                auto sph_comp = sph_comp_coef.first;
                auto sph_coef = sph_comp_coef.second;

                auto sph_ind = sph_comp * _orb_indices[0] + _orb_indices[ind];

                cart_sph_f[cart_ind].push_back(std::pair<int, double>({sph_ind, sph_coef}));
            }
        }
    }

    return cart_sph_f;
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
