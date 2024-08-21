#include "GtoPairBlock.hpp"

#include <algorithm>
#include <cmath>

#include "CustomViews.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"

CGtoPairBlock::CGtoPairBlock(const std::vector<TPoint<double>>& bra_coordinates,
                             const std::vector<TPoint<double>>& ket_coordinates,
                             const std::vector<double>&         bra_exponents,
                             const std::vector<double>&         ket_exponents,
                             const std::vector<double>&         norms,
                             const std::vector<double>&         overlaps,
                             const std::vector<size_t>&         bra_orb_indices,
                             const std::vector<size_t>&         ket_orb_indices,
                             const std::vector<int>&            bra_atm_indices,
                             const std::vector<int>&            ket_atm_indices,
                             const std::pair<int, int>&         angular_momentums,
                             const int                          nppairs)

    : _bra_coordinates(bra_coordinates)

    , _ket_coordinates(ket_coordinates)

    , _bra_exponents(bra_exponents)

    , _ket_exponents(ket_exponents)

    , _norms(norms)

    , _overlaps(overlaps)

    , _bra_orb_indices(bra_orb_indices)

    , _ket_orb_indices(ket_orb_indices)

    , _bra_atm_indices(bra_atm_indices)

    , _ket_atm_indices(ket_atm_indices)

    , _angular_momentums(angular_momentums)

    , _nppairs(nppairs)
{
}

CGtoPairBlock::CGtoPairBlock(const CGtoBlock& gto_block)

    : _bra_coordinates{}

    , _ket_coordinates{}

    , _bra_exponents{}

    , _ket_exponents{}

    , _norms{}

    , _overlaps{}

    , _bra_orb_indices{}

    , _ket_orb_indices{}

    , _bra_atm_indices{}

    , _ket_atm_indices{}

    , _angular_momentums({-1, -1})

    , _nppairs(-1)
{
    // set up pi value

    const auto fpi = mathconst::pi_value();

    // fetch basis functions block data

    const auto coords = gto_block.coordinates();

    const auto gexps = gto_block.exponents();

    const auto gnorms = gto_block.normalization_factors();

    const auto orbidx = gto_block.orbital_indices();

    const auto atmidx = gto_block.atomic_indices();

    // set up dimensions

    auto ncgtos = gto_block.number_of_basis_functions();

    auto npgtos = gto_block.number_of_primitives();

    // allocate memory

    const auto cdim = ncgtos * (ncgtos + 1) / 2;

    _bra_coordinates = std::vector<TPoint<double>>(cdim, TPoint<double>({0.0, 0.0, 0.0}));

    _ket_coordinates = std::vector<TPoint<double>>(cdim, TPoint<double>({0.0, 0.0, 0.0}));

    _bra_orb_indices = std::vector<size_t>(cdim + 1, 0);

    _ket_orb_indices = std::vector<size_t>(cdim + 1, 0);

    _bra_atm_indices = std::vector<int>(cdim, -1);

    _ket_atm_indices = std::vector<int>(cdim, -1);

    const auto pdim = cdim * npgtos * npgtos;

    _bra_exponents = std::vector<double>(pdim, 0.0);

    _ket_exponents = std::vector<double>(pdim, 0.0);

    _norms = std::vector<double>(pdim, 0.0);

    _overlaps = std::vector<double>(pdim, 0.0);

    // set up GTO pairs data

    const auto angmom = gto_block.angular_momentum();

    _bra_orb_indices[0] = orbidx[0];

    _ket_orb_indices[0] = orbidx[0];

    _angular_momentums = {angmom, angmom};

    _nppairs = npgtos * npgtos;

    std::ranges::for_each(views::triangular(ncgtos), [&](const auto& index_ij) {
        const auto [i, j] = index_ij;
        const auto r_a    = coords[i];
        const auto r_b    = coords[j];
        const auto ij_idx = mathfunc::uplo_index(i, j);

        _bra_coordinates[ij_idx]     = r_a;
        _ket_coordinates[ij_idx]     = r_b;
        _bra_orb_indices[ij_idx + 1] = orbidx[i + 1];
        _ket_orb_indices[ij_idx + 1] = orbidx[j + 1];
        _bra_atm_indices[ij_idx]     = atmidx[i];
        _ket_atm_indices[ij_idx]     = atmidx[j];

        const auto r2ab = r_a.distance_square(r_b);
        std::ranges::for_each(views::rectangular(npgtos, npgtos), [&](const auto& index_kl) {
            const auto [i, j] = index_ij;
            const auto [k, l] = index_kl;
            const auto koff   = k * ncgtos + i;
            const auto loff   = l * ncgtos + j;
            const auto ijoff  = (k * npgtos + l) * cdim + ij_idx;

            _bra_exponents[ijoff] = gexps[koff];
            _ket_exponents[ijoff] = gexps[loff];
            _norms[ijoff]         = gnorms[koff] * gnorms[loff];
            const auto fe_ab      = 1.0 / (gexps[koff] + gexps[loff]);
            const auto fz_ab      = gexps[koff] * gexps[loff] * fe_ab;
            const auto fact       = fpi * fe_ab;
            _overlaps[ijoff]      = fact * std::sqrt(fact);
            if (i != j)
            {
                _overlaps[ijoff] *= std::exp(-fz_ab * r2ab);
            }
        });
    });
}

CGtoPairBlock::CGtoPairBlock(const CGtoBlock& bra_gto_block, const CGtoBlock& ket_gto_block)

    : _bra_coordinates{}

    , _ket_coordinates{}

    , _bra_exponents{}

    , _ket_exponents{}

    , _norms{}

    , _overlaps{}

    , _bra_orb_indices{}

    , _ket_orb_indices{}

    , _bra_atm_indices{}

    , _ket_atm_indices{}

    , _angular_momentums({-1, -1})

    , _nppairs(-1)
{
    // set up pi value

    const auto fpi = mathconst::pi_value();

    // fetch GTOs block data on bra side

    const auto bcoords = bra_gto_block.coordinates();

    const auto bexps = bra_gto_block.exponents();

    const auto bnorms = bra_gto_block.normalization_factors();

    const auto borbidx = bra_gto_block.orbital_indices();

    const auto batmidx = bra_gto_block.atomic_indices();

    // fetch GTOs block data on ket side

    const auto kcoords = ket_gto_block.coordinates();

    const auto kexps = ket_gto_block.exponents();

    const auto knorms = ket_gto_block.normalization_factors();

    const auto korbidx = ket_gto_block.orbital_indices();

    const auto katmidx = ket_gto_block.atomic_indices();

    // set up dimensions

    auto bcgtos = bra_gto_block.number_of_basis_functions();

    auto bpgtos = bra_gto_block.number_of_primitives();

    auto kcgtos = ket_gto_block.number_of_basis_functions();

    auto kpgtos = ket_gto_block.number_of_primitives();

    // reserve memory

    const auto cdim = bcgtos * kcgtos;

    _bra_coordinates = std::vector<TPoint<double>>(cdim, TPoint<double>({0.0, 0.0, 0.0}));

    _ket_coordinates = std::vector<TPoint<double>>(cdim, TPoint<double>({0.0, 0.0, 0.0}));

    _bra_orb_indices = std::vector<size_t>(cdim + 1, 0);

    _ket_orb_indices = std::vector<size_t>(cdim + 1, 0);

    _bra_atm_indices = std::vector<int>(cdim, -1);

    _ket_atm_indices = std::vector<int>(cdim, -1);

    const auto pdim = cdim * bpgtos * kpgtos;

    _bra_exponents = std::vector<double>(pdim, 0.0);

    _ket_exponents = std::vector<double>(pdim, 0.0);

    _norms = std::vector<double>(pdim, 0.0);

    _overlaps = std::vector<double>(pdim, 0.0);

    // set up GTO pairs data

    const auto bangmom = bra_gto_block.angular_momentum();

    const auto kangmom = ket_gto_block.angular_momentum();

    _bra_orb_indices[0] = borbidx[0];

    _ket_orb_indices[0] = korbidx[0];

    _angular_momentums = {bangmom, kangmom};

    _nppairs = bpgtos * kpgtos;

    std::ranges::for_each(views::rectangular(bcgtos, kcgtos), [&](const auto& index_ij) {
        const auto [i, j] = index_ij;
        const auto r_a    = bcoords[i];
        const auto r_b    = kcoords[j];
        const auto ij_idx = i * kcgtos + j;

        _bra_coordinates[ij_idx]     = r_a;
        _ket_coordinates[ij_idx]     = r_b;
        _bra_orb_indices[ij_idx + 1] = borbidx[i + 1];
        _ket_orb_indices[ij_idx + 1] = korbidx[j + 1];
        _bra_atm_indices[ij_idx]     = batmidx[i];
        _ket_atm_indices[ij_idx]     = katmidx[j];

        const auto r2ab = r_a.distance_square(r_b);
        std::ranges::for_each(views::rectangular(bpgtos, kpgtos), [&](const auto& index_kl) {
            const auto [i, j] = index_ij;
            const auto [k, l] = index_kl;
            const auto koff   = k * bcgtos + i;
            const auto loff   = l * kcgtos + j;
            const auto ijoff  = (k * kpgtos + l) * cdim + ij_idx;

            _bra_exponents[ijoff] = bexps[koff];
            _ket_exponents[ijoff] = kexps[loff];
            _norms[ijoff]         = bnorms[koff] * knorms[loff];
            const auto fe_ab      = 1.0 / (bexps[koff] + kexps[loff]);
            const auto fz_ab      = bexps[koff] * kexps[loff] * fe_ab;
            const auto fact       = fpi * fe_ab;
            _overlaps[ijoff]      = fact * std::sqrt(fact) * std::exp(-fz_ab * r2ab);
        });
    });
}

auto
CGtoPairBlock::operator==(const CGtoPairBlock& other) const -> bool
{
    if (_angular_momentums != other._angular_momentums)
    {
        return false;
    }
    else if (_nppairs != other._nppairs)
    {
        return false;
    }
    else if (_bra_atm_indices != other._bra_atm_indices)
    {
        return false;
    }
    else if (_ket_atm_indices != other._ket_atm_indices)
    {
        return false;
    }
    else if (_bra_orb_indices != other._bra_orb_indices)
    {
        return false;
    }
    else if (_ket_orb_indices != other._ket_orb_indices)
    {
        return false;
    }
    else if (!std::ranges::equal(_overlaps, other._overlaps, [](auto lhs, auto rhs) -> bool { return mathfunc::equal(lhs, rhs, 1.0e-12, 1.0e-12); }))
    {
        return false;
    }
    else if (!std::ranges::equal(_norms, other._norms, [](auto lhs, auto rhs) -> bool { return mathfunc::equal(lhs, rhs, 1.0e-12, 1.0e-12); }))
    {
        return false;
    }
    else if (!std::ranges::equal(
                 _bra_exponents, other._bra_exponents, [](auto lhs, auto rhs) -> bool { return mathfunc::equal(lhs, rhs, 1.0e-12, 1.0e-12); }))
    {
        return false;
    }
    else if (!std::ranges::equal(
                 _ket_exponents, other._ket_exponents, [](auto lhs, auto rhs) -> bool { return mathfunc::equal(lhs, rhs, 1.0e-12, 1.0e-12); }))
    {
        return false;
    }
    else if (!std::ranges::equal(_bra_coordinates, other._bra_coordinates))
    {
        return false;
    }
    else
    {
        return std::ranges::equal(_ket_coordinates, other._ket_coordinates);
    }
}

auto
CGtoPairBlock::bra_coordinates() const -> std::vector<TPoint<double>>
{
    return _bra_coordinates;
}

auto
CGtoPairBlock::ket_coordinates() const -> std::vector<TPoint<double>>
{
    return _ket_coordinates;
}

auto
CGtoPairBlock::bra_exponents() const -> std::vector<double>
{
    return _bra_exponents;
}

auto
CGtoPairBlock::ket_exponents() const -> std::vector<double>
{
    return _ket_exponents;
}

auto
CGtoPairBlock::normalization_factors() const -> std::vector<double>
{
    return _norms;
}

auto
CGtoPairBlock::overlap_factors() const -> std::vector<double>
{
    return _overlaps;
}

auto
CGtoPairBlock::bra_orbital_indices() const -> std::vector<size_t>
{
    return _bra_orb_indices;
}

auto
CGtoPairBlock::ket_orbital_indices() const -> std::vector<size_t>
{
    return _ket_orb_indices;
}

auto
CGtoPairBlock::bra_atomic_indices() const -> std::vector<int>
{
    return _bra_atm_indices;
}

auto
CGtoPairBlock::ket_atomic_indices() const -> std::vector<int>
{
    return _ket_atm_indices;
}

auto
CGtoPairBlock::angular_momentums() const -> std::pair<int, int>
{
    return _angular_momentums;
}

auto
CGtoPairBlock::number_of_primitive_pairs() const -> int
{
    return _nppairs;
}

auto
CGtoPairBlock::number_of_contracted_pairs() const -> size_t
{
    return _bra_coordinates.size();
}
