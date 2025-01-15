#include "BaseCorePotential.hpp"

#include <algorithm>
#include <cmath>
#include <ranges>

#include "MathFunc.hpp"

CBaseCorePotential::CBaseCorePotential()

    : _exponents{}

    , _factors{}

    , _radial_orders{}
{
}

CBaseCorePotential::CBaseCorePotential(const std::vector<double> &exponents,
                                       const std::vector<double> &factors,
                                       const std::vector<int>    &radial_orders)

    : _exponents(exponents)

    , _factors(factors)

    , _radial_orders(radial_orders)
{
}

CBaseCorePotential::CBaseCorePotential(const CBaseCorePotential &other)

    : _exponents(other._exponents)

    , _factors(other._factors)

    , _radial_orders(other._radial_orders)
{
}

CBaseCorePotential::CBaseCorePotential(CBaseCorePotential &&other) noexcept

    : _exponents{}

    , _factors{}

    , _radial_orders{}
{
    std::swap(_exponents, other._exponents);

    std::swap(_factors, other._factors);

    std::swap(_radial_orders, other._radial_orders);
}

auto
CBaseCorePotential::operator=(const CBaseCorePotential &other) -> CBaseCorePotential &
{
    _exponents = other._exponents;

    _factors = other._factors;

    _radial_orders = other._radial_orders;

    return *this;
}

auto
CBaseCorePotential::operator=(CBaseCorePotential &&other) noexcept -> CBaseCorePotential &
{
    std::swap(_exponents, other._exponents);

    std::swap(_factors, other._factors);

    std::swap(_radial_orders, other._radial_orders);

    return *this;
}

auto
CBaseCorePotential::operator==(const CBaseCorePotential &other) const -> bool
{
    if (_radial_orders != other._radial_orders)
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
        return std::ranges::equal(_factors, other._factors, [](auto lhs, auto rhs) -> bool { return mathfunc::equal(lhs, rhs, 1.0e-12, 1.0e-12); });
    }
}

auto
CBaseCorePotential::operator!=(const CBaseCorePotential &other) const -> bool
{
    return !(*this == other);
}

auto
CBaseCorePotential::set_exponents(const std::vector<double> &exponents) -> void
{
    _exponents = exponents;
}

auto
CBaseCorePotential::set_factors(const std::vector<double> &factors) -> void
{
    _factors = factors;
}

auto
CBaseCorePotential::set_radial_orders(const std::vector<int> &radial_orders) -> void
{
    _radial_orders = radial_orders; 
}

auto
CBaseCorePotential::add(const double exponent, const double factor, const int radial_order) -> void
{
    _exponents.push_back(exponent);

    _factors.push_back(factor);
    
    _radial_orders.push_back(radial_order);
}

auto
CBaseCorePotential::get_exponents() const -> std::vector<double>
{
    return _exponents;
}

auto
CBaseCorePotential::get_factors() const -> std::vector<double>
{
    return _factors;
}

auto
CBaseCorePotential::get_radial_orders() const -> std::vector<int>
{
    return _radial_orders;
}

auto
CBaseCorePotential::number_of_primitive_potentials() const -> size_t
{
    return _exponents.size();
}
