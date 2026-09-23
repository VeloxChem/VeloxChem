//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#include "PolarizableSite.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "EmbeddingError.hpp"

CPolarizableSite::CPolarizableSite(const int                    order,
                                   const double                 charge,
                                   const std::array<double, 3> &dipole,
                                   const std::array<double, 6> &quadrupole)

    : _order(order)
    , _charge(charge)
    , _dipole(dipole)
    , _quadrupole(quadrupole)
{
    embedding::require((order >= 0) && (order <= 2),
                              std::string("PolarizableSite: The highest moment of a site is the quadrupole"));
}

auto
CPolarizableSite::get_order() const -> int
{
    return _order;
}

auto
CPolarizableSite::get_charge() const -> double
{
    return _charge;
}

auto
CPolarizableSite::get_dipole() const -> const std::array<double, 3> &
{
    return _dipole;
}

auto
CPolarizableSite::get_quadrupole() const -> const std::array<double, 6> &
{
    return _quadrupole;
}

auto
CPolarizableSite::quadrupole_trace() const -> double
{
    // the diagonal of the packed upper triangle: xx, yy and zz

    return _quadrupole[0] + _quadrupole[3] + _quadrupole[5];
}

auto
CPolarizableSite::set_isotropic_polarizability(const double alpha) -> void
{
    _isotropic = true;

    _polarizable = true;

    // NOTE: written on the diagonal as well, so that a caller which takes the six
    // components does not have to ask which kind of site it has.

    _polarizability = {alpha, 0.0, 0.0, alpha, 0.0, alpha};
}

auto
CPolarizableSite::set_polarizability(const std::array<double, 6> &alpha) -> void
{
    _isotropic = false;

    _polarizable = true;

    _polarizability = alpha;
}

auto
CPolarizableSite::is_isotropic() const -> bool
{
    return _isotropic;
}

auto
CPolarizableSite::is_polarizable() const -> bool
{
    return _polarizable;
}

auto
CPolarizableSite::get_polarizability() const -> const std::array<double, 6> &
{
    return _polarizability;
}

auto
CPolarizableSite::get_isotropic_polarizability() const -> double
{
    embedding::require(_isotropic,
                              std::string("PolarizableSite: The polarizability of this site is anisotropic and is "
                                          "not one number"));

    return _polarizability[0];
}

auto
CPolarizableSite::set_gaussian(const double multipole_width, const double polarizability_width) -> void
{
    embedding::require((multipole_width > 0.0) && (polarizability_width > 0.0),
                              std::string("PolarizableSite: The width of a Gaussian site is positive"));

    _form = pesite::gaussian;

    _multipole_width = multipole_width;

    _polarizability_width = polarizability_width;
}

auto
CPolarizableSite::get_form() const -> pesite
{
    return _form;
}

auto
CPolarizableSite::get_multipole_width() const -> double
{
    embedding::require(_form == pesite::gaussian,
                              std::string("PolarizableSite: A point site has no width"));

    return _multipole_width;
}

auto
CPolarizableSite::get_polarizability_width() const -> double
{
    embedding::require(_form == pesite::gaussian,
                              std::string("PolarizableSite: A point site has no width"));

    return _polarizability_width;
}

auto
CPolarizableSite::matches(const CPolarizableSite &other, const double tolerance) const -> bool
{
    auto close = [&](const double lhs, const double rhs) -> bool {
        return std::fabs(lhs - rhs) <= tolerance * std::max({1.0, std::fabs(lhs), std::fabs(rhs)});
    };

    if (_order != other._order) return false;

    if (_form != other._form) return false;

    if (_polarizable != other._polarizable) return false;

    if (_isotropic != other._isotropic) return false;

    if (!close(_charge, other._charge)) return false;

    for (size_t i = 0; i < _dipole.size(); i++)
    {
        if (!close(_dipole[i], other._dipole[i])) return false;
    }

    for (size_t i = 0; i < _quadrupole.size(); i++)
    {
        if (!close(_quadrupole[i], other._quadrupole[i])) return false;
    }

    for (size_t i = 0; i < _polarizability.size(); i++)
    {
        if (!close(_polarizability[i], other._polarizability[i])) return false;
    }

    if (_form == pesite::gaussian)
    {
        if (!close(_multipole_width, other._multipole_width)) return false;

        if (!close(_polarizability_width, other._polarizability_width)) return false;
    }

    return true;
}
