//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#include "PolarizableForceField.hpp"

#include "ErrorHandler.hpp"

CPolarizableForceField::CPolarizableForceField(const std::string &name)

    : _name(name)
{
}

auto
CPolarizableForceField::get_name() const -> const std::string &
{
    return _name;
}

auto
CPolarizableForceField::add_site(const CPolarizableSite &site) -> void
{
    _sites.push_back(site);
}

auto
CPolarizableForceField::number_of_sites() const -> size_t
{
    return _sites.size();
}

auto
CPolarizableForceField::get_site(const size_t index) const -> const CPolarizableSite &
{
    errors::assertMsgCritical(index < _sites.size(),
                              std::string("PolarizableForceField: There is no site of this index"));

    return _sites[index];
}

auto
CPolarizableForceField::get_sites() const -> const std::vector<CPolarizableSite> &
{
    return _sites;
}

auto
CPolarizableForceField::total_charge() const -> double
{
    double total = 0.0;

    for (const auto &site : _sites) total += site.get_charge();

    return total;
}

auto
CPolarizableForceField::is_polarizable() const -> bool
{
    for (const auto &site : _sites)
    {
        if (site.is_polarizable()) return true;
    }

    return false;
}
