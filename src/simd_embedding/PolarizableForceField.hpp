//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#ifndef PolarizableForceField_hpp
#define PolarizableForceField_hpp

#include <string>
#include <vector>

#include "PolarizableSite.hpp"

/// @brief Class CPolarizableForceField holds the force field of one kind of
/// molecule of a polarizable environment.
///
/// @note One site for each atom, in the order the atoms of the molecule are given.
/// Applying the force field to an instance is a walk of the two together, which is
/// why the order is the whole of the correspondence: there is nothing here which
/// says which atom a site belongs to, and a force field handed a molecule whose
/// atoms come in another order answers parameters for the wrong atoms without
/// anything looking amiss.
///
/// @note **One of these serves every molecule of its kind.** The solvent geometry
/// is frozen, so the five hundred waters of a box are five hundred coordinates and
/// one force field, and the parameters are held once.
class CPolarizableForceField
{
   public:
    /// @brief The default constructor.
    CPolarizableForceField() = default;

    /// @brief The constructor with a name.
    /// @param name The name of the kind of molecule, as the structure calls it.
    CPolarizableForceField(const std::string &name);

    /// @brief Gets the name of the kind of molecule.
    auto get_name() const -> const std::string &;

    /// @brief Adds a site to the force field.
    /// @param site The site, which belongs to the atom after the ones already
    /// added.
    auto add_site(const CPolarizableSite &site) -> void;

    /// @brief Gets the number of sites.
    auto number_of_sites() const -> size_t;

    /// @brief Gets a site.
    /// @param index The index of the site, which is the index of its atom.
    auto get_site(const size_t index) const -> const CPolarizableSite &;

    /// @brief Gets the sites.
    auto get_sites() const -> const std::vector<CPolarizableSite> &;

    /// @brief Gets the sum of the permanent charges of the sites.
    /// @return The total charge.
    /// @note A neutral molecule answers zero here to whatever the parameters were
    /// fitted to, and a caller which expects a neutral solvent can see whether it
    /// has one rather than discover it in an energy.
    auto total_charge() const -> double;

    /// @brief Checks whether any site carries a polarizability.
    /// @return True if one does.
    auto is_polarizable() const -> bool;

    /// @brief Checks whether another force field carries the same parameters.
    /// @param other The other force field.
    /// @param tolerance The tolerance the numbers are compared within.
    /// @return True if it does, site for site and in the same order.
    /// @note The name is not compared. This answers whether two force fields say
    /// the same thing, which is what settles whether one name stands for one set
    /// of parameters.
    auto matches(const CPolarizableForceField &other, const double tolerance = 1.0e-12) const -> bool;

   private:
    /// @brief The name of the kind of molecule.
    std::string _name;

    /// @brief The sites, one for each atom and in the order of the atoms.
    std::vector<CPolarizableSite> _sites;
};

#endif /* PolarizableForceField_hpp */
