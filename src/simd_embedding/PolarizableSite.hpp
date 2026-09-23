//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#ifndef PolarizableSite_hpp
#define PolarizableSite_hpp

#include <array>
#include <cstddef>

/// @brief The form a site's charge and polarizability take.
/// @note A point site carries its moments at a point and is singular where the
/// distance goes to zero; a Gaussian one spreads them over a normalized spherical
/// Gaussian of a given width and is finite everywhere. The widths belong to the
/// Gaussian form alone and mean nothing for a point.
enum class pesite
{
    point,
    gaussian
};

/// @brief Class CPolarizableSite holds the force field parameters of one site of
/// a polarizable environment.
///
/// @note **There are no coordinates here.** These are the parameters of a kind of
/// site and not of one instance of it: the solvent geometry is frozen, so the
/// water molecules of a box share one set of them and differ only in where they
/// sit. A site is paired with a coordinate when the force field is applied.
///
/// @note Which follows from that, and is the caller's to get right: a charge is a
/// number and transfers as it stands, but a dipole is a vector and a quadrupole is
/// a tensor, and both are written in the frame the parameters were made in. An
/// instance turned some other way needs them turned with it. Nothing here rotates
/// anything, and nothing here can tell whether a caller forgot to: the charges
/// stay right while the higher moments point the wrong way, which reads as a small
/// error in the energy rather than as a mistake.
///
/// @note Atomic units throughout, as the rest of this layer.
class CPolarizableSite
{
   public:
    /// @brief The default constructor, which leaves a site with no moments and no
    /// polarizability.
    CPolarizableSite() = default;

    /// @brief The constructor of a point site.
    /// @param order The highest moment the site carries: zero for a charge, one
    /// for a charge and a dipole, two for those and a quadrupole.
    /// @param charge The permanent charge.
    /// @param dipole The permanent dipole, x y z.
    /// @param quadrupole The permanent quadrupole, xx xy xz yy yz zz.
    CPolarizableSite(const int                    order,
                     const double                 charge,
                     const std::array<double, 3> &dipole,
                     const std::array<double, 6> &quadrupole);

    /// @brief Gets the highest moment the site carries.
    /// @return Zero, one or two.
    /// @note A moment above this is absent and not zero. The two are different
    /// things to whoever contracts them: an absent one is skipped where a zero one
    /// is multiplied.
    auto get_order() const -> int;

    /// @brief Gets the permanent charge.
    auto get_charge() const -> double;

    /// @brief Gets the permanent dipole, x y z.
    auto get_dipole() const -> const std::array<double, 3> &;

    /// @brief Gets the permanent quadrupole, xx xy xz yy yz zz.
    auto get_quadrupole() const -> const std::array<double, 6> &;

    /// @brief Gets the trace of the quadrupole.
    /// @return The sum of its diagonal.
    /// @note The moments are kept as they were given and are not made traceless.
    /// A potential which supplies them traceless answers zero here and one which
    /// does not answers what it carries, so a caller which needs one convention
    /// can see which it has rather than assume.
    auto quadrupole_trace() const -> double;

    /// @brief Sets the polarizability to an isotropic one.
    /// @param alpha The polarizability.
    auto set_isotropic_polarizability(const double alpha) -> void;

    /// @brief Sets the polarizability to an anisotropic one.
    /// @param alpha The polarizability, xx xy xz yy yz zz.
    auto set_polarizability(const std::array<double, 6> &alpha) -> void;

    /// @brief Checks whether the polarizability is isotropic.
    /// @return True if it is.
    /// @note Kept although the six components are stored either way: an isotropic
    /// site's induced dipole is a number times a field where an anisotropic one is
    /// a three by three solve, and the distinction is worth more than the five
    /// doubles it saves.
    auto is_isotropic() const -> bool;

    /// @brief Checks whether the site is polarizable at all.
    /// @return True if it carries a polarizability.
    auto is_polarizable() const -> bool;

    /// @brief Gets the polarizability, xx xy xz yy yz zz.
    /// @note An isotropic polarizability answers here as well, on the diagonal.
    auto get_polarizability() const -> const std::array<double, 6> &;

    /// @brief Gets the isotropic polarizability.
    /// @return The polarizability.
    /// @note Refused for an anisotropic site rather than answering a third of the
    /// trace, which is a number nobody asked for.
    auto get_isotropic_polarizability() const -> double;

    /// @brief Sets the form of the site to a Gaussian of the given widths.
    /// @param multipole_width The width of the charge distribution.
    /// @param polarizability_width The width of the polarizability.
    auto set_gaussian(const double multipole_width, const double polarizability_width) -> void;

    /// @brief Gets the form of the site.
    auto get_form() const -> pesite;

    /// @brief Gets the width of the charge distribution.
    /// @return The width.
    /// @note Refused for a point site, which has no width: answering a zero would
    /// be taken for a width of zero and divided by.
    auto get_multipole_width() const -> double;

    /// @brief Gets the width of the polarizability.
    /// @return The width.
    /// @note Refused for a point site, for the reason above.
    auto get_polarizability_width() const -> double;

   private:
    /// @brief The highest moment the site carries.
    int _order = 0;

    /// @brief The permanent charge.
    double _charge = 0.0;

    /// @brief The permanent dipole, x y z.
    std::array<double, 3> _dipole = {0.0, 0.0, 0.0};

    /// @brief The permanent quadrupole, xx xy xz yy yz zz.
    std::array<double, 6> _quadrupole = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    /// @brief Whether the polarizability is isotropic.
    bool _isotropic = true;

    /// @brief Whether the site carries a polarizability.
    bool _polarizable = false;

    /// @brief The polarizability, xx xy xz yy yz zz.
    std::array<double, 6> _polarizability = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    /// @brief The form of the site.
    pesite _form = pesite::point;

    /// @brief The width of the charge distribution, for a Gaussian site.
    double _multipole_width = 0.0;

    /// @brief The width of the polarizability, for a Gaussian site.
    double _polarizability_width = 0.0;
};

#endif /* PolarizableSite_hpp */
