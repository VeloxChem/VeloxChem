//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "ExportEmbedding.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "EmbeddingRegion.hpp"
#include "ExportGeneral.hpp"
#include "PermanentMultipoles.hpp"
#include "PolarizableEmbedding.hpp"
#include "PolarizableForceField.hpp"
#include "PolarizableSite.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_embedding {  // vlx_embedding namespace

/// @brief An array of doubles which a list of lists is accepted for.
using CArrayOfDoubles = py::array_t<double, py::array::c_style | py::array::forcecast>;

/// @brief Refuses what the caller has asked for.
/// @param message The reason.
/// @note This layer throws where the classes underneath it assert, because what
/// is refused here is what somebody typed. The classes are asserted against a
/// mistake in the code which calls them, which is a reason to stop; a misspelled
/// parameter is a reason to say so and let the caller try again.
static auto
_refuse(const std::string &message) -> void
{
    throw std::runtime_error(message);
}

/// @brief Reads a vector of three from a Python object.
/// @param obj The object.
/// @param name The name of the keyword, for the message if it is not one.
/// @return The vector, x y z.
static auto
_vector_from(const py::object &obj, const std::string &name) -> std::array<double, 3>
{
    const auto arr = obj.cast<CArrayOfDoubles>();

    if ((arr.ndim() != 1) || (arr.shape(0) != 3))
    {
        _refuse(std::string("PolarizableSite: The ") + name + std::string(" must have three components"));
    }

    const auto val = arr.unchecked<1>();

    return std::array<double, 3>({val(0), val(1), val(2)});
}

/// @brief Reads the six components of a symmetric tensor from a three by three
/// Python object.
/// @param obj The object.
/// @param name The name of the keyword, for the message if it is not one.
/// @return The components, xx xy xz yy yz zz.
/// @note The tensor must be symmetric, and one which is not is refused rather
/// than halved or read from one triangle. An array whose layout is not the one
/// the caller thinks it is arrives here as an unsymmetric one, and that is the
/// only sign of it: every component is a plausible number.
static auto
_symmetric_from(const py::object &obj, const std::string &name) -> std::array<double, 6>
{
    const auto arr = obj.cast<CArrayOfDoubles>();

    if ((arr.ndim() != 2) || (arr.shape(0) != 3) || (arr.shape(1) != 3))
    {
        _refuse(std::string("PolarizableSite: The ") + name + std::string(" must be a three by three array"));
    }

    const auto val = arr.unchecked<2>();

    for (int i = 0; i < 3; i++)
    {
        for (int j = i + 1; j < 3; j++)
        {
            const double scale = std::max({1.0, std::fabs(val(i, j)), std::fabs(val(j, i))});

            if (std::fabs(val(i, j) - val(j, i)) > 1.0e-10 * scale)
            {
                _refuse(std::string("PolarizableSite: The ") + name + std::string(" must be symmetric"));
            }
        }
    }

    return std::array<double, 6>({val(0, 0), val(0, 1), val(0, 2), val(1, 1), val(1, 2), val(2, 2)});
}

/// @brief Makes a three by three array of a symmetric tensor.
/// @param tensor The components, xx xy xz yy yz zz.
/// @return The array.
static auto
_array_of_symmetric(const std::array<double, 6> &tensor) -> py::array_t<double>
{
    const std::array<double, 9> full({tensor[0],
                                      tensor[1],
                                      tensor[2],
                                      tensor[1],
                                      tensor[3],
                                      tensor[4],
                                      tensor[2],
                                      tensor[4],
                                      tensor[5]});

    return vlx_general::pointer_to_numpy(full.data(), {3, 3});
}

/// @brief Sets the polarizability of a site from a Python object.
/// @param site The site.
/// @param alpha A number for an isotropic polarizability, a three by three array
/// for an anisotropic one.
static auto
_set_polarizability(CPolarizableSite &site, const py::object &alpha) -> void
{
    if (py::isinstance<py::float_>(alpha) || py::isinstance<py::int_>(alpha))
    {
        const double value = alpha.cast<double>();

        if (value < 0.0) _refuse(std::string("PolarizableSite: The polarizability is negative"));

        site.set_isotropic_polarizability(value);
    }
    else
    {
        site.set_polarizability(_symmetric_from(alpha, "polarizability"));
    }
}

/// @brief Makes a site of the keywords which describe it.
/// @param charge The permanent charge.
/// @param dipole The permanent dipole, or none.
/// @param quadrupole The permanent quadrupole, or none.
/// @param polarizability A number for an isotropic polarizability, a three by
/// three array for an anisotropic one, or none.
/// @param multipole_width The width of the Gaussian charges, or none.
/// @param polarizability_width The width of the Gaussian polarizabilities, or none.
/// @return The site.
/// @note The order of the site follows from the moments which are given: a
/// quadrupole makes it two, a dipole alone makes it one, neither makes it zero.
static auto
_make_site(const double      charge,
           const py::object &dipole,
           const py::object &quadrupole,
           const py::object &polarizability,
           const py::object &multipole_width,
           const py::object &polarizability_width) -> std::shared_ptr<CPolarizableSite>
{
    const bool has_dipole = !dipole.is_none();

    const bool has_quadrupole = !quadrupole.is_none();

    const int order = has_quadrupole ? 2 : (has_dipole ? 1 : 0);

    const auto mu = has_dipole ? _vector_from(dipole, "dipole") : std::array<double, 3>({0.0, 0.0, 0.0});

    const auto theta = has_quadrupole ? _symmetric_from(quadrupole, "quadrupole")
                                      : std::array<double, 6>({0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

    auto site = std::make_shared<CPolarizableSite>(order, charge, mu, theta);

    if (!polarizability.is_none()) _set_polarizability(*site, polarizability);

    const bool has_multipole_width = !multipole_width.is_none();

    const bool has_polarizability_width = !polarizability_width.is_none();

    if (has_multipole_width != has_polarizability_width)
    {
        _refuse(std::string("PolarizableSite: A Gaussian site needs both of multipole_width and polarizability_width"));
    }

    if (has_multipole_width)
    {
        const double width = multipole_width.cast<double>();

        const double alpha_width = polarizability_width.cast<double>();

        if ((width <= 0.0) || (alpha_width <= 0.0))
        {
            _refuse(std::string("PolarizableSite: The width of a Gaussian site is not positive"));
        }

        site->set_gaussian(width, alpha_width);
    }

    return site;
}

/// @brief Makes a site of a description of one.
/// @param desc The description, which is a site or the keywords of one.
/// @return The site.
/// @note A key which is not one of the keywords is refused rather than passed
/// over, so that a misspelled parameter is not silently absent from the site.
static auto
_site_of(const py::handle &desc) -> std::shared_ptr<CPolarizableSite>
{
    if (py::isinstance<CPolarizableSite>(desc))
    {
        return std::make_shared<CPolarizableSite>(desc.cast<const CPolarizableSite &>());
    }

    if (!py::isinstance<py::dict>(desc))
    {
        _refuse(std::string("PolarizableForceField: A site is given as a dictionary of parameters or as a PolarizableSite"));
    }

    const auto keywords = desc.cast<py::dict>();

    const std::vector<std::string> allowed(
        {"charge", "dipole", "quadrupole", "polarizability", "multipole_width", "polarizability_width"});

    for (const auto &item : keywords)
    {
        const auto key = item.first.cast<std::string>();

        if (std::find(allowed.begin(), allowed.end(), key) == allowed.end())
        {
            _refuse(std::string("PolarizableForceField: A site has no parameter named ") + key);
        }
    }

    auto value_of = [&](const char *key) -> py::object {
        return keywords.contains(key) ? keywords[key].cast<py::object>() : py::none();
    };

    const auto charge = keywords.contains("charge") ? keywords["charge"].cast<double>() : 0.0;

    return _make_site(charge,
                      value_of("dipole"),
                      value_of("quadrupole"),
                      value_of("polarizability"),
                      value_of("multipole_width"),
                      value_of("polarizability_width"));
}

/// @brief Makes a force field of a name and the descriptions of its sites.
/// @param name The name of the kind of molecule.
/// @param sites The descriptions, one for each atom and in the order of the atoms.
/// @return The force field.
static auto
_make_force_field(const std::string &name, const py::object &sites) -> std::shared_ptr<CPolarizableForceField>
{
    auto force_field = std::make_shared<CPolarizableForceField>(name);

    if (!sites.is_none())
    {
        for (const auto &desc : sites)
        {
            force_field->add_site(*_site_of(desc));
        }
    }

    return force_field;
}

/// @brief Gets the index of a site, counted from the end if it is negative.
/// @param self The force field.
/// @param index The index.
/// @return The index counted from the start.
static auto
_index_of_site(const CPolarizableForceField &self, const py::ssize_t index) -> size_t
{
    const auto nsites = static_cast<py::ssize_t>(self.number_of_sites());

    const auto at = (index < 0) ? index + nsites : index;

    if ((at < 0) || (at >= nsites)) throw py::index_error("PolarizableForceField: There is no site of this index");

    return static_cast<size_t>(at);
}

/// @brief Makes the arrays of a set of permanent multipoles.
/// @param gathered The multipoles.
/// @param owner The Python object which holds them.
/// @return The coordinates, of shape (n, 3), and the values, of shape (n,) for
/// the charges and (n, 3) or (n, 6) above them.
/// @note Views onto what the embedding holds and not copies of it, so asking
/// for them costs nothing however often it is done. Each array takes the object
/// it came from as its base, so the environment cannot be collected while a
/// view of it is still held.
/// @note They are still only as good as the arrays underneath: adding a
/// molecule to the environment builds those again and leaves a view which was
/// taken before it pointing at what was there then. Take them again after
/// adding.
/// @brief Makes an array read only.
/// @param array The array.
/// @return The array, which can no longer be written through.
/// @note These are views onto what the environment holds, so writing through
/// one would change the environment itself and every later reader of it, with
/// nothing to say that it had happened. A caller which wants to change the
/// numbers takes a copy.
static auto
_read_only(py::array_t<double> array) -> py::array_t<double>
{
    py::detail::array_proxy(array.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;

    return array;
}

static auto
_arrays_of(const TPermanentMultipoles &gathered, const py::object &owner) -> py::tuple
{
    const auto nsites = static_cast<py::ssize_t>(gathered.number_of_sites());

    const auto ncomponents = static_cast<py::ssize_t>(gathered.components());

    const auto stride = static_cast<py::ssize_t>(sizeof(double));

    auto coordinates = py::array_t<double>(
        {nsites, static_cast<py::ssize_t>(3)}, {3 * stride, stride}, gathered.coordinates.data(), owner);

    auto values = (ncomponents == 1)
                      ? py::array_t<double>({nsites}, {stride}, gathered.values.data(), owner)
                      : py::array_t<double>({nsites, ncomponents}, {ncomponents * stride, stride}, gathered.values.data(), owner);

    return py::make_tuple(_read_only(coordinates), _read_only(values));
}

auto
export_embedding(py::module &m) -> void
{
    // pesite enumeration

    py::enum_<pesite>(m, "pesite").value("point", pesite::point).value("gaussian", pesite::gaussian);

    // CPolarizableSite class

    PyClass<CPolarizableSite>(m, "PolarizableSite")
        .def(py::init<>())
        .def(py::init(&_make_site),
             "charge"_a               = 0.0,
             "dipole"_a               = py::none(),
             "quadrupole"_a           = py::none(),
             "polarizability"_a       = py::none(),
             "multipole_width"_a      = py::none(),
             "polarizability_width"_a = py::none())
        .def("get_order", &CPolarizableSite::get_order, "Gets the highest moment the site carries.")
        .def("get_charge", &CPolarizableSite::get_charge, "Gets the permanent charge.")
        .def(
            "get_dipole",
            [](const CPolarizableSite &self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.get_dipole().data(), {3});
            },
            "Gets the permanent dipole as a vector of three.")
        .def(
            "get_quadrupole",
            [](const CPolarizableSite &self) -> py::array_t<double> { return _array_of_symmetric(self.get_quadrupole()); },
            "Gets the permanent quadrupole as a three by three array.")
        .def("quadrupole_trace", &CPolarizableSite::quadrupole_trace, "Gets the trace of the quadrupole.")
        .def(
            "set_polarizability",
            [](CPolarizableSite &self, const py::object &alpha) -> void { _set_polarizability(self, alpha); },
            "Sets a polarizability, isotropic if it is a number and anisotropic if it is a three by three array.",
            "alpha"_a)
        .def("is_isotropic", &CPolarizableSite::is_isotropic, "Checks whether the polarizability is isotropic.")
        .def("is_polarizable", &CPolarizableSite::is_polarizable, "Checks whether the site carries a polarizability.")
        .def(
            "get_polarizability",
            [](const CPolarizableSite &self) -> py::array_t<double> { return _array_of_symmetric(self.get_polarizability()); },
            "Gets the polarizability as a three by three array.")
        .def("get_isotropic_polarizability",
             &CPolarizableSite::get_isotropic_polarizability,
             "Gets the isotropic polarizability, which an anisotropic site has none of.")
        .def("set_gaussian",
             &CPolarizableSite::set_gaussian,
             "Makes the site a Gaussian one of these widths.",
             "multipole_width"_a,
             "polarizability_width"_a)
        .def("get_form", &CPolarizableSite::get_form, "Gets the form of the site, point or Gaussian.")
        .def("get_multipole_width",
             &CPolarizableSite::get_multipole_width,
             "Gets the width of the Gaussian charges, which a point site has none of.")
        .def("get_polarizability_width",
             &CPolarizableSite::get_polarizability_width,
             "Gets the width of the Gaussian polarizabilities, which a point site has none of.")
        .def("matches",
             &CPolarizableSite::matches,
             "Checks whether another site carries the same parameters.",
             "other"_a,
             "tolerance"_a = 1.0e-12);

    // CPolarizableForceField class

    PyClass<CPolarizableForceField>(m, "PolarizableForceField")
        .def(py::init<>())
        .def(py::init(&_make_force_field), "name"_a = std::string(), "sites"_a = py::none())
        .def("get_name", &CPolarizableForceField::get_name, "Gets the name of the kind of molecule.")
        .def(
            "add_site",
            [](CPolarizableForceField &self, const py::handle &desc) -> void { self.add_site(*_site_of(desc)); },
            "Adds the site of the atom after the ones already added.",
            "site"_a)
        .def("number_of_sites", &CPolarizableForceField::number_of_sites, "Gets the number of sites.")
        .def(
            "get_site",
            [](const CPolarizableForceField &self, const py::ssize_t index) -> const CPolarizableSite & {
                return self.get_site(_index_of_site(self, index));
            },
            "Gets the site of the atom of this index.",
            "index"_a,
            py::return_value_policy::reference_internal)
        .def("get_sites", &CPolarizableForceField::get_sites, "Gets the sites, in the order of the atoms.")
        .def("total_charge", &CPolarizableForceField::total_charge, "Gets the sum of the permanent charges.")
        .def("is_polarizable", &CPolarizableForceField::is_polarizable, "Checks whether any site carries a polarizability.")
        .def("__len__", &CPolarizableForceField::number_of_sites)
        .def(
            "__getitem__",
            [](const CPolarizableForceField &self, const py::ssize_t index) -> const CPolarizableSite & {
                return self.get_site(_index_of_site(self, index));
            },
            py::return_value_policy::reference_internal)
        .def("matches",
             &CPolarizableForceField::matches,
             "Checks whether another force field carries the same parameters.",
             "other"_a,
             "tolerance"_a = 1.0e-12);

    // CEmbeddingRegion class

    PyClass<CEmbeddingRegion>(m, "EmbeddingRegion")
        .def(py::init<>())
        .def(py::init<const bool>(), "allows_polarizabilities"_a)
        .def("allows_polarizabilities",
             &CEmbeddingRegion::allows_polarizabilities,
             "Checks whether a force field of this region may carry a polarizability.")
        .def("add_force_field",
             &CEmbeddingRegion::add_force_field,
             "Adds a force field, or finds the one already held, and gets its identifier.",
             "force_field"_a)
        .def("add_molecule",
             py::overload_cast<const CMolecule &, const int>(&CEmbeddingRegion::add_molecule),
             "Adds a molecule of a force field already held.",
             "molecule"_a,
             "identifier"_a)
        .def("add_molecule",
             py::overload_cast<const CMolecule &, const CPolarizableForceField &>(&CEmbeddingRegion::add_molecule),
             "Adds a molecule and the force field which describes it.",
             "molecule"_a,
             "force_field"_a)
        .def("index_of_force_field",
             &CEmbeddingRegion::index_of_force_field,
             "Gets the index of the force field of this name, or minus one.",
             "name"_a)
        .def("number_of_molecules", &CEmbeddingRegion::number_of_molecules, "Gets the number of molecules.")
        .def("number_of_force_fields",
             &CEmbeddingRegion::number_of_force_fields,
             "Gets the number of kinds of molecule the region holds.")
        .def("get_molecule",
             &CEmbeddingRegion::get_molecule,
             "Gets a molecule.",
             "index"_a,
             py::return_value_policy::reference_internal)
        .def("get_identifier", &CEmbeddingRegion::get_identifier, "Gets the identifier of a molecule.", "index"_a)
        .def("get_force_field",
             &CEmbeddingRegion::get_force_field,
             "Gets a force field.",
             "index"_a,
             py::return_value_policy::reference_internal)
        .def("force_field_of",
             &CEmbeddingRegion::force_field_of,
             "Gets the force field of a molecule.",
             "index"_a,
             py::return_value_policy::reference_internal)
        .def("get_molecules", &CEmbeddingRegion::get_molecules, "Gets the molecules.")
        .def("get_identifiers",
             [](const CEmbeddingRegion &self) -> std::vector<int> { return self.get_identifiers(); },
             "Gets the identifiers, one for each molecule.")
        .def("get_force_fields", &CEmbeddingRegion::get_force_fields, "Gets the force fields.")
        .def("number_of_sites", &CEmbeddingRegion::number_of_sites, "Gets the number of sites the region carries.")
        .def("number_of_polarizable_sites",
             &CEmbeddingRegion::number_of_polarizable_sites,
             "Gets the number of sites which carry a polarizability.")
        .def("is_polarizable", &CEmbeddingRegion::is_polarizable, "Checks whether any molecule of the region is polarizable.")
        .def(
            "permanent_multipoles",
            [](const py::object &self, const int order) -> py::tuple {
                return _arrays_of(self.cast<const CEmbeddingRegion &>().permanent_multipoles(order), self);
            },
            "Gets the coordinates and the values of the permanent multipoles of an order.",
            "order"_a)
        .def("version", &CEmbeddingRegion::version, "Gets how many times the region has been added to.")
        .def("__len__", &CEmbeddingRegion::number_of_molecules);

    // CPolarizableEmbedding class

    PyClass<CPolarizableEmbedding>(m, "PolarizableEmbedding")
        .def(py::init<>())
        .def("get_polarizable_region",
             &CPolarizableEmbedding::get_polarizable_region,
             "Gets the polarizable region.",
             py::return_value_policy::reference_internal)
        .def("polarizable_region",
             &CPolarizableEmbedding::polarizable_region,
             "Gets the polarizable region, to add to.",
             py::return_value_policy::reference_internal)
        .def("get_nonpolarizable_region",
             &CPolarizableEmbedding::get_nonpolarizable_region,
             "Gets the nonpolarizable region.",
             py::return_value_policy::reference_internal)
        .def("nonpolarizable_region",
             &CPolarizableEmbedding::nonpolarizable_region,
             "Gets the nonpolarizable region, to add to.",
             py::return_value_policy::reference_internal)
        .def("number_of_molecules", &CPolarizableEmbedding::number_of_molecules, "Gets the number of molecules of both regions.")
        .def("number_of_sites", &CPolarizableEmbedding::number_of_sites, "Gets the number of sites of both regions.")
        .def("number_of_polarizable_sites",
             &CPolarizableEmbedding::number_of_polarizable_sites,
             "Gets the number of sites which carry a polarizability.")
        .def("is_polarizable", &CPolarizableEmbedding::is_polarizable, "Checks whether the environment polarizes at all.")
        .def(
            "permanent_multipoles",
            [](const py::object &self, const int order) -> py::tuple {
                return _arrays_of(self.cast<const CPolarizableEmbedding &>().permanent_multipoles(order), self);
            },
            "Gets the coordinates and the values of the permanent multipoles of an order, "
            "of both regions and the polarizable one first.",
            "order"_a)
        .def("permanent_nuclear_energy",
             &CPolarizableEmbedding::permanent_nuclear_energy,
             "Computes what the permanent charges of the environment do to the nuclei of the quantum region.",
             "molecule"_a,
             "basis"_a);
}

}  // namespace vlx_embedding
