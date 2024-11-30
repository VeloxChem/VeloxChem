//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "ExportMoldata.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "ChemicalElement.hpp"
#include "Codata.hpp"
#include "DispersionModel.hpp"
#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "Molecule.hpp"
#include "PartialCharges.hpp"
#include "Point.hpp"
#include "StringFormat.hpp"

namespace vlx_moldata {

static std::shared_ptr<CMolecule>
CMolecule_from_labels_and_array(const std::vector<std::string>& labels,
                                const py::array_t<double>&      py_coords,
                                const std::string&              units,
                                const std::vector<std::pair<std::string, std::string>>& atom_basis_labels)
{
    // sanity check

    std::string errmol("Molecule: Inconsistent sizes");

    errors::assertMsgCritical(py_coords.shape(0) == static_cast<py::ssize_t>(labels.size()), errmol);

    errors::assertMsgCritical(py_coords.shape(0) == static_cast<py::ssize_t>(atom_basis_labels.size()), errmol);

    errors::assertMsgCritical(py_coords.shape(1) == 3, errmol);

    std::string errstyle("Molecule: Expecting C-style contiguous array");

    auto c_style = py::detail::check_flags(py_coords.ptr(), py::array::c_style);

    errors::assertMsgCritical(c_style, errstyle);

    // scaling factor

    auto scale = 1.0;

    if ((units.length() >= 3) && (format::upper_case(units) == std::string("ANGSTROM").substr(0, units.length())))
    {
        scale = 1.0 / units::bohr_in_angstrom();
    }
    else if ((format::upper_case(units) == std::string("AU")) || (format::upper_case(units) == std::string("BOHR")))
    {
        scale = 1.0;
    }
    else
    {
        std::string errunit("Molecule: Invalid unit for coordinates");

        errors::assertMsgCritical(false, errunit);
    }

    // coordinates

    std::vector<TPoint<double>> coords_au;

    for (size_t i = 0; i < labels.size(); i++)
    {
        TPoint<double> pxyz({
            py_coords.data()[i * 3 + 0] * scale,
            py_coords.data()[i * 3 + 1] * scale,
            py_coords.data()[i * 3 + 2] * scale
        });

        coords_au.push_back(pxyz);
    }

    // elemental IDs

    std::vector<int> elem_ids;

    for (const auto& label : labels)
    {
        auto elem_id = chem_elem::identifier(format::upper_case(label));

        errors::assertMsgCritical(elem_id != -1, std::string("Molecule: Invalid element label"));

        elem_ids.push_back(elem_id);
    }

    // molecule

    return std::make_shared<CMolecule>(elem_ids, coords_au, std::string("BOHR"), atom_basis_labels);
}

static std::shared_ptr<CMolecule>
CMolecule_from_labels_and_array_2(const std::vector<std::string>& labels,
                                  const py::array_t<double>&      py_coords,
                                  const std::string&              units = std::string("angstrom"))
{
    const auto natoms = static_cast<int>(labels.size());

    std::vector<std::pair<std::string, std::string>> atom_basis_labels;

    for (int a = 0; a < natoms; a++)
    {
        atom_basis_labels.push_back(std::make_pair(std::string(""), format::upper_case(labels[a])));
    }

    return CMolecule_from_labels_and_array(labels, py_coords, units, atom_basis_labels);
}

static std::shared_ptr<CMolecule>
CMolecule_from_elemids_and_array(const std::vector<int>& elem_ids,
                                 const py::array_t<double>& py_coords,
                                 const std::string& units,
                                 const std::vector<std::pair<std::string, std::string>>& atom_basis_labels)
{
    std::vector<std::string> labels;

    for (const auto elem_id : elem_ids)
    {
        const auto label = chem_elem::label(elem_id);

        labels.push_back(label);
    }

    return CMolecule_from_labels_and_array(labels, py_coords, units, atom_basis_labels);
}

static std::shared_ptr<CMolecule>
CMolecule_from_elemids_and_array_2(const std::vector<int>& elem_ids, const py::array_t<double>& py_coords, const std::string& units = std::string("angstrom"))
{
    const auto natoms = static_cast<int>(elem_ids.size());

    std::vector<std::pair<std::string, std::string>> atom_basis_labels;

    for (int a = 0; a < natoms; a++)
    {
        const auto label = chem_elem::label(elem_ids[a]);

        atom_basis_labels.push_back(std::make_pair(std::string(""), format::upper_case(label)));
    }

    return CMolecule_from_elemids_and_array(elem_ids, py_coords, units, atom_basis_labels);
}

void
export_moldata(py::module &m)
{
    // exposing functions from ChemicalElement.hpp
    m.def("is_chemical_element", &chem_elem::valid_identifier, "Checks if identifier is chemical element number.");
    m.def("chemical_element_name", &chem_elem::name, "Gets chemical element name.");
    m.def("chemical_element_label", &chem_elem::label, "Gets chemical element label.");
    m.def("chemical_element_identifier", &chem_elem::identifier, "Gets chemical element identifier.");
    m.def("chemical_element_mass", &chem_elem::mass, "Gets chemical element mass.");
    m.def("chemical_element_max_angular_momentum",
          &chem_elem::max_angular_momentum,
          "Gets maximum angular momentum of atomic shell in chemical element.");
    m.def("chemical_element_max_identifier", &chem_elem::max_identifier, "Gets maximum value of chemical element number.");

    // CMolecule class
    PyClass<CMolecule>(m, "Molecule")
        .def(py::init<>())
        .def(py::init<const std::vector<int> &, const std::vector<TPoint<double>> &, const std::string &>())
        .def(py::init<const std::vector<int>&, const std::vector<TPoint<double>>&, const std::string&, const std::vector<std::pair<std::string, std::string>>&>())
        .def(py::init<const CMolecule &>())
        .def(py::init<const CMolecule &, const CMolecule &>())
        .def(py::init(&CMolecule_from_labels_and_array), "symbols"_a, "coordinates"_a, "units"_a, "atom_basis_labels"_a)
        .def(py::init(&CMolecule_from_labels_and_array_2), "symbols"_a, "coordinates"_a, "units"_a = std::string("angstrom"))
        .def(py::init(&CMolecule_from_elemids_and_array), "Zs"_a, "coordinates"_a, "units"_a, "atom_basis_labels"_a)
        .def(py::init(&CMolecule_from_elemids_and_array_2), "Zs"_a, "coordinates"_a, "units"_a = std::string("angstrom"))
        .def(py::pickle(
            [](const CMolecule &mol) { return py::make_tuple(mol.identifiers(), mol.coordinates("au"), mol.atom_basis_labels(), mol.get_charge(), mol.get_multiplicity()); },
            [](py::tuple t) {
                auto mol = CMolecule(t[0].cast<std::vector<int>>(), t[1].cast<std::vector<TPoint<double>>>(), std::string("au"), t[2].cast<std::vector<std::pair<std::string, std::string>>>());
                mol.set_charge(t[3].cast<double>());
                mol.set_multiplicity(t[4].cast<int>());
                return mol;
            }))
        .def("add_atom", py::overload_cast<const int, const TPoint<double>&, const std::string&>(&CMolecule::add_atom),
          "Adds atom to molecule.")
        .def("add_atom", py::overload_cast<const int, const TPoint<double>&, const std::string&, const std::pair<std::string, std::string>&>(&CMolecule::add_atom),
          "Adds atom to molecule.")
        .def("slice", &CMolecule::slice, "Creates a new molecule by slicing selected atoms from molecule.")
        .def("set_charge", &CMolecule::set_charge, "Sets charge of molecule.")
        .def("set_multiplicity", &CMolecule::set_multiplicity, "Sets spin multiplicity of molecule.")
        .def("get_charge", &CMolecule::get_charge, "Gets charge of molecule.")
        .def("get_multiplicity", &CMolecule::get_multiplicity, "Gets spin multiplicity of molecule.")
        .def("number_of_atoms", py::overload_cast<>(&CMolecule::number_of_atoms, py::const_), "Gets total number of atoms in molecule.")
        .def("number_of_atoms",
             py::overload_cast<const int>(&CMolecule::number_of_atoms, py::const_),
             "Gets number of atoms belonging to specific chemical element in "
             "molecule.")
        .def("number_of_atoms",
             py::overload_cast<const int, const int, const int>(&CMolecule::number_of_atoms, py::const_),
             "Gets number of atoms belonging to specific chemical element in "
             "range of atoms in molecule.")
        .def("get_elemental_composition", &CMolecule::elemental_composition, "Gets set of unique chemical elements in molecule.")
        .def("number_of_electrons", &CMolecule::number_of_electrons, "Gets a number of electrons in molecule.")
        .def("get_identifiers", &CMolecule::identifiers, "Gets a vector of elemental identidiers in molecule.")
        .def("get_atom_basis_labels", &CMolecule::atom_basis_labels, "Gets a vector of atom basis set labels.")
        .def("get_coordinates", &CMolecule::coordinates, py::arg("units") = "au", "Gets coordinates of atoms in molecules")
        .def("get_charges", &CMolecule::charges, "Gets a vector of atomic charges in molecule.")
        .def(
            "get_element_ids",
            [](const CMolecule& self) -> py::array_t<double> {
                const auto elem_ids = self.charges();
                const auto n_atoms  = self.number_of_atoms();
                return vlx_general::pointer_to_numpy(elem_ids.data(), {n_atoms});
            },
            "Gets nuclear charges for molecule.")
        .def(
            "get_partial_charges",
            [](const CMolecule& self, const double net_charge) -> std::vector<double> { return parchg::getPartialCharges(self, net_charge); },
            "Gets partial charges for molecule.")
        .def(
            "vdw_radii_to_numpy",
            [](const CMolecule& self) -> py::array_t<double> {
                auto atomradii = self.get_vdw_radii();
                return vlx_general::pointer_to_numpy(atomradii.data(), {static_cast<int>(atomradii.size())});
            },
            "Gets VDW radii for molecule.")
        .def(
            "mk_radii_to_numpy",
            [](const CMolecule& self) -> py::array_t<double> {
                auto atomradii = self.get_mk_radii();
                return vlx_general::pointer_to_numpy(atomradii.data(), {static_cast<int>(atomradii.size())});
            },
            "Gets MK radii for molecule.")
        .def(
            "chelpg_radii_to_numpy",
            [](const CMolecule& self) -> py::array_t<double> {
                auto atomradii = self.get_chelpg_radii();
                return vlx_general::pointer_to_numpy(atomradii.data(), {static_cast<int>(atomradii.size())});
            },
            "Gets CHELPG radii for molecule.")
        .def(
            "covalent_radii_to_numpy",
            [](const CMolecule& self) -> py::array_t<double> {
                auto atomradii = self.get_covalent_radii();
                return vlx_general::pointer_to_numpy(atomradii.data(), {static_cast<int>(atomradii.size())});
            },
            "Gets covalent radii for molecule.")
        .def("get_masses", &CMolecule::masses, "Gets a vector of atomic masses in molecule.")
        .def("get_labels", &CMolecule::labels, "Gets a vector of atomic labels in molecule.")
        .def("get_label", &CMolecule::label, "Gets an atomic labels of specific atom in molecule.")
        .def(
            "get_atom_coordinates",
            [](const CMolecule& self, const int iatom, const std::string& unit) -> std::array<double, 3> {
                const auto atom_coords = self.atom_coordinates(iatom, unit);
                return atom_coords.coordinates();
            },
            "Gets coordinates of a given atom.",
             "iatom"_a,
             "unit"_a = std::string("au"))
        .def("set_atom_coordinates", &CMolecule::set_atom_coordinates, "Sets coordinates (x,y,z) of atom.", "iatom"_a, "xyz"_a)
        .def("atom_indices", &CMolecule::atom_indices, "Gets indices of atoms with requested atomic label.")
        .def("nuclear_repulsion_energy",
             &CMolecule::nuclear_repulsion_energy,
             "Gets nuclear repulsion energy for molecule assuming point charge "
             "model for nucleus.")
        .def("check_proximity",
             &CMolecule::check_proximity,
             "Checks if proximity requirement is satisfied by all pairs of atoms "
             "in molecule.")
        .def("__eq__", [](const CMolecule &self, const CMolecule &other) { return self == other; })
        .def("__copy__", [](const CMolecule &self) { return CMolecule(self); })
        .def("__deepcopy__", [](const CMolecule &self, py::dict) { return CMolecule(self); });

    // CDispersionModel class

    PyClass<CDispersionModel>(m, "DispersionModel")
        .def(py::init<>())
        .def("compute",
             &CDispersionModel::compute,
             "Computes dispersion energy and gradient for a given molecule and a given density functional.",
             "molecule"_a,
             "xcLabel"_a)
        .def("get_energy", &CDispersionModel::getEnergy, "Gets dispersion energy.")
        .def("get_gradient", &CDispersionModel::getGradient, "Gets dispersion gradient.");
}

}  // namespace vlx_moldata
