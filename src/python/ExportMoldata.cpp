#include "ExportMoldata.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "ChemicalElement.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

namespace vlx_moldata {

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
    m.def("chemical_element_max_identifier",
          &chem_elem::max_identifier,
          "Gets maximum value of chemical element number.");

    // CMolecule class
    PyClass<CMolecule>(m, "Molecule")
        .def(py::init<>())
        .def(py::init<const std::vector<int> &, const std::vector<TPoint<double>> &, const std::string &>())
        .def(py::init<const CMolecule &>())
        .def(py::init<const CMolecule &, const CMolecule &>())
        .def(py::pickle(
            [](const CMolecule &mol) {
                return py::make_tuple(
                    mol.identifiers(), mol.coordinates("au"), mol.get_charge(), mol.get_multiplicity());
            },
            [](py::tuple t) {
                auto mol = CMolecule(
                    t[0].cast<std::vector<int>>(), t[1].cast<std::vector<TPoint<double>>>(), std::string("au"));
                mol.set_charge(t[2].cast<double>());
                mol.set_multiplicity(t[3].cast<int>());
                return mol;
            }))
        .def("add_atom", &CMolecule::add_atom, "Adds atom to molecule.")
        .def("slice", &CMolecule::slice, "Creates a new molecule by slicing selected atoms from molecule.")
        .def("set_charge", &CMolecule::set_charge, "Sets charge of molecule.")
        .def("set_multiplicity", &CMolecule::set_multiplicity, "Sets spin multiplicity of molecule.")
        .def("get_charge", &CMolecule::get_charge, "Gets charge of molecule.")
        .def("get_multiplicity", &CMolecule::get_multiplicity, "Gets spin multiplicity of molecule.")
        .def("number_of_atoms",
             py::overload_cast<>(&CMolecule::number_of_atoms, py::const_),
             "Gets total number of atoms in molecule.")
        .def("number_of_atoms",
             py::overload_cast<const int>(&CMolecule::number_of_atoms, py::const_),
             "Gets number of atoms belonging to specific chemical element in "
             "molecule.")
        .def("number_of_atoms",
             py::overload_cast<const int, const int, const int>(&CMolecule::number_of_atoms, py::const_),
             "Gets number of atoms belonging to specific chemical element in "
             "range of atoms in molecule.")
        .def("get_elemental_composition",
             &CMolecule::elemental_composition,
             "Gets set of unique chemical elements in molecule.")
        .def("number_of_electrons", &CMolecule::number_of_electrons, "Gets a number of electrons in molecule.")
        .def("get_identifiers", &CMolecule::identifiers, "Gets a vector of elemental identidiers in molecule.")
        .def("get_coordinates",
             &CMolecule::coordinates,
             py::arg("units") = "au",
             "Gets coordinates of atoms in molecules")
        .def("get_charges", &CMolecule::charges, "Gets a vector of atomic charges in molecule.")
        .def("get_masses", &CMolecule::masses, "Gets a vector of atomic masses in molecule.")
        .def("get_labels", &CMolecule::labels, "Gets a vector of atomic labels in molecule.")
        .def("get_label", &CMolecule::label, "Gets an atomic labels of specific atom in molecule.")
        .def("get_atom_coordinates", &CMolecule::atom_coordinates, "Gets coordinates [x,y,z] of atom.")
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

    // ... add class bindings here ...
}

}  // namespace vlx_moldata
