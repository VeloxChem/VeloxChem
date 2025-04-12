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

#include "ExportMoldata.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <algorithm>

#include "AtomicPartialChargesModel.hpp"
#include "ChemicalElement.hpp"
#include "Codata.hpp"
#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "Molecule.hpp"
#include "StringFormat.hpp"

namespace vlx_moldata {  // vlx_moldata namespace

// Exports classes/functions in src/moldata to python

void
export_moldata(py::module& m)
{
    // CChemicalElement class

    PyClass<CChemicalElement>(m, "ChemicalElement")
        .def(py::init<>())
        .def("set_atom_type",
             py::overload_cast<const std::string&>(&CChemicalElement::setAtomType),
             "Sets chemical element properties using name of chemical element.")
        .def("set_atom_type",
             py::overload_cast<const int64_t>(&CChemicalElement::setAtomType),
             "Sets chemical element properties using chemical element number.")
        .def("set_isotope", &CChemicalElement::setIsotope, "Sets chemical element isotope using isotope number.")
        .def("get_name", &CChemicalElement::getName, "Gets name of chemical element.")
        .def("get_identifier", &CChemicalElement::getIdentifier, "Gets identifier of chemical element.")
        .def("get_mass", &CChemicalElement::getAtomicMass, "Gets atomic mass of chemical element.")
        .def("get_charge", &CChemicalElement::getAtomicCharge, "Gets atomic charge of chemical element.")
        .def("get_max_angular_momentum",
             &CChemicalElement::getMaxAngularMomentum,
             "Gets maximum angular momentum of occupied electron shell in chemical element.")
        .def("get_max_identifier", &CChemicalElement::getMaxIdentifier, "Gets maximum elemental number of supported chemical elements.");

    // CMolecule class

    PyClass<CMolecule>(m, "Molecule")
        .def(py::init<>())
        .def(py::init<const std::vector<int64_t>&, const std::vector<TPoint3D>&, const std::string&, const std::vector<std::string>&>())
        .def(py::init<const std::vector<std::string>&, const std::vector<TPoint3D>&, const std::string&, const std::vector<std::string>&>())
        .def(py::init<const CMolecule&>())
        .def(py::init<const CMolecule&, const CMolecule&>())
        .def("add_atom", py::overload_cast<const std::string&, const TPoint3D&, const std::string&, const std::string&>(&CMolecule::addAtom), "Adds atom to molecule.")
        .def("add_atom", py::overload_cast<const int64_t, const TPoint3D&, const std::string&, const std::string&>(&CMolecule::addAtom), "Adds atom to molecule.")
        .def("set_charge", &CMolecule::setCharge, "Sets charge of molecule.")
        .def("set_multiplicity", &CMolecule::setMultiplicity, "Sets spin multiplicity of molecule.")
        .def(
            "check_multiplicity",
            [](const CMolecule& self) -> void {
                auto multip = self.getMultiplicity() % 2;
                auto nelec  = self.getNumberOfElectrons() % 2;
                bool flag   = true;
                if ((multip == 0) && (nelec != 1)) flag = false;
                if ((multip == 1) && (nelec != 0)) flag = false;
                std::string errmult("Molecule.check_multiplicity: Incompatible multiplicity and number of electrons");
                errors::assertMsgCritical(flag, errmult);
            },
            "Checks spin multiplicity of molecule.")
        .def("get_charge", &CMolecule::getCharge, "Gets charge of molecule.")
        .def("get_multiplicity", &CMolecule::getMultiplicity, "Gets spin multiplicity of molecule.")
        .def("number_of_atoms", py::overload_cast<>(&CMolecule::getNumberOfAtoms, py::const_), "Gets total number of atoms in molecule.")
        .def("number_of_atoms",
             py::overload_cast<const int64_t>(&CMolecule::getNumberOfAtoms, py::const_),
             "Gets number of atoms belonging to specific chemical element in molecule.")
        .def("number_of_atoms",
             py::overload_cast<const int64_t, const int64_t, const int64_t>(&CMolecule::getNumberOfAtoms, py::const_),
             "Gets number of atoms belonging to specific chemical element in list of atoms in molecule.")
        .def("get_elemental_composition", &CMolecule::getElementalComposition, "Gets set of unique chemical elements in molecule.")
        .def("number_of_electrons", &CMolecule::getNumberOfElectrons, "Gets a number of electrons in molecule.")
        .def("get_identifiers", &CMolecule::getIdsElemental, "Gets a vector of elemental identidiers in molecule.")
        .def("get_coordinates", &CMolecule::getCoordinates, py::arg("units") = "au", "Gets coordinates of atoms in molecules")
        .def("get_charges", &CMolecule::getCharges, "Gets a vector of atomic charges in molecule.")
        .def(
            "get_element_ids",
            [](const CMolecule& self) -> py::array_t<double> {
                const auto elem_ids = self.getCharges();
                const auto n_atoms  = self.getNumberOfAtoms();
                return vlx_general::pointer_to_numpy(elem_ids.data(), {n_atoms});
            },
            "Gets nuclear charges for molecule.")
        .def(
            "get_partial_charges",
            [](const CMolecule& self, const double net_charge) -> std::vector<double> { return atmparchg::getPartialCharges(self, net_charge); },
            "Gets partial charges for molecule.")
        .def("get_masses", &CMolecule::getMasses, "Gets a vector of atomic masses in molecule.")
        .def("get_labels", &CMolecule::getLabels, "Gets a vector of atomic labels in molecule.")
        .def("get_label", &CMolecule::getLabel, "Gets an atomic labels of specific atom in molecule.")
        .def("get_basis_set_labels", &CMolecule::getBasisSetLabels, "Gets a vector of names of atom basis sets in molecule.")
        .def("get_atom_basis_set_label", &CMolecule::getAtomBasisSetLabel, "Gets basis set name of specific atom in molecule.")
        .def("get_atom_coordinates", &CMolecule::getAtomCoordinates, "Gets coordinates [x,y,z] of atom.")
        .def("atom_indexes", &CMolecule::getAtomIndexes, "Gets indexes of atoms with requested atomic label")
        .def(
            "number_of_alpha_electrons",
            [](const CMolecule& self) -> int32_t {
                int32_t nelec  = self.getNumberOfElectrons();
                int32_t mult_1 = self.getMultiplicity() - 1;
                return (nelec + mult_1) / 2;
            },
            "Gets number of alpha electrons.")
        .def(
            "number_of_beta_electrons",
            [](const CMolecule& self) -> int32_t {
                int32_t nelec  = self.getNumberOfElectrons();
                int32_t mult_1 = self.getMultiplicity() - 1;
                return (nelec - mult_1) / 2;
            },
            "Gets number of beta electrons.")
        .def("nuclear_repulsion_energy",
             &CMolecule::getNuclearRepulsionEnergy,
             "Gets nuclear repulsion energy for molecule assuming point charge model for nucleus.")
        .def("check_proximity", &CMolecule::checkProximity, "Checks if proximity requirement is satisfied by all pairs of atoms in molecule..")
        .def("get_string", &CMolecule::printGeometry, "Creates string representation of molecule.");
}

}  // namespace vlx_moldata
