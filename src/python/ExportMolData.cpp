//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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

#include "ExportMolData.hpp"

#include <mpi.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <string>
#include <vector>

#include "ChemicalElement.hpp"
#include "Codata.hpp"
#include "CoordinationNumber.hpp"
#include "DispersionModel.hpp"
#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "Molecule.hpp"
#include "PartialCharges.hpp"
#include "StringFormat.hpp"
#include "CommonNeighbors.hpp"
#include "AtomicRadii.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_moldata {  // vlx_moldata namespace

// Helper function for CMolecule constructor

static std::shared_ptr<CMolecule>
CMolecule_from_coords(const std::vector<std::string>& labels, const std::vector<double>& coords_raw, const std::string& units)
{
    // NOTE:
    // The C++ Molecule constructor expects the coordinates to be arranged as 3 x natoms,
    // namely {x1, x2, x3, x4, ..., y1, y2, y3, y4, ..., z1, z2, z3, z4, ...}

    // sanity check

    std::string errmol("Molecule: Inconsistent lengths of lists");

    errors::assertMsgCritical(coords_raw.size() == labels.size() * 3, errmol);

    // scaling factor

    auto scale = 1.0;

    if ((units.length() >= 3) && (fstr::upcase(units) == std::string("ANGSTROM").substr(0, units.length())))
    {
        scale = 1.0 / units::getBohrValueInAngstroms();
    }
    else if ((fstr::upcase(units) == std::string("AU")) || (fstr::upcase(units) == std::string("BOHR")))
    {
        scale = 1.0;
    }
    else
    {
        std::string errunit("Molecule: Invalid unit for coordinates");

        errors::assertMsgCritical(false, errunit);
    }

    std::vector<double> coords_au(coords_raw.size());

    for (size_t i = 0; i < coords_au.size(); i++)
    {
        coords_au[i] = coords_raw[i] * scale;
    }

    // form charge, mass, label and elemental ID vectors

    std::string errelm("Molecule: Unsupported chemical element");

    std::vector<double> charges;

    std::vector<double> masses;

    std::vector<int32_t> idselem;

    const int32_t natoms = static_cast<int32_t>(labels.size());

    for (int32_t i = 0; i < natoms; i++)
    {
        CChemicalElement chemelm;

        auto err = chemelm.setAtomType(fstr::upcase(labels[i]));

        errors::assertMsgCritical(err, errelm);

        charges.push_back(chemelm.getAtomicCharge());

        masses.push_back(chemelm.getAtomicMass());

        idselem.push_back(chemelm.getIdentifier());
    }

    // form molecule

    return std::make_shared<CMolecule>(coords_au, charges, masses, labels, idselem);
}

static std::shared_ptr<CMolecule>
CMolecule_from_array(const std::vector<std::string>&                labels,
                     const py::array_t<double, py::array::f_style>& py_coords,
                     const std::string&                             units = std::string("angstrom"))
{
    // NOTE:
    // The Python Molecule constructor expects the coordinates as a 2d numpy array,
    // namely np.array([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3], [x4, y4, z4], ...])

    // sanity check

    std::string errmol("Molecule: Inconsistent size");

    errors::assertMsgCritical(py_coords.shape(0) == static_cast<ssize_t>(labels.size()), errmol);

    errors::assertMsgCritical(py_coords.shape(1) == 3, errmol);

    // form coordinate vector

    std::vector<double> coords(py_coords.size());

    std::memcpy(coords.data(), py_coords.data(), py_coords.size() * sizeof(double));

    return CMolecule_from_coords(labels, coords, units);
}

static std::shared_ptr<CMolecule>
CMolecule_from_array_2(const std::vector<int32_t>& idselem, const py::array_t<double>& py_coords, const std::string& units = std::string("angstrom"))
{
    std::vector<std::string> labels;

    std::string errelm("Molecule: Unsupported element id");

    for (size_t i = 0; i < idselem.size(); i++)
    {
        CChemicalElement chemelm;

        auto err = chemelm.setAtomType(idselem[i]);

        errors::assertMsgCritical(err, errelm);

        labels.push_back(chemelm.getName());
    }

    return CMolecule_from_array(labels, py_coords, units);
}

// Helper function for checking multiplicity of molecule

static void
CMolecule_check_multiplicity(const CMolecule& self)
{
    auto multip = self.getMultiplicity() % 2;

    auto nelec = self.getNumberOfElectrons() % 2;

    bool flag = true;

    if ((multip == 0) && (nelec != 1)) flag = false;

    if ((multip == 1) && (nelec != 0)) flag = false;

    std::string errmult("Molecule.check_multiplicity: Incompatible multiplicity and number of electrons");

    errors::assertMsgCritical(flag, errmult);
}

// Helper function for checking proximity of atoms

static void
CMolecule_check_proximity(const CMolecule& self, const double minDistance)
{
    std::string errproxi("Molecule.check_proximity: Atoms too close");

    errors::assertMsgCritical(self.checkProximity(minDistance), errproxi);
}

// Exports classes/functions in src/moldata to python

void
export_moldata(py::module& m)
{
    // CMolecule class

    PyClass<CMolecule>(m, "Molecule")
        .def(py::init<>())
        .def(py::init<const CMolecule&>())
        .def(py::init<const CMolecule&, const CMolecule&>())
        .def(py::init(&CMolecule_from_array), "symbols"_a, "coordinates"_a, "units"_a = std::string("angstrom"))
        .def(py::init(&CMolecule_from_array_2), "Zs"_a, "coordinates"_a, "units"_a = std::string("angstrom"))
        .def("get_string", &CMolecule::printGeometry, "Prints geometry of molecule as table to output stream.")
        .def("set_charge", &CMolecule::setCharge, "Sets charge of molecule object.", "charge"_a)
        .def("get_charge", &CMolecule::getCharge, "Gets charge of molecule.")
        .def("set_multiplicity", &CMolecule::setMultiplicity, "Sets spin multiplicity of molecule object.", "multiplicity"_a)
        .def("get_multiplicity", &CMolecule::getMultiplicity, "Gets spin multiplicity of molecule.")
        .def("check_multiplicity", &CMolecule_check_multiplicity, "Checks multiplicity of molecule.")
        .def("check_proximity", &CMolecule_check_proximity, "Checks proximity of atoms.", "minDistance"_a)
        .def("get_elemental_composition", &CMolecule::getElementalComposition, "Gets set of unique chemical elements in molecule.")
        .def("nuclear_repulsion_energy",
             &CMolecule::getNuclearRepulsionEnergy,
             "Gets nuclear repulsion energy for molecule assuming point charge model for nucleus.")
        .def("get_sub_molecule",
             &CMolecule::getSubMolecule,
             "Creates a sub-molecule object by slicing the molecule object.",
             "startIndex"_a,
             "numAtoms"_a)
        .def("number_of_atoms", vlx_general::overload_cast_<>()(&CMolecule::getNumberOfAtoms, py::const_), "Gets total number of atoms in molecule.")
        .def("number_of_atoms",
             vlx_general::overload_cast_<const int32_t>()(&CMolecule::getNumberOfAtoms, py::const_),
             "Gets number of atoms belonging to specific chemical element in molecule.",
             "idElemental"_a)
        .def("number_of_atoms",
             vlx_general::overload_cast_<const int32_t, const int32_t, const int32_t>()(&CMolecule::getNumberOfAtoms, py::const_),
             "Gets number of atoms belonging to specific chemical element in list of atoms in molecule.",
             "iAtom"_a,
             "nAtoms"_a,
             "idElemental"_a)
        .def("number_of_electrons", &CMolecule::getNumberOfElectrons, "Gets a number of electrons in molecule.")
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
        .def(
            "x_to_numpy",
            [](const CMolecule& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoordinatesX(), self.getNumberOfAtoms());
            },
            "Gets X coordinates as numpy array.")
        .def(
            "y_to_numpy",
            [](const CMolecule& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoordinatesY(), self.getNumberOfAtoms());
            },
            "Gets Y coordinates as numpy array.")
        .def(
            "z_to_numpy",
            [](const CMolecule& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoordinatesZ(), self.getNumberOfAtoms());
            },
            "Gets Z coordinates as numpy array.")
        .def(
            "partial_charges",
            [](const CMolecule& self) -> py::array_t<double> {
                auto chg = parchg::getPartialCharges(self, self.getCharge());
                return vlx_general::pointer_to_numpy(chg.data(), static_cast<int32_t>(chg.size()));
            },
            "Gets partial charges for molecule.")
        .def(
            "vdw_radii_to_numpy",
            [](const CMolecule& self) -> py::array_t<double> {
                auto atomradii = self.getVdwRadii();
                return vlx_general::pointer_to_numpy(atomradii.data(), static_cast<int32_t>(atomradii.size()));
            },
            "Gets VDW radii for molecule.")
        .def(
            "mk_radii_to_numpy",
            [](const CMolecule& self) ->py::array_t<double> {
                auto atomradii = self.getMkRadii();
                return vlx_general::pointer_to_numpy(atomradii.data(), static_cast<int32_t>(atomradii.size()));
            },
            "Gets MK radii for molecule.")
        .def(
            "covalent_radii_to_numpy",
            [](const CMolecule& self) ->py::array_t<double> {
                auto atomradii = self.getCovalentRadii();
                return vlx_general::pointer_to_numpy(atomradii.data(), static_cast<int32_t>(atomradii.size()));
            },
            "Gets covalent radii for molecule.")
        .def(
            "elem_ids_to_numpy",
            [](const CMolecule& self) -> py::array_t<int32_t> {
                return vlx_general::pointer_to_numpy(self.getIdsElemental(), self.getNumberOfAtoms());
            },
            "Gets nuclear charges for molecule.")
        .def(
            "masses_to_numpy",
            [](const CMolecule& self) -> py::array_t<double> {
                auto masses = self.getMasses();
                return vlx_general::pointer_to_numpy(masses.data(), masses.size());
            },
            "Gets masses for molecule.")
        .def(
            "broadcast",
            [](CMolecule& self, int32_t rank, py::object py_comm) -> void {
                auto comm = vlx_general::get_mpi_comm(py_comm);
                self.broadcast(rank, *comm);
            },
            "Broadcasts Molecule object.",
            "rank"_a,
            "py_comm"_a)
        .def(py::self == py::self);

    // CChemicalElement class

    PyClass<CChemicalElement>(m, "ChemicalElement")
        .def(py::init<>())
        .def("set_atom_type",
             vlx_general::overload_cast_<const std::string&>()(&CChemicalElement::setAtomType),
             "Sets chemical element properties using name of chemical element.",
             "atomLabel"_a)
        .def("set_atom_type",
             vlx_general::overload_cast_<const int32_t>()(&CChemicalElement::setAtomType),
             "Sets chemical element object properties using chemical element number.",
             "idElemental"_a)
        .def("get_name", &CChemicalElement::getName, "Gets name of chemical element.")
        .def(py::self == py::self);

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
    
    // CCommonNeighbors class

    PyClass<CCommonNeighbors>(m, "CommonNeighbors")
        .def(py::init<>())
        .def(py::init<const CMolecule&, const double>())
        .def("generate", &CCommonNeighbors::generate)
        .def("comp_cna", &CCommonNeighbors::compJaccardIndex)
        .def("__repr__", &CCommonNeighbors::getSignaturesRepr);
}

}  // namespace vlx_moldata
