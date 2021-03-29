//
//                           VELOXCHEM 1.0-RC
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

#include "CoordinationNumber.hpp"
#include "ChemicalElement.hpp"
#include "Codata.hpp"
#include "DispersionModel.hpp"
#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "Molecule.hpp"
#include "PartialCharges.hpp"
#include "StringFormat.hpp"
#include "VdwRadii.hpp"

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

// Helper function for getting number of alpha/beta electrons

static int32_t
CMolecule_alpha_elec(const CMolecule& self)
{
    int32_t nelec = self.getNumberOfElectrons();

    int32_t mult_1 = self.getMultiplicity() - 1;

    return (nelec + mult_1) / 2;
}

static int32_t
CMolecule_beta_elec(const CMolecule& self)
{
    int32_t nelec = self.getNumberOfElectrons();

    int32_t mult_1 = self.getMultiplicity() - 1;

    return (nelec - mult_1) / 2;
}

// Helper function for getting coordinates as numpy array

static py::array_t<double>
CMolecule_x_to_numpy(const CMolecule& self)
{
    return vlx_general::pointer_to_numpy(self.getCoordinatesX(), self.getNumberOfAtoms());
}

static py::array_t<double>
CMolecule_y_to_numpy(const CMolecule& self)
{
    return vlx_general::pointer_to_numpy(self.getCoordinatesY(), self.getNumberOfAtoms());
}

static py::array_t<double>
CMolecule_z_to_numpy(const CMolecule& self)
{
    return vlx_general::pointer_to_numpy(self.getCoordinatesZ(), self.getNumberOfAtoms());
}

// Helper function for getting coodination number for molecule

static py::array_t<double>
CMolecule_coordination_numbers(const CMolecule& self)
{
    auto cn = coordnum::getCoordinationNumber(self);

    return vlx_general::pointer_to_numpy(cn.data(), static_cast<int32_t>(cn.size()));
}

// Helper function for getting partial charges for molecule

static py::array_t<double>
CMolecule_partial_charges(const CMolecule& self)
{
    auto chg = parchg::getPartialCharges(self, self.getCharge());

    return vlx_general::pointer_to_numpy(chg.data(), static_cast<int32_t>(chg.size()));
}

// Helper function for getting VDW radii for molecule

static py::array_t<double>
CMolecule_vdw_radii_to_numpy(const CMolecule& self)
{
    auto atomradii = self.getVdwRadii();

    return vlx_general::pointer_to_numpy(atomradii.data(), static_cast<int32_t>(atomradii.size()));
}

// Helper function for getting nuclear charges for molecule

static py::array_t<int32_t>
CMolecule_elem_ids_to_numpy(const CMolecule& self)
{
    return vlx_general::pointer_to_numpy(self.getIdsElemental(), self.getNumberOfAtoms());
}

// Helper function for getting masses for molecule

static py::array_t<double>
CMolecule_masses_to_numpy(const CMolecule& self)
{
    auto masses = self.getMasses();

    return vlx_general::pointer_to_numpy(masses.data(), masses.size());
}

// Helper function for getting elemental composition

static py::list
CMolecule_get_elem_comp(const CMolecule& self)
{
    py::list elemcomp;

    auto elmlist = self.getElementalComposition();

    for (auto p = elmlist.cbegin(); p != elmlist.cend(); ++p)
    {
        elemcomp.append(*p);
    }

    return elemcomp;
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

// Helper function for broadcasting CMolecule object

static void
CMolecule_broadcast(CMolecule& self, int32_t rank, py::object py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    self.broadcast(rank, *comm_ptr);
}

// Exports classes/functions in src/moldata to python

void
export_moldata(py::module& m)
{
    // CMolecule class

    py::class_<CMolecule, std::shared_ptr<CMolecule>>(m, "Molecule")
        .def(py::init<>())
        .def(py::init<const CMolecule&>())
        .def(py::init<const CMolecule&, const CMolecule&>())
        .def(py::init(&CMolecule_from_array), "symbols"_a, "coordinates"_a, "units"_a = std::string("angstrom"))
        .def(py::init(&CMolecule_from_array_2), "Zs"_a, "coordinates"_a, "units"_a = std::string("angstrom"))
        .def("set_charge", &CMolecule::setCharge)
        .def("get_charge", &CMolecule::getCharge)
        .def("set_multiplicity", &CMolecule::setMultiplicity)
        .def("get_multiplicity", &CMolecule::getMultiplicity)
        .def("check_multiplicity", &CMolecule_check_multiplicity)
        .def("get_string", &CMolecule::printGeometry)
        .def("check_proximity", &CMolecule_check_proximity)
        .def("get_sub_molecule", &CMolecule::getSubMolecule)
        .def("number_of_atoms", vlx_general::overload_cast_<>()(&CMolecule::getNumberOfAtoms, py::const_))
        .def("number_of_atoms", vlx_general::overload_cast_<const int32_t>()(&CMolecule::getNumberOfAtoms, py::const_), "idElemental"_a)
        .def("number_of_atoms",
             vlx_general::overload_cast_<const int32_t, const int32_t, const int32_t>()(&CMolecule::getNumberOfAtoms, py::const_),
             "iatom"_a,
             "natoms"_a,
             "idElemental"_a)
        .def("number_of_electrons", &CMolecule::getNumberOfElectrons)
        .def("number_of_alpha_electrons", &CMolecule_alpha_elec)
        .def("number_of_beta_electrons", &CMolecule_beta_elec)
        .def("nuclear_repulsion_energy", &CMolecule::getNuclearRepulsionEnergy)
        .def("x_to_numpy", &CMolecule_x_to_numpy)
        .def("y_to_numpy", &CMolecule_y_to_numpy)
        .def("z_to_numpy", &CMolecule_z_to_numpy)
        .def("coordination_numbers", &CMolecule_coordination_numbers)
        .def("partial_charges", &CMolecule_partial_charges)
        .def("vdw_radii_to_numpy", &CMolecule_vdw_radii_to_numpy)
        .def("elem_ids_to_numpy", &CMolecule_elem_ids_to_numpy)
        .def("masses_to_numpy", &CMolecule_masses_to_numpy)
        .def("get_elemental_composition", &CMolecule_get_elem_comp)
        .def("broadcast", &CMolecule_broadcast)
        .def(py::self == py::self);

    // CChemicalElement class

    py::class_<CChemicalElement, std::shared_ptr<CChemicalElement>>(m, "ChemicalElement")
        .def(py::init<>())
        .def("set_atom_type", vlx_general::overload_cast_<const std::string&>()(&CChemicalElement::setAtomType), "label"_a)
        .def("set_atom_type", vlx_general::overload_cast_<const int32_t>()(&CChemicalElement::setAtomType), "idElemental"_a)
        .def("get_name", &CChemicalElement::getName)
        .def(py::self == py::self);

    // CDispersionModel class

    py::class_<CDispersionModel, std::shared_ptr<CDispersionModel>>(m, "DispersionModel")
        .def(py::init<>())
        .def("compute", &CDispersionModel::compute)
        .def("get_energy", &CDispersionModel::getEnergy)
        .def("get_gradient", &CDispersionModel::getGradient);
}

}  // namespace vlx_moldata
