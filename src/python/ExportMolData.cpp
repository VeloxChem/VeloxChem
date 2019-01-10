//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <mpi.h>
#include <memory>
#include <vector>
#include <string>

#include "Molecule.hpp"
#include "VdwRadii.hpp"
#include "ErrorHandler.hpp"
#include "StringFormat.hpp"
#include "ChemicalElement.hpp"
#include "ExportGeneral.hpp"
#include "ExportMolData.hpp"

namespace py = pybind11;

namespace vlx_moldata { // vlx_moldata namespace

// Helper function for CMolecule constructor

static std::shared_ptr<CMolecule>
CMolecule_from_xyz(const std::vector<std::string>& labels,
                   const std::vector<double>& x_coords,
                   const std::vector<double>& y_coords,
                   const std::vector<double>& z_coords)
{
    // form coordinate vector

    std::string errmol("Molecule.from_xyz: Inconsistent lengths of lists");

    errors::assertMsgCritical(x_coords.size() == labels.size(), errmol);

    errors::assertMsgCritical(y_coords.size() == labels.size(), errmol);

    errors::assertMsgCritical(z_coords.size() == labels.size(), errmol);

    std::vector<double> coords;

    coords.insert(coords.end(), x_coords.begin(), x_coords.end());

    coords.insert(coords.end(), y_coords.begin(), y_coords.end());

    coords.insert(coords.end(), z_coords.begin(), z_coords.end());

    // form charge, mass, label and elemental ID vectors

    std::string errelm("Molecule.from_xyz: Unsupported chemical element");

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

    return std::shared_ptr<CMolecule>(
            new CMolecule(coords, charges, masses, labels, idselem)
            );
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

static py::array
CMolecule_x_to_numpy(const CMolecule& self)
{
    py::list rx;

    for (int32_t i = 0; i < self.getNumberOfAtoms(); i++)
    {
        rx.append(self.getCoordinatesX()[i]);
    }

    return py::array(rx);
}

static py::array
CMolecule_y_to_numpy(const CMolecule& self)
{
    py::list ry;

    for (int32_t i = 0; i < self.getNumberOfAtoms(); i++)
    {
        ry.append(self.getCoordinatesY()[i]);
    }

    return py::array(ry);
}

static py::array
CMolecule_z_to_numpy(const CMolecule& self)
{
    py::list rz;

    for (int32_t i = 0; i < self.getNumberOfAtoms(); i++)
    {
        rz.append(self.getCoordinatesZ()[i]);
    }

    return py::array(rz);
}

// Helper function for getting VDW radii for molecule

static py::array
CMolecule_vdw_radii_to_numpy(const CMolecule& self)
{
    auto natoms = self.getNumberOfAtoms();

    auto atomradii = vdwradii::getRadii(self);

    py::list radii;

    for (int32_t i = 0; i < natoms; i++)
    {
        radii.append(atomradii[i]);
    }

    return py::array(radii);
}

// Helper function for getting nuclear charges for molecule

static py::array
CMolecule_elem_ids_to_numpy(const CMolecule& self)
{
    auto natoms = self.getNumberOfAtoms();

    auto idselem = self.getIdsElemental();

    py::list ids;

    for (int32_t i = 0; i < natoms; i++)
    {
        ids.append(idselem[i]);
    }

    return py::array(ids);
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

    std::string errmult("Molecule.check_multiplicity: ");

    errmult += "Incompatble multiplicity & number of electrons";

    errors::assertMsgCritical(flag, errmult);
}

// Helper function for checking proximity of atoms

static void
CMolecule_check_proximity(const CMolecule& self,
                          const double     minDistance)
{
    std::string errproxi("Molecule.check_proximity: Atoms too close");

    errors::assertMsgCritical(self.checkProximity(minDistance), errproxi);
}

// Helper function for broadcasting CMolecule object

static void
CMolecule_broadcast(CMolecule& self,
                    int32_t    rank,
                    py::object py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    self.broadcast(rank, *comm_ptr);
}

// Exports classes/functions in src/moldata to python

void export_moldata(py::module& m)
{
    // CMolecule class

    py::class_< CMolecule, std::shared_ptr<CMolecule> >
        (
            m, "Molecule"
        )
        .def(py::init<>())
        .def(py::init<const CMolecule&>())
        .def(py::init<const CMolecule&, const CMolecule&>())
        .def(py::init(&CMolecule_from_xyz))
        .def("set_charge", &CMolecule::setCharge)
        .def("get_charge", &CMolecule::getCharge)
        .def("set_multiplicity", &CMolecule::setMultiplicity)
        .def("get_multiplicity", &CMolecule::getMultiplicity)
        .def("check_multiplicity", &CMolecule_check_multiplicity)
        .def("get_string", &CMolecule::printGeometry)
        .def("check_proximity", &CMolecule_check_proximity)
        .def("get_sub_molecule", &CMolecule::getSubMolecule)
        .def("number_of_atoms",
             (int32_t (CMolecule::*)() const)
             &CMolecule::getNumberOfAtoms)
        .def("number_of_atoms",
             (int32_t (CMolecule::*)(const int32_t) const)
             &CMolecule::getNumberOfAtoms)
        .def("number_of_atoms",
             (int32_t (CMolecule::*)(const int32_t,
                                     const int32_t,
                                     const int32_t) const)
             &CMolecule::getNumberOfAtoms)
        .def("number_of_electrons", &CMolecule::getNumberOfElectrons)
        .def("number_of_alpha_electrons", &CMolecule_alpha_elec)
        .def("number_of_beta_electrons", &CMolecule_beta_elec)
        .def("nuclear_repulsion_energy", &CMolecule::getNuclearRepulsionEnergy)
        .def("x_to_numpy", &CMolecule_x_to_numpy)
        .def("y_to_numpy", &CMolecule_y_to_numpy)
        .def("z_to_numpy", &CMolecule_z_to_numpy)
        .def("vdw_radii_to_numpy", &CMolecule_vdw_radii_to_numpy)
        .def("elem_ids_to_numpy", &CMolecule_elem_ids_to_numpy)
        .def("get_elemental_composition", &CMolecule_get_elem_comp)
        .def("broadcast", &CMolecule_broadcast)
        .def(py::self == py::self)
    ;

    // CChemicalElement class

    py::class_< CChemicalElement, std::shared_ptr<CChemicalElement> >
        (
            m, "ChemicalElement"
        )
        .def(py::init<>())
        .def("set_atom_type",
             (bool (CChemicalElement::*)(const std::string&))
             &CChemicalElement::setAtomType)
        .def("set_atom_type",
             (bool (CChemicalElement::*)(const int32_t))
             &CChemicalElement::setAtomType)
        .def("get_name", &CChemicalElement::getName)
        .def(py::self == py::self)
    ;
}

} // vlx_moldata namespace
