//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <mpi.h>
#include <memory>
#include <vector>
#include <string>

#include "Codata.hpp"
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
CMolecule_from_coords(const std::vector<std::string>& labels,
                      const std::vector<double>&      coords_raw,
                      const std::string&              units)
{
    // NOTE:
    // The C++ Molecule constructor expects the coordinates to be arranged as 3 x natoms,
    // namely {x1, x2, x3, x4, ..., y1, y2, y3, y4, ..., z1, z2, z3, z4, ...}

    // sanity check

    std::string errmol("CMolecule_from_coords: Inconsistent lengths of lists");

    errors::assertMsgCritical(coords_raw.size() == labels.size() * 3, errmol);

    // scaling factor

    auto scale = 1.0 / units::getBohrValueInAngstroms();

    if (fstr::upcase(units) == "AU" || fstr::upcase(units) == "BOHR" ||
        fstr::upcase(units) == "BOHRS")
    {
        scale = 1.0;
    }

    std::vector<double> coords_au(coords_raw.size());

    for (size_t i = 0; i < coords_au.size(); i++)
    {
        coords_au[i] = coords_raw[i] * scale;
    }

    // form charge, mass, label and elemental ID vectors

    std::string errelm("CMolecule_from_coords: Unsupported chemical element");

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
            new CMolecule(coords_au, charges, masses, labels, idselem)
            );
}

static std::shared_ptr<CMolecule>
CMolecule_from_array(const std::vector<std::string>& labels,
                     const py::array_t<double>&      py_coords,
                     const std::string&              units=std::string("angs"))
{
    // NOTE:
    // The Python Molecule constructor expects the coordinates as a 2d numpy array,
    // namely np.array([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3], [x4, y4, z4], ...])

    // sanity check

    std::string errsrc("CMolecule_from_array: need a contiguous numpy array");

    auto c_style = py::detail::check_flags(py_coords.ptr(), py::array::c_style);

    auto f_style = py::detail::check_flags(py_coords.ptr(), py::array::f_style);

    errors::assertMsgCritical(c_style ^ f_style, errsrc);

    std::string errmol("CMolecule_from_array: Inconsistent size");

    errors::assertMsgCritical(
        py_coords.shape(0) == static_cast<ssize_t>(labels.size()), errmol);

    errors::assertMsgCritical(py_coords.shape(1) == 3, errmol);

    // form coordinate vector

    std::vector<double> coords(py_coords.size());

    if (c_style)
    {
        for (size_t d = 0; d < 3; d++)
        {
            for (size_t a = 0; a < labels.size(); a++)
            {
                // need to transpose py_coords for the C++ Molecule contructor
                coords[d * labels.size() + a] = py_coords.data()[a * 3 + d];
            }
        }
    }
    else if (f_style)
    {
        // no need to transpose py_coords for fortran style numpy array
        std::memcpy(coords.data(), py_coords.data(), py_coords.size() * sizeof(double));
    }

    return CMolecule_from_coords(labels, coords, units);
}

static std::shared_ptr<CMolecule>
CMolecule_from_array_2(const std::vector<int32_t>& idselem,
                       const py::array_t<double>&  py_coords,
                       const std::string&          units=std::string("angs"))
{
    std::vector<std::string> labels;

    std::string errelm("CMolecule_from_array: Unsupported element id");

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
    return vlx_general::pointer_to_numpy(self.getCoordinatesX(),
                                         self.getNumberOfAtoms());
}

static py::array_t<double>
CMolecule_y_to_numpy(const CMolecule& self)
{
    return vlx_general::pointer_to_numpy(self.getCoordinatesY(),
                                         self.getNumberOfAtoms());
}

static py::array_t<double>
CMolecule_z_to_numpy(const CMolecule& self)
{
    return vlx_general::pointer_to_numpy(self.getCoordinatesZ(),
                                         self.getNumberOfAtoms());
}

// Helper function for getting VDW radii for molecule

static py::array_t<double>
CMolecule_vdw_radii_to_numpy(const CMolecule& self)
{
    auto atomradii = vdwradii::getRadii(self);

    return vlx_general::pointer_to_numpy(atomradii.data(),
                                         self.getNumberOfAtoms());
}

// Helper function for getting nuclear charges for molecule

static py::array
CMolecule_elem_ids_to_numpy(const CMolecule& self)
{
    py::list ids;

    for (int32_t i = 0; i < self.getNumberOfAtoms(); i++)
    {
        ids.append(self.getIdsElemental()[i]);
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
        .def(py::init(&CMolecule_from_array),
             py::arg(), py::arg(), py::arg("units")=std::string("angs"))
        .def(py::init(&CMolecule_from_array_2),
             py::arg(), py::arg(), py::arg("units")=std::string("angs"))
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
