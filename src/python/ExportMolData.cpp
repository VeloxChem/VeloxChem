//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>
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

namespace bp = boost::python;

namespace bp_moldata { // bp_moldata namespace

// Helper function for CMolecule constructor

static std::shared_ptr<CMolecule>
CMolecule_from_xyz(const bp::list& label_list,
                   const bp::list& x_list,
                   const bp::list& y_list,
                   const bp::list& z_list)
{
    // form coordinate vector

    std::string errmol("Molecule.from_xyz - Inconsistent lengths of lists");

    const int32_t natoms = (int32_t)bp::len(label_list);

    errors::assertMsgCritical(bp::len(x_list) == natoms, errmol);

    errors::assertMsgCritical(bp::len(y_list) == natoms, errmol);

    errors::assertMsgCritical(bp::len(z_list) == natoms, errmol);

    std::vector<double> coords;

    for (int32_t i = 0; i < natoms; i++)
    {
        coords.push_back(bp::extract<double>(x_list[i]));
    }

    for (int32_t i = 0; i < natoms; i++)
    {
        coords.push_back(bp::extract<double>(y_list[i]));
    }

    for (int32_t i = 0; i < natoms; i++)
    {
        coords.push_back(bp::extract<double>(z_list[i]));
    }

    // form charge, mass, label and elemental ID vectors

    std::string errelm("Molecule.from_xyz - Unsupported chemical element");

    std::vector<double> charges;

    std::vector<double> masses;

    std::vector<std::string> labels;

    std::vector<int32_t> idselem;

    for (int32_t i = 0; i < natoms; i++)
    {
        std::string lbl = bp::extract<std::string>(label_list[i]);

        CChemicalElement chemelm;

        auto err = chemelm.setAtomType(fstr::upcase(lbl));

        errors::assertMsgCritical(err, errelm);

        charges.push_back(chemelm.getAtomicCharge());

        masses.push_back(chemelm.getAtomicMass());

        labels.push_back(lbl);

        idselem.push_back(chemelm.getIdentifier());
    }

    // form molecule

    return std::shared_ptr<CMolecule>(
            new CMolecule(coords, charges, masses, labels, idselem)
            );
}

// Helper function for getting coordinates as numpy array

static np::ndarray
CMolecule_coordinates_to_numpy(const CMolecule& self)
{
    bp::list coords, rx, ry, rz;

    for (int32_t i = 0; i < self.getNumberOfAtoms(); i++)
    {
        rx.append(self.getCoordinatesX()[i]);

        ry.append(self.getCoordinatesY()[i]);

        rz.append(self.getCoordinatesZ()[i]);
    }

    coords.append(rx);

    coords.append(ry);

    coords.append(rz);

    return np::array(coords);
}

// Helper function for getting VDW radii for molecule

static np::ndarray
CMolecule_vdw_radii_to_numpy(const CMolecule& self)
{
    auto natoms = self.getNumberOfAtoms();

    auto atomradii = vdwradii::getRadii(self);

    bp::list radii;

    for (int32_t i = 0; i < natoms; i++)
    {
        radii.append(atomradii[i]);
    }

    return np::array(radii);
}

// Helper function for getting nuclear charges for molecule

static bp::list
CMolecule_get_ids_elem(const CMolecule& self)
{
    auto natoms = self.getNumberOfAtoms();

    auto idselem = self.getIdsElemental();

    bp::list ids;

    for (int32_t i = 0; i < natoms; i++)
    {
        ids.append(idselem[i]);
    }

    return ids;
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

// Helper function for broadcasting CMolecule object

static void
CMolecule_broadcast(CMolecule& self,
                    int32_t    rank,
                    bp::object py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    self.broadcast(rank, *comm_ptr);
}

// Exports classes/functions in src/moldata to python

void export_moldata()
{
    // CMolecule class
    // Note: Need member function pointers for proper overloading

    int32_t (CMolecule::*number_of_atoms_1)(
            ) const
        = &CMolecule::getNumberOfAtoms;

    int32_t (CMolecule::*number_of_atoms_2)(
            const int32_t idElemental) const
        = &CMolecule::getNumberOfAtoms;

    int32_t (CMolecule::*number_of_atoms_3)(
            const int32_t iAtom,
            const int32_t nAtoms,
            const int32_t idElemental) const
        = &CMolecule::getNumberOfAtoms;

    bp::class_< CMolecule, std::shared_ptr<CMolecule> >
        (
            "Molecule",
            bp::init<>()
        )
        .def(bp::init<const CMolecule&>())
        .def(bp::init<const CMolecule&, const CMolecule&>())
        .def("from_xyz", &CMolecule_from_xyz)
        .staticmethod("from_xyz")
        .def("set_charge", &CMolecule::setCharge)
        .def("get_charge", &CMolecule::getCharge)
        .def("set_multiplicity", &CMolecule::setMultiplicity)
        .def("get_multiplicity", &CMolecule::getMultiplicity)
        .def("check_multiplicity", &CMolecule_check_multiplicity)
        .def("print_geometry", &CMolecule::printGeometry)
        .def("check_proximity", &CMolecule::checkProximity)
        .def("get_sub_molecule", &CMolecule::getSubMolecule)
        .def("number_of_atoms", number_of_atoms_1)
        .def("number_of_atoms", number_of_atoms_2)
        .def("number_of_atoms", number_of_atoms_3)
        .def("number_of_electrons", &CMolecule::getNumberOfElectrons)
        .def("nuclear_repulsion_energy", &CMolecule::getNuclearRepulsionEnergy)
        .def("coordinates_to_numpy", &CMolecule_coordinates_to_numpy)
        .def("vdw_radii_to_numpy", &CMolecule_vdw_radii_to_numpy)
        .def("get_ids_elem", &CMolecule_get_ids_elem)
        .def("broadcast", &CMolecule_broadcast)
        .def(bp::self == bp::other<CMolecule>())
    ;
}

} // bp_moldata namespace
