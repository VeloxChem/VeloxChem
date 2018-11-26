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
#include "ExportGeneral.hpp"
#include "ExportMolData.hpp"

namespace bp = boost::python;

namespace bp_moldata { // bp_moldata namespace

// Helper function for CMolecule constructor

static std::shared_ptr<CMolecule>
CMolecule_from_list(const bp::list& coord_list,
                    const bp::list& charge_list,
                    const bp::list& mass_list,
                    const bp::list& label_list,
                    const bp::list& idelem_list)
{
    std::string errmol("Molecule.from_numpy_list: Inconsistent lengths of lists!");

    const int natoms = bp::len(charge_list);
    errors::assertMsgCritical(bp::len(coord_list)  == natoms * 3, errmol);
    errors::assertMsgCritical(bp::len(mass_list)   == natoms, errmol);
    errors::assertMsgCritical(bp::len(label_list)  == natoms, errmol);
    errors::assertMsgCritical(bp::len(idelem_list) == natoms, errmol);

    std::vector<double> coords;
    std::vector<double> charges;
    std::vector<double> masses;
    std::vector<std::string> labels;
    std::vector<int32_t> idselem;

    for (int i = 0; i < natoms * 3; i++)
    {
        double x = bp::extract<double>(coord_list[i]);
        coords.push_back(x);
    }

    for (int i = 0; i < natoms; i++)
    {
        double      chg  = bp::extract<double>(charge_list[i]);
        double      mass = bp::extract<double>(mass_list[i]);
        std::string lbl  = bp::extract<std::string>(label_list[i]);
        int32_t     id   = bp::extract<int32_t>(idelem_list[i]);

        charges.push_back(chg);
        masses.push_back(mass);
        labels.push_back(lbl);
        idselem.push_back(id);
    }

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
CMolecule_vdw_radii_to_numpy(CMolecule& self)
{
    auto atomradii = vdwradii::getRadii(self);

    bp::list radii;

    for (size_t i = 0; i < atomradii.size(); i++)
    {
        radii.append(atomradii[i]);
    }

    return np::array(radii);
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
        .def("from_list", &CMolecule_from_list)
        .staticmethod("from_list")
        .def("set_charge", &CMolecule::setCharge)
        .def("get_charge", &CMolecule::getCharge)
        .def("set_multiplicity", &CMolecule::setMultiplicity)
        .def("get_multiplicity", &CMolecule::getMultiplicity)
        .def("print_geometry", &CMolecule::printGeometry)
        .def("get_sub_molecule", &CMolecule::getSubMolecule)
        .def("number_of_atoms", number_of_atoms_1)
        .def("number_of_atoms", number_of_atoms_2)
        .def("number_of_atoms", number_of_atoms_3)
        .def("number_of_electrons", &CMolecule::getNumberOfElectrons)
        .def("nuclear_repulsion_energy", &CMolecule::getNuclearRepulsionEnergy)
        .def("coordinates_to_numpy", &CMolecule_coordinates_to_numpy)
        .def("vdw_radii_to_numpy", &CMolecule_vdw_radii_to_numpy)
        .def("broadcast", &CMolecule_broadcast)
        .def(bp::self == bp::other<CMolecule>())
    ;
}

} // bp_moldata namespace
