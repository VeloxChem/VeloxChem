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

#include "ExportGeneral.hpp"

#include <mpi.h>
// see here: https://github.com/mpi4py/mpi4py/issues/19#issuecomment-768143143
#ifdef MSMPI_VER
#define PyMPI_HAVE_MPI_Message 1
#endif
#include <mpi4py/mpi4py.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <string>

#ifdef ENABLE_MKL
#include "ConfigMKL.hpp"
#endif
#include "Codata.hpp"
#include "ErrorHandler.hpp"
#include "MpiFunc.hpp"
#include "SpinBlock.hpp"
#include "StringFormat.hpp"
#include "TwoIndexes.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_general {  // vlx_general namespace

// Gets MPI_Comm pointer from a mpi4py communicator object
// Not a static function; used in other files

MPI_Comm*
get_mpi_comm(py::object py_comm)
{
    auto comm_ptr = PyMPIComm_Get(py_comm.ptr());

    if (!comm_ptr) throw py::error_already_set();

    return comm_ptr;
}

// Exports classes/functions in src/general to python

void
export_general(py::module& m)
{
    // configure MKL single dynamic library

#ifdef ENABLE_MKL
    configure_mkl_rt();
#endif

    // initialize mpi4py's C-API

    if (import_mpi4py() < 0)
    {
        // mpi4py calls the Python C API
        // we let pybind11 give us the detailed traceback
        throw py::error_already_set();
    }

    // szblock enum class

    // clang-format off
    py::enum_<szblock>(m, "szblock")
        .value("aa", szblock::aa)
        .value("ab", szblock::ab)
        .value("ba", szblock::ba)
        .value("bb", szblock::bb);
    // clang-format on

    // CTwoIndexes class

    PyClass<CTwoIndexes>(m, "TwoIndexes")
        .def(py::init<>())
        .def(py::init<const int32_t, const int32_t>())
        .def("first", &CTwoIndexes::first, "Gets first index from pair of indexes.")
        .def("second", &CTwoIndexes::second, "Gets second index from pair of indexes.");

    // exposing functions

    m.def("mpi_master", &mpi::master, "Gets default rank of master MPI process.");
    m.def("mpi_initialized", &mpi::initialized, "Check if MPI has been initialized.");
    m.def(
        "mpi_size_limit",
        []() -> int32_t { return static_cast<int32_t>(1 << 30) / 5 * 9; },
        "Gets the size limit in MPI communication (below 2^31-1).");
    m.def(
        "is_mpi_master",
        [](py::object py_comm) -> bool {
            if (py_comm.is_none())
                return (mpi::rank(MPI_COMM_WORLD) == mpi::master());
            else
                return (mpi::rank(*get_mpi_comm(py_comm)) == mpi::master());
        },
        "Checks if a MPI process is the master process.",
        "py_comm"_a = py::none());
    m.def(
        "mpi_barrier",
        [](py::object py_comm) {
            if (py_comm.is_none())
                MPI_Barrier(MPI_COMM_WORLD);
            else
                MPI_Barrier(*get_mpi_comm(py_comm));
        },
        "Synchronize all MPI processes using barrier.",
        "py_comm"_a = py::none());
    m.def(
        "is_single_node",
        [](py::object py_comm) -> bool {
            if (py_comm.is_none())
                return (mpi::nodes(MPI_COMM_WORLD) == 1);
            else
                return (mpi::nodes(*get_mpi_comm(py_comm)) == 1);
        },
        "Checks if there is only one MPI process in the communicator.",
        "py_comm"_a = py::none());

    m.def("assert_msg_critical", &errors::assertMsgCritical, "Prints message and aborts in case of a critical error.", "condition"_a, "message"_a);

    m.def("bohr_in_angstroms", &units::getBohrValueInAngstroms, "Gets Bohr value in Angstroms.");
    m.def("hartree_in_ev", &units::getHartreeValueInElectronVolts, "Gets Hartree value in electronvolts.");
    m.def("hartree_in_kcalpermol", &units::getHartreeValueInKiloCaloriePerMole, "Gets Hartree value in kcal/mol.");
    m.def("hartree_in_inverse_nm", &units::getHartreeValueInInverseNanometer, "Gets Hartree value in inverse nanometer.");
    m.def("hartree_in_wavenumbers", &units::getHartreeValueInWavenumbers, "Gets Hartree value in reciprocal cm.");
    m.def("electron_mass_in_amu", &units::getElectronMassInAtomicMassUnit, "Gets electron mass in amu.");
    m.def("amu_in_electron_masses", &units::getAtomicMassUnitInElectronMasses, "Gets atomic mass unit in electron masses.");
    m.def("amu_in_kg", &units::getAtomicMassUnitInKg, "Gets atomic mass unit in kg.");
    m.def("speed_of_light_in_vacuum_in_SI", &units::getSpeedOfLightInVacuumInSI, "Gets speed of light in vacuum in SI.");
    m.def("avogadro_constant", &units::getAvogadroConstant, "Gets Avogadro constant.");
    m.def("boltzmann_in_evperkelvin", &units::getBoltzmannConstantInElectronVoltsPerKelvin, "Gets Boltzmann constant in eV/K.");
    m.def("boltzmann_in_hartreeperkelvin", &units::getBoltzmannConstantInHartreePerKelvin, "Gets Boltzmann constant in Hartree/K.");

    m.def("dipole_in_debye", &units::getDipoleInDebye, "Gets convertion factor for dipole moment (a.u. -> Debye).");
    m.def("rotatory_strength_in_cgs", &units::getRotatoryStrengthInCGS, "Gets convertion factor for rotatory strength (a.u. -> 10^-40 cgs).");
    m.def("extinction_coefficient_from_beta",
          &units::getExtinctionCoefficientFromBeta,
          "Gets factor needed for the calculation of the extinction coefficent from the electric-dipole magnetic-dipole polarizability beta.");
    m.def("fine_structure_constant", &units::getFineStructureConstant, "Gets fine-structure constant.");

    m.def(
        "to_angular_momentum",
        [](const int32_t angl) -> std::string { return fstr::to_AngularMomentum(angl); },
        "Converts angular momentum integer to string.",
        "angl"_a);
    m.def(
        "to_angular_momentum",
        [](const std::string& label) -> int32_t { return fstr::to_AngularMomentum(label); },
        "Converts angular momentum string to integer.",
        "label"_a);
}

}  // namespace vlx_general
