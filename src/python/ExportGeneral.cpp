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

#include "Codata.hpp"
#include "ErrorHandler.hpp"
#include "MpiFunc.hpp"
#include "SpinBlock.hpp"
#include "StringFormat.hpp"
#ifdef ENABLE_MKL
#include "ConfigMKL.hpp"
#endif

namespace py = pybind11;

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

// Helper function for checking master node

static bool
is_mpi_master(py::object py_comm)
{
    if (py_comm.is_none())
    {
        return (mpi::rank(MPI_COMM_WORLD) == mpi::master());
    }
    else
    {
        auto comm = get_mpi_comm(py_comm);
        return (mpi::rank(*comm) == mpi::master());
    }
}

// Helper function for checking number of nodes

static bool
is_single_node(py::object py_comm)
{
    if (py_comm.is_none())
    {
        return (mpi::nodes(MPI_COMM_WORLD) == 1);
    }
    else
    {
        auto comm = get_mpi_comm(py_comm);
        return (mpi::nodes(*comm) == 1);
    }
}

// Helper function for checking number of nodes

// Helper function for getting the size limit in MPI communication

static int32_t
mpi_size_limit()
{
    return static_cast<int32_t>(1 << 30) / 5 * 9;
}

// Helper functions for getting shape and strides from int32_t dimension

static std::vector<ssize_t>
dimension_to_shape(const std::vector<int32_t>& dimension)
{
    std::vector<ssize_t> shape;

    for (size_t i = 0; i < dimension.size(); i++)
    {
        shape.push_back(static_cast<ssize_t>(dimension[i]));
    }

    return shape;
}

static std::vector<ssize_t>
dimension_to_strides(const std::vector<int32_t>& dimension, size_t sizeoftype)
{
    std::vector<ssize_t> strides;

    for (size_t i = 0; i < dimension.size(); i++)
    {
        size_t strd = 1;

        for (size_t j = i + 1; j < dimension.size(); j++)
        {
            strd *= static_cast<size_t>(dimension[j]);
        }

        strides.push_back(static_cast<ssize_t>(strd * sizeoftype));
    }

    return strides;
}

// Gets numpy array from double pointer and int32_t dimensions
// Not static functions; used in other files

py::array_t<double>
pointer_to_numpy(const double* ptr, const std::vector<int32_t>& dimension)
{
    if (ptr == nullptr || dimension.size() == 0)
    {
        return py::array_t<double>();
    }
    else
    {
        return py::array_t<double>(dimension_to_shape(dimension), dimension_to_strides(dimension, sizeof(double)), ptr);
    }
}

py::array_t<double>
pointer_to_numpy(const double* ptr, int32_t nElements)
{
    return pointer_to_numpy(ptr, std::vector<int32_t>({nElements}));
}

py::array_t<double>
pointer_to_numpy(const double* ptr, int32_t nRows, int32_t nColumns)
{
    return pointer_to_numpy(ptr, std::vector<int32_t>({nRows, nColumns}));
}

// Gets numpy array from int32_t pointer and dimensions
// Not static functions; used in other files

py::array_t<int32_t>
pointer_to_numpy(const int32_t* ptr, const std::vector<int32_t>& dimension)
{
    if (ptr == nullptr || dimension.size() == 0)
    {
        return py::array_t<int32_t>();
    }
    else
    {
        return py::array_t<int32_t>(dimension_to_shape(dimension), dimension_to_strides(dimension, sizeof(int32_t)), ptr);
    }
}

py::array_t<int32_t>
pointer_to_numpy(const int32_t* ptr, int32_t nElements)
{
    return pointer_to_numpy(ptr, std::vector<int32_t>({nElements}));
}

py::array_t<int32_t>
pointer_to_numpy(const int32_t* ptr, int32_t nRows, int32_t nColumns)
{
    return pointer_to_numpy(ptr, std::vector<int32_t>({nRows, nColumns}));
}

// Helper function for converting angular momentum

static std::string
string_to_angular_momentum(const int32_t angl)
{
    return fstr::to_AngularMomentum(angl);
}

static int32_t
integer_to_angular_momentum(const std::string& label)
{
    return fstr::to_AngularMomentum(label);
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

    // exposing functions

    m.def("mpi_master", &mpi::master);

    m.def("is_mpi_master", &is_mpi_master, py::arg("py_comm") = py::none());

    m.def("is_single_node", &is_single_node, py::arg("py_comm") = py::none());

    m.def("mpi_size_limit", &mpi_size_limit);

    m.def("mpi_initialized", &mpi::initialized);

    m.def("assert_msg_critical", &errors::assertMsgCritical);

    m.def("bohr_in_angstroms", &units::getBohrValueInAngstroms);

    m.def("hartree_in_ev", &units::getHartreeValueInElectronVolts);

    m.def("hartree_in_kcalpermol", &units::getHartreeValueInKiloCaloriePerMole);

    m.def("dipole_in_debye", &units::getDipoleInDebye);

    m.def("rotatory_strength_in_cgs", &units::getRotatoryStrengthInCGS);

    m.def("to_angular_momentum", &string_to_angular_momentum);

    m.def("to_angular_momentum", &integer_to_angular_momentum);
}

}  // namespace vlx_general
