//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <mpi.h>
#include <mpi4py/mpi4py.h>
#include <string>

#include "Codata.hpp"
#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "MpiFunc.hpp"
#include "SpinBlock.hpp"
#include "StringFormat.hpp"

namespace py = pybind11;

namespace vlx_general {  // vlx_general namespace

// Gets MPI_Comm pointer from a mpi4py communicator object
// Not a static function; used in other files

MPI_Comm*
get_mpi_comm(py::object py_comm)
{
    PyObject* py_obj = py_comm.ptr();

    MPI_Comm* comm_ptr = PyMPIComm_Get(py_obj);

    if (comm_ptr == NULL) throw py::error_already_set();

    return comm_ptr;
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
            strd *= dimension[j];
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
        return py::array_t<int32_t>(
            dimension_to_shape(dimension), dimension_to_strides(dimension, sizeof(int32_t)), ptr);
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
    // initialize mpi4py's C-API

    auto err = import_mpi4py();

    std::string errmpi4py("mpi4py: failed to import mpi4py");

    errors::assertMsgCritical(err == 0, errmpi4py);

    // szblock enum class

    py::enum_<szblock>(m, "szblock")
        .value("aa", szblock::aa)
        .value("ab", szblock::ab)
        .value("ba", szblock::ba)
        .value("bb", szblock::bb);

    // exposing functions

    m.def("mpi_master", &mpi::master);

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
