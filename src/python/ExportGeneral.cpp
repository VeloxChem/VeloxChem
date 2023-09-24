#include "ExportGeneral.hpp"

#include <mpi.h>
// see here: https://github.com/mpi4py/mpi4py/issues/19#issuecomment-768143143
#ifdef MSMPI_VER
#define PyMPI_HAVE_MPI_Message 1
#endif
#include <mpi4py/mpi4py.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "BatchFunc.hpp"
#include "Codata.hpp"
#include "GtoBlock.hpp"
#include "MpiFunc.hpp"
#include "OpenMPFunc.hpp"
#include "StringFormat.hpp"

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

// Gets shape and strides from dimension

static auto
get_shape_and_strides(const std::vector<int64_t>& dimension) -> std::tuple<std::vector<py::ssize_t>, std::vector<py::ssize_t>>
{
    std::vector<py::ssize_t> shape, strides;

    for (size_t i = 0; i < dimension.size(); i++)
    {
        shape.push_back(static_cast<py::ssize_t>(dimension[i]));

        size_t strd = 1;

        for (size_t j = i + 1; j < dimension.size(); j++)
        {
            strd *= dimension[j];
        }

        strides.push_back(static_cast<py::ssize_t>(strd * sizeof(double)));
    }

    return {shape, strides};
}

// Creates numpy ndarray from pointer and dimension
// Not a static function; used in other files

auto
pointer_to_numpy(const double* ptr, const std::vector<int64_t>& dimension) -> py::array_t<double>
{
    if ((ptr == nullptr) || (dimension.size() == 0))
    {
        return py::array_t<double>();
    }
    else
    {
        const auto [shape, strides] = get_shape_and_strides(dimension);

        return py::array_t<double>(shape, strides, ptr);
    }
}

// Exports classes/functions in src/general to python

auto
export_general(py::module& m) -> void
{
    // initialize mpi4py's C-API

    if (import_mpi4py() < 0)
    {
        // mpi4py calls the Python C API
        // we let pybind11 give us the detailed traceback
        throw py::error_already_set();
    }

    // exposing enum from FmtType.hpp

    py::enum_<fmt_t>(m, "fmt_t").value("center", fmt_t::center).value("left", fmt_t::left).value("right", fmt_t::right);

    // exposing functions from Codata.hpp

    m.def("bohr_in_angstroms", &units::getBohrValueInAngstroms, "Gets Bohr value in Angstroms.");

    m.def("hartree_in_ev", &units::getHartreeValueInElectronVolts, "Gets Hartree value in electronvolts.");

    // exposing functions from BatchFunc.hpp

    m.def("get_batch_index", &batch::getBatchIndex, "Gets starting index of batch.");

    m.def("number_of_batches", &batch::getNumberOfBatches, "Gets number of batches for given vector size.");

    m.def("get_batch_range", &batch::getBatchRange, "Gets range of specific batch.");

    // exposing functions from OpenMPFunc.hpp

    m.def("set_number_of_threads", &omp::setNumberOfThreads, "Sets number of OMP threads to requested value.");

    m.def("get_number_of_threads", &omp::getNumberOfThreads, "Gets number of OMP threads available.");

    m.def("make_workgroup",
          py::overload_cast<const std::vector<CGtoBlock>&>(&omp::makeWorkGroup),
          "Gets work group for given vector of basis function blocks.");

    m.def("make_workgroup",
          py::overload_cast<const std::vector<CGtoBlock>&, const std::vector<CGtoBlock>&>(&omp::makeWorkGroup),
          "Gets work group for given two vectors of basis function blocks.");

    // exposing functions from StringFormat.hpp

    m.def("upcase", &fstr::upcase, "Convers string to upper case string.");

    m.def("format", &fstr::format, "Formats string to string with specified alignment.");

    m.def(
        "to_string",
        [](const double source, const size_t presicion, const size_t width, const fmt_t aligment) -> std::string {
            return fstr::to_string(source, presicion, width, aligment);
        },
        "Formats double precision number to string with specified alignment.");

    m.def(
        "to_string",
        [](const double source, const size_t presicion) -> std::string { return fstr::to_string(source, presicion); },
        "Formats double precision number to string with specified alignment.");

    m.def(
        "to_string",
        [](const int64_t source, const size_t width, const fmt_t aligment) -> std::string { return fstr::to_string(source, width, aligment); },
        "Formats integer number to string with specified alignment.");

    m.def(
        "to_string", [](const bool source) -> std::string { return fstr::to_string(source); }, "Formats bool to string.");

    m.def(
        "to_angular_momentum",
        [](const int64_t angmom) -> std::string { return fstr::to_AngularMomentum(angmom); },
        "Converts angular momentum integer to string.");

    m.def(
        "to_angular_momentum",
        [](const std::string& label) -> int64_t { return fstr::to_AngularMomentum(label); },
        "Converts angular momentum string to integer.");

    m.def("mpi_master", &mpi::master, "Returns rank of MPI master process.");

    m.def(
        "bcast_scalar",
        [](const int64_t val, py::object py_comm) -> int64_t { return mpi::bcastScalar(val, *get_mpi_comm(py_comm)); },
        "Broadcasts scalar.");

    m.def(
        "bcast_scalar",
        [](const double val, py::object py_comm) -> double { return mpi::bcastScalar(val, *get_mpi_comm(py_comm)); },
        "Broadcasts scalar.");

    m.def(
        "bcast_dense_matrix",
        [](const CDenseMatrix& matrix, py::object py_comm) -> CDenseMatrix { return mpi::bcastDenseMatrix(matrix, *get_mpi_comm(py_comm)); },
        "Broadcasts dense matrix.");

    m.def(
        "scatter_vector",
        [](const std::vector<int64_t>& vec, py::object py_comm) -> std::vector<int64_t> {
            return mpi::scatterStdVector(vec, *get_mpi_comm(py_comm));
        },
        "Scatters vector.");

    m.def(
        "gather_dense_matrices_by_columns",
        [](const CDenseMatrix& matrix, py::object py_comm) -> CDenseMatrix {
            return mpi::gatherDenseMatricesByColumns(matrix, *get_mpi_comm(py_comm));
        },
        "Gathers dense matrices by columns.");
}

}  // namespace vlx_general
