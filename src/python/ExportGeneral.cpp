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
#include "MpiFunc.hpp"
#include "OpenMPFunc.hpp"
#include "Point.hpp"
#include "SphericalMomentum.hpp"
#include "StringFormat.hpp"
#include "TensorComponents.hpp"
#include "TensorLabels.hpp"

namespace vlx_general {

// Gets MPI_Comm pointer from a mpi4py communicator object
// Not a static function; used in other files

auto
get_mpi_comm(py::object py_comm) -> MPI_Comm*
{
    auto comm_ptr = PyMPIComm_Get(py_comm.ptr());

    if (!comm_ptr) throw py::error_already_set();

    return comm_ptr;
}

// Gets shape and strides from dimension

static auto
get_shape_and_strides(const std::vector<int>& dimension) -> std::tuple<std::vector<py::ssize_t>, std::vector<py::ssize_t>>
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
pointer_to_numpy(const double* ptr, const std::vector<int>& dimension) -> py::array_t<double>
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

auto
export_general(py::module &m) -> void
{
    /// exposing functions from StringFormat.hpp
    m.def("upper_case", &format::upper_case, "Creates upper cased copy of given string.");
    m.def("lower_case", &format::lower_case, "Creates lower cased copy of given string.");

    // exposing functions from MpiFunc.hpp
    m.def("mpi_master", &mpi::master, "Returns rank of MPI master process.");
    m.def("mpi_initialized", &mpi::initialized, "Check if MPI has been initialized.");
    m.def(
        "mpi_size_limit",
        []() -> int { return static_cast<int>(1 << 30) / 5 * 9; },
        "Gets the size limit in MPI communication (below 2^31-1).");

    // exposing functions from Codata.hpp
    m.def("bohr_in_angstrom", &units::bohr_in_angstrom, "Gets Bohr value in Angstroms.");
    m.def("hartree_in_ev", &units::hartree_in_ev, "Gets Hartree value in electronvolts.");
    m.def("hartree_in_kcalpermol", &units::getHartreeValueInKiloCaloriePerMole, "Gets Hartree value in kcal/mol.");
    m.def("hartree_in_inverse_nm", &units::getHartreeValueInInverseNanometer, "Gets Hartree value in inverse nanometer.");
    m.def("hartree_in_wavenumber", &units::getHartreeValueInWavenumbers, "Gets Hartree value in reciprocal cm.");
    m.def("hartree_in_wavenumbers", &units::getHartreeValueInWavenumbers, "Gets Hartree value in reciprocal cm.");
    m.def("electron_mass_in_amu", &units::getElectronMassInAtomicMassUnit, "Gets electron mass in amu.");
    m.def("amu_in_electron_mass", &units::getAtomicMassUnitInElectronMasses, "Gets atomic mass unit in electron masses.");
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

    // exposing functions from TensorLabels.hpp
    m.def("tensor_cartesian_labels", &tensor::cartesian_labels, "Gets all Cartesian component labels of tensor.");
    m.def("tensor_spherical_labels", &tensor::spherical_labels, "Gets all spherical component labels of tensor.");
    m.def("tensor_cartesian_index", &tensor::cartesian_index, "Gets index of Cartesian tensor component.");
    m.def("tensor_label", &tensor::label, "Gets label of tensor.");
    m.def("tensor_order", &tensor::order, "Gets order of tensor.");

    // exposing functions from TensorComponents.hpp
    m.def(
        "number_of_cartesian_components",
        [](const int order) -> int { return tensor::number_of_cartesian_components(std::array<int, 1>{order}); },
        "Gets number of Cartesian components in tensor.");
    m.def(
        "number_of_cartesian_components",
        [](const std::array<int, 2> orders) -> int { return tensor::number_of_cartesian_components(orders); },
        "Gets number of Cartesian components in array of tensors.");
    m.def(
        "number_of_cartesian_components",
        [](const std::array<int, 3> orders) -> int { return tensor::number_of_cartesian_components(orders); },
        "Gets number of Cartesian components in array of tensors.");
    m.def(
        "number_of_cartesian_components",
        [](const std::array<int, 4> orders) -> int { return tensor::number_of_cartesian_components(orders); },
        "Gets number of Cartesian components in array of tensors.");
    m.def(
        "number_of_spherical_components",
        [](const int order) -> int { return tensor::number_of_spherical_components(std::array<int, 1>{order}); },
        "Gets number of spherical components in tensor.");
    m.def(
        "number_of_spherical_components",
        [](const std::array<int, 2> orders) -> int { return tensor::number_of_spherical_components(orders); },
        "Gets number of spherical components in array of tensors.");
    m.def(
        "number_of_spherical_components",
        [](const std::array<int, 3> orders) -> int { return tensor::number_of_spherical_components(orders); },
        "Gets number of spherical components in array of tensors.");
    m.def(
        "number_of_spherical_components",
        [](const std::array<int, 4> orders) -> int { return tensor::number_of_spherical_components(orders); },
        "Gets number of spherical components in array of tensors.");

    // exposing functions from BatchFunc.hpp
    m.def("number_of_batches", &batch::number_of_batches<size_t>, "Gets number of batches.");
    m.def("batch_range",
          py::overload_cast<const size_t, const size_t, const size_t>(&batch::batch_range<size_t>),
          "Gets [first, last) range for requested batch.");
    m.def("batch_range",
          py::overload_cast<const size_t, const size_t, const size_t, const size_t>(&batch::batch_range<size_t>),
          "Gets [first, last) range for requested batch.");

    // exposing functions from OpenMPFunc.hpp
    m.def("set_number_of_threads", &omp::set_number_of_threads, "Sets number of OMP threads to requested value.");
    m.def("get_number_of_threads", &omp::get_number_of_threads, "Gets number of OMP threads available.");
    m.def("make_work_tasks",
          py::overload_cast<const std::vector<CGtoBlock> &>(&omp::make_work_tasks),
          "Gets work tasks for given vector of basis function blocks.");
    m.def("make_work_tasks",
          py::overload_cast<const std::vector<CGtoBlock> &, const std::vector<CGtoBlock> &>(&omp::make_work_tasks),
          "Gets work tasks for given two vectors of basis function blocks.");
    m.def("make_diag_work_tasks", &omp::make_diag_work_group);
    m.def("make_work_group", &omp::make_work_group);

    // exposing functions from SphericalMomentum.hpp
    m.def("spherical_momentum_s_factors", spher_mom::transformation_factors<0>, "Gets transformation factors for S type spherical momentum.");
    m.def("spherical_momentum_p_factors", spher_mom::transformation_factors<1>, "Gets transformation factors for P type spherical momentum.");
    m.def("spherical_momentum_d_factors", spher_mom::transformation_factors<2>, "Gets transformation factors for D type spherical momentum.");
    m.def("spherical_momentum_f_factors", spher_mom::transformation_factors<3>, "Gets transformation factors for F type spherical momentum.");
    m.def("spherical_momentum_g_factors", spher_mom::transformation_factors<4>, "Gets transformation factors for G type spherical momentum.");

    // TPoint class
    PyClass<TPoint<double>>(m, "Point")
        .def(py::init<>())
        .def(py::init<const std::array<double, 3> &>())
        .def(py::init<const TPoint<double> &>())
        .def(py::pickle([](const TPoint<double> &pnt) { return py::make_tuple(pnt.coordinates()); },
                        [](py::tuple t) { return TPoint<double>(t[0].cast<std::array<double, 3>>()); }))
        .def("coordinates", &TPoint<double>::coordinates, "Getter for Cartesian coordinates.")
        .def("scale", &TPoint<double>::scale, "Scales Cartesian coordinates by factor.")
        .def("length_square", &TPoint<double>::length_square, "Computes square of length for vector given by point.")
        .def("length", &TPoint<double>::length, "Computes length for vector given by point.")
        .def("distance_square", &TPoint<double>::distance_square, "Computes square of distance between two points.")
        .def("distance", &TPoint<double>::distance, "Computes distance between two points.")
        .def("__eq__", [](const TPoint<double> &self, const TPoint<double> &other) { return self == other; })
        .def("__ne__", [](const TPoint<double> &self, const TPoint<double> &other) { return self != other; })
        .def("__copy__", [](const TPoint<double> &self) { return TPoint<double>(self); })
        .def("__deepcopy__", [](const TPoint<double> &self, py::dict) { return TPoint<double>(self); });
}

}  // namespace vlx_general
