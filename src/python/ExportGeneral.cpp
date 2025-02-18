#include "ExportGeneral.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "BatchFunc.hpp"
#include "Codata.hpp"
#include "OpenMPFunc.hpp"
#include "Point.hpp"
#include "SphericalMomentum.hpp"
#include "StringFormat.hpp"
#include "TensorComponents.hpp"
#include "TensorLabels.hpp"
#include "T3FlatBuffer.hpp"
#include "T3RectFlatBuffer.hpp"
#include "MathFunc.hpp"

namespace vlx_general {

auto
export_general(py::module &m) -> void
{
    /// exposing functions from StringFormat.hpp
    m.def("upper_case", &format::upper_case, "Creates upper cased copy of given string.");
    m.def("lower_case", &format::lower_case, "Creates lower cased copy of given string.");

    // exposing functions from Codata.hpp
    m.def("bohr_in_angstrom", &units::bohr_in_angstrom, "Gets Bohr value in Angstrom.");
    m.def("hartree_in_ev", &units::hartree_in_ev, "Gets Hartree value in electronvolt.");

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
    m.def("make_diag_work_tasks", &omp::make_diag_work_group, "Gets work tasks for diagonal integrals.");
    m.def("make_work_group",
          py::overload_cast<const std::vector<CBlockedGtoPairBlock>&, const int>(&omp::make_work_group),
          "Gets work group for ERIs.");
    m.def("make_work_group",
          py::overload_cast<const std::vector<CBlockedGtoPairBlock>&, const int, const int>(&omp::make_work_group),
          "Gets work group for ERIs.");
    m.def("partition_atoms", &omp::partition_atoms, "Get atomic indices of partitioned atoms list.");

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
    
    // CT3FlatBuffer class
    PyClass<CT3FlatBuffer<double>>(m, "T3FlatBuffer")
        .def(py::init<>())
        .def(py::init<const std::vector<size_t>&, const size_t>())
        .def(py::init<const CT3FlatBuffer<double> &>())
        .def("indices", &CT3FlatBuffer<double>::indices, "Gets indices vector along x axis of tensor.")
        .def(
            "values",
             [](const CT3FlatBuffer<double> &self, const size_t index) -> py::array_t<double> {
                 const auto ndim = self.width();
                 const auto nelems = static_cast<py::ssize_t>(ndim * (ndim + 1) / 2);
                 const auto tdim  = static_cast<py::ssize_t>(sizeof(double));
                 return py::array_t<double>(
                                            std::vector<py::ssize_t>({nelems,}), std::vector<py::ssize_t>({tdim, }), self.data(index));
            },
            "Gets slice of tensor values along y,z axes.")
        .def("value",
             [](const CT3FlatBuffer<double> &self, const size_t index, const size_t i, const size_t j) -> double {
                if (i <= j)
                {
                    return self.data(index)[mathfunc::uplo_rm_index(i, j, self.width())];
                }
                else
                {
                    return self.data(index)[mathfunc::uplo_rm_index(j, i, self.width())];
                }
             },
             "Gets tensor element value.")
        .def("__copy__", [](CT3FlatBuffer<double> &self) { return CT3FlatBuffer<double>(self); })
        .def("__deepcopy__", [](const CT3FlatBuffer<double> &self, py::dict) { return CT3FlatBuffer<double>(self); });
    
    // CT3RectFlatBuffer class
    PyClass<CT3RectFlatBuffer<double>>(m, "T3RectFlatBuffer")
        .def(py::init<>())
        .def(py::init<const std::vector<size_t>&, const std::vector<size_t>&, const std::vector<size_t>&>())
        .def(py::init<const CT3RectFlatBuffer<double> &>())
        .def("indices", &CT3RectFlatBuffer<double>::indices, "Gets indices vector along x axis of tensor.")
        .def("row_indices", &CT3RectFlatBuffer<double>::row_indices, "Gets indices vector along y axis of tensor.")
        .def("col_indices", &CT3RectFlatBuffer<double>::col_indices, "Gets indices vector along z axis of tensor.")
        .def("number_of_rows", &CT3RectFlatBuffer<double>::number_of_rows, "Gets number of rows along y axis of tensor.")
        .def("number_of_columns", &CT3RectFlatBuffer<double>::number_of_columns, "Gets number of columns along z axis of tensor.")
        .def(
            "values",
             [](const CT3RectFlatBuffer<double> &self, const size_t index) -> py::array_t<double> {
                 const auto nrows = self.number_of_rows();
                 const auto ncols = self.number_of_columns();
                 const auto nelems = static_cast<py::ssize_t>(nrows * ncols);
                 const auto tdim  = static_cast<py::ssize_t>(sizeof(double));
                 return py::array_t<double>(
                                            std::vector<py::ssize_t>({nelems,}), std::vector<py::ssize_t>({tdim, }), self.data(index));
            },
            "Gets slice of tensor values along y,z axes.")
        .def("value",
             [](const CT3RectFlatBuffer<double> &self, const size_t index, const size_t i, const size_t j) -> double {
                return self.data(index)[i * self.number_of_columns() + j];
             },
             "Gets tensor element value.")
        .def("__copy__", [](CT3RectFlatBuffer<double> &self) { return CT3RectFlatBuffer<double>(self); })
        .def("__deepcopy__", [](const CT3RectFlatBuffer<double> &self, py::dict) { return CT3RectFlatBuffer<double>(self); });
}

}  // namespace vlx_general
