//
//  Tabula — custom-recursion molecular-integral machinery.
//  Python (pybind11) bindings.
//

#include "ExportTabula.hpp"

#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <utility>

#include "TabulaDenseMatrix.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace tabula {  // tabula namespace

auto
export_tabula(py::module& m) -> void
{
    // tabula::Symmetry enum class

    py::enum_<Symmetry>(m, "TabulaSymmetry", "The symmetry of a Tabula dense matrix.")
        .value("general", Symmetry::general)
        .value("symmetric", Symmetry::symmetric)
        .value("antisymmetric", Symmetry::antisymmetric);

    // tabula::DenseMatrix class

    py::class_<DenseMatrix, std::shared_ptr<DenseMatrix>>(m, "TabulaDenseMatrix", "A dense, row-major Tabula matrix.")
        .def(py::init<>())
        .def(py::init<std::size_t, std::size_t>(), "rows"_a, "columns"_a)
        .def(py::init<std::size_t, std::size_t, Symmetry>(), "rows"_a, "columns"_a, "symmetry"_a)
        .def("rows", &DenseMatrix::rows, "Gets the number of rows.")
        .def("columns", &DenseMatrix::columns, "Gets the number of columns.")
        .def("symmetry", &DenseMatrix::symmetry, "Gets the matrix symmetry.")
        .def("zero", &DenseMatrix::zero, "Sets all elements to zero.")
        .def("scale", &DenseMatrix::scale, "Scales all elements by a factor.", "factor"_a)
        .def("symmetrize", &DenseMatrix::symmetrize, "Mirrors the upper triangle into the lower one per the matrix symmetry.")
        .def(
            "__getitem__",
            [](const DenseMatrix& self, const std::pair<std::size_t, std::size_t>& index) -> double {
                return self(index.first, index.second);
            },
            "index"_a)
        .def(
            "__setitem__",
            [](DenseMatrix& self, const std::pair<std::size_t, std::size_t>& index, const double value) -> void {
                self(index.first, index.second) = value;
            },
            "index"_a,
            "value"_a)
        .def(
            "to_numpy",
            [](py::object obj) -> py::array_t<double> {
                auto& self = obj.cast<DenseMatrix&>();
                return py::array_t<double>({self.rows(), self.columns()},
                                           {self.columns() * sizeof(double), sizeof(double)},
                                           self.values(),
                                           obj);
            },
            "Returns a zero-copy NumPy view of the matrix.");
}

}  // namespace tabula
