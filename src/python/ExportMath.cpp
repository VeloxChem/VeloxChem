#include "ExportMath.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <array>
#include <utility>
#include <vector>

#include "CustomViews.hpp"
#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "MathConst.hpp"
#include "Matrices.hpp"
#include "MatricesFunc.hpp"
#include "Matrix.hpp"
#include "MatrixFunc.hpp"
#include "MolecularBasis.hpp"
#include "SubMatrix.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_math {

// Helper function for CDenseMatrix constructor

auto
CDenseMatrix_from_numpy(const py::array_t<double>& arr) -> std::shared_ptr<CDenseMatrix>
{
    // check dimension

    std::string errdim("DenseMatrix: Expecting a 2D numpy array");

    errors::assertMsgCritical(arr.ndim() == 2, errdim);

    if (arr.data() == nullptr || arr.size() == 0)
    {
        return std::make_shared<CDenseMatrix>();
    }

    // check that the numpy array is c-style contiguous

    std::string errsrc("DenseMatrix: Expecting a contiguous numpy array");

    auto c_style = py::detail::check_flags(arr.ptr(), py::array::c_style);

    auto f_style = py::detail::check_flags(arr.ptr(), py::array::f_style);

    errors::assertMsgCritical(c_style | f_style, errsrc);

    // create CDenseMatrix from numpy array

    auto nrows = static_cast<int>(arr.shape(0));

    auto ncols = static_cast<int>(arr.shape(1));

    CDenseMatrix submat(nrows, ncols);

    auto submat_ptr = submat.values();

    if (c_style)
    {
        std::memcpy(submat_ptr, arr.data(), arr.size() * sizeof(double));
    }
    else if (f_style)
    {
        for (py::ssize_t i = 0; i < arr.shape(0); i++)
        {
            for (py::ssize_t j = 0; j < arr.shape(1); j++)
            {
                submat_ptr[i * arr.shape(1) + j] = arr.data()[j * arr.shape(0) + i];
            }
        }
    }

    return std::make_shared<CDenseMatrix>(submat);
}

auto
export_math(py::module &m) -> void
{
    // exposing enum from Matrix.hpp
    // clang-format off
    py::enum_<mat_t>(m, "mat_t")
        .value("symmetric", mat_t::symmetric)
        .value("antisymmetric", mat_t::antisymmetric)
        .value("general", mat_t::general);
    // clang-format on

    // exposing functions from MathConst.hpp
    m.def("pi_value", &mathconst::pi_value, "Gets PI value.");

    // exposing functions from MatrixFunc.hpp
    m.def(
        "make_matrix",
        [](const CMolecularBasis &basis, const mat_t mat_type) -> std::shared_ptr<CMatrix> {
            return std::make_shared<CMatrix>(matfunc::make_matrix(basis, mat_type));
        },
        "Creates matrix for given basis.");
    m.def(
        "make_matrix",
        [](const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis) -> std::shared_ptr<CMatrix> {
            return std::make_shared<CMatrix>(matfunc::make_matrix(bra_basis, ket_basis));
        },
        "Creates matrix for given pair of bases.");

    // exposing functions from MatrixFunc.hpp
    m.def(
        "make_matrices",
        [](const int order, const CMolecularBasis &basis, const mat_t mtype) -> std::shared_ptr<CMatrices> {
            return std::make_shared<CMatrices>(matfunc::make_matrices(std::array<int, 1>{order}, basis, mtype));
        },
        "Creates matrices for given basis.");
    m.def(
        "make_matrices",
        [](const std::array<int, 2> orders, const CMolecularBasis &basis, const mat_t mtype) -> std::shared_ptr<CMatrices> {
            return std::make_shared<CMatrices>(matfunc::make_matrices(orders, basis, mtype));
        },
        "Creates matrices for given basis.");
    m.def(
        "make_matrices",
        [](const std::array<int, 3> orders, const CMolecularBasis &basis, const mat_t mtype) -> std::shared_ptr<CMatrices> {
            return std::make_shared<CMatrices>(matfunc::make_matrices(orders, basis, mtype));
        },
        "Creates matrices for given basis.");
    m.def(
        "make_matrices",
        [](const std::array<int, 4> orders, const CMolecularBasis &basis, const mat_t mtype) -> std::shared_ptr<CMatrices> {
            return std::make_shared<CMatrices>(matfunc::make_matrices(orders, basis, mtype));
        },
        "Creates matrices for given basis.");
    m.def(
        "make_matrices",
        [](const int order, const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis) -> std::shared_ptr<CMatrices> {
            return std::make_shared<CMatrices>(matfunc::make_matrices(std::array<int, 1>{order}, bra_basis, ket_basis));
        },
        "Creates matrices for given basis.");
    m.def(
        "make_matrices",
        [](const std::array<int, 2> orders, const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis) -> std::shared_ptr<CMatrices> {
            return std::make_shared<CMatrices>(matfunc::make_matrices(orders, bra_basis, ket_basis));
        },
        "Creates matrices for given basis.");
    m.def(
        "make_matrices",
        [](const std::array<int, 3> orders, const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis) -> std::shared_ptr<CMatrices> {
            return std::make_shared<CMatrices>(matfunc::make_matrices(orders, bra_basis, ket_basis));
        },
        "Creates matrices for given basis.");
    m.def(
        "make_matrices",
        [](const std::array<int, 4> orders, const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis) -> std::shared_ptr<CMatrices> {
            return std::make_shared<CMatrices>(matfunc::make_matrices(orders, bra_basis, ket_basis));
        },
        "Creates matrices for given basis.");

    // CSubMatrix class
    PyClass<CSubMatrix>(m, "SubMatrix")
        .def(py::init<>())
        .def(py::init<const CSubMatrix &>())
        .def(py::init<const std::array<size_t, 4> &>())
        .def(py::init<const std::array<size_t, 4> &, const double>())
        .def(py::init<const std::vector<double> &, const std::array<size_t, 4> &>())
        .def(py::pickle([](const CSubMatrix &submat) { return py::make_tuple(submat.get_values(), submat.get_dimensions()); },
                        [](py::tuple t) { return CSubMatrix(t[0].cast<std::vector<double>>(), t[1].cast<std::array<size_t, 4>>()); }))
        .def(
            "at",
            [](CSubMatrix &self, const std::pair<size_t, size_t> &index, const double value) -> void { self.at(index) = value; },
            "Sets value of submatrix element using local indexing scheme.")
        .def(
            "at",
            [](const CSubMatrix &self, const std::pair<size_t, size_t> &index) -> double { return self.at(index); },
            "Gets value of submatrix element using local indexing scheme.")
        .def("set_offsets", &CSubMatrix::set_offsets, "Sets offsets for rows and columns of submatrix in supermatrix.")
        .def(
            "set_values",
            [](CSubMatrix &self, const py::array_t<double> &values) -> void {
                if (values.ndim() == 2)
                {
                    const auto nrows = self.number_of_rows();
                    const auto ncols = self.number_of_columns();
                    if ((nrows == values.shape(0)) && (ncols == values.shape(1)))
                    {
                        std::ranges::for_each(views::rectangular(nrows, ncols),
                                              [&](const auto &index) { self.at(index) = values.at(index.first, index.second); });
                    }
                }
            },
            "Sets values of submatrix using numpy array.")
        .def("zero", &CSubMatrix::zero, "Sets values of submatrix to zero.")
        .def("scale", &CSubMatrix::scale, "Scales values of submatrix by given factor.")
        .def("symmetrize", &CSubMatrix::symmetrize, "Symmetrizes values of square submatrix.")
        .def(
            "to_numpy",
            [](const CSubMatrix &self) -> py::array_t<double> {
                const auto nrows = static_cast<py::ssize_t>(self.number_of_rows());
                const auto ncols = static_cast<py::ssize_t>(self.number_of_columns());
                const auto tdim  = static_cast<py::ssize_t>(sizeof(double));
                return py::array_t<double>(
                    std::vector<py::ssize_t>({nrows, ncols}), std::vector<py::ssize_t>({ncols * tdim, tdim}), self.get_values().data());
            },
            "Converts submatrix to numpy array.")
        .def("get_dimensions", &CSubMatrix::get_dimensions, "Gets dimensions of submatrix.")
        .def("offset_of_rows", &CSubMatrix::offset_of_rows, "Gets offset of rows in submatrix.")
        .def("offset_of_columns", &CSubMatrix::offset_of_columns, "Gets offset of columns in submatrix.")
        .def("number_of_rows", &CSubMatrix::number_of_rows, "Gets number of rows in submatrix.")
        .def("number_of_columns", &CSubMatrix::number_of_columns, "Gets number of columns in submatrix.")
        .def("number_of_elements", &CSubMatrix::number_of_elements, "Gets number of elements in submatrix.")
        .def("is_square", &CSubMatrix::is_square, "Checks if submatrix is square matrix")
        .def("__setitem__", [](CSubMatrix &self, const std::pair<size_t, size_t> &index, const double value) { self[index] = value; })
        .def("__getitem__", [](CSubMatrix &self, const std::pair<size_t, size_t> &index) { return self[index]; })
        .def("__add__", [](const CSubMatrix &self, const CSubMatrix &other) { return self + other; })
        .def("__eq__", [](const CSubMatrix &self, const CSubMatrix &other) { return self == other; })
        .def("__ne__", [](const CSubMatrix &self, const CSubMatrix &other) { return self != other; })
        .def("__copy__", [](const CSubMatrix &self) { return CSubMatrix(self); })
        .def("__deepcopy__", [](const CSubMatrix &self, py::dict) { return CSubMatrix(self); });

    // CMatrix class
    PyClass<CMatrix>(m, "Matrix")
        .def(py::init<>())
        .def(py::init<const std::map<std::pair<int, int>, CSubMatrix> &, const mat_t>())
        .def(py::init<const std::vector<std::pair<int, int>> &, const std::vector<CSubMatrix> &, const mat_t>())
        .def(py::init<const CMatrix &>())
        .def(py::pickle([](const CMatrix &mat) { return py::make_tuple(mat.get_type(), mat.angular_pairs(), mat.sub_matrices()); },
                        [](py::tuple t) {
                            return CMatrix(t[1].cast<std::vector<std::pair<int, int>>>(), t[2].cast<std::vector<CSubMatrix>>(), t[0].cast<mat_t>());
                        }))
        .def("add", py::overload_cast<const CSubMatrix &, const std::pair<int, int> &>(&CMatrix::add), "Adds submatrix to matrix.")
        .def("add", py::overload_cast<const std::array<size_t, 4> &, const std::pair<int, int> &>(&CMatrix::add), "Adds submatrix to matrix.")
        .def("set_type", &CMatrix::set_type, "Sets matrix type.")
        .def("zero", &CMatrix::zero, "Sets values of matrix to zero.")
        .def("scale", &CMatrix::scale, "Scales matrix values by factor.")
        .def("symmetrize", &CMatrix::symmetrize, "Symmetrizes values of diagonal blocks of symmetric matrix.")
        .def(
            "set_values",
            [](CMatrix &self, const py::array_t<double> &values) -> void {
                if (values.ndim() == 2)
                {
                    const auto nrows = self.number_of_rows();
                    const auto ncols = self.number_of_columns();
                    if ((nrows == values.shape(0)) && (ncols == values.shape(1)))
                    {
                        std::ranges::for_each(self.angular_pairs(), [&](const auto &key) {
                            auto       submat = self.sub_matrix(key);
                            const auto dims   = submat->get_dimensions();
                            std::ranges::for_each(views::rectangular(dims[2], dims[3]), [&](const auto &index) {
                                submat->at(index) = values.at(index.first + dims[0], index.second + dims[1]);
                            });
                        });
                    }
                }
            },
            "Sets values of matrix using numpy array.")
        .def("get_type", &CMatrix::get_type, "Gets matrix type.")
        .def("angular_pairs", &CMatrix::angular_pairs, "Gets vector of angular pairs for stored submatrices.")
        .def(
            "submatrix",
            [](const CMatrix &self, const std::pair<int, int> &angpair) -> std::shared_ptr<CSubMatrix> {
                if (auto submat = self.sub_matrix(angpair); submat != nullptr)
                {
                    return std::make_shared<CSubMatrix>(*submat);
                }
                else
                {
                    return std::make_shared<CSubMatrix>();
                }
            },
            "Gets specific submatrix from matrix.")
        .def("is_angular_order", &CMatrix::is_angular_order, "Checks if submatrix with this angular pair is stored in matrix.")
        .def("number_of_rows", &CMatrix::number_of_rows, "Number of rows in matrix.")
        .def("number_of_columns", &CMatrix::number_of_columns, "Number of columns in matrix.")
        .def(
            "full_matrix",
            [](const CMatrix &self) -> std::shared_ptr<CSubMatrix> { return std::make_shared<CSubMatrix>(self.full_matrix()); },
            "Creates full matrix representation of matrix.")
        .def("__add__", [](const CMatrix &self, const CMatrix &other) { return self + other; })
        .def("__eq__", [](const CMatrix &self, const CMatrix &other) { return self == other; })
        .def("__ne__", [](const CMatrix &self, const CMatrix &other) { return self != other; })
        .def("__copy__", [](const CMatrix &self) { return CMatrix(self); })
        .def("__deepcopy__", [](const CMatrix &self, py::dict) { return CMatrix(self); });

    // CMatrices class
    PyClass<CMatrices>(m, "Matrices")
        .def(py::init<>())
        .def(py::init<const std::map<std::string, CMatrix> &>())
        .def(py::init<const CMatrices &>())
        .def("add", py::overload_cast<const CMatrix &, const std::string &>(&CMatrices::add), "Adds matrix to matrix.")
        .def("add", py::overload_cast<const CMatrix &, const int>(&CMatrices::add), "Adds matrix to matrix.")
        .def("zero", &CMatrices::zero, "Sets values of matrices to zero.")
        .def("scale", &CMatrices::scale, "Scales matrices values by factor.")
        .def("symmetrize", &CMatrices::symmetrize, "Symmetrizes values of diagonal blocks of symmetric matrices.")
        .def("keys", &CMatrices::keys, "Gets vector of keys for stored matrices.")
        .def(
            "matrix",
            [](const CMatrices &self, const std::string &label) -> std::shared_ptr<CMatrix> {
                if (auto mat = self.matrix(label); mat != nullptr)
                {
                    return std::make_shared<CMatrix>(*mat);
                }
                else
                {
                    return std::make_shared<CMatrix>();
                }
            },
            "Gets specific matrix from matrices.")
        .def(
            "matrix",
            [](const CMatrices &self, const int key) -> std::shared_ptr<CMatrix> {
                if (auto mat = self.matrix(key); mat != nullptr)
                {
                    return std::make_shared<CMatrix>(*mat);
                }
                else
                {
                    return std::make_shared<CMatrix>();
                }
            },
            "Gets specific matrix from matrices.")
        .def("__eq__", [](const CMatrices &self, const CMatrices &other) { return self == other; })
        .def("__ne__", [](const CMatrices &self, const CMatrices &other) { return self != other; })
        .def("__copy__", [](const CMatrices &self) { return CMatrices(self); })
        .def("__deepcopy__", [](const CMatrices &self, py::dict) { return CMatrices(self); });
}

}  // namespace vlx_math
