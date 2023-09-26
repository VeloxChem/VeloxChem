#include "ExportMath.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

#include "AngularMomentum.hpp"
#include "CantorFunc.hpp"
#include "DenseMatrix.hpp"
#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "MathConst.hpp"
#include "Matrices.hpp"
#include "Matrix.hpp"
#include "MatrixFunc.hpp"
#include "MatrixIndex.hpp"
#include "MatrixType.hpp"
#include "SubMatrix.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_math {  // vlx_math namespace

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

    std::vector<double> vec(arr.size());

    if (c_style)
    {
        std::memcpy(vec.data(), arr.data(), arr.size() * sizeof(double));
    }
    else if (f_style)
    {
        for (py::ssize_t i = 0; i < arr.shape(0); i++)
        {
            for (py::ssize_t j = 0; j < arr.shape(1); j++)
            {
                vec.data()[i * arr.shape(1) + j] = arr.data()[j * arr.shape(0) + i];
            }
        }
    }

    int64_t nrows = static_cast<int64_t>(arr.shape(0));

    int64_t ncols = static_cast<int64_t>(arr.shape(1));

    return std::make_shared<CDenseMatrix>(vec, nrows, ncols);
}

// Exports classes/functions in src/math to python

auto
export_math(py::module& m) -> void
{
    // exposing enum from MatrixType.hpp

    py::enum_<mat_t>(m, "mat_t").value("symm", mat_t::symm).value("antisymm", mat_t::antisymm).value("gen", mat_t::gen);

    // exposing functions from MathConst.hpp

    m.def("get_pi", &mathconst::getPiValue, "Gets PI value.");

    // exposing functions from MathIndex.hpp

    m.def("uplo_index", &mathfunc::uplo_index, "Gets index of upper triangular matrix.");

    // exposing functions from CantorFunc.hpp

    m.def("cantor_index", &mathfunc::getCantorIndex, "Computes Cantor index for pair of non-negative numbers.");

    m.def("cantor_pair", &mathfunc::getCantorPair, "Extractes pair of non-negative numbers from Cantor index.");

    // exposing functions from AngularMomentum.hpp

    m.def("to_spherical_components",
          py::overload_cast<int64_t>(angmom::to_SphericalComponents),
          "Gets number of spherical components of given angular momentum.");

    m.def("to_spherical_components",
          py::overload_cast<int64_t, int64_t>(angmom::to_SphericalComponents),
          "Gets number of spherical components of given pair of angular momentums.");

    m.def("to_cartesian_components",
          py::overload_cast<int64_t>(angmom::to_CartesianComponents),
          "Gets number of Cartesian components of given angular momentum.");

    m.def("to_cartesian_components",
          py::overload_cast<int64_t, int64_t>(angmom::to_CartesianComponents),
          "Gets number of Cartesian components of given pair of angular momentums.");

    m.def("angular_component_to_str", &angmom::getStringOfAngularMomentum, "Gets string of angular momentum component.");

    // exposing functions from MatrixFunc.hpp

    m.def(
        "make_matrix",
        [](const CMolecularBasis& basis, const mat_t mat_type) -> std::shared_ptr<CMatrix> {
            return std::make_shared<CMatrix>(matfunc::makeMatrix(basis, mat_type));
        },
        "Creates matrix for given basis.");

    m.def(
        "make_matrix",
        [](const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> std::shared_ptr<CMatrix> {
            return std::make_shared<CMatrix>(matfunc::makeMatrix(bra_basis, ket_basis));
        },
        "Creates matrix for given pair of bases.");

    // CSubMatrix class

    PyClass<CSubMatrix>(m, "SubMatrix")
        .def(py::init<>())
        .def(py::init<const CSubMatrix&>())
        .def(py::init<const T4Index&>())
        .def(py::init<const std::vector<double>&, const T4Index&>())
        .def(
            "at",
            [](CSubMatrix& self, const int64_t irow, const int64_t icol, const double value) -> void { self.at(irow, icol, true) = value; },
            "Sets value of submatrix element using supermatrix indexing scheme.")
        .def(
            "at",
            [](const CSubMatrix& self, const int64_t irow, const int64_t icol) -> double { return self.at(irow, icol, true); },
            "Gets value of submatrix element using supermatrix indexing scheme.")
        .def(
            "set_value",
            [](CSubMatrix& self, const int64_t irow, const int64_t icol, const double value) -> void { self.at(irow, icol, false) = value; },
            "Sets value of submatrix element.")
        .def(
            "get_value",
            [](const CSubMatrix& self, const int64_t irow, const int64_t icol) -> double { return self.at(irow, icol, false); },
            "Gets value of submatrix element.")
        .def("set_offsets", &CSubMatrix::setOffsets, "Sets offsets for rows and columns of submatrix in supermatrix.")
        .def(
            "set_values",
            [](CSubMatrix& self, const py::array_t<double>& values) -> void {
                if (values.ndim() == 2)
                {
                    const auto nrows = static_cast<py::ssize_t>(self.getNumberOfRows());

                    const auto ncols = static_cast<py::ssize_t>(self.getNumberOfColumns());

                    if ((nrows == values.shape(0)) && (ncols == values.shape(1)))
                    {
                        std::memcpy(self.getData(), values.data(), nrows * ncols * sizeof(double));
                    }
                }
            },
            "Sets values of submatrix using numpy array.")
        .def("zero", &CSubMatrix::zero, "Sets values of submatrix to zero.")
        .def(
            "to_numpy",
            [](const CSubMatrix& self) -> py::array_t<double> {
                const auto nrows = static_cast<py::ssize_t>(self.getNumberOfRows());

                const auto ncols = static_cast<py::ssize_t>(self.getNumberOfColumns());

                const auto tdim = static_cast<py::ssize_t>(sizeof(double));

                return py::array_t<double>(std::vector<py::ssize_t>({nrows, ncols}), std::vector<py::ssize_t>({ncols * tdim, tdim}), self.getData());
            },
            "Converts submatrix to numpy array.")
        .def("get_dimensions", &CSubMatrix::getDimensions, "Gets dimensions of submatrix.")
        .def("offset_of_rows", &CSubMatrix::getOffsetOfRows, "Gets offset of rows in submatrix.")
        .def("offset_of_columns", &CSubMatrix::getOffsetOfColumns, "Gets offset of columns in submatrix.")
        .def("number_of_rows", &CSubMatrix::getNumberOfRows, "Gets number of rows in submatrix.")
        .def("number_of_columns", &CSubMatrix::getNumberOfColumns, "Gets number of columns in submatrix.");

    // CMatrix class

    PyClass<CMatrix>(m, "Matrix")
        .def(py::init<>())
        .def(py::init<const std::map<T2Pair, CSubMatrix>&, const mat_t>())
        .def(py::init<const CMatrix&>())
        .def("add", py::overload_cast<const CSubMatrix&, const T2Pair&>(&CMatrix::add), "Adds submatrix to matrix.")
        .def("add", py::overload_cast<const T4Index&, const T2Pair&>(&CMatrix::add), "Adds submatrix to matrix.")
        .def("set_type", &CMatrix::setType, "Sets matrix type.")
        .def("zero", &CMatrix::zero, "Sets values of matrix to zero.")
        .def("get_type", &CMatrix::getType, "Gets matrix type.")
        .def("get_angular_pairs", &CMatrix::getAngularPairs, "Gets vector of angular pairs for stored submatrices.")
        .def(
            "get_submatrix",
            [](const CMatrix& self, const T2Pair& angpair) -> std::shared_ptr<CSubMatrix> {
                if (auto submat = self.getSubMatrix(angpair); submat != nullptr)
                {
                    return std::make_shared<CSubMatrix>(*submat);
                }
                else
                {
                    return std::make_shared<CSubMatrix>();
                }
            },
            "Gets specific submatrix from matrix.")
        .def("is_angular_order", &CMatrix::isAngularOrder, "Checks if submatrix with this angular pair is stored in matrix.")
        .def("number_of_rows", &CMatrix::getNumberOfRows, "Number of rows in matrix.")
        .def("number_of_columns", &CMatrix::getNumberOfColumns, "Number of columns in matrix.")
        .def(
            "get_full_matrix",
            [](const CMatrix& self) -> std::shared_ptr<CSubMatrix> { return std::make_shared<CSubMatrix>(self.getFullMatrix()); },
            "Creates full matrix representation of matrix.");

    // CMatrices class

    PyClass<CMatrices>(m, "Matrices")
        .def(py::init<>())
        .def(py::init<const std::map<int64_t, CMatrix>&>())
        .def(py::init<const CMatrix&, const std::vector<int64_t>&>())
        .def(py::init<const CMatrices&>())
        .def("add", py::overload_cast<const CMatrix&, const int64_t>(&CMatrices::add), "Adds matrix to matrices.")
        .def("add", py::overload_cast<const CMatrix&, const std::string&>(&CMatrices::add), "Adds matrix to matrices.")
        .def("zero", &CMatrices::zero, "Sets values of matrices to zero.")
        .def("get_keys", &CMatrices::getKeys, "Gets keys of each matrix in matrices.")
        .def(
            "get_matrix",
            [](const CMatrices& self, const int64_t key) -> std::shared_ptr<CMatrix> {
                if (auto mat = self.getMatrix(key); mat != nullptr)
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
            "get_matrix",
            [](const CMatrices& self, const std::string& label) -> std::shared_ptr<CMatrix> {
                if (auto mat = self.getMatrix(label); mat != nullptr)
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
            "get_matrix",
            [](const CMatrices& self, const int64_t atom, const std::string& label) -> std::shared_ptr<CMatrix> {
                if (auto mat = self.getMatrix(atom, label); mat != nullptr)
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
            "get_matrix",
            [](const CMatrices& self, const int64_t atom_a, const std::string& label_a, const int64_t atom_b, const std::string& label_b)
                -> std::shared_ptr<CMatrix> {
                if (auto mat = self.getMatrix(atom_a, label_a, atom_b, label_b); mat != nullptr)
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
            "get_matrix",
            [](const CMatrices&   self,
               const int64_t      atom_a,
               const std::string& label_a,
               const int64_t      atom_b,
               const std::string& label_b,
               const int64_t      atom_c,
               const std::string& label_c) -> std::shared_ptr<CMatrix> {
                if (auto mat = self.getMatrix(atom_a, label_a, atom_b, label_b, atom_c, label_c); mat != nullptr)
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
            "get_matrix",
            [](const CMatrices&   self,
               const int64_t      atom_a,
               const std::string& label_a,
               const int64_t      atom_b,
               const std::string& label_b,
               const int64_t      atom_c,
               const std::string& label_c,
               const int64_t      atom_d,
               const std::string& label_d) -> std::shared_ptr<CMatrix> {
                if (auto mat = self.getMatrix(atom_a, label_a, atom_b, label_b, atom_c, label_c, atom_d, label_d); mat != nullptr)
                {
                    return std::make_shared<CMatrix>(*mat);
                }
                else
                {
                    return std::make_shared<CMatrix>();
                }
            },
            "Gets specific matrix from matrices.");

    // CDenseMatrix class

    PyClass<CDenseMatrix>(m, "DenseMatrix")
        .def(py::init<>())
        .def(py::init<const int64_t, const int64_t>())
        .def(py::init<const CDenseMatrix&>())
        .def(py::init(&CDenseMatrix_from_numpy))
        .def("number_of_rows", &CDenseMatrix::getNumberOfRows, "Gets number of rows in dense matrix.")
        .def("number_of_columns", &CDenseMatrix::getNumberOfColumns, "Gets number of columns in dense matrix.")
        .def("symmetrize", &CDenseMatrix::symmetrize, "Symmetrizes elements of square matrix: a_ij = a_ji = (a_ij + a_ji).")
        .def("slice",
             &CDenseMatrix::slice,
             "Creates dense matrix object by slicing columns at selected position from this dense matrix object.",
             "i_column"_a,
             "n_columns"_a)
        .def(
            "to_numpy",
            [](const CDenseMatrix& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.values(), std::vector<int64_t>{self.getNumberOfRows(), self.getNumberOfColumns()});
            },
            "Converts DenseMatrix to numpy array.")
        .def(py::self == py::self);
}

}  // namespace vlx_math
