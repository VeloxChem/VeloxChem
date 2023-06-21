#include "ExportMath.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "MathConst.hpp"
#include "AngularMomentum.hpp"
#include "SubMatrix.hpp"
#include "MatrixType.hpp"
#include "Matrix.hpp"
#include "MatrixFunc.hpp"

namespace vlx_math {  // vlx_math namespace

// Exports classes/functions in src/math to python

auto
export_math(py::module& m) -> void
{
    // exposing enum from MatrixType.hpp
    
    py::enum_<mat_t>(m, "mat_t")
        .value("symm",     mat_t::symm)
        .value("antisymm", mat_t::antisymm)
        .value("gen",      mat_t::gen);
    
    // exposing functions from MathConst.hpp

    m.def("get_pi",
          &mathconst::getPiValue,
          "Gets PI value.");
    
    // exposing functions from AngularMomentum.hpp
    
    m.def("to_spherical_components",
          py::overload_cast<int64_t>(angmom::to_SphericalComponents),
          "Gets number of spherical components of given angular momentum.");
    
    m.def("to_spherical_components",
          py::overload_cast<int64_t,int64_t>(angmom::to_SphericalComponents),
          "Gets number of spherical components of given pair of angular momentums.");
    
    m.def("to_cartesian_components",
          py::overload_cast<int64_t>(angmom::to_CartesianComponents),
          "Gets number of Cartesian components of given angular momentum.");
    
    m.def("to_cartesian_components",
          py::overload_cast<int64_t,int64_t>(angmom::to_CartesianComponents),
          "Gets number of Cartesian components of given pair of angular momentums.");
    
    m.def("angular_component_to_str",
          &angmom::getStringOfAngularMomentum,
          "Gets string of angular momentum component.");
    
    // exposing functions from MatrixFunc.hpp

    m.def("make_matrix",
          [](const CMolecularBasis& basis,
             const mat_t            mat_type) -> std::shared_ptr<CMatrix>
          {
             return std::make_shared<CMatrix>(matfunc::makeMatrix(basis, mat_type));
          },
          "Creates matrix for given basis.");
    
    m.def("make_matrix",
          [](const CMolecularBasis& bra_basis,
             const CMolecularBasis& ket_basis) -> std::shared_ptr<CMatrix>
          {
              return std::make_shared<CMatrix>(matfunc::makeMatrix(bra_basis, ket_basis));
          },
          "Creates matrix for given pair of bases.");
    
    // CSubMatrix class

    PyClass<CSubMatrix>(m, "SubMatrix")
        .def(py::init<>())
        .def(py::init<const CSubMatrix&>())
        .def(py::init<const T4Index&>())
        .def(py::init<const std::vector<double>&,
                      const T4Index&>())
        .def("at",
             [](      CSubMatrix& self,
                const int64_t     irow,
                const int64_t     icol,
                const double      value) -> void
             {
                self.at(irow, icol, true) = value;
             },
             "Sets value of submatrix element using supermatrix indexing scheme.")
        .def("at",
             [](const CSubMatrix& self,
                const int64_t     irow,
                const int64_t     icol) -> double
             {
                return self.at(irow, icol, true);
             },
             "Gets value of submatrix element using supermatrix indexing scheme.")
        .def("set_value",
             [](      CSubMatrix& self,
                const int64_t     irow,
                const int64_t     icol,
                const double      value) -> void
             {
                self.at(irow, icol, false) = value;
             },
             "Sets value of submatrix element.")
        .def("get_value",
             [](const CSubMatrix& self,
                const int64_t     irow,
                const int64_t     icol) -> double
             {
                return self.at(irow, icol, false);
             },
             "Gets value of submatrix element.")
        .def("set_offsets",
             &CSubMatrix::setOffsets,
             "Sets offsets for rows and columns of submatrix in supermatrix.")
        .def("set_values",
             [](      CSubMatrix& self,
                const py::array_t<double>& values) -> void
             {
                if (values.ndim() == 2)
                {
                    const auto nrows = static_cast<py::ssize_t>(self.getNumberOfRows());
                
                    const auto ncols = static_cast<py::ssize_t>(self.getNumberOfColumns());
                
                    if ((nrows == values.shape(0)) &&  (ncols == values.shape(1)))
                    {
                        std::memcpy(self.getData(), values.data(), nrows * ncols * sizeof(double));
                    }
                }
             },
             "Sets values of submatrix using numpy array.")
        .def("to_numpy",
             [](const CSubMatrix& self) -> py::array_t<double>
             {
                const auto nrows = static_cast<py::ssize_t>(self.getNumberOfRows());
                
                const auto ncols = static_cast<py::ssize_t>(self.getNumberOfColumns());
            
                const auto tdim = static_cast<py::ssize_t>(sizeof(double));
            
                return py::array_t<double>(std::vector<py::ssize_t>({nrows, ncols}),
                                           std::vector<py::ssize_t>({ncols * tdim, tdim}),
                                           self.getData());
             },
            "Converts submatrix to numpy array.")
        .def("get_dimensions",
             &CSubMatrix::getDimensions,
             "Gets dimensions ofsubmatrix.")
        .def("offset_of_rows",
             &CSubMatrix::getOffsetOfRows,
             "Gets offset of rows in submatrix.")
        .def("offset_of_columns",
             &CSubMatrix::getOffsetOfColumns,
             "Gets offset of columns in submatrix.")
        .def("number_of_rows",
             &CSubMatrix::getNumberOfRows,
             "Gets number of rows in submatrix.")
        .def("number_of_columns",
             &CSubMatrix::getNumberOfColumns,
             "Gets number of columns in submatrix.");
    
    // CMatrix class

    PyClass<CMatrix>(m, "Matrix")
        .def(py::init<>())
        .def(py::init<const std::map<T2Pair, CSubMatrix>&,
                      const mat_t>())
        .def(py::init<CMatrix&>())
        .def("add",
             py::overload_cast<const CSubMatrix&,
                               const T2Pair&>(&CMatrix::add),
             "Adds submatrix to matrix.")
        .def("add",
             py::overload_cast<const T4Index&,
                               const T2Pair&>(&CMatrix::add),
             "Adds submatrix to matrix.")
        .def("set_type",
             &CMatrix::setType,
             "Sets matrix type.")
        .def("get_type",
             &CMatrix::getType,
             "Gets matrix type.")
        .def("get_angular_pairs",
              &CMatrix::getAngularPairs,
              "Gets vector of angular pairs for stored submatrices.")
        .def("get_submatrix",
             [](const CMatrix& self,
                const T2Pair&  angpair) -> std::shared_ptr<CSubMatrix>
             {
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
        .def("is_angular_order",
              &CMatrix::isAngularOrder,
             "Checks if submatrix with this angular pair is stored in matrix.")
        .def("number_of_rows",
             &CMatrix::getNumberOfRows,
             "Number of rows in matrix.")
        .def("number_of_columns",
             &CMatrix::getNumberOfColumns,
             "Number of columns in matrix.")
        .def("get_full_matrix",
             [](const CMatrix& self) -> std::shared_ptr<CSubMatrix>
             {
                return std::make_shared<CSubMatrix>(self.getFullMatrix());
             },
             "Creates full matrix representation of matrix.");

    // ...
}

}  // namespace vlx_math 
