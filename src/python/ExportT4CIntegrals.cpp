#include "ExportT4CIntegrals.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "FockMatrices.hpp"
#include "FockMatrix.hpp"
#include "FockType.hpp"

namespace vlx_t4cintegrals {  // vlx_t4cintegrals namespace

// Exports classes/functions in src/t4c_* to python

void
export_t4cintegrals(py::module& m)
{
    // exposing enum from FockType.hpp

    py::enum_<fock_t>(m, "fock_t")
        .value("restjk", fock_t::restjk)
        .value("restjkx", fock_t::restjkx)
        .value("restj", fock_t::restj)
        .value("restk", fock_t::restk)
        .value("restkx", fock_t::restkx)
        .value("rgenjk", fock_t::rgenjk)
        .value("rgenjkx", fock_t::rgenjkx)
        .value("rgenj", fock_t::rgenj)
        .value("rgenk", fock_t::rgenk)
        .value("rgenkx", fock_t::rgenkx)
        .value("unrestjk", fock_t::unrestjk)
        .value("unrestj", fock_t::unrestj)
        .value("unrestjkx", fock_t::unrestjkx);

    // CFockMatrix class

    PyClass<CFockMatrix>(m, "FockMatrix")
        .def(py::init<>())
        .def(py::init<const CMatrix&, const double, const fock_t>())
        .def(py::init<const CFockMatrix&>())
        .def("set_exchange_scale", &CFockMatrix::setExchangeScale, "Sets scaling factor of exchange contribution to Fock matrix.")
        .def("set_fock_type", &CFockMatrix::setFockType, "Sets Fock matrix type.")
        .def("zero", &CFockMatrix::zero, "Sets values of Fock matrix to zero.")
        .def("get_exchange_scale", &CFockMatrix::getExchangeScale, "Gets scaling factor of exchange contribution to Fock matrix.")
        .def("get_fock_type", &CFockMatrix::getFockType, "Gets Fock matrix type.")
        .def("get_storage_type", &CFockMatrix::getStorageType, "Gets Fock matrix storage type.")
        .def("get_angular_pairs", &CFockMatrix::getAngularPairs, "Gets vector of angular pairs for stored submatrices.")
        .def(
            "get_submatrix",
            [](const CFockMatrix& self, const T2Pair& angpair) -> std::shared_ptr<CSubMatrix> {
                if (auto submat = self.getSubMatrix(angpair); submat != nullptr)
                {
                    return std::make_shared<CSubMatrix>(*submat);
                }
                else
                {
                    return std::make_shared<CSubMatrix>();
                }
            },
            "Gets specific submatrix from Fock matrix.")
        .def(
            "get_matrix",
            [](const CFockMatrix& self) -> std::shared_ptr<CMatrix> {
                if (auto mat = self.getMatrix(); mat != nullptr)
                {
                    return std::make_shared<CMatrix>(*mat);
                }
                else
                {
                    return std::make_shared<CMatrix>();
                }
            },
            "Gets values matrix from Fock matrix.")
        .def("is_angular_order", &CFockMatrix::isAngularOrder, "Checks if submatrix with this angular pair is stored in Fock matrix.")
        .def("number_of_rows", &CFockMatrix::getNumberOfRows, "Number of rows in Fock matrix.")
        .def("number_of_columns", &CFockMatrix::getNumberOfColumns, "Number of columns in Fock matrix.")
        .def(
            "get_full_matrix",
            [](const CFockMatrix& self) -> std::shared_ptr<CSubMatrix> { return std::make_shared<CSubMatrix>(self.getFullMatrix()); },
            "Creates full matrix representation of Fock matrix.");

    // CFockMatrices class

    PyClass<CFockMatrices>(m, "FockMatrices")
        .def(py::init<>())
        .def(py::init<const std::vector<CFockMatrix>&>())
        .def(py::init<const CFockMatrices&>())
        .def("add", &CFockMatrices::add, "Adds Fock matrix to Fock matrices.")
        .def("zero", &CFockMatrices::zero, "Sets values of Fock matrices to zero.")
        .def(
            "get_matrix",
            [](const CFockMatrices& self, const int64_t index) -> std::shared_ptr<CFockMatrix> {
                if (auto mat = self.getMatrix(index); mat != nullptr)
                {
                    return std::make_shared<CFockMatrix>(*mat);
                }
                else
                {
                    return std::make_shared<CFockMatrix>();
                }
            },
            "Gets specific Fock matrix from Fock matrices.")
        .def("number_of_matrices", &CFockMatrices::getNumberOfMatrices, "Gets number of Fock matrices.");
}

}  // namespace vlx_t4cintegrals
