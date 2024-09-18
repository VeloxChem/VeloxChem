#include "ExportT4CIntegrals.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "FockDriver.hpp"
#include "FockGeomX000Driver.hpp"
#include "T4CScreener.hpp"
#include "T4CUtils.hpp"

namespace vlx_t4cintegrals {  // vlx_t4cintegrals namespace

// Exports classes/functions in src/t4c_* to python

void
export_t4cintegrals(py::module& m)
{
    // FockDriver class
    PyClass<CFockDriver>(m, "FockDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CFockDriver&     fock_drv,
               const CMolecularBasis& basis,
               const CMolecule&       molecule,
               const CMatrix&         density,
               const std::string&     label,
               const double           exchange_factor,
               const double           omega) -> std::shared_ptr<CMatrix> {
                return std::make_shared<CMatrix>(fock_drv.compute(basis, molecule, density, label, exchange_factor, omega));
            },
            "Computes single Fock matrix of requested type for given molecule and basis.")
        .def(
            "compute",
            [](const CFockDriver&  fock_drv,
               const CT4CScreener& screener,
               const CMatrix&      density,
               const std::string&  label,
               const double        exchange_factor,
               const double        omega,
               const int           ithreshold) -> std::shared_ptr<CMatrix> {
                return std::make_shared<CMatrix>(fock_drv.compute(screener, density, label, exchange_factor, omega, ithreshold));
            },
            "Computes single Fock matrix of requested type for two-electron integrals screener.")
        .def(
            "compute",
            [](const CFockDriver&              fock_drv,
               const CT4CScreener&             screener,
               const CMatrices&                densities,
               const std::vector<std::string>& labels,
               const double                    exchange_factor,
               const double                    omega,
               const int                       ithreshold) -> std::shared_ptr<CMatrices> {
                return std::make_shared<CMatrices>(fock_drv.compute(screener, densities, labels, exchange_factor, omega, ithreshold));
            },
            "Computes Fock matrices of requested type for two-electron integrals screener.")
        .def(
            "compute",
            [](const CFockDriver&  fock_drv,
               const CT4CScreener& screener,
               const int           rank,
               const int           nodes,
               const CMatrix&      density,
               const std::string&  label,
               const double        exchange_factor,
               const double        omega,
               const int           ithreshold) -> std::shared_ptr<CMatrix> {
                return std::make_shared<CMatrix>(fock_drv.compute(screener, rank, nodes, density, label, exchange_factor, omega, ithreshold));
            },
            "Computes single Fock matrix of requested type for two-electron integrals screener.")
        .def(
            "compute",
            [](const CFockDriver&              fock_drv,
               const CT4CScreener&             screener,
               const int                       rank,
               const int                       nodes,
               const CMatrices&                densities,
               const std::vector<std::string>& labels,
               const double                    exchange_factor,
               const double                    omega,
               const int                       ithreshold) -> std::shared_ptr<CMatrices> {
                return std::make_shared<CMatrices>(fock_drv.compute(screener, rank, nodes, densities, labels, exchange_factor, omega, ithreshold));
            },
            "Computes Fock matrices of requested type for two-electron integrals screener.");

    // CT4CScreener class
    PyClass<CT4CScreener>(m, "T4CScreener")
        .def(py::init<>())
        .def(py::init<const CT4CScreener&>())
        .def(py::init<const std::vector<CBlockedGtoPairBlock>&>())
        .def(py::pickle([](const CT4CScreener& screener) { return py::make_tuple(screener.gto_pair_blocks()); },
                        [](py::tuple t) { return CT4CScreener(t[0].cast<std::vector<CBlockedGtoPairBlock>>()); }))
        .def("partition", &CT4CScreener::partition, "Partition basis funtion pairs blocks for given molecule and basis.")
        .def("gto_pair_blocks", &CT4CScreener::gto_pair_blocks, "Gets vector of blocked basis function pairs blocks.")
        .def("__eq__", [](const CT4CScreener& self, const CT4CScreener& other) { return self == other; })
        .def("__ne__", [](const CT4CScreener& self, const CT4CScreener& other) { return self != other; })
        .def("__copy__", [](const CT4CScreener& self) { return CT4CScreener(self); })
        .def("__deepcopy__", [](const CT4CScreener& self, py::dict) { return CT4CScreener(self); });

    // CFockGeom1000Driver class
    PyClass<CFockGeomX000Driver<1>>(m, "FockGeom1000Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CFockGeomX000Driver<1>& fock_drv,
               const CMolecularBasis&        basis,
               const CMolecule&              molecule,
               const CMatrix&                density,
               const int                     iatom,
               const std::string&            label,
               const double                  exchange_factor,
               const double                  omega) -> std::shared_ptr<CMatrices> {
                return std::make_shared<CMatrices>(fock_drv.compute(basis, molecule, density, iatom, label, exchange_factor, omega));
            },
            "Computes gradient of Fock matrix of requested type for given molecule and basis.");
}

}  // namespace vlx_t4cintegrals
