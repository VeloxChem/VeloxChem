#include "ExportOrbdata.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>
#include <vector>

#include "AtomBasis.hpp"
#include "BasisFunction.hpp"
#include "GtoBlock.hpp"
#include "GtoFunc.hpp"
#include "MolecularBasis.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_orbdata {  // vlx_orbdata namespace

void
export_orbdata(py::module &m)
{
    // exposing functions from GtoFunc.hpp
    m.def("make_gto_blocks",
          py::overload_cast<const CMolecularBasis &, const CMolecule &>(&gtofunc::make_gto_blocks),
          "Creates vector of basis functions blocks for given basis and "
          "molecule.");
    m.def("make_gto_blocks",
          py::overload_cast<const CMolecularBasis &, const CMolecule &, const std::vector<int> &>(&gtofunc::make_gto_blocks),
          "Creates vector of basis functions blocks for selected atoms in given "
          "basis and molecule.");

    // CBasisFunction class
    PyClass<CBasisFunction>(m, "BasisFunction")
        .def(py::init<>())
        .def(py::init<const CBasisFunction &>())
        .def(py::init<const std::vector<double> &, const std::vector<double> &, const int>())
        .def(py::pickle(
            [](const CBasisFunction &bf) { return py::make_tuple(bf.get_exponents(), bf.get_normalization_factors(), bf.get_angular_momentum()); },
            [](py::tuple t) { return CBasisFunction(t[0].cast<std::vector<double>>(), t[1].cast<std::vector<double>>(), t[2].cast<int>()); }))
        .def("set_exponents", &CBasisFunction::set_exponents, "Sets exponents of basis function.")
        .def("set_normalization_factors", &CBasisFunction::set_normalization_factors, "Gets name of chemical element.")
        .def("set_angular_momentum", &CBasisFunction::set_angular_momentum, "Sets angular momentum of basis function.")
        .def("add",
             &CBasisFunction::add,
             "Add primitive i.e. exponent and normalization factor to basis "
             "function.")
        .def("normalize", &CBasisFunction::normalize, "Normalizes primitive GTOs in basis function.")
        .def("get_exponents", &CBasisFunction::get_exponents, "Gets vector of exponents in basis function.")
        .def("get_normalization_factors", &CBasisFunction::get_normalization_factors, "Gets vector of normalization factors in basis function.")
        .def("get_angular_momentum", &CBasisFunction::get_angular_momentum, "Gets angular momentum of basis function.")
        .def("number_of_primitives", &CBasisFunction::number_of_primitive_functions, "Gets number of primitives in basis function.")
        .def("__eq__", [](const CBasisFunction &self, const CBasisFunction &other) { return self == other; })
        .def("__copy__", [](const CBasisFunction &self) { return CBasisFunction(self); })
        .def("__deepcopy__", [](const CBasisFunction &self, py::dict) { return CBasisFunction(self); });

    // CAtomBasis class
    PyClass<CAtomBasis>(m, "AtomBasis")
        .def(py::init<>())
        .def(py::init<const CAtomBasis &>())
        .def(py::init<const std::vector<CBasisFunction> &, const std::string &, const std::string &, const int>())
        .def(py::pickle(
            [](const CAtomBasis &abas) {
                return py::make_tuple(abas.basis_functions(), abas.get_name(), abas.get_ecp_label(), abas.get_identifier());
            },
            [](py::tuple t) {
                return CAtomBasis(t[0].cast<std::vector<CBasisFunction>>(), t[1].cast<std::string>(), t[2].cast<std::string>(), t[3].cast<int>());
            }))
        .def("set_identifier", &CAtomBasis::set_identifier, "Sets identifier of atom basis.")
        .def("set_name", &CAtomBasis::set_name, "Sets name of atom basis.")
        .def("set_ecp_label", &CAtomBasis::set_ecp_label, "Sets effective core potential label of atom basis.")
        .def("add", &CAtomBasis::add, "Adds basis function to atom basis.")
        .def("reduce_to_valence_basis", &CAtomBasis::reduce_to_valence_basis, "Reduces atom basis to it's valence only form.")
        .def("get_basis_functions", py::overload_cast<>(&CAtomBasis::basis_functions, py::const_), "Gets GTOs.")
        .def("get_basis_functions",
             py::overload_cast<const int>(&CAtomBasis::basis_functions, py::const_),
             "Gets GTOs with specific angular momentum.")
        .def("get_basis_functions",
             py::overload_cast<const int, const size_t>(&CAtomBasis::basis_functions, py::const_),
             "Gets GTOs with specific angular momentum and number of primitives.")
        .def("get_identifier", &CAtomBasis::get_identifier, "Gets identifier of atom basis.")
        .def("get_name", &CAtomBasis::get_name, "Gets name of atom basis.")
        .def("get_ecp_label", &CAtomBasis::get_ecp_label, "Gets effective core potential label of atom basis.")
        .def("need_ecp", &CAtomBasis::need_ecp, "Checks if atom basis requires effective core potential.")
        .def("max_angular_momentum", &CAtomBasis::max_angular_momentum, "Gets maximum angular momentum in atom basis.")
        .def("number_of_basis_functions",
             py::overload_cast<const int>(&CAtomBasis::number_of_basis_functions, py::const_),
             "Gets number of GTOs with specific angular momentum.")
        .def("number_of_basis_functions",
             py::overload_cast<const int, const size_t>(&CAtomBasis::number_of_basis_functions, py::const_),
             "Gets number of GTOs with specific angular momentum and number of "
             "primitives.")
        .def("number_of_primitive_basis_functions",
             &CAtomBasis::number_of_primitive_functions,
             "Gets number of primitive GTOs with specific angular momentum.")
        .def("contraction_depths", &CAtomBasis::contraction_depths, "Gets contraction depths of GTOs with specific angular momentum.")
        .def("contraction_str", &CAtomBasis::contraction_string, "Gets contraction string of atom basis.")
        .def("primitives_str", &CAtomBasis::primitives_string, "Gets primitive GTOs string of atom basis.")
        .def("__eq__", [](const CAtomBasis &self, const CAtomBasis &other) { return self == other; })
        .def("__copy__", [](const CAtomBasis &self) { return CAtomBasis(self); })
        .def("__deepcopy__", [](const CAtomBasis &self, py::dict) { return CAtomBasis(self); });

    // CMolecularBasis class
    PyClass<CMolecularBasis>(m, "MolecularBasis")
        .def(py::init<>())
        .def(py::init<const CMolecularBasis &>())
        .def(py::init<const std::vector<CAtomBasis> &, const std::vector<int> &>())
        .def(py::pickle([](const CMolecularBasis &mbas) { return py::make_tuple(mbas.basis_sets(), mbas.basis_sets_indices()); },
                        [](py::tuple t) { return CMolecularBasis(t[0].cast<std::vector<CAtomBasis>>(), t[1].cast<std::vector<int>>()); }))
        .def("add", &CMolecularBasis::add, "Adds atomic basis to molecular basis.")
        .def("slice", &CMolecularBasis::slice, "Slices fraction of molecular basis for specific atoms.")
        .def("reduce_to_valence_basis", &CMolecularBasis::reduce_to_valence_basis, "Reduces molecular basis to it's valence only form.")
        .def("basis_sets", &CMolecularBasis::basis_sets, "Gets unique atomic basis sets in molecular basis")
        .def("basis_sets_indices", &CMolecularBasis::basis_sets_indices, "Gets vector of basis sets indices.")
        .def("max_angular_momentum",
             py::overload_cast<>(&CMolecularBasis::max_angular_momentum, py::const_),
             "Gets maximum angular momentum of molecular basis.")
        .def("max_angular_momentum",
             py::overload_cast<const std::vector<int> &>(&CMolecularBasis::max_angular_momentum, py::const_),
             "Gets maximum angular momentum of molecular basis for list of "
             "specific atoms.")
        .def("basis_functions", py::overload_cast<>(&CMolecularBasis::basis_functions, py::const_), "Gets vector of GTOs from molecular basis.")
        .def("basis_functions",
             py::overload_cast<const int>(&CMolecularBasis::basis_functions, py::const_),
             "Gets vector of GTOs with specific angular momentum from molecular "
             "basis.")
        .def("basis_functions",
             py::overload_cast<const int, const size_t>(&CMolecularBasis::basis_functions, py::const_),
             "Gets vector of GTOs with specific angular momentum and number of "
             "primitive GTOs from "
             "molecular basis.")
        .def("basis_functions",
             py::overload_cast<const std::vector<int> &>(&CMolecularBasis::basis_functions, py::const_),
             "Gets vector of GTOs from molecular basis for list of specific "
             "atoms.")
        .def("basis_functions",
             py::overload_cast<const std::vector<int> &, const int>(&CMolecularBasis::basis_functions, py::const_),
             "Gets vector of GTOs with specific angular momentum from molecular "
             "basis for list of specific atoms.")
        .def("basis_functions",
             py::overload_cast<const std::vector<int> &, const int, const size_t>(&CMolecularBasis::basis_functions, py::const_),
             "Gets vector of GTOs with specific angular momentum and number of "
             "primitive GTOs from molecular basis for list of specific atoms.")
        .def("atomic_indices",
             py::overload_cast<>(&CMolecularBasis::atomic_indices, py::const_),
             "Gets vector of atomic indices for GTOs from molecular basis.")
        .def("atomic_indices",
             py::overload_cast<const int>(&CMolecularBasis::atomic_indices, py::const_),
             "Gets vector of atomic indices for GTOs with specific angular "
             "momentum from molecular basis.")
        .def("atomic_indices",
             py::overload_cast<const int, const size_t>(&CMolecularBasis::atomic_indices, py::const_),
             "Gets vector of atomic indices for GTOs with specific angular "
             "momentum and number of primitive GTOs from molecular basis.")
        .def("atomic_indices",
             py::overload_cast<const std::vector<int> &>(&CMolecularBasis::atomic_indices, py::const_),
             "Gets vector of atomic indices for GTOs from molecular basis for "
             "list of specific "
             "atoms.")
        .def("atomic_indices",
             py::overload_cast<const std::vector<int> &, const int>(&CMolecularBasis::atomic_indices, py::const_),
             "Gets vector of atomic indices for GTOs with specific angular "
             "momentum from molecular basis for list of specific atoms.")
        .def("atomic_indices",
             py::overload_cast<const std::vector<int> &, const int, const size_t>(&CMolecularBasis::atomic_indices, py::const_),
             "Gets vector of atomic indices for GTOs with specific angular "
             "momentum and number of primitive GTOs from molecular basis for "
             "list of "
             "specific atoms.")
        .def("number_of_basis_functions",
             py::overload_cast<const int>(&CMolecularBasis::number_of_basis_functions, py::const_),
             "Gets number of GTOs with specific angular momentum from molecular "
             "basis.")
        .def("number_of_basis_functions",
             py::overload_cast<const int, const size_t>(&CMolecularBasis::number_of_basis_functions, py::const_),
             "Gets of GTOs with specific angular momentum and number of "
             "primitive GTOs from molecular basis.")
        .def("number_of_basis_functions",
             py::overload_cast<const std::vector<int> &, const int>(&CMolecularBasis::number_of_basis_functions, py::const_),
             "Gets number of GTOs with specific angular momentum from molecular "
             "basis for list of specific atoms.")
        .def("number_of_basis_functions",
             py::overload_cast<const std::vector<int> &, const int, const size_t>(&CMolecularBasis::number_of_basis_functions, py::const_),
             "Gets number of GTOs with specific angular momentum and number of "
             "primitive GTOs from molecular basis for list of specific atoms.")
        .def("number_of_primitive_basis_functions",
             py::overload_cast<const int>(&CMolecularBasis::number_of_primitive_functions, py::const_),
             "Gets number of primitive GTOs with specific angular momentum from "
             "molecular basis.")
        .def("number_of_primitive_basis_functions",
             py::overload_cast<const std::vector<int> &, const int>(&CMolecularBasis::number_of_primitive_functions, py::const_),
             "Gets number of primitive GTOs with specific angular momentum from "
             "molecular basis for list of specific atoms.")
        .def("contraction_depths",
             py::overload_cast<const int>(&CMolecularBasis::contraction_depths, py::const_),
             "Gets contraction depths of GTOs with specific angular momentum "
             "from molecular basis.")
        .def("contraction_depths",
             py::overload_cast<const std::vector<int> &, const int>(&CMolecularBasis::contraction_depths, py::const_),
             "Gets contraction depths of GTOs with specific angular momentum "
             "from molecular basis for list of specific atoms.")
        .def("get_dimensions_of_basis",
             py::overload_cast<>(&CMolecularBasis::dimensions_of_basis, py::const_),
             "Gets full dimensions of basis of molecular basis.")
        .def("get_dimensions_of_basis",
             py::overload_cast<const int>(&CMolecularBasis::dimensions_of_basis, py::const_),
             "Gets partial dimensions of basis of molecular basis up to specific "
             "angular momentum.")
        .def("get_dimensions_of_primitive_basis",
             py::overload_cast<>(&CMolecularBasis::dimensions_of_primitive_basis, py::const_),
             "Gets full dimensions of primitive basis of molecular basis.")
        .def("get_index_map",
             py::overload_cast<const int, const size_t>(&CMolecularBasis::index_map, py::const_),
             "Gets compressed global index map for basis sets of specific "
             "angular "
             "momentum and number of primitive GTOs from molecular basis.")
        .def("get_index_map",
             py::overload_cast<const std::vector<int> &, const int, const size_t>(&CMolecularBasis::index_map, py::const_),
             "Gets compressed global index map for basis sets of specific "
             "angular "
             "momentum and number of primitive GTOs from molecular basis for "
             "list of specific atoms.")
        .def("get_main_basis_label", &CMolecularBasis::main_basis_label, "Gets main basis set label in molecular basis.")
        .def("__eq__", [](const CMolecularBasis &self, const CMolecularBasis &other) { return self == other; })
        .def("__copy__", [](const CMolecularBasis &self) { return CMolecularBasis(self); })
        .def("__deepcopy__", [](const CMolecularBasis &self, py::dict) { return CMolecularBasis(self); });

    // CGtoBlock class
    PyClass<CGtoBlock>(m, "GtoBlock")
        .def(py::init<>())
        .def(py::init<const CGtoBlock &>())
        .def(py::init<const CMolecularBasis &, const CMolecule &, const int, const int>())
        .def(py::init<const CMolecularBasis &, const CMolecule &, const std::vector<int> &, const int, const int>())
        .def("coordinates", &CGtoBlock::coordinates, "Gets vector of basis function Cartesian coordinates.")
        .def("exponents", &CGtoBlock::exponents, "Gets vector of basis function exponents.")
        .def("normalization_factors", &CGtoBlock::normalization_factors, "Gets vector of basis function normalization factors.")
        .def("orbital_indices", &CGtoBlock::orbital_indices, "Gets vector of orbital indices of GTOs.")
        .def("atomic_indices", &CGtoBlock::atomic_indices, "Gets vector of atomic indices of GTOs.")
        .def("angular_momentum", &CGtoBlock::angular_momentum, "Gets angular momentum of GTOs block.")
        .def("number_of_primitives", &CGtoBlock::number_of_primitives, "Gets number of primitive GTOs in basis function.")
        .def("number_of_basis_functions", &CGtoBlock::number_of_basis_functions, "Gets number of GTOs in basis function block.")
        .def("__eq__", [](const CGtoBlock &self, const CGtoBlock &other) { return self == other; })
        .def("__copy__", [](const CGtoBlock &self) { return CGtoBlock(self); })
        .def("__deepcopy__", [](const CGtoBlock &self, py::dict) { return CGtoBlock(self); });
}

}  // namespace vlx_orbdata
