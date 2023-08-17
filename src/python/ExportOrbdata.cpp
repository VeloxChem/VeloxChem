#include "ExportOrbdata.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "AtomBasis.hpp"
#include "BasisFunction.hpp"
#include "GtoBlock.hpp"
#include "GtoPairBlock.hpp"
#include "GtoFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

namespace vlx_orbdata {  // vlx_orbdata namespace

// Exports classes/functions in src/orbdata to python

void
export_orbdata(py::module &m)
{
    // exposing functions from GtoFunc.hpp

    m.def("make_gto_blocks",
          py::overload_cast<const CMolecularBasis &, const CMolecule &>(&gtofunc::makeGtoBlocks),
          "Creates vector of GtoBlocks for given basis and molecule.");
    m.def("make_gto_blocks",
          py::overload_cast<const CMolecularBasis &, const CMolecule &, const std::vector<int64_t> &>(&gtofunc::makeGtoBlocks),
          "Creates vector of GtoBlocks for selected atoms in given basis and molecule.");

    // CBasisFunction class

    PyClass<CBasisFunction>(m, "BasisFunction")
        .def(py::init<>())
        .def(py::init<const CBasisFunction &>())
        .def(py::init<const std::vector<double> &, const std::vector<double> &, const int64_t>())
        .def("set_exponents", &CBasisFunction::setExponents, "Sets exponents of basis function.")
        .def("set_normalization_factors", &CBasisFunction::setNormalizationFactors, "Gets name of chemical element.")
        .def("set_angular_momentum", &CBasisFunction::setAngularMomentum, "Sets angular momentum of basis function.")
        .def("add", &CBasisFunction::add, "Add primitive i.e. exponent and normalization factor to basis function.")
        .def("normalize", &CBasisFunction::normalize, "Normalizes primitive GTOs in basis function.")
        .def("normalize", &CBasisFunction::normalize, "Normalizes primitive GTOs in basis function.")
        .def("get_exponents", &CBasisFunction::getExponents, "Gets vector of exponents in basis function.")
        .def("get_normalization_factors", &CBasisFunction::getNormalizationFactors, "Gets vector of normalization factors in basis function.")
        .def("get_angular_momentum", &CBasisFunction::getAngularMomentum, "Gets angular momentum of basis function.")
        .def("number_of_primitives", &CBasisFunction::getNumberOfPrimitiveFunctions, "Gets number of primitives in basis function.");

    // CAtomBasis class

    PyClass<CAtomBasis>(m, "AtomBasis")
        .def(py::init<>())
        .def(py::init<const CAtomBasis &>())
        .def(py::init<const std::vector<CBasisFunction> &, const int64_t, const std::string &, const std::string &>())
        .def("set_identifier", &CAtomBasis::setIdentifier, "Sets identifier of atom basis.")
        .def("set_name", &CAtomBasis::setName, "Sets name of atom basis.")
        .def("set_ecp_label", &CAtomBasis::setEffectiveCorePotentialLabel, "Sets effective core potential label of atom basis.")
        .def("add", &CAtomBasis::add, "Adds basis function to atom basis.")
        .def("reduce_to_valence_basis", &CAtomBasis::reduceToValenceBasis, "Reduces atom basis to it's valence only form.")
        .def("get_basis_functions", py::overload_cast<>(&CAtomBasis::getBasisFunctions, py::const_), "Gets GTOs.")
        .def("get_basis_functions",
             py::overload_cast<const int64_t>(&CAtomBasis::getBasisFunctions, py::const_),
             "Gets GTOs with specific angular momentum.")
        .def("get_basis_functions",
             py::overload_cast<const int64_t, const int64_t>(&CAtomBasis::getBasisFunctions, py::const_),
             "Gets GTOs with specific angular momentum and number of primitives.")
        .def("get_identifier", &CAtomBasis::getIdentifier, "Gets identifier of atom basis.")
        .def("get_name", &CAtomBasis::getName, "Gets name of atom basis.")
        .def("get_ecp_label", &CAtomBasis::getEffectiveCorePotentialLabel, "Gets effective core potential label of atom basis.")
        .def("need_ecp", &CAtomBasis::needEffectiveCorePotential, "Checks if atom basis requires effective core potential.")
        .def("max_angular_momentum", &CAtomBasis::getMaxAngularMomentum, "Gets maximum angular momentum in atom basis.")
        .def("number_of_basis_functions",
             py::overload_cast<const int64_t>(&CAtomBasis::getNumberOfBasisFunctions, py::const_),
             "Gets number of GTOs with specific angular momentum.")
        .def("number_of_basis_functions",
             py::overload_cast<const int64_t, const int64_t>(&CAtomBasis::getNumberOfBasisFunctions, py::const_),
             "Gets number of GTOs with specific angular momentum and number of primitives.")
        .def("number_of_primitive_basis_functions",
             &CAtomBasis::getNumberOfPrimitiveFunctions,
             "Gets number of primitive GTOs with specific angular momentum.")
        .def("contraction_depths", &CAtomBasis::getContractionDepths, "Gets contraction depths of GTOs with specific angular momentum.")
        .def("contraction_str", &CAtomBasis::getContractionString, "Gets contraction string of atom basis.")
        .def("primitives_str", &CAtomBasis::getPrimitivesString, "Gets primitive GTOs string of atom basis.");

    // CMolecularBasis class

    PyClass<CMolecularBasis>(m, "MolecularBasis")
        .def(py::init<>())
        .def(py::init<const CMolecularBasis &>())
        .def(py::init<const std::vector<CAtomBasis> &, const std::vector<int64_t> &>())
        .def("add", &CMolecularBasis::add, "Adds atomic basis to molecular basis.")
        .def("slice", &CMolecularBasis::slice, "Slices fraction of molecular basis for specific atoms.")
        .def("reduce_to_valence_basis", &CMolecularBasis::reduceToValenceBasis, "Reduces molecular basis to it's valence only form.")
        .def("get_basis_sets", &CMolecularBasis::getBasisSets, "Gets unique atomic basis sets in molecular basis")
        .def("get_basis_sets_indexes", &CMolecularBasis::getBasisSetsIndexes, "Gets vector of basis sets indexes.")
        .def("max_angular_momentum",
             py::overload_cast<>(&CMolecularBasis::getMaxAngularMomentum, py::const_),
             "Gets maximum angular momentum of molecular basis.")
        .def("max_angular_momentum",
             py::overload_cast<const std::vector<int64_t> &>(&CMolecularBasis::getMaxAngularMomentum, py::const_),
             "Gets maximum angular momentum of molecular basis for list of specific atoms.")
        .def("get_basis_functions", py::overload_cast<>(&CMolecularBasis::getBasisFunctions, py::const_), "Gets vector of GTOs from molecular basis.")
        .def("get_basis_functions",
             py::overload_cast<const int64_t>(&CMolecularBasis::getBasisFunctions, py::const_),
             "Gets vector of GTOs with specific angular momentum from molecular basis.")
        .def("get_basis_functions",
             py::overload_cast<const int64_t, const int64_t>(&CMolecularBasis::getBasisFunctions, py::const_),
             "Gets vector of GTOs with specific angular momentum and number of primitive GTOs from molecular basis.")
        .def("get_basis_functions",
             py::overload_cast<const std::vector<int64_t> &>(&CMolecularBasis::getBasisFunctions, py::const_),
             "Gets vector of GTOs from molecular basis for list of specific atoms.")
        .def("get_basis_functions",
             py::overload_cast<const std::vector<int64_t> &, const int64_t>(&CMolecularBasis::getBasisFunctions, py::const_),
             "Gets vector of GTOs with specific angular momentum from molecular basis for list of specific atoms.")
        .def("get_basis_functions",
             py::overload_cast<const std::vector<int64_t> &, const int64_t, const int64_t>(&CMolecularBasis::getBasisFunctions, py::const_),
             "Gets vector of GTOs with specific angular momentum and number of primitive GTOs from molecular basis for list of specific atoms.")
        .def("get_atomic_indexes",
             py::overload_cast<>(&CMolecularBasis::getAtomicIndexes, py::const_),
             "Gets vector of atomic indexes for GTOs from molecular basis.")
        .def("get_atomic_indexes",
             py::overload_cast<const int64_t>(&CMolecularBasis::getAtomicIndexes, py::const_),
             "Gets vector of atomic indexes for GTOs with specific angular momentum from molecular basis.")
        .def("get_atomic_indexes",
             py::overload_cast<const int64_t, const int64_t>(&CMolecularBasis::getAtomicIndexes, py::const_),
             "Gets vector of atomic indexes for GTOs with specific angular momentum and number of primitive GTOs from molecular basis.")
        .def("get_atomic_indexes",
             py::overload_cast<const std::vector<int64_t> &>(&CMolecularBasis::getAtomicIndexes, py::const_),
             "Gets vector of atomic indexes for GTOs from molecular basis for list of specific atoms.")
        .def("get_atomic_indexes",
             py::overload_cast<const std::vector<int64_t> &, const int64_t>(&CMolecularBasis::getAtomicIndexes, py::const_),
             "Gets vector of atomic indexes for GTOs with specific angular momentum from molecular basis for list of specific atoms.")
        .def("get_atomic_indexes",
             py::overload_cast<const std::vector<int64_t> &, const int64_t, const int64_t>(&CMolecularBasis::getAtomicIndexes, py::const_),
             "Gets vector of atomic indexes for GTOs with specific angular momentum and number of primitive GTOs from molecular basis for list of "
             "specific atoms.")
        .def("number_of_basis_functions",
             py::overload_cast<const int64_t>(&CMolecularBasis::getNumberOfBasisFunctions, py::const_),
             "Gets number of GTOs with specific angular momentum from molecular basis.")
        .def("number_of_basis_functions",
             py::overload_cast<const int64_t, const int64_t>(&CMolecularBasis::getNumberOfBasisFunctions, py::const_),
             "Gets of GTOs with specific angular momentum and number of primitive GTOs from molecular basis.")
        .def("number_of_basis_functions",
             py::overload_cast<const std::vector<int64_t> &, const int64_t>(&CMolecularBasis::getNumberOfBasisFunctions, py::const_),
             "Gets number of GTOs with specific angular momentum from molecular basis for list of specific atoms.")
        .def("number_of_basis_functions",
             py::overload_cast<const std::vector<int64_t> &, const int64_t, const int64_t>(&CMolecularBasis::getNumberOfBasisFunctions, py::const_),
             "Gets number of GTOs with specific angular momentum and number of primitive GTOs from molecular basis for list of specific atoms.")
        .def("number_of_primitive_basis_functions",
             py::overload_cast<const int64_t>(&CMolecularBasis::getNumberOfPrimitiveFunctions, py::const_),
             "Gets number of primitive GTOs with specific angular momentum from molecular basis.")
        .def("number_of_primitive_basis_functions",
             py::overload_cast<const std::vector<int64_t> &, const int64_t>(&CMolecularBasis::getNumberOfPrimitiveFunctions, py::const_),
             "Gets number of primitive GTOs with specific angular momentum from molecular basis for list of specific atoms.")
        .def("contraction_depths",
             py::overload_cast<const int64_t>(&CMolecularBasis::getContractionDepths, py::const_),
             "Gets contraction depths of GTOs with specific angular momentum from molecular basis.")
        .def("contraction_depths",
             py::overload_cast<const std::vector<int64_t> &, const int64_t>(&CMolecularBasis::getContractionDepths, py::const_),
             "Gets contraction depths of GTOs with specific angular momentum from molecular basis for list of specific atoms.")
        .def("get_dimensions_of_basis",
             py::overload_cast<>(&CMolecularBasis::getDimensionsOfBasis, py::const_),
             "Gets full dimensions of basis of molecular basis.")
        .def("get_dimensions_of_basis",
             py::overload_cast<const int64_t>(&CMolecularBasis::getDimensionsOfBasis, py::const_),
             "Gets partial dimensions of basis of molecular basis up to specific angular momentum.")
        .def("get_dimensions_of_primitive_basis",
             py::overload_cast<>(&CMolecularBasis::getDimensionsOfPrimitiveBasis, py::const_),
             "Gets full dimensions of primitive basis of molecular basis.")
        .def("get_index_map",
             py::overload_cast<const int64_t, const int64_t>(&CMolecularBasis::getIndexMap, py::const_),
             "Gets compressed global index map for basis sets of specific angular momentum and number of primitive GTOs from molecular basis.")
        .def("get_index_map",
             py::overload_cast<const std::vector<int64_t> &, const int64_t, const int64_t>(&CMolecularBasis::getIndexMap, py::const_),
             "Gets compressed global index map for basis sets of specific angular momentum and number of primitive GTOs from molecular basis for "
             "list of specific atoms.")
        .def("info_str",
             py::overload_cast<const std::string &>(&CMolecularBasis::printBasis, py::const_),
             "Generates basis set information string with specific title for molecular basis.")
        .def(
            "info_str", py::overload_cast<>(&CMolecularBasis::printBasis, py::const_), "Generates basis set information string for molecular basis.");

    // CGtoBlock class

    PyClass<CGtoBlock>(m, "GtoBlock")
        .def(py::init<>())
        .def(py::init<const CGtoBlock &>())
        .def(py::init<const std::vector<TPoint3D> &,
                      const std::vector<double> &,
                      const std::vector<double> &,
                      const std::vector<int64_t> &,
                      const std::vector<int64_t> &,
                      const int64_t,
                      const int64_t>())
        .def(py::init<const CMolecularBasis &, const CMolecule &, const int64_t, const int64_t>())
        .def(py::init<const CMolecularBasis &, const CMolecule &, const std::vector<int64_t> &, const int64_t, const int64_t>())
        .def("get_coordinates", &CGtoBlock::getCoordinates, "Gets vector of basis function coordinates.")
        .def("get_exponents", &CGtoBlock::getExponents, "Gets vector of basis function exponents.")
        .def("get_normalization_factors", &CGtoBlock::getNormalizationFactors, "Gets vector of basis function normalization factors.")
        .def("get_orbital_indexes", &CGtoBlock::getOrbitalIndexes, "Gets vector of orbital indexes of GTOs.")
        .def("get_atomic_indexes", &CGtoBlock::getAtomicIndexes, "Gets vector of atomic indexes of GTOs.")
        .def("get_angular_momentum", &CGtoBlock::getAngularMomentum, "Gets angular momentum of GTOs block.")
        .def("number_of_primitives", &CGtoBlock::getNumberOfPrimitives, "Gets number of primitive GTOs in basis function.")
        .def("number_of_basis_functions", &CGtoBlock::getNumberOfBasisFunctions, "Gets number of GTOs in basis function block.");
    
    // CGtoPairBlock class

    PyClass<CGtoPairBlock>(m, "GtoPairBlock")
        .def(py::init<>())
        .def(py::init<const CGtoPairBlock&>())
        .def(py::init<const std::vector<TPairOfPoints3D>&,
                      const std::vector<TPoint2D>&,
                      const std::vector<TPoint2D>&,
                      const std::vector<T2Index>&,
                      const std::vector<T2Index>&,
                      const T2Index&,
                      const int64_t>())
        .def(py::init<const CGtoBlock&>())
        .def(py::init<const CGtoBlock&, const CGtoBlock&>())
        .def("get_coordinates", &CGtoPairBlock::getCoordinates, "Gets vector of atoms pair coordinates.")
        .def("get_exponents", &CGtoPairBlock::getExponents, "Gets vector of basis function exponents pairs.")
        .def("get_normalization_factors", &CGtoPairBlock::getNormalizationFactors, "Gets vector of basis function normalization factors pairs.")
        .def("get_orbital_indexes", &CGtoPairBlock::getOrbitalIndexes, "Gets vector of orbital index pairs.")
        .def("get_atomic_indexes", &CGtoPairBlock::getAtomicIndexes, "Gets vector of atomic index pairs.")
        .def("get_angular_momentum", &CGtoPairBlock::getAngularMomentums, "Gets angular momentum of GTOs pair.")
        .def("number_of_primitives", &CGtoPairBlock::getNumberOfPrimitivePairs, "Gets number of primitive GTO pairs in GTO pair.")
        .def("number_of_basis_functions", &CGtoPairBlock::getNumberOfContractedPairs, "Gets number of contracted GTO pairs in basis function pairs block.");

    // ...
}

}  // namespace vlx_orbdata
