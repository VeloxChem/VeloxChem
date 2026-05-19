//
//  Tabula — custom-recursion molecular-integral machinery.
//  Python (pybind11) bindings.
//

#include "ExportTabula.hpp"

#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <utility>

#include "TabulaBlockSparseMatrix.hpp"
#include "TabulaDenseMatrix.hpp"
#include "TabulaMixedPrecisionBlockSparseMatrix.hpp"
#include "TabulaOverlapDriver.hpp"

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

    // tabula::BlockSparseMatrix::Block descriptor

    py::class_<BlockSparseMatrix::Block>(m, "TabulaBlockSparseBlock", "A block descriptor of a Tabula block-sparse matrix.")
        .def_readonly("group_a", &BlockSparseMatrix::Block::groupA)
        .def_readonly("group_b", &BlockSparseMatrix::Block::groupB)
        .def_readonly("row_count", &BlockSparseMatrix::Block::rowCount)
        .def_readonly("column_count", &BlockSparseMatrix::Block::columnCount)
        .def_readonly("offset", &BlockSparseMatrix::Block::offset);

    // tabula::BlockSparseMatrix class

    py::class_<BlockSparseMatrix, std::shared_ptr<BlockSparseMatrix>>(
        m, "TabulaBlockSparseMatrix", "A symmetric block-sparse Tabula matrix.")
        .def(py::init<>())
        .def(py::init<std::size_t,
                      const std::vector<std::vector<std::size_t>>&,
                      const std::vector<std::pair<std::size_t, std::size_t>>&>(),
             "dimension"_a,
             "group_global_ao"_a,
             "group_pairs"_a)
        .def("dimension", &BlockSparseMatrix::dimension, "Gets the full-matrix dimension.")
        .def("number_of_groups", &BlockSparseMatrix::number_of_groups, "Gets the number of AO groups.")
        .def("number_of_blocks", &BlockSparseMatrix::number_of_blocks, "Gets the number of stored blocks.")
        .def("block", &BlockSparseMatrix::block, "Gets a stored block descriptor.", "index"_a)
        .def("group_global_ao", &BlockSparseMatrix::group_global_ao, "Gets a group's global AO indices.", "group"_a)
        .def("value", &BlockSparseMatrix::value, "Gets a block-local element.", "block_index"_a, "row"_a, "column"_a)
        .def("set_value", &BlockSparseMatrix::set_value, "Sets a block-local element.", "block_index"_a, "row"_a, "column"_a, "value"_a)
        .def("stored_element_count", &BlockSparseMatrix::stored_element_count, "Gets the stored scalar count.")
        .def("to_dense", &BlockSparseMatrix::to_dense, "Reconstructs the dense symmetric matrix.");

    // tabula::MixedPrecisionBlockSparseMatrix::Block descriptor

    py::class_<MixedPrecisionBlockSparseMatrix::Block>(
        m, "TabulaMixedPrecisionBlock", "A block descriptor of a Tabula mixed-precision block-sparse matrix.")
        .def_readonly("group_a", &MixedPrecisionBlockSparseMatrix::Block::groupA)
        .def_readonly("group_b", &MixedPrecisionBlockSparseMatrix::Block::groupB)
        .def_readonly("row_count", &MixedPrecisionBlockSparseMatrix::Block::rowCount)
        .def_readonly("column_count", &MixedPrecisionBlockSparseMatrix::Block::columnCount)
        .def_readonly("offset", &MixedPrecisionBlockSparseMatrix::Block::offset)
        .def_readonly("is_single_precision", &MixedPrecisionBlockSparseMatrix::Block::isSinglePrecision);

    // tabula::MixedPrecisionBlockSparseMatrix class

    py::class_<MixedPrecisionBlockSparseMatrix, std::shared_ptr<MixedPrecisionBlockSparseMatrix>>(
        m, "TabulaMixedPrecisionBlockSparseMatrix", "A mixed-precision symmetric block-sparse Tabula matrix.")
        .def(py::init<>())
        .def(py::init<const BlockSparseMatrix&, double>(), "source"_a, "precision_threshold"_a)
        .def("dimension", &MixedPrecisionBlockSparseMatrix::dimension, "Gets the full-matrix dimension.")
        .def("number_of_blocks", &MixedPrecisionBlockSparseMatrix::number_of_blocks, "Gets the number of stored blocks.")
        .def("block", &MixedPrecisionBlockSparseMatrix::block, "Gets a stored block descriptor.", "index"_a)
        .def("group_global_ao", &MixedPrecisionBlockSparseMatrix::group_global_ao, "Gets a group's global AO indices.", "group"_a)
        .def("value", &MixedPrecisionBlockSparseMatrix::value, "Gets a block-local element.", "block_index"_a, "row"_a, "column"_a)
        .def("precision_threshold", &MixedPrecisionBlockSparseMatrix::precision_threshold, "Gets the precision threshold.")
        .def("single_block_count", &MixedPrecisionBlockSparseMatrix::single_block_count, "Gets the single-precision block count.")
        .def("double_block_count", &MixedPrecisionBlockSparseMatrix::double_block_count, "Gets the double-precision block count.")
        .def("stored_element_count", &MixedPrecisionBlockSparseMatrix::stored_element_count, "Gets the stored scalar count.")
        .def("stored_byte_count", &MixedPrecisionBlockSparseMatrix::stored_byte_count, "Gets the stored footprint in bytes.")
        .def("to_dense", &MixedPrecisionBlockSparseMatrix::to_dense, "Reconstructs the dense symmetric matrix.");

    // tabula::OverlapDriver class

    py::class_<OverlapDriver, std::shared_ptr<OverlapDriver>>(
        m, "TabulaOverlapDriver", "Driver for the Tabula two-center overlap integral.")
        .def(py::init<>())
        .def(
            "compute",
            [](const OverlapDriver& self, const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) {
                return self.compute(molecule, basis, threshold);
            },
            "Computes the overlap matrix.", "molecule"_a, "basis"_a, "threshold"_a = 0.0);

    // overlap driver — per-phase wall-time profile of a compute run

    m.def(
        "tabula_overlap_profile",
        [](const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) -> py::dict {
            OverlapProfile profile;
            OverlapDriver().compute(molecule, basis, threshold, &profile);

            py::dict result;
            result["make_blocks"] = profile.make_blocks;
            result["pair_setup"]  = profile.pair_setup;
            result["screen"]      = profile.screen;
            result["kernel"]      = profile.kernel;
            result["scatter"]     = profile.scatter;
            result["symmetrize"]  = profile.symmetrize;
            return result;
        },
        "Computes the overlap matrix and returns the per-phase wall-time breakdown.",
        "molecule"_a, "basis"_a, "threshold"_a = 0.0);

    // block-pair parallel loop — per-thread load balance of the last run

    m.def(
        "overlap_thread_balance",
        []() -> py::dict {
            const auto balance = overlap_thread_balance();

            py::dict result;
            result["wall"]  = balance.wall;
            result["busy"]  = balance.busy;
            result["pairs"] = balance.pairs;
            return result;
        },
        "Gets the per-thread load balance of the most recent overlap compute.");

}

}  // namespace tabula

PYBIND11_MODULE(tabulalib, m)
{
    tabula::export_tabula(m);
}
