//
//  Tabula — custom-recursion molecular-integral machinery.
//  Python (pybind11) bindings.
//

#include "ExportTabula.hpp"

#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include "TabulaBlockSparseMatrix.hpp"
#include "TabulaBoys.hpp"
#include "TabulaChargeDipoleDriver.hpp"
#include "TabulaCoulombDriver.hpp"
#include "TabulaDenseMatrix.hpp"
#include "TabulaKineticDriver.hpp"
#include "TabulaMixedPrecisionBlockSparseMatrix.hpp"
#include "TabulaNuclearAttractionDriver.hpp"
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
            "Computes the overlap matrix.", "molecule"_a, "basis"_a, "threshold"_a = 0.0)
        .def(
            "compute_sparse",
            [](const OverlapDriver& self, const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) {
                return self.computeSparse(molecule, basis, threshold);
            },
            "Computes the overlap matrix in block-sparse storage.", "molecule"_a, "basis"_a, "threshold"_a = 0.0);

    // overlap driver — per-phase wall-time profile of a compute run

    m.def(
        "tabula_overlap_profile",
        [](const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) -> py::dict {
            KernelProfile profile;
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

    m.def(
        "tabula_overlap_profile_sparse",
        [](const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) -> py::dict {
            KernelProfile profile;
            OverlapDriver().computeSparse(molecule, basis, threshold, &profile);

            py::dict result;
            result["make_blocks"] = profile.make_blocks;
            result["pair_setup"]  = profile.pair_setup;
            result["screen"]      = profile.screen;
            result["kernel"]      = profile.kernel;
            result["scatter"]     = profile.scatter;
            result["symmetrize"]  = profile.symmetrize;
            return result;
        },
        "Computes the block-sparse overlap matrix and returns the per-phase wall-time breakdown.",
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

    // tabula::KineticDriver class

    py::class_<KineticDriver, std::shared_ptr<KineticDriver>>(
        m, "TabulaKineticDriver", "Driver for the Tabula two-center kinetic-energy integral.")
        .def(py::init<>())
        .def(
            "compute",
            [](const KineticDriver& self, const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) {
                return self.compute(molecule, basis, threshold);
            },
            "Computes the kinetic-energy matrix.", "molecule"_a, "basis"_a, "threshold"_a = 0.0)
        .def(
            "compute_sparse",
            [](const KineticDriver& self, const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) {
                return self.computeSparse(molecule, basis, threshold);
            },
            "Computes the kinetic-energy matrix in block-sparse storage.", "molecule"_a, "basis"_a, "threshold"_a = 0.0);

    // kinetic driver — per-phase wall-time profile of a compute run

    m.def(
        "tabula_kinetic_profile",
        [](const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) -> py::dict {
            KernelProfile profile;
            KineticDriver().compute(molecule, basis, threshold, &profile);

            py::dict result;
            result["make_blocks"] = profile.make_blocks;
            result["pair_setup"]  = profile.pair_setup;
            result["screen"]      = profile.screen;
            result["kernel"]      = profile.kernel;
            result["scatter"]     = profile.scatter;
            result["symmetrize"]  = profile.symmetrize;
            return result;
        },
        "Computes the kinetic-energy matrix and returns the per-phase wall-time breakdown.",
        "molecule"_a, "basis"_a, "threshold"_a = 0.0);

    m.def(
        "tabula_kinetic_profile_sparse",
        [](const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) -> py::dict {
            KernelProfile profile;
            KineticDriver().computeSparse(molecule, basis, threshold, &profile);

            py::dict result;
            result["make_blocks"] = profile.make_blocks;
            result["pair_setup"]  = profile.pair_setup;
            result["screen"]      = profile.screen;
            result["kernel"]      = profile.kernel;
            result["scatter"]     = profile.scatter;
            result["symmetrize"]  = profile.symmetrize;
            return result;
        },
        "Computes the block-sparse kinetic-energy matrix and returns the per-phase wall-time breakdown.",
        "molecule"_a, "basis"_a, "threshold"_a = 0.0);

    m.def(
        "kinetic_thread_balance",
        []() -> py::dict {
            const auto balance = kinetic_thread_balance();

            py::dict result;
            result["wall"]  = balance.wall;
            result["busy"]  = balance.busy;
            result["pairs"] = balance.pairs;
            return result;
        },
        "Gets the per-thread load balance of the most recent kinetic compute.");

    // tabula::CoulombDriver class

    py::class_<CoulombDriver, std::shared_ptr<CoulombDriver>>(
        m, "TabulaCoulombDriver", "Driver for the Tabula two-center Coulomb integral.")
        .def(py::init<>())
        .def(
            "compute",
            [](const CoulombDriver& self, const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) {
                return self.compute(molecule, basis, threshold);
            },
            "Computes the Coulomb matrix.", "molecule"_a, "basis"_a, "threshold"_a = 0.0)
        .def(
            "compute_sparse",
            [](const CoulombDriver& self, const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) {
                return self.computeSparse(molecule, basis, threshold);
            },
            "Computes the Coulomb matrix in block-sparse storage.", "molecule"_a, "basis"_a, "threshold"_a = 0.0);

    // Coulomb driver — per-phase wall-time profile of a compute run

    m.def(
        "tabula_coulomb_profile",
        [](const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) -> py::dict {
            KernelProfile profile;
            CoulombDriver().compute(molecule, basis, threshold, &profile);

            py::dict result;
            result["make_blocks"] = profile.make_blocks;
            result["pair_setup"]  = profile.pair_setup;
            result["screen"]      = profile.screen;
            result["kernel"]      = profile.kernel;
            result["scatter"]     = profile.scatter;
            result["symmetrize"]  = profile.symmetrize;
            return result;
        },
        "Computes the Coulomb matrix and returns the per-phase wall-time breakdown.",
        "molecule"_a, "basis"_a, "threshold"_a = 0.0);

    m.def(
        "tabula_coulomb_profile_sparse",
        [](const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) -> py::dict {
            KernelProfile profile;
            CoulombDriver().computeSparse(molecule, basis, threshold, &profile);

            py::dict result;
            result["make_blocks"] = profile.make_blocks;
            result["pair_setup"]  = profile.pair_setup;
            result["screen"]      = profile.screen;
            result["kernel"]      = profile.kernel;
            result["scatter"]     = profile.scatter;
            result["symmetrize"]  = profile.symmetrize;
            return result;
        },
        "Computes the block-sparse Coulomb matrix and returns the per-phase wall-time breakdown.",
        "molecule"_a, "basis"_a, "threshold"_a = 0.0);

    m.def(
        "coulomb_thread_balance",
        []() -> py::dict {
            const auto balance = coulomb_thread_balance();

            py::dict result;
            result["wall"]  = balance.wall;
            result["busy"]  = balance.busy;
            result["pairs"] = balance.pairs;
            return result;
        },
        "Gets the per-thread load balance of the most recent Coulomb compute.");

    // tabula::NuclearAttractionDriver class

    py::class_<NuclearAttractionDriver, std::shared_ptr<NuclearAttractionDriver>>(
        m, "TabulaNuclearAttractionDriver", "Driver for the Tabula two-center nuclear-attraction integral.")
        .def(py::init<>())
        .def(
            "compute",
            [](const NuclearAttractionDriver& self, const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) {
                return self.compute(molecule, basis, threshold);
            },
            "Computes the nuclear-attraction matrix over the molecule's nuclei. "
            "threshold < 0 (default) auto-screens large molecules; 0 is exact dense.",
            "molecule"_a, "basis"_a, "threshold"_a = -1.0)
        .def(
            "compute_sparse",
            [](const NuclearAttractionDriver& self, const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) {
                return self.computeSparse(molecule, basis, threshold);
            },
            "Computes the nuclear-attraction matrix in block-sparse storage, over the nuclei.",
            "molecule"_a, "basis"_a, "threshold"_a = 0.0)
        .def(
            "compute_external",
            [](const NuclearAttractionDriver&            self,
               const CMolecule&                          molecule,
               const CMolecularBasis&                    basis,
               const std::vector<double>&                magnitudes,
               const std::vector<std::array<double, 3>>& coordinates,
               const double                              threshold) {
                return self.compute(molecule, basis, magnitudes, coordinates, threshold);
            },
            "Computes the matrix over external point charges (coordinates in au). "
            "threshold < 0 (default) auto-screens large molecules; 0 is exact dense.",
            "molecule"_a, "basis"_a, "magnitudes"_a, "coordinates"_a, "threshold"_a = -1.0)
        .def(
            "compute_external_sparse",
            [](const NuclearAttractionDriver&            self,
               const CMolecule&                          molecule,
               const CMolecularBasis&                    basis,
               const std::vector<double>&                magnitudes,
               const std::vector<std::array<double, 3>>& coordinates,
               const double                              threshold) {
                return self.computeSparse(molecule, basis, magnitudes, coordinates, threshold);
            },
            "Block-sparse variant over external point charges.",
            "molecule"_a, "basis"_a, "magnitudes"_a, "coordinates"_a, "threshold"_a = 0.0);

    m.def(
        "tabula_nuclear_attraction_profile",
        [](const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) -> py::dict {
            KernelProfile profile;
            NuclearAttractionDriver().compute(molecule, basis, threshold, &profile);

            py::dict result;
            result["make_blocks"] = profile.make_blocks;
            result["pair_setup"]  = profile.pair_setup;
            result["screen"]      = profile.screen;
            result["kernel"]      = profile.kernel;
            result["scatter"]     = profile.scatter;
            result["symmetrize"]  = profile.symmetrize;
            return result;
        },
        "Computes the nuclear-attraction matrix and returns the per-phase wall-time breakdown.",
        "molecule"_a, "basis"_a, "threshold"_a = 0.0);

    m.def(
        "nuclear_attraction_thread_balance",
        []() -> py::dict {
            const auto balance = nuclear_attraction_thread_balance();

            py::dict result;
            result["wall"]  = balance.wall;
            result["busy"]  = balance.busy;
            result["pairs"] = balance.pairs;
            return result;
        },
        "Gets the per-thread load balance of the most recent nuclear-attraction compute.");

    // tabula::ChargeDipoleDriver class

    py::class_<ChargeDipoleDriver, std::shared_ptr<ChargeDipoleDriver>>(
        m, "TabulaChargeDipoleDriver", "Driver for the Tabula two-center charge-dipole integral.")
        .def(py::init<>())
        .def(
            "compute",
            [](const ChargeDipoleDriver&                 self,
               const CMolecule&                          molecule,
               const CMolecularBasis&                    basis,
               const std::vector<std::array<double, 3>>& moments,
               const std::vector<std::array<double, 3>>& coordinates,
               const double                              threshold) {
                return self.compute(molecule, basis, moments, coordinates, threshold);
            },
            "Computes the charge-dipole matrix Sum_N d_N . (a|(r-N)/|r-N|^3|c) over point dipoles "
            "(moments and coordinates in au). threshold < 0 (default) auto-screens large molecules; 0 is exact dense.",
            "molecule"_a, "basis"_a, "moments"_a, "coordinates"_a, "threshold"_a = -1.0)
        .def(
            "compute_sparse",
            [](const ChargeDipoleDriver&                 self,
               const CMolecule&                          molecule,
               const CMolecularBasis&                    basis,
               const std::vector<std::array<double, 3>>& moments,
               const std::vector<std::array<double, 3>>& coordinates,
               const double                              threshold) {
                return self.computeSparse(molecule, basis, moments, coordinates, threshold);
            },
            "Computes the charge-dipole matrix in block-sparse storage, over point dipoles.",
            "molecule"_a, "basis"_a, "moments"_a, "coordinates"_a, "threshold"_a = 0.0);

    m.def(
        "charge_dipole_thread_balance",
        []() -> py::dict {
            const auto balance = charge_dipole_thread_balance();

            py::dict result;
            result["wall"]  = balance.wall;
            result["busy"]  = balance.busy;
            result["pairs"] = balance.pairs;
            return result;
        },
        "Gets the per-thread load balance of the most recent charge-dipole compute.");

    // tabula::boys — the Coulomb / nuclear-attraction integral seed

    m.def(
        "tabula_boys",
        [](const int order, const double x) -> std::vector<double> {
            std::vector<double> results(static_cast<std::size_t>(order) + 1, 0.0);
            boys(order, x, results.data());
            return results;
        },
        "Evaluates the Boys function F_0 … F_order at x.",
        "order"_a, "x"_a);

}

}  // namespace tabula

PYBIND11_MODULE(tabulalib, m)
{
    tabula::export_tabula(m);
}
