//
//  Tabula — custom-recursion molecular-integral machinery.
//  Python (pybind11) bindings.
//

#include "ExportTabula.hpp"

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <utility>

#include "GtoBlock.hpp"
#include "TabulaBlockSparseMatrix.hpp"
#include "TabulaContraction.hpp"
#include "TabulaDenseMatrix.hpp"
#include "TabulaGtoPairBlock.hpp"
#include "TabulaMDRecursion.hpp"
#include "TabulaMixedPrecisionBlockSparseMatrix.hpp"
#include "TabulaOverlapDriver.hpp"
#include "TabulaOverlapRecursion.hpp"
#include "TabulaOverlapTransform.hpp"

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

    // tabula::GtoPairBlock class

    py::class_<GtoPairBlock>(m, "TabulaGtoPairBlock", "A Tabula-owned basis-function-pair block.")
        .def(py::init<const CGtoBlock&, const CGtoBlock&>(), "bra"_a, "ket"_a)
        .def(py::init<const CGtoBlock&, const CGtoBlock&, const ScreeningEstimator&, double>(),
             "bra"_a,
             "ket"_a,
             "estimator"_a,
             "threshold"_a)
        .def("angular_momentums", &GtoPairBlock::angular_momentums, "Gets the (l_a, l_c) angular momentums.")
        .def("number_of_contracted_pairs", &GtoPairBlock::number_of_contracted_pairs, "Gets the contracted-pair count.")
        .def("number_of_primitive_pairs", &GtoPairBlock::number_of_primitive_pairs, "Gets the primitive-pair count.")
        .def("bra_coordinates", &GtoPairBlock::bra_coordinates, "Gets the bra-center coordinates.")
        .def("ket_coordinates", &GtoPairBlock::ket_coordinates, "Gets the ket-center coordinates.")
        .def("bra_orbital_indices", &GtoPairBlock::bra_orbital_indices, "Gets the bra-side AO indices.")
        .def("ket_orbital_indices", &GtoPairBlock::ket_orbital_indices, "Gets the ket-side AO indices.")
        .def(
            "bra_exponents",
            [](const GtoPairBlock& self) -> py::array_t<double> {
                const auto pdim = self.number_of_contracted_pairs() * static_cast<std::size_t>(self.number_of_primitive_pairs());
                return py::array_t<double>({pdim}, {sizeof(double)}, self.bra_exponents());
            },
            "Gets the primitive-pair bra exponents.")
        .def(
            "ket_exponents",
            [](const GtoPairBlock& self) -> py::array_t<double> {
                const auto pdim = self.number_of_contracted_pairs() * static_cast<std::size_t>(self.number_of_primitive_pairs());
                return py::array_t<double>({pdim}, {sizeof(double)}, self.ket_exponents());
            },
            "Gets the primitive-pair ket exponents.")
        .def(
            "weights",
            [](const GtoPairBlock& self) -> py::array_t<double> {
                const auto pdim = self.number_of_contracted_pairs() * static_cast<std::size_t>(self.number_of_primitive_pairs());
                return py::array_t<double>({pdim}, {sizeof(double)}, self.weights());
            },
            "Gets the primitive-pair weights (normalization × overlap factor).");

    // overlap recursion — step (a): the seed ladder [0]^m

    m.def(
        "tabula_overlap_seed",
        [](const GtoPairBlock& pair_block) -> py::array_t<double> {
            const auto seed = compute_overlap_seed(pair_block);

            const auto [l_a, l_c] = pair_block.angular_momentums();
            const auto order      = l_a + l_c;

            const auto cdim    = pair_block.number_of_contracted_pairs();
            const auto nppairs = static_cast<std::size_t>(pair_block.number_of_primitive_pairs());
            const auto pdim    = cdim * nppairs;
            const auto stride  = ((pdim + 7) / 8) * 8;

            py::array_t<double> result({static_cast<std::size_t>(order + 1), pdim});
            auto                view = result.mutable_unchecked<2>();

            for (int m = 0; m <= order; m++)
            {
                for (std::size_t k = 0; k < pdim; k++)
                {
                    view(m, k) = seed[static_cast<std::size_t>(m) * stride + k];
                }
            }

            return result;
        },
        "Computes the overlap seed ladder [0]^m of a basis-function-pair block.",
        "pair_block"_a);

    // overlap recursion — step (b): the seed ladder contracted over primitives

    m.def(
        "tabula_overlap_contracted",
        [](const GtoPairBlock& pair_block) -> py::array_t<double> {
            const auto seed = compute_overlap_seed(pair_block);

            const auto angular_momentums = pair_block.angular_momentums();
            const auto rows              = static_cast<std::size_t>(
                angular_momentums.first + angular_momentums.second + 1);

            const auto cdim    = pair_block.number_of_contracted_pairs();
            const auto nppairs = static_cast<std::size_t>(pair_block.number_of_primitive_pairs());

            const auto contracted = contract_primitive_pairs(seed, rows, cdim, nppairs);

            const auto stride = ((cdim + 7) / 8) * 8;

            py::array_t<double> result({rows, cdim});
            auto                view = result.mutable_unchecked<2>();

            for (std::size_t m = 0; m < rows; m++)
            {
                for (std::size_t ij = 0; ij < cdim; ij++)
                {
                    view(m, ij) = contracted[m * stride + ij];
                }
            }

            return result;
        },
        "Computes the contracted overlap seed ladder [0]^m of a basis-function-pair block.",
        "pair_block"_a);

    // overlap recursion — step (c): the single-centre MD recursion [r]^0

    m.def(
        "tabula_overlap_rterms",
        [](const GtoPairBlock& pair_block) -> py::array_t<double> {
            const auto angular_momentums = pair_block.angular_momentums();
            const auto order             = static_cast<std::size_t>(
                angular_momentums.first + angular_momentums.second);

            const auto cdim    = pair_block.number_of_contracted_pairs();
            const auto nppairs = static_cast<std::size_t>(pair_block.number_of_primitive_pairs());

            // steps (a) + (b) — the contracted seed ladder
            const auto seed       = compute_overlap_seed(pair_block);
            const auto contracted = contract_primitive_pairs(seed, order + 1, cdim, nppairs);

            // AC = A − C, per contracted pair
            const auto bra_coords = pair_block.bra_coordinates();
            const auto ket_coords = pair_block.ket_coordinates();

            std::vector<double> ac_x(cdim), ac_y(cdim), ac_z(cdim);
            for (std::size_t ij = 0; ij < cdim; ij++)
            {
                const auto a = bra_coords[ij].coordinates();
                const auto c = ket_coords[ij].coordinates();
                ac_x[ij]     = a[0] - c[0];
                ac_y[ij]     = a[1] - c[1];
                ac_z[ij]     = a[2] - c[2];
            }

            // step (c) — the single-centre MD recursion
            const auto rterms = compute_one_center_md(contracted, order, cdim, ac_x, ac_y, ac_z);

            const auto monomials = (order + 1) * (order + 2) / 2;
            const auto stride    = ((cdim + 7) / 8) * 8;

            py::array_t<double> result({monomials, cdim});
            auto                view = result.mutable_unchecked<2>();

            for (std::size_t r = 0; r < monomials; r++)
            {
                for (std::size_t ij = 0; ij < cdim; ij++)
                {
                    view(r, ij) = rterms[r * stride + ij];
                }
            }

            return result;
        },
        "Computes the single-centre MD recursion terms [r]^0, |r| = l_a+l_c, of a pair block.",
        "pair_block"_a);

    // overlap recursion — step (d): the assembled spherical overlap block

    m.def(
        "tabula_overlap_spherical",
        [](const GtoPairBlock& pair_block) -> py::array_t<double> {
            const auto angular_momentums = pair_block.angular_momentums();
            const auto l_a               = angular_momentums.first;
            const auto l_c               = angular_momentums.second;
            const auto order             = static_cast<std::size_t>(l_a + l_c);

            const auto cdim    = pair_block.number_of_contracted_pairs();
            const auto nppairs = static_cast<std::size_t>(pair_block.number_of_primitive_pairs());

            // steps (a) + (b) — the contracted seed ladder
            const auto seed       = compute_overlap_seed(pair_block);
            const auto contracted = contract_primitive_pairs(seed, order + 1, cdim, nppairs);

            // AC = A − C, per contracted pair
            const auto bra_coords = pair_block.bra_coordinates();
            const auto ket_coords = pair_block.ket_coordinates();

            std::vector<double> ac_x(cdim), ac_y(cdim), ac_z(cdim);
            for (std::size_t ij = 0; ij < cdim; ij++)
            {
                const auto a = bra_coords[ij].coordinates();
                const auto c = ket_coords[ij].coordinates();
                ac_x[ij]     = a[0] - c[0];
                ac_y[ij]     = a[1] - c[1];
                ac_z[ij]     = a[2] - c[2];
            }

            // step (c) — the single-centre MD recursion
            const auto rterms = compute_one_center_md(contracted, order, cdim, ac_x, ac_y, ac_z);

            // step (d) — the Cartesian-to-spherical assembly
            const auto stride         = ((cdim + 7) / 8) * 8;
            const auto component_rows = static_cast<std::size_t>((2 * l_a + 1) * (2 * l_c + 1));

            std::vector<double> spherical(component_rows * stride, 0.0);
            overlap_transform(l_a, l_c, rterms.data(), cdim, spherical.data());

            py::array_t<double> result({component_rows, cdim});
            auto                view = result.mutable_unchecked<2>();

            for (std::size_t s = 0; s < component_rows; s++)
            {
                for (std::size_t ij = 0; ij < cdim; ij++)
                {
                    view(s, ij) = spherical[s * stride + ij];
                }
            }

            return result;
        },
        "Computes the assembled spherical overlap block of a basis-function-pair block.",
        "pair_block"_a);

    // overlap driver — per-phase wall-time profile of a compute run

    m.def(
        "tabula_overlap_profile",
        [](const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) -> py::dict {
            OverlapProfile profile;
            OverlapDriver().compute(molecule, basis, threshold, &profile);

            py::dict result;
            result["make_blocks"] = profile.make_blocks;
            result["pair_setup"]  = profile.pair_setup;
            result["seed"]        = profile.seed;
            result["contract"]    = profile.contract;
            result["md"]          = profile.md;
            result["transform"]   = profile.transform;
            result["scatter"]     = profile.scatter;
            result["symmetrize"]  = profile.symmetrize;
            return result;
        },
        "Computes the overlap matrix and returns the per-phase wall-time breakdown.",
        "molecule"_a, "basis"_a, "threshold"_a = 0.0);

    // overlap seed — accumulated allocate / row0 / ladder profile

    m.def(
        "seed_profile",
        []() -> py::dict {
            const auto profile = seed_profile();

            py::dict result;
            result["allocate"] = profile.allocate;
            result["row0"]     = profile.row0;
            result["ladder"]   = profile.ladder;
            return result;
        },
        "Gets the accumulated compute_overlap_seed allocate / row0 / ladder profile.");

    m.def("reset_seed_profile", &reset_seed_profile, "Resets the accumulated compute_overlap_seed profile.");

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
