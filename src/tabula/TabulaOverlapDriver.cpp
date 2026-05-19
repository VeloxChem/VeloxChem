//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center overlap driver.
//

#include "TabulaOverlapDriver.hpp"

#include <cstddef>
#include <utility>
#include <vector>

#include "GtoBlock.hpp"
#include "GtoFunc.hpp"
#include "GtoPairBlock.hpp"
#include "Point.hpp"
#include "TabulaContraction.hpp"
#include "TabulaMDRecursion.hpp"
#include "TabulaOverlapRecursion.hpp"
#include "TabulaOverlapScreener.hpp"
#include "TabulaOverlapTransform.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

/// @brief Evaluates a screened `CGtoPairBlock` into the overlap matrix — the
/// late-contraction recursion end to end: the seed ladder (a), the
/// primitive-pair contraction (b), the single-centre MD recursion (c), the
/// Cartesian-to-spherical assembly (d), and the scatter into the matrix (e).
///
/// An off-diagonal block pair also fills its transpose.
auto
evaluate_pair_block(DenseMatrix &matrix, const CGtoPairBlock &pair_block, const bool diagonal_pair) -> void
{
    const auto angular_momentums = pair_block.angular_momentums();
    const auto l_a               = angular_momentums.first;
    const auto l_c               = angular_momentums.second;
    const auto order             = static_cast<std::size_t>(l_a + l_c);

    const auto cdim = pair_block.number_of_contracted_pairs();
    if (cdim == 0) return;

    const auto nppairs = static_cast<std::size_t>(pair_block.number_of_primitive_pairs());

    // (a) + (b) — the contracted seed ladder [0]^m
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

    // (c) — the single-centre MD recursion [r]^0
    const auto rterms = compute_one_center_md(contracted, order, cdim, ac_x, ac_y, ac_z);

    // (d) — the Cartesian-to-spherical assembly
    const auto stride          = ((cdim + 7) / 8) * 8;
    const auto bra_components  = 2 * l_a + 1;
    const auto ket_components  = 2 * l_c + 1;

    std::vector<double> spherical(static_cast<std::size_t>(bra_components * ket_components) * stride, 0.0);
    overlap_transform(l_a, l_c, rterms.data(), cdim, spherical.data());

    // (e) — scatter the spherical block into the matrix; the orbital indices
    // carry the AO component stride ([0]) and the per-pair offset ([ij+1])
    const auto bra_orbitals = pair_block.bra_orbital_indices();
    const auto ket_orbitals = pair_block.ket_orbital_indices();

    for (int ca = 0; ca < bra_components; ca++)
    {
        for (int cc = 0; cc < ket_components; cc++)
        {
            const auto *row = spherical.data() + static_cast<std::size_t>(ca * ket_components + cc) * stride;

            for (std::size_t ij = 0; ij < cdim; ij++)
            {
                const auto r = static_cast<std::size_t>(ca) * bra_orbitals[0] + bra_orbitals[ij + 1];
                const auto c = static_cast<std::size_t>(cc) * ket_orbitals[0] + ket_orbitals[ij + 1];

                matrix(r, c) = row[ij];

                if (!diagonal_pair)
                {
                    matrix(c, r) = row[ij];
                }
            }
        }
    }
}

}  // namespace

auto
OverlapDriver::compute(const CMolecule &molecule, const CMolecularBasis &basis, const double threshold) const -> DenseMatrix
{
    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto dimension = static_cast<std::size_t>(gtofunc::getNumberOfAtomicOrbitals(gto_blocks));

    DenseMatrix matrix(dimension, dimension, Symmetry::symmetric);

    // the triangular basis-function-block pairs are the unit of parallel work
    // — each pair writes a disjoint AO region of the matrix, so the concurrent
    // writes do not race
    std::vector<std::pair<std::size_t, std::size_t>> block_pairs;
    for (std::size_t i = 0; i < gto_blocks.size(); i++)
    {
        for (std::size_t j = 0; j <= i; j++)
        {
            block_pairs.push_back({i, j});
        }
    }

#pragma omp parallel for schedule(dynamic)
    for (int p = 0; p < static_cast<int>(block_pairs.size()); p++)
    {
        const auto i = block_pairs[static_cast<std::size_t>(p)].first;
        const auto j = block_pairs[static_cast<std::size_t>(p)].second;

        const auto estimator = make_overlap_screening_estimator(gto_blocks[i].angular_momentum(),
                                                                gto_blocks[j].angular_momentum());

        // a screened pair block — the negligible contracted-GTO pairs are
        // dropped at construction
        const CGtoPairBlock pair_block(gto_blocks[i], gto_blocks[j], estimator, threshold);

        evaluate_pair_block(matrix, pair_block, i == j);
    }

    return matrix;
}

}  // namespace tabula
