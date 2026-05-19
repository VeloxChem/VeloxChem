//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center overlap driver.
//

#include "TabulaOverlapDriver.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "GtoBlock.hpp"
#include "GtoFunc.hpp"
#include "GtoPairBlock.hpp"
#include "Point.hpp"
#include "TabulaOverlapScreener.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

// ====================================================================
//  RECURSION SEAM
//
//  This is where Tabula's custom recursion plugs in. The scaffold below
//  owns the contracted-pair loop, the late contraction, and the AO
//  scatter; it calls this routine once per primitive pair and expects the
//  *primitive* overlap integrals back (no contraction coefficients — the
//  scaffold applies those).
// ====================================================================

/// @brief Fills `primitive` with the primitive overlap integrals of one
/// primitive pair, for the `(la, lc)` shell pair.
///
/// `primitive` has length `(2*la+1)*(2*lc+1)`, row-major over the bra/ket
/// solid-harmonic components: element `(ma, mc)` is at `ma*(2*lc+1) + mc`.
///
/// STUB — currently zero-filled. The custom recursion replaces this body.
auto
compute_primitive_overlap(double*                              primitive,
                          const int                            la,
                          const int                            lc,
                          [[maybe_unused]] const double         alpha,
                          [[maybe_unused]] const double         beta,
                          [[maybe_unused]] const TPoint<double>& rA,
                          [[maybe_unused]] const TPoint<double>& rC) -> void
{
    const auto count = (2 * la + 1) * (2 * lc + 1);

    // TODO: custom recursion — primitive overlap for the (la, lc) shell pair.
    for (int k = 0; k < count; k++)
    {
        primitive[k] = 0.0;
    }
}

/// @brief Evaluates a screened `CGtoPairBlock` into the overlap matrix — the
/// late-contraction scaffold.
///
/// For each surviving contracted-GTO pair the primitive integrals are summed
/// over the primitive pairs, weighted by the pair's contraction coefficient
/// (`normalization_factors`, the combined bra/ket primitive weight), and the
/// contracted block is scattered into the global matrix. An off-diagonal
/// block pair also fills its transpose.
auto
compute_pair_block(DenseMatrix& matrix, const CGtoPairBlock& pairBlock, const bool diagonalPair) -> void
{
    const auto [la, lc] = pairBlock.angular_momentums();

    const auto cdim    = pairBlock.number_of_contracted_pairs();
    const auto nppairs = pairBlock.number_of_primitive_pairs();

    const auto ma = 2 * la + 1;
    const auto mc = 2 * lc + 1;

    // basis-function-pair-block SoA data
    const auto braCoords = pairBlock.bra_coordinates();
    const auto ketCoords = pairBlock.ket_coordinates();
    const auto braExps   = pairBlock.bra_exponents();
    const auto ketExps   = pairBlock.ket_exponents();
    const auto norms     = pairBlock.normalization_factors();

    // orbital indices — [0] is the AO component stride, [ij+1] the per-pair
    // global orbital offset
    const auto braOrb = pairBlock.bra_orbital_indices();
    const auto ketOrb = pairBlock.ket_orbital_indices();

    std::vector<double> primitive(static_cast<std::size_t>(ma * mc), 0.0);
    std::vector<double> contracted(static_cast<std::size_t>(ma * mc), 0.0);

    for (std::size_t ij = 0; ij < cdim; ij++)
    {
        std::fill(contracted.begin(), contracted.end(), 0.0);

        const auto rA = braCoords[ij];
        const auto rC = ketCoords[ij];

        // late contraction — sum the primitive integrals over the primitive
        // pairs, weighted by the contraction coefficients
        for (int pp = 0; pp < nppairs; pp++)
        {
            const auto ijoff = static_cast<std::size_t>(pp) * cdim + ij;

            const auto alpha  = braExps[ijoff];
            const auto beta   = ketExps[ijoff];
            const auto weight = norms[ijoff];

            compute_primitive_overlap(primitive.data(), la, lc, alpha, beta, rA, rC);

            for (int k = 0; k < ma * mc; k++)
            {
                contracted[k] += weight * primitive[k];
            }
        }

        // scatter the contracted block into the global matrix
        for (int ia = 0; ia < ma; ia++)
        {
            const auto row = static_cast<std::size_t>(ia) * braOrb[0] + braOrb[ij + 1];

            for (int ic = 0; ic < mc; ic++)
            {
                const auto col = static_cast<std::size_t>(ic) * ketOrb[0] + ketOrb[ij + 1];

                const auto value = contracted[ia * mc + ic];

                matrix(row, col) = value;

                // an off-diagonal block pair also fills its transpose; a
                // diagonal pair already visits both (a,c) and (c,a)
                if (!diagonalPair)
                {
                    matrix(col, row) = value;
                }
            }
        }
    }
}

}  // namespace

auto
OverlapDriver::compute(const CMolecule& molecule, const CMolecularBasis& basis, const double threshold) const -> DenseMatrix
{
    const auto gtoBlocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto dimension = static_cast<std::size_t>(gtofunc::getNumberOfAtomicOrbitals(gtoBlocks));

    DenseMatrix matrix(dimension, dimension, Symmetry::symmetric);

    // triangular loop over the basis-function-block pairs
    for (std::size_t i = 0; i < gtoBlocks.size(); i++)
    {
        for (std::size_t j = 0; j <= i; j++)
        {
            const auto estimator = make_overlap_screening_estimator(gtoBlocks[i].angular_momentum(),
                                                                    gtoBlocks[j].angular_momentum());

            // a screened pair block — the negligible contracted-GTO pairs are
            // dropped at construction
            const CGtoPairBlock pairBlock(gtoBlocks[i], gtoBlocks[j], estimator, threshold);

            compute_pair_block(matrix, pairBlock, i == j);
        }
    }

    return matrix;
}

}  // namespace tabula
