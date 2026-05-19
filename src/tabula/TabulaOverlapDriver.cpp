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
#include "Point.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

// ====================================================================
//  RECURSION SEAM
//
//  This is where Tabula's custom recursion plugs in. The scaffold below
//  owns the block-pair loop, the late contraction, and the AO scatter; it
//  calls this routine once per primitive pair and expects the *primitive*
//  overlap integrals back (no contraction coefficients — the scaffold
//  applies those).
// ====================================================================

/// @brief Fills `primitive` with the primitive overlap integrals of one
/// primitive pair, for the `(la, lc)` shell pair.
///
/// `primitive` has length `(2*la+1) * (2*lc+1)`, row-major over the bra/ket
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

/// @brief Evaluates one bra/ket basis-function-block pair into the overlap
/// matrix — the late-contraction scaffold.
///
/// For each contracted-GTO pair the primitive integrals are summed over the
/// primitive pairs, weighted by the contraction coefficients (the block's
/// normalization factors), and the contracted block is scattered into the
/// global matrix. An off-diagonal block pair also fills its transpose.
auto
compute_block_pair(DenseMatrix&     matrix,
                   const CGtoBlock& braBlock,
                   const CGtoBlock& ketBlock,
                   const bool       diagonalPair) -> void
{
    const auto la = braBlock.angular_momentum();
    const auto lc = ketBlock.angular_momentum();

    const auto nca = braBlock.number_of_basis_functions();
    const auto ncc = ketBlock.number_of_basis_functions();
    const auto npa = braBlock.number_of_primitives();
    const auto npc = ketBlock.number_of_primitives();

    const auto ma = 2 * la + 1;
    const auto mc = 2 * lc + 1;

    // basis-function-block SoA data
    const auto braCoords = braBlock.coordinates();
    const auto ketCoords = ketBlock.coordinates();
    const auto braExps   = braBlock.exponents();
    const auto ketExps   = ketBlock.exponents();
    const auto braNorms  = braBlock.normalization_factors();
    const auto ketNorms  = ketBlock.normalization_factors();

    // global AO indices — component-major: (component, contracted GTO)
    const auto braAO = braBlock.getAtomicOrbitalsIndexes();
    const auto ketAO = ketBlock.getAtomicOrbitalsIndexes();

    std::vector<double> primitive(static_cast<std::size_t>(ma * mc), 0.0);
    std::vector<double> contracted(static_cast<std::size_t>(ma * mc), 0.0);

    for (int a = 0; a < nca; a++)
    {
        for (int c = 0; c < ncc; c++)
        {
            std::fill(contracted.begin(), contracted.end(), 0.0);

            // late contraction — sum the primitive integrals over the
            // primitive pairs, weighted by the contraction coefficients
            for (int p = 0; p < npa; p++)
            {
                const auto alpha = braExps[p * nca + a];
                const auto cA    = braNorms[p * nca + a];

                for (int q = 0; q < npc; q++)
                {
                    const auto beta = ketExps[q * ncc + c];
                    const auto cB   = ketNorms[q * ncc + c];

                    compute_primitive_overlap(primitive.data(), la, lc, alpha, beta, braCoords[a], ketCoords[c]);

                    const auto weight = cA * cB;

                    for (int k = 0; k < ma * mc; k++)
                    {
                        contracted[k] += weight * primitive[k];
                    }
                }
            }

            // scatter the contracted block into the global matrix
            for (int ia = 0; ia < ma; ia++)
            {
                const auto row = static_cast<std::size_t>(braAO[ia * nca + a]);

                for (int ic = 0; ic < mc; ic++)
                {
                    const auto col = static_cast<std::size_t>(ketAO[ic * ncc + c]);

                    const auto value = contracted[ia * mc + ic];

                    matrix(row, col) = value;

                    // an off-diagonal block pair also fills its transpose;
                    // a diagonal pair already visits both (a,c) and (c,a)
                    if (!diagonalPair)
                    {
                        matrix(col, row) = value;
                    }
                }
            }
        }
    }
}

}  // namespace

auto
OverlapDriver::compute(const CMolecule& molecule, const CMolecularBasis& basis) const -> DenseMatrix
{
    const auto gtoBlocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto dimension = static_cast<std::size_t>(gtofunc::getNumberOfAtomicOrbitals(gtoBlocks));

    DenseMatrix matrix(dimension, dimension, Symmetry::symmetric);

    // triangular loop over the basis-function-block pairs
    for (std::size_t i = 0; i < gtoBlocks.size(); i++)
    {
        for (std::size_t j = 0; j <= i; j++)
        {
            compute_block_pair(matrix, gtoBlocks[i], gtoBlocks[j], i == j);
        }
    }

    return matrix;
}

}  // namespace tabula
