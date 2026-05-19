//
//  Tabula — custom-recursion molecular-integral machinery.
//  Primitive-pair contraction.
//

#include "TabulaContraction.hpp"

namespace tabula {  // tabula namespace

auto
contract_primitive_pairs(const std::vector<double> &primitive,
                         const std::size_t          rows,
                         const std::size_t          cdim,
                         const std::size_t          nppairs) -> std::vector<double>
{
    // a single primitive pair per contracted pair — the contraction is the
    // identity (in/out row strides coincide), so the primitive buffer is
    // already the contracted buffer; skip the summation
    if (nppairs == 1) return primitive;

    // input / output row strides — padded to a multiple of 8 for aligned
    // SIMD, matching the seed-buffer convention (see compute_overlap_seed)
    const auto pdim       = cdim * nppairs;
    const auto in_stride  = ((pdim + 7) / 8) * 8;
    const auto out_stride = ((cdim + 7) / 8) * 8;

    std::vector<double> contracted(rows * out_stride, 0.0);

    if (cdim == 0 || nppairs == 0) return contracted;

    for (std::size_t m = 0; m < rows; m++)
    {
        const auto *in_row  = primitive.data() + m * in_stride;
        auto       *out_row = contracted.data() + m * out_stride;

        // primitive pair 0 — initialize the contracted row
#pragma omp simd
        for (std::size_t ij = 0; ij < cdim; ij++)
        {
            out_row[ij] = in_row[ij];
        }

        // primitive pairs 1 … nppairs-1 — accumulate
        for (std::size_t pp = 1; pp < nppairs; pp++)
        {
            const auto *block = in_row + pp * cdim;

#pragma omp simd
            for (std::size_t ij = 0; ij < cdim; ij++)
            {
                out_row[ij] += block[ij];
            }
        }
    }

    return contracted;
}

}  // namespace tabula
