//
//  Tabula — custom-recursion molecular-integral machinery.
//  Tabula-owned basis-function-pair block.
//

#include "TabulaGtoPairBlock.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace tabula {  // tabula namespace

GtoPairBlock::GtoPairBlock(const CGtoBlock &bra_gto_block, const CGtoBlock &ket_gto_block)
{
    _build(bra_gto_block, ket_gto_block, nullptr, 0.0, 0, bra_gto_block.number_of_basis_functions());
}

GtoPairBlock::GtoPairBlock(const CGtoBlock          &bra_gto_block,
                           const CGtoBlock          &ket_gto_block,
                           const ScreeningEstimator &estimator,
                           const double              threshold)
{
    _build(bra_gto_block, ket_gto_block, &estimator, threshold, 0, bra_gto_block.number_of_basis_functions());
}

GtoPairBlock::GtoPairBlock(const CGtoBlock &bra_gto_block, const CGtoBlock &ket_gto_block, const int bra_begin, const int bra_end)
{
    _build(bra_gto_block, ket_gto_block, nullptr, 0.0, bra_begin, bra_end);
}

auto
GtoPairBlock::_build(const CGtoBlock          &bra_gto_block,
                     const CGtoBlock          &ket_gto_block,
                     const ScreeningEstimator *estimator,
                     const double              threshold,
                     const int                 bra_begin,
                     const int                 bra_end) -> void
{
    const auto fpi = mathconst::pi_value();

    // fetch the bra/ket block data

    const auto bcoords = bra_gto_block.coordinates();
    const auto bexps   = bra_gto_block.exponents();
    const auto bnorms  = bra_gto_block.normalization_factors();
    const auto borbidx = bra_gto_block.orbital_indices();

    const auto kcoords = ket_gto_block.coordinates();
    const auto kexps   = ket_gto_block.exponents();
    const auto knorms  = ket_gto_block.normalization_factors();
    const auto korbidx = ket_gto_block.orbital_indices();

    const int bcgtos = bra_gto_block.number_of_basis_functions();
    const int bpgtos = bra_gto_block.number_of_primitives();
    const int kcgtos = ket_gto_block.number_of_basis_functions();
    const int kpgtos = ket_gto_block.number_of_primitives();

    // the bra contracted-GTO range this pair block covers
    const int bra_count = bra_end - bra_begin;

    _angular_momentums = {bra_gto_block.angular_momentum(), ket_gto_block.angular_momentum()};
    _nppairs           = bpgtos * kpgtos;

    // a threshold at or below 0 keeps every pair — the estimator pre-pass is
    // skipped entirely

    const bool screened = (estimator != nullptr) && (threshold > 0.0);

    // the screening pre-pass — collect the surviving contracted-GTO pairs

    std::vector<std::pair<int, int>> pairs;
    if (screened)
    {
        for (int i = bra_begin; i < bra_end; i++)
        {
            for (int j = 0; j < kcgtos; j++)
            {
                const auto r = std::sqrt(bcoords[i].distance_square(kcoords[j]));

                const auto estimate = (*estimator)(bra_gto_block.screening_data(static_cast<std::size_t>(i)),
                                                   ket_gto_block.screening_data(static_cast<std::size_t>(j)),
                                                   r);

                if (estimate >= threshold) pairs.push_back({i, j});
            }
        }
    }

    const std::size_t cdim = screened ? pairs.size() : static_cast<std::size_t>(bra_count) * static_cast<std::size_t>(kcgtos);
    const std::size_t pdim = cdim * static_cast<std::size_t>(_nppairs);

    // the contracted-pair data — coordinates and AO indices; element 0 of the
    // orbital-index vectors is the AO component stride

    _bra_coordinates.assign(cdim, TPoint<double>({0.0, 0.0, 0.0}));
    _ket_coordinates.assign(cdim, TPoint<double>({0.0, 0.0, 0.0}));
    _bra_orb_indices.assign(cdim + 1, 0);
    _ket_orb_indices.assign(cdim + 1, 0);
    _bra_orb_indices[0] = borbidx[0];
    _ket_orb_indices[0] = korbidx[0];

    // the primitive-pair arrays — allocated uninitialized; the fill loop
    // writes every element

    _bra_exponents = std::make_unique_for_overwrite<double[]>(pdim);
    _ket_exponents = std::make_unique_for_overwrite<double[]>(pdim);
    _weights       = std::make_unique_for_overwrite<double[]>(pdim);

    if (cdim == 0) return;

    // the per-contracted-pair coordinates, AO indices, and squared center
    // distance

    std::vector<double> r2ab(cdim);
    for (std::size_t ij = 0; ij < cdim; ij++)
    {
        const int i = screened ? pairs[ij].first : bra_begin + static_cast<int>(ij) / kcgtos;
        const int j = screened ? pairs[ij].second : static_cast<int>(ij) % kcgtos;

        _bra_coordinates[ij]     = bcoords[i];
        _ket_coordinates[ij]     = kcoords[j];
        _bra_orb_indices[ij + 1] = borbidx[i + 1];
        _ket_orb_indices[ij + 1] = korbidx[j + 1];
        r2ab[ij]                 = bcoords[i].distance_square(kcoords[j]);
    }

    // the primitive-pair fill — primitive-pair-major, the value for primitive
    // pair pp of contracted pair ij at pp·cdim + ij. The overlap factor
    // (π/(α+γ))^(3/2)·exp(−ρR²) folds the normalization factor in, so the
    // seed step is a plain product.

    double *bra_exp = _bra_exponents.get();
    double *ket_exp = _ket_exponents.get();
    double *weight  = _weights.get();

    if (!screened)
    {
        // unscreened — ij = i·kcgtos + j, so the inner j sweep writes a
        // contiguous run and the ket data loads contiguously
        for (int k = 0; k < bpgtos; k++)
        {
            for (int l = 0; l < kpgtos; l++)
            {
                const std::size_t pp = (static_cast<std::size_t>(k) * static_cast<std::size_t>(kpgtos) + static_cast<std::size_t>(l)) * cdim;

                for (int i = bra_begin; i < bra_end; i++)
                {
                    const double      bexp = bexps[k * bcgtos + i];
                    const double      bnrm = bnorms[k * bcgtos + i];
                    const std::size_t base = pp + static_cast<std::size_t>(i - bra_begin) * static_cast<std::size_t>(kcgtos);

                    const double *kexp_row = kexps.data() + l * kcgtos;
                    const double *knrm_row = knorms.data() + l * kcgtos;
                    const double *r2_row   = r2ab.data() + static_cast<std::size_t>(i - bra_begin) * static_cast<std::size_t>(kcgtos);

#pragma omp simd
                    for (int j = 0; j < kcgtos; j++)
                    {
                        const double g    = kexp_row[j];
                        const double fe   = 1.0 / (bexp + g);
                        const double fz   = bexp * g * fe;
                        const double fact = fpi * fe;

                        bra_exp[base + j] = bexp;
                        ket_exp[base + j] = g;
                        weight[base + j]  = bnrm * knrm_row[j] * fact * std::sqrt(fact) * std::exp(-fz * r2_row[j]);
                    }
                }
            }
        }
    }
    else
    {
        // screened — the surviving contracted pairs are an arbitrary subset
        for (int k = 0; k < bpgtos; k++)
        {
            for (int l = 0; l < kpgtos; l++)
            {
                const std::size_t pp = (static_cast<std::size_t>(k) * static_cast<std::size_t>(kpgtos) + static_cast<std::size_t>(l)) * cdim;

                for (std::size_t ij = 0; ij < cdim; ij++)
                {
                    const int i = pairs[ij].first;
                    const int j = pairs[ij].second;

                    const double bexp = bexps[k * bcgtos + i];
                    const double g    = kexps[l * kcgtos + j];
                    const double fe   = 1.0 / (bexp + g);
                    const double fz   = bexp * g * fe;
                    const double fact = fpi * fe;

                    bra_exp[pp + ij] = bexp;
                    ket_exp[pp + ij] = g;
                    weight[pp + ij]  = bnorms[k * bcgtos + i] * knorms[l * kcgtos + j] * fact * std::sqrt(fact) * std::exp(-fz * r2ab[ij]);
                }
            }
        }
    }
}

auto
GtoPairBlock::angular_momentums() const -> std::pair<int, int>
{
    return _angular_momentums;
}

auto
GtoPairBlock::number_of_contracted_pairs() const -> std::size_t
{
    return _bra_coordinates.size();
}

auto
GtoPairBlock::number_of_primitive_pairs() const -> int
{
    return _nppairs;
}

auto
GtoPairBlock::bra_coordinates() const -> const std::vector<TPoint<double>> &
{
    return _bra_coordinates;
}

auto
GtoPairBlock::ket_coordinates() const -> const std::vector<TPoint<double>> &
{
    return _ket_coordinates;
}

auto
GtoPairBlock::bra_orbital_indices() const -> const std::vector<std::size_t> &
{
    return _bra_orb_indices;
}

auto
GtoPairBlock::ket_orbital_indices() const -> const std::vector<std::size_t> &
{
    return _ket_orb_indices;
}

auto
GtoPairBlock::bra_exponents() const -> const double *
{
    return _bra_exponents.get();
}

auto
GtoPairBlock::ket_exponents() const -> const double *
{
    return _ket_exponents.get();
}

auto
GtoPairBlock::weights() const -> const double *
{
    return _weights.get();
}

}  // namespace tabula
