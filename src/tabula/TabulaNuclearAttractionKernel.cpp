//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center nuclear-attraction kernel — the shared, table-driven
//  three-phase engine that interprets the generated `(l_a, l_c)` tables.
//

#include "TabulaNuclearAttractionKernel.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "MathConst.hpp"
#include "TabulaBoys.hpp"
#include "TabulaNuclearAttractionTable.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

/// @brief The ket contracted-GTO tile width. Caps the expanded-V scratch
/// (`ev_count · tile` doubles, thread-local, grow-only) and is the phase-3
/// SIMD axis.
constexpr int tile = 64;

/// @brief `x` raised to a small non-negative integer power.
inline auto
ipow(const double x, const int n) -> double
{
    double result = 1.0;
    for (int i = 0; i < n; i++) result *= x;
    return result;
}

}  // namespace

auto
nuclear_attraction_kernel(const int              l_a,
                          const int              l_c,
                          const KernelBlockData &bra,
                          const int              bra_begin,
                          const int              bra_end,
                          const KernelBlockData &ket,
                          const ChargeSet       &charges,
                          double                *spherical) -> void
{
    const detail::NuclearTable &t = detail::nuclear_table(l_a, l_c);

    const double      two_pi  = 2.0 * mathconst::pi_value();
    const int         ev      = t.ev_count;
    const int         max_rac = l_a + l_c;       // largest R_AC power
    // The v- and u-recurrences each lift a Q per reduction, so the Q-monomial
    // degree reaches a + c = l_a + l_c (not just the bra's l_a).
    const int         max_q   = l_a + l_c;

    // largest µ / (γ/ζ) powers over the expanded V — for the per-(prim,lane)
    // power tables that replace the V-fold's per-ev ipow calls
    int max_mu = 0, max_goz = 0;
    for (int e = 0; e < ev; e++)
    {
        if (t.ev_mu_power[e] > max_mu) max_mu = t.ev_mu_power[e];
        if (t.ev_goz_power[e] > max_goz) max_goz = t.ev_goz_power[e];
    }

    // charge data and per-(prim,lane) charge-batch buffers, laid out
    // [order][charge] / [k][charge] so the per-ev reduction over charges
    // vectorizes. Thread-local, grow-only.
    const int     nc  = charges.count;
    const double *chm = charges.magnitudes;
    const double *chx = charges.x;
    const double *chy = charges.y;
    const double *chz = charges.z;

    thread_local std::vector<double> cb_qx, cb_qy, cb_qz, cb_T, cb_boys, cb_qxp, cb_qyp, cb_qzp;
    const auto grow = [](std::vector<double> &v, std::size_t n) { if (v.size() < n) v.resize(n); };
    grow(cb_qx, nc);
    grow(cb_qy, nc);
    grow(cb_qz, nc);
    grow(cb_T, nc);
    grow(cb_boys, static_cast<std::size_t>(t.max_order + 1) * nc);
    grow(cb_qxp, static_cast<std::size_t>(max_q + 1) * nc);
    grow(cb_qyp, static_cast<std::size_t>(max_q + 1) * nc);
    grow(cb_qzp, static_cast<std::size_t>(max_q + 1) * nc);
    const std::size_t cdim    = static_cast<std::size_t>(bra_end - bra_begin) * static_cast<std::size_t>(ket.ncgtos);
    const std::size_t stride  = ((cdim + 7) / 8) * 8;
    const int         ket_components = 2 * l_c + 1;

    // expanded-V grid (ev × tile) and the per-lane-charge accumulator (ev) —
    // thread-local, grow-only
    thread_local std::vector<double> v_grid;
    thread_local std::vector<double> acc;
    if (v_grid.size() < static_cast<std::size_t>(ev) * tile) v_grid.resize(static_cast<std::size_t>(ev) * tile);
    if (acc.size() < static_cast<std::size_t>(ev)) acc.resize(static_cast<std::size_t>(ev));
    double *const vg = v_grid.data();
    double *const ac = acc.data();

    // per-tile lane geometry and R_AC power tables
    double acx[tile], acy[tile], acz[tile], r2[tile];
    double powx[9 * tile], powy[9 * tile], powz[9 * tile];

    double boys_v[16];

    for (int i = bra_begin; i < bra_end; i++)
    {
        const double      bx      = bra.x[i];
        const double      by      = bra.y[i];
        const double      bz      = bra.z[i];
        const std::size_t ij_base = static_cast<std::size_t>(i - bra_begin) * static_cast<std::size_t>(ket.ncgtos);

        for (int j0 = 0; j0 < ket.ncgtos; j0 += tile)
        {
            const int w = (j0 + tile < ket.ncgtos) ? tile : (ket.ncgtos - j0);

            // AC = A − C, per ket contracted GTO of the tile, and its powers
            for (int jj = 0; jj < w; jj++)
            {
                const double dx = bx - ket.x[j0 + jj];
                const double dy = by - ket.y[j0 + jj];
                const double dz = bz - ket.z[j0 + jj];
                acx[jj] = dx;
                acy[jj] = dy;
                acz[jj] = dz;
                r2[jj]  = dx * dx + dy * dy + dz * dz;
                powx[jj] = 1.0;
                powy[jj] = 1.0;
                powz[jj] = 1.0;
            }
            for (int k = 1; k <= max_rac; k++)
            {
                for (int jj = 0; jj < w; jj++)
                {
                    powx[k * tile + jj] = powx[(k - 1) * tile + jj] * acx[jj];
                    powy[k * tile + jj] = powy[(k - 1) * tile + jj] * acy[jj];
                    powz[k * tile + jj] = powz[(k - 1) * tile + jj] * acz[jj];
                }
            }

            // ── Phase 1 — charge sum into the contracted expanded-V grid ──
            for (int e = 0; e < ev; e++)
                for (int jj = 0; jj < w; jj++) vg[e * tile + jj] = 0.0;

            for (int kb = 0; kb < bra.nprims; kb++)
            {
                const double a      = bra.exponents[static_cast<std::size_t>(kb) * bra.ncgtos + i];
                const double bn     = bra.norms[static_cast<std::size_t>(kb) * bra.ncgtos + i];
                const double ia_pow = ipow(0.5 / a, l_a);  // (1/2α)^l_a

                for (int lk = 0; lk < ket.nprims; lk++)
                {
                    const double *ke = ket.exponents + static_cast<std::size_t>(lk) * ket.ncgtos + j0;
                    const double *kn = ket.norms + static_cast<std::size_t>(lk) * ket.ncgtos + j0;

                    // scalar over ket lanes — the Boys function is scalar
                    for (int jj = 0; jj < w; jj++)
                    {
                        const double g    = ke[jj];
                        const double zeta = a + g;
                        const double mu   = a * g / zeta;     // reduced exponent aγ/ζ
                        const double goz  = g / zeta;         // γ/ζ
                        const double sac  = std::exp(-mu * r2[jj]);
                        const double pref = (two_pi / zeta) * ia_pow * ipow(-0.5 / g, l_c);
                        const double w0   = bn * kn[jj] * pref * sac;

                        const double axj = acx[jj], ayj = acy[jj], azj = acz[jj];

                        // Q and T per charge
                        for (int c = 0; c < nc; c++)
                        {
                            const double qx = (bx - chx[c]) - goz * axj;
                            const double qy = (by - chy[c]) - goz * ayj;
                            const double qz = (bz - chz[c]) - goz * azj;
                            cb_qx[c] = qx;
                            cb_qy[c] = qy;
                            cb_qz[c] = qz;
                            cb_T[c]  = zeta * (qx * qx + qy * qy + qz * qz);
                        }

                        // Boys per charge, laid out [order][charge]
                        for (int c = 0; c < nc; c++)
                        {
                            boys(t.max_order, cb_T[c], boys_v);
                            for (int o = 0; o <= t.max_order; o++) cb_boys[static_cast<std::size_t>(o) * nc + c] = boys_v[o];
                        }

                        // Q-powers, laid out [k][charge]
                        for (int c = 0; c < nc; c++)
                        {
                            cb_qxp[c] = 1.0;
                            cb_qyp[c] = 1.0;
                            cb_qzp[c] = 1.0;
                        }
                        for (int k = 1; k <= max_q; k++)
                        {
                            const std::size_t cur = static_cast<std::size_t>(k) * nc;
                            const std::size_t prv = static_cast<std::size_t>(k - 1) * nc;
                            for (int c = 0; c < nc; c++)
                            {
                                cb_qxp[cur + c] = cb_qxp[prv + c] * cb_qx[c];
                                cb_qyp[cur + c] = cb_qyp[prv + c] * cb_qy[c];
                                cb_qzp[cur + c] = cb_qzp[prv + c] * cb_qz[c];
                            }
                        }

                        // charge sum per ev — vectorized reduction over charges
                        for (int e = 0; e < ev; e++)
                        {
                            const double *Bo = &cb_boys[static_cast<std::size_t>(t.ev_order[e]) * nc];
                            const double *Qx = &cb_qxp[static_cast<std::size_t>(t.ev_qx[e]) * nc];
                            const double *Qy = &cb_qyp[static_cast<std::size_t>(t.ev_qy[e]) * nc];
                            const double *Qz = &cb_qzp[static_cast<std::size_t>(t.ev_qz[e]) * nc];
                            double s = 0.0;
#pragma omp simd reduction(+ : s)
                            for (int c = 0; c < nc; c++) s += chm[c] * Bo[c] * Qx[c] * Qy[c] * Qz[c];
                            ac[e] = s;
                        }

                        // power tables — replace the per-ev ipow in the fold
                        double       mupow[17], gozpow[17], seedpow[17];
                        const double neg2zeta = -2.0 * zeta;
                        mupow[0] = gozpow[0] = seedpow[0] = 1.0;
                        for (int k = 1; k <= max_mu; k++) mupow[k] = mupow[k - 1] * mu;
                        for (int k = 1; k <= max_goz; k++) gozpow[k] = gozpow[k - 1] * goz;
                        for (int k = 1; k <= t.max_order; k++) seedpow[k] = seedpow[k - 1] * neg2zeta;

                        for (int e = 0; e < ev; e++)
                        {
                            const double scale = mupow[t.ev_mu_power[e]] * gozpow[t.ev_goz_power[e]] * seedpow[t.ev_order[e]];
                            vg[e * tile + jj] += w0 * scale * ac[e];
                        }
                    }
                }
            }

            // ── Phase 3 — pure-R_AC M·V per contracted pair → spherical ──
            int r = 0;
            for (int k = 0; k < t.component_count; k++)
            {
                const int     out_row = (t.component_bra_m[k] + l_a) * ket_components + (t.component_ket_m[k] + l_c);
                double *const out     = spherical + static_cast<std::size_t>(out_row) * stride + ij_base +
                                    static_cast<std::size_t>(j0);

                for (int jj = 0; jj < w; jj++) out[jj] = 0.0;

                while (r < t.m_row_count && t.m_fields[6 * r] == k)
                {
                    const int     ev_index = t.m_fields[6 * r + 1] | (static_cast<int>(t.m_fields[6 * r + 2]) << 8);
                    const int     rx       = t.m_fields[6 * r + 3];
                    const int     ry       = t.m_fields[6 * r + 4];
                    const int     rz       = t.m_fields[6 * r + 5];
                    const double  coef     = t.m_coef[r];
                    const double *vev      = vg + static_cast<std::size_t>(ev_index) * tile;
                    const double *px       = powx + static_cast<std::size_t>(rx) * tile;
                    const double *py       = powy + static_cast<std::size_t>(ry) * tile;
                    const double *pz       = powz + static_cast<std::size_t>(rz) * tile;

#pragma omp simd
                    for (int jj = 0; jj < w; jj++) out[jj] += coef * px[jj] * py[jj] * pz[jj] * vev[jj];

                    r++;
                }
            }
        }
    }
}

}  // namespace tabula
