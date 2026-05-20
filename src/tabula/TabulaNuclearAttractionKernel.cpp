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
    // Only l_a ≤ l_c tables are stored; for l_a > l_c use the (l_c, l_a) table
    // with the bra/ket roles swapped (the block is the transpose).
    const bool                  swap      = l_a > l_c;
    const detail::NuclearTable &t         = detail::nuclear_table(swap ? l_c : l_a, swap ? l_a : l_c);
    const double *const         coef_dict = detail::nuclear_coef_dict();

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

    thread_local std::vector<double> cb_qx, cb_qy, cb_qz, cb_T, cb_boys, cb_qxp, cb_qyp, cb_qzp, cb_chmbo, cb_qm;
    const auto grow = [](std::vector<double> &v, std::size_t n) { if (v.size() < n) v.resize(n); };
    grow(cb_qx, nc);
    grow(cb_qy, nc);
    grow(cb_qz, nc);
    grow(cb_T, nc);
    grow(cb_boys, static_cast<std::size_t>(t.max_order + 1) * nc);
    grow(cb_qxp, static_cast<std::size_t>(max_q + 1) * nc);
    grow(cb_qyp, static_cast<std::size_t>(max_q + 1) * nc);
    grow(cb_qzp, static_cast<std::size_t>(max_q + 1) * nc);

    // Distinct Q-monomial map (a per-(l_a,l_c) property; built once per call).
    // Many expanded-V entries share a Q-monomial (differing only in seed
    // order), so the charge sum precomputes the shared `chm·Boys` and Q-product
    // factors and reduces each ev with one multiply instead of four.
    const int                qstride = max_q + 1;
    thread_local std::vector<int> qmono_of_ev, qm_qx, qm_qy, qm_qz, qm_seen;
    qmono_of_ev.resize(ev);
    qm_qx.clear();
    qm_qy.clear();
    qm_qz.clear();
    qm_seen.assign(static_cast<std::size_t>(qstride) * qstride * qstride, -1);
    for (int e = 0; e < ev; e++)
    {
        const int key = (t.ev_qx[e] * qstride + t.ev_qy[e]) * qstride + t.ev_qz[e];
        int       idx = qm_seen[static_cast<std::size_t>(key)];
        if (idx < 0)
        {
            idx               = static_cast<int>(qm_qx.size());
            qm_seen[static_cast<std::size_t>(key)] = idx;
            qm_qx.push_back(t.ev_qx[e]);
            qm_qy.push_back(t.ev_qy[e]);
            qm_qz.push_back(t.ev_qz[e]);
        }
        qmono_of_ev[e] = idx;
    }
    const int dqm = static_cast<int>(qm_qx.size());
    grow(cb_chmbo, static_cast<std::size_t>(t.max_order + 1) * nc);
    grow(cb_qm, static_cast<std::size_t>(dqm) * nc);
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

    // per-tile lane geometry and R_AC power tables. `acx` holds the *table*
    // R_AC (flipped when swapped); `ranx` the table-bra centre for R_AN.
    double acx[tile], acy[tile], acz[tile], r2[tile];
    double ranx[tile], rany[tile], ranz[tile];
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
                const double dx = bx - ket.x[j0 + jj];  // A − C
                const double dy = by - ket.y[j0 + jj];
                const double dz = bz - ket.z[j0 + jj];
                r2[jj]   = dx * dx + dy * dy + dz * dz;
                acx[jj]  = swap ? -dx : dx;              // table R_AC
                acy[jj]  = swap ? -dy : dy;
                acz[jj]  = swap ? -dz : dz;
                ranx[jj] = swap ? ket.x[j0 + jj] : bx;  // table-bra centre (R_AN)
                rany[jj] = swap ? ket.y[j0 + jj] : by;
                ranz[jj] = swap ? ket.z[j0 + jj] : bz;
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
                const double a  = bra.exponents[static_cast<std::size_t>(kb) * bra.ncgtos + i];
                const double bn = bra.norms[static_cast<std::size_t>(kb) * bra.ncgtos + i];

                for (int lk = 0; lk < ket.nprims; lk++)
                {
                    const double *ke = ket.exponents + static_cast<std::size_t>(lk) * ket.ncgtos + j0;
                    const double *kn = ket.norms + static_cast<std::size_t>(lk) * ket.ncgtos + j0;

                    // scalar over ket lanes — the Boys function is scalar
                    for (int jj = 0; jj < w; jj++)
                    {
                        const double g     = ke[jj];
                        const double zeta  = a + g;
                        const double mu    = a * g / zeta;     // reduced exponent aγ/ζ (symmetric)
                        // bra/ket exponents in the table's convention — swapped
                        // for l_a > l_c, where the (l_c, l_a) table is used
                        const double alpha = swap ? g : a;
                        const double gamma = swap ? a : g;
                        const double goz   = gamma / zeta;     // (table γ)/ζ
                        const double sac   = std::exp(-mu * r2[jj]);
                        const double pref  = (two_pi / zeta) * ipow(0.5 / alpha, t.l_a) * ipow(-0.5 / gamma, t.l_c);
                        const double w0    = bn * kn[jj] * pref * sac;

                        const double axj = acx[jj], ayj = acy[jj], azj = acz[jj];
                        const double rnx = ranx[jj], rny = rany[jj], rnz = ranz[jj];

                        // Q and T per charge — Q = R_AN − (γ/ζ)·R_AC (table convention)
                        for (int c = 0; c < nc; c++)
                        {
                            const double qx = (rnx - chx[c]) - goz * axj;
                            const double qy = (rny - chy[c]) - goz * ayj;
                            const double qz = (rnz - chz[c]) - goz * azj;
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

                        // fold the charge magnitude into Boys, once per order
                        for (int o = 0; o <= t.max_order; o++)
                        {
                            const double     *Bo  = &cb_boys[static_cast<std::size_t>(o) * nc];
                            double           *cBo = &cb_chmbo[static_cast<std::size_t>(o) * nc];
#pragma omp simd
                            for (int c = 0; c < nc; c++) cBo[c] = chm[c] * Bo[c];
                        }

                        // the Q-product, once per distinct Q-monomial
                        for (int d = 0; d < dqm; d++)
                        {
                            const double *Qx  = &cb_qxp[static_cast<std::size_t>(qm_qx[d]) * nc];
                            const double *Qy  = &cb_qyp[static_cast<std::size_t>(qm_qy[d]) * nc];
                            const double *Qz  = &cb_qzp[static_cast<std::size_t>(qm_qz[d]) * nc];
                            double       *QMd = &cb_qm[static_cast<std::size_t>(d) * nc];
#pragma omp simd
                            for (int c = 0; c < nc; c++) QMd[c] = Qx[c] * Qy[c] * Qz[c];
                        }

                        // charge sum per ev — one multiply per term
                        for (int e = 0; e < ev; e++)
                        {
                            const double *Bo  = &cb_chmbo[static_cast<std::size_t>(t.ev_order[e]) * nc];
                            const double *QMe = &cb_qm[static_cast<std::size_t>(qmono_of_ev[e]) * nc];
                            double        s   = 0.0;
#pragma omp simd reduction(+ : s)
                            for (int c = 0; c < nc; c++) s += Bo[c] * QMe[c];
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
                // map the table component to our (l_a, l_c) output; swapped, the
                // table's (bra, ket) projections are our (l_c, l_a) → transpose
                const int     braM    = t.component_bra_m[k];
                const int     ketM    = t.component_ket_m[k];
                const int     ma      = swap ? ketM : braM;
                const int     mc      = swap ? braM : ketM;
                const int     out_row = (ma + l_a) * ket_components + (mc + l_c);
                double *const out     = spherical + static_cast<std::size_t>(out_row) * stride + ij_base +
                                    static_cast<std::size_t>(j0);

                for (int jj = 0; jj < w; jj++) out[jj] = 0.0;

                while (r < t.m_row_count && static_cast<int>(t.m_fields[r] & 0x7Fu) == k)
                {
                    const auto    f        = t.m_fields[r];  // component | ev<<7 | racx<<19 | racy<<23 | racz<<27
                    const int     ev_index = static_cast<int>((f >> 7) & 0xFFFu);
                    const int     rx       = static_cast<int>((f >> 19) & 0xFu);
                    const int     ry       = static_cast<int>((f >> 23) & 0xFu);
                    const int     rz       = static_cast<int>((f >> 27) & 0xFu);
                    const double  coef     = coef_dict[t.m_coef_idx[r]];
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
