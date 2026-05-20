//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center charge-dipole kernel — the three field axes evaluated in a single
//  pass: the per-charge geometry, Boys, Q-powers, and seed scales are shared
//  across the axes (they do not depend on the field axis), and only the
//  per-axis table fold and the M·V contraction differ. Separate from the
//  nuclear engine; reuses its three-phase shape, coef dictionary, and symmetry.
//

#include "TabulaChargeDipoleKernel.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "MathConst.hpp"
#include "TabulaBoys.hpp"
#include "TabulaChargeDipoleTable.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

/// @brief The ket contracted-GTO tile width.
constexpr int tile = 64;

/// @brief The field is `−∇_Q` of the charge-dipole table (the `∇_{R_N} = −∇_Q`
/// sign applied downstream of the recursion's `+1` Q-shift); pinned by
/// validation against `d·∇_R V`.
constexpr double field_sign = -1.0;

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
charge_dipole_kernel(const int              l_a,
                     const int              l_c,
                     const KernelBlockData &bra,
                     const int              bra_begin,
                     const int              bra_end,
                     const KernelBlockData &ket,
                     const DipoleSet       &dipoles,
                     double                *spherical) -> void
{
    // Only l_a ≤ l_c tables are stored; for l_a > l_c use the (l_c, l_a) tables
    // transposed (swapped bra/ket roles), as for nuclear attraction. The three
    // field axes share the same (l_a, l_c) angular structure — only the M rows
    // and expanded-V differ — so the per-charge work is computed once.
    const bool                       swap      = l_a > l_c;
    const int                        la_t      = swap ? l_c : l_a;
    const int                        lc_t      = swap ? l_a : l_c;
    const detail::ChargeDipoleTable *t[3]      = {&detail::chargedipole_table(0, la_t, lc_t),
                                                  &detail::chargedipole_table(1, la_t, lc_t),
                                                  &detail::chargedipole_table(2, la_t, lc_t)};
    const double *const              coef_dict = detail::chargedipole_coef_dict();
    const double *const              moment[3] = {dipoles.mx, dipoles.my, dipoles.mz};
    const int                        nc        = dipoles.count;

    const double two_pi  = 2.0 * mathconst::pi_value();
    const int    max_rac = l_a + l_c;  // shared R_AC degree (the ∇_Q shift is Q-only)

    // global maxima over the three tables — size the shared per-charge buffers
    int max_order = 0, max_mu = 0, max_goz = 0, max_q = 0;
    for (int a = 0; a < 3; a++)
    {
        if (t[a]->max_order > max_order) max_order = t[a]->max_order;
        for (int e = 0; e < t[a]->ev_count; e++)
        {
            if (t[a]->ev_mu_power[e] > max_mu) max_mu = t[a]->ev_mu_power[e];
            if (t[a]->ev_goz_power[e] > max_goz) max_goz = t[a]->ev_goz_power[e];
            if (t[a]->ev_qx[e] > max_q) max_q = t[a]->ev_qx[e];
            if (t[a]->ev_qy[e] > max_q) max_q = t[a]->ev_qy[e];
            if (t[a]->ev_qz[e] > max_q) max_q = t[a]->ev_qz[e];
        }
    }

    // shared per-charge buffers
    thread_local std::vector<double> cb_qx, cb_qy, cb_qz, cb_T, cb_boys, cb_qxp, cb_qyp, cb_qzp;
    const auto grow = [](std::vector<double> &v, std::size_t n) { if (v.size() < n) v.resize(n); };
    grow(cb_qx, nc);
    grow(cb_qy, nc);
    grow(cb_qz, nc);
    grow(cb_T, nc);
    grow(cb_boys, static_cast<std::size_t>(max_order + 1) * nc);
    grow(cb_qxp, static_cast<std::size_t>(max_q + 1) * nc);
    grow(cb_qyp, static_cast<std::size_t>(max_q + 1) * nc);
    grow(cb_qzp, static_cast<std::size_t>(max_q + 1) * nc);

    // per-axis state: the distinct Q-monomial map, and the fold scratch
    const int                                     qstride = max_q + 1;
    thread_local std::array<std::vector<int>, 3>  qmono_of_ev, qm_qx, qm_qy, qm_qz;
    thread_local std::array<std::vector<double>, 3> vg, ac, wbo, qmprod;
    thread_local std::vector<int>                 qseen;
    int                                           dqm[3];
    qseen.assign(static_cast<std::size_t>(qstride) * qstride * qstride, -1);
    for (int a = 0; a < 3; a++)
    {
        const int ev = t[a]->ev_count;
        qmono_of_ev[a].resize(ev);
        qm_qx[a].clear();
        qm_qy[a].clear();
        qm_qz[a].clear();
        for (int e = 0; e < ev; e++)
        {
            const int key = (t[a]->ev_qx[e] * qstride + t[a]->ev_qy[e]) * qstride + t[a]->ev_qz[e];
            int       idx = qseen[static_cast<std::size_t>(key)];
            if (idx < 0)
            {
                idx                            = static_cast<int>(qm_qx[a].size());
                qseen[static_cast<std::size_t>(key)] = idx;
                qm_qx[a].push_back(t[a]->ev_qx[e]);
                qm_qy[a].push_back(t[a]->ev_qy[e]);
                qm_qz[a].push_back(t[a]->ev_qz[e]);
            }
            qmono_of_ev[a][e] = idx;
        }
        dqm[a] = static_cast<int>(qm_qx[a].size());
        // reset the seen entries this axis touched (cheap vs re-clearing all)
        for (int d = 0; d < dqm[a]; d++)
        {
            const int key                       = (qm_qx[a][d] * qstride + qm_qy[a][d]) * qstride + qm_qz[a][d];
            qseen[static_cast<std::size_t>(key)] = -1;
        }
        grow(wbo[a], static_cast<std::size_t>(max_order + 1) * nc);
        grow(qmprod[a], static_cast<std::size_t>(dqm[a]) * nc);
        if (vg[a].size() < static_cast<std::size_t>(ev) * tile) vg[a].resize(static_cast<std::size_t>(ev) * tile);
        if (ac[a].size() < static_cast<std::size_t>(ev)) ac[a].resize(static_cast<std::size_t>(ev));
    }

    const std::size_t cdim           = static_cast<std::size_t>(bra_end - bra_begin) * static_cast<std::size_t>(ket.ncgtos);
    const std::size_t stride         = ((cdim + 7) / 8) * 8;
    const int         ket_components = 2 * l_c + 1;
    const std::size_t components     = static_cast<std::size_t>((2 * l_a + 1) * ket_components);

    // zero the whole output block once — the three axes accumulate into it
    for (std::size_t k = 0; k < components * stride; k++) spherical[k] = 0.0;

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

            // shared per-tile geometry (axis-independent)
            for (int jj = 0; jj < w; jj++)
            {
                const double dx = bx - ket.x[j0 + jj];
                const double dy = by - ket.y[j0 + jj];
                const double dz = bz - ket.z[j0 + jj];
                r2[jj]   = dx * dx + dy * dy + dz * dz;
                acx[jj]  = swap ? -dx : dx;
                acy[jj]  = swap ? -dy : dy;
                acz[jj]  = swap ? -dz : dz;
                ranx[jj] = swap ? ket.x[j0 + jj] : bx;
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

            // ── Phase 1 — charge sum into the per-axis expanded-V grids ──
            for (int a = 0; a < 3; a++)
                for (int e = 0; e < t[a]->ev_count; e++)
                    for (int jj = 0; jj < w; jj++) vg[a][e * tile + jj] = 0.0;

            for (int kb = 0; kb < bra.nprims; kb++)
            {
                const double a_exp = bra.exponents[static_cast<std::size_t>(kb) * bra.ncgtos + i];
                const double bn    = bra.norms[static_cast<std::size_t>(kb) * bra.ncgtos + i];

                for (int lk = 0; lk < ket.nprims; lk++)
                {
                    const double *ke = ket.exponents + static_cast<std::size_t>(lk) * ket.ncgtos + j0;
                    const double *kn = ket.norms + static_cast<std::size_t>(lk) * ket.ncgtos + j0;

                    for (int jj = 0; jj < w; jj++)
                    {
                        const double g     = ke[jj];
                        const double zeta  = a_exp + g;
                        const double mu    = a_exp * g / zeta;
                        const double alpha = swap ? g : a_exp;
                        const double gamma = swap ? a_exp : g;
                        const double goz   = gamma / zeta;
                        const double sac   = std::exp(-mu * r2[jj]);
                        const double pref  = (two_pi / zeta) * ipow(0.5 / alpha, la_t) * ipow(-0.5 / gamma, lc_t);
                        const double w0    = field_sign * bn * kn[jj] * pref * sac;

                        const double axj = acx[jj], ayj = acy[jj], azj = acz[jj];
                        const double rnx = ranx[jj], rny = rany[jj], rnz = ranz[jj];

                        // shared per-charge geometry: Q, T (axis-independent)
                        for (int c = 0; c < nc; c++)
                        {
                            const double qx = (rnx - dipoles.x[c]) - goz * axj;
                            const double qy = (rny - dipoles.y[c]) - goz * ayj;
                            const double qz = (rnz - dipoles.z[c]) - goz * azj;
                            cb_qx[c] = qx;
                            cb_qy[c] = qy;
                            cb_qz[c] = qz;
                            cb_T[c]  = zeta * (qx * qx + qy * qy + qz * qz);
                        }

                        // shared Boys (the expensive scalar) — once for all axes
                        for (int c = 0; c < nc; c++)
                        {
                            boys(max_order, cb_T[c], boys_v);
                            for (int o = 0; o <= max_order; o++) cb_boys[static_cast<std::size_t>(o) * nc + c] = boys_v[o];
                        }

                        // shared Q-powers
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

                        // shared seed/scale tables (mu, goz, zeta are axis-independent)
                        double       mupow[17], gozpow[17], seedpow[17];
                        const double neg2zeta = -2.0 * zeta;
                        mupow[0] = gozpow[0] = seedpow[0] = 1.0;
                        for (int k = 1; k <= max_mu; k++) mupow[k] = mupow[k - 1] * mu;
                        for (int k = 1; k <= max_goz; k++) gozpow[k] = gozpow[k - 1] * goz;
                        for (int k = 1; k <= max_order; k++) seedpow[k] = seedpow[k - 1] * neg2zeta;

                        // per-axis fold — distinct table, distinct dipole moment
                        for (int a = 0; a < 3; a++)
                        {
                            const detail::ChargeDipoleTable &ta = *t[a];
                            const double *const              wt = moment[a];

                            for (int o = 0; o <= ta.max_order; o++)
                            {
                                const double *Bo  = &cb_boys[static_cast<std::size_t>(o) * nc];
                                double       *wBo = &wbo[a][static_cast<std::size_t>(o) * nc];
#pragma omp simd
                                for (int c = 0; c < nc; c++) wBo[c] = wt[c] * Bo[c];
                            }

                            for (int d = 0; d < dqm[a]; d++)
                            {
                                const double *Qx  = &cb_qxp[static_cast<std::size_t>(qm_qx[a][d]) * nc];
                                const double *Qy  = &cb_qyp[static_cast<std::size_t>(qm_qy[a][d]) * nc];
                                const double *Qz  = &cb_qzp[static_cast<std::size_t>(qm_qz[a][d]) * nc];
                                double       *QMd = &qmprod[a][static_cast<std::size_t>(d) * nc];
#pragma omp simd
                                for (int c = 0; c < nc; c++) QMd[c] = Qx[c] * Qy[c] * Qz[c];
                            }

                            const int ev = ta.ev_count;
                            for (int e = 0; e < ev; e++)
                            {
                                const double *Bo  = &wbo[a][static_cast<std::size_t>(ta.ev_order[e]) * nc];
                                const double *QMe = &qmprod[a][static_cast<std::size_t>(qmono_of_ev[a][e]) * nc];
                                double        s   = 0.0;
#pragma omp simd reduction(+ : s)
                                for (int c = 0; c < nc; c++) s += Bo[c] * QMe[c];
                                ac[a][e] = s;
                            }

                            for (int e = 0; e < ev; e++)
                            {
                                const double scale =
                                    mupow[ta.ev_mu_power[e]] * gozpow[ta.ev_goz_power[e]] * seedpow[ta.ev_order[e]];
                                vg[a][e * tile + jj] += w0 * scale * ac[a][e];
                            }
                        }
                    }
                }
            }

            // ── Phase 3 — pure-R_AC M·V per axis, accumulated into spherical ──
            for (int a = 0; a < 3; a++)
            {
                const detail::ChargeDipoleTable &ta = *t[a];
                int                              r  = 0;
                for (int k = 0; k < ta.component_count; k++)
                {
                    const int     braM    = ta.component_bra_m[k];
                    const int     ketM    = ta.component_ket_m[k];
                    const int     ma      = swap ? ketM : braM;
                    const int     mc      = swap ? braM : ketM;
                    const int     out_row = (ma + l_a) * ket_components + (mc + l_c);
                    double *const out     = spherical + static_cast<std::size_t>(out_row) * stride + ij_base +
                                        static_cast<std::size_t>(j0);

                    while (r < ta.m_row_count && ta.m_fields[6 * r] == k)
                    {
                        const int     ev_index = ta.m_fields[6 * r + 1] | (static_cast<int>(ta.m_fields[6 * r + 2]) << 8);
                        const int     rx       = ta.m_fields[6 * r + 3];
                        const int     ry       = ta.m_fields[6 * r + 4];
                        const int     rz       = ta.m_fields[6 * r + 5];
                        const double  coef     = coef_dict[ta.m_coef_idx[r]];
                        const double *vev      = vg[a].data() + static_cast<std::size_t>(ev_index) * tile;
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
}

}  // namespace tabula
