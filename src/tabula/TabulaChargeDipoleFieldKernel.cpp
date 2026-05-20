//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center charge-dipole field kernel — the transpose of the matrix kernel.
//  Same per-axis field integral and per-point seed; the density is folded into
//  the M-rows so the result is the electric field at each external point.
//

#include "TabulaChargeDipoleFieldKernel.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "MathConst.hpp"
#include "TabulaBoys.hpp"
#include "TabulaChargeDipoleTable.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

/// @brief The field is `−∇_Q` of the charge-dipole table (the same sign the
/// matrix kernel is validated with).
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
charge_dipole_field_kernel(const int              l_a,
                           const int              l_c,
                           const KernelBlockData &bra,
                           const int              bra_begin,
                           const int              bra_end,
                           const KernelBlockData &ket,
                           const double          *point_x,
                           const double          *point_y,
                           const double          *point_z,
                           const int              n_points,
                           const double          *density_block,
                           const double           weight,
                           double                *field) -> void
{
    const bool                       swap      = l_a > l_c;
    const int                        la_t      = swap ? l_c : l_a;
    const int                        lc_t      = swap ? l_a : l_c;
    const detail::ChargeDipoleTable *t[3]      = {&detail::chargedipole_table(0, la_t, lc_t),
                                                  &detail::chargedipole_table(1, la_t, lc_t),
                                                  &detail::chargedipole_table(2, la_t, lc_t)};
    const double *const              coef_dict = detail::chargedipole_coef_dict();
    const int                        np        = n_points;

    const double two_pi  = 2.0 * mathconst::pi_value();
    const int    max_rac = l_a + l_c;

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

    // shared per-point buffers
    thread_local std::vector<double> cb_qx, cb_qy, cb_qz, cb_T, cb_boys, cb_qxp, cb_qyp, cb_qzp;
    const auto grow = [](std::vector<double> &v, std::size_t n) { if (v.size() < n) v.resize(n); };
    grow(cb_qx, np);
    grow(cb_qy, np);
    grow(cb_qz, np);
    grow(cb_T, np);
    grow(cb_boys, static_cast<std::size_t>(max_order + 1) * np);
    grow(cb_qxp, static_cast<std::size_t>(max_q + 1) * np);
    grow(cb_qyp, static_cast<std::size_t>(max_q + 1) * np);
    grow(cb_qzp, static_cast<std::size_t>(max_q + 1) * np);

    // per-axis state: distinct Q-monomial map, per-point expanded-V, per-ev weight
    const int                                       qstride = max_q + 1;
    thread_local std::array<std::vector<int>, 3>    qmono_of_ev, qm_qx, qm_qy, qm_qz;
    thread_local std::array<std::vector<double>, 3> vn, qmprod, dweight;
    thread_local std::vector<int>                   qseen;
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
                idx                                  = static_cast<int>(qm_qx[a].size());
                qseen[static_cast<std::size_t>(key)] = idx;
                qm_qx[a].push_back(t[a]->ev_qx[e]);
                qm_qy[a].push_back(t[a]->ev_qy[e]);
                qm_qz[a].push_back(t[a]->ev_qz[e]);
            }
            qmono_of_ev[a][e] = idx;
        }
        const int dqm = static_cast<int>(qm_qx[a].size());
        for (int d = 0; d < dqm; d++)
        {
            const int key                        = (qm_qx[a][d] * qstride + qm_qy[a][d]) * qstride + qm_qz[a][d];
            qseen[static_cast<std::size_t>(key)] = -1;
        }
        grow(qmprod[a], static_cast<std::size_t>(dqm) * np);
        grow(vn[a], static_cast<std::size_t>(ev) * np);
        if (dweight[a].size() < static_cast<std::size_t>(ev)) dweight[a].resize(ev);
    }

    const int         ket_components = 2 * l_c + 1;
    const std::size_t cdim           = static_cast<std::size_t>(bra_end - bra_begin) * static_cast<std::size_t>(ket.ncgtos);

    for (int i = bra_begin; i < bra_end; i++)
    {
        const std::size_t i_local = static_cast<std::size_t>(i - bra_begin);
        const double      bx      = bra.x[i];
        const double      by      = bra.y[i];
        const double      bz      = bra.z[i];

        for (int jj = 0; jj < ket.ncgtos; jj++)
        {
            const double dx  = bx - ket.x[jj];
            const double dy  = by - ket.y[jj];
            const double dz  = bz - ket.z[jj];
            const double r2  = dx * dx + dy * dy + dz * dz;
            const double acx = swap ? -dx : dx;
            const double acy = swap ? -dy : dy;
            const double acz = swap ? -dz : dz;
            const double rnx = swap ? ket.x[jj] : bx;
            const double rny = swap ? ket.y[jj] : by;
            const double rnz = swap ? ket.z[jj] : bz;

            double powx[9], powy[9], powz[9];
            powx[0] = powy[0] = powz[0] = 1.0;
            for (int k = 1; k <= max_rac; k++)
            {
                powx[k] = powx[k - 1] * acx;
                powy[k] = powy[k - 1] * acy;
                powz[k] = powz[k - 1] * acz;
            }

            // ── Phase 1 — per-point expanded-V (kept per point) ──
            for (int a = 0; a < 3; a++)
                for (std::size_t e = 0; e < static_cast<std::size_t>(t[a]->ev_count) * np; e++) vn[a][e] = 0.0;

            for (int kb = 0; kb < bra.nprims; kb++)
            {
                const double a_exp = bra.exponents[static_cast<std::size_t>(kb) * bra.ncgtos + i];
                const double bn    = bra.norms[static_cast<std::size_t>(kb) * bra.ncgtos + i];

                for (int lk = 0; lk < ket.nprims; lk++)
                {
                    const double g    = ket.exponents[static_cast<std::size_t>(lk) * ket.ncgtos + jj];
                    const double kn   = ket.norms[static_cast<std::size_t>(lk) * ket.ncgtos + jj];
                    const double zeta  = a_exp + g;
                    const double mu    = a_exp * g / zeta;
                    const double alpha = swap ? g : a_exp;
                    const double gamma = swap ? a_exp : g;
                    const double goz   = gamma / zeta;
                    const double sac   = std::exp(-mu * r2);
                    const double pref  = (two_pi / zeta) * ipow(0.5 / alpha, la_t) * ipow(-0.5 / gamma, lc_t);
                    const double w0    = field_sign * bn * kn * pref * sac;

                    for (int p = 0; p < np; p++)
                    {
                        const double qx = (rnx - point_x[p]) - goz * acx;
                        const double qy = (rny - point_y[p]) - goz * acy;
                        const double qz = (rnz - point_z[p]) - goz * acz;
                        cb_qx[p] = qx;
                        cb_qy[p] = qy;
                        cb_qz[p] = qz;
                        cb_T[p]  = zeta * (qx * qx + qy * qy + qz * qz);
                    }

                    double boys_v[16];
                    for (int p = 0; p < np; p++)
                    {
                        boys(max_order, cb_T[p], boys_v);
                        for (int o = 0; o <= max_order; o++) cb_boys[static_cast<std::size_t>(o) * np + p] = boys_v[o];
                    }

                    for (int p = 0; p < np; p++)
                    {
                        cb_qxp[p] = 1.0;
                        cb_qyp[p] = 1.0;
                        cb_qzp[p] = 1.0;
                    }
                    for (int k = 1; k <= max_q; k++)
                    {
                        const std::size_t cur = static_cast<std::size_t>(k) * np;
                        const std::size_t prv = static_cast<std::size_t>(k - 1) * np;
                        for (int p = 0; p < np; p++)
                        {
                            cb_qxp[cur + p] = cb_qxp[prv + p] * cb_qx[p];
                            cb_qyp[cur + p] = cb_qyp[prv + p] * cb_qy[p];
                            cb_qzp[cur + p] = cb_qzp[prv + p] * cb_qz[p];
                        }
                    }

                    double       mupow[17], gozpow[17], seedpow[17];
                    const double neg2zeta = -2.0 * zeta;
                    mupow[0] = gozpow[0] = seedpow[0] = 1.0;
                    for (int k = 1; k <= max_mu; k++) mupow[k] = mupow[k - 1] * mu;
                    for (int k = 1; k <= max_goz; k++) gozpow[k] = gozpow[k - 1] * goz;
                    for (int k = 1; k <= max_order; k++) seedpow[k] = seedpow[k - 1] * neg2zeta;

                    for (int a = 0; a < 3; a++)
                    {
                        const detail::ChargeDipoleTable &ta  = *t[a];
                        const int                        dqm = static_cast<int>(qm_qx[a].size());
                        for (int d = 0; d < dqm; d++)
                        {
                            const double *Qx  = &cb_qxp[static_cast<std::size_t>(qm_qx[a][d]) * np];
                            const double *Qy  = &cb_qyp[static_cast<std::size_t>(qm_qy[a][d]) * np];
                            const double *Qz  = &cb_qzp[static_cast<std::size_t>(qm_qz[a][d]) * np];
                            double       *QMd = &qmprod[a][static_cast<std::size_t>(d) * np];
#pragma omp simd
                            for (int p = 0; p < np; p++) QMd[p] = Qx[p] * Qy[p] * Qz[p];
                        }

                        for (int e = 0; e < ta.ev_count; e++)
                        {
                            const double scale =
                                mupow[ta.ev_mu_power[e]] * gozpow[ta.ev_goz_power[e]] * seedpow[ta.ev_order[e]];
                            const double  c0  = w0 * scale;
                            const double *Bo  = &cb_boys[static_cast<std::size_t>(ta.ev_order[e]) * np];
                            const double *QMe = &qmprod[a][static_cast<std::size_t>(qmono_of_ev[a][e]) * np];
                            double       *VNe = &vn[a][static_cast<std::size_t>(e) * np];
#pragma omp simd
                            for (int p = 0; p < np; p++) VNe[p] += c0 * Bo[p] * QMe[p];
                        }
                    }
                }
            }

            // ── Phase 3 — fold the density into the M-rows, accumulate E ──
            for (int a = 0; a < 3; a++)
            {
                const detail::ChargeDipoleTable &ta = *t[a];
                for (int e = 0; e < ta.ev_count; e++) dweight[a][e] = 0.0;

                int r = 0;
                for (int k = 0; k < ta.component_count; k++)
                {
                    const int     braM    = ta.component_bra_m[k];
                    const int     ketM    = ta.component_ket_m[k];
                    const int     ma      = swap ? ketM : braM;
                    const int     mc      = swap ? braM : ketM;
                    const int     out_row = (ma + l_a) * ket_components + (mc + l_c);
                    const double  dval    = density_block[static_cast<std::size_t>(out_row) * cdim +
                                                       i_local * static_cast<std::size_t>(ket.ncgtos) +
                                                       static_cast<std::size_t>(jj)];

                    while (r < ta.m_row_count && ta.m_fields[6 * r] == k)
                    {
                        const int    ev_index = ta.m_fields[6 * r + 1] | (static_cast<int>(ta.m_fields[6 * r + 2]) << 8);
                        const int    rx       = ta.m_fields[6 * r + 3];
                        const int    ry       = ta.m_fields[6 * r + 4];
                        const int    rz       = ta.m_fields[6 * r + 5];
                        const double coef     = coef_dict[ta.m_coef_idx[r]];
                        dweight[a][ev_index] += coef * powx[rx] * powy[ry] * powz[rz] * dval;
                        r++;
                    }
                }

                double *const Fa = &field[static_cast<std::size_t>(a) * np];
                for (int e = 0; e < ta.ev_count; e++)
                {
                    const double dw = weight * dweight[a][e];
                    if (dw == 0.0) continue;
                    const double *VNe = &vn[a][static_cast<std::size_t>(e) * np];
#pragma omp simd
                    for (int p = 0; p < np; p++) Fa[p] += dw * VNe[p];
                }
            }
        }
    }
}

}  // namespace tabula
