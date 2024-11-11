#include "NuclearPotentialPrimRecGS.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_gs(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_gs,
                               const size_t              idx_npot_0_ds,
                               const size_t              idx_npot_1_ds,
                               const size_t              idx_npot_0_fs,
                               const size_t              idx_npot_1_fs,
                               const CSimdArray<double>& factors,
                               const size_t              idx_rpa,
                               const size_t              idx_rpc,
                               const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : DS

    auto ta_xx_0_0 = pbuffer.data(idx_npot_0_ds);

    auto ta_yy_0_0 = pbuffer.data(idx_npot_0_ds + 3);

    auto ta_zz_0_0 = pbuffer.data(idx_npot_0_ds + 5);

    // Set up components of auxiliary buffer : DS

    auto ta_xx_0_1 = pbuffer.data(idx_npot_1_ds);

    auto ta_yy_0_1 = pbuffer.data(idx_npot_1_ds + 3);

    auto ta_zz_0_1 = pbuffer.data(idx_npot_1_ds + 5);

    // Set up components of auxiliary buffer : FS

    auto ta_xxx_0_0 = pbuffer.data(idx_npot_0_fs);

    auto ta_xxz_0_0 = pbuffer.data(idx_npot_0_fs + 2);

    auto ta_xyy_0_0 = pbuffer.data(idx_npot_0_fs + 3);

    auto ta_xzz_0_0 = pbuffer.data(idx_npot_0_fs + 5);

    auto ta_yyy_0_0 = pbuffer.data(idx_npot_0_fs + 6);

    auto ta_yyz_0_0 = pbuffer.data(idx_npot_0_fs + 7);

    auto ta_yzz_0_0 = pbuffer.data(idx_npot_0_fs + 8);

    auto ta_zzz_0_0 = pbuffer.data(idx_npot_0_fs + 9);

    // Set up components of auxiliary buffer : FS

    auto ta_xxx_0_1 = pbuffer.data(idx_npot_1_fs);

    auto ta_xxz_0_1 = pbuffer.data(idx_npot_1_fs + 2);

    auto ta_xyy_0_1 = pbuffer.data(idx_npot_1_fs + 3);

    auto ta_xzz_0_1 = pbuffer.data(idx_npot_1_fs + 5);

    auto ta_yyy_0_1 = pbuffer.data(idx_npot_1_fs + 6);

    auto ta_yyz_0_1 = pbuffer.data(idx_npot_1_fs + 7);

    auto ta_yzz_0_1 = pbuffer.data(idx_npot_1_fs + 8);

    auto ta_zzz_0_1 = pbuffer.data(idx_npot_1_fs + 9);

    // Set up components of targeted buffer : GS

    auto ta_xxxx_0_0 = pbuffer.data(idx_npot_0_gs);

    auto ta_xxxy_0_0 = pbuffer.data(idx_npot_0_gs + 1);

    auto ta_xxxz_0_0 = pbuffer.data(idx_npot_0_gs + 2);

    auto ta_xxyy_0_0 = pbuffer.data(idx_npot_0_gs + 3);

    auto ta_xxyz_0_0 = pbuffer.data(idx_npot_0_gs + 4);

    auto ta_xxzz_0_0 = pbuffer.data(idx_npot_0_gs + 5);

    auto ta_xyyy_0_0 = pbuffer.data(idx_npot_0_gs + 6);

    auto ta_xyyz_0_0 = pbuffer.data(idx_npot_0_gs + 7);

    auto ta_xyzz_0_0 = pbuffer.data(idx_npot_0_gs + 8);

    auto ta_xzzz_0_0 = pbuffer.data(idx_npot_0_gs + 9);

    auto ta_yyyy_0_0 = pbuffer.data(idx_npot_0_gs + 10);

    auto ta_yyyz_0_0 = pbuffer.data(idx_npot_0_gs + 11);

    auto ta_yyzz_0_0 = pbuffer.data(idx_npot_0_gs + 12);

    auto ta_yzzz_0_0 = pbuffer.data(idx_npot_0_gs + 13);

    auto ta_zzzz_0_0 = pbuffer.data(idx_npot_0_gs + 14);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             pa_z,        \
                             pc_x,        \
                             pc_y,        \
                             pc_z,        \
                             ta_xx_0_0,   \
                             ta_xx_0_1,   \
                             ta_xxx_0_0,  \
                             ta_xxx_0_1,  \
                             ta_xxxx_0_0, \
                             ta_xxxy_0_0, \
                             ta_xxxz_0_0, \
                             ta_xxyy_0_0, \
                             ta_xxyz_0_0, \
                             ta_xxz_0_0,  \
                             ta_xxz_0_1,  \
                             ta_xxzz_0_0, \
                             ta_xyy_0_0,  \
                             ta_xyy_0_1,  \
                             ta_xyyy_0_0, \
                             ta_xyyz_0_0, \
                             ta_xyzz_0_0, \
                             ta_xzz_0_0,  \
                             ta_xzz_0_1,  \
                             ta_xzzz_0_0, \
                             ta_yy_0_0,   \
                             ta_yy_0_1,   \
                             ta_yyy_0_0,  \
                             ta_yyy_0_1,  \
                             ta_yyyy_0_0, \
                             ta_yyyz_0_0, \
                             ta_yyz_0_0,  \
                             ta_yyz_0_1,  \
                             ta_yyzz_0_0, \
                             ta_yzz_0_0,  \
                             ta_yzz_0_1,  \
                             ta_yzzz_0_0, \
                             ta_zz_0_0,   \
                             ta_zz_0_1,   \
                             ta_zzz_0_0,  \
                             ta_zzz_0_1,  \
                             ta_zzzz_0_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxx_0_0[i] = 3.0 * ta_xx_0_0[i] * fe_0 - 3.0 * ta_xx_0_1[i] * fe_0 + ta_xxx_0_0[i] * pa_x[i] - ta_xxx_0_1[i] * pc_x[i];

        ta_xxxy_0_0[i] = ta_xxx_0_0[i] * pa_y[i] - ta_xxx_0_1[i] * pc_y[i];

        ta_xxxz_0_0[i] = ta_xxx_0_0[i] * pa_z[i] - ta_xxx_0_1[i] * pc_z[i];

        ta_xxyy_0_0[i] = ta_yy_0_0[i] * fe_0 - ta_yy_0_1[i] * fe_0 + ta_xyy_0_0[i] * pa_x[i] - ta_xyy_0_1[i] * pc_x[i];

        ta_xxyz_0_0[i] = ta_xxz_0_0[i] * pa_y[i] - ta_xxz_0_1[i] * pc_y[i];

        ta_xxzz_0_0[i] = ta_zz_0_0[i] * fe_0 - ta_zz_0_1[i] * fe_0 + ta_xzz_0_0[i] * pa_x[i] - ta_xzz_0_1[i] * pc_x[i];

        ta_xyyy_0_0[i] = ta_yyy_0_0[i] * pa_x[i] - ta_yyy_0_1[i] * pc_x[i];

        ta_xyyz_0_0[i] = ta_yyz_0_0[i] * pa_x[i] - ta_yyz_0_1[i] * pc_x[i];

        ta_xyzz_0_0[i] = ta_yzz_0_0[i] * pa_x[i] - ta_yzz_0_1[i] * pc_x[i];

        ta_xzzz_0_0[i] = ta_zzz_0_0[i] * pa_x[i] - ta_zzz_0_1[i] * pc_x[i];

        ta_yyyy_0_0[i] = 3.0 * ta_yy_0_0[i] * fe_0 - 3.0 * ta_yy_0_1[i] * fe_0 + ta_yyy_0_0[i] * pa_y[i] - ta_yyy_0_1[i] * pc_y[i];

        ta_yyyz_0_0[i] = ta_yyy_0_0[i] * pa_z[i] - ta_yyy_0_1[i] * pc_z[i];

        ta_yyzz_0_0[i] = ta_zz_0_0[i] * fe_0 - ta_zz_0_1[i] * fe_0 + ta_yzz_0_0[i] * pa_y[i] - ta_yzz_0_1[i] * pc_y[i];

        ta_yzzz_0_0[i] = ta_zzz_0_0[i] * pa_y[i] - ta_zzz_0_1[i] * pc_y[i];

        ta_zzzz_0_0[i] = 3.0 * ta_zz_0_0[i] * fe_0 - 3.0 * ta_zz_0_1[i] * fe_0 + ta_zzz_0_0[i] * pa_z[i] - ta_zzz_0_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
