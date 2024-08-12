#include "NuclearPotentialGeom020PrimRecSS.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_ss(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_020_0_ss,
                                        const size_t              idx_npot_1_ss,
                                        const size_t              idx_npot_2_ss,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpc,
                                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : SS

    auto ta_0_0_1 = pbuffer.data(idx_npot_1_ss);

    // Set up components of auxiliary buffer : SS

    auto ta_0_0_2 = pbuffer.data(idx_npot_2_ss);

    // Set up components of auxiliary buffer : SS

    auto ta2_xx_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss);

    auto ta2_xy_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss + 1);

    auto ta2_xz_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss + 2);

    auto ta2_yy_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss + 3);

    auto ta2_yz_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss + 4);

    auto ta2_zz_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss + 5);

#pragma omp simd aligned( \
        pc_x, pc_y, pc_z, ta2_xx_0_0_0, ta2_xy_0_0_0, ta2_xz_0_0_0, ta2_yy_0_0_0, ta2_yz_0_0_0, ta2_zz_0_0_0, ta_0_0_1, ta_0_0_2, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fact = 2.0 * (a_exp + b_exps[i]);

        ta2_xx_0_0_0[i] = fact * fact * pc_x[i] * pc_x[i] * ta_0_0_2[i] - fact * ta_0_0_1[i];

        ta2_xy_0_0_0[i] = fact * fact * pc_x[i] * pc_y[i] * ta_0_0_2[i];

        ta2_xz_0_0_0[i] = fact * fact * pc_x[i] * pc_z[i] * ta_0_0_2[i];

        ta2_yy_0_0_0[i] = fact * fact * pc_y[i] * pc_y[i] * ta_0_0_2[i] - fact * ta_0_0_1[i];

        ta2_yz_0_0_0[i] = fact * fact * pc_y[i] * pc_z[i] * ta_0_0_2[i];

        ta2_zz_0_0_0[i] = fact * fact * pc_z[i] * pc_z[i] * ta_0_0_2[i] - fact * ta_0_0_1[i];
    }
}

}  // namespace npotrec
