#include "NuclearPotentialGeom010PrimRecSS.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_ss(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_ss,
                                        const size_t              idx_npot_1_ss,
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

    auto ta1_x_0_0_0 = pbuffer.data(idx_npot_geom_010_0_ss);

    auto ta1_y_0_0_0 = pbuffer.data(idx_npot_geom_010_0_ss + 1);

    auto ta1_z_0_0_0 = pbuffer.data(idx_npot_geom_010_0_ss + 2);

#pragma omp simd aligned(pc_x, pc_y, pc_z, ta1_x_0_0_0, ta1_y_0_0_0, ta1_z_0_0_0, ta_0_0_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fact = 2.0 * (a_exp + b_exps[i]) * ta_0_0_1[i];

        ta1_x_0_0_0[i] = fact * pc_x[i];

        ta1_y_0_0_0[i] = fact * pc_y[i];

        ta1_z_0_0_0[i] = fact * pc_z[i];
    }
}

}  // namespace npotrec
