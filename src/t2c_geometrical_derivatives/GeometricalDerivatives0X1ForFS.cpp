#include "GeometricalDerivatives0X1ForFS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_fs(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_fs,
                       const int idx_op_fp,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : FP

        auto to_xxx_x = pbuffer.data(idx_op_fp + i * 30 + 0);

        auto to_xxx_y = pbuffer.data(idx_op_fp + i * 30 + 1);

        auto to_xxx_z = pbuffer.data(idx_op_fp + i * 30 + 2);

        auto to_xxy_x = pbuffer.data(idx_op_fp + i * 30 + 3);

        auto to_xxy_y = pbuffer.data(idx_op_fp + i * 30 + 4);

        auto to_xxy_z = pbuffer.data(idx_op_fp + i * 30 + 5);

        auto to_xxz_x = pbuffer.data(idx_op_fp + i * 30 + 6);

        auto to_xxz_y = pbuffer.data(idx_op_fp + i * 30 + 7);

        auto to_xxz_z = pbuffer.data(idx_op_fp + i * 30 + 8);

        auto to_xyy_x = pbuffer.data(idx_op_fp + i * 30 + 9);

        auto to_xyy_y = pbuffer.data(idx_op_fp + i * 30 + 10);

        auto to_xyy_z = pbuffer.data(idx_op_fp + i * 30 + 11);

        auto to_xyz_x = pbuffer.data(idx_op_fp + i * 30 + 12);

        auto to_xyz_y = pbuffer.data(idx_op_fp + i * 30 + 13);

        auto to_xyz_z = pbuffer.data(idx_op_fp + i * 30 + 14);

        auto to_xzz_x = pbuffer.data(idx_op_fp + i * 30 + 15);

        auto to_xzz_y = pbuffer.data(idx_op_fp + i * 30 + 16);

        auto to_xzz_z = pbuffer.data(idx_op_fp + i * 30 + 17);

        auto to_yyy_x = pbuffer.data(idx_op_fp + i * 30 + 18);

        auto to_yyy_y = pbuffer.data(idx_op_fp + i * 30 + 19);

        auto to_yyy_z = pbuffer.data(idx_op_fp + i * 30 + 20);

        auto to_yyz_x = pbuffer.data(idx_op_fp + i * 30 + 21);

        auto to_yyz_y = pbuffer.data(idx_op_fp + i * 30 + 22);

        auto to_yyz_z = pbuffer.data(idx_op_fp + i * 30 + 23);

        auto to_yzz_x = pbuffer.data(idx_op_fp + i * 30 + 24);

        auto to_yzz_y = pbuffer.data(idx_op_fp + i * 30 + 25);

        auto to_yzz_z = pbuffer.data(idx_op_fp + i * 30 + 26);

        auto to_zzz_x = pbuffer.data(idx_op_fp + i * 30 + 27);

        auto to_zzz_y = pbuffer.data(idx_op_fp + i * 30 + 28);

        auto to_zzz_z = pbuffer.data(idx_op_fp + i * 30 + 29);

        // Set up 0-10 components of targeted buffer : FS

        auto to_0_x_xxx_0 = pbuffer.data(idx_op_geom_001_fs + 0 * op_comps * 10 + i * 10 + 0);

        auto to_0_x_xxy_0 = pbuffer.data(idx_op_geom_001_fs + 0 * op_comps * 10 + i * 10 + 1);

        auto to_0_x_xxz_0 = pbuffer.data(idx_op_geom_001_fs + 0 * op_comps * 10 + i * 10 + 2);

        auto to_0_x_xyy_0 = pbuffer.data(idx_op_geom_001_fs + 0 * op_comps * 10 + i * 10 + 3);

        auto to_0_x_xyz_0 = pbuffer.data(idx_op_geom_001_fs + 0 * op_comps * 10 + i * 10 + 4);

        auto to_0_x_xzz_0 = pbuffer.data(idx_op_geom_001_fs + 0 * op_comps * 10 + i * 10 + 5);

        auto to_0_x_yyy_0 = pbuffer.data(idx_op_geom_001_fs + 0 * op_comps * 10 + i * 10 + 6);

        auto to_0_x_yyz_0 = pbuffer.data(idx_op_geom_001_fs + 0 * op_comps * 10 + i * 10 + 7);

        auto to_0_x_yzz_0 = pbuffer.data(idx_op_geom_001_fs + 0 * op_comps * 10 + i * 10 + 8);

        auto to_0_x_zzz_0 = pbuffer.data(idx_op_geom_001_fs + 0 * op_comps * 10 + i * 10 + 9);

        #pragma omp simd aligned(to_0_x_xxx_0, to_0_x_xxy_0, to_0_x_xxz_0, to_0_x_xyy_0, to_0_x_xyz_0, to_0_x_xzz_0, to_0_x_yyy_0, to_0_x_yyz_0, to_0_x_yzz_0, to_0_x_zzz_0, to_xxx_x, to_xxy_x, to_xxz_x, to_xyy_x, to_xyz_x, to_xzz_x, to_yyy_x, to_yyz_x, to_yzz_x, to_zzz_x, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxx_0[k] = 2.0 * to_xxx_x[k] * tke_0;

            to_0_x_xxy_0[k] = 2.0 * to_xxy_x[k] * tke_0;

            to_0_x_xxz_0[k] = 2.0 * to_xxz_x[k] * tke_0;

            to_0_x_xyy_0[k] = 2.0 * to_xyy_x[k] * tke_0;

            to_0_x_xyz_0[k] = 2.0 * to_xyz_x[k] * tke_0;

            to_0_x_xzz_0[k] = 2.0 * to_xzz_x[k] * tke_0;

            to_0_x_yyy_0[k] = 2.0 * to_yyy_x[k] * tke_0;

            to_0_x_yyz_0[k] = 2.0 * to_yyz_x[k] * tke_0;

            to_0_x_yzz_0[k] = 2.0 * to_yzz_x[k] * tke_0;

            to_0_x_zzz_0[k] = 2.0 * to_zzz_x[k] * tke_0;
        }

        // Set up 10-20 components of targeted buffer : FS

        auto to_0_y_xxx_0 = pbuffer.data(idx_op_geom_001_fs + 1 * op_comps * 10 + i * 10 + 0);

        auto to_0_y_xxy_0 = pbuffer.data(idx_op_geom_001_fs + 1 * op_comps * 10 + i * 10 + 1);

        auto to_0_y_xxz_0 = pbuffer.data(idx_op_geom_001_fs + 1 * op_comps * 10 + i * 10 + 2);

        auto to_0_y_xyy_0 = pbuffer.data(idx_op_geom_001_fs + 1 * op_comps * 10 + i * 10 + 3);

        auto to_0_y_xyz_0 = pbuffer.data(idx_op_geom_001_fs + 1 * op_comps * 10 + i * 10 + 4);

        auto to_0_y_xzz_0 = pbuffer.data(idx_op_geom_001_fs + 1 * op_comps * 10 + i * 10 + 5);

        auto to_0_y_yyy_0 = pbuffer.data(idx_op_geom_001_fs + 1 * op_comps * 10 + i * 10 + 6);

        auto to_0_y_yyz_0 = pbuffer.data(idx_op_geom_001_fs + 1 * op_comps * 10 + i * 10 + 7);

        auto to_0_y_yzz_0 = pbuffer.data(idx_op_geom_001_fs + 1 * op_comps * 10 + i * 10 + 8);

        auto to_0_y_zzz_0 = pbuffer.data(idx_op_geom_001_fs + 1 * op_comps * 10 + i * 10 + 9);

        #pragma omp simd aligned(to_0_y_xxx_0, to_0_y_xxy_0, to_0_y_xxz_0, to_0_y_xyy_0, to_0_y_xyz_0, to_0_y_xzz_0, to_0_y_yyy_0, to_0_y_yyz_0, to_0_y_yzz_0, to_0_y_zzz_0, to_xxx_y, to_xxy_y, to_xxz_y, to_xyy_y, to_xyz_y, to_xzz_y, to_yyy_y, to_yyz_y, to_yzz_y, to_zzz_y, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxx_0[k] = 2.0 * to_xxx_y[k] * tke_0;

            to_0_y_xxy_0[k] = 2.0 * to_xxy_y[k] * tke_0;

            to_0_y_xxz_0[k] = 2.0 * to_xxz_y[k] * tke_0;

            to_0_y_xyy_0[k] = 2.0 * to_xyy_y[k] * tke_0;

            to_0_y_xyz_0[k] = 2.0 * to_xyz_y[k] * tke_0;

            to_0_y_xzz_0[k] = 2.0 * to_xzz_y[k] * tke_0;

            to_0_y_yyy_0[k] = 2.0 * to_yyy_y[k] * tke_0;

            to_0_y_yyz_0[k] = 2.0 * to_yyz_y[k] * tke_0;

            to_0_y_yzz_0[k] = 2.0 * to_yzz_y[k] * tke_0;

            to_0_y_zzz_0[k] = 2.0 * to_zzz_y[k] * tke_0;
        }

        // Set up 20-30 components of targeted buffer : FS

        auto to_0_z_xxx_0 = pbuffer.data(idx_op_geom_001_fs + 2 * op_comps * 10 + i * 10 + 0);

        auto to_0_z_xxy_0 = pbuffer.data(idx_op_geom_001_fs + 2 * op_comps * 10 + i * 10 + 1);

        auto to_0_z_xxz_0 = pbuffer.data(idx_op_geom_001_fs + 2 * op_comps * 10 + i * 10 + 2);

        auto to_0_z_xyy_0 = pbuffer.data(idx_op_geom_001_fs + 2 * op_comps * 10 + i * 10 + 3);

        auto to_0_z_xyz_0 = pbuffer.data(idx_op_geom_001_fs + 2 * op_comps * 10 + i * 10 + 4);

        auto to_0_z_xzz_0 = pbuffer.data(idx_op_geom_001_fs + 2 * op_comps * 10 + i * 10 + 5);

        auto to_0_z_yyy_0 = pbuffer.data(idx_op_geom_001_fs + 2 * op_comps * 10 + i * 10 + 6);

        auto to_0_z_yyz_0 = pbuffer.data(idx_op_geom_001_fs + 2 * op_comps * 10 + i * 10 + 7);

        auto to_0_z_yzz_0 = pbuffer.data(idx_op_geom_001_fs + 2 * op_comps * 10 + i * 10 + 8);

        auto to_0_z_zzz_0 = pbuffer.data(idx_op_geom_001_fs + 2 * op_comps * 10 + i * 10 + 9);

        #pragma omp simd aligned(to_0_z_xxx_0, to_0_z_xxy_0, to_0_z_xxz_0, to_0_z_xyy_0, to_0_z_xyz_0, to_0_z_xzz_0, to_0_z_yyy_0, to_0_z_yyz_0, to_0_z_yzz_0, to_0_z_zzz_0, to_xxx_z, to_xxy_z, to_xxz_z, to_xyy_z, to_xyz_z, to_xzz_z, to_yyy_z, to_yyz_z, to_yzz_z, to_zzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxx_0[k] = 2.0 * to_xxx_z[k] * tke_0;

            to_0_z_xxy_0[k] = 2.0 * to_xxy_z[k] * tke_0;

            to_0_z_xxz_0[k] = 2.0 * to_xxz_z[k] * tke_0;

            to_0_z_xyy_0[k] = 2.0 * to_xyy_z[k] * tke_0;

            to_0_z_xyz_0[k] = 2.0 * to_xyz_z[k] * tke_0;

            to_0_z_xzz_0[k] = 2.0 * to_xzz_z[k] * tke_0;

            to_0_z_yyy_0[k] = 2.0 * to_yyy_z[k] * tke_0;

            to_0_z_yyz_0[k] = 2.0 * to_yyz_z[k] * tke_0;

            to_0_z_yzz_0[k] = 2.0 * to_yzz_z[k] * tke_0;

            to_0_z_zzz_0[k] = 2.0 * to_zzz_z[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

