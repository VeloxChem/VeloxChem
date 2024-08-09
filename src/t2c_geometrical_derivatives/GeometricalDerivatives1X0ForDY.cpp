#include "GeometricalDerivatives1X0ForDY.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_10_dx(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_100_ds,
                        const size_t idx_op_ps,
                        const size_t idx_op_fs,
                        const size_t op_comps,
                        const size_t ket_comps,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    for (size_t i = 0; i < op_comps; i++)
    {
        for (size_t j = 0; j < ket_comps; j++)
        {
            // Set up components of auxiliary buffer : PS

            auto to_x_0 = pbuffer.data(idx_op_ps + i * 3 * ket_comps + 0 * ket_comps + j);

            auto to_y_0 = pbuffer.data(idx_op_ps + i * 3 * ket_comps + 1 * ket_comps + j);

            auto to_z_0 = pbuffer.data(idx_op_ps + i * 3 * ket_comps + 2 * ket_comps + j);

            // Set up components of auxiliary buffer : FS

            auto to_xxx_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 0 * ket_comps + j);

            auto to_xxy_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 1 * ket_comps + j);

            auto to_xxz_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 2 * ket_comps + j);

            auto to_xyy_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 3 * ket_comps + j);

            auto to_xyz_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 4 * ket_comps + j);

            auto to_xzz_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 5 * ket_comps + j);

            auto to_yyy_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 6 * ket_comps + j);

            auto to_yyz_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 7 * ket_comps + j);

            auto to_yzz_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 8 * ket_comps + j);

            auto to_zzz_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 9 * ket_comps + j);

            // Set up 0-6 components of targeted buffer : DS

            auto to_x_0_xx_0 = pbuffer.data(idx_op_geom_100_ds + 0 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 0 * ket_comps + j);

            auto to_x_0_xy_0 = pbuffer.data(idx_op_geom_100_ds + 0 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 1 * ket_comps + j);

            auto to_x_0_xz_0 = pbuffer.data(idx_op_geom_100_ds + 0 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 2 * ket_comps + j);

            auto to_x_0_yy_0 = pbuffer.data(idx_op_geom_100_ds + 0 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 3 * ket_comps + j);

            auto to_x_0_yz_0 = pbuffer.data(idx_op_geom_100_ds + 0 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 4 * ket_comps + j);

            auto to_x_0_zz_0 = pbuffer.data(idx_op_geom_100_ds + 0 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 5 * ket_comps + j);

            #pragma omp simd aligned(to_x_0, to_x_0_xx_0, to_x_0_xy_0, to_x_0_xz_0, to_x_0_yy_0, to_x_0_yz_0, to_x_0_zz_0, to_xxx_0, to_xxy_0, to_xxz_0, to_xyy_0, to_xyz_0, to_xzz_0, to_y_0, to_z_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_x_0_xx_0[k] = -2.0 * to_x_0[k] + 2.0 * to_xxx_0[k] * tbe_0;

                to_x_0_xy_0[k] = -to_y_0[k] + 2.0 * to_xxy_0[k] * tbe_0;

                to_x_0_xz_0[k] = -to_z_0[k] + 2.0 * to_xxz_0[k] * tbe_0;

                to_x_0_yy_0[k] = 2.0 * to_xyy_0[k] * tbe_0;

                to_x_0_yz_0[k] = 2.0 * to_xyz_0[k] * tbe_0;

                to_x_0_zz_0[k] = 2.0 * to_xzz_0[k] * tbe_0;
            }

            // Set up 6-12 components of targeted buffer : DS

            auto to_y_0_xx_0 = pbuffer.data(idx_op_geom_100_ds + 1 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 0 * ket_comps + j);

            auto to_y_0_xy_0 = pbuffer.data(idx_op_geom_100_ds + 1 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 1 * ket_comps + j);

            auto to_y_0_xz_0 = pbuffer.data(idx_op_geom_100_ds + 1 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 2 * ket_comps + j);

            auto to_y_0_yy_0 = pbuffer.data(idx_op_geom_100_ds + 1 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 3 * ket_comps + j);

            auto to_y_0_yz_0 = pbuffer.data(idx_op_geom_100_ds + 1 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 4 * ket_comps + j);

            auto to_y_0_zz_0 = pbuffer.data(idx_op_geom_100_ds + 1 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 5 * ket_comps + j);

            #pragma omp simd aligned(to_x_0, to_xxy_0, to_xyy_0, to_xyz_0, to_y_0, to_y_0_xx_0, to_y_0_xy_0, to_y_0_xz_0, to_y_0_yy_0, to_y_0_yz_0, to_y_0_zz_0, to_yyy_0, to_yyz_0, to_yzz_0, to_z_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_y_0_xx_0[k] = 2.0 * to_xxy_0[k] * tbe_0;

                to_y_0_xy_0[k] = -to_x_0[k] + 2.0 * to_xyy_0[k] * tbe_0;

                to_y_0_xz_0[k] = 2.0 * to_xyz_0[k] * tbe_0;

                to_y_0_yy_0[k] = -2.0 * to_y_0[k] + 2.0 * to_yyy_0[k] * tbe_0;

                to_y_0_yz_0[k] = -to_z_0[k] + 2.0 * to_yyz_0[k] * tbe_0;

                to_y_0_zz_0[k] = 2.0 * to_yzz_0[k] * tbe_0;
            }

            // Set up 12-18 components of targeted buffer : DS

            auto to_z_0_xx_0 = pbuffer.data(idx_op_geom_100_ds + 2 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 0 * ket_comps + j);

            auto to_z_0_xy_0 = pbuffer.data(idx_op_geom_100_ds + 2 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 1 * ket_comps + j);

            auto to_z_0_xz_0 = pbuffer.data(idx_op_geom_100_ds + 2 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 2 * ket_comps + j);

            auto to_z_0_yy_0 = pbuffer.data(idx_op_geom_100_ds + 2 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 3 * ket_comps + j);

            auto to_z_0_yz_0 = pbuffer.data(idx_op_geom_100_ds + 2 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 4 * ket_comps + j);

            auto to_z_0_zz_0 = pbuffer.data(idx_op_geom_100_ds + 2 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 5 * ket_comps + j);

            #pragma omp simd aligned(to_x_0, to_xxz_0, to_xyz_0, to_xzz_0, to_y_0, to_yyz_0, to_yzz_0, to_z_0, to_z_0_xx_0, to_z_0_xy_0, to_z_0_xz_0, to_z_0_yy_0, to_z_0_yz_0, to_z_0_zz_0, to_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_z_0_xx_0[k] = 2.0 * to_xxz_0[k] * tbe_0;

                to_z_0_xy_0[k] = 2.0 * to_xyz_0[k] * tbe_0;

                to_z_0_xz_0[k] = -to_x_0[k] + 2.0 * to_xzz_0[k] * tbe_0;

                to_z_0_yy_0[k] = 2.0 * to_yyz_0[k] * tbe_0;

                to_z_0_yz_0[k] = -to_y_0[k] + 2.0 * to_yzz_0[k] * tbe_0;

                to_z_0_zz_0[k] = -2.0 * to_z_0[k] + 2.0 * to_zzz_0[k] * tbe_0;
            }

        }
    }

}

} // t2cgeom namespace

