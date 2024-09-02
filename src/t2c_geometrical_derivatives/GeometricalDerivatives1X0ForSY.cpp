#include "GeometricalDerivatives1X0ForSY.hpp"

namespace t2cgeom {  // t2cgeom namespace

auto
comp_prim_op_geom_10_sx(CSimdArray<double>& pbuffer,
                        const size_t        idx_op_geom_100_ss,
                        const size_t        idx_op_ps,
                        const size_t        op_comps,
                        const size_t        ket_comps,
                        const double        a_exp) -> void
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

            // Set up components of targeted buffer : SS

            auto to_x_0_0_0 = pbuffer.data(idx_op_geom_100_ss + 0 * op_comps * 1 * ket_comps + i * 1 * ket_comps + 0 * ket_comps + j);

            auto to_y_0_0_0 = pbuffer.data(idx_op_geom_100_ss + 1 * op_comps * 1 * ket_comps + i * 1 * ket_comps + 0 * ket_comps + j);

            auto to_z_0_0_0 = pbuffer.data(idx_op_geom_100_ss + 2 * op_comps * 1 * ket_comps + i * 1 * ket_comps + 0 * ket_comps + j);

#pragma omp simd aligned(to_x_0, to_x_0_0_0, to_y_0, to_y_0_0_0, to_z_0, to_z_0_0_0 : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_x_0_0_0[k] = 2.0 * to_x_0[k] * tbe_0;

                to_y_0_0_0[k] = 2.0 * to_y_0[k] * tbe_0;

                to_z_0_0_0[k] = 2.0 * to_z_0[k] * tbe_0;
            }
        }
    }
}

}  // namespace t2cgeom
