#include "T2CHrrABRecFP.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_fp(CSimdArray<double>& cbuffer, 
            const size_t idx_fp,
            const size_t idx_fs,
            const size_t idx_gs,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : FS

    auto t_xxx_0 = cbuffer.data(idx_fs);

    auto t_xxy_0 = cbuffer.data(idx_fs + 1);

    auto t_xxz_0 = cbuffer.data(idx_fs + 2);

    auto t_xyy_0 = cbuffer.data(idx_fs + 3);

    auto t_xyz_0 = cbuffer.data(idx_fs + 4);

    auto t_xzz_0 = cbuffer.data(idx_fs + 5);

    auto t_yyy_0 = cbuffer.data(idx_fs + 6);

    auto t_yyz_0 = cbuffer.data(idx_fs + 7);

    auto t_yzz_0 = cbuffer.data(idx_fs + 8);

    auto t_zzz_0 = cbuffer.data(idx_fs + 9);

    // Set up components of auxiliary buffer : GS

    auto t_xxxx_0 = cbuffer.data(idx_gs);

    auto t_xxxy_0 = cbuffer.data(idx_gs + 1);

    auto t_xxxz_0 = cbuffer.data(idx_gs + 2);

    auto t_xxyy_0 = cbuffer.data(idx_gs + 3);

    auto t_xxyz_0 = cbuffer.data(idx_gs + 4);

    auto t_xxzz_0 = cbuffer.data(idx_gs + 5);

    auto t_xyyy_0 = cbuffer.data(idx_gs + 6);

    auto t_xyyz_0 = cbuffer.data(idx_gs + 7);

    auto t_xyzz_0 = cbuffer.data(idx_gs + 8);

    auto t_xzzz_0 = cbuffer.data(idx_gs + 9);

    auto t_yyyy_0 = cbuffer.data(idx_gs + 10);

    auto t_yyyz_0 = cbuffer.data(idx_gs + 11);

    auto t_yyzz_0 = cbuffer.data(idx_gs + 12);

    auto t_yzzz_0 = cbuffer.data(idx_gs + 13);

    auto t_zzzz_0 = cbuffer.data(idx_gs + 14);

    // Set up components of targeted buffer : FP

    auto t_xxx_x = cbuffer.data(idx_fp);

    auto t_xxx_y = cbuffer.data(idx_fp + 1);

    auto t_xxx_z = cbuffer.data(idx_fp + 2);

    auto t_xxy_x = cbuffer.data(idx_fp + 3);

    auto t_xxy_y = cbuffer.data(idx_fp + 4);

    auto t_xxy_z = cbuffer.data(idx_fp + 5);

    auto t_xxz_x = cbuffer.data(idx_fp + 6);

    auto t_xxz_y = cbuffer.data(idx_fp + 7);

    auto t_xxz_z = cbuffer.data(idx_fp + 8);

    auto t_xyy_x = cbuffer.data(idx_fp + 9);

    auto t_xyy_y = cbuffer.data(idx_fp + 10);

    auto t_xyy_z = cbuffer.data(idx_fp + 11);

    auto t_xyz_x = cbuffer.data(idx_fp + 12);

    auto t_xyz_y = cbuffer.data(idx_fp + 13);

    auto t_xyz_z = cbuffer.data(idx_fp + 14);

    auto t_xzz_x = cbuffer.data(idx_fp + 15);

    auto t_xzz_y = cbuffer.data(idx_fp + 16);

    auto t_xzz_z = cbuffer.data(idx_fp + 17);

    auto t_yyy_x = cbuffer.data(idx_fp + 18);

    auto t_yyy_y = cbuffer.data(idx_fp + 19);

    auto t_yyy_z = cbuffer.data(idx_fp + 20);

    auto t_yyz_x = cbuffer.data(idx_fp + 21);

    auto t_yyz_y = cbuffer.data(idx_fp + 22);

    auto t_yyz_z = cbuffer.data(idx_fp + 23);

    auto t_yzz_x = cbuffer.data(idx_fp + 24);

    auto t_yzz_y = cbuffer.data(idx_fp + 25);

    auto t_yzz_z = cbuffer.data(idx_fp + 26);

    auto t_zzz_x = cbuffer.data(idx_fp + 27);

    auto t_zzz_y = cbuffer.data(idx_fp + 28);

    auto t_zzz_z = cbuffer.data(idx_fp + 29);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_xxx_0, t_xxx_x, t_xxx_y, t_xxx_z, t_xxxx_0, t_xxxy_0, t_xxxz_0, t_xxy_0, t_xxy_x, t_xxy_y, t_xxy_z, t_xxyy_0, t_xxyz_0, t_xxz_0, t_xxz_x, t_xxz_y, t_xxz_z, t_xxzz_0, t_xyy_0, t_xyy_x, t_xyy_y, t_xyy_z, t_xyyy_0, t_xyyz_0, t_xyz_0, t_xyz_x, t_xyz_y, t_xyz_z, t_xyzz_0, t_xzz_0, t_xzz_x, t_xzz_y, t_xzz_z, t_xzzz_0, t_yyy_0, t_yyy_x, t_yyy_y, t_yyy_z, t_yyyy_0, t_yyyz_0, t_yyz_0, t_yyz_x, t_yyz_y, t_yyz_z, t_yyzz_0, t_yzz_0, t_yzz_x, t_yzz_y, t_yzz_z, t_yzzz_0, t_zzz_0, t_zzz_x, t_zzz_y, t_zzz_z, t_zzzz_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xxx_x[i] = t_xxx_0[i] * ab_x[i] + t_xxxx_0[i];

        t_xxx_y[i] = t_xxx_0[i] * ab_y[i] + t_xxxy_0[i];

        t_xxx_z[i] = t_xxx_0[i] * ab_z[i] + t_xxxz_0[i];

        t_xxy_x[i] = t_xxy_0[i] * ab_x[i] + t_xxxy_0[i];

        t_xxy_y[i] = t_xxy_0[i] * ab_y[i] + t_xxyy_0[i];

        t_xxy_z[i] = t_xxy_0[i] * ab_z[i] + t_xxyz_0[i];

        t_xxz_x[i] = t_xxz_0[i] * ab_x[i] + t_xxxz_0[i];

        t_xxz_y[i] = t_xxz_0[i] * ab_y[i] + t_xxyz_0[i];

        t_xxz_z[i] = t_xxz_0[i] * ab_z[i] + t_xxzz_0[i];

        t_xyy_x[i] = t_xyy_0[i] * ab_x[i] + t_xxyy_0[i];

        t_xyy_y[i] = t_xyy_0[i] * ab_y[i] + t_xyyy_0[i];

        t_xyy_z[i] = t_xyy_0[i] * ab_z[i] + t_xyyz_0[i];

        t_xyz_x[i] = t_xyz_0[i] * ab_x[i] + t_xxyz_0[i];

        t_xyz_y[i] = t_xyz_0[i] * ab_y[i] + t_xyyz_0[i];

        t_xyz_z[i] = t_xyz_0[i] * ab_z[i] + t_xyzz_0[i];

        t_xzz_x[i] = t_xzz_0[i] * ab_x[i] + t_xxzz_0[i];

        t_xzz_y[i] = t_xzz_0[i] * ab_y[i] + t_xyzz_0[i];

        t_xzz_z[i] = t_xzz_0[i] * ab_z[i] + t_xzzz_0[i];

        t_yyy_x[i] = t_yyy_0[i] * ab_x[i] + t_xyyy_0[i];

        t_yyy_y[i] = t_yyy_0[i] * ab_y[i] + t_yyyy_0[i];

        t_yyy_z[i] = t_yyy_0[i] * ab_z[i] + t_yyyz_0[i];

        t_yyz_x[i] = t_yyz_0[i] * ab_x[i] + t_xyyz_0[i];

        t_yyz_y[i] = t_yyz_0[i] * ab_y[i] + t_yyyz_0[i];

        t_yyz_z[i] = t_yyz_0[i] * ab_z[i] + t_yyzz_0[i];

        t_yzz_x[i] = t_yzz_0[i] * ab_x[i] + t_xyzz_0[i];

        t_yzz_y[i] = t_yzz_0[i] * ab_y[i] + t_yyzz_0[i];

        t_yzz_z[i] = t_yzz_0[i] * ab_z[i] + t_yzzz_0[i];

        t_zzz_x[i] = t_zzz_0[i] * ab_x[i] + t_xzzz_0[i];

        t_zzz_y[i] = t_zzz_0[i] * ab_y[i] + t_yzzz_0[i];

        t_zzz_z[i] = t_zzz_0[i] * ab_z[i] + t_zzzz_0[i];
    }
}

} // t2chrr namespace

