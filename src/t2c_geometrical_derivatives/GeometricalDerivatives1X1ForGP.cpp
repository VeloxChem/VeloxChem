#include "GeometricalDerivatives1X1ForGP.hpp"

namespace t2cgeom {  // t2cgeom namespace

auto
comp_prim_op_geom_11_gp(CSimdArray<double>&       pbuffer,
                        const size_t              idx_op_geom_101_gp,
                        const size_t              idx_op_fs,
                        const size_t              idx_op_fd,
                        const size_t              idx_op_hs,
                        const size_t              idx_op_hd,
                        const size_t              op_comps,
                        const CSimdArray<double>& factors,
                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : FS

        auto to_xxx_0 = pbuffer.data(idx_op_fs + i * 10 + 0);

        auto to_xxy_0 = pbuffer.data(idx_op_fs + i * 10 + 1);

        auto to_xxz_0 = pbuffer.data(idx_op_fs + i * 10 + 2);

        auto to_xyy_0 = pbuffer.data(idx_op_fs + i * 10 + 3);

        auto to_xyz_0 = pbuffer.data(idx_op_fs + i * 10 + 4);

        auto to_xzz_0 = pbuffer.data(idx_op_fs + i * 10 + 5);

        auto to_yyy_0 = pbuffer.data(idx_op_fs + i * 10 + 6);

        auto to_yyz_0 = pbuffer.data(idx_op_fs + i * 10 + 7);

        auto to_yzz_0 = pbuffer.data(idx_op_fs + i * 10 + 8);

        auto to_zzz_0 = pbuffer.data(idx_op_fs + i * 10 + 9);

        // Set up components of auxiliary buffer : FD

        auto to_xxx_xx = pbuffer.data(idx_op_fd + i * 60 + 0);

        auto to_xxx_xy = pbuffer.data(idx_op_fd + i * 60 + 1);

        auto to_xxx_xz = pbuffer.data(idx_op_fd + i * 60 + 2);

        auto to_xxx_yy = pbuffer.data(idx_op_fd + i * 60 + 3);

        auto to_xxx_yz = pbuffer.data(idx_op_fd + i * 60 + 4);

        auto to_xxx_zz = pbuffer.data(idx_op_fd + i * 60 + 5);

        auto to_xxy_xx = pbuffer.data(idx_op_fd + i * 60 + 6);

        auto to_xxy_xy = pbuffer.data(idx_op_fd + i * 60 + 7);

        auto to_xxy_xz = pbuffer.data(idx_op_fd + i * 60 + 8);

        auto to_xxy_yy = pbuffer.data(idx_op_fd + i * 60 + 9);

        auto to_xxy_yz = pbuffer.data(idx_op_fd + i * 60 + 10);

        auto to_xxy_zz = pbuffer.data(idx_op_fd + i * 60 + 11);

        auto to_xxz_xx = pbuffer.data(idx_op_fd + i * 60 + 12);

        auto to_xxz_xy = pbuffer.data(idx_op_fd + i * 60 + 13);

        auto to_xxz_xz = pbuffer.data(idx_op_fd + i * 60 + 14);

        auto to_xxz_yy = pbuffer.data(idx_op_fd + i * 60 + 15);

        auto to_xxz_yz = pbuffer.data(idx_op_fd + i * 60 + 16);

        auto to_xxz_zz = pbuffer.data(idx_op_fd + i * 60 + 17);

        auto to_xyy_xx = pbuffer.data(idx_op_fd + i * 60 + 18);

        auto to_xyy_xy = pbuffer.data(idx_op_fd + i * 60 + 19);

        auto to_xyy_xz = pbuffer.data(idx_op_fd + i * 60 + 20);

        auto to_xyy_yy = pbuffer.data(idx_op_fd + i * 60 + 21);

        auto to_xyy_yz = pbuffer.data(idx_op_fd + i * 60 + 22);

        auto to_xyy_zz = pbuffer.data(idx_op_fd + i * 60 + 23);

        auto to_xyz_xx = pbuffer.data(idx_op_fd + i * 60 + 24);

        auto to_xyz_xy = pbuffer.data(idx_op_fd + i * 60 + 25);

        auto to_xyz_xz = pbuffer.data(idx_op_fd + i * 60 + 26);

        auto to_xyz_yy = pbuffer.data(idx_op_fd + i * 60 + 27);

        auto to_xyz_yz = pbuffer.data(idx_op_fd + i * 60 + 28);

        auto to_xyz_zz = pbuffer.data(idx_op_fd + i * 60 + 29);

        auto to_xzz_xx = pbuffer.data(idx_op_fd + i * 60 + 30);

        auto to_xzz_xy = pbuffer.data(idx_op_fd + i * 60 + 31);

        auto to_xzz_xz = pbuffer.data(idx_op_fd + i * 60 + 32);

        auto to_xzz_yy = pbuffer.data(idx_op_fd + i * 60 + 33);

        auto to_xzz_yz = pbuffer.data(idx_op_fd + i * 60 + 34);

        auto to_xzz_zz = pbuffer.data(idx_op_fd + i * 60 + 35);

        auto to_yyy_xx = pbuffer.data(idx_op_fd + i * 60 + 36);

        auto to_yyy_xy = pbuffer.data(idx_op_fd + i * 60 + 37);

        auto to_yyy_xz = pbuffer.data(idx_op_fd + i * 60 + 38);

        auto to_yyy_yy = pbuffer.data(idx_op_fd + i * 60 + 39);

        auto to_yyy_yz = pbuffer.data(idx_op_fd + i * 60 + 40);

        auto to_yyy_zz = pbuffer.data(idx_op_fd + i * 60 + 41);

        auto to_yyz_xx = pbuffer.data(idx_op_fd + i * 60 + 42);

        auto to_yyz_xy = pbuffer.data(idx_op_fd + i * 60 + 43);

        auto to_yyz_xz = pbuffer.data(idx_op_fd + i * 60 + 44);

        auto to_yyz_yy = pbuffer.data(idx_op_fd + i * 60 + 45);

        auto to_yyz_yz = pbuffer.data(idx_op_fd + i * 60 + 46);

        auto to_yyz_zz = pbuffer.data(idx_op_fd + i * 60 + 47);

        auto to_yzz_xx = pbuffer.data(idx_op_fd + i * 60 + 48);

        auto to_yzz_xy = pbuffer.data(idx_op_fd + i * 60 + 49);

        auto to_yzz_xz = pbuffer.data(idx_op_fd + i * 60 + 50);

        auto to_yzz_yy = pbuffer.data(idx_op_fd + i * 60 + 51);

        auto to_yzz_yz = pbuffer.data(idx_op_fd + i * 60 + 52);

        auto to_yzz_zz = pbuffer.data(idx_op_fd + i * 60 + 53);

        auto to_zzz_xx = pbuffer.data(idx_op_fd + i * 60 + 54);

        auto to_zzz_xy = pbuffer.data(idx_op_fd + i * 60 + 55);

        auto to_zzz_xz = pbuffer.data(idx_op_fd + i * 60 + 56);

        auto to_zzz_yy = pbuffer.data(idx_op_fd + i * 60 + 57);

        auto to_zzz_yz = pbuffer.data(idx_op_fd + i * 60 + 58);

        auto to_zzz_zz = pbuffer.data(idx_op_fd + i * 60 + 59);

        // Set up components of auxiliary buffer : HS

        auto to_xxxxx_0 = pbuffer.data(idx_op_hs + i * 21 + 0);

        auto to_xxxxy_0 = pbuffer.data(idx_op_hs + i * 21 + 1);

        auto to_xxxxz_0 = pbuffer.data(idx_op_hs + i * 21 + 2);

        auto to_xxxyy_0 = pbuffer.data(idx_op_hs + i * 21 + 3);

        auto to_xxxyz_0 = pbuffer.data(idx_op_hs + i * 21 + 4);

        auto to_xxxzz_0 = pbuffer.data(idx_op_hs + i * 21 + 5);

        auto to_xxyyy_0 = pbuffer.data(idx_op_hs + i * 21 + 6);

        auto to_xxyyz_0 = pbuffer.data(idx_op_hs + i * 21 + 7);

        auto to_xxyzz_0 = pbuffer.data(idx_op_hs + i * 21 + 8);

        auto to_xxzzz_0 = pbuffer.data(idx_op_hs + i * 21 + 9);

        auto to_xyyyy_0 = pbuffer.data(idx_op_hs + i * 21 + 10);

        auto to_xyyyz_0 = pbuffer.data(idx_op_hs + i * 21 + 11);

        auto to_xyyzz_0 = pbuffer.data(idx_op_hs + i * 21 + 12);

        auto to_xyzzz_0 = pbuffer.data(idx_op_hs + i * 21 + 13);

        auto to_xzzzz_0 = pbuffer.data(idx_op_hs + i * 21 + 14);

        auto to_yyyyy_0 = pbuffer.data(idx_op_hs + i * 21 + 15);

        auto to_yyyyz_0 = pbuffer.data(idx_op_hs + i * 21 + 16);

        auto to_yyyzz_0 = pbuffer.data(idx_op_hs + i * 21 + 17);

        auto to_yyzzz_0 = pbuffer.data(idx_op_hs + i * 21 + 18);

        auto to_yzzzz_0 = pbuffer.data(idx_op_hs + i * 21 + 19);

        auto to_zzzzz_0 = pbuffer.data(idx_op_hs + i * 21 + 20);

        // Set up components of auxiliary buffer : HD

        auto to_xxxxx_xx = pbuffer.data(idx_op_hd + i * 126 + 0);

        auto to_xxxxx_xy = pbuffer.data(idx_op_hd + i * 126 + 1);

        auto to_xxxxx_xz = pbuffer.data(idx_op_hd + i * 126 + 2);

        auto to_xxxxx_yy = pbuffer.data(idx_op_hd + i * 126 + 3);

        auto to_xxxxx_yz = pbuffer.data(idx_op_hd + i * 126 + 4);

        auto to_xxxxx_zz = pbuffer.data(idx_op_hd + i * 126 + 5);

        auto to_xxxxy_xx = pbuffer.data(idx_op_hd + i * 126 + 6);

        auto to_xxxxy_xy = pbuffer.data(idx_op_hd + i * 126 + 7);

        auto to_xxxxy_xz = pbuffer.data(idx_op_hd + i * 126 + 8);

        auto to_xxxxy_yy = pbuffer.data(idx_op_hd + i * 126 + 9);

        auto to_xxxxy_yz = pbuffer.data(idx_op_hd + i * 126 + 10);

        auto to_xxxxy_zz = pbuffer.data(idx_op_hd + i * 126 + 11);

        auto to_xxxxz_xx = pbuffer.data(idx_op_hd + i * 126 + 12);

        auto to_xxxxz_xy = pbuffer.data(idx_op_hd + i * 126 + 13);

        auto to_xxxxz_xz = pbuffer.data(idx_op_hd + i * 126 + 14);

        auto to_xxxxz_yy = pbuffer.data(idx_op_hd + i * 126 + 15);

        auto to_xxxxz_yz = pbuffer.data(idx_op_hd + i * 126 + 16);

        auto to_xxxxz_zz = pbuffer.data(idx_op_hd + i * 126 + 17);

        auto to_xxxyy_xx = pbuffer.data(idx_op_hd + i * 126 + 18);

        auto to_xxxyy_xy = pbuffer.data(idx_op_hd + i * 126 + 19);

        auto to_xxxyy_xz = pbuffer.data(idx_op_hd + i * 126 + 20);

        auto to_xxxyy_yy = pbuffer.data(idx_op_hd + i * 126 + 21);

        auto to_xxxyy_yz = pbuffer.data(idx_op_hd + i * 126 + 22);

        auto to_xxxyy_zz = pbuffer.data(idx_op_hd + i * 126 + 23);

        auto to_xxxyz_xx = pbuffer.data(idx_op_hd + i * 126 + 24);

        auto to_xxxyz_xy = pbuffer.data(idx_op_hd + i * 126 + 25);

        auto to_xxxyz_xz = pbuffer.data(idx_op_hd + i * 126 + 26);

        auto to_xxxyz_yy = pbuffer.data(idx_op_hd + i * 126 + 27);

        auto to_xxxyz_yz = pbuffer.data(idx_op_hd + i * 126 + 28);

        auto to_xxxyz_zz = pbuffer.data(idx_op_hd + i * 126 + 29);

        auto to_xxxzz_xx = pbuffer.data(idx_op_hd + i * 126 + 30);

        auto to_xxxzz_xy = pbuffer.data(idx_op_hd + i * 126 + 31);

        auto to_xxxzz_xz = pbuffer.data(idx_op_hd + i * 126 + 32);

        auto to_xxxzz_yy = pbuffer.data(idx_op_hd + i * 126 + 33);

        auto to_xxxzz_yz = pbuffer.data(idx_op_hd + i * 126 + 34);

        auto to_xxxzz_zz = pbuffer.data(idx_op_hd + i * 126 + 35);

        auto to_xxyyy_xx = pbuffer.data(idx_op_hd + i * 126 + 36);

        auto to_xxyyy_xy = pbuffer.data(idx_op_hd + i * 126 + 37);

        auto to_xxyyy_xz = pbuffer.data(idx_op_hd + i * 126 + 38);

        auto to_xxyyy_yy = pbuffer.data(idx_op_hd + i * 126 + 39);

        auto to_xxyyy_yz = pbuffer.data(idx_op_hd + i * 126 + 40);

        auto to_xxyyy_zz = pbuffer.data(idx_op_hd + i * 126 + 41);

        auto to_xxyyz_xx = pbuffer.data(idx_op_hd + i * 126 + 42);

        auto to_xxyyz_xy = pbuffer.data(idx_op_hd + i * 126 + 43);

        auto to_xxyyz_xz = pbuffer.data(idx_op_hd + i * 126 + 44);

        auto to_xxyyz_yy = pbuffer.data(idx_op_hd + i * 126 + 45);

        auto to_xxyyz_yz = pbuffer.data(idx_op_hd + i * 126 + 46);

        auto to_xxyyz_zz = pbuffer.data(idx_op_hd + i * 126 + 47);

        auto to_xxyzz_xx = pbuffer.data(idx_op_hd + i * 126 + 48);

        auto to_xxyzz_xy = pbuffer.data(idx_op_hd + i * 126 + 49);

        auto to_xxyzz_xz = pbuffer.data(idx_op_hd + i * 126 + 50);

        auto to_xxyzz_yy = pbuffer.data(idx_op_hd + i * 126 + 51);

        auto to_xxyzz_yz = pbuffer.data(idx_op_hd + i * 126 + 52);

        auto to_xxyzz_zz = pbuffer.data(idx_op_hd + i * 126 + 53);

        auto to_xxzzz_xx = pbuffer.data(idx_op_hd + i * 126 + 54);

        auto to_xxzzz_xy = pbuffer.data(idx_op_hd + i * 126 + 55);

        auto to_xxzzz_xz = pbuffer.data(idx_op_hd + i * 126 + 56);

        auto to_xxzzz_yy = pbuffer.data(idx_op_hd + i * 126 + 57);

        auto to_xxzzz_yz = pbuffer.data(idx_op_hd + i * 126 + 58);

        auto to_xxzzz_zz = pbuffer.data(idx_op_hd + i * 126 + 59);

        auto to_xyyyy_xx = pbuffer.data(idx_op_hd + i * 126 + 60);

        auto to_xyyyy_xy = pbuffer.data(idx_op_hd + i * 126 + 61);

        auto to_xyyyy_xz = pbuffer.data(idx_op_hd + i * 126 + 62);

        auto to_xyyyy_yy = pbuffer.data(idx_op_hd + i * 126 + 63);

        auto to_xyyyy_yz = pbuffer.data(idx_op_hd + i * 126 + 64);

        auto to_xyyyy_zz = pbuffer.data(idx_op_hd + i * 126 + 65);

        auto to_xyyyz_xx = pbuffer.data(idx_op_hd + i * 126 + 66);

        auto to_xyyyz_xy = pbuffer.data(idx_op_hd + i * 126 + 67);

        auto to_xyyyz_xz = pbuffer.data(idx_op_hd + i * 126 + 68);

        auto to_xyyyz_yy = pbuffer.data(idx_op_hd + i * 126 + 69);

        auto to_xyyyz_yz = pbuffer.data(idx_op_hd + i * 126 + 70);

        auto to_xyyyz_zz = pbuffer.data(idx_op_hd + i * 126 + 71);

        auto to_xyyzz_xx = pbuffer.data(idx_op_hd + i * 126 + 72);

        auto to_xyyzz_xy = pbuffer.data(idx_op_hd + i * 126 + 73);

        auto to_xyyzz_xz = pbuffer.data(idx_op_hd + i * 126 + 74);

        auto to_xyyzz_yy = pbuffer.data(idx_op_hd + i * 126 + 75);

        auto to_xyyzz_yz = pbuffer.data(idx_op_hd + i * 126 + 76);

        auto to_xyyzz_zz = pbuffer.data(idx_op_hd + i * 126 + 77);

        auto to_xyzzz_xx = pbuffer.data(idx_op_hd + i * 126 + 78);

        auto to_xyzzz_xy = pbuffer.data(idx_op_hd + i * 126 + 79);

        auto to_xyzzz_xz = pbuffer.data(idx_op_hd + i * 126 + 80);

        auto to_xyzzz_yy = pbuffer.data(idx_op_hd + i * 126 + 81);

        auto to_xyzzz_yz = pbuffer.data(idx_op_hd + i * 126 + 82);

        auto to_xyzzz_zz = pbuffer.data(idx_op_hd + i * 126 + 83);

        auto to_xzzzz_xx = pbuffer.data(idx_op_hd + i * 126 + 84);

        auto to_xzzzz_xy = pbuffer.data(idx_op_hd + i * 126 + 85);

        auto to_xzzzz_xz = pbuffer.data(idx_op_hd + i * 126 + 86);

        auto to_xzzzz_yy = pbuffer.data(idx_op_hd + i * 126 + 87);

        auto to_xzzzz_yz = pbuffer.data(idx_op_hd + i * 126 + 88);

        auto to_xzzzz_zz = pbuffer.data(idx_op_hd + i * 126 + 89);

        auto to_yyyyy_xx = pbuffer.data(idx_op_hd + i * 126 + 90);

        auto to_yyyyy_xy = pbuffer.data(idx_op_hd + i * 126 + 91);

        auto to_yyyyy_xz = pbuffer.data(idx_op_hd + i * 126 + 92);

        auto to_yyyyy_yy = pbuffer.data(idx_op_hd + i * 126 + 93);

        auto to_yyyyy_yz = pbuffer.data(idx_op_hd + i * 126 + 94);

        auto to_yyyyy_zz = pbuffer.data(idx_op_hd + i * 126 + 95);

        auto to_yyyyz_xx = pbuffer.data(idx_op_hd + i * 126 + 96);

        auto to_yyyyz_xy = pbuffer.data(idx_op_hd + i * 126 + 97);

        auto to_yyyyz_xz = pbuffer.data(idx_op_hd + i * 126 + 98);

        auto to_yyyyz_yy = pbuffer.data(idx_op_hd + i * 126 + 99);

        auto to_yyyyz_yz = pbuffer.data(idx_op_hd + i * 126 + 100);

        auto to_yyyyz_zz = pbuffer.data(idx_op_hd + i * 126 + 101);

        auto to_yyyzz_xx = pbuffer.data(idx_op_hd + i * 126 + 102);

        auto to_yyyzz_xy = pbuffer.data(idx_op_hd + i * 126 + 103);

        auto to_yyyzz_xz = pbuffer.data(idx_op_hd + i * 126 + 104);

        auto to_yyyzz_yy = pbuffer.data(idx_op_hd + i * 126 + 105);

        auto to_yyyzz_yz = pbuffer.data(idx_op_hd + i * 126 + 106);

        auto to_yyyzz_zz = pbuffer.data(idx_op_hd + i * 126 + 107);

        auto to_yyzzz_xx = pbuffer.data(idx_op_hd + i * 126 + 108);

        auto to_yyzzz_xy = pbuffer.data(idx_op_hd + i * 126 + 109);

        auto to_yyzzz_xz = pbuffer.data(idx_op_hd + i * 126 + 110);

        auto to_yyzzz_yy = pbuffer.data(idx_op_hd + i * 126 + 111);

        auto to_yyzzz_yz = pbuffer.data(idx_op_hd + i * 126 + 112);

        auto to_yyzzz_zz = pbuffer.data(idx_op_hd + i * 126 + 113);

        auto to_yzzzz_xx = pbuffer.data(idx_op_hd + i * 126 + 114);

        auto to_yzzzz_xy = pbuffer.data(idx_op_hd + i * 126 + 115);

        auto to_yzzzz_xz = pbuffer.data(idx_op_hd + i * 126 + 116);

        auto to_yzzzz_yy = pbuffer.data(idx_op_hd + i * 126 + 117);

        auto to_yzzzz_yz = pbuffer.data(idx_op_hd + i * 126 + 118);

        auto to_yzzzz_zz = pbuffer.data(idx_op_hd + i * 126 + 119);

        auto to_zzzzz_xx = pbuffer.data(idx_op_hd + i * 126 + 120);

        auto to_zzzzz_xy = pbuffer.data(idx_op_hd + i * 126 + 121);

        auto to_zzzzz_xz = pbuffer.data(idx_op_hd + i * 126 + 122);

        auto to_zzzzz_yy = pbuffer.data(idx_op_hd + i * 126 + 123);

        auto to_zzzzz_yz = pbuffer.data(idx_op_hd + i * 126 + 124);

        auto to_zzzzz_zz = pbuffer.data(idx_op_hd + i * 126 + 125);

        // Set up 0-3 components of targeted buffer : GP

        auto to_x_x_xxxx_x = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 0);

        auto to_x_x_xxxx_y = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 1);

        auto to_x_x_xxxx_z = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 2);

#pragma omp simd aligned(to_x_x_xxxx_x,     \
                             to_x_x_xxxx_y, \
                             to_x_x_xxxx_z, \
                             to_xxx_0,      \
                             to_xxx_xx,     \
                             to_xxx_xy,     \
                             to_xxx_xz,     \
                             to_xxxxx_0,    \
                             to_xxxxx_xx,   \
                             to_xxxxx_xy,   \
                             to_xxxxx_xz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxxx_x[k] = 4.0 * to_xxx_0[k] - 8.0 * to_xxx_xx[k] * tke_0 - 2.0 * to_xxxxx_0[k] * tbe_0 + 4.0 * to_xxxxx_xx[k] * tbe_0 * tke_0;

            to_x_x_xxxx_y[k] = -8.0 * to_xxx_xy[k] * tke_0 + 4.0 * to_xxxxx_xy[k] * tbe_0 * tke_0;

            to_x_x_xxxx_z[k] = -8.0 * to_xxx_xz[k] * tke_0 + 4.0 * to_xxxxx_xz[k] * tbe_0 * tke_0;
        }

        // Set up 3-6 components of targeted buffer : GP

        auto to_x_x_xxxy_x = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 3);

        auto to_x_x_xxxy_y = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 4);

        auto to_x_x_xxxy_z = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 5);

#pragma omp simd aligned(to_x_x_xxxy_x,     \
                             to_x_x_xxxy_y, \
                             to_x_x_xxxy_z, \
                             to_xxxxy_0,    \
                             to_xxxxy_xx,   \
                             to_xxxxy_xy,   \
                             to_xxxxy_xz,   \
                             to_xxy_0,      \
                             to_xxy_xx,     \
                             to_xxy_xy,     \
                             to_xxy_xz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxxy_x[k] = 3.0 * to_xxy_0[k] - 6.0 * to_xxy_xx[k] * tke_0 - 2.0 * to_xxxxy_0[k] * tbe_0 + 4.0 * to_xxxxy_xx[k] * tbe_0 * tke_0;

            to_x_x_xxxy_y[k] = -6.0 * to_xxy_xy[k] * tke_0 + 4.0 * to_xxxxy_xy[k] * tbe_0 * tke_0;

            to_x_x_xxxy_z[k] = -6.0 * to_xxy_xz[k] * tke_0 + 4.0 * to_xxxxy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 6-9 components of targeted buffer : GP

        auto to_x_x_xxxz_x = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 6);

        auto to_x_x_xxxz_y = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 7);

        auto to_x_x_xxxz_z = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 8);

#pragma omp simd aligned(to_x_x_xxxz_x,     \
                             to_x_x_xxxz_y, \
                             to_x_x_xxxz_z, \
                             to_xxxxz_0,    \
                             to_xxxxz_xx,   \
                             to_xxxxz_xy,   \
                             to_xxxxz_xz,   \
                             to_xxz_0,      \
                             to_xxz_xx,     \
                             to_xxz_xy,     \
                             to_xxz_xz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxxz_x[k] = 3.0 * to_xxz_0[k] - 6.0 * to_xxz_xx[k] * tke_0 - 2.0 * to_xxxxz_0[k] * tbe_0 + 4.0 * to_xxxxz_xx[k] * tbe_0 * tke_0;

            to_x_x_xxxz_y[k] = -6.0 * to_xxz_xy[k] * tke_0 + 4.0 * to_xxxxz_xy[k] * tbe_0 * tke_0;

            to_x_x_xxxz_z[k] = -6.0 * to_xxz_xz[k] * tke_0 + 4.0 * to_xxxxz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 9-12 components of targeted buffer : GP

        auto to_x_x_xxyy_x = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 9);

        auto to_x_x_xxyy_y = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 10);

        auto to_x_x_xxyy_z = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 11);

#pragma omp simd aligned(to_x_x_xxyy_x,     \
                             to_x_x_xxyy_y, \
                             to_x_x_xxyy_z, \
                             to_xxxyy_0,    \
                             to_xxxyy_xx,   \
                             to_xxxyy_xy,   \
                             to_xxxyy_xz,   \
                             to_xyy_0,      \
                             to_xyy_xx,     \
                             to_xyy_xy,     \
                             to_xyy_xz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxyy_x[k] = 2.0 * to_xyy_0[k] - 4.0 * to_xyy_xx[k] * tke_0 - 2.0 * to_xxxyy_0[k] * tbe_0 + 4.0 * to_xxxyy_xx[k] * tbe_0 * tke_0;

            to_x_x_xxyy_y[k] = -4.0 * to_xyy_xy[k] * tke_0 + 4.0 * to_xxxyy_xy[k] * tbe_0 * tke_0;

            to_x_x_xxyy_z[k] = -4.0 * to_xyy_xz[k] * tke_0 + 4.0 * to_xxxyy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 12-15 components of targeted buffer : GP

        auto to_x_x_xxyz_x = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 12);

        auto to_x_x_xxyz_y = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 13);

        auto to_x_x_xxyz_z = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 14);

#pragma omp simd aligned(to_x_x_xxyz_x,     \
                             to_x_x_xxyz_y, \
                             to_x_x_xxyz_z, \
                             to_xxxyz_0,    \
                             to_xxxyz_xx,   \
                             to_xxxyz_xy,   \
                             to_xxxyz_xz,   \
                             to_xyz_0,      \
                             to_xyz_xx,     \
                             to_xyz_xy,     \
                             to_xyz_xz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxyz_x[k] = 2.0 * to_xyz_0[k] - 4.0 * to_xyz_xx[k] * tke_0 - 2.0 * to_xxxyz_0[k] * tbe_0 + 4.0 * to_xxxyz_xx[k] * tbe_0 * tke_0;

            to_x_x_xxyz_y[k] = -4.0 * to_xyz_xy[k] * tke_0 + 4.0 * to_xxxyz_xy[k] * tbe_0 * tke_0;

            to_x_x_xxyz_z[k] = -4.0 * to_xyz_xz[k] * tke_0 + 4.0 * to_xxxyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 15-18 components of targeted buffer : GP

        auto to_x_x_xxzz_x = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 15);

        auto to_x_x_xxzz_y = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 16);

        auto to_x_x_xxzz_z = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 17);

#pragma omp simd aligned(to_x_x_xxzz_x,     \
                             to_x_x_xxzz_y, \
                             to_x_x_xxzz_z, \
                             to_xxxzz_0,    \
                             to_xxxzz_xx,   \
                             to_xxxzz_xy,   \
                             to_xxxzz_xz,   \
                             to_xzz_0,      \
                             to_xzz_xx,     \
                             to_xzz_xy,     \
                             to_xzz_xz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxzz_x[k] = 2.0 * to_xzz_0[k] - 4.0 * to_xzz_xx[k] * tke_0 - 2.0 * to_xxxzz_0[k] * tbe_0 + 4.0 * to_xxxzz_xx[k] * tbe_0 * tke_0;

            to_x_x_xxzz_y[k] = -4.0 * to_xzz_xy[k] * tke_0 + 4.0 * to_xxxzz_xy[k] * tbe_0 * tke_0;

            to_x_x_xxzz_z[k] = -4.0 * to_xzz_xz[k] * tke_0 + 4.0 * to_xxxzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 18-21 components of targeted buffer : GP

        auto to_x_x_xyyy_x = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 18);

        auto to_x_x_xyyy_y = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 19);

        auto to_x_x_xyyy_z = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 20);

#pragma omp simd aligned(to_x_x_xyyy_x,     \
                             to_x_x_xyyy_y, \
                             to_x_x_xyyy_z, \
                             to_xxyyy_0,    \
                             to_xxyyy_xx,   \
                             to_xxyyy_xy,   \
                             to_xxyyy_xz,   \
                             to_yyy_0,      \
                             to_yyy_xx,     \
                             to_yyy_xy,     \
                             to_yyy_xz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyyy_x[k] = to_yyy_0[k] - 2.0 * to_yyy_xx[k] * tke_0 - 2.0 * to_xxyyy_0[k] * tbe_0 + 4.0 * to_xxyyy_xx[k] * tbe_0 * tke_0;

            to_x_x_xyyy_y[k] = -2.0 * to_yyy_xy[k] * tke_0 + 4.0 * to_xxyyy_xy[k] * tbe_0 * tke_0;

            to_x_x_xyyy_z[k] = -2.0 * to_yyy_xz[k] * tke_0 + 4.0 * to_xxyyy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 21-24 components of targeted buffer : GP

        auto to_x_x_xyyz_x = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 21);

        auto to_x_x_xyyz_y = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 22);

        auto to_x_x_xyyz_z = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 23);

#pragma omp simd aligned(to_x_x_xyyz_x,     \
                             to_x_x_xyyz_y, \
                             to_x_x_xyyz_z, \
                             to_xxyyz_0,    \
                             to_xxyyz_xx,   \
                             to_xxyyz_xy,   \
                             to_xxyyz_xz,   \
                             to_yyz_0,      \
                             to_yyz_xx,     \
                             to_yyz_xy,     \
                             to_yyz_xz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyyz_x[k] = to_yyz_0[k] - 2.0 * to_yyz_xx[k] * tke_0 - 2.0 * to_xxyyz_0[k] * tbe_0 + 4.0 * to_xxyyz_xx[k] * tbe_0 * tke_0;

            to_x_x_xyyz_y[k] = -2.0 * to_yyz_xy[k] * tke_0 + 4.0 * to_xxyyz_xy[k] * tbe_0 * tke_0;

            to_x_x_xyyz_z[k] = -2.0 * to_yyz_xz[k] * tke_0 + 4.0 * to_xxyyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 24-27 components of targeted buffer : GP

        auto to_x_x_xyzz_x = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 24);

        auto to_x_x_xyzz_y = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 25);

        auto to_x_x_xyzz_z = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 26);

#pragma omp simd aligned(to_x_x_xyzz_x,     \
                             to_x_x_xyzz_y, \
                             to_x_x_xyzz_z, \
                             to_xxyzz_0,    \
                             to_xxyzz_xx,   \
                             to_xxyzz_xy,   \
                             to_xxyzz_xz,   \
                             to_yzz_0,      \
                             to_yzz_xx,     \
                             to_yzz_xy,     \
                             to_yzz_xz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyzz_x[k] = to_yzz_0[k] - 2.0 * to_yzz_xx[k] * tke_0 - 2.0 * to_xxyzz_0[k] * tbe_0 + 4.0 * to_xxyzz_xx[k] * tbe_0 * tke_0;

            to_x_x_xyzz_y[k] = -2.0 * to_yzz_xy[k] * tke_0 + 4.0 * to_xxyzz_xy[k] * tbe_0 * tke_0;

            to_x_x_xyzz_z[k] = -2.0 * to_yzz_xz[k] * tke_0 + 4.0 * to_xxyzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 27-30 components of targeted buffer : GP

        auto to_x_x_xzzz_x = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 27);

        auto to_x_x_xzzz_y = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 28);

        auto to_x_x_xzzz_z = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 29);

#pragma omp simd aligned(to_x_x_xzzz_x,     \
                             to_x_x_xzzz_y, \
                             to_x_x_xzzz_z, \
                             to_xxzzz_0,    \
                             to_xxzzz_xx,   \
                             to_xxzzz_xy,   \
                             to_xxzzz_xz,   \
                             to_zzz_0,      \
                             to_zzz_xx,     \
                             to_zzz_xy,     \
                             to_zzz_xz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xzzz_x[k] = to_zzz_0[k] - 2.0 * to_zzz_xx[k] * tke_0 - 2.0 * to_xxzzz_0[k] * tbe_0 + 4.0 * to_xxzzz_xx[k] * tbe_0 * tke_0;

            to_x_x_xzzz_y[k] = -2.0 * to_zzz_xy[k] * tke_0 + 4.0 * to_xxzzz_xy[k] * tbe_0 * tke_0;

            to_x_x_xzzz_z[k] = -2.0 * to_zzz_xz[k] * tke_0 + 4.0 * to_xxzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 30-33 components of targeted buffer : GP

        auto to_x_x_yyyy_x = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 30);

        auto to_x_x_yyyy_y = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 31);

        auto to_x_x_yyyy_z = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 32);

#pragma omp simd aligned(to_x_x_yyyy_x, to_x_x_yyyy_y, to_x_x_yyyy_z, to_xyyyy_0, to_xyyyy_xx, to_xyyyy_xy, to_xyyyy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyyy_x[k] = -2.0 * to_xyyyy_0[k] * tbe_0 + 4.0 * to_xyyyy_xx[k] * tbe_0 * tke_0;

            to_x_x_yyyy_y[k] = 4.0 * to_xyyyy_xy[k] * tbe_0 * tke_0;

            to_x_x_yyyy_z[k] = 4.0 * to_xyyyy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 33-36 components of targeted buffer : GP

        auto to_x_x_yyyz_x = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 33);

        auto to_x_x_yyyz_y = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 34);

        auto to_x_x_yyyz_z = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 35);

#pragma omp simd aligned(to_x_x_yyyz_x, to_x_x_yyyz_y, to_x_x_yyyz_z, to_xyyyz_0, to_xyyyz_xx, to_xyyyz_xy, to_xyyyz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyyz_x[k] = -2.0 * to_xyyyz_0[k] * tbe_0 + 4.0 * to_xyyyz_xx[k] * tbe_0 * tke_0;

            to_x_x_yyyz_y[k] = 4.0 * to_xyyyz_xy[k] * tbe_0 * tke_0;

            to_x_x_yyyz_z[k] = 4.0 * to_xyyyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 36-39 components of targeted buffer : GP

        auto to_x_x_yyzz_x = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 36);

        auto to_x_x_yyzz_y = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 37);

        auto to_x_x_yyzz_z = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 38);

#pragma omp simd aligned(to_x_x_yyzz_x, to_x_x_yyzz_y, to_x_x_yyzz_z, to_xyyzz_0, to_xyyzz_xx, to_xyyzz_xy, to_xyyzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyzz_x[k] = -2.0 * to_xyyzz_0[k] * tbe_0 + 4.0 * to_xyyzz_xx[k] * tbe_0 * tke_0;

            to_x_x_yyzz_y[k] = 4.0 * to_xyyzz_xy[k] * tbe_0 * tke_0;

            to_x_x_yyzz_z[k] = 4.0 * to_xyyzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 39-42 components of targeted buffer : GP

        auto to_x_x_yzzz_x = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 39);

        auto to_x_x_yzzz_y = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 40);

        auto to_x_x_yzzz_z = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 41);

#pragma omp simd aligned(to_x_x_yzzz_x, to_x_x_yzzz_y, to_x_x_yzzz_z, to_xyzzz_0, to_xyzzz_xx, to_xyzzz_xy, to_xyzzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yzzz_x[k] = -2.0 * to_xyzzz_0[k] * tbe_0 + 4.0 * to_xyzzz_xx[k] * tbe_0 * tke_0;

            to_x_x_yzzz_y[k] = 4.0 * to_xyzzz_xy[k] * tbe_0 * tke_0;

            to_x_x_yzzz_z[k] = 4.0 * to_xyzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 42-45 components of targeted buffer : GP

        auto to_x_x_zzzz_x = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 42);

        auto to_x_x_zzzz_y = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 43);

        auto to_x_x_zzzz_z = pbuffer.data(idx_op_geom_101_gp + 0 * op_comps * 45 + i * 45 + 44);

#pragma omp simd aligned(to_x_x_zzzz_x, to_x_x_zzzz_y, to_x_x_zzzz_z, to_xzzzz_0, to_xzzzz_xx, to_xzzzz_xy, to_xzzzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_zzzz_x[k] = -2.0 * to_xzzzz_0[k] * tbe_0 + 4.0 * to_xzzzz_xx[k] * tbe_0 * tke_0;

            to_x_x_zzzz_y[k] = 4.0 * to_xzzzz_xy[k] * tbe_0 * tke_0;

            to_x_x_zzzz_z[k] = 4.0 * to_xzzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 45-48 components of targeted buffer : GP

        auto to_x_y_xxxx_x = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 0);

        auto to_x_y_xxxx_y = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 1);

        auto to_x_y_xxxx_z = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 2);

#pragma omp simd aligned(to_x_y_xxxx_x,     \
                             to_x_y_xxxx_y, \
                             to_x_y_xxxx_z, \
                             to_xxx_0,      \
                             to_xxx_xy,     \
                             to_xxx_yy,     \
                             to_xxx_yz,     \
                             to_xxxxx_0,    \
                             to_xxxxx_xy,   \
                             to_xxxxx_yy,   \
                             to_xxxxx_yz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxxx_x[k] = -8.0 * to_xxx_xy[k] * tke_0 + 4.0 * to_xxxxx_xy[k] * tbe_0 * tke_0;

            to_x_y_xxxx_y[k] = 4.0 * to_xxx_0[k] - 8.0 * to_xxx_yy[k] * tke_0 - 2.0 * to_xxxxx_0[k] * tbe_0 + 4.0 * to_xxxxx_yy[k] * tbe_0 * tke_0;

            to_x_y_xxxx_z[k] = -8.0 * to_xxx_yz[k] * tke_0 + 4.0 * to_xxxxx_yz[k] * tbe_0 * tke_0;
        }

        // Set up 48-51 components of targeted buffer : GP

        auto to_x_y_xxxy_x = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 3);

        auto to_x_y_xxxy_y = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 4);

        auto to_x_y_xxxy_z = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 5);

#pragma omp simd aligned(to_x_y_xxxy_x,     \
                             to_x_y_xxxy_y, \
                             to_x_y_xxxy_z, \
                             to_xxxxy_0,    \
                             to_xxxxy_xy,   \
                             to_xxxxy_yy,   \
                             to_xxxxy_yz,   \
                             to_xxy_0,      \
                             to_xxy_xy,     \
                             to_xxy_yy,     \
                             to_xxy_yz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxxy_x[k] = -6.0 * to_xxy_xy[k] * tke_0 + 4.0 * to_xxxxy_xy[k] * tbe_0 * tke_0;

            to_x_y_xxxy_y[k] = 3.0 * to_xxy_0[k] - 6.0 * to_xxy_yy[k] * tke_0 - 2.0 * to_xxxxy_0[k] * tbe_0 + 4.0 * to_xxxxy_yy[k] * tbe_0 * tke_0;

            to_x_y_xxxy_z[k] = -6.0 * to_xxy_yz[k] * tke_0 + 4.0 * to_xxxxy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 51-54 components of targeted buffer : GP

        auto to_x_y_xxxz_x = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 6);

        auto to_x_y_xxxz_y = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 7);

        auto to_x_y_xxxz_z = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 8);

#pragma omp simd aligned(to_x_y_xxxz_x,     \
                             to_x_y_xxxz_y, \
                             to_x_y_xxxz_z, \
                             to_xxxxz_0,    \
                             to_xxxxz_xy,   \
                             to_xxxxz_yy,   \
                             to_xxxxz_yz,   \
                             to_xxz_0,      \
                             to_xxz_xy,     \
                             to_xxz_yy,     \
                             to_xxz_yz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxxz_x[k] = -6.0 * to_xxz_xy[k] * tke_0 + 4.0 * to_xxxxz_xy[k] * tbe_0 * tke_0;

            to_x_y_xxxz_y[k] = 3.0 * to_xxz_0[k] - 6.0 * to_xxz_yy[k] * tke_0 - 2.0 * to_xxxxz_0[k] * tbe_0 + 4.0 * to_xxxxz_yy[k] * tbe_0 * tke_0;

            to_x_y_xxxz_z[k] = -6.0 * to_xxz_yz[k] * tke_0 + 4.0 * to_xxxxz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 54-57 components of targeted buffer : GP

        auto to_x_y_xxyy_x = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 9);

        auto to_x_y_xxyy_y = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 10);

        auto to_x_y_xxyy_z = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 11);

#pragma omp simd aligned(to_x_y_xxyy_x,     \
                             to_x_y_xxyy_y, \
                             to_x_y_xxyy_z, \
                             to_xxxyy_0,    \
                             to_xxxyy_xy,   \
                             to_xxxyy_yy,   \
                             to_xxxyy_yz,   \
                             to_xyy_0,      \
                             to_xyy_xy,     \
                             to_xyy_yy,     \
                             to_xyy_yz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxyy_x[k] = -4.0 * to_xyy_xy[k] * tke_0 + 4.0 * to_xxxyy_xy[k] * tbe_0 * tke_0;

            to_x_y_xxyy_y[k] = 2.0 * to_xyy_0[k] - 4.0 * to_xyy_yy[k] * tke_0 - 2.0 * to_xxxyy_0[k] * tbe_0 + 4.0 * to_xxxyy_yy[k] * tbe_0 * tke_0;

            to_x_y_xxyy_z[k] = -4.0 * to_xyy_yz[k] * tke_0 + 4.0 * to_xxxyy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 57-60 components of targeted buffer : GP

        auto to_x_y_xxyz_x = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 12);

        auto to_x_y_xxyz_y = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 13);

        auto to_x_y_xxyz_z = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 14);

#pragma omp simd aligned(to_x_y_xxyz_x,     \
                             to_x_y_xxyz_y, \
                             to_x_y_xxyz_z, \
                             to_xxxyz_0,    \
                             to_xxxyz_xy,   \
                             to_xxxyz_yy,   \
                             to_xxxyz_yz,   \
                             to_xyz_0,      \
                             to_xyz_xy,     \
                             to_xyz_yy,     \
                             to_xyz_yz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxyz_x[k] = -4.0 * to_xyz_xy[k] * tke_0 + 4.0 * to_xxxyz_xy[k] * tbe_0 * tke_0;

            to_x_y_xxyz_y[k] = 2.0 * to_xyz_0[k] - 4.0 * to_xyz_yy[k] * tke_0 - 2.0 * to_xxxyz_0[k] * tbe_0 + 4.0 * to_xxxyz_yy[k] * tbe_0 * tke_0;

            to_x_y_xxyz_z[k] = -4.0 * to_xyz_yz[k] * tke_0 + 4.0 * to_xxxyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 60-63 components of targeted buffer : GP

        auto to_x_y_xxzz_x = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 15);

        auto to_x_y_xxzz_y = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 16);

        auto to_x_y_xxzz_z = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 17);

#pragma omp simd aligned(to_x_y_xxzz_x,     \
                             to_x_y_xxzz_y, \
                             to_x_y_xxzz_z, \
                             to_xxxzz_0,    \
                             to_xxxzz_xy,   \
                             to_xxxzz_yy,   \
                             to_xxxzz_yz,   \
                             to_xzz_0,      \
                             to_xzz_xy,     \
                             to_xzz_yy,     \
                             to_xzz_yz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxzz_x[k] = -4.0 * to_xzz_xy[k] * tke_0 + 4.0 * to_xxxzz_xy[k] * tbe_0 * tke_0;

            to_x_y_xxzz_y[k] = 2.0 * to_xzz_0[k] - 4.0 * to_xzz_yy[k] * tke_0 - 2.0 * to_xxxzz_0[k] * tbe_0 + 4.0 * to_xxxzz_yy[k] * tbe_0 * tke_0;

            to_x_y_xxzz_z[k] = -4.0 * to_xzz_yz[k] * tke_0 + 4.0 * to_xxxzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 63-66 components of targeted buffer : GP

        auto to_x_y_xyyy_x = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 18);

        auto to_x_y_xyyy_y = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 19);

        auto to_x_y_xyyy_z = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 20);

#pragma omp simd aligned(to_x_y_xyyy_x,     \
                             to_x_y_xyyy_y, \
                             to_x_y_xyyy_z, \
                             to_xxyyy_0,    \
                             to_xxyyy_xy,   \
                             to_xxyyy_yy,   \
                             to_xxyyy_yz,   \
                             to_yyy_0,      \
                             to_yyy_xy,     \
                             to_yyy_yy,     \
                             to_yyy_yz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyyy_x[k] = -2.0 * to_yyy_xy[k] * tke_0 + 4.0 * to_xxyyy_xy[k] * tbe_0 * tke_0;

            to_x_y_xyyy_y[k] = to_yyy_0[k] - 2.0 * to_yyy_yy[k] * tke_0 - 2.0 * to_xxyyy_0[k] * tbe_0 + 4.0 * to_xxyyy_yy[k] * tbe_0 * tke_0;

            to_x_y_xyyy_z[k] = -2.0 * to_yyy_yz[k] * tke_0 + 4.0 * to_xxyyy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 66-69 components of targeted buffer : GP

        auto to_x_y_xyyz_x = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 21);

        auto to_x_y_xyyz_y = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 22);

        auto to_x_y_xyyz_z = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 23);

#pragma omp simd aligned(to_x_y_xyyz_x,     \
                             to_x_y_xyyz_y, \
                             to_x_y_xyyz_z, \
                             to_xxyyz_0,    \
                             to_xxyyz_xy,   \
                             to_xxyyz_yy,   \
                             to_xxyyz_yz,   \
                             to_yyz_0,      \
                             to_yyz_xy,     \
                             to_yyz_yy,     \
                             to_yyz_yz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyyz_x[k] = -2.0 * to_yyz_xy[k] * tke_0 + 4.0 * to_xxyyz_xy[k] * tbe_0 * tke_0;

            to_x_y_xyyz_y[k] = to_yyz_0[k] - 2.0 * to_yyz_yy[k] * tke_0 - 2.0 * to_xxyyz_0[k] * tbe_0 + 4.0 * to_xxyyz_yy[k] * tbe_0 * tke_0;

            to_x_y_xyyz_z[k] = -2.0 * to_yyz_yz[k] * tke_0 + 4.0 * to_xxyyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 69-72 components of targeted buffer : GP

        auto to_x_y_xyzz_x = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 24);

        auto to_x_y_xyzz_y = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 25);

        auto to_x_y_xyzz_z = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 26);

#pragma omp simd aligned(to_x_y_xyzz_x,     \
                             to_x_y_xyzz_y, \
                             to_x_y_xyzz_z, \
                             to_xxyzz_0,    \
                             to_xxyzz_xy,   \
                             to_xxyzz_yy,   \
                             to_xxyzz_yz,   \
                             to_yzz_0,      \
                             to_yzz_xy,     \
                             to_yzz_yy,     \
                             to_yzz_yz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyzz_x[k] = -2.0 * to_yzz_xy[k] * tke_0 + 4.0 * to_xxyzz_xy[k] * tbe_0 * tke_0;

            to_x_y_xyzz_y[k] = to_yzz_0[k] - 2.0 * to_yzz_yy[k] * tke_0 - 2.0 * to_xxyzz_0[k] * tbe_0 + 4.0 * to_xxyzz_yy[k] * tbe_0 * tke_0;

            to_x_y_xyzz_z[k] = -2.0 * to_yzz_yz[k] * tke_0 + 4.0 * to_xxyzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 72-75 components of targeted buffer : GP

        auto to_x_y_xzzz_x = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 27);

        auto to_x_y_xzzz_y = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 28);

        auto to_x_y_xzzz_z = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 29);

#pragma omp simd aligned(to_x_y_xzzz_x,     \
                             to_x_y_xzzz_y, \
                             to_x_y_xzzz_z, \
                             to_xxzzz_0,    \
                             to_xxzzz_xy,   \
                             to_xxzzz_yy,   \
                             to_xxzzz_yz,   \
                             to_zzz_0,      \
                             to_zzz_xy,     \
                             to_zzz_yy,     \
                             to_zzz_yz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xzzz_x[k] = -2.0 * to_zzz_xy[k] * tke_0 + 4.0 * to_xxzzz_xy[k] * tbe_0 * tke_0;

            to_x_y_xzzz_y[k] = to_zzz_0[k] - 2.0 * to_zzz_yy[k] * tke_0 - 2.0 * to_xxzzz_0[k] * tbe_0 + 4.0 * to_xxzzz_yy[k] * tbe_0 * tke_0;

            to_x_y_xzzz_z[k] = -2.0 * to_zzz_yz[k] * tke_0 + 4.0 * to_xxzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 75-78 components of targeted buffer : GP

        auto to_x_y_yyyy_x = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 30);

        auto to_x_y_yyyy_y = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 31);

        auto to_x_y_yyyy_z = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 32);

#pragma omp simd aligned(to_x_y_yyyy_x, to_x_y_yyyy_y, to_x_y_yyyy_z, to_xyyyy_0, to_xyyyy_xy, to_xyyyy_yy, to_xyyyy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyyy_x[k] = 4.0 * to_xyyyy_xy[k] * tbe_0 * tke_0;

            to_x_y_yyyy_y[k] = -2.0 * to_xyyyy_0[k] * tbe_0 + 4.0 * to_xyyyy_yy[k] * tbe_0 * tke_0;

            to_x_y_yyyy_z[k] = 4.0 * to_xyyyy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 78-81 components of targeted buffer : GP

        auto to_x_y_yyyz_x = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 33);

        auto to_x_y_yyyz_y = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 34);

        auto to_x_y_yyyz_z = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 35);

#pragma omp simd aligned(to_x_y_yyyz_x, to_x_y_yyyz_y, to_x_y_yyyz_z, to_xyyyz_0, to_xyyyz_xy, to_xyyyz_yy, to_xyyyz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyyz_x[k] = 4.0 * to_xyyyz_xy[k] * tbe_0 * tke_0;

            to_x_y_yyyz_y[k] = -2.0 * to_xyyyz_0[k] * tbe_0 + 4.0 * to_xyyyz_yy[k] * tbe_0 * tke_0;

            to_x_y_yyyz_z[k] = 4.0 * to_xyyyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 81-84 components of targeted buffer : GP

        auto to_x_y_yyzz_x = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 36);

        auto to_x_y_yyzz_y = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 37);

        auto to_x_y_yyzz_z = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 38);

#pragma omp simd aligned(to_x_y_yyzz_x, to_x_y_yyzz_y, to_x_y_yyzz_z, to_xyyzz_0, to_xyyzz_xy, to_xyyzz_yy, to_xyyzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyzz_x[k] = 4.0 * to_xyyzz_xy[k] * tbe_0 * tke_0;

            to_x_y_yyzz_y[k] = -2.0 * to_xyyzz_0[k] * tbe_0 + 4.0 * to_xyyzz_yy[k] * tbe_0 * tke_0;

            to_x_y_yyzz_z[k] = 4.0 * to_xyyzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 84-87 components of targeted buffer : GP

        auto to_x_y_yzzz_x = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 39);

        auto to_x_y_yzzz_y = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 40);

        auto to_x_y_yzzz_z = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 41);

#pragma omp simd aligned(to_x_y_yzzz_x, to_x_y_yzzz_y, to_x_y_yzzz_z, to_xyzzz_0, to_xyzzz_xy, to_xyzzz_yy, to_xyzzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yzzz_x[k] = 4.0 * to_xyzzz_xy[k] * tbe_0 * tke_0;

            to_x_y_yzzz_y[k] = -2.0 * to_xyzzz_0[k] * tbe_0 + 4.0 * to_xyzzz_yy[k] * tbe_0 * tke_0;

            to_x_y_yzzz_z[k] = 4.0 * to_xyzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 87-90 components of targeted buffer : GP

        auto to_x_y_zzzz_x = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 42);

        auto to_x_y_zzzz_y = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 43);

        auto to_x_y_zzzz_z = pbuffer.data(idx_op_geom_101_gp + 1 * op_comps * 45 + i * 45 + 44);

#pragma omp simd aligned(to_x_y_zzzz_x, to_x_y_zzzz_y, to_x_y_zzzz_z, to_xzzzz_0, to_xzzzz_xy, to_xzzzz_yy, to_xzzzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_zzzz_x[k] = 4.0 * to_xzzzz_xy[k] * tbe_0 * tke_0;

            to_x_y_zzzz_y[k] = -2.0 * to_xzzzz_0[k] * tbe_0 + 4.0 * to_xzzzz_yy[k] * tbe_0 * tke_0;

            to_x_y_zzzz_z[k] = 4.0 * to_xzzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 90-93 components of targeted buffer : GP

        auto to_x_z_xxxx_x = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 0);

        auto to_x_z_xxxx_y = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 1);

        auto to_x_z_xxxx_z = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 2);

#pragma omp simd aligned(to_x_z_xxxx_x,     \
                             to_x_z_xxxx_y, \
                             to_x_z_xxxx_z, \
                             to_xxx_0,      \
                             to_xxx_xz,     \
                             to_xxx_yz,     \
                             to_xxx_zz,     \
                             to_xxxxx_0,    \
                             to_xxxxx_xz,   \
                             to_xxxxx_yz,   \
                             to_xxxxx_zz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxxx_x[k] = -8.0 * to_xxx_xz[k] * tke_0 + 4.0 * to_xxxxx_xz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_y[k] = -8.0 * to_xxx_yz[k] * tke_0 + 4.0 * to_xxxxx_yz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_z[k] = 4.0 * to_xxx_0[k] - 8.0 * to_xxx_zz[k] * tke_0 - 2.0 * to_xxxxx_0[k] * tbe_0 + 4.0 * to_xxxxx_zz[k] * tbe_0 * tke_0;
        }

        // Set up 93-96 components of targeted buffer : GP

        auto to_x_z_xxxy_x = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 3);

        auto to_x_z_xxxy_y = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 4);

        auto to_x_z_xxxy_z = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 5);

#pragma omp simd aligned(to_x_z_xxxy_x,     \
                             to_x_z_xxxy_y, \
                             to_x_z_xxxy_z, \
                             to_xxxxy_0,    \
                             to_xxxxy_xz,   \
                             to_xxxxy_yz,   \
                             to_xxxxy_zz,   \
                             to_xxy_0,      \
                             to_xxy_xz,     \
                             to_xxy_yz,     \
                             to_xxy_zz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxxy_x[k] = -6.0 * to_xxy_xz[k] * tke_0 + 4.0 * to_xxxxy_xz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_y[k] = -6.0 * to_xxy_yz[k] * tke_0 + 4.0 * to_xxxxy_yz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_z[k] = 3.0 * to_xxy_0[k] - 6.0 * to_xxy_zz[k] * tke_0 - 2.0 * to_xxxxy_0[k] * tbe_0 + 4.0 * to_xxxxy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 96-99 components of targeted buffer : GP

        auto to_x_z_xxxz_x = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 6);

        auto to_x_z_xxxz_y = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 7);

        auto to_x_z_xxxz_z = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 8);

#pragma omp simd aligned(to_x_z_xxxz_x,     \
                             to_x_z_xxxz_y, \
                             to_x_z_xxxz_z, \
                             to_xxxxz_0,    \
                             to_xxxxz_xz,   \
                             to_xxxxz_yz,   \
                             to_xxxxz_zz,   \
                             to_xxz_0,      \
                             to_xxz_xz,     \
                             to_xxz_yz,     \
                             to_xxz_zz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxxz_x[k] = -6.0 * to_xxz_xz[k] * tke_0 + 4.0 * to_xxxxz_xz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_y[k] = -6.0 * to_xxz_yz[k] * tke_0 + 4.0 * to_xxxxz_yz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_z[k] = 3.0 * to_xxz_0[k] - 6.0 * to_xxz_zz[k] * tke_0 - 2.0 * to_xxxxz_0[k] * tbe_0 + 4.0 * to_xxxxz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 99-102 components of targeted buffer : GP

        auto to_x_z_xxyy_x = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 9);

        auto to_x_z_xxyy_y = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 10);

        auto to_x_z_xxyy_z = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 11);

#pragma omp simd aligned(to_x_z_xxyy_x,     \
                             to_x_z_xxyy_y, \
                             to_x_z_xxyy_z, \
                             to_xxxyy_0,    \
                             to_xxxyy_xz,   \
                             to_xxxyy_yz,   \
                             to_xxxyy_zz,   \
                             to_xyy_0,      \
                             to_xyy_xz,     \
                             to_xyy_yz,     \
                             to_xyy_zz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxyy_x[k] = -4.0 * to_xyy_xz[k] * tke_0 + 4.0 * to_xxxyy_xz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_y[k] = -4.0 * to_xyy_yz[k] * tke_0 + 4.0 * to_xxxyy_yz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_z[k] = 2.0 * to_xyy_0[k] - 4.0 * to_xyy_zz[k] * tke_0 - 2.0 * to_xxxyy_0[k] * tbe_0 + 4.0 * to_xxxyy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 102-105 components of targeted buffer : GP

        auto to_x_z_xxyz_x = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 12);

        auto to_x_z_xxyz_y = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 13);

        auto to_x_z_xxyz_z = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 14);

#pragma omp simd aligned(to_x_z_xxyz_x,     \
                             to_x_z_xxyz_y, \
                             to_x_z_xxyz_z, \
                             to_xxxyz_0,    \
                             to_xxxyz_xz,   \
                             to_xxxyz_yz,   \
                             to_xxxyz_zz,   \
                             to_xyz_0,      \
                             to_xyz_xz,     \
                             to_xyz_yz,     \
                             to_xyz_zz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxyz_x[k] = -4.0 * to_xyz_xz[k] * tke_0 + 4.0 * to_xxxyz_xz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_y[k] = -4.0 * to_xyz_yz[k] * tke_0 + 4.0 * to_xxxyz_yz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_z[k] = 2.0 * to_xyz_0[k] - 4.0 * to_xyz_zz[k] * tke_0 - 2.0 * to_xxxyz_0[k] * tbe_0 + 4.0 * to_xxxyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 105-108 components of targeted buffer : GP

        auto to_x_z_xxzz_x = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 15);

        auto to_x_z_xxzz_y = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 16);

        auto to_x_z_xxzz_z = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 17);

#pragma omp simd aligned(to_x_z_xxzz_x,     \
                             to_x_z_xxzz_y, \
                             to_x_z_xxzz_z, \
                             to_xxxzz_0,    \
                             to_xxxzz_xz,   \
                             to_xxxzz_yz,   \
                             to_xxxzz_zz,   \
                             to_xzz_0,      \
                             to_xzz_xz,     \
                             to_xzz_yz,     \
                             to_xzz_zz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxzz_x[k] = -4.0 * to_xzz_xz[k] * tke_0 + 4.0 * to_xxxzz_xz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_y[k] = -4.0 * to_xzz_yz[k] * tke_0 + 4.0 * to_xxxzz_yz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_z[k] = 2.0 * to_xzz_0[k] - 4.0 * to_xzz_zz[k] * tke_0 - 2.0 * to_xxxzz_0[k] * tbe_0 + 4.0 * to_xxxzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 108-111 components of targeted buffer : GP

        auto to_x_z_xyyy_x = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 18);

        auto to_x_z_xyyy_y = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 19);

        auto to_x_z_xyyy_z = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 20);

#pragma omp simd aligned(to_x_z_xyyy_x,     \
                             to_x_z_xyyy_y, \
                             to_x_z_xyyy_z, \
                             to_xxyyy_0,    \
                             to_xxyyy_xz,   \
                             to_xxyyy_yz,   \
                             to_xxyyy_zz,   \
                             to_yyy_0,      \
                             to_yyy_xz,     \
                             to_yyy_yz,     \
                             to_yyy_zz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyyy_x[k] = -2.0 * to_yyy_xz[k] * tke_0 + 4.0 * to_xxyyy_xz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_y[k] = -2.0 * to_yyy_yz[k] * tke_0 + 4.0 * to_xxyyy_yz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_z[k] = to_yyy_0[k] - 2.0 * to_yyy_zz[k] * tke_0 - 2.0 * to_xxyyy_0[k] * tbe_0 + 4.0 * to_xxyyy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 111-114 components of targeted buffer : GP

        auto to_x_z_xyyz_x = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 21);

        auto to_x_z_xyyz_y = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 22);

        auto to_x_z_xyyz_z = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 23);

#pragma omp simd aligned(to_x_z_xyyz_x,     \
                             to_x_z_xyyz_y, \
                             to_x_z_xyyz_z, \
                             to_xxyyz_0,    \
                             to_xxyyz_xz,   \
                             to_xxyyz_yz,   \
                             to_xxyyz_zz,   \
                             to_yyz_0,      \
                             to_yyz_xz,     \
                             to_yyz_yz,     \
                             to_yyz_zz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyyz_x[k] = -2.0 * to_yyz_xz[k] * tke_0 + 4.0 * to_xxyyz_xz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_y[k] = -2.0 * to_yyz_yz[k] * tke_0 + 4.0 * to_xxyyz_yz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_z[k] = to_yyz_0[k] - 2.0 * to_yyz_zz[k] * tke_0 - 2.0 * to_xxyyz_0[k] * tbe_0 + 4.0 * to_xxyyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 114-117 components of targeted buffer : GP

        auto to_x_z_xyzz_x = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 24);

        auto to_x_z_xyzz_y = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 25);

        auto to_x_z_xyzz_z = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 26);

#pragma omp simd aligned(to_x_z_xyzz_x,     \
                             to_x_z_xyzz_y, \
                             to_x_z_xyzz_z, \
                             to_xxyzz_0,    \
                             to_xxyzz_xz,   \
                             to_xxyzz_yz,   \
                             to_xxyzz_zz,   \
                             to_yzz_0,      \
                             to_yzz_xz,     \
                             to_yzz_yz,     \
                             to_yzz_zz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyzz_x[k] = -2.0 * to_yzz_xz[k] * tke_0 + 4.0 * to_xxyzz_xz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_y[k] = -2.0 * to_yzz_yz[k] * tke_0 + 4.0 * to_xxyzz_yz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_z[k] = to_yzz_0[k] - 2.0 * to_yzz_zz[k] * tke_0 - 2.0 * to_xxyzz_0[k] * tbe_0 + 4.0 * to_xxyzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 117-120 components of targeted buffer : GP

        auto to_x_z_xzzz_x = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 27);

        auto to_x_z_xzzz_y = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 28);

        auto to_x_z_xzzz_z = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 29);

#pragma omp simd aligned(to_x_z_xzzz_x,     \
                             to_x_z_xzzz_y, \
                             to_x_z_xzzz_z, \
                             to_xxzzz_0,    \
                             to_xxzzz_xz,   \
                             to_xxzzz_yz,   \
                             to_xxzzz_zz,   \
                             to_zzz_0,      \
                             to_zzz_xz,     \
                             to_zzz_yz,     \
                             to_zzz_zz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xzzz_x[k] = -2.0 * to_zzz_xz[k] * tke_0 + 4.0 * to_xxzzz_xz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_y[k] = -2.0 * to_zzz_yz[k] * tke_0 + 4.0 * to_xxzzz_yz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_z[k] = to_zzz_0[k] - 2.0 * to_zzz_zz[k] * tke_0 - 2.0 * to_xxzzz_0[k] * tbe_0 + 4.0 * to_xxzzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 120-123 components of targeted buffer : GP

        auto to_x_z_yyyy_x = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 30);

        auto to_x_z_yyyy_y = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 31);

        auto to_x_z_yyyy_z = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 32);

#pragma omp simd aligned(to_x_z_yyyy_x, to_x_z_yyyy_y, to_x_z_yyyy_z, to_xyyyy_0, to_xyyyy_xz, to_xyyyy_yz, to_xyyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyyy_x[k] = 4.0 * to_xyyyy_xz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_y[k] = 4.0 * to_xyyyy_yz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_z[k] = -2.0 * to_xyyyy_0[k] * tbe_0 + 4.0 * to_xyyyy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 123-126 components of targeted buffer : GP

        auto to_x_z_yyyz_x = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 33);

        auto to_x_z_yyyz_y = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 34);

        auto to_x_z_yyyz_z = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 35);

#pragma omp simd aligned(to_x_z_yyyz_x, to_x_z_yyyz_y, to_x_z_yyyz_z, to_xyyyz_0, to_xyyyz_xz, to_xyyyz_yz, to_xyyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyyz_x[k] = 4.0 * to_xyyyz_xz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_y[k] = 4.0 * to_xyyyz_yz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_z[k] = -2.0 * to_xyyyz_0[k] * tbe_0 + 4.0 * to_xyyyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 126-129 components of targeted buffer : GP

        auto to_x_z_yyzz_x = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 36);

        auto to_x_z_yyzz_y = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 37);

        auto to_x_z_yyzz_z = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 38);

#pragma omp simd aligned(to_x_z_yyzz_x, to_x_z_yyzz_y, to_x_z_yyzz_z, to_xyyzz_0, to_xyyzz_xz, to_xyyzz_yz, to_xyyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyzz_x[k] = 4.0 * to_xyyzz_xz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_y[k] = 4.0 * to_xyyzz_yz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_z[k] = -2.0 * to_xyyzz_0[k] * tbe_0 + 4.0 * to_xyyzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 129-132 components of targeted buffer : GP

        auto to_x_z_yzzz_x = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 39);

        auto to_x_z_yzzz_y = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 40);

        auto to_x_z_yzzz_z = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 41);

#pragma omp simd aligned(to_x_z_yzzz_x, to_x_z_yzzz_y, to_x_z_yzzz_z, to_xyzzz_0, to_xyzzz_xz, to_xyzzz_yz, to_xyzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yzzz_x[k] = 4.0 * to_xyzzz_xz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_y[k] = 4.0 * to_xyzzz_yz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_z[k] = -2.0 * to_xyzzz_0[k] * tbe_0 + 4.0 * to_xyzzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 132-135 components of targeted buffer : GP

        auto to_x_z_zzzz_x = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 42);

        auto to_x_z_zzzz_y = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 43);

        auto to_x_z_zzzz_z = pbuffer.data(idx_op_geom_101_gp + 2 * op_comps * 45 + i * 45 + 44);

#pragma omp simd aligned(to_x_z_zzzz_x, to_x_z_zzzz_y, to_x_z_zzzz_z, to_xzzzz_0, to_xzzzz_xz, to_xzzzz_yz, to_xzzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_zzzz_x[k] = 4.0 * to_xzzzz_xz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_y[k] = 4.0 * to_xzzzz_yz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_z[k] = -2.0 * to_xzzzz_0[k] * tbe_0 + 4.0 * to_xzzzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 135-138 components of targeted buffer : GP

        auto to_y_x_xxxx_x = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 0);

        auto to_y_x_xxxx_y = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 1);

        auto to_y_x_xxxx_z = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 2);

#pragma omp simd aligned(to_xxxxy_0, to_xxxxy_xx, to_xxxxy_xy, to_xxxxy_xz, to_y_x_xxxx_x, to_y_x_xxxx_y, to_y_x_xxxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxxx_x[k] = -2.0 * to_xxxxy_0[k] * tbe_0 + 4.0 * to_xxxxy_xx[k] * tbe_0 * tke_0;

            to_y_x_xxxx_y[k] = 4.0 * to_xxxxy_xy[k] * tbe_0 * tke_0;

            to_y_x_xxxx_z[k] = 4.0 * to_xxxxy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 138-141 components of targeted buffer : GP

        auto to_y_x_xxxy_x = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 3);

        auto to_y_x_xxxy_y = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 4);

        auto to_y_x_xxxy_z = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 5);

#pragma omp simd aligned(to_xxx_0,          \
                             to_xxx_xx,     \
                             to_xxx_xy,     \
                             to_xxx_xz,     \
                             to_xxxyy_0,    \
                             to_xxxyy_xx,   \
                             to_xxxyy_xy,   \
                             to_xxxyy_xz,   \
                             to_y_x_xxxy_x, \
                             to_y_x_xxxy_y, \
                             to_y_x_xxxy_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxxy_x[k] = to_xxx_0[k] - 2.0 * to_xxx_xx[k] * tke_0 - 2.0 * to_xxxyy_0[k] * tbe_0 + 4.0 * to_xxxyy_xx[k] * tbe_0 * tke_0;

            to_y_x_xxxy_y[k] = -2.0 * to_xxx_xy[k] * tke_0 + 4.0 * to_xxxyy_xy[k] * tbe_0 * tke_0;

            to_y_x_xxxy_z[k] = -2.0 * to_xxx_xz[k] * tke_0 + 4.0 * to_xxxyy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 141-144 components of targeted buffer : GP

        auto to_y_x_xxxz_x = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 6);

        auto to_y_x_xxxz_y = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 7);

        auto to_y_x_xxxz_z = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 8);

#pragma omp simd aligned(to_xxxyz_0, to_xxxyz_xx, to_xxxyz_xy, to_xxxyz_xz, to_y_x_xxxz_x, to_y_x_xxxz_y, to_y_x_xxxz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxxz_x[k] = -2.0 * to_xxxyz_0[k] * tbe_0 + 4.0 * to_xxxyz_xx[k] * tbe_0 * tke_0;

            to_y_x_xxxz_y[k] = 4.0 * to_xxxyz_xy[k] * tbe_0 * tke_0;

            to_y_x_xxxz_z[k] = 4.0 * to_xxxyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 144-147 components of targeted buffer : GP

        auto to_y_x_xxyy_x = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 9);

        auto to_y_x_xxyy_y = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 10);

        auto to_y_x_xxyy_z = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 11);

#pragma omp simd aligned(to_xxy_0,          \
                             to_xxy_xx,     \
                             to_xxy_xy,     \
                             to_xxy_xz,     \
                             to_xxyyy_0,    \
                             to_xxyyy_xx,   \
                             to_xxyyy_xy,   \
                             to_xxyyy_xz,   \
                             to_y_x_xxyy_x, \
                             to_y_x_xxyy_y, \
                             to_y_x_xxyy_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxyy_x[k] = 2.0 * to_xxy_0[k] - 4.0 * to_xxy_xx[k] * tke_0 - 2.0 * to_xxyyy_0[k] * tbe_0 + 4.0 * to_xxyyy_xx[k] * tbe_0 * tke_0;

            to_y_x_xxyy_y[k] = -4.0 * to_xxy_xy[k] * tke_0 + 4.0 * to_xxyyy_xy[k] * tbe_0 * tke_0;

            to_y_x_xxyy_z[k] = -4.0 * to_xxy_xz[k] * tke_0 + 4.0 * to_xxyyy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 147-150 components of targeted buffer : GP

        auto to_y_x_xxyz_x = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 12);

        auto to_y_x_xxyz_y = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 13);

        auto to_y_x_xxyz_z = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 14);

#pragma omp simd aligned(to_xxyyz_0,        \
                             to_xxyyz_xx,   \
                             to_xxyyz_xy,   \
                             to_xxyyz_xz,   \
                             to_xxz_0,      \
                             to_xxz_xx,     \
                             to_xxz_xy,     \
                             to_xxz_xz,     \
                             to_y_x_xxyz_x, \
                             to_y_x_xxyz_y, \
                             to_y_x_xxyz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxyz_x[k] = to_xxz_0[k] - 2.0 * to_xxz_xx[k] * tke_0 - 2.0 * to_xxyyz_0[k] * tbe_0 + 4.0 * to_xxyyz_xx[k] * tbe_0 * tke_0;

            to_y_x_xxyz_y[k] = -2.0 * to_xxz_xy[k] * tke_0 + 4.0 * to_xxyyz_xy[k] * tbe_0 * tke_0;

            to_y_x_xxyz_z[k] = -2.0 * to_xxz_xz[k] * tke_0 + 4.0 * to_xxyyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 150-153 components of targeted buffer : GP

        auto to_y_x_xxzz_x = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 15);

        auto to_y_x_xxzz_y = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 16);

        auto to_y_x_xxzz_z = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 17);

#pragma omp simd aligned(to_xxyzz_0, to_xxyzz_xx, to_xxyzz_xy, to_xxyzz_xz, to_y_x_xxzz_x, to_y_x_xxzz_y, to_y_x_xxzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxzz_x[k] = -2.0 * to_xxyzz_0[k] * tbe_0 + 4.0 * to_xxyzz_xx[k] * tbe_0 * tke_0;

            to_y_x_xxzz_y[k] = 4.0 * to_xxyzz_xy[k] * tbe_0 * tke_0;

            to_y_x_xxzz_z[k] = 4.0 * to_xxyzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 153-156 components of targeted buffer : GP

        auto to_y_x_xyyy_x = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 18);

        auto to_y_x_xyyy_y = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 19);

        auto to_y_x_xyyy_z = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 20);

#pragma omp simd aligned(to_xyy_0,          \
                             to_xyy_xx,     \
                             to_xyy_xy,     \
                             to_xyy_xz,     \
                             to_xyyyy_0,    \
                             to_xyyyy_xx,   \
                             to_xyyyy_xy,   \
                             to_xyyyy_xz,   \
                             to_y_x_xyyy_x, \
                             to_y_x_xyyy_y, \
                             to_y_x_xyyy_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyyy_x[k] = 3.0 * to_xyy_0[k] - 6.0 * to_xyy_xx[k] * tke_0 - 2.0 * to_xyyyy_0[k] * tbe_0 + 4.0 * to_xyyyy_xx[k] * tbe_0 * tke_0;

            to_y_x_xyyy_y[k] = -6.0 * to_xyy_xy[k] * tke_0 + 4.0 * to_xyyyy_xy[k] * tbe_0 * tke_0;

            to_y_x_xyyy_z[k] = -6.0 * to_xyy_xz[k] * tke_0 + 4.0 * to_xyyyy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 156-159 components of targeted buffer : GP

        auto to_y_x_xyyz_x = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 21);

        auto to_y_x_xyyz_y = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 22);

        auto to_y_x_xyyz_z = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 23);

#pragma omp simd aligned(to_xyyyz_0,        \
                             to_xyyyz_xx,   \
                             to_xyyyz_xy,   \
                             to_xyyyz_xz,   \
                             to_xyz_0,      \
                             to_xyz_xx,     \
                             to_xyz_xy,     \
                             to_xyz_xz,     \
                             to_y_x_xyyz_x, \
                             to_y_x_xyyz_y, \
                             to_y_x_xyyz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyyz_x[k] = 2.0 * to_xyz_0[k] - 4.0 * to_xyz_xx[k] * tke_0 - 2.0 * to_xyyyz_0[k] * tbe_0 + 4.0 * to_xyyyz_xx[k] * tbe_0 * tke_0;

            to_y_x_xyyz_y[k] = -4.0 * to_xyz_xy[k] * tke_0 + 4.0 * to_xyyyz_xy[k] * tbe_0 * tke_0;

            to_y_x_xyyz_z[k] = -4.0 * to_xyz_xz[k] * tke_0 + 4.0 * to_xyyyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 159-162 components of targeted buffer : GP

        auto to_y_x_xyzz_x = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 24);

        auto to_y_x_xyzz_y = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 25);

        auto to_y_x_xyzz_z = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 26);

#pragma omp simd aligned(to_xyyzz_0,        \
                             to_xyyzz_xx,   \
                             to_xyyzz_xy,   \
                             to_xyyzz_xz,   \
                             to_xzz_0,      \
                             to_xzz_xx,     \
                             to_xzz_xy,     \
                             to_xzz_xz,     \
                             to_y_x_xyzz_x, \
                             to_y_x_xyzz_y, \
                             to_y_x_xyzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyzz_x[k] = to_xzz_0[k] - 2.0 * to_xzz_xx[k] * tke_0 - 2.0 * to_xyyzz_0[k] * tbe_0 + 4.0 * to_xyyzz_xx[k] * tbe_0 * tke_0;

            to_y_x_xyzz_y[k] = -2.0 * to_xzz_xy[k] * tke_0 + 4.0 * to_xyyzz_xy[k] * tbe_0 * tke_0;

            to_y_x_xyzz_z[k] = -2.0 * to_xzz_xz[k] * tke_0 + 4.0 * to_xyyzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 162-165 components of targeted buffer : GP

        auto to_y_x_xzzz_x = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 27);

        auto to_y_x_xzzz_y = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 28);

        auto to_y_x_xzzz_z = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 29);

#pragma omp simd aligned(to_xyzzz_0, to_xyzzz_xx, to_xyzzz_xy, to_xyzzz_xz, to_y_x_xzzz_x, to_y_x_xzzz_y, to_y_x_xzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xzzz_x[k] = -2.0 * to_xyzzz_0[k] * tbe_0 + 4.0 * to_xyzzz_xx[k] * tbe_0 * tke_0;

            to_y_x_xzzz_y[k] = 4.0 * to_xyzzz_xy[k] * tbe_0 * tke_0;

            to_y_x_xzzz_z[k] = 4.0 * to_xyzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 165-168 components of targeted buffer : GP

        auto to_y_x_yyyy_x = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 30);

        auto to_y_x_yyyy_y = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 31);

        auto to_y_x_yyyy_z = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 32);

#pragma omp simd aligned(to_y_x_yyyy_x,     \
                             to_y_x_yyyy_y, \
                             to_y_x_yyyy_z, \
                             to_yyy_0,      \
                             to_yyy_xx,     \
                             to_yyy_xy,     \
                             to_yyy_xz,     \
                             to_yyyyy_0,    \
                             to_yyyyy_xx,   \
                             to_yyyyy_xy,   \
                             to_yyyyy_xz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyyy_x[k] = 4.0 * to_yyy_0[k] - 8.0 * to_yyy_xx[k] * tke_0 - 2.0 * to_yyyyy_0[k] * tbe_0 + 4.0 * to_yyyyy_xx[k] * tbe_0 * tke_0;

            to_y_x_yyyy_y[k] = -8.0 * to_yyy_xy[k] * tke_0 + 4.0 * to_yyyyy_xy[k] * tbe_0 * tke_0;

            to_y_x_yyyy_z[k] = -8.0 * to_yyy_xz[k] * tke_0 + 4.0 * to_yyyyy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 168-171 components of targeted buffer : GP

        auto to_y_x_yyyz_x = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 33);

        auto to_y_x_yyyz_y = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 34);

        auto to_y_x_yyyz_z = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 35);

#pragma omp simd aligned(to_y_x_yyyz_x,     \
                             to_y_x_yyyz_y, \
                             to_y_x_yyyz_z, \
                             to_yyyyz_0,    \
                             to_yyyyz_xx,   \
                             to_yyyyz_xy,   \
                             to_yyyyz_xz,   \
                             to_yyz_0,      \
                             to_yyz_xx,     \
                             to_yyz_xy,     \
                             to_yyz_xz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyyz_x[k] = 3.0 * to_yyz_0[k] - 6.0 * to_yyz_xx[k] * tke_0 - 2.0 * to_yyyyz_0[k] * tbe_0 + 4.0 * to_yyyyz_xx[k] * tbe_0 * tke_0;

            to_y_x_yyyz_y[k] = -6.0 * to_yyz_xy[k] * tke_0 + 4.0 * to_yyyyz_xy[k] * tbe_0 * tke_0;

            to_y_x_yyyz_z[k] = -6.0 * to_yyz_xz[k] * tke_0 + 4.0 * to_yyyyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 171-174 components of targeted buffer : GP

        auto to_y_x_yyzz_x = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 36);

        auto to_y_x_yyzz_y = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 37);

        auto to_y_x_yyzz_z = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 38);

#pragma omp simd aligned(to_y_x_yyzz_x,     \
                             to_y_x_yyzz_y, \
                             to_y_x_yyzz_z, \
                             to_yyyzz_0,    \
                             to_yyyzz_xx,   \
                             to_yyyzz_xy,   \
                             to_yyyzz_xz,   \
                             to_yzz_0,      \
                             to_yzz_xx,     \
                             to_yzz_xy,     \
                             to_yzz_xz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyzz_x[k] = 2.0 * to_yzz_0[k] - 4.0 * to_yzz_xx[k] * tke_0 - 2.0 * to_yyyzz_0[k] * tbe_0 + 4.0 * to_yyyzz_xx[k] * tbe_0 * tke_0;

            to_y_x_yyzz_y[k] = -4.0 * to_yzz_xy[k] * tke_0 + 4.0 * to_yyyzz_xy[k] * tbe_0 * tke_0;

            to_y_x_yyzz_z[k] = -4.0 * to_yzz_xz[k] * tke_0 + 4.0 * to_yyyzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 174-177 components of targeted buffer : GP

        auto to_y_x_yzzz_x = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 39);

        auto to_y_x_yzzz_y = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 40);

        auto to_y_x_yzzz_z = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 41);

#pragma omp simd aligned(to_y_x_yzzz_x,     \
                             to_y_x_yzzz_y, \
                             to_y_x_yzzz_z, \
                             to_yyzzz_0,    \
                             to_yyzzz_xx,   \
                             to_yyzzz_xy,   \
                             to_yyzzz_xz,   \
                             to_zzz_0,      \
                             to_zzz_xx,     \
                             to_zzz_xy,     \
                             to_zzz_xz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yzzz_x[k] = to_zzz_0[k] - 2.0 * to_zzz_xx[k] * tke_0 - 2.0 * to_yyzzz_0[k] * tbe_0 + 4.0 * to_yyzzz_xx[k] * tbe_0 * tke_0;

            to_y_x_yzzz_y[k] = -2.0 * to_zzz_xy[k] * tke_0 + 4.0 * to_yyzzz_xy[k] * tbe_0 * tke_0;

            to_y_x_yzzz_z[k] = -2.0 * to_zzz_xz[k] * tke_0 + 4.0 * to_yyzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 177-180 components of targeted buffer : GP

        auto to_y_x_zzzz_x = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 42);

        auto to_y_x_zzzz_y = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 43);

        auto to_y_x_zzzz_z = pbuffer.data(idx_op_geom_101_gp + 3 * op_comps * 45 + i * 45 + 44);

#pragma omp simd aligned(to_y_x_zzzz_x, to_y_x_zzzz_y, to_y_x_zzzz_z, to_yzzzz_0, to_yzzzz_xx, to_yzzzz_xy, to_yzzzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_zzzz_x[k] = -2.0 * to_yzzzz_0[k] * tbe_0 + 4.0 * to_yzzzz_xx[k] * tbe_0 * tke_0;

            to_y_x_zzzz_y[k] = 4.0 * to_yzzzz_xy[k] * tbe_0 * tke_0;

            to_y_x_zzzz_z[k] = 4.0 * to_yzzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 180-183 components of targeted buffer : GP

        auto to_y_y_xxxx_x = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 0);

        auto to_y_y_xxxx_y = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 1);

        auto to_y_y_xxxx_z = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 2);

#pragma omp simd aligned(to_xxxxy_0, to_xxxxy_xy, to_xxxxy_yy, to_xxxxy_yz, to_y_y_xxxx_x, to_y_y_xxxx_y, to_y_y_xxxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxxx_x[k] = 4.0 * to_xxxxy_xy[k] * tbe_0 * tke_0;

            to_y_y_xxxx_y[k] = -2.0 * to_xxxxy_0[k] * tbe_0 + 4.0 * to_xxxxy_yy[k] * tbe_0 * tke_0;

            to_y_y_xxxx_z[k] = 4.0 * to_xxxxy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 183-186 components of targeted buffer : GP

        auto to_y_y_xxxy_x = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 3);

        auto to_y_y_xxxy_y = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 4);

        auto to_y_y_xxxy_z = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 5);

#pragma omp simd aligned(to_xxx_0,          \
                             to_xxx_xy,     \
                             to_xxx_yy,     \
                             to_xxx_yz,     \
                             to_xxxyy_0,    \
                             to_xxxyy_xy,   \
                             to_xxxyy_yy,   \
                             to_xxxyy_yz,   \
                             to_y_y_xxxy_x, \
                             to_y_y_xxxy_y, \
                             to_y_y_xxxy_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxxy_x[k] = -2.0 * to_xxx_xy[k] * tke_0 + 4.0 * to_xxxyy_xy[k] * tbe_0 * tke_0;

            to_y_y_xxxy_y[k] = to_xxx_0[k] - 2.0 * to_xxx_yy[k] * tke_0 - 2.0 * to_xxxyy_0[k] * tbe_0 + 4.0 * to_xxxyy_yy[k] * tbe_0 * tke_0;

            to_y_y_xxxy_z[k] = -2.0 * to_xxx_yz[k] * tke_0 + 4.0 * to_xxxyy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 186-189 components of targeted buffer : GP

        auto to_y_y_xxxz_x = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 6);

        auto to_y_y_xxxz_y = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 7);

        auto to_y_y_xxxz_z = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 8);

#pragma omp simd aligned(to_xxxyz_0, to_xxxyz_xy, to_xxxyz_yy, to_xxxyz_yz, to_y_y_xxxz_x, to_y_y_xxxz_y, to_y_y_xxxz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxxz_x[k] = 4.0 * to_xxxyz_xy[k] * tbe_0 * tke_0;

            to_y_y_xxxz_y[k] = -2.0 * to_xxxyz_0[k] * tbe_0 + 4.0 * to_xxxyz_yy[k] * tbe_0 * tke_0;

            to_y_y_xxxz_z[k] = 4.0 * to_xxxyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 189-192 components of targeted buffer : GP

        auto to_y_y_xxyy_x = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 9);

        auto to_y_y_xxyy_y = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 10);

        auto to_y_y_xxyy_z = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 11);

#pragma omp simd aligned(to_xxy_0,          \
                             to_xxy_xy,     \
                             to_xxy_yy,     \
                             to_xxy_yz,     \
                             to_xxyyy_0,    \
                             to_xxyyy_xy,   \
                             to_xxyyy_yy,   \
                             to_xxyyy_yz,   \
                             to_y_y_xxyy_x, \
                             to_y_y_xxyy_y, \
                             to_y_y_xxyy_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxyy_x[k] = -4.0 * to_xxy_xy[k] * tke_0 + 4.0 * to_xxyyy_xy[k] * tbe_0 * tke_0;

            to_y_y_xxyy_y[k] = 2.0 * to_xxy_0[k] - 4.0 * to_xxy_yy[k] * tke_0 - 2.0 * to_xxyyy_0[k] * tbe_0 + 4.0 * to_xxyyy_yy[k] * tbe_0 * tke_0;

            to_y_y_xxyy_z[k] = -4.0 * to_xxy_yz[k] * tke_0 + 4.0 * to_xxyyy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 192-195 components of targeted buffer : GP

        auto to_y_y_xxyz_x = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 12);

        auto to_y_y_xxyz_y = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 13);

        auto to_y_y_xxyz_z = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 14);

#pragma omp simd aligned(to_xxyyz_0,        \
                             to_xxyyz_xy,   \
                             to_xxyyz_yy,   \
                             to_xxyyz_yz,   \
                             to_xxz_0,      \
                             to_xxz_xy,     \
                             to_xxz_yy,     \
                             to_xxz_yz,     \
                             to_y_y_xxyz_x, \
                             to_y_y_xxyz_y, \
                             to_y_y_xxyz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxyz_x[k] = -2.0 * to_xxz_xy[k] * tke_0 + 4.0 * to_xxyyz_xy[k] * tbe_0 * tke_0;

            to_y_y_xxyz_y[k] = to_xxz_0[k] - 2.0 * to_xxz_yy[k] * tke_0 - 2.0 * to_xxyyz_0[k] * tbe_0 + 4.0 * to_xxyyz_yy[k] * tbe_0 * tke_0;

            to_y_y_xxyz_z[k] = -2.0 * to_xxz_yz[k] * tke_0 + 4.0 * to_xxyyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 195-198 components of targeted buffer : GP

        auto to_y_y_xxzz_x = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 15);

        auto to_y_y_xxzz_y = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 16);

        auto to_y_y_xxzz_z = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 17);

#pragma omp simd aligned(to_xxyzz_0, to_xxyzz_xy, to_xxyzz_yy, to_xxyzz_yz, to_y_y_xxzz_x, to_y_y_xxzz_y, to_y_y_xxzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxzz_x[k] = 4.0 * to_xxyzz_xy[k] * tbe_0 * tke_0;

            to_y_y_xxzz_y[k] = -2.0 * to_xxyzz_0[k] * tbe_0 + 4.0 * to_xxyzz_yy[k] * tbe_0 * tke_0;

            to_y_y_xxzz_z[k] = 4.0 * to_xxyzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 198-201 components of targeted buffer : GP

        auto to_y_y_xyyy_x = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 18);

        auto to_y_y_xyyy_y = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 19);

        auto to_y_y_xyyy_z = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 20);

#pragma omp simd aligned(to_xyy_0,          \
                             to_xyy_xy,     \
                             to_xyy_yy,     \
                             to_xyy_yz,     \
                             to_xyyyy_0,    \
                             to_xyyyy_xy,   \
                             to_xyyyy_yy,   \
                             to_xyyyy_yz,   \
                             to_y_y_xyyy_x, \
                             to_y_y_xyyy_y, \
                             to_y_y_xyyy_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyyy_x[k] = -6.0 * to_xyy_xy[k] * tke_0 + 4.0 * to_xyyyy_xy[k] * tbe_0 * tke_0;

            to_y_y_xyyy_y[k] = 3.0 * to_xyy_0[k] - 6.0 * to_xyy_yy[k] * tke_0 - 2.0 * to_xyyyy_0[k] * tbe_0 + 4.0 * to_xyyyy_yy[k] * tbe_0 * tke_0;

            to_y_y_xyyy_z[k] = -6.0 * to_xyy_yz[k] * tke_0 + 4.0 * to_xyyyy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 201-204 components of targeted buffer : GP

        auto to_y_y_xyyz_x = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 21);

        auto to_y_y_xyyz_y = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 22);

        auto to_y_y_xyyz_z = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 23);

#pragma omp simd aligned(to_xyyyz_0,        \
                             to_xyyyz_xy,   \
                             to_xyyyz_yy,   \
                             to_xyyyz_yz,   \
                             to_xyz_0,      \
                             to_xyz_xy,     \
                             to_xyz_yy,     \
                             to_xyz_yz,     \
                             to_y_y_xyyz_x, \
                             to_y_y_xyyz_y, \
                             to_y_y_xyyz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyyz_x[k] = -4.0 * to_xyz_xy[k] * tke_0 + 4.0 * to_xyyyz_xy[k] * tbe_0 * tke_0;

            to_y_y_xyyz_y[k] = 2.0 * to_xyz_0[k] - 4.0 * to_xyz_yy[k] * tke_0 - 2.0 * to_xyyyz_0[k] * tbe_0 + 4.0 * to_xyyyz_yy[k] * tbe_0 * tke_0;

            to_y_y_xyyz_z[k] = -4.0 * to_xyz_yz[k] * tke_0 + 4.0 * to_xyyyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 204-207 components of targeted buffer : GP

        auto to_y_y_xyzz_x = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 24);

        auto to_y_y_xyzz_y = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 25);

        auto to_y_y_xyzz_z = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 26);

#pragma omp simd aligned(to_xyyzz_0,        \
                             to_xyyzz_xy,   \
                             to_xyyzz_yy,   \
                             to_xyyzz_yz,   \
                             to_xzz_0,      \
                             to_xzz_xy,     \
                             to_xzz_yy,     \
                             to_xzz_yz,     \
                             to_y_y_xyzz_x, \
                             to_y_y_xyzz_y, \
                             to_y_y_xyzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyzz_x[k] = -2.0 * to_xzz_xy[k] * tke_0 + 4.0 * to_xyyzz_xy[k] * tbe_0 * tke_0;

            to_y_y_xyzz_y[k] = to_xzz_0[k] - 2.0 * to_xzz_yy[k] * tke_0 - 2.0 * to_xyyzz_0[k] * tbe_0 + 4.0 * to_xyyzz_yy[k] * tbe_0 * tke_0;

            to_y_y_xyzz_z[k] = -2.0 * to_xzz_yz[k] * tke_0 + 4.0 * to_xyyzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 207-210 components of targeted buffer : GP

        auto to_y_y_xzzz_x = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 27);

        auto to_y_y_xzzz_y = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 28);

        auto to_y_y_xzzz_z = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 29);

#pragma omp simd aligned(to_xyzzz_0, to_xyzzz_xy, to_xyzzz_yy, to_xyzzz_yz, to_y_y_xzzz_x, to_y_y_xzzz_y, to_y_y_xzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xzzz_x[k] = 4.0 * to_xyzzz_xy[k] * tbe_0 * tke_0;

            to_y_y_xzzz_y[k] = -2.0 * to_xyzzz_0[k] * tbe_0 + 4.0 * to_xyzzz_yy[k] * tbe_0 * tke_0;

            to_y_y_xzzz_z[k] = 4.0 * to_xyzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 210-213 components of targeted buffer : GP

        auto to_y_y_yyyy_x = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 30);

        auto to_y_y_yyyy_y = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 31);

        auto to_y_y_yyyy_z = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 32);

#pragma omp simd aligned(to_y_y_yyyy_x,     \
                             to_y_y_yyyy_y, \
                             to_y_y_yyyy_z, \
                             to_yyy_0,      \
                             to_yyy_xy,     \
                             to_yyy_yy,     \
                             to_yyy_yz,     \
                             to_yyyyy_0,    \
                             to_yyyyy_xy,   \
                             to_yyyyy_yy,   \
                             to_yyyyy_yz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyyy_x[k] = -8.0 * to_yyy_xy[k] * tke_0 + 4.0 * to_yyyyy_xy[k] * tbe_0 * tke_0;

            to_y_y_yyyy_y[k] = 4.0 * to_yyy_0[k] - 8.0 * to_yyy_yy[k] * tke_0 - 2.0 * to_yyyyy_0[k] * tbe_0 + 4.0 * to_yyyyy_yy[k] * tbe_0 * tke_0;

            to_y_y_yyyy_z[k] = -8.0 * to_yyy_yz[k] * tke_0 + 4.0 * to_yyyyy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 213-216 components of targeted buffer : GP

        auto to_y_y_yyyz_x = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 33);

        auto to_y_y_yyyz_y = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 34);

        auto to_y_y_yyyz_z = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 35);

#pragma omp simd aligned(to_y_y_yyyz_x,     \
                             to_y_y_yyyz_y, \
                             to_y_y_yyyz_z, \
                             to_yyyyz_0,    \
                             to_yyyyz_xy,   \
                             to_yyyyz_yy,   \
                             to_yyyyz_yz,   \
                             to_yyz_0,      \
                             to_yyz_xy,     \
                             to_yyz_yy,     \
                             to_yyz_yz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyyz_x[k] = -6.0 * to_yyz_xy[k] * tke_0 + 4.0 * to_yyyyz_xy[k] * tbe_0 * tke_0;

            to_y_y_yyyz_y[k] = 3.0 * to_yyz_0[k] - 6.0 * to_yyz_yy[k] * tke_0 - 2.0 * to_yyyyz_0[k] * tbe_0 + 4.0 * to_yyyyz_yy[k] * tbe_0 * tke_0;

            to_y_y_yyyz_z[k] = -6.0 * to_yyz_yz[k] * tke_0 + 4.0 * to_yyyyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 216-219 components of targeted buffer : GP

        auto to_y_y_yyzz_x = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 36);

        auto to_y_y_yyzz_y = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 37);

        auto to_y_y_yyzz_z = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 38);

#pragma omp simd aligned(to_y_y_yyzz_x,     \
                             to_y_y_yyzz_y, \
                             to_y_y_yyzz_z, \
                             to_yyyzz_0,    \
                             to_yyyzz_xy,   \
                             to_yyyzz_yy,   \
                             to_yyyzz_yz,   \
                             to_yzz_0,      \
                             to_yzz_xy,     \
                             to_yzz_yy,     \
                             to_yzz_yz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyzz_x[k] = -4.0 * to_yzz_xy[k] * tke_0 + 4.0 * to_yyyzz_xy[k] * tbe_0 * tke_0;

            to_y_y_yyzz_y[k] = 2.0 * to_yzz_0[k] - 4.0 * to_yzz_yy[k] * tke_0 - 2.0 * to_yyyzz_0[k] * tbe_0 + 4.0 * to_yyyzz_yy[k] * tbe_0 * tke_0;

            to_y_y_yyzz_z[k] = -4.0 * to_yzz_yz[k] * tke_0 + 4.0 * to_yyyzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 219-222 components of targeted buffer : GP

        auto to_y_y_yzzz_x = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 39);

        auto to_y_y_yzzz_y = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 40);

        auto to_y_y_yzzz_z = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 41);

#pragma omp simd aligned(to_y_y_yzzz_x,     \
                             to_y_y_yzzz_y, \
                             to_y_y_yzzz_z, \
                             to_yyzzz_0,    \
                             to_yyzzz_xy,   \
                             to_yyzzz_yy,   \
                             to_yyzzz_yz,   \
                             to_zzz_0,      \
                             to_zzz_xy,     \
                             to_zzz_yy,     \
                             to_zzz_yz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yzzz_x[k] = -2.0 * to_zzz_xy[k] * tke_0 + 4.0 * to_yyzzz_xy[k] * tbe_0 * tke_0;

            to_y_y_yzzz_y[k] = to_zzz_0[k] - 2.0 * to_zzz_yy[k] * tke_0 - 2.0 * to_yyzzz_0[k] * tbe_0 + 4.0 * to_yyzzz_yy[k] * tbe_0 * tke_0;

            to_y_y_yzzz_z[k] = -2.0 * to_zzz_yz[k] * tke_0 + 4.0 * to_yyzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 222-225 components of targeted buffer : GP

        auto to_y_y_zzzz_x = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 42);

        auto to_y_y_zzzz_y = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 43);

        auto to_y_y_zzzz_z = pbuffer.data(idx_op_geom_101_gp + 4 * op_comps * 45 + i * 45 + 44);

#pragma omp simd aligned(to_y_y_zzzz_x, to_y_y_zzzz_y, to_y_y_zzzz_z, to_yzzzz_0, to_yzzzz_xy, to_yzzzz_yy, to_yzzzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_zzzz_x[k] = 4.0 * to_yzzzz_xy[k] * tbe_0 * tke_0;

            to_y_y_zzzz_y[k] = -2.0 * to_yzzzz_0[k] * tbe_0 + 4.0 * to_yzzzz_yy[k] * tbe_0 * tke_0;

            to_y_y_zzzz_z[k] = 4.0 * to_yzzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 225-228 components of targeted buffer : GP

        auto to_y_z_xxxx_x = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 0);

        auto to_y_z_xxxx_y = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 1);

        auto to_y_z_xxxx_z = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 2);

#pragma omp simd aligned(to_xxxxy_0, to_xxxxy_xz, to_xxxxy_yz, to_xxxxy_zz, to_y_z_xxxx_x, to_y_z_xxxx_y, to_y_z_xxxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxxx_x[k] = 4.0 * to_xxxxy_xz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_y[k] = 4.0 * to_xxxxy_yz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_z[k] = -2.0 * to_xxxxy_0[k] * tbe_0 + 4.0 * to_xxxxy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 228-231 components of targeted buffer : GP

        auto to_y_z_xxxy_x = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 3);

        auto to_y_z_xxxy_y = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 4);

        auto to_y_z_xxxy_z = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 5);

#pragma omp simd aligned(to_xxx_0,          \
                             to_xxx_xz,     \
                             to_xxx_yz,     \
                             to_xxx_zz,     \
                             to_xxxyy_0,    \
                             to_xxxyy_xz,   \
                             to_xxxyy_yz,   \
                             to_xxxyy_zz,   \
                             to_y_z_xxxy_x, \
                             to_y_z_xxxy_y, \
                             to_y_z_xxxy_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxxy_x[k] = -2.0 * to_xxx_xz[k] * tke_0 + 4.0 * to_xxxyy_xz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_y[k] = -2.0 * to_xxx_yz[k] * tke_0 + 4.0 * to_xxxyy_yz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_z[k] = to_xxx_0[k] - 2.0 * to_xxx_zz[k] * tke_0 - 2.0 * to_xxxyy_0[k] * tbe_0 + 4.0 * to_xxxyy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 231-234 components of targeted buffer : GP

        auto to_y_z_xxxz_x = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 6);

        auto to_y_z_xxxz_y = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 7);

        auto to_y_z_xxxz_z = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 8);

#pragma omp simd aligned(to_xxxyz_0, to_xxxyz_xz, to_xxxyz_yz, to_xxxyz_zz, to_y_z_xxxz_x, to_y_z_xxxz_y, to_y_z_xxxz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxxz_x[k] = 4.0 * to_xxxyz_xz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_y[k] = 4.0 * to_xxxyz_yz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_z[k] = -2.0 * to_xxxyz_0[k] * tbe_0 + 4.0 * to_xxxyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 234-237 components of targeted buffer : GP

        auto to_y_z_xxyy_x = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 9);

        auto to_y_z_xxyy_y = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 10);

        auto to_y_z_xxyy_z = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 11);

#pragma omp simd aligned(to_xxy_0,          \
                             to_xxy_xz,     \
                             to_xxy_yz,     \
                             to_xxy_zz,     \
                             to_xxyyy_0,    \
                             to_xxyyy_xz,   \
                             to_xxyyy_yz,   \
                             to_xxyyy_zz,   \
                             to_y_z_xxyy_x, \
                             to_y_z_xxyy_y, \
                             to_y_z_xxyy_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxyy_x[k] = -4.0 * to_xxy_xz[k] * tke_0 + 4.0 * to_xxyyy_xz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_y[k] = -4.0 * to_xxy_yz[k] * tke_0 + 4.0 * to_xxyyy_yz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_z[k] = 2.0 * to_xxy_0[k] - 4.0 * to_xxy_zz[k] * tke_0 - 2.0 * to_xxyyy_0[k] * tbe_0 + 4.0 * to_xxyyy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 237-240 components of targeted buffer : GP

        auto to_y_z_xxyz_x = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 12);

        auto to_y_z_xxyz_y = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 13);

        auto to_y_z_xxyz_z = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 14);

#pragma omp simd aligned(to_xxyyz_0,        \
                             to_xxyyz_xz,   \
                             to_xxyyz_yz,   \
                             to_xxyyz_zz,   \
                             to_xxz_0,      \
                             to_xxz_xz,     \
                             to_xxz_yz,     \
                             to_xxz_zz,     \
                             to_y_z_xxyz_x, \
                             to_y_z_xxyz_y, \
                             to_y_z_xxyz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxyz_x[k] = -2.0 * to_xxz_xz[k] * tke_0 + 4.0 * to_xxyyz_xz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_y[k] = -2.0 * to_xxz_yz[k] * tke_0 + 4.0 * to_xxyyz_yz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_z[k] = to_xxz_0[k] - 2.0 * to_xxz_zz[k] * tke_0 - 2.0 * to_xxyyz_0[k] * tbe_0 + 4.0 * to_xxyyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 240-243 components of targeted buffer : GP

        auto to_y_z_xxzz_x = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 15);

        auto to_y_z_xxzz_y = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 16);

        auto to_y_z_xxzz_z = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 17);

#pragma omp simd aligned(to_xxyzz_0, to_xxyzz_xz, to_xxyzz_yz, to_xxyzz_zz, to_y_z_xxzz_x, to_y_z_xxzz_y, to_y_z_xxzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxzz_x[k] = 4.0 * to_xxyzz_xz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_y[k] = 4.0 * to_xxyzz_yz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_z[k] = -2.0 * to_xxyzz_0[k] * tbe_0 + 4.0 * to_xxyzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 243-246 components of targeted buffer : GP

        auto to_y_z_xyyy_x = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 18);

        auto to_y_z_xyyy_y = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 19);

        auto to_y_z_xyyy_z = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 20);

#pragma omp simd aligned(to_xyy_0,          \
                             to_xyy_xz,     \
                             to_xyy_yz,     \
                             to_xyy_zz,     \
                             to_xyyyy_0,    \
                             to_xyyyy_xz,   \
                             to_xyyyy_yz,   \
                             to_xyyyy_zz,   \
                             to_y_z_xyyy_x, \
                             to_y_z_xyyy_y, \
                             to_y_z_xyyy_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyyy_x[k] = -6.0 * to_xyy_xz[k] * tke_0 + 4.0 * to_xyyyy_xz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_y[k] = -6.0 * to_xyy_yz[k] * tke_0 + 4.0 * to_xyyyy_yz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_z[k] = 3.0 * to_xyy_0[k] - 6.0 * to_xyy_zz[k] * tke_0 - 2.0 * to_xyyyy_0[k] * tbe_0 + 4.0 * to_xyyyy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 246-249 components of targeted buffer : GP

        auto to_y_z_xyyz_x = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 21);

        auto to_y_z_xyyz_y = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 22);

        auto to_y_z_xyyz_z = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 23);

#pragma omp simd aligned(to_xyyyz_0,        \
                             to_xyyyz_xz,   \
                             to_xyyyz_yz,   \
                             to_xyyyz_zz,   \
                             to_xyz_0,      \
                             to_xyz_xz,     \
                             to_xyz_yz,     \
                             to_xyz_zz,     \
                             to_y_z_xyyz_x, \
                             to_y_z_xyyz_y, \
                             to_y_z_xyyz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyyz_x[k] = -4.0 * to_xyz_xz[k] * tke_0 + 4.0 * to_xyyyz_xz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_y[k] = -4.0 * to_xyz_yz[k] * tke_0 + 4.0 * to_xyyyz_yz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_z[k] = 2.0 * to_xyz_0[k] - 4.0 * to_xyz_zz[k] * tke_0 - 2.0 * to_xyyyz_0[k] * tbe_0 + 4.0 * to_xyyyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 249-252 components of targeted buffer : GP

        auto to_y_z_xyzz_x = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 24);

        auto to_y_z_xyzz_y = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 25);

        auto to_y_z_xyzz_z = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 26);

#pragma omp simd aligned(to_xyyzz_0,        \
                             to_xyyzz_xz,   \
                             to_xyyzz_yz,   \
                             to_xyyzz_zz,   \
                             to_xzz_0,      \
                             to_xzz_xz,     \
                             to_xzz_yz,     \
                             to_xzz_zz,     \
                             to_y_z_xyzz_x, \
                             to_y_z_xyzz_y, \
                             to_y_z_xyzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyzz_x[k] = -2.0 * to_xzz_xz[k] * tke_0 + 4.0 * to_xyyzz_xz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_y[k] = -2.0 * to_xzz_yz[k] * tke_0 + 4.0 * to_xyyzz_yz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_z[k] = to_xzz_0[k] - 2.0 * to_xzz_zz[k] * tke_0 - 2.0 * to_xyyzz_0[k] * tbe_0 + 4.0 * to_xyyzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 252-255 components of targeted buffer : GP

        auto to_y_z_xzzz_x = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 27);

        auto to_y_z_xzzz_y = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 28);

        auto to_y_z_xzzz_z = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 29);

#pragma omp simd aligned(to_xyzzz_0, to_xyzzz_xz, to_xyzzz_yz, to_xyzzz_zz, to_y_z_xzzz_x, to_y_z_xzzz_y, to_y_z_xzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xzzz_x[k] = 4.0 * to_xyzzz_xz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_y[k] = 4.0 * to_xyzzz_yz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_z[k] = -2.0 * to_xyzzz_0[k] * tbe_0 + 4.0 * to_xyzzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 255-258 components of targeted buffer : GP

        auto to_y_z_yyyy_x = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 30);

        auto to_y_z_yyyy_y = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 31);

        auto to_y_z_yyyy_z = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 32);

#pragma omp simd aligned(to_y_z_yyyy_x,     \
                             to_y_z_yyyy_y, \
                             to_y_z_yyyy_z, \
                             to_yyy_0,      \
                             to_yyy_xz,     \
                             to_yyy_yz,     \
                             to_yyy_zz,     \
                             to_yyyyy_0,    \
                             to_yyyyy_xz,   \
                             to_yyyyy_yz,   \
                             to_yyyyy_zz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyyy_x[k] = -8.0 * to_yyy_xz[k] * tke_0 + 4.0 * to_yyyyy_xz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_y[k] = -8.0 * to_yyy_yz[k] * tke_0 + 4.0 * to_yyyyy_yz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_z[k] = 4.0 * to_yyy_0[k] - 8.0 * to_yyy_zz[k] * tke_0 - 2.0 * to_yyyyy_0[k] * tbe_0 + 4.0 * to_yyyyy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 258-261 components of targeted buffer : GP

        auto to_y_z_yyyz_x = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 33);

        auto to_y_z_yyyz_y = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 34);

        auto to_y_z_yyyz_z = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 35);

#pragma omp simd aligned(to_y_z_yyyz_x,     \
                             to_y_z_yyyz_y, \
                             to_y_z_yyyz_z, \
                             to_yyyyz_0,    \
                             to_yyyyz_xz,   \
                             to_yyyyz_yz,   \
                             to_yyyyz_zz,   \
                             to_yyz_0,      \
                             to_yyz_xz,     \
                             to_yyz_yz,     \
                             to_yyz_zz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyyz_x[k] = -6.0 * to_yyz_xz[k] * tke_0 + 4.0 * to_yyyyz_xz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_y[k] = -6.0 * to_yyz_yz[k] * tke_0 + 4.0 * to_yyyyz_yz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_z[k] = 3.0 * to_yyz_0[k] - 6.0 * to_yyz_zz[k] * tke_0 - 2.0 * to_yyyyz_0[k] * tbe_0 + 4.0 * to_yyyyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 261-264 components of targeted buffer : GP

        auto to_y_z_yyzz_x = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 36);

        auto to_y_z_yyzz_y = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 37);

        auto to_y_z_yyzz_z = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 38);

#pragma omp simd aligned(to_y_z_yyzz_x,     \
                             to_y_z_yyzz_y, \
                             to_y_z_yyzz_z, \
                             to_yyyzz_0,    \
                             to_yyyzz_xz,   \
                             to_yyyzz_yz,   \
                             to_yyyzz_zz,   \
                             to_yzz_0,      \
                             to_yzz_xz,     \
                             to_yzz_yz,     \
                             to_yzz_zz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyzz_x[k] = -4.0 * to_yzz_xz[k] * tke_0 + 4.0 * to_yyyzz_xz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_y[k] = -4.0 * to_yzz_yz[k] * tke_0 + 4.0 * to_yyyzz_yz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_z[k] = 2.0 * to_yzz_0[k] - 4.0 * to_yzz_zz[k] * tke_0 - 2.0 * to_yyyzz_0[k] * tbe_0 + 4.0 * to_yyyzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 264-267 components of targeted buffer : GP

        auto to_y_z_yzzz_x = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 39);

        auto to_y_z_yzzz_y = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 40);

        auto to_y_z_yzzz_z = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 41);

#pragma omp simd aligned(to_y_z_yzzz_x,     \
                             to_y_z_yzzz_y, \
                             to_y_z_yzzz_z, \
                             to_yyzzz_0,    \
                             to_yyzzz_xz,   \
                             to_yyzzz_yz,   \
                             to_yyzzz_zz,   \
                             to_zzz_0,      \
                             to_zzz_xz,     \
                             to_zzz_yz,     \
                             to_zzz_zz,     \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yzzz_x[k] = -2.0 * to_zzz_xz[k] * tke_0 + 4.0 * to_yyzzz_xz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_y[k] = -2.0 * to_zzz_yz[k] * tke_0 + 4.0 * to_yyzzz_yz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_z[k] = to_zzz_0[k] - 2.0 * to_zzz_zz[k] * tke_0 - 2.0 * to_yyzzz_0[k] * tbe_0 + 4.0 * to_yyzzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 267-270 components of targeted buffer : GP

        auto to_y_z_zzzz_x = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 42);

        auto to_y_z_zzzz_y = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 43);

        auto to_y_z_zzzz_z = pbuffer.data(idx_op_geom_101_gp + 5 * op_comps * 45 + i * 45 + 44);

#pragma omp simd aligned(to_y_z_zzzz_x, to_y_z_zzzz_y, to_y_z_zzzz_z, to_yzzzz_0, to_yzzzz_xz, to_yzzzz_yz, to_yzzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_zzzz_x[k] = 4.0 * to_yzzzz_xz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_y[k] = 4.0 * to_yzzzz_yz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_z[k] = -2.0 * to_yzzzz_0[k] * tbe_0 + 4.0 * to_yzzzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 270-273 components of targeted buffer : GP

        auto to_z_x_xxxx_x = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 0);

        auto to_z_x_xxxx_y = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 1);

        auto to_z_x_xxxx_z = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 2);

#pragma omp simd aligned(to_xxxxz_0, to_xxxxz_xx, to_xxxxz_xy, to_xxxxz_xz, to_z_x_xxxx_x, to_z_x_xxxx_y, to_z_x_xxxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxxx_x[k] = -2.0 * to_xxxxz_0[k] * tbe_0 + 4.0 * to_xxxxz_xx[k] * tbe_0 * tke_0;

            to_z_x_xxxx_y[k] = 4.0 * to_xxxxz_xy[k] * tbe_0 * tke_0;

            to_z_x_xxxx_z[k] = 4.0 * to_xxxxz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 273-276 components of targeted buffer : GP

        auto to_z_x_xxxy_x = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 3);

        auto to_z_x_xxxy_y = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 4);

        auto to_z_x_xxxy_z = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 5);

#pragma omp simd aligned(to_xxxyz_0, to_xxxyz_xx, to_xxxyz_xy, to_xxxyz_xz, to_z_x_xxxy_x, to_z_x_xxxy_y, to_z_x_xxxy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxxy_x[k] = -2.0 * to_xxxyz_0[k] * tbe_0 + 4.0 * to_xxxyz_xx[k] * tbe_0 * tke_0;

            to_z_x_xxxy_y[k] = 4.0 * to_xxxyz_xy[k] * tbe_0 * tke_0;

            to_z_x_xxxy_z[k] = 4.0 * to_xxxyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 276-279 components of targeted buffer : GP

        auto to_z_x_xxxz_x = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 6);

        auto to_z_x_xxxz_y = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 7);

        auto to_z_x_xxxz_z = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 8);

#pragma omp simd aligned(to_xxx_0,          \
                             to_xxx_xx,     \
                             to_xxx_xy,     \
                             to_xxx_xz,     \
                             to_xxxzz_0,    \
                             to_xxxzz_xx,   \
                             to_xxxzz_xy,   \
                             to_xxxzz_xz,   \
                             to_z_x_xxxz_x, \
                             to_z_x_xxxz_y, \
                             to_z_x_xxxz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxxz_x[k] = to_xxx_0[k] - 2.0 * to_xxx_xx[k] * tke_0 - 2.0 * to_xxxzz_0[k] * tbe_0 + 4.0 * to_xxxzz_xx[k] * tbe_0 * tke_0;

            to_z_x_xxxz_y[k] = -2.0 * to_xxx_xy[k] * tke_0 + 4.0 * to_xxxzz_xy[k] * tbe_0 * tke_0;

            to_z_x_xxxz_z[k] = -2.0 * to_xxx_xz[k] * tke_0 + 4.0 * to_xxxzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 279-282 components of targeted buffer : GP

        auto to_z_x_xxyy_x = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 9);

        auto to_z_x_xxyy_y = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 10);

        auto to_z_x_xxyy_z = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 11);

#pragma omp simd aligned(to_xxyyz_0, to_xxyyz_xx, to_xxyyz_xy, to_xxyyz_xz, to_z_x_xxyy_x, to_z_x_xxyy_y, to_z_x_xxyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxyy_x[k] = -2.0 * to_xxyyz_0[k] * tbe_0 + 4.0 * to_xxyyz_xx[k] * tbe_0 * tke_0;

            to_z_x_xxyy_y[k] = 4.0 * to_xxyyz_xy[k] * tbe_0 * tke_0;

            to_z_x_xxyy_z[k] = 4.0 * to_xxyyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 282-285 components of targeted buffer : GP

        auto to_z_x_xxyz_x = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 12);

        auto to_z_x_xxyz_y = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 13);

        auto to_z_x_xxyz_z = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 14);

#pragma omp simd aligned(to_xxy_0,          \
                             to_xxy_xx,     \
                             to_xxy_xy,     \
                             to_xxy_xz,     \
                             to_xxyzz_0,    \
                             to_xxyzz_xx,   \
                             to_xxyzz_xy,   \
                             to_xxyzz_xz,   \
                             to_z_x_xxyz_x, \
                             to_z_x_xxyz_y, \
                             to_z_x_xxyz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxyz_x[k] = to_xxy_0[k] - 2.0 * to_xxy_xx[k] * tke_0 - 2.0 * to_xxyzz_0[k] * tbe_0 + 4.0 * to_xxyzz_xx[k] * tbe_0 * tke_0;

            to_z_x_xxyz_y[k] = -2.0 * to_xxy_xy[k] * tke_0 + 4.0 * to_xxyzz_xy[k] * tbe_0 * tke_0;

            to_z_x_xxyz_z[k] = -2.0 * to_xxy_xz[k] * tke_0 + 4.0 * to_xxyzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 285-288 components of targeted buffer : GP

        auto to_z_x_xxzz_x = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 15);

        auto to_z_x_xxzz_y = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 16);

        auto to_z_x_xxzz_z = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 17);

#pragma omp simd aligned(to_xxz_0,          \
                             to_xxz_xx,     \
                             to_xxz_xy,     \
                             to_xxz_xz,     \
                             to_xxzzz_0,    \
                             to_xxzzz_xx,   \
                             to_xxzzz_xy,   \
                             to_xxzzz_xz,   \
                             to_z_x_xxzz_x, \
                             to_z_x_xxzz_y, \
                             to_z_x_xxzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxzz_x[k] = 2.0 * to_xxz_0[k] - 4.0 * to_xxz_xx[k] * tke_0 - 2.0 * to_xxzzz_0[k] * tbe_0 + 4.0 * to_xxzzz_xx[k] * tbe_0 * tke_0;

            to_z_x_xxzz_y[k] = -4.0 * to_xxz_xy[k] * tke_0 + 4.0 * to_xxzzz_xy[k] * tbe_0 * tke_0;

            to_z_x_xxzz_z[k] = -4.0 * to_xxz_xz[k] * tke_0 + 4.0 * to_xxzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 288-291 components of targeted buffer : GP

        auto to_z_x_xyyy_x = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 18);

        auto to_z_x_xyyy_y = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 19);

        auto to_z_x_xyyy_z = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 20);

#pragma omp simd aligned(to_xyyyz_0, to_xyyyz_xx, to_xyyyz_xy, to_xyyyz_xz, to_z_x_xyyy_x, to_z_x_xyyy_y, to_z_x_xyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyyy_x[k] = -2.0 * to_xyyyz_0[k] * tbe_0 + 4.0 * to_xyyyz_xx[k] * tbe_0 * tke_0;

            to_z_x_xyyy_y[k] = 4.0 * to_xyyyz_xy[k] * tbe_0 * tke_0;

            to_z_x_xyyy_z[k] = 4.0 * to_xyyyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 291-294 components of targeted buffer : GP

        auto to_z_x_xyyz_x = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 21);

        auto to_z_x_xyyz_y = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 22);

        auto to_z_x_xyyz_z = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 23);

#pragma omp simd aligned(to_xyy_0,          \
                             to_xyy_xx,     \
                             to_xyy_xy,     \
                             to_xyy_xz,     \
                             to_xyyzz_0,    \
                             to_xyyzz_xx,   \
                             to_xyyzz_xy,   \
                             to_xyyzz_xz,   \
                             to_z_x_xyyz_x, \
                             to_z_x_xyyz_y, \
                             to_z_x_xyyz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyyz_x[k] = to_xyy_0[k] - 2.0 * to_xyy_xx[k] * tke_0 - 2.0 * to_xyyzz_0[k] * tbe_0 + 4.0 * to_xyyzz_xx[k] * tbe_0 * tke_0;

            to_z_x_xyyz_y[k] = -2.0 * to_xyy_xy[k] * tke_0 + 4.0 * to_xyyzz_xy[k] * tbe_0 * tke_0;

            to_z_x_xyyz_z[k] = -2.0 * to_xyy_xz[k] * tke_0 + 4.0 * to_xyyzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 294-297 components of targeted buffer : GP

        auto to_z_x_xyzz_x = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 24);

        auto to_z_x_xyzz_y = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 25);

        auto to_z_x_xyzz_z = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 26);

#pragma omp simd aligned(to_xyz_0,          \
                             to_xyz_xx,     \
                             to_xyz_xy,     \
                             to_xyz_xz,     \
                             to_xyzzz_0,    \
                             to_xyzzz_xx,   \
                             to_xyzzz_xy,   \
                             to_xyzzz_xz,   \
                             to_z_x_xyzz_x, \
                             to_z_x_xyzz_y, \
                             to_z_x_xyzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyzz_x[k] = 2.0 * to_xyz_0[k] - 4.0 * to_xyz_xx[k] * tke_0 - 2.0 * to_xyzzz_0[k] * tbe_0 + 4.0 * to_xyzzz_xx[k] * tbe_0 * tke_0;

            to_z_x_xyzz_y[k] = -4.0 * to_xyz_xy[k] * tke_0 + 4.0 * to_xyzzz_xy[k] * tbe_0 * tke_0;

            to_z_x_xyzz_z[k] = -4.0 * to_xyz_xz[k] * tke_0 + 4.0 * to_xyzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 297-300 components of targeted buffer : GP

        auto to_z_x_xzzz_x = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 27);

        auto to_z_x_xzzz_y = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 28);

        auto to_z_x_xzzz_z = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 29);

#pragma omp simd aligned(to_xzz_0,          \
                             to_xzz_xx,     \
                             to_xzz_xy,     \
                             to_xzz_xz,     \
                             to_xzzzz_0,    \
                             to_xzzzz_xx,   \
                             to_xzzzz_xy,   \
                             to_xzzzz_xz,   \
                             to_z_x_xzzz_x, \
                             to_z_x_xzzz_y, \
                             to_z_x_xzzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xzzz_x[k] = 3.0 * to_xzz_0[k] - 6.0 * to_xzz_xx[k] * tke_0 - 2.0 * to_xzzzz_0[k] * tbe_0 + 4.0 * to_xzzzz_xx[k] * tbe_0 * tke_0;

            to_z_x_xzzz_y[k] = -6.0 * to_xzz_xy[k] * tke_0 + 4.0 * to_xzzzz_xy[k] * tbe_0 * tke_0;

            to_z_x_xzzz_z[k] = -6.0 * to_xzz_xz[k] * tke_0 + 4.0 * to_xzzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 300-303 components of targeted buffer : GP

        auto to_z_x_yyyy_x = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 30);

        auto to_z_x_yyyy_y = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 31);

        auto to_z_x_yyyy_z = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 32);

#pragma omp simd aligned(to_yyyyz_0, to_yyyyz_xx, to_yyyyz_xy, to_yyyyz_xz, to_z_x_yyyy_x, to_z_x_yyyy_y, to_z_x_yyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyyy_x[k] = -2.0 * to_yyyyz_0[k] * tbe_0 + 4.0 * to_yyyyz_xx[k] * tbe_0 * tke_0;

            to_z_x_yyyy_y[k] = 4.0 * to_yyyyz_xy[k] * tbe_0 * tke_0;

            to_z_x_yyyy_z[k] = 4.0 * to_yyyyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 303-306 components of targeted buffer : GP

        auto to_z_x_yyyz_x = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 33);

        auto to_z_x_yyyz_y = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 34);

        auto to_z_x_yyyz_z = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 35);

#pragma omp simd aligned(to_yyy_0,          \
                             to_yyy_xx,     \
                             to_yyy_xy,     \
                             to_yyy_xz,     \
                             to_yyyzz_0,    \
                             to_yyyzz_xx,   \
                             to_yyyzz_xy,   \
                             to_yyyzz_xz,   \
                             to_z_x_yyyz_x, \
                             to_z_x_yyyz_y, \
                             to_z_x_yyyz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyyz_x[k] = to_yyy_0[k] - 2.0 * to_yyy_xx[k] * tke_0 - 2.0 * to_yyyzz_0[k] * tbe_0 + 4.0 * to_yyyzz_xx[k] * tbe_0 * tke_0;

            to_z_x_yyyz_y[k] = -2.0 * to_yyy_xy[k] * tke_0 + 4.0 * to_yyyzz_xy[k] * tbe_0 * tke_0;

            to_z_x_yyyz_z[k] = -2.0 * to_yyy_xz[k] * tke_0 + 4.0 * to_yyyzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 306-309 components of targeted buffer : GP

        auto to_z_x_yyzz_x = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 36);

        auto to_z_x_yyzz_y = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 37);

        auto to_z_x_yyzz_z = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 38);

#pragma omp simd aligned(to_yyz_0,          \
                             to_yyz_xx,     \
                             to_yyz_xy,     \
                             to_yyz_xz,     \
                             to_yyzzz_0,    \
                             to_yyzzz_xx,   \
                             to_yyzzz_xy,   \
                             to_yyzzz_xz,   \
                             to_z_x_yyzz_x, \
                             to_z_x_yyzz_y, \
                             to_z_x_yyzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyzz_x[k] = 2.0 * to_yyz_0[k] - 4.0 * to_yyz_xx[k] * tke_0 - 2.0 * to_yyzzz_0[k] * tbe_0 + 4.0 * to_yyzzz_xx[k] * tbe_0 * tke_0;

            to_z_x_yyzz_y[k] = -4.0 * to_yyz_xy[k] * tke_0 + 4.0 * to_yyzzz_xy[k] * tbe_0 * tke_0;

            to_z_x_yyzz_z[k] = -4.0 * to_yyz_xz[k] * tke_0 + 4.0 * to_yyzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 309-312 components of targeted buffer : GP

        auto to_z_x_yzzz_x = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 39);

        auto to_z_x_yzzz_y = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 40);

        auto to_z_x_yzzz_z = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 41);

#pragma omp simd aligned(to_yzz_0,          \
                             to_yzz_xx,     \
                             to_yzz_xy,     \
                             to_yzz_xz,     \
                             to_yzzzz_0,    \
                             to_yzzzz_xx,   \
                             to_yzzzz_xy,   \
                             to_yzzzz_xz,   \
                             to_z_x_yzzz_x, \
                             to_z_x_yzzz_y, \
                             to_z_x_yzzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yzzz_x[k] = 3.0 * to_yzz_0[k] - 6.0 * to_yzz_xx[k] * tke_0 - 2.0 * to_yzzzz_0[k] * tbe_0 + 4.0 * to_yzzzz_xx[k] * tbe_0 * tke_0;

            to_z_x_yzzz_y[k] = -6.0 * to_yzz_xy[k] * tke_0 + 4.0 * to_yzzzz_xy[k] * tbe_0 * tke_0;

            to_z_x_yzzz_z[k] = -6.0 * to_yzz_xz[k] * tke_0 + 4.0 * to_yzzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 312-315 components of targeted buffer : GP

        auto to_z_x_zzzz_x = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 42);

        auto to_z_x_zzzz_y = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 43);

        auto to_z_x_zzzz_z = pbuffer.data(idx_op_geom_101_gp + 6 * op_comps * 45 + i * 45 + 44);

#pragma omp simd aligned(to_z_x_zzzz_x,     \
                             to_z_x_zzzz_y, \
                             to_z_x_zzzz_z, \
                             to_zzz_0,      \
                             to_zzz_xx,     \
                             to_zzz_xy,     \
                             to_zzz_xz,     \
                             to_zzzzz_0,    \
                             to_zzzzz_xx,   \
                             to_zzzzz_xy,   \
                             to_zzzzz_xz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_zzzz_x[k] = 4.0 * to_zzz_0[k] - 8.0 * to_zzz_xx[k] * tke_0 - 2.0 * to_zzzzz_0[k] * tbe_0 + 4.0 * to_zzzzz_xx[k] * tbe_0 * tke_0;

            to_z_x_zzzz_y[k] = -8.0 * to_zzz_xy[k] * tke_0 + 4.0 * to_zzzzz_xy[k] * tbe_0 * tke_0;

            to_z_x_zzzz_z[k] = -8.0 * to_zzz_xz[k] * tke_0 + 4.0 * to_zzzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 315-318 components of targeted buffer : GP

        auto to_z_y_xxxx_x = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 0);

        auto to_z_y_xxxx_y = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 1);

        auto to_z_y_xxxx_z = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 2);

#pragma omp simd aligned(to_xxxxz_0, to_xxxxz_xy, to_xxxxz_yy, to_xxxxz_yz, to_z_y_xxxx_x, to_z_y_xxxx_y, to_z_y_xxxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxxx_x[k] = 4.0 * to_xxxxz_xy[k] * tbe_0 * tke_0;

            to_z_y_xxxx_y[k] = -2.0 * to_xxxxz_0[k] * tbe_0 + 4.0 * to_xxxxz_yy[k] * tbe_0 * tke_0;

            to_z_y_xxxx_z[k] = 4.0 * to_xxxxz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 318-321 components of targeted buffer : GP

        auto to_z_y_xxxy_x = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 3);

        auto to_z_y_xxxy_y = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 4);

        auto to_z_y_xxxy_z = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 5);

#pragma omp simd aligned(to_xxxyz_0, to_xxxyz_xy, to_xxxyz_yy, to_xxxyz_yz, to_z_y_xxxy_x, to_z_y_xxxy_y, to_z_y_xxxy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxxy_x[k] = 4.0 * to_xxxyz_xy[k] * tbe_0 * tke_0;

            to_z_y_xxxy_y[k] = -2.0 * to_xxxyz_0[k] * tbe_0 + 4.0 * to_xxxyz_yy[k] * tbe_0 * tke_0;

            to_z_y_xxxy_z[k] = 4.0 * to_xxxyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 321-324 components of targeted buffer : GP

        auto to_z_y_xxxz_x = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 6);

        auto to_z_y_xxxz_y = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 7);

        auto to_z_y_xxxz_z = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 8);

#pragma omp simd aligned(to_xxx_0,          \
                             to_xxx_xy,     \
                             to_xxx_yy,     \
                             to_xxx_yz,     \
                             to_xxxzz_0,    \
                             to_xxxzz_xy,   \
                             to_xxxzz_yy,   \
                             to_xxxzz_yz,   \
                             to_z_y_xxxz_x, \
                             to_z_y_xxxz_y, \
                             to_z_y_xxxz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxxz_x[k] = -2.0 * to_xxx_xy[k] * tke_0 + 4.0 * to_xxxzz_xy[k] * tbe_0 * tke_0;

            to_z_y_xxxz_y[k] = to_xxx_0[k] - 2.0 * to_xxx_yy[k] * tke_0 - 2.0 * to_xxxzz_0[k] * tbe_0 + 4.0 * to_xxxzz_yy[k] * tbe_0 * tke_0;

            to_z_y_xxxz_z[k] = -2.0 * to_xxx_yz[k] * tke_0 + 4.0 * to_xxxzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 324-327 components of targeted buffer : GP

        auto to_z_y_xxyy_x = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 9);

        auto to_z_y_xxyy_y = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 10);

        auto to_z_y_xxyy_z = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 11);

#pragma omp simd aligned(to_xxyyz_0, to_xxyyz_xy, to_xxyyz_yy, to_xxyyz_yz, to_z_y_xxyy_x, to_z_y_xxyy_y, to_z_y_xxyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxyy_x[k] = 4.0 * to_xxyyz_xy[k] * tbe_0 * tke_0;

            to_z_y_xxyy_y[k] = -2.0 * to_xxyyz_0[k] * tbe_0 + 4.0 * to_xxyyz_yy[k] * tbe_0 * tke_0;

            to_z_y_xxyy_z[k] = 4.0 * to_xxyyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 327-330 components of targeted buffer : GP

        auto to_z_y_xxyz_x = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 12);

        auto to_z_y_xxyz_y = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 13);

        auto to_z_y_xxyz_z = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 14);

#pragma omp simd aligned(to_xxy_0,          \
                             to_xxy_xy,     \
                             to_xxy_yy,     \
                             to_xxy_yz,     \
                             to_xxyzz_0,    \
                             to_xxyzz_xy,   \
                             to_xxyzz_yy,   \
                             to_xxyzz_yz,   \
                             to_z_y_xxyz_x, \
                             to_z_y_xxyz_y, \
                             to_z_y_xxyz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxyz_x[k] = -2.0 * to_xxy_xy[k] * tke_0 + 4.0 * to_xxyzz_xy[k] * tbe_0 * tke_0;

            to_z_y_xxyz_y[k] = to_xxy_0[k] - 2.0 * to_xxy_yy[k] * tke_0 - 2.0 * to_xxyzz_0[k] * tbe_0 + 4.0 * to_xxyzz_yy[k] * tbe_0 * tke_0;

            to_z_y_xxyz_z[k] = -2.0 * to_xxy_yz[k] * tke_0 + 4.0 * to_xxyzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 330-333 components of targeted buffer : GP

        auto to_z_y_xxzz_x = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 15);

        auto to_z_y_xxzz_y = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 16);

        auto to_z_y_xxzz_z = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 17);

#pragma omp simd aligned(to_xxz_0,          \
                             to_xxz_xy,     \
                             to_xxz_yy,     \
                             to_xxz_yz,     \
                             to_xxzzz_0,    \
                             to_xxzzz_xy,   \
                             to_xxzzz_yy,   \
                             to_xxzzz_yz,   \
                             to_z_y_xxzz_x, \
                             to_z_y_xxzz_y, \
                             to_z_y_xxzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxzz_x[k] = -4.0 * to_xxz_xy[k] * tke_0 + 4.0 * to_xxzzz_xy[k] * tbe_0 * tke_0;

            to_z_y_xxzz_y[k] = 2.0 * to_xxz_0[k] - 4.0 * to_xxz_yy[k] * tke_0 - 2.0 * to_xxzzz_0[k] * tbe_0 + 4.0 * to_xxzzz_yy[k] * tbe_0 * tke_0;

            to_z_y_xxzz_z[k] = -4.0 * to_xxz_yz[k] * tke_0 + 4.0 * to_xxzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 333-336 components of targeted buffer : GP

        auto to_z_y_xyyy_x = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 18);

        auto to_z_y_xyyy_y = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 19);

        auto to_z_y_xyyy_z = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 20);

#pragma omp simd aligned(to_xyyyz_0, to_xyyyz_xy, to_xyyyz_yy, to_xyyyz_yz, to_z_y_xyyy_x, to_z_y_xyyy_y, to_z_y_xyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyyy_x[k] = 4.0 * to_xyyyz_xy[k] * tbe_0 * tke_0;

            to_z_y_xyyy_y[k] = -2.0 * to_xyyyz_0[k] * tbe_0 + 4.0 * to_xyyyz_yy[k] * tbe_0 * tke_0;

            to_z_y_xyyy_z[k] = 4.0 * to_xyyyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 336-339 components of targeted buffer : GP

        auto to_z_y_xyyz_x = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 21);

        auto to_z_y_xyyz_y = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 22);

        auto to_z_y_xyyz_z = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 23);

#pragma omp simd aligned(to_xyy_0,          \
                             to_xyy_xy,     \
                             to_xyy_yy,     \
                             to_xyy_yz,     \
                             to_xyyzz_0,    \
                             to_xyyzz_xy,   \
                             to_xyyzz_yy,   \
                             to_xyyzz_yz,   \
                             to_z_y_xyyz_x, \
                             to_z_y_xyyz_y, \
                             to_z_y_xyyz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyyz_x[k] = -2.0 * to_xyy_xy[k] * tke_0 + 4.0 * to_xyyzz_xy[k] * tbe_0 * tke_0;

            to_z_y_xyyz_y[k] = to_xyy_0[k] - 2.0 * to_xyy_yy[k] * tke_0 - 2.0 * to_xyyzz_0[k] * tbe_0 + 4.0 * to_xyyzz_yy[k] * tbe_0 * tke_0;

            to_z_y_xyyz_z[k] = -2.0 * to_xyy_yz[k] * tke_0 + 4.0 * to_xyyzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 339-342 components of targeted buffer : GP

        auto to_z_y_xyzz_x = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 24);

        auto to_z_y_xyzz_y = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 25);

        auto to_z_y_xyzz_z = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 26);

#pragma omp simd aligned(to_xyz_0,          \
                             to_xyz_xy,     \
                             to_xyz_yy,     \
                             to_xyz_yz,     \
                             to_xyzzz_0,    \
                             to_xyzzz_xy,   \
                             to_xyzzz_yy,   \
                             to_xyzzz_yz,   \
                             to_z_y_xyzz_x, \
                             to_z_y_xyzz_y, \
                             to_z_y_xyzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyzz_x[k] = -4.0 * to_xyz_xy[k] * tke_0 + 4.0 * to_xyzzz_xy[k] * tbe_0 * tke_0;

            to_z_y_xyzz_y[k] = 2.0 * to_xyz_0[k] - 4.0 * to_xyz_yy[k] * tke_0 - 2.0 * to_xyzzz_0[k] * tbe_0 + 4.0 * to_xyzzz_yy[k] * tbe_0 * tke_0;

            to_z_y_xyzz_z[k] = -4.0 * to_xyz_yz[k] * tke_0 + 4.0 * to_xyzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 342-345 components of targeted buffer : GP

        auto to_z_y_xzzz_x = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 27);

        auto to_z_y_xzzz_y = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 28);

        auto to_z_y_xzzz_z = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 29);

#pragma omp simd aligned(to_xzz_0,          \
                             to_xzz_xy,     \
                             to_xzz_yy,     \
                             to_xzz_yz,     \
                             to_xzzzz_0,    \
                             to_xzzzz_xy,   \
                             to_xzzzz_yy,   \
                             to_xzzzz_yz,   \
                             to_z_y_xzzz_x, \
                             to_z_y_xzzz_y, \
                             to_z_y_xzzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xzzz_x[k] = -6.0 * to_xzz_xy[k] * tke_0 + 4.0 * to_xzzzz_xy[k] * tbe_0 * tke_0;

            to_z_y_xzzz_y[k] = 3.0 * to_xzz_0[k] - 6.0 * to_xzz_yy[k] * tke_0 - 2.0 * to_xzzzz_0[k] * tbe_0 + 4.0 * to_xzzzz_yy[k] * tbe_0 * tke_0;

            to_z_y_xzzz_z[k] = -6.0 * to_xzz_yz[k] * tke_0 + 4.0 * to_xzzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 345-348 components of targeted buffer : GP

        auto to_z_y_yyyy_x = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 30);

        auto to_z_y_yyyy_y = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 31);

        auto to_z_y_yyyy_z = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 32);

#pragma omp simd aligned(to_yyyyz_0, to_yyyyz_xy, to_yyyyz_yy, to_yyyyz_yz, to_z_y_yyyy_x, to_z_y_yyyy_y, to_z_y_yyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyyy_x[k] = 4.0 * to_yyyyz_xy[k] * tbe_0 * tke_0;

            to_z_y_yyyy_y[k] = -2.0 * to_yyyyz_0[k] * tbe_0 + 4.0 * to_yyyyz_yy[k] * tbe_0 * tke_0;

            to_z_y_yyyy_z[k] = 4.0 * to_yyyyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 348-351 components of targeted buffer : GP

        auto to_z_y_yyyz_x = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 33);

        auto to_z_y_yyyz_y = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 34);

        auto to_z_y_yyyz_z = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 35);

#pragma omp simd aligned(to_yyy_0,          \
                             to_yyy_xy,     \
                             to_yyy_yy,     \
                             to_yyy_yz,     \
                             to_yyyzz_0,    \
                             to_yyyzz_xy,   \
                             to_yyyzz_yy,   \
                             to_yyyzz_yz,   \
                             to_z_y_yyyz_x, \
                             to_z_y_yyyz_y, \
                             to_z_y_yyyz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyyz_x[k] = -2.0 * to_yyy_xy[k] * tke_0 + 4.0 * to_yyyzz_xy[k] * tbe_0 * tke_0;

            to_z_y_yyyz_y[k] = to_yyy_0[k] - 2.0 * to_yyy_yy[k] * tke_0 - 2.0 * to_yyyzz_0[k] * tbe_0 + 4.0 * to_yyyzz_yy[k] * tbe_0 * tke_0;

            to_z_y_yyyz_z[k] = -2.0 * to_yyy_yz[k] * tke_0 + 4.0 * to_yyyzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 351-354 components of targeted buffer : GP

        auto to_z_y_yyzz_x = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 36);

        auto to_z_y_yyzz_y = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 37);

        auto to_z_y_yyzz_z = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 38);

#pragma omp simd aligned(to_yyz_0,          \
                             to_yyz_xy,     \
                             to_yyz_yy,     \
                             to_yyz_yz,     \
                             to_yyzzz_0,    \
                             to_yyzzz_xy,   \
                             to_yyzzz_yy,   \
                             to_yyzzz_yz,   \
                             to_z_y_yyzz_x, \
                             to_z_y_yyzz_y, \
                             to_z_y_yyzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyzz_x[k] = -4.0 * to_yyz_xy[k] * tke_0 + 4.0 * to_yyzzz_xy[k] * tbe_0 * tke_0;

            to_z_y_yyzz_y[k] = 2.0 * to_yyz_0[k] - 4.0 * to_yyz_yy[k] * tke_0 - 2.0 * to_yyzzz_0[k] * tbe_0 + 4.0 * to_yyzzz_yy[k] * tbe_0 * tke_0;

            to_z_y_yyzz_z[k] = -4.0 * to_yyz_yz[k] * tke_0 + 4.0 * to_yyzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 354-357 components of targeted buffer : GP

        auto to_z_y_yzzz_x = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 39);

        auto to_z_y_yzzz_y = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 40);

        auto to_z_y_yzzz_z = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 41);

#pragma omp simd aligned(to_yzz_0,          \
                             to_yzz_xy,     \
                             to_yzz_yy,     \
                             to_yzz_yz,     \
                             to_yzzzz_0,    \
                             to_yzzzz_xy,   \
                             to_yzzzz_yy,   \
                             to_yzzzz_yz,   \
                             to_z_y_yzzz_x, \
                             to_z_y_yzzz_y, \
                             to_z_y_yzzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yzzz_x[k] = -6.0 * to_yzz_xy[k] * tke_0 + 4.0 * to_yzzzz_xy[k] * tbe_0 * tke_0;

            to_z_y_yzzz_y[k] = 3.0 * to_yzz_0[k] - 6.0 * to_yzz_yy[k] * tke_0 - 2.0 * to_yzzzz_0[k] * tbe_0 + 4.0 * to_yzzzz_yy[k] * tbe_0 * tke_0;

            to_z_y_yzzz_z[k] = -6.0 * to_yzz_yz[k] * tke_0 + 4.0 * to_yzzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 357-360 components of targeted buffer : GP

        auto to_z_y_zzzz_x = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 42);

        auto to_z_y_zzzz_y = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 43);

        auto to_z_y_zzzz_z = pbuffer.data(idx_op_geom_101_gp + 7 * op_comps * 45 + i * 45 + 44);

#pragma omp simd aligned(to_z_y_zzzz_x,     \
                             to_z_y_zzzz_y, \
                             to_z_y_zzzz_z, \
                             to_zzz_0,      \
                             to_zzz_xy,     \
                             to_zzz_yy,     \
                             to_zzz_yz,     \
                             to_zzzzz_0,    \
                             to_zzzzz_xy,   \
                             to_zzzzz_yy,   \
                             to_zzzzz_yz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_zzzz_x[k] = -8.0 * to_zzz_xy[k] * tke_0 + 4.0 * to_zzzzz_xy[k] * tbe_0 * tke_0;

            to_z_y_zzzz_y[k] = 4.0 * to_zzz_0[k] - 8.0 * to_zzz_yy[k] * tke_0 - 2.0 * to_zzzzz_0[k] * tbe_0 + 4.0 * to_zzzzz_yy[k] * tbe_0 * tke_0;

            to_z_y_zzzz_z[k] = -8.0 * to_zzz_yz[k] * tke_0 + 4.0 * to_zzzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 360-363 components of targeted buffer : GP

        auto to_z_z_xxxx_x = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 0);

        auto to_z_z_xxxx_y = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 1);

        auto to_z_z_xxxx_z = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 2);

#pragma omp simd aligned(to_xxxxz_0, to_xxxxz_xz, to_xxxxz_yz, to_xxxxz_zz, to_z_z_xxxx_x, to_z_z_xxxx_y, to_z_z_xxxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxxx_x[k] = 4.0 * to_xxxxz_xz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_y[k] = 4.0 * to_xxxxz_yz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_z[k] = -2.0 * to_xxxxz_0[k] * tbe_0 + 4.0 * to_xxxxz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 363-366 components of targeted buffer : GP

        auto to_z_z_xxxy_x = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 3);

        auto to_z_z_xxxy_y = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 4);

        auto to_z_z_xxxy_z = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 5);

#pragma omp simd aligned(to_xxxyz_0, to_xxxyz_xz, to_xxxyz_yz, to_xxxyz_zz, to_z_z_xxxy_x, to_z_z_xxxy_y, to_z_z_xxxy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxxy_x[k] = 4.0 * to_xxxyz_xz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_y[k] = 4.0 * to_xxxyz_yz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_z[k] = -2.0 * to_xxxyz_0[k] * tbe_0 + 4.0 * to_xxxyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 366-369 components of targeted buffer : GP

        auto to_z_z_xxxz_x = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 6);

        auto to_z_z_xxxz_y = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 7);

        auto to_z_z_xxxz_z = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 8);

#pragma omp simd aligned(to_xxx_0,          \
                             to_xxx_xz,     \
                             to_xxx_yz,     \
                             to_xxx_zz,     \
                             to_xxxzz_0,    \
                             to_xxxzz_xz,   \
                             to_xxxzz_yz,   \
                             to_xxxzz_zz,   \
                             to_z_z_xxxz_x, \
                             to_z_z_xxxz_y, \
                             to_z_z_xxxz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxxz_x[k] = -2.0 * to_xxx_xz[k] * tke_0 + 4.0 * to_xxxzz_xz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_y[k] = -2.0 * to_xxx_yz[k] * tke_0 + 4.0 * to_xxxzz_yz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_z[k] = to_xxx_0[k] - 2.0 * to_xxx_zz[k] * tke_0 - 2.0 * to_xxxzz_0[k] * tbe_0 + 4.0 * to_xxxzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 369-372 components of targeted buffer : GP

        auto to_z_z_xxyy_x = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 9);

        auto to_z_z_xxyy_y = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 10);

        auto to_z_z_xxyy_z = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 11);

#pragma omp simd aligned(to_xxyyz_0, to_xxyyz_xz, to_xxyyz_yz, to_xxyyz_zz, to_z_z_xxyy_x, to_z_z_xxyy_y, to_z_z_xxyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxyy_x[k] = 4.0 * to_xxyyz_xz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_y[k] = 4.0 * to_xxyyz_yz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_z[k] = -2.0 * to_xxyyz_0[k] * tbe_0 + 4.0 * to_xxyyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 372-375 components of targeted buffer : GP

        auto to_z_z_xxyz_x = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 12);

        auto to_z_z_xxyz_y = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 13);

        auto to_z_z_xxyz_z = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 14);

#pragma omp simd aligned(to_xxy_0,          \
                             to_xxy_xz,     \
                             to_xxy_yz,     \
                             to_xxy_zz,     \
                             to_xxyzz_0,    \
                             to_xxyzz_xz,   \
                             to_xxyzz_yz,   \
                             to_xxyzz_zz,   \
                             to_z_z_xxyz_x, \
                             to_z_z_xxyz_y, \
                             to_z_z_xxyz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxyz_x[k] = -2.0 * to_xxy_xz[k] * tke_0 + 4.0 * to_xxyzz_xz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_y[k] = -2.0 * to_xxy_yz[k] * tke_0 + 4.0 * to_xxyzz_yz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_z[k] = to_xxy_0[k] - 2.0 * to_xxy_zz[k] * tke_0 - 2.0 * to_xxyzz_0[k] * tbe_0 + 4.0 * to_xxyzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 375-378 components of targeted buffer : GP

        auto to_z_z_xxzz_x = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 15);

        auto to_z_z_xxzz_y = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 16);

        auto to_z_z_xxzz_z = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 17);

#pragma omp simd aligned(to_xxz_0,          \
                             to_xxz_xz,     \
                             to_xxz_yz,     \
                             to_xxz_zz,     \
                             to_xxzzz_0,    \
                             to_xxzzz_xz,   \
                             to_xxzzz_yz,   \
                             to_xxzzz_zz,   \
                             to_z_z_xxzz_x, \
                             to_z_z_xxzz_y, \
                             to_z_z_xxzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxzz_x[k] = -4.0 * to_xxz_xz[k] * tke_0 + 4.0 * to_xxzzz_xz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_y[k] = -4.0 * to_xxz_yz[k] * tke_0 + 4.0 * to_xxzzz_yz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_z[k] = 2.0 * to_xxz_0[k] - 4.0 * to_xxz_zz[k] * tke_0 - 2.0 * to_xxzzz_0[k] * tbe_0 + 4.0 * to_xxzzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 378-381 components of targeted buffer : GP

        auto to_z_z_xyyy_x = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 18);

        auto to_z_z_xyyy_y = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 19);

        auto to_z_z_xyyy_z = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 20);

#pragma omp simd aligned(to_xyyyz_0, to_xyyyz_xz, to_xyyyz_yz, to_xyyyz_zz, to_z_z_xyyy_x, to_z_z_xyyy_y, to_z_z_xyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyyy_x[k] = 4.0 * to_xyyyz_xz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_y[k] = 4.0 * to_xyyyz_yz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_z[k] = -2.0 * to_xyyyz_0[k] * tbe_0 + 4.0 * to_xyyyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 381-384 components of targeted buffer : GP

        auto to_z_z_xyyz_x = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 21);

        auto to_z_z_xyyz_y = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 22);

        auto to_z_z_xyyz_z = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 23);

#pragma omp simd aligned(to_xyy_0,          \
                             to_xyy_xz,     \
                             to_xyy_yz,     \
                             to_xyy_zz,     \
                             to_xyyzz_0,    \
                             to_xyyzz_xz,   \
                             to_xyyzz_yz,   \
                             to_xyyzz_zz,   \
                             to_z_z_xyyz_x, \
                             to_z_z_xyyz_y, \
                             to_z_z_xyyz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyyz_x[k] = -2.0 * to_xyy_xz[k] * tke_0 + 4.0 * to_xyyzz_xz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_y[k] = -2.0 * to_xyy_yz[k] * tke_0 + 4.0 * to_xyyzz_yz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_z[k] = to_xyy_0[k] - 2.0 * to_xyy_zz[k] * tke_0 - 2.0 * to_xyyzz_0[k] * tbe_0 + 4.0 * to_xyyzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 384-387 components of targeted buffer : GP

        auto to_z_z_xyzz_x = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 24);

        auto to_z_z_xyzz_y = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 25);

        auto to_z_z_xyzz_z = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 26);

#pragma omp simd aligned(to_xyz_0,          \
                             to_xyz_xz,     \
                             to_xyz_yz,     \
                             to_xyz_zz,     \
                             to_xyzzz_0,    \
                             to_xyzzz_xz,   \
                             to_xyzzz_yz,   \
                             to_xyzzz_zz,   \
                             to_z_z_xyzz_x, \
                             to_z_z_xyzz_y, \
                             to_z_z_xyzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyzz_x[k] = -4.0 * to_xyz_xz[k] * tke_0 + 4.0 * to_xyzzz_xz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_y[k] = -4.0 * to_xyz_yz[k] * tke_0 + 4.0 * to_xyzzz_yz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_z[k] = 2.0 * to_xyz_0[k] - 4.0 * to_xyz_zz[k] * tke_0 - 2.0 * to_xyzzz_0[k] * tbe_0 + 4.0 * to_xyzzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 387-390 components of targeted buffer : GP

        auto to_z_z_xzzz_x = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 27);

        auto to_z_z_xzzz_y = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 28);

        auto to_z_z_xzzz_z = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 29);

#pragma omp simd aligned(to_xzz_0,          \
                             to_xzz_xz,     \
                             to_xzz_yz,     \
                             to_xzz_zz,     \
                             to_xzzzz_0,    \
                             to_xzzzz_xz,   \
                             to_xzzzz_yz,   \
                             to_xzzzz_zz,   \
                             to_z_z_xzzz_x, \
                             to_z_z_xzzz_y, \
                             to_z_z_xzzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xzzz_x[k] = -6.0 * to_xzz_xz[k] * tke_0 + 4.0 * to_xzzzz_xz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_y[k] = -6.0 * to_xzz_yz[k] * tke_0 + 4.0 * to_xzzzz_yz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_z[k] = 3.0 * to_xzz_0[k] - 6.0 * to_xzz_zz[k] * tke_0 - 2.0 * to_xzzzz_0[k] * tbe_0 + 4.0 * to_xzzzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 390-393 components of targeted buffer : GP

        auto to_z_z_yyyy_x = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 30);

        auto to_z_z_yyyy_y = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 31);

        auto to_z_z_yyyy_z = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 32);

#pragma omp simd aligned(to_yyyyz_0, to_yyyyz_xz, to_yyyyz_yz, to_yyyyz_zz, to_z_z_yyyy_x, to_z_z_yyyy_y, to_z_z_yyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyyy_x[k] = 4.0 * to_yyyyz_xz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_y[k] = 4.0 * to_yyyyz_yz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_z[k] = -2.0 * to_yyyyz_0[k] * tbe_0 + 4.0 * to_yyyyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 393-396 components of targeted buffer : GP

        auto to_z_z_yyyz_x = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 33);

        auto to_z_z_yyyz_y = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 34);

        auto to_z_z_yyyz_z = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 35);

#pragma omp simd aligned(to_yyy_0,          \
                             to_yyy_xz,     \
                             to_yyy_yz,     \
                             to_yyy_zz,     \
                             to_yyyzz_0,    \
                             to_yyyzz_xz,   \
                             to_yyyzz_yz,   \
                             to_yyyzz_zz,   \
                             to_z_z_yyyz_x, \
                             to_z_z_yyyz_y, \
                             to_z_z_yyyz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyyz_x[k] = -2.0 * to_yyy_xz[k] * tke_0 + 4.0 * to_yyyzz_xz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_y[k] = -2.0 * to_yyy_yz[k] * tke_0 + 4.0 * to_yyyzz_yz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_z[k] = to_yyy_0[k] - 2.0 * to_yyy_zz[k] * tke_0 - 2.0 * to_yyyzz_0[k] * tbe_0 + 4.0 * to_yyyzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 396-399 components of targeted buffer : GP

        auto to_z_z_yyzz_x = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 36);

        auto to_z_z_yyzz_y = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 37);

        auto to_z_z_yyzz_z = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 38);

#pragma omp simd aligned(to_yyz_0,          \
                             to_yyz_xz,     \
                             to_yyz_yz,     \
                             to_yyz_zz,     \
                             to_yyzzz_0,    \
                             to_yyzzz_xz,   \
                             to_yyzzz_yz,   \
                             to_yyzzz_zz,   \
                             to_z_z_yyzz_x, \
                             to_z_z_yyzz_y, \
                             to_z_z_yyzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyzz_x[k] = -4.0 * to_yyz_xz[k] * tke_0 + 4.0 * to_yyzzz_xz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_y[k] = -4.0 * to_yyz_yz[k] * tke_0 + 4.0 * to_yyzzz_yz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_z[k] = 2.0 * to_yyz_0[k] - 4.0 * to_yyz_zz[k] * tke_0 - 2.0 * to_yyzzz_0[k] * tbe_0 + 4.0 * to_yyzzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 399-402 components of targeted buffer : GP

        auto to_z_z_yzzz_x = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 39);

        auto to_z_z_yzzz_y = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 40);

        auto to_z_z_yzzz_z = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 41);

#pragma omp simd aligned(to_yzz_0,          \
                             to_yzz_xz,     \
                             to_yzz_yz,     \
                             to_yzz_zz,     \
                             to_yzzzz_0,    \
                             to_yzzzz_xz,   \
                             to_yzzzz_yz,   \
                             to_yzzzz_zz,   \
                             to_z_z_yzzz_x, \
                             to_z_z_yzzz_y, \
                             to_z_z_yzzz_z, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yzzz_x[k] = -6.0 * to_yzz_xz[k] * tke_0 + 4.0 * to_yzzzz_xz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_y[k] = -6.0 * to_yzz_yz[k] * tke_0 + 4.0 * to_yzzzz_yz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_z[k] = 3.0 * to_yzz_0[k] - 6.0 * to_yzz_zz[k] * tke_0 - 2.0 * to_yzzzz_0[k] * tbe_0 + 4.0 * to_yzzzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 402-405 components of targeted buffer : GP

        auto to_z_z_zzzz_x = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 42);

        auto to_z_z_zzzz_y = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 43);

        auto to_z_z_zzzz_z = pbuffer.data(idx_op_geom_101_gp + 8 * op_comps * 45 + i * 45 + 44);

#pragma omp simd aligned(to_z_z_zzzz_x,     \
                             to_z_z_zzzz_y, \
                             to_z_z_zzzz_z, \
                             to_zzz_0,      \
                             to_zzz_xz,     \
                             to_zzz_yz,     \
                             to_zzz_zz,     \
                             to_zzzzz_0,    \
                             to_zzzzz_xz,   \
                             to_zzzzz_yz,   \
                             to_zzzzz_zz,   \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_zzzz_x[k] = -8.0 * to_zzz_xz[k] * tke_0 + 4.0 * to_zzzzz_xz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_y[k] = -8.0 * to_zzz_yz[k] * tke_0 + 4.0 * to_zzzzz_yz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_z[k] = 4.0 * to_zzz_0[k] - 8.0 * to_zzz_zz[k] * tke_0 - 2.0 * to_zzzzz_0[k] * tbe_0 + 4.0 * to_zzzzz_zz[k] * tbe_0 * tke_0;
        }
    }
}

}  // namespace t2cgeom
