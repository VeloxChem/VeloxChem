//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

namespace kinvecfunc { // kinvecfunc namespace

    // SIMD elementary functions for (D|T|D) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xx_s_0(double fx,
                                 double pa_x,
                                 double pa_xx,
                                 double pb_x,
                                 double pb_xx,
                                 double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx

                     + 0.5 * pa_xx * fx

                     + 2.0 * pa_x * fx * pb_x

                     + 0.5 * fx * pb_xx

                     + pa_xx * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xx_r_0(double fga,
                                 double fgb,
                                 double fx,
                                 double fz,
                                 double pa_x,
                                 double pa_xx,
                                 double pb_x,
                                 double pb_xx,
                                 double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fx * fz

                     - 0.5 * fx * fz * fgb

                     - 0.5 * fz * fga * fx

                     - pa_xx * fz * fgb

                     + 3.0 * pa_xx * fz * fx

                     + 12.0 * pa_x * fz * fx * pb_x

                     - fz * fga * pb_xx

                     + 3.0 * fz * fx * pb_xx

                     + 8.0 * pa_xx * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xy_s_0(double fx,
                                 double pa_x,
                                 double pa_xx,
                                 double pb_xy,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (pa_x * fx * pb_y

                     + 0.5 * fx * pb_xy

                     + pa_xx * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xy_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_x,
                                 double pa_xx,
                                 double pb_xy,
                                 double pb_y,
                                 double r_0_0)
    {
        return r_0_0 * (6.0 * pa_x * fz * fx * pb_y

                     - fz * fga * pb_xy

                     + 3.0 * fz * fx * pb_xy

                     + 8.0 * pa_xx * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xz_s_0(double fx,
                                 double pa_x,
                                 double pa_xx,
                                 double pb_xz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (pa_x * fx * pb_z

                     + 0.5 * fx * pb_xz

                     + pa_xx * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xz_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_x,
                                 double pa_xx,
                                 double pb_xz,
                                 double pb_z,
                                 double r_0_0)
    {
        return r_0_0 * (6.0 * pa_x * fz * fx * pb_z

                     - fz * fga * pb_xz

                     + 3.0 * fz * fx * pb_xz

                     + 8.0 * pa_xx * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yy_s_0(double fx,
                                 double pa_xx,
                                 double pb_yy,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_xx * fx

                     + 0.5 * fx * pb_yy

                     + pa_xx * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yy_r_0(double fga,
                                 double fgb,
                                 double fx,
                                 double fz,
                                 double pa_xx,
                                 double pb_yy,
                                 double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fz * fgb

                     - 0.5 * fz * fga * fx

                     - pa_xx * fz * fgb

                     + fz * fx * fx

                     + 3.0 * pa_xx * fz * fx

                     - fz * fga * pb_yy

                     + 3.0 * fz * fx * pb_yy

                     + 8.0 * pa_xx * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yz_s_0(double fx,
                                 double pa_xx,
                                 double pb_yz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_yz

                     + pa_xx * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yz_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_xx,
                                 double pb_yz,
                                 double r_0_0)
    {
        return r_0_0 * (- fz * fga * pb_yz

                     + 3.0 * fz * fx * pb_yz

                     + 8.0 * pa_xx * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_zz_s_0(double fx,
                                 double pa_xx,
                                 double pb_zz,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_xx * fx

                     + 0.5 * fx * pb_zz

                     + pa_xx * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_zz_r_0(double fga,
                                 double fgb,
                                 double fx,
                                 double fz,
                                 double pa_xx,
                                 double pb_zz,
                                 double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fz * fgb

                     - 0.5 * fz * fga * fx

                     - pa_xx * fz * fgb

                     + fz * fx * fx

                     + 3.0 * pa_xx * fz * fx

                     - fz * fga * pb_zz

                     + 3.0 * fz * fx * pb_zz

                     + 8.0 * pa_xx * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xx_s_0(double fx,
                                 double pa_xy,
                                 double pa_y,
                                 double pb_x,
                                 double pb_xx,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xy * fx

                     + fx * pa_y * pb_x

                     + pa_xy * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xx_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pa_xy,
                                 double pa_y,
                                 double pb_x,
                                 double pb_xx,
                                 double r_0_0)
    {
        return r_0_0 * (- pa_xy * fz * fgb

                     + 3.0 * pa_xy * fz * fx

                     + 6.0 * fx * fz * pa_y * pb_x

                     + 8.0 * pa_xy * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xy_s_0(double fx,
                                 double pa_x,
                                 double pa_xy,
                                 double pa_y,
                                 double pb_x,
                                 double pb_xy,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_x * fx * pb_x

                     + 0.5 * fx * pa_y * pb_y

                     + pa_xy * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xy_r_0(double fx,
                                 double fz,
                                 double pa_x,
                                 double pa_xy,
                                 double pa_y,
                                 double pb_x,
                                 double pb_xy,
                                 double pb_y,
                                 double r_0_0)
    {
        return r_0_0 * (fx * fx * fz

                     + 3.0 * pa_x * fz * fx * pb_x

                     + 3.0 * fx * fz * pa_y * pb_y

                     + 8.0 * pa_xy * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xz_s_0(double fx,
                                 double pa_xy,
                                 double pa_y,
                                 double pb_xz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pa_y * pb_z

                     + pa_xy * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xz_r_0(double fx,
                                 double fz,
                                 double pa_xy,
                                 double pa_y,
                                 double pb_xz,
                                 double pb_z,
                                 double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fz * pa_y * pb_z

                     + 8.0 * pa_xy * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yy_s_0(double fx,
                                 double pa_x,
                                 double pa_xy,
                                 double pb_y,
                                 double pb_yy,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xy * fx

                     + pa_x * fx * pb_y

                     + pa_xy * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yy_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pa_x,
                                 double pa_xy,
                                 double pb_y,
                                 double pb_yy,
                                 double r_0_0)
    {
        return r_0_0 * (- pa_xy * fz * fgb

                     + 3.0 * pa_xy * fz * fx

                     + 6.0 * pa_x * fz * fx * pb_y

                     + 8.0 * pa_xy * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yz_s_0(double fx,
                                 double pa_x,
                                 double pa_xy,
                                 double pb_yz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * pb_z

                     + pa_xy * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yz_r_0(double fx,
                                 double fz,
                                 double pa_x,
                                 double pa_xy,
                                 double pb_yz,
                                 double pb_z,
                                 double r_0_0)
    {
        return r_0_0 * (3.0 * pa_x * fz * fx * pb_z

                     + 8.0 * pa_xy * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_zz_s_0(double fx,
                                 double pa_xy,
                                 double pb_zz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xy * fx

                     + pa_xy * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_zz_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pa_xy,
                                 double pb_zz,
                                 double r_0_0)
    {
        return r_0_0 * (- pa_xy * fz * fgb

                     + 3.0 * pa_xy * fz * fx

                     + 8.0 * pa_xy * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xx_s_0(double fx,
                                 double pa_xz,
                                 double pa_z,
                                 double pb_x,
                                 double pb_xx,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xz * fx

                     + fx * pa_z * pb_x

                     + pa_xz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xx_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pa_xz,
                                 double pa_z,
                                 double pb_x,
                                 double pb_xx,
                                 double r_0_0)
    {
        return r_0_0 * (- pa_xz * fz * fgb

                     + 3.0 * pa_xz * fz * fx

                     + 6.0 * fx * fz * pa_z * pb_x

                     + 8.0 * pa_xz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xy_s_0(double fx,
                                 double pa_xz,
                                 double pa_z,
                                 double pb_xy,
                                 double pb_y,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pa_z * pb_y

                     + pa_xz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xy_r_0(double fx,
                                 double fz,
                                 double pa_xz,
                                 double pa_z,
                                 double pb_xy,
                                 double pb_y,
                                 double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fz * pa_z * pb_y

                     + 8.0 * pa_xz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xz_s_0(double fx,
                                 double pa_x,
                                 double pa_xz,
                                 double pa_z,
                                 double pb_x,
                                 double pb_xz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_x * fx * pb_x

                     + 0.5 * fx * pa_z * pb_z

                     + pa_xz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xz_r_0(double fx,
                                 double fz,
                                 double pa_x,
                                 double pa_xz,
                                 double pa_z,
                                 double pb_x,
                                 double pb_xz,
                                 double pb_z,
                                 double r_0_0)
    {
        return r_0_0 * (fx * fx * fz

                     + 3.0 * pa_x * fz * fx * pb_x

                     + 3.0 * fx * fz * pa_z * pb_z

                     + 8.0 * pa_xz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yy_s_0(double fx,
                                 double pa_xz,
                                 double pb_yy,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xz * fx

                     + pa_xz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yy_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pa_xz,
                                 double pb_yy,
                                 double r_0_0)
    {
        return r_0_0 * (- pa_xz * fz * fgb

                     + 3.0 * pa_xz * fz * fx

                     + 8.0 * pa_xz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yz_s_0(double fx,
                                 double pa_x,
                                 double pa_xz,
                                 double pb_y,
                                 double pb_yz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * pb_y

                     + pa_xz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yz_r_0(double fx,
                                 double fz,
                                 double pa_x,
                                 double pa_xz,
                                 double pb_y,
                                 double pb_yz,
                                 double r_0_0)
    {
        return r_0_0 * (3.0 * pa_x * fz * fx * pb_y

                     + 8.0 * pa_xz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_zz_s_0(double fx,
                                 double pa_x,
                                 double pa_xz,
                                 double pb_z,
                                 double pb_zz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xz * fx

                     + pa_x * fx * pb_z

                     + pa_xz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_zz_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pa_x,
                                 double pa_xz,
                                 double pb_z,
                                 double pb_zz,
                                 double r_0_0)
    {
        return r_0_0 * (- pa_xz * fz * fgb

                     + 3.0 * pa_xz * fz * fx

                     + 6.0 * pa_x * fz * fx * pb_z

                     + 8.0 * pa_xz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xx_s_0(double fx,
                                 double pa_yy,
                                 double pb_xx,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_yy * fx

                     + 0.5 * fx * pb_xx

                     + pa_yy * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xx_r_0(double fga,
                                 double fgb,
                                 double fx,
                                 double fz,
                                 double pa_yy,
                                 double pb_xx,
                                 double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fz * fgb

                     - 0.5 * fz * fga * fx

                     - pa_yy * fz * fgb

                     + fz * fx * fx

                     + 3.0 * pa_yy * fz * fx

                     - fz * fga * pb_xx

                     + 3.0 * fz * fx * pb_xx

                     + 8.0 * pa_yy * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xy_s_0(double fx,
                                 double pa_y,
                                 double pa_yy,
                                 double pb_x,
                                 double pb_xy,
                                 double s_0_0)
    {
        return s_0_0 * (pa_y * fx * pb_x

                     + 0.5 * fx * pb_xy

                     + pa_yy * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xy_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_y,
                                 double pa_yy,
                                 double pb_x,
                                 double pb_xy,
                                 double r_0_0)
    {
        return r_0_0 * (6.0 * pa_y * fz * fx * pb_x

                     - fz * fga * pb_xy

                     + 3.0 * fz * fx * pb_xy

                     + 8.0 * pa_yy * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xz_s_0(double fx,
                                 double pa_yy,
                                 double pb_xz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_xz

                     + pa_yy * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xz_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_yy,
                                 double pb_xz,
                                 double r_0_0)
    {
        return r_0_0 * (- fz * fga * pb_xz

                     + 3.0 * fz * fx * pb_xz

                     + 8.0 * pa_yy * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yy_s_0(double fx,
                                 double pa_y,
                                 double pa_yy,
                                 double pb_y,
                                 double pb_yy,
                                 double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx

                     + 0.5 * pa_yy * fx

                     + 2.0 * pa_y * fx * pb_y

                     + 0.5 * fx * pb_yy

                     + pa_yy * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yy_r_0(double fga,
                                 double fgb,
                                 double fx,
                                 double fz,
                                 double pa_y,
                                 double pa_yy,
                                 double pb_y,
                                 double pb_yy,
                                 double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fx * fz

                     - 0.5 * fx * fz * fgb

                     - 0.5 * fz * fga * fx

                     - pa_yy * fz * fgb

                     + 3.0 * pa_yy * fz * fx

                     + 12.0 * pa_y * fz * fx * pb_y

                     - fz * fga * pb_yy

                     + 3.0 * fz * fx * pb_yy

                     + 8.0 * pa_yy * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yz_s_0(double fx,
                                 double pa_y,
                                 double pa_yy,
                                 double pb_yz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (pa_y * fx * pb_z

                     + 0.5 * fx * pb_yz

                     + pa_yy * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yz_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_y,
                                 double pa_yy,
                                 double pb_yz,
                                 double pb_z,
                                 double r_0_0)
    {
        return r_0_0 * (6.0 * pa_y * fz * fx * pb_z

                     - fz * fga * pb_yz

                     + 3.0 * fz * fx * pb_yz

                     + 8.0 * pa_yy * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_zz_s_0(double fx,
                                 double pa_yy,
                                 double pb_zz,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_yy * fx

                     + 0.5 * fx * pb_zz

                     + pa_yy * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_zz_r_0(double fga,
                                 double fgb,
                                 double fx,
                                 double fz,
                                 double pa_yy,
                                 double pb_zz,
                                 double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fz * fgb

                     - 0.5 * fz * fga * fx

                     - pa_yy * fz * fgb

                     + fz * fx * fx

                     + 3.0 * pa_yy * fz * fx

                     - fz * fga * pb_zz

                     + 3.0 * fz * fx * pb_zz

                     + 8.0 * pa_yy * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xx_s_0(double fx,
                                 double pa_yz,
                                 double pb_xx,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_yz * fx

                     + pa_yz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xx_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pa_yz,
                                 double pb_xx,
                                 double r_0_0)
    {
        return r_0_0 * (- pa_yz * fz * fgb

                     + 3.0 * pa_yz * fz * fx

                     + 8.0 * pa_yz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xy_s_0(double fx,
                                 double pa_yz,
                                 double pa_z,
                                 double pb_x,
                                 double pb_xy,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pa_z * pb_x

                     + pa_yz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xy_r_0(double fx,
                                 double fz,
                                 double pa_yz,
                                 double pa_z,
                                 double pb_x,
                                 double pb_xy,
                                 double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fz * pa_z * pb_x

                     + 8.0 * pa_yz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xz_s_0(double fx,
                                 double pa_y,
                                 double pa_yz,
                                 double pb_x,
                                 double pb_xz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx * pb_x

                     + pa_yz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xz_r_0(double fx,
                                 double fz,
                                 double pa_y,
                                 double pa_yz,
                                 double pb_x,
                                 double pb_xz,
                                 double r_0_0)
    {
        return r_0_0 * (3.0 * pa_y * fz * fx * pb_x

                     + 8.0 * pa_yz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yy_s_0(double fx,
                                 double pa_yz,
                                 double pa_z,
                                 double pb_y,
                                 double pb_yy,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_yz * fx

                     + fx * pa_z * pb_y

                     + pa_yz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yy_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pa_yz,
                                 double pa_z,
                                 double pb_y,
                                 double pb_yy,
                                 double r_0_0)
    {
        return r_0_0 * (- pa_yz * fz * fgb

                     + 3.0 * pa_yz * fz * fx

                     + 6.0 * fx * fz * pa_z * pb_y

                     + 8.0 * pa_yz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yz_s_0(double fx,
                                 double pa_y,
                                 double pa_yz,
                                 double pa_z,
                                 double pb_y,
                                 double pb_yz,
                                 double pb_z,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_y * fx * pb_y

                     + 0.5 * fx * pa_z * pb_z

                     + pa_yz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yz_r_0(double fx,
                                 double fz,
                                 double pa_y,
                                 double pa_yz,
                                 double pa_z,
                                 double pb_y,
                                 double pb_yz,
                                 double pb_z,
                                 double r_0_0)
    {
        return r_0_0 * (fx * fx * fz

                     + 3.0 * pa_y * fz * fx * pb_y

                     + 3.0 * fx * fz * pa_z * pb_z

                     + 8.0 * pa_yz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_zz_s_0(double fx,
                                 double pa_y,
                                 double pa_yz,
                                 double pb_z,
                                 double pb_zz,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * pa_yz * fx

                     + pa_y * fx * pb_z

                     + pa_yz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_zz_r_0(double fgb,
                                 double fx,
                                 double fz,
                                 double pa_y,
                                 double pa_yz,
                                 double pb_z,
                                 double pb_zz,
                                 double r_0_0)
    {
        return r_0_0 * (- pa_yz * fz * fgb

                     + 3.0 * pa_yz * fz * fx

                     + 6.0 * pa_y * fz * fx * pb_z

                     + 8.0 * pa_yz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xx_s_0(double fx,
                                 double pa_zz,
                                 double pb_xx,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_zz * fx

                     + 0.5 * fx * pb_xx

                     + pa_zz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xx_r_0(double fga,
                                 double fgb,
                                 double fx,
                                 double fz,
                                 double pa_zz,
                                 double pb_xx,
                                 double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fz * fgb

                     - 0.5 * fz * fga * fx

                     - pa_zz * fz * fgb

                     + fz * fx * fx

                     + 3.0 * pa_zz * fz * fx

                     - fz * fga * pb_xx

                     + 3.0 * fz * fx * pb_xx

                     + 8.0 * pa_zz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xy_s_0(double fx,
                                 double pa_zz,
                                 double pb_xy,
                                 double s_0_0)
    {
        return s_0_0 * (0.5 * fx * pb_xy

                     + pa_zz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xy_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_zz,
                                 double pb_xy,
                                 double r_0_0)
    {
        return r_0_0 * (- fz * fga * pb_xy

                     + 3.0 * fz * fx * pb_xy

                     + 8.0 * pa_zz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xz_s_0(double fx,
                                 double pa_z,
                                 double pa_zz,
                                 double pb_x,
                                 double pb_xz,
                                 double s_0_0)
    {
        return s_0_0 * (pa_z * fx * pb_x

                     + 0.5 * fx * pb_xz

                     + pa_zz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xz_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_z,
                                 double pa_zz,
                                 double pb_x,
                                 double pb_xz,
                                 double r_0_0)
    {
        return r_0_0 * (6.0 * pa_z * fz * fx * pb_x

                     - fz * fga * pb_xz

                     + 3.0 * fz * fx * pb_xz

                     + 8.0 * pa_zz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yy_s_0(double fx,
                                 double pa_zz,
                                 double pb_yy,
                                 double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx

                     + 0.5 * pa_zz * fx

                     + 0.5 * fx * pb_yy

                     + pa_zz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yy_r_0(double fga,
                                 double fgb,
                                 double fx,
                                 double fz,
                                 double pa_zz,
                                 double pb_yy,
                                 double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fz * fgb

                     - 0.5 * fz * fga * fx

                     - pa_zz * fz * fgb

                     + fz * fx * fx

                     + 3.0 * pa_zz * fz * fx

                     - fz * fga * pb_yy

                     + 3.0 * fz * fx * pb_yy

                     + 8.0 * pa_zz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yz_s_0(double fx,
                                 double pa_z,
                                 double pa_zz,
                                 double pb_y,
                                 double pb_yz,
                                 double s_0_0)
    {
        return s_0_0 * (pa_z * fx * pb_y

                     + 0.5 * fx * pb_yz

                     + pa_zz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yz_r_0(double fga,
                                 double fx,
                                 double fz,
                                 double pa_z,
                                 double pa_zz,
                                 double pb_y,
                                 double pb_yz,
                                 double r_0_0)
    {
        return r_0_0 * (6.0 * pa_z * fz * fx * pb_y

                     - fz * fga * pb_yz

                     + 3.0 * fz * fx * pb_yz

                     + 8.0 * pa_zz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_zz_s_0(double fx,
                                 double pa_z,
                                 double pa_zz,
                                 double pb_z,
                                 double pb_zz,
                                 double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx

                     + 0.5 * pa_zz * fx

                     + 2.0 * pa_z * fx * pb_z

                     + 0.5 * fx * pb_zz

                     + pa_zz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_zz_r_0(double fga,
                                 double fgb,
                                 double fx,
                                 double fz,
                                 double pa_z,
                                 double pa_zz,
                                 double pb_z,
                                 double pb_zz,
                                 double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fx * fz

                     - 0.5 * fx * fz * fgb

                     - 0.5 * fz * fga * fx

                     - pa_zz * fz * fgb

                     + 3.0 * pa_zz * fz * fx

                     + 12.0 * pa_z * fz * fx * pb_z

                     - fz * fga * pb_zz

                     + 3.0 * fz * fx * pb_zz

                     + 8.0 * pa_zz * fz * pb_zz);

    }

    // SIMD elementary functions for (D|T|F) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxx_s_0(double fx,
                                  double pa_x,
                                  double pa_xx,
                                  double pb_x,
                                  double pb_xx,
                                  double pb_xxx,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * fx * fx

                     + 2.25 * fx * fx * pb_x

                     + 1.5 * pa_xx * pb_x * fx

                     + 3.0 * pa_x * fx * pb_xx

                     + 0.5 * fx * pb_xxx

                     + pa_xx * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxx_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xx,
                                  double pb_x,
                                  double pb_xx,
                                  double pb_xxx,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_x * fx * fz * fgb

                     + 9.0 * pa_x * fz * fx * fx

                     + 13.5 * fx * fx * fz * pb_x

                     - 1.5 * fx * pb_x * fz * fgb

                     - 1.5 * fz * fga * pb_x * fx

                     - 3.0 * pa_xx * pb_x * fz * fgb

                     + 12.0 * pa_xx * fz * pb_x * fx

                     + 24.0 * pa_x * fz * fx * pb_xx

                     - fz * fga * pb_xxx

                     + 4.0 * fz * fx * pb_xxx

                     + 10.0 * pa_xx * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxy_s_0(double fx,
                                  double pa_x,
                                  double pa_xx,
                                  double pb_xxy,
                                  double pb_xy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_y

                     + 0.5 * pa_xx * fx * pb_y

                     + 2.0 * pa_x * fx * pb_xy

                     + 0.5 * fx * pb_xxy

                     + pa_xx * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xx,
                                  double pb_xxy,
                                  double pb_xy,
                                  double pb_y,
                                  double r_0_0)
    {
        return r_0_0 * (4.5 * fx * fx * fz * pb_y

                     - 0.5 * fx * fz * fgb * pb_y

                     - 0.5 * fz * fga * fx * pb_y

                     - pa_xx * fz * fgb * pb_y

                     + 4.0 * pa_xx * fz * fx * pb_y

                     + 16.0 * pa_x * fz * fx * pb_xy

                     - fz * fga * pb_xxy

                     + 4.0 * fz * fx * pb_xxy

                     + 10.0 * pa_xx * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxz_s_0(double fx,
                                  double pa_x,
                                  double pa_xx,
                                  double pb_xxz,
                                  double pb_xz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_z

                     + 0.5 * pa_xx * fx * pb_z

                     + 2.0 * pa_x * fx * pb_xz

                     + 0.5 * fx * pb_xxz

                     + pa_xx * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xx,
                                  double pb_xxz,
                                  double pb_xz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (4.5 * fx * fx * fz * pb_z

                     - 0.5 * fx * fz * fgb * pb_z

                     - 0.5 * fz * fga * fx * pb_z

                     - pa_xx * fz * fgb * pb_z

                     + 4.0 * pa_xx * fz * fx * pb_z

                     + 16.0 * pa_x * fz * fx * pb_xz

                     - fz * fga * pb_xxz

                     + 4.0 * fz * fx * pb_xxz

                     + 10.0 * pa_xx * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xyy_s_0(double fx,
                                  double pa_x,
                                  double pa_xx,
                                  double pb_x,
                                  double pb_xyy,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * fx

                     + 0.25 * fx * fx * pb_x

                     + 0.5 * pa_xx * pb_x * fx

                     + pa_x * fx * pb_yy

                     + 0.5 * fx * pb_xyy

                     + pa_xx * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xyy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xx,
                                  double pb_x,
                                  double pb_xyy,
                                  double pb_yy,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_x * fx * fz * fgb

                     + 3.0 * pa_x * fz * fx * fx

                     - 0.5 * fx * pb_x * fz * fgb

                     - 0.5 * fz * fga * pb_x * fx

                     - pa_xx * pb_x * fz * fgb

                     + 1.5 * fz * fx * fx * pb_x

                     + 4.0 * pa_xx * fz * pb_x * fx

                     + 8.0 * pa_x * fz * fx * pb_yy

                     - fz * fga * pb_xyy

                     + 4.0 * fz * fx * pb_xyy

                     + 10.0 * pa_xx * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xyz_s_0(double fx,
                                  double pa_x,
                                  double pa_xx,
                                  double pb_xyz,
                                  double pb_yz,
                                  double s_0_0)
    {
        return s_0_0 * (pa_x * fx * pb_yz

                     + 0.5 * fx * pb_xyz

                     + pa_xx * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xyz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xx,
                                  double pb_xyz,
                                  double pb_yz,
                                  double r_0_0)
    {
        return r_0_0 * (8.0 * pa_x * fz * fx * pb_yz

                     - fz * fga * pb_xyz

                     + 4.0 * fz * fx * pb_xyz

                     + 10.0 * pa_xx * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xzz_s_0(double fx,
                                  double pa_x,
                                  double pa_xx,
                                  double pb_x,
                                  double pb_xzz,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * fx

                     + 0.25 * fx * fx * pb_x

                     + 0.5 * pa_xx * pb_x * fx

                     + pa_x * fx * pb_zz

                     + 0.5 * fx * pb_xzz

                     + pa_xx * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xzz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xx,
                                  double pb_x,
                                  double pb_xzz,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_x * fx * fz * fgb

                     + 3.0 * pa_x * fz * fx * fx

                     - 0.5 * fx * pb_x * fz * fgb

                     - 0.5 * fz * fga * pb_x * fx

                     - pa_xx * pb_x * fz * fgb

                     + 1.5 * fz * fx * fx * pb_x

                     + 4.0 * pa_xx * fz * pb_x * fx

                     + 8.0 * pa_x * fz * fx * pb_zz

                     - fz * fga * pb_xzz

                     + 4.0 * fz * fx * pb_xzz

                     + 10.0 * pa_xx * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yyy_s_0(double fx,
                                  double pa_xx,
                                  double pb_y,
                                  double pb_yyy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_y

                     + 1.5 * pa_xx * pb_y * fx

                     + 0.5 * fx * pb_yyy

                     + pa_xx * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yyy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xx,
                                  double pb_y,
                                  double pb_yyy,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pb_y * fz * fgb

                     - 1.5 * fz * fga * pb_y * fx

                     - 3.0 * pa_xx * pb_y * fz * fgb

                     + 4.5 * fz * fx * fx * pb_y

                     + 12.0 * pa_xx * fz * pb_y * fx

                     - fz * fga * pb_yyy

                     + 4.0 * fz * fx * pb_yyy

                     + 10.0 * pa_xx * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yyz_s_0(double fx,
                                  double pa_xx,
                                  double pb_yyz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_z

                     + 0.5 * pa_xx * fx * pb_z

                     + 0.5 * fx * pb_yyz

                     + pa_xx * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yyz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xx,
                                  double pb_yyz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fz * fgb * pb_z

                     - 0.5 * fz * fga * fx * pb_z

                     - pa_xx * fz * fgb * pb_z

                     + 1.5 * fz * fx * fx * pb_z

                     + 4.0 * pa_xx * fz * fx * pb_z

                     - fz * fga * pb_yyz

                     + 4.0 * fz * fx * pb_yyz

                     + 10.0 * pa_xx * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yzz_s_0(double fx,
                                  double pa_xx,
                                  double pb_y,
                                  double pb_yzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_y

                     + 0.5 * pa_xx * pb_y * fx

                     + 0.5 * fx * pb_yzz

                     + pa_xx * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yzz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xx,
                                  double pb_y,
                                  double pb_yzz,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pb_y * fz * fgb

                     - 0.5 * fz * fga * pb_y * fx

                     - pa_xx * pb_y * fz * fgb

                     + 1.5 * fz * fx * fx * pb_y

                     + 4.0 * pa_xx * fz * pb_y * fx

                     - fz * fga * pb_yzz

                     + 4.0 * fz * fx * pb_yzz

                     + 10.0 * pa_xx * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_zzz_s_0(double fx,
                                  double pa_xx,
                                  double pb_z,
                                  double pb_zzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_z

                     + 1.5 * pa_xx * pb_z * fx

                     + 0.5 * fx * pb_zzz

                     + pa_xx * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_zzz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xx,
                                  double pb_z,
                                  double pb_zzz,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pb_z * fz * fgb

                     - 1.5 * fz * fga * pb_z * fx

                     - 3.0 * pa_xx * pb_z * fz * fgb

                     + 4.5 * fz * fx * fx * pb_z

                     + 12.0 * pa_xx * fz * pb_z * fx

                     - fz * fga * pb_zzz

                     + 4.0 * fz * fx * pb_zzz

                     + 10.0 * pa_xx * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxx_s_0(double fx,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_x,
                                  double pb_xx,
                                  double pb_xxx,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_y

                     + 1.5 * pa_xy * pb_x * fx

                     + 1.5 * fx * pa_y * pb_xx

                     + pa_xy * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxx_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_x,
                                  double pb_xx,
                                  double pb_xxx,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_y * fz * fgb

                     + 4.5 * fx * fx * fz * pa_y

                     - 3.0 * pa_xy * pb_x * fz * fgb

                     + 12.0 * pa_xy * fz * pb_x * fx

                     + 12.0 * fx * fz * pa_y * pb_xx

                     + 10.0 * pa_xy * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxy_s_0(double fx,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_x,
                                  double pb_xx,
                                  double pb_xxy,
                                  double pb_xy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx

                     + 0.5 * fx * fx * pb_x

                     + 0.5 * pa_xy * fx * pb_y

                     + 0.5 * pa_x * fx * pb_xx

                     + fx * pa_y * pb_xy

                     + pa_xy * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxy_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_x,
                                  double pb_xx,
                                  double pb_xxy,
                                  double pb_xy,
                                  double pb_y,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fz * fgb

                     + 1.5 * pa_x * fz * fx * fx

                     + 3.0 * fx * fx * fz * pb_x

                     - pa_xy * fz * fgb * pb_y

                     + 4.0 * pa_xy * fz * fx * pb_y

                     + 4.0 * pa_x * fz * fx * pb_xx

                     + 8.0 * fx * fz * pa_y * pb_xy

                     + 10.0 * pa_xy * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxz_s_0(double fx,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_xxz,
                                  double pb_xz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xy * fx * pb_z

                     + fx * pa_y * pb_xz

                     + pa_xy * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_xxz,
                                  double pb_xz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_xy * fz * fgb * pb_z

                     + 4.0 * pa_xy * fz * fx * pb_z

                     + 8.0 * fx * fz * pa_y * pb_xz

                     + 10.0 * pa_xy * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xyy_s_0(double fx,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_xyy,
                                  double pb_y,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_y

                     + 0.5 * fx * fx * pb_y

                     + 0.5 * pa_xy * pb_x * fx

                     + pa_x * fx * pb_xy

                     + 0.5 * fx * pa_y * pb_yy

                     + pa_xy * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xyy_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_xyy,
                                  double pb_y,
                                  double pb_yy,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_y * fz * fgb

                     + 1.5 * fx * fx * fz * pa_y

                     + 3.0 * fx * fx * fz * pb_y

                     - pa_xy * pb_x * fz * fgb

                     + 4.0 * pa_xy * fz * pb_x * fx

                     + 8.0 * pa_x * fz * fx * pb_xy

                     + 4.0 * fx * fz * pa_y * pb_yy

                     + 10.0 * pa_xy * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xyz_s_0(double fx,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_xyz,
                                  double pb_xz,
                                  double pb_yz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_z

                     + 0.5 * pa_x * fx * pb_xz

                     + 0.5 * fx * pa_y * pb_yz

                     + pa_xy * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xyz_r_0(double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_xyz,
                                  double pb_xz,
                                  double pb_yz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (1.5 * fx * fx * fz * pb_z

                     + 4.0 * pa_x * fz * fx * pb_xz

                     + 4.0 * fx * fz * pa_y * pb_yz

                     + 10.0 * pa_xy * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xzz_s_0(double fx,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_x,
                                  double pb_xzz,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_y

                     + 0.5 * pa_xy * pb_x * fx

                     + 0.5 * fx * pa_y * pb_zz

                     + pa_xy * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xzz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_x,
                                  double pb_xzz,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_y * fz * fgb

                     + 1.5 * fx * fx * fz * pa_y

                     - pa_xy * pb_x * fz * fgb

                     + 4.0 * pa_xy * fz * pb_x * fx

                     + 4.0 * fx * fz * pa_y * pb_zz

                     + 10.0 * pa_xy * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yyy_s_0(double fx,
                                  double pa_x,
                                  double pa_xy,
                                  double pb_y,
                                  double pb_yy,
                                  double pb_yyy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx

                     + 1.5 * pa_xy * pb_y * fx

                     + 1.5 * pa_x * fx * pb_yy

                     + pa_xy * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yyy_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xy,
                                  double pb_y,
                                  double pb_yy,
                                  double pb_yyy,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fz * fgb

                     + 4.5 * pa_x * fz * fx * fx

                     - 3.0 * pa_xy * pb_y * fz * fgb

                     + 12.0 * pa_xy * fz * pb_y * fx

                     + 12.0 * pa_x * fz * fx * pb_yy

                     + 10.0 * pa_xy * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yyz_s_0(double fx,
                                  double pa_x,
                                  double pa_xy,
                                  double pb_yyz,
                                  double pb_yz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xy * fx * pb_z

                     + pa_x * fx * pb_yz

                     + pa_xy * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yyz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xy,
                                  double pb_yyz,
                                  double pb_yz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_xy * fz * fgb * pb_z

                     + 4.0 * pa_xy * fz * fx * pb_z

                     + 8.0 * pa_x * fz * fx * pb_yz

                     + 10.0 * pa_xy * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yzz_s_0(double fx,
                                  double pa_x,
                                  double pa_xy,
                                  double pb_y,
                                  double pb_yzz,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx

                     + 0.5 * pa_xy * pb_y * fx

                     + 0.5 * pa_x * fx * pb_zz

                     + pa_xy * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yzz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xy,
                                  double pb_y,
                                  double pb_yzz,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fz * fgb

                     + 1.5 * pa_x * fz * fx * fx

                     - pa_xy * pb_y * fz * fgb

                     + 4.0 * pa_xy * fz * pb_y * fx

                     + 4.0 * pa_x * fz * fx * pb_zz

                     + 10.0 * pa_xy * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_zzz_s_0(double fx,
                                  double pa_xy,
                                  double pb_z,
                                  double pb_zzz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xy * pb_z * fx

                     + pa_xy * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_zzz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xy,
                                  double pb_z,
                                  double pb_zzz,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * pb_z * fz * fgb

                     + 12.0 * pa_xy * fz * pb_z * fx

                     + 10.0 * pa_xy * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxx_s_0(double fx,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xx,
                                  double pb_xxx,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z

                     + 1.5 * pa_xz * pb_x * fx

                     + 1.5 * fx * pa_z * pb_xx

                     + pa_xz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxx_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xx,
                                  double pb_xxx,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_z * fz * fgb

                     + 4.5 * fx * fx * fz * pa_z

                     - 3.0 * pa_xz * pb_x * fz * fgb

                     + 12.0 * pa_xz * fz * pb_x * fx

                     + 12.0 * fx * fz * pa_z * pb_xx

                     + 10.0 * pa_xz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxy_s_0(double fx,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_xxy,
                                  double pb_xy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xz * fx * pb_y

                     + fx * pa_z * pb_xy

                     + pa_xz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxy_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_xxy,
                                  double pb_xy,
                                  double pb_y,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_xz * fz * fgb * pb_y

                     + 4.0 * pa_xz * fz * fx * pb_y

                     + 8.0 * fx * fz * pa_z * pb_xy

                     + 10.0 * pa_xz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxz_s_0(double fx,
                                  double pa_x,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xx,
                                  double pb_xxz,
                                  double pb_xz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx

                     + 0.5 * fx * fx * pb_x

                     + 0.5 * pa_xz * fx * pb_z

                     + 0.5 * pa_x * fx * pb_xx

                     + fx * pa_z * pb_xz

                     + pa_xz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xx,
                                  double pb_xxz,
                                  double pb_xz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fz * fgb

                     + 1.5 * pa_x * fz * fx * fx

                     + 3.0 * fx * fx * fz * pb_x

                     - pa_xz * fz * fgb * pb_z

                     + 4.0 * pa_xz * fz * fx * pb_z

                     + 4.0 * pa_x * fz * fx * pb_xx

                     + 8.0 * fx * fz * pa_z * pb_xz

                     + 10.0 * pa_xz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xyy_s_0(double fx,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xyy,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_z

                     + 0.5 * pa_xz * pb_x * fx

                     + 0.5 * fx * pa_z * pb_yy

                     + pa_xz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xyy_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xyy,
                                  double pb_yy,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_z * fz * fgb

                     + 1.5 * fx * fx * fz * pa_z

                     - pa_xz * pb_x * fz * fgb

                     + 4.0 * pa_xz * fz * pb_x * fx

                     + 4.0 * fx * fz * pa_z * pb_yy

                     + 10.0 * pa_xz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xyz_s_0(double fx,
                                  double pa_x,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_xy,
                                  double pb_xyz,
                                  double pb_y,
                                  double pb_yz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_y

                     + 0.5 * pa_x * fx * pb_xy

                     + 0.5 * fx * pa_z * pb_yz

                     + pa_xz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xyz_r_0(double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_xy,
                                  double pb_xyz,
                                  double pb_y,
                                  double pb_yz,
                                  double r_0_0)
    {
        return r_0_0 * (1.5 * fx * fx * fz * pb_y

                     + 4.0 * pa_x * fz * fx * pb_xy

                     + 4.0 * fx * fz * pa_z * pb_yz

                     + 10.0 * pa_xz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xzz_s_0(double fx,
                                  double pa_x,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xz,
                                  double pb_xzz,
                                  double pb_z,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_z

                     + 0.5 * fx * fx * pb_z

                     + 0.5 * pa_xz * pb_x * fx

                     + pa_x * fx * pb_xz

                     + 0.5 * fx * pa_z * pb_zz

                     + pa_xz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xzz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xz,
                                  double pb_xzz,
                                  double pb_z,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_z * fz * fgb

                     + 1.5 * fx * fx * fz * pa_z

                     + 3.0 * fx * fx * fz * pb_z

                     - pa_xz * pb_x * fz * fgb

                     + 4.0 * pa_xz * fz * pb_x * fx

                     + 8.0 * pa_x * fz * fx * pb_xz

                     + 4.0 * fx * fz * pa_z * pb_zz

                     + 10.0 * pa_xz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yyy_s_0(double fx,
                                  double pa_xz,
                                  double pb_y,
                                  double pb_yyy,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_xz * pb_y * fx

                     + pa_xz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yyy_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xz,
                                  double pb_y,
                                  double pb_yyy,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * pb_y * fz * fgb

                     + 12.0 * pa_xz * fz * pb_y * fx

                     + 10.0 * pa_xz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yyz_s_0(double fx,
                                  double pa_x,
                                  double pa_xz,
                                  double pb_yy,
                                  double pb_yyz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx

                     + 0.5 * pa_xz * fx * pb_z

                     + 0.5 * pa_x * fx * pb_yy

                     + pa_xz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yyz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xz,
                                  double pb_yy,
                                  double pb_yyz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fz * fgb

                     + 1.5 * pa_x * fz * fx * fx

                     - pa_xz * fz * fgb * pb_z

                     + 4.0 * pa_xz * fz * fx * pb_z

                     + 4.0 * pa_x * fz * fx * pb_yy

                     + 10.0 * pa_xz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yzz_s_0(double fx,
                                  double pa_x,
                                  double pa_xz,
                                  double pb_y,
                                  double pb_yz,
                                  double pb_yzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xz * pb_y * fx

                     + pa_x * fx * pb_yz

                     + pa_xz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yzz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xz,
                                  double pb_y,
                                  double pb_yz,
                                  double pb_yzz,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_xz * pb_y * fz * fgb

                     + 4.0 * pa_xz * fz * pb_y * fx

                     + 8.0 * pa_x * fz * fx * pb_yz

                     + 10.0 * pa_xz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_zzz_s_0(double fx,
                                  double pa_x,
                                  double pa_xz,
                                  double pb_z,
                                  double pb_zz,
                                  double pb_zzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx

                     + 1.5 * pa_xz * pb_z * fx

                     + 1.5 * pa_x * fx * pb_zz

                     + pa_xz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_zzz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xz,
                                  double pb_z,
                                  double pb_zz,
                                  double pb_zzz,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fz * fgb

                     + 4.5 * pa_x * fz * fx * fx

                     - 3.0 * pa_xz * pb_z * fz * fgb

                     + 12.0 * pa_xz * fz * pb_z * fx

                     + 12.0 * pa_x * fz * fx * pb_zz

                     + 10.0 * pa_xz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxx_s_0(double fx,
                                  double pa_yy,
                                  double pb_x,
                                  double pb_xxx,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_x

                     + 1.5 * pa_yy * pb_x * fx

                     + 0.5 * fx * pb_xxx

                     + pa_yy * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxx_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_yy,
                                  double pb_x,
                                  double pb_xxx,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pb_x * fz * fgb

                     - 1.5 * fz * fga * pb_x * fx

                     - 3.0 * pa_yy * pb_x * fz * fgb

                     + 4.5 * fz * fx * fx * pb_x

                     + 12.0 * pa_yy * fz * pb_x * fx

                     - fz * fga * pb_xxx

                     + 4.0 * fz * fx * pb_xxx

                     + 10.0 * pa_yy * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxy_s_0(double fx,
                                  double pa_y,
                                  double pa_yy,
                                  double pb_xx,
                                  double pb_xxy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx * fx

                     + 0.25 * fx * fx * pb_y

                     + 0.5 * pa_yy * fx * pb_y

                     + pa_y * fx * pb_xx

                     + 0.5 * fx * pb_xxy

                     + pa_yy * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yy,
                                  double pb_xx,
                                  double pb_xxy,
                                  double pb_y,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_y * fx * fz * fgb

                     + 3.0 * pa_y * fz * fx * fx

                     - 0.5 * fx * fz * fgb * pb_y

                     - 0.5 * fz * fga * fx * pb_y

                     - pa_yy * fz * fgb * pb_y

                     + 1.5 * fz * fx * fx * pb_y

                     + 4.0 * pa_yy * fz * fx * pb_y

                     + 8.0 * pa_y * fz * fx * pb_xx

                     - fz * fga * pb_xxy

                     + 4.0 * fz * fx * pb_xxy

                     + 10.0 * pa_yy * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxz_s_0(double fx,
                                  double pa_yy,
                                  double pb_xxz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_z

                     + 0.5 * pa_yy * fx * pb_z

                     + 0.5 * fx * pb_xxz

                     + pa_yy * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_yy,
                                  double pb_xxz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fz * fgb * pb_z

                     - 0.5 * fz * fga * fx * pb_z

                     - pa_yy * fz * fgb * pb_z

                     + 1.5 * fz * fx * fx * pb_z

                     + 4.0 * pa_yy * fz * fx * pb_z

                     - fz * fga * pb_xxz

                     + 4.0 * fz * fx * pb_xxz

                     + 10.0 * pa_yy * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xyy_s_0(double fx,
                                  double pa_y,
                                  double pa_yy,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_xyy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_x

                     + 0.5 * pa_yy * pb_x * fx

                     + 2.0 * pa_y * fx * pb_xy

                     + 0.5 * fx * pb_xyy

                     + pa_yy * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xyy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yy,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_xyy,
                                  double r_0_0)
    {
        return r_0_0 * (4.5 * fx * fx * fz * pb_x

                     - 0.5 * fx * pb_x * fz * fgb

                     - 0.5 * fz * fga * pb_x * fx

                     - pa_yy * pb_x * fz * fgb

                     + 4.0 * pa_yy * fz * pb_x * fx

                     + 16.0 * pa_y * fz * fx * pb_xy

                     - fz * fga * pb_xyy

                     + 4.0 * fz * fx * pb_xyy

                     + 10.0 * pa_yy * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xyz_s_0(double fx,
                                  double pa_y,
                                  double pa_yy,
                                  double pb_xyz,
                                  double pb_xz,
                                  double s_0_0)
    {
        return s_0_0 * (pa_y * fx * pb_xz

                     + 0.5 * fx * pb_xyz

                     + pa_yy * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xyz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yy,
                                  double pb_xyz,
                                  double pb_xz,
                                  double r_0_0)
    {
        return r_0_0 * (8.0 * pa_y * fz * fx * pb_xz

                     - fz * fga * pb_xyz

                     + 4.0 * fz * fx * pb_xyz

                     + 10.0 * pa_yy * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xzz_s_0(double fx,
                                  double pa_yy,
                                  double pb_x,
                                  double pb_xzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_x

                     + 0.5 * pa_yy * pb_x * fx

                     + 0.5 * fx * pb_xzz

                     + pa_yy * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xzz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_yy,
                                  double pb_x,
                                  double pb_xzz,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pb_x * fz * fgb

                     - 0.5 * fz * fga * pb_x * fx

                     - pa_yy * pb_x * fz * fgb

                     + 1.5 * fz * fx * fx * pb_x

                     + 4.0 * pa_yy * fz * pb_x * fx

                     - fz * fga * pb_xzz

                     + 4.0 * fz * fx * pb_xzz

                     + 10.0 * pa_yy * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yyy_s_0(double fx,
                                  double pa_y,
                                  double pa_yy,
                                  double pb_y,
                                  double pb_yy,
                                  double pb_yyy,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * fx * fx

                     + 2.25 * fx * fx * pb_y

                     + 1.5 * pa_yy * pb_y * fx

                     + 3.0 * pa_y * fx * pb_yy

                     + 0.5 * fx * pb_yyy

                     + pa_yy * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yyy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yy,
                                  double pb_y,
                                  double pb_yy,
                                  double pb_yyy,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_y * fx * fz * fgb

                     + 9.0 * pa_y * fz * fx * fx

                     + 13.5 * fx * fx * fz * pb_y

                     - 1.5 * fx * pb_y * fz * fgb

                     - 1.5 * fz * fga * pb_y * fx

                     - 3.0 * pa_yy * pb_y * fz * fgb

                     + 12.0 * pa_yy * fz * pb_y * fx

                     + 24.0 * pa_y * fz * fx * pb_yy

                     - fz * fga * pb_yyy

                     + 4.0 * fz * fx * pb_yyy

                     + 10.0 * pa_yy * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yyz_s_0(double fx,
                                  double pa_y,
                                  double pa_yy,
                                  double pb_yyz,
                                  double pb_yz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_z

                     + 0.5 * pa_yy * fx * pb_z

                     + 2.0 * pa_y * fx * pb_yz

                     + 0.5 * fx * pb_yyz

                     + pa_yy * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yyz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yy,
                                  double pb_yyz,
                                  double pb_yz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (4.5 * fx * fx * fz * pb_z

                     - 0.5 * fx * fz * fgb * pb_z

                     - 0.5 * fz * fga * fx * pb_z

                     - pa_yy * fz * fgb * pb_z

                     + 4.0 * pa_yy * fz * fx * pb_z

                     + 16.0 * pa_y * fz * fx * pb_yz

                     - fz * fga * pb_yyz

                     + 4.0 * fz * fx * pb_yyz

                     + 10.0 * pa_yy * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yzz_s_0(double fx,
                                  double pa_y,
                                  double pa_yy,
                                  double pb_y,
                                  double pb_yzz,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx * fx

                     + 0.25 * fx * fx * pb_y

                     + 0.5 * pa_yy * pb_y * fx

                     + pa_y * fx * pb_zz

                     + 0.5 * fx * pb_yzz

                     + pa_yy * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yzz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yy,
                                  double pb_y,
                                  double pb_yzz,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_y * fx * fz * fgb

                     + 3.0 * pa_y * fz * fx * fx

                     - 0.5 * fx * pb_y * fz * fgb

                     - 0.5 * fz * fga * pb_y * fx

                     - pa_yy * pb_y * fz * fgb

                     + 1.5 * fz * fx * fx * pb_y

                     + 4.0 * pa_yy * fz * pb_y * fx

                     + 8.0 * pa_y * fz * fx * pb_zz

                     - fz * fga * pb_yzz

                     + 4.0 * fz * fx * pb_yzz

                     + 10.0 * pa_yy * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_zzz_s_0(double fx,
                                  double pa_yy,
                                  double pb_z,
                                  double pb_zzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_z

                     + 1.5 * pa_yy * pb_z * fx

                     + 0.5 * fx * pb_zzz

                     + pa_yy * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_zzz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_yy,
                                  double pb_z,
                                  double pb_zzz,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pb_z * fz * fgb

                     - 1.5 * fz * fga * pb_z * fx

                     - 3.0 * pa_yy * pb_z * fz * fgb

                     + 4.5 * fz * fx * fx * pb_z

                     + 12.0 * pa_yy * fz * pb_z * fx

                     - fz * fga * pb_zzz

                     + 4.0 * fz * fx * pb_zzz

                     + 10.0 * pa_yy * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxx_s_0(double fx,
                                  double pa_yz,
                                  double pb_x,
                                  double pb_xxx,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_yz * pb_x * fx

                     + pa_yz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxx_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_yz,
                                  double pb_x,
                                  double pb_xxx,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * pb_x * fz * fgb

                     + 12.0 * pa_yz * fz * pb_x * fx

                     + 10.0 * pa_yz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxy_s_0(double fx,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_xx,
                                  double pb_xxy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_z

                     + 0.5 * pa_yz * fx * pb_y

                     + 0.5 * fx * pa_z * pb_xx

                     + pa_yz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxy_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_xx,
                                  double pb_xxy,
                                  double pb_y,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_z * fz * fgb

                     + 1.5 * fx * fx * fz * pa_z

                     - pa_yz * fz * fgb * pb_y

                     + 4.0 * pa_yz * fz * fx * pb_y

                     + 4.0 * fx * fz * pa_z * pb_xx

                     + 10.0 * pa_yz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxz_s_0(double fx,
                                  double pa_y,
                                  double pa_yz,
                                  double pb_xx,
                                  double pb_xxz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_y * fx * fx

                     + 0.5 * pa_yz * fx * pb_z

                     + 0.5 * pa_y * fx * pb_xx

                     + pa_yz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yz,
                                  double pb_xx,
                                  double pb_xxz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_y * fx * fz * fgb

                     + 1.5 * pa_y * fz * fx * fx

                     - pa_yz * fz * fgb * pb_z

                     + 4.0 * pa_yz * fz * fx * pb_z

                     + 4.0 * pa_y * fz * fx * pb_xx

                     + 10.0 * pa_yz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xyy_s_0(double fx,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_xyy,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_yz * pb_x * fx

                     + fx * pa_z * pb_xy

                     + pa_yz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xyy_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_xyy,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_yz * pb_x * fz * fgb

                     + 4.0 * pa_yz * fz * pb_x * fx

                     + 8.0 * fx * fz * pa_z * pb_xy

                     + 10.0 * pa_yz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xyz_s_0(double fx,
                                  double pa_y,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_xyz,
                                  double pb_xz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_x

                     + 0.5 * pa_y * fx * pb_xy

                     + 0.5 * fx * pa_z * pb_xz

                     + pa_yz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xyz_r_0(double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_xyz,
                                  double pb_xz,
                                  double r_0_0)
    {
        return r_0_0 * (1.5 * fx * fx * fz * pb_x

                     + 4.0 * pa_y * fz * fx * pb_xy

                     + 4.0 * fx * fz * pa_z * pb_xz

                     + 10.0 * pa_yz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xzz_s_0(double fx,
                                  double pa_y,
                                  double pa_yz,
                                  double pb_x,
                                  double pb_xz,
                                  double pb_xzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_yz * pb_x * fx

                     + pa_y * fx * pb_xz

                     + pa_yz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xzz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yz,
                                  double pb_x,
                                  double pb_xz,
                                  double pb_xzz,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_yz * pb_x * fz * fgb

                     + 4.0 * pa_yz * fz * pb_x * fx

                     + 8.0 * pa_y * fz * fx * pb_xz

                     + 10.0 * pa_yz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yyy_s_0(double fx,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_y,
                                  double pb_yy,
                                  double pb_yyy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z

                     + 1.5 * pa_yz * pb_y * fx

                     + 1.5 * fx * pa_z * pb_yy

                     + pa_yz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yyy_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_y,
                                  double pb_yy,
                                  double pb_yyy,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_z * fz * fgb

                     + 4.5 * fx * fx * fz * pa_z

                     - 3.0 * pa_yz * pb_y * fz * fgb

                     + 12.0 * pa_yz * fz * pb_y * fx

                     + 12.0 * fx * fz * pa_z * pb_yy

                     + 10.0 * pa_yz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yyz_s_0(double fx,
                                  double pa_y,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_y,
                                  double pb_yy,
                                  double pb_yyz,
                                  double pb_yz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_y * fx * fx

                     + 0.5 * fx * fx * pb_y

                     + 0.5 * pa_yz * fx * pb_z

                     + 0.5 * pa_y * fx * pb_yy

                     + fx * pa_z * pb_yz

                     + pa_yz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yyz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_y,
                                  double pb_yy,
                                  double pb_yyz,
                                  double pb_yz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_y * fx * fz * fgb

                     + 1.5 * pa_y * fz * fx * fx

                     + 3.0 * fx * fx * fz * pb_y

                     - pa_yz * fz * fgb * pb_z

                     + 4.0 * pa_yz * fz * fx * pb_z

                     + 4.0 * pa_y * fz * fx * pb_yy

                     + 8.0 * fx * fz * pa_z * pb_yz

                     + 10.0 * pa_yz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yzz_s_0(double fx,
                                  double pa_y,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_y,
                                  double pb_yz,
                                  double pb_yzz,
                                  double pb_z,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_z

                     + 0.5 * fx * fx * pb_z

                     + 0.5 * pa_yz * pb_y * fx

                     + pa_y * fx * pb_yz

                     + 0.5 * fx * pa_z * pb_zz

                     + pa_yz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yzz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_y,
                                  double pb_yz,
                                  double pb_yzz,
                                  double pb_z,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_z * fz * fgb

                     + 1.5 * fx * fx * fz * pa_z

                     + 3.0 * fx * fx * fz * pb_z

                     - pa_yz * pb_y * fz * fgb

                     + 4.0 * pa_yz * fz * pb_y * fx

                     + 8.0 * pa_y * fz * fx * pb_yz

                     + 4.0 * fx * fz * pa_z * pb_zz

                     + 10.0 * pa_yz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_zzz_s_0(double fx,
                                  double pa_y,
                                  double pa_yz,
                                  double pb_z,
                                  double pb_zz,
                                  double pb_zzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx

                     + 1.5 * pa_yz * pb_z * fx

                     + 1.5 * pa_y * fx * pb_zz

                     + pa_yz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_zzz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yz,
                                  double pb_z,
                                  double pb_zz,
                                  double pb_zzz,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * fz * fgb

                     + 4.5 * pa_y * fz * fx * fx

                     - 3.0 * pa_yz * pb_z * fz * fgb

                     + 12.0 * pa_yz * fz * pb_z * fx

                     + 12.0 * pa_y * fz * fx * pb_zz

                     + 10.0 * pa_yz * fz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxx_s_0(double fx,
                                  double pa_zz,
                                  double pb_x,
                                  double pb_xxx,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_x

                     + 1.5 * pa_zz * pb_x * fx

                     + 0.5 * fx * pb_xxx

                     + pa_zz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxx_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_zz,
                                  double pb_x,
                                  double pb_xxx,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pb_x * fz * fgb

                     - 1.5 * fz * fga * pb_x * fx

                     - 3.0 * pa_zz * pb_x * fz * fgb

                     + 4.5 * fz * fx * fx * pb_x

                     + 12.0 * pa_zz * fz * pb_x * fx

                     - fz * fga * pb_xxx

                     + 4.0 * fz * fx * pb_xxx

                     + 10.0 * pa_zz * fz * pb_xxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxy_s_0(double fx,
                                  double pa_zz,
                                  double pb_xxy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_y

                     + 0.5 * pa_zz * fx * pb_y

                     + 0.5 * fx * pb_xxy

                     + pa_zz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_zz,
                                  double pb_xxy,
                                  double pb_y,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fz * fgb * pb_y

                     - 0.5 * fz * fga * fx * pb_y

                     - pa_zz * fz * fgb * pb_y

                     + 1.5 * fz * fx * fx * pb_y

                     + 4.0 * pa_zz * fz * fx * pb_y

                     - fz * fga * pb_xxy

                     + 4.0 * fz * fx * pb_xxy

                     + 10.0 * pa_zz * fz * pb_xxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxz_s_0(double fx,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_xx,
                                  double pb_xxz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_z * fx * fx

                     + 0.25 * fx * fx * pb_z

                     + 0.5 * pa_zz * fx * pb_z

                     + pa_z * fx * pb_xx

                     + 0.5 * fx * pb_xxz

                     + pa_zz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_xx,
                                  double pb_xxz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_z * fx * fz * fgb

                     + 3.0 * pa_z * fz * fx * fx

                     - 0.5 * fx * fz * fgb * pb_z

                     - 0.5 * fz * fga * fx * pb_z

                     - pa_zz * fz * fgb * pb_z

                     + 1.5 * fz * fx * fx * pb_z

                     + 4.0 * pa_zz * fz * fx * pb_z

                     + 8.0 * pa_z * fz * fx * pb_xx

                     - fz * fga * pb_xxz

                     + 4.0 * fz * fx * pb_xxz

                     + 10.0 * pa_zz * fz * pb_xxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xyy_s_0(double fx,
                                  double pa_zz,
                                  double pb_x,
                                  double pb_xyy,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_x

                     + 0.5 * pa_zz * pb_x * fx

                     + 0.5 * fx * pb_xyy

                     + pa_zz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xyy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_zz,
                                  double pb_x,
                                  double pb_xyy,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pb_x * fz * fgb

                     - 0.5 * fz * fga * pb_x * fx

                     - pa_zz * pb_x * fz * fgb

                     + 1.5 * fz * fx * fx * pb_x

                     + 4.0 * pa_zz * fz * pb_x * fx

                     - fz * fga * pb_xyy

                     + 4.0 * fz * fx * pb_xyy

                     + 10.0 * pa_zz * fz * pb_xyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xyz_s_0(double fx,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_xy,
                                  double pb_xyz,
                                  double s_0_0)
    {
        return s_0_0 * (pa_z * fx * pb_xy

                     + 0.5 * fx * pb_xyz

                     + pa_zz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xyz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_xy,
                                  double pb_xyz,
                                  double r_0_0)
    {
        return r_0_0 * (8.0 * pa_z * fz * fx * pb_xy

                     - fz * fga * pb_xyz

                     + 4.0 * fz * fx * pb_xyz

                     + 10.0 * pa_zz * fz * pb_xyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xzz_s_0(double fx,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_x,
                                  double pb_xz,
                                  double pb_xzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_x

                     + 0.5 * pa_zz * pb_x * fx

                     + 2.0 * pa_z * fx * pb_xz

                     + 0.5 * fx * pb_xzz

                     + pa_zz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xzz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_x,
                                  double pb_xz,
                                  double pb_xzz,
                                  double r_0_0)
    {
        return r_0_0 * (4.5 * fx * fx * fz * pb_x

                     - 0.5 * fx * pb_x * fz * fgb

                     - 0.5 * fz * fga * pb_x * fx

                     - pa_zz * pb_x * fz * fgb

                     + 4.0 * pa_zz * fz * pb_x * fx

                     + 16.0 * pa_z * fz * fx * pb_xz

                     - fz * fga * pb_xzz

                     + 4.0 * fz * fx * pb_xzz

                     + 10.0 * pa_zz * fz * pb_xzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yyy_s_0(double fx,
                                  double pa_zz,
                                  double pb_y,
                                  double pb_yyy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_y

                     + 1.5 * pa_zz * pb_y * fx

                     + 0.5 * fx * pb_yyy

                     + pa_zz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yyy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_zz,
                                  double pb_y,
                                  double pb_yyy,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pb_y * fz * fgb

                     - 1.5 * fz * fga * pb_y * fx

                     - 3.0 * pa_zz * pb_y * fz * fgb

                     + 4.5 * fz * fx * fx * pb_y

                     + 12.0 * pa_zz * fz * pb_y * fx

                     - fz * fga * pb_yyy

                     + 4.0 * fz * fx * pb_yyy

                     + 10.0 * pa_zz * fz * pb_yyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yyz_s_0(double fx,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_yy,
                                  double pb_yyz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_z * fx * fx

                     + 0.25 * fx * fx * pb_z

                     + 0.5 * pa_zz * fx * pb_z

                     + pa_z * fx * pb_yy

                     + 0.5 * fx * pb_yyz

                     + pa_zz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yyz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_yy,
                                  double pb_yyz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_z * fx * fz * fgb

                     + 3.0 * pa_z * fz * fx * fx

                     - 0.5 * fx * fz * fgb * pb_z

                     - 0.5 * fz * fga * fx * pb_z

                     - pa_zz * fz * fgb * pb_z

                     + 1.5 * fz * fx * fx * pb_z

                     + 4.0 * pa_zz * fz * fx * pb_z

                     + 8.0 * pa_z * fz * fx * pb_yy

                     - fz * fga * pb_yyz

                     + 4.0 * fz * fx * pb_yyz

                     + 10.0 * pa_zz * fz * pb_yyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yzz_s_0(double fx,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_y,
                                  double pb_yz,
                                  double pb_yzz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_y

                     + 0.5 * pa_zz * pb_y * fx

                     + 2.0 * pa_z * fx * pb_yz

                     + 0.5 * fx * pb_yzz

                     + pa_zz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yzz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_y,
                                  double pb_yz,
                                  double pb_yzz,
                                  double r_0_0)
    {
        return r_0_0 * (4.5 * fx * fx * fz * pb_y

                     - 0.5 * fx * pb_y * fz * fgb

                     - 0.5 * fz * fga * pb_y * fx

                     - pa_zz * pb_y * fz * fgb

                     + 4.0 * pa_zz * fz * pb_y * fx

                     + 16.0 * pa_z * fz * fx * pb_yz

                     - fz * fga * pb_yzz

                     + 4.0 * fz * fx * pb_yzz

                     + 10.0 * pa_zz * fz * pb_yzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_zzz_s_0(double fx,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_z,
                                  double pb_zz,
                                  double pb_zzz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * fx * fx

                     + 2.25 * fx * fx * pb_z

                     + 1.5 * pa_zz * pb_z * fx

                     + 3.0 * pa_z * fx * pb_zz

                     + 0.5 * fx * pb_zzz

                     + pa_zz * pb_zzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_zzz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_z,
                                  double pb_zz,
                                  double pb_zzz,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_z * fx * fz * fgb

                     + 9.0 * pa_z * fz * fx * fx

                     + 13.5 * fx * fx * fz * pb_z

                     - 1.5 * fx * pb_z * fz * fgb

                     - 1.5 * fz * fga * pb_z * fx

                     - 3.0 * pa_zz * pb_z * fz * fgb

                     + 12.0 * pa_zz * fz * pb_z * fx

                     + 24.0 * pa_z * fz * fx * pb_zz

                     - fz * fga * pb_zzz

                     + 4.0 * fz * fx * pb_zzz

                     + 10.0 * pa_zz * fz * pb_zzz);

    }

    // SIMD elementary functions for (F|T|D) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xx_s_0(double fx,
                                  double pa_x,
                                  double pa_xx,
                                  double pa_xxx,
                                  double pb_x,
                                  double pb_xx,
                                  double s_0_0)
    {
        return s_0_0 * (2.25 * pa_x * fx * fx

                     + 1.5 * fx * fx * pb_x

                     + 0.5 * pa_xxx * fx

                     + 3.0 * pa_xx * fx * pb_x

                     + 1.5 * pa_x * fx * pb_xx

                     + pa_xxx * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xx_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xx,
                                  double pa_xxx,
                                  double pb_x,
                                  double pb_xx,
                                  double r_0_0)
    {
        return r_0_0 * (13.5 * pa_x * fx * fx * fz

                     - 1.5 * pa_x * fx * fz * fgb

                     - 1.5 * pa_x * fz * fga * fx

                     - 3.0 * fx * fz * fga * pb_x

                     - pa_xxx * fz * fgb

                     + 9.0 * fx * fx * fz * pb_x

                     + 4.0 * pa_xxx * fz * fx

                     + 24.0 * pa_xx * fz * fx * pb_x

                     - 3.0 * pa_x * fz * fga * pb_xx

                     + 12.0 * pa_x * fz * fx * pb_xx

                     + 10.0 * pa_xxx * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xy_s_0(double fx,
                                  double pa_x,
                                  double pa_xx,
                                  double pa_xxx,
                                  double pb_xy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_y

                     + 1.5 * pa_xx * fx * pb_y

                     + 1.5 * pa_x * fx * pb_xy

                     + pa_xxx * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xy_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xx,
                                  double pa_xxx,
                                  double pb_xy,
                                  double pb_y,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fz * fga * pb_y

                     + 4.5 * fx * fx * fz * pb_y

                     + 12.0 * pa_xx * fz * fx * pb_y

                     - 3.0 * pa_x * fz * fga * pb_xy

                     + 12.0 * pa_x * fz * fx * pb_xy

                     + 10.0 * pa_xxx * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xz_s_0(double fx,
                                  double pa_x,
                                  double pa_xx,
                                  double pa_xxx,
                                  double pb_xz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_z

                     + 1.5 * pa_xx * fx * pb_z

                     + 1.5 * pa_x * fx * pb_xz

                     + pa_xxx * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_xz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xx,
                                  double pa_xxx,
                                  double pb_xz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fz * fga * pb_z

                     + 4.5 * fx * fx * fz * pb_z

                     + 12.0 * pa_xx * fz * fx * pb_z

                     - 3.0 * pa_x * fz * fga * pb_xz

                     + 12.0 * pa_x * fz * fx * pb_xz

                     + 10.0 * pa_xxx * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_yy_s_0(double fx,
                                  double pa_x,
                                  double pa_xxx,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx

                     + 0.5 * pa_xxx * fx

                     + 1.5 * pa_x * fx * pb_yy

                     + pa_xxx * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_yy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xxx,
                                  double pb_yy,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fz * fgb

                     - 1.5 * pa_x * fz * fga * fx

                     - pa_xxx * fz * fgb

                     + 4.5 * pa_x * fz * fx * fx

                     + 4.0 * pa_xxx * fz * fx

                     - 3.0 * pa_x * fz * fga * pb_yy

                     + 12.0 * pa_x * fz * fx * pb_yy

                     + 10.0 * pa_xxx * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_yz_s_0(double fx,
                                  double pa_x,
                                  double pa_xxx,
                                  double pb_yz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * fx * pb_yz

                     + pa_xxx * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_yz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xxx,
                                  double pb_yz,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_x * fz * fga * pb_yz

                     + 12.0 * pa_x * fz * fx * pb_yz

                     + 10.0 * pa_xxx * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_zz_s_0(double fx,
                                  double pa_x,
                                  double pa_xxx,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx

                     + 0.5 * pa_xxx * fx

                     + 1.5 * pa_x * fx * pb_zz

                     + pa_xxx * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxx_zz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xxx,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fz * fgb

                     - 1.5 * pa_x * fz * fga * fx

                     - pa_xxx * fz * fgb

                     + 4.5 * pa_x * fz * fx * fx

                     + 4.0 * pa_xxx * fz * fx

                     - 3.0 * pa_x * fz * fga * pb_zz

                     + 12.0 * pa_x * fz * fx * pb_zz

                     + 10.0 * pa_xxx * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xx_s_0(double fx,
                                  double pa_xxy,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_x,
                                  double pb_xx,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_y

                     + 0.5 * pa_xxy * fx

                     + 2.0 * pa_xy * fx * pb_x

                     + 0.5 * fx * pa_y * pb_xx

                     + pa_xxy * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xx_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xxy,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_x,
                                  double pb_xx,
                                  double r_0_0)
    {
        return r_0_0 * (4.5 * fx * fx * fz * pa_y

                     - 0.5 * fx * pa_y * fz * fgb

                     - 0.5 * fz * fga * pa_y * fx

                     - pa_xxy * fz * fgb

                     + 4.0 * pa_xxy * fz * fx

                     + 16.0 * pa_xy * fx * fz * pb_x

                     - fz * fga * pa_y * pb_xx

                     + 4.0 * fx * fz * pa_y * pb_xx

                     + 10.0 * pa_xxy * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xy_s_0(double fx,
                                  double pa_x,
                                  double pa_xx,
                                  double pa_xxy,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * fx

                     + 0.25 * fx * fx * pb_x

                     + 0.5 * pa_xx * fx * pb_x

                     + pa_xy * fx * pb_y

                     + 0.5 * fx * pa_y * pb_xy

                     + pa_xxy * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xy_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xx,
                                  double pa_xxy,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_y,
                                  double r_0_0)
    {
        return r_0_0 * (3.0 * pa_x * fx * fx * fz

                     - 0.5 * fz * fga * fx * pb_x

                     + 1.5 * fx * fx * fz * pb_x

                     + 4.0 * pa_xx * fz * fx * pb_x

                     + 8.0 * pa_xy * fx * fz * pb_y

                     - fz * fga * pa_y * pb_xy

                     + 4.0 * fx * fz * pa_y * pb_xy

                     + 10.0 * pa_xxy * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xz_s_0(double fx,
                                  double pa_xxy,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_xz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (pa_xy * fx * pb_z

                     + 0.5 * fx * pa_y * pb_xz

                     + pa_xxy * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_xz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_xxy,
                                  double pa_xy,
                                  double pa_y,
                                  double pb_xz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (8.0 * pa_xy * fx * fz * pb_z

                     - fz * fga * pa_y * pb_xz

                     + 4.0 * fx * fz * pa_y * pb_xz

                     + 10.0 * pa_xxy * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_yy_s_0(double fx,
                                  double pa_xx,
                                  double pa_xxy,
                                  double pa_y,
                                  double pb_y,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_y

                     + 0.5 * fx * fx * pb_y

                     + 0.5 * pa_xxy * fx

                     + pa_xx * fx * pb_y

                     + 0.5 * fx * pa_y * pb_yy

                     + pa_xxy * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_yy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xx,
                                  double pa_xxy,
                                  double pa_y,
                                  double pb_y,
                                  double pb_yy,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_y * fz * fgb

                     - 0.5 * fz * fga * pa_y * fx

                     - fz * fga * fx * pb_y

                     - pa_xxy * fz * fgb

                     + 1.5 * fx * fx * fz * pa_y

                     + 3.0 * fx * fx * fz * pb_y

                     + 4.0 * pa_xxy * fz * fx

                     + 8.0 * pa_xx * fz * fx * pb_y

                     - fz * fga * pa_y * pb_yy

                     + 4.0 * fx * fz * pa_y * pb_yy

                     + 10.0 * pa_xxy * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_yz_s_0(double fx,
                                  double pa_xx,
                                  double pa_xxy,
                                  double pa_y,
                                  double pb_yz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_z

                     + 0.5 * pa_xx * fx * pb_z

                     + 0.5 * fx * pa_y * pb_yz

                     + pa_xxy * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_yz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_xx,
                                  double pa_xxy,
                                  double pa_y,
                                  double pb_yz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fz * fga * fx * pb_z

                     + 1.5 * fx * fx * fz * pb_z

                     + 4.0 * pa_xx * fz * fx * pb_z

                     - fz * fga * pa_y * pb_yz

                     + 4.0 * fx * fz * pa_y * pb_yz

                     + 10.0 * pa_xxy * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_zz_s_0(double fx,
                                  double pa_xxy,
                                  double pa_y,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_y

                     + 0.5 * pa_xxy * fx

                     + 0.5 * fx * pa_y * pb_zz

                     + pa_xxy * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxy_zz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xxy,
                                  double pa_y,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_y * fz * fgb

                     - 0.5 * fz * fga * pa_y * fx

                     - pa_xxy * fz * fgb

                     + 1.5 * fx * fx * fz * pa_y

                     + 4.0 * pa_xxy * fz * fx

                     - fz * fga * pa_y * pb_zz

                     + 4.0 * fx * fz * pa_y * pb_zz

                     + 10.0 * pa_xxy * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xx_s_0(double fx,
                                  double pa_xxz,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xx,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z

                     + 0.5 * pa_xxz * fx

                     + 2.0 * pa_xz * fx * pb_x

                     + 0.5 * fx * pa_z * pb_xx

                     + pa_xxz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xx_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xxz,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xx,
                                  double r_0_0)
    {
        return r_0_0 * (4.5 * fx * fx * fz * pa_z

                     - 0.5 * fx * pa_z * fz * fgb

                     - 0.5 * fz * fga * pa_z * fx

                     - pa_xxz * fz * fgb

                     + 4.0 * pa_xxz * fz * fx

                     + 16.0 * pa_xz * fx * fz * pb_x

                     - fz * fga * pa_z * pb_xx

                     + 4.0 * fx * fz * pa_z * pb_xx

                     + 10.0 * pa_xxz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xy_s_0(double fx,
                                  double pa_xxz,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_xy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (pa_xz * fx * pb_y

                     + 0.5 * fx * pa_z * pb_xy

                     + pa_xxz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xy_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_xxz,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_xy,
                                  double pb_y,
                                  double r_0_0)
    {
        return r_0_0 * (8.0 * pa_xz * fx * fz * pb_y

                     - fz * fga * pa_z * pb_xy

                     + 4.0 * fx * fz * pa_z * pb_xy

                     + 10.0 * pa_xxz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xz_s_0(double fx,
                                  double pa_x,
                                  double pa_xx,
                                  double pa_xxz,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * fx

                     + 0.25 * fx * fx * pb_x

                     + 0.5 * pa_xx * fx * pb_x

                     + pa_xz * fx * pb_z

                     + 0.5 * fx * pa_z * pb_xz

                     + pa_xxz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_xz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xx,
                                  double pa_xxz,
                                  double pa_xz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (3.0 * pa_x * fx * fx * fz

                     - 0.5 * fz * fga * fx * pb_x

                     + 1.5 * fx * fx * fz * pb_x

                     + 4.0 * pa_xx * fz * fx * pb_x

                     + 8.0 * pa_xz * fx * fz * pb_z

                     - fz * fga * pa_z * pb_xz

                     + 4.0 * fx * fz * pa_z * pb_xz

                     + 10.0 * pa_xxz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_yy_s_0(double fx,
                                  double pa_xxz,
                                  double pa_z,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_z

                     + 0.5 * pa_xxz * fx

                     + 0.5 * fx * pa_z * pb_yy

                     + pa_xxz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_yy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xxz,
                                  double pa_z,
                                  double pb_yy,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_z * fz * fgb

                     - 0.5 * fz * fga * pa_z * fx

                     - pa_xxz * fz * fgb

                     + 1.5 * fx * fx * fz * pa_z

                     + 4.0 * pa_xxz * fz * fx

                     - fz * fga * pa_z * pb_yy

                     + 4.0 * fx * fz * pa_z * pb_yy

                     + 10.0 * pa_xxz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_yz_s_0(double fx,
                                  double pa_xx,
                                  double pa_xxz,
                                  double pa_z,
                                  double pb_y,
                                  double pb_yz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_y

                     + 0.5 * pa_xx * fx * pb_y

                     + 0.5 * fx * pa_z * pb_yz

                     + pa_xxz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_yz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_xx,
                                  double pa_xxz,
                                  double pa_z,
                                  double pb_y,
                                  double pb_yz,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fz * fga * fx * pb_y

                     + 1.5 * fx * fx * fz * pb_y

                     + 4.0 * pa_xx * fz * fx * pb_y

                     - fz * fga * pa_z * pb_yz

                     + 4.0 * fx * fz * pa_z * pb_yz

                     + 10.0 * pa_xxz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_zz_s_0(double fx,
                                  double pa_xx,
                                  double pa_xxz,
                                  double pa_z,
                                  double pb_z,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_z

                     + 0.5 * fx * fx * pb_z

                     + 0.5 * pa_xxz * fx

                     + pa_xx * fx * pb_z

                     + 0.5 * fx * pa_z * pb_zz

                     + pa_xxz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxz_zz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xx,
                                  double pa_xxz,
                                  double pa_z,
                                  double pb_z,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_z * fz * fgb

                     - 0.5 * fz * fga * pa_z * fx

                     - fz * fga * fx * pb_z

                     - pa_xxz * fz * fgb

                     + 1.5 * fx * fx * fz * pa_z

                     + 3.0 * fx * fx * fz * pb_z

                     + 4.0 * pa_xxz * fz * fx

                     + 8.0 * pa_xx * fz * fx * pb_z

                     - fz * fga * pa_z * pb_zz

                     + 4.0 * fx * fz * pa_z * pb_zz

                     + 10.0 * pa_xxz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xx_s_0(double fx,
                                  double pa_x,
                                  double pa_xyy,
                                  double pa_yy,
                                  double pb_x,
                                  double pb_xx,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx

                     + 0.5 * fx * fx * pb_x

                     + 0.5 * pa_xyy * fx

                     + fx * pa_yy * pb_x

                     + 0.5 * pa_x * fx * pb_xx

                     + pa_xyy * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xx_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xyy,
                                  double pa_yy,
                                  double pb_x,
                                  double pb_xx,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fz * fgb

                     - 0.5 * pa_x * fz * fga * fx

                     - fx * fz * fga * pb_x

                     - pa_xyy * fz * fgb

                     + 1.5 * pa_x * fz * fx * fx

                     + 3.0 * fx * fx * fz * pb_x

                     + 4.0 * pa_xyy * fz * fx

                     + 8.0 * fx * pa_yy * fz * pb_x

                     - pa_x * fz * fga * pb_xx

                     + 4.0 * pa_x * fz * fx * pb_xx

                     + 10.0 * pa_xyy * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xy_s_0(double fx,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_xyy,
                                  double pa_y,
                                  double pa_yy,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * fx * fx * pa_y

                     + 0.25 * fx * fx * pb_y

                     + pa_xy * fx * pb_x

                     + 0.5 * fx * pa_yy * pb_y

                     + 0.5 * pa_x * fx * pb_xy

                     + pa_xyy * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xy_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_xyy,
                                  double pa_y,
                                  double pa_yy,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_y,
                                  double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fx * pa_y * fz

                     - 0.5 * fx * fz * fga * pb_y

                     + 1.5 * fx * fx * fz * pb_y

                     + 8.0 * pa_xy * fz * fx * pb_x

                     + 4.0 * fx * pa_yy * fz * pb_y

                     - pa_x * fz * fga * pb_xy

                     + 4.0 * pa_x * fz * fx * pb_xy

                     + 10.0 * pa_xyy * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xz_s_0(double fx,
                                  double pa_x,
                                  double pa_xyy,
                                  double pa_yy,
                                  double pb_xz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_z

                     + 0.5 * fx * pa_yy * pb_z

                     + 0.5 * pa_x * fx * pb_xz

                     + pa_xyy * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_xz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xyy,
                                  double pa_yy,
                                  double pb_xz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fz * fga * pb_z

                     + 1.5 * fx * fx * fz * pb_z

                     + 4.0 * fx * pa_yy * fz * pb_z

                     - pa_x * fz * fga * pb_xz

                     + 4.0 * pa_x * fz * fx * pb_xz

                     + 10.0 * pa_xyy * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_yy_s_0(double fx,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_xyy,
                                  double pb_y,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx

                     + 0.5 * pa_xyy * fx

                     + 2.0 * pa_xy * fx * pb_y

                     + 0.5 * pa_x * fx * pb_yy

                     + pa_xyy * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_yy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_xyy,
                                  double pb_y,
                                  double pb_yy,
                                  double r_0_0)
    {
        return r_0_0 * (4.5 * pa_x * fx * fx * fz

                     - 0.5 * pa_x * fx * fz * fgb

                     - 0.5 * pa_x * fz * fga * fx

                     - pa_xyy * fz * fgb

                     + 4.0 * pa_xyy * fz * fx

                     + 16.0 * pa_xy * fz * fx * pb_y

                     - pa_x * fz * fga * pb_yy

                     + 4.0 * pa_x * fz * fx * pb_yy

                     + 10.0 * pa_xyy * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_yz_s_0(double fx,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_xyy,
                                  double pb_yz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (pa_xy * fx * pb_z

                     + 0.5 * pa_x * fx * pb_yz

                     + pa_xyy * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_yz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_xyy,
                                  double pb_yz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (8.0 * pa_xy * fz * fx * pb_z

                     - pa_x * fz * fga * pb_yz

                     + 4.0 * pa_x * fz * fx * pb_yz

                     + 10.0 * pa_xyy * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_zz_s_0(double fx,
                                  double pa_x,
                                  double pa_xyy,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx

                     + 0.5 * pa_xyy * fx

                     + 0.5 * pa_x * fx * pb_zz

                     + pa_xyy * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyy_zz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xyy,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fz * fgb

                     - 0.5 * pa_x * fz * fga * fx

                     - pa_xyy * fz * fgb

                     + 1.5 * pa_x * fz * fx * fx

                     + 4.0 * pa_xyy * fz * fx

                     - pa_x * fz * fga * pb_zz

                     + 4.0 * pa_x * fz * fx * pb_zz

                     + 10.0 * pa_xyy * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xx_s_0(double fx,
                                  double pa_xyz,
                                  double pa_yz,
                                  double pb_x,
                                  double pb_xx,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xyz * fx

                     + fx * pa_yz * pb_x

                     + pa_xyz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xx_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xyz,
                                  double pa_yz,
                                  double pb_x,
                                  double pb_xx,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_xyz * fz * fgb

                     + 4.0 * pa_xyz * fz * fx

                     + 8.0 * fx * pa_yz * fz * pb_x

                     + 10.0 * pa_xyz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xy_s_0(double fx,
                                  double pa_xyz,
                                  double pa_xz,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_z

                     + 0.5 * pa_xz * fx * pb_x

                     + 0.5 * fx * pa_yz * pb_y

                     + pa_xyz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xy_r_0(double fx,
                                  double fz,
                                  double pa_xyz,
                                  double pa_xz,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xy,
                                  double pb_y,
                                  double r_0_0)
    {
        return r_0_0 * (1.5 * fx * fx * fz * pa_z

                     + 4.0 * pa_xz * fx * fz * pb_x

                     + 4.0 * fx * pa_yz * fz * pb_y

                     + 10.0 * pa_xyz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xz_s_0(double fx,
                                  double pa_xy,
                                  double pa_xyz,
                                  double pa_y,
                                  double pa_yz,
                                  double pb_x,
                                  double pb_xz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_y

                     + 0.5 * pa_xy * fx * pb_x

                     + 0.5 * fx * pa_yz * pb_z

                     + pa_xyz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_xz_r_0(double fx,
                                  double fz,
                                  double pa_xy,
                                  double pa_xyz,
                                  double pa_y,
                                  double pa_yz,
                                  double pb_x,
                                  double pb_xz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (1.5 * fx * fx * pa_y * fz

                     + 4.0 * pa_xy * fz * fx * pb_x

                     + 4.0 * fx * pa_yz * fz * pb_z

                     + 10.0 * pa_xyz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_yy_s_0(double fx,
                                  double pa_xyz,
                                  double pa_xz,
                                  double pb_y,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xyz * fx

                     + pa_xz * fx * pb_y

                     + pa_xyz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_yy_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xyz,
                                  double pa_xz,
                                  double pb_y,
                                  double pb_yy,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_xyz * fz * fgb

                     + 4.0 * pa_xyz * fz * fx

                     + 8.0 * pa_xz * fx * fz * pb_y

                     + 10.0 * pa_xyz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_yz_s_0(double fx,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_xyz,
                                  double pa_xz,
                                  double pb_y,
                                  double pb_yz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx

                     + 0.5 * pa_xy * fx * pb_y

                     + 0.5 * pa_xz * fx * pb_z

                     + pa_xyz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_yz_r_0(double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xy,
                                  double pa_xyz,
                                  double pa_xz,
                                  double pb_y,
                                  double pb_yz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (1.5 * pa_x * fx * fx * fz

                     + 4.0 * pa_xy * fz * fx * pb_y

                     + 4.0 * pa_xz * fx * fz * pb_z

                     + 10.0 * pa_xyz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_zz_s_0(double fx,
                                  double pa_xy,
                                  double pa_xyz,
                                  double pb_z,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xyz * fx

                     + pa_xy * fx * pb_z

                     + pa_xyz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyz_zz_r_0(double fgb,
                                  double fx,
                                  double fz,
                                  double pa_xy,
                                  double pa_xyz,
                                  double pb_z,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- pa_xyz * fz * fgb

                     + 4.0 * pa_xyz * fz * fx

                     + 8.0 * pa_xy * fz * fx * pb_z

                     + 10.0 * pa_xyz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xx_s_0(double fx,
                                  double pa_x,
                                  double pa_xzz,
                                  double pa_zz,
                                  double pb_x,
                                  double pb_xx,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx

                     + 0.5 * fx * fx * pb_x

                     + 0.5 * pa_xzz * fx

                     + fx * pa_zz * pb_x

                     + 0.5 * pa_x * fx * pb_xx

                     + pa_xzz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xx_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xzz,
                                  double pa_zz,
                                  double pb_x,
                                  double pb_xx,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fz * fgb

                     - 0.5 * pa_x * fz * fga * fx

                     - fx * fz * fga * pb_x

                     - pa_xzz * fz * fgb

                     + 1.5 * pa_x * fz * fx * fx

                     + 3.0 * fx * fx * fz * pb_x

                     + 4.0 * pa_xzz * fz * fx

                     + 8.0 * fx * pa_zz * fz * pb_x

                     - pa_x * fz * fga * pb_xx

                     + 4.0 * pa_x * fz * fx * pb_xx

                     + 10.0 * pa_xzz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xy_s_0(double fx,
                                  double pa_x,
                                  double pa_xzz,
                                  double pa_zz,
                                  double pb_xy,
                                  double pb_y,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_y

                     + 0.5 * fx * pa_zz * pb_y

                     + 0.5 * pa_x * fx * pb_xy

                     + pa_xzz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xy_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xzz,
                                  double pa_zz,
                                  double pb_xy,
                                  double pb_y,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fz * fga * pb_y

                     + 1.5 * fx * fx * fz * pb_y

                     + 4.0 * fx * pa_zz * fz * pb_y

                     - pa_x * fz * fga * pb_xy

                     + 4.0 * pa_x * fz * fx * pb_xy

                     + 10.0 * pa_xzz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xz_s_0(double fx,
                                  double pa_x,
                                  double pa_xz,
                                  double pa_xzz,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_x,
                                  double pb_xz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * fx * fx * pa_z

                     + 0.25 * fx * fx * pb_z

                     + pa_xz * fx * pb_x

                     + 0.5 * fx * pa_zz * pb_z

                     + 0.5 * pa_x * fx * pb_xz

                     + pa_xzz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_xz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xz,
                                  double pa_xzz,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_x,
                                  double pb_xz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fx * pa_z * fz

                     - 0.5 * fx * fz * fga * pb_z

                     + 1.5 * fx * fx * fz * pb_z

                     + 8.0 * pa_xz * fz * fx * pb_x

                     + 4.0 * fx * pa_zz * fz * pb_z

                     - pa_x * fz * fga * pb_xz

                     + 4.0 * pa_x * fz * fx * pb_xz

                     + 10.0 * pa_xzz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_yy_s_0(double fx,
                                  double pa_x,
                                  double pa_xzz,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx

                     + 0.5 * pa_xzz * fx

                     + 0.5 * pa_x * fx * pb_yy

                     + pa_xzz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_yy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xzz,
                                  double pb_yy,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fz * fgb

                     - 0.5 * pa_x * fz * fga * fx

                     - pa_xzz * fz * fgb

                     + 1.5 * pa_x * fz * fx * fx

                     + 4.0 * pa_xzz * fz * fx

                     - pa_x * fz * fga * pb_yy

                     + 4.0 * pa_x * fz * fx * pb_yy

                     + 10.0 * pa_xzz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_yz_s_0(double fx,
                                  double pa_x,
                                  double pa_xz,
                                  double pa_xzz,
                                  double pb_y,
                                  double pb_yz,
                                  double s_0_0)
    {
        return s_0_0 * (pa_xz * fx * pb_y

                     + 0.5 * pa_x * fx * pb_yz

                     + pa_xzz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_yz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xz,
                                  double pa_xzz,
                                  double pb_y,
                                  double pb_yz,
                                  double r_0_0)
    {
        return r_0_0 * (8.0 * pa_xz * fz * fx * pb_y

                     - pa_x * fz * fga * pb_yz

                     + 4.0 * pa_x * fz * fx * pb_yz

                     + 10.0 * pa_xzz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_zz_s_0(double fx,
                                  double pa_x,
                                  double pa_xz,
                                  double pa_xzz,
                                  double pb_z,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx

                     + 0.5 * pa_xzz * fx

                     + 2.0 * pa_xz * fx * pb_z

                     + 0.5 * pa_x * fx * pb_zz

                     + pa_xzz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzz_zz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_x,
                                  double pa_xz,
                                  double pa_xzz,
                                  double pb_z,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (4.5 * pa_x * fx * fx * fz

                     - 0.5 * pa_x * fx * fz * fgb

                     - 0.5 * pa_x * fz * fga * fx

                     - pa_xzz * fz * fgb

                     + 4.0 * pa_xzz * fz * fx

                     + 16.0 * pa_xz * fz * fx * pb_z

                     - pa_x * fz * fga * pb_zz

                     + 4.0 * pa_x * fz * fx * pb_zz

                     + 10.0 * pa_xzz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xx_s_0(double fx,
                                  double pa_y,
                                  double pa_yyy,
                                  double pb_xx,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx

                     + 0.5 * pa_yyy * fx

                     + 1.5 * pa_y * fx * pb_xx

                     + pa_yyy * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xx_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yyy,
                                  double pb_xx,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * fz * fgb

                     - 1.5 * pa_y * fz * fga * fx

                     - pa_yyy * fz * fgb

                     + 4.5 * pa_y * fz * fx * fx

                     + 4.0 * pa_yyy * fz * fx

                     - 3.0 * pa_y * fz * fga * pb_xx

                     + 12.0 * pa_y * fz * fx * pb_xx

                     + 10.0 * pa_yyy * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xy_s_0(double fx,
                                  double pa_y,
                                  double pa_yy,
                                  double pa_yyy,
                                  double pb_x,
                                  double pb_xy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_x

                     + 1.5 * pa_yy * fx * pb_x

                     + 1.5 * pa_y * fx * pb_xy

                     + pa_yyy * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xy_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yy,
                                  double pa_yyy,
                                  double pb_x,
                                  double pb_xy,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fz * fga * pb_x

                     + 4.5 * fx * fx * fz * pb_x

                     + 12.0 * pa_yy * fz * fx * pb_x

                     - 3.0 * pa_y * fz * fga * pb_xy

                     + 12.0 * pa_y * fz * fx * pb_xy

                     + 10.0 * pa_yyy * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xz_s_0(double fx,
                                  double pa_y,
                                  double pa_yyy,
                                  double pb_xz,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * fx * pb_xz

                     + pa_yyy * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_xz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yyy,
                                  double pb_xz,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_y * fz * fga * pb_xz

                     + 12.0 * pa_y * fz * fx * pb_xz

                     + 10.0 * pa_yyy * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_yy_s_0(double fx,
                                  double pa_y,
                                  double pa_yy,
                                  double pa_yyy,
                                  double pb_y,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (2.25 * pa_y * fx * fx

                     + 1.5 * fx * fx * pb_y

                     + 0.5 * pa_yyy * fx

                     + 3.0 * pa_yy * fx * pb_y

                     + 1.5 * pa_y * fx * pb_yy

                     + pa_yyy * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_yy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yy,
                                  double pa_yyy,
                                  double pb_y,
                                  double pb_yy,
                                  double r_0_0)
    {
        return r_0_0 * (13.5 * pa_y * fx * fx * fz

                     - 1.5 * pa_y * fx * fz * fgb

                     - 1.5 * pa_y * fz * fga * fx

                     - 3.0 * fx * fz * fga * pb_y

                     - pa_yyy * fz * fgb

                     + 9.0 * fx * fx * fz * pb_y

                     + 4.0 * pa_yyy * fz * fx

                     + 24.0 * pa_yy * fz * fx * pb_y

                     - 3.0 * pa_y * fz * fga * pb_yy

                     + 12.0 * pa_y * fz * fx * pb_yy

                     + 10.0 * pa_yyy * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_yz_s_0(double fx,
                                  double pa_y,
                                  double pa_yy,
                                  double pa_yyy,
                                  double pb_yz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_z

                     + 1.5 * pa_yy * fx * pb_z

                     + 1.5 * pa_y * fx * pb_yz

                     + pa_yyy * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_yz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yy,
                                  double pa_yyy,
                                  double pb_yz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fz * fga * pb_z

                     + 4.5 * fx * fx * fz * pb_z

                     + 12.0 * pa_yy * fz * fx * pb_z

                     - 3.0 * pa_y * fz * fga * pb_yz

                     + 12.0 * pa_y * fz * fx * pb_yz

                     + 10.0 * pa_yyy * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_zz_s_0(double fx,
                                  double pa_y,
                                  double pa_yyy,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx

                     + 0.5 * pa_yyy * fx

                     + 1.5 * pa_y * fx * pb_zz

                     + pa_yyy * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyy_zz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yyy,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * fz * fgb

                     - 1.5 * pa_y * fz * fga * fx

                     - pa_yyy * fz * fgb

                     + 4.5 * pa_y * fz * fx * fx

                     + 4.0 * pa_yyy * fz * fx

                     - 3.0 * pa_y * fz * fga * pb_zz

                     + 12.0 * pa_y * fz * fx * pb_zz

                     + 10.0 * pa_yyy * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xx_s_0(double fx,
                                  double pa_yyz,
                                  double pa_z,
                                  double pb_xx,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_z

                     + 0.5 * pa_yyz * fx

                     + 0.5 * fx * pa_z * pb_xx

                     + pa_yyz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xx_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_yyz,
                                  double pa_z,
                                  double pb_xx,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_z * fz * fgb

                     - 0.5 * fz * fga * pa_z * fx

                     - pa_yyz * fz * fgb

                     + 1.5 * fx * fx * fz * pa_z

                     + 4.0 * pa_yyz * fz * fx

                     - fz * fga * pa_z * pb_xx

                     + 4.0 * fx * fz * pa_z * pb_xx

                     + 10.0 * pa_yyz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xy_s_0(double fx,
                                  double pa_yyz,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xy,
                                  double s_0_0)
    {
        return s_0_0 * (pa_yz * fx * pb_x

                     + 0.5 * fx * pa_z * pb_xy

                     + pa_yyz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xy_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_yyz,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xy,
                                  double r_0_0)
    {
        return r_0_0 * (8.0 * pa_yz * fx * fz * pb_x

                     - fz * fga * pa_z * pb_xy

                     + 4.0 * fx * fz * pa_z * pb_xy

                     + 10.0 * pa_yyz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xz_s_0(double fx,
                                  double pa_yy,
                                  double pa_yyz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_x

                     + 0.5 * pa_yy * fx * pb_x

                     + 0.5 * fx * pa_z * pb_xz

                     + pa_yyz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_xz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_yy,
                                  double pa_yyz,
                                  double pa_z,
                                  double pb_x,
                                  double pb_xz,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fz * fga * fx * pb_x

                     + 1.5 * fx * fx * fz * pb_x

                     + 4.0 * pa_yy * fz * fx * pb_x

                     - fz * fga * pa_z * pb_xz

                     + 4.0 * fx * fz * pa_z * pb_xz

                     + 10.0 * pa_yyz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_yy_s_0(double fx,
                                  double pa_yyz,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_y,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z

                     + 0.5 * pa_yyz * fx

                     + 2.0 * pa_yz * fx * pb_y

                     + 0.5 * fx * pa_z * pb_yy

                     + pa_yyz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_yy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_yyz,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_y,
                                  double pb_yy,
                                  double r_0_0)
    {
        return r_0_0 * (4.5 * fx * fx * fz * pa_z

                     - 0.5 * fx * pa_z * fz * fgb

                     - 0.5 * fz * fga * pa_z * fx

                     - pa_yyz * fz * fgb

                     + 4.0 * pa_yyz * fz * fx

                     + 16.0 * pa_yz * fx * fz * pb_y

                     - fz * fga * pa_z * pb_yy

                     + 4.0 * fx * fz * pa_z * pb_yy

                     + 10.0 * pa_yyz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_yz_s_0(double fx,
                                  double pa_y,
                                  double pa_yy,
                                  double pa_yyz,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_y,
                                  double pb_yz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx * fx

                     + 0.25 * fx * fx * pb_y

                     + 0.5 * pa_yy * fx * pb_y

                     + pa_yz * fx * pb_z

                     + 0.5 * fx * pa_z * pb_yz

                     + pa_yyz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_yz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yy,
                                  double pa_yyz,
                                  double pa_yz,
                                  double pa_z,
                                  double pb_y,
                                  double pb_yz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (3.0 * pa_y * fx * fx * fz

                     - 0.5 * fz * fga * fx * pb_y

                     + 1.5 * fx * fx * fz * pb_y

                     + 4.0 * pa_yy * fz * fx * pb_y

                     + 8.0 * pa_yz * fx * fz * pb_z

                     - fz * fga * pa_z * pb_yz

                     + 4.0 * fx * fz * pa_z * pb_yz

                     + 10.0 * pa_yyz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_zz_s_0(double fx,
                                  double pa_yy,
                                  double pa_yyz,
                                  double pa_z,
                                  double pb_z,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_z

                     + 0.5 * fx * fx * pb_z

                     + 0.5 * pa_yyz * fx

                     + pa_yy * fx * pb_z

                     + 0.5 * fx * pa_z * pb_zz

                     + pa_yyz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyz_zz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_yy,
                                  double pa_yyz,
                                  double pa_z,
                                  double pb_z,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_z * fz * fgb

                     - 0.5 * fz * fga * pa_z * fx

                     - fz * fga * fx * pb_z

                     - pa_yyz * fz * fgb

                     + 1.5 * fx * fx * fz * pa_z

                     + 3.0 * fx * fx * fz * pb_z

                     + 4.0 * pa_yyz * fz * fx

                     + 8.0 * pa_yy * fz * fx * pb_z

                     - fz * fga * pa_z * pb_zz

                     + 4.0 * fx * fz * pa_z * pb_zz

                     + 10.0 * pa_yyz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xx_s_0(double fx,
                                  double pa_y,
                                  double pa_yzz,
                                  double pb_xx,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_y * fx * fx

                     + 0.5 * pa_yzz * fx

                     + 0.5 * pa_y * fx * pb_xx

                     + pa_yzz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xx_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yzz,
                                  double pb_xx,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_y * fx * fz * fgb

                     - 0.5 * pa_y * fz * fga * fx

                     - pa_yzz * fz * fgb

                     + 1.5 * pa_y * fz * fx * fx

                     + 4.0 * pa_yzz * fz * fx

                     - pa_y * fz * fga * pb_xx

                     + 4.0 * pa_y * fz * fx * pb_xx

                     + 10.0 * pa_yzz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xy_s_0(double fx,
                                  double pa_y,
                                  double pa_yzz,
                                  double pa_zz,
                                  double pb_x,
                                  double pb_xy,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pb_x

                     + 0.5 * fx * pa_zz * pb_x

                     + 0.5 * pa_y * fx * pb_xy

                     + pa_yzz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xy_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yzz,
                                  double pa_zz,
                                  double pb_x,
                                  double pb_xy,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fz * fga * pb_x

                     + 1.5 * fx * fx * fz * pb_x

                     + 4.0 * fx * pa_zz * fz * pb_x

                     - pa_y * fz * fga * pb_xy

                     + 4.0 * pa_y * fz * fx * pb_xy

                     + 10.0 * pa_yzz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xz_s_0(double fx,
                                  double pa_y,
                                  double pa_yz,
                                  double pa_yzz,
                                  double pb_x,
                                  double pb_xz,
                                  double s_0_0)
    {
        return s_0_0 * (pa_yz * fx * pb_x

                     + 0.5 * pa_y * fx * pb_xz

                     + pa_yzz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_xz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yz,
                                  double pa_yzz,
                                  double pb_x,
                                  double pb_xz,
                                  double r_0_0)
    {
        return r_0_0 * (8.0 * pa_yz * fz * fx * pb_x

                     - pa_y * fz * fga * pb_xz

                     + 4.0 * pa_y * fz * fx * pb_xz

                     + 10.0 * pa_yzz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_yy_s_0(double fx,
                                  double pa_y,
                                  double pa_yzz,
                                  double pa_zz,
                                  double pb_y,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.25 * pa_y * fx * fx

                     + 0.5 * fx * fx * pb_y

                     + 0.5 * pa_yzz * fx

                     + fx * pa_zz * pb_y

                     + 0.5 * pa_y * fx * pb_yy

                     + pa_yzz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_yy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yzz,
                                  double pa_zz,
                                  double pb_y,
                                  double pb_yy,
                                  double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_y * fx * fz * fgb

                     - 0.5 * pa_y * fz * fga * fx

                     - fx * fz * fga * pb_y

                     - pa_yzz * fz * fgb

                     + 1.5 * pa_y * fz * fx * fx

                     + 3.0 * fx * fx * fz * pb_y

                     + 4.0 * pa_yzz * fz * fx

                     + 8.0 * fx * pa_zz * fz * pb_y

                     - pa_y * fz * fga * pb_yy

                     + 4.0 * pa_y * fz * fx * pb_yy

                     + 10.0 * pa_yzz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_yz_s_0(double fx,
                                  double pa_y,
                                  double pa_yz,
                                  double pa_yzz,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_y,
                                  double pb_yz,
                                  double pb_z,
                                  double s_0_0)
    {
        return s_0_0 * (0.5 * fx * fx * pa_z

                     + 0.25 * fx * fx * pb_z

                     + pa_yz * fx * pb_y

                     + 0.5 * fx * pa_zz * pb_z

                     + 0.5 * pa_y * fx * pb_yz

                     + pa_yzz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_yz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yz,
                                  double pa_yzz,
                                  double pa_z,
                                  double pa_zz,
                                  double pb_y,
                                  double pb_yz,
                                  double pb_z,
                                  double r_0_0)
    {
        return r_0_0 * (3.0 * fx * fx * pa_z * fz

                     - 0.5 * fx * fz * fga * pb_z

                     + 1.5 * fx * fx * fz * pb_z

                     + 8.0 * pa_yz * fz * fx * pb_y

                     + 4.0 * fx * pa_zz * fz * pb_z

                     - pa_y * fz * fga * pb_yz

                     + 4.0 * pa_y * fz * fx * pb_yz

                     + 10.0 * pa_yzz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_zz_s_0(double fx,
                                  double pa_y,
                                  double pa_yz,
                                  double pa_yzz,
                                  double pb_z,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx

                     + 0.5 * pa_yzz * fx

                     + 2.0 * pa_yz * fx * pb_z

                     + 0.5 * pa_y * fx * pb_zz

                     + pa_yzz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzz_zz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_y,
                                  double pa_yz,
                                  double pa_yzz,
                                  double pb_z,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (4.5 * pa_y * fx * fx * fz

                     - 0.5 * pa_y * fx * fz * fgb

                     - 0.5 * pa_y * fz * fga * fx

                     - pa_yzz * fz * fgb

                     + 4.0 * pa_yzz * fz * fx

                     + 16.0 * pa_yz * fz * fx * pb_z

                     - pa_y * fz * fga * pb_zz

                     + 4.0 * pa_y * fz * fx * pb_zz

                     + 10.0 * pa_yzz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xx_s_0(double fx,
                                  double pa_z,
                                  double pa_zzz,
                                  double pb_xx,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_z * fx * fx

                     + 0.5 * pa_zzz * fx

                     + 1.5 * pa_z * fx * pb_xx

                     + pa_zzz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xx_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_z,
                                  double pa_zzz,
                                  double pb_xx,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_z * fx * fz * fgb

                     - 1.5 * pa_z * fz * fga * fx

                     - pa_zzz * fz * fgb

                     + 4.5 * pa_z * fz * fx * fx

                     + 4.0 * pa_zzz * fz * fx

                     - 3.0 * pa_z * fz * fga * pb_xx

                     + 12.0 * pa_z * fz * fx * pb_xx

                     + 10.0 * pa_zzz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xy_s_0(double fx,
                                  double pa_z,
                                  double pa_zzz,
                                  double pb_xy,
                                  double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * fx * pb_xy

                     + pa_zzz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xy_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_z,
                                  double pa_zzz,
                                  double pb_xy,
                                  double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_z * fz * fga * pb_xy

                     + 12.0 * pa_z * fz * fx * pb_xy

                     + 10.0 * pa_zzz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xz_s_0(double fx,
                                  double pa_z,
                                  double pa_zz,
                                  double pa_zzz,
                                  double pb_x,
                                  double pb_xz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_x

                     + 1.5 * pa_zz * fx * pb_x

                     + 1.5 * pa_z * fx * pb_xz

                     + pa_zzz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_xz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_z,
                                  double pa_zz,
                                  double pa_zzz,
                                  double pb_x,
                                  double pb_xz,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fz * fga * pb_x

                     + 4.5 * fx * fx * fz * pb_x

                     + 12.0 * pa_zz * fz * fx * pb_x

                     - 3.0 * pa_z * fz * fga * pb_xz

                     + 12.0 * pa_z * fz * fx * pb_xz

                     + 10.0 * pa_zzz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_yy_s_0(double fx,
                                  double pa_z,
                                  double pa_zzz,
                                  double pb_yy,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * pa_z * fx * fx

                     + 0.5 * pa_zzz * fx

                     + 1.5 * pa_z * fx * pb_yy

                     + pa_zzz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_yy_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_z,
                                  double pa_zzz,
                                  double pb_yy,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_z * fx * fz * fgb

                     - 1.5 * pa_z * fz * fga * fx

                     - pa_zzz * fz * fgb

                     + 4.5 * pa_z * fz * fx * fx

                     + 4.0 * pa_zzz * fz * fx

                     - 3.0 * pa_z * fz * fga * pb_yy

                     + 12.0 * pa_z * fz * fx * pb_yy

                     + 10.0 * pa_zzz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_yz_s_0(double fx,
                                  double pa_z,
                                  double pa_zz,
                                  double pa_zzz,
                                  double pb_y,
                                  double pb_yz,
                                  double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_y

                     + 1.5 * pa_zz * fx * pb_y

                     + 1.5 * pa_z * fx * pb_yz

                     + pa_zzz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_yz_r_0(double fga,
                                  double fx,
                                  double fz,
                                  double pa_z,
                                  double pa_zz,
                                  double pa_zzz,
                                  double pb_y,
                                  double pb_yz,
                                  double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fz * fga * pb_y

                     + 4.5 * fx * fx * fz * pb_y

                     + 12.0 * pa_zz * fz * fx * pb_y

                     - 3.0 * pa_z * fz * fga * pb_yz

                     + 12.0 * pa_z * fz * fx * pb_yz

                     + 10.0 * pa_zzz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_zz_s_0(double fx,
                                  double pa_z,
                                  double pa_zz,
                                  double pa_zzz,
                                  double pb_z,
                                  double pb_zz,
                                  double s_0_0)
    {
        return s_0_0 * (2.25 * pa_z * fx * fx

                     + 1.5 * fx * fx * pb_z

                     + 0.5 * pa_zzz * fx

                     + 3.0 * pa_zz * fx * pb_z

                     + 1.5 * pa_z * fx * pb_zz

                     + pa_zzz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzz_zz_r_0(double fga,
                                  double fgb,
                                  double fx,
                                  double fz,
                                  double pa_z,
                                  double pa_zz,
                                  double pa_zzz,
                                  double pb_z,
                                  double pb_zz,
                                  double r_0_0)
    {
        return r_0_0 * (13.5 * pa_z * fx * fx * fz

                     - 1.5 * pa_z * fx * fz * fgb

                     - 1.5 * pa_z * fz * fga * fx

                     - 3.0 * fx * fz * fga * pb_z

                     - pa_zzz * fz * fgb

                     + 9.0 * fx * fx * fz * pb_z

                     + 4.0 * pa_zzz * fz * fx

                     + 24.0 * pa_zz * fz * fx * pb_z

                     - 3.0 * pa_z * fz * fga * pb_zz

                     + 12.0 * pa_z * fz * fx * pb_zz

                     + 10.0 * pa_zzz * fz * pb_zz);

    }

    // SIMD elementary functions for (D|T|G) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxxx_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double pb_xxxx,
                                   double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx

                     + 0.75 * pa_xx * fx * fx

                     + 6.0 * pa_x * fx * fx * pb_x

                     + 4.5 * fx * fx * pb_xx

                     + 3.0 * pa_xx * pb_xx * fx

                     + 4.0 * pa_x * fx * pb_xxx

                     + 0.5 * fx * pb_xxxx

                     + pa_xx * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxxx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double pb_xxxx,
                                   double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * fz * fgb

                     + 11.25 * fx * fx * fx * fz

                     - 0.75 * fz * fga * fx * fx

                     - 3.0 * pa_xx * fx * fz * fgb

                     - 12.0 * pa_x * fx * pb_x * fz * fgb

                     + 6.0 * pa_xx * fz * fx * fx

                     + 48.0 * pa_x * fz * fx * fx * pb_x

                     + 36.0 * fx * fx * fz * pb_xx

                     - 3.0 * fx * pb_xx * fz * fgb

                     - 3.0 * fz * fga * pb_xx * fx

                     - 6.0 * pa_xx * pb_xx * fz * fgb

                     + 30.0 * pa_xx * fz * pb_xx * fx

                     + 40.0 * pa_x * fz * fx * pb_xxx

                     - fz * fga * pb_xxxx

                     + 5.0 * fz * fx * pb_xxxx

                     + 12.0 * pa_xx * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxxy_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_xxxy,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * fx * fx * pb_y

                     + 2.25 * fx * fx * pb_xy

                     + 1.5 * pa_xx * pb_xy * fx

                     + 3.0 * pa_x * fx * pb_xxy

                     + 0.5 * fx * pb_xxxy

                     + pa_xx * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxxy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_xxxy,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_x * fx * fz * fgb * pb_y

                     + 12.0 * pa_x * fz * fx * fx * pb_y

                     + 18.0 * fx * fx * fz * pb_xy

                     - 1.5 * fx * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pb_xy * fx

                     - 3.0 * pa_xx * pb_xy * fz * fgb

                     + 15.0 * pa_xx * fz * pb_xy * fx

                     + 30.0 * pa_x * fz * fx * pb_xxy

                     - fz * fga * pb_xxxy

                     + 5.0 * fz * fx * pb_xxxy

                     + 12.0 * pa_xx * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxxz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_xxxz,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * fx * fx * pb_z

                     + 2.25 * fx * fx * pb_xz

                     + 1.5 * pa_xx * pb_xz * fx

                     + 3.0 * pa_x * fx * pb_xxz

                     + 0.5 * fx * pb_xxxz

                     + pa_xx * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxxz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_xxxz,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_x * fx * fz * fgb * pb_z

                     + 12.0 * pa_x * fz * fx * fx * pb_z

                     + 18.0 * fx * fx * fz * pb_xz

                     - 1.5 * fx * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pb_xz * fx

                     - 3.0 * pa_xx * pb_xz * fz * fgb

                     + 15.0 * pa_xx * fz * pb_xz * fx

                     + 30.0 * pa_x * fz * fx * pb_xxz

                     - fz * fga * pb_xxxz

                     + 5.0 * fz * fx * pb_xxxz

                     + 12.0 * pa_xx * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxyy_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxyy,
                                   double pb_xyy,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.25 * pa_xx * fx * fx

                     + pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pb_yy

                     + 0.25 * fx * fx * pb_xx

                     + 0.5 * pa_xx * pb_xx * fx

                     + 0.5 * pa_xx * fx * pb_yy

                     + 2.0 * pa_x * fx * pb_xyy

                     + 0.5 * fx * pb_xxyy

                     + pa_xx * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxyy,
                                   double pb_xyy,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- fx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     - 0.25 * fz * fga * fx * fx

                     - pa_xx * fx * fz * fgb

                     - 2.0 * pa_x * fx * pb_x * fz * fgb

                     + 2.0 * pa_xx * fz * fx * fx

                     + 8.0 * pa_x * fz * fx * fx * pb_x

                     + 6.0 * fx * fx * fz * pb_yy

                     - 0.5 * fx * pb_xx * fz * fgb

                     - 0.5 * fx * fz * fgb * pb_yy

                     - 0.5 * fz * fga * pb_xx * fx

                     - 0.5 * fz * fga * fx * pb_yy

                     - pa_xx * pb_xx * fz * fgb

                     - pa_xx * fz * fgb * pb_yy

                     + 2.0 * fz * fx * fx * pb_xx

                     + 5.0 * pa_xx * fz * pb_xx * fx

                     + 5.0 * pa_xx * fz * fx * pb_yy

                     + 20.0 * pa_x * fz * fx * pb_xyy

                     - fz * fga * pb_xxyy

                     + 5.0 * fz * fx * pb_xxyy

                     + 12.0 * pa_xx * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_xxyz,
                                   double pb_xyz,
                                   double pb_yz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_yz

                     + 0.5 * pa_xx * fx * pb_yz

                     + 2.0 * pa_x * fx * pb_xyz

                     + 0.5 * fx * pb_xxyz

                     + pa_xx * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_xxyz,
                                   double pb_xyz,
                                   double pb_yz,
                                   double r_0_0)
    {
        return r_0_0 * (6.0 * fx * fx * fz * pb_yz

                     - 0.5 * fx * fz * fgb * pb_yz

                     - 0.5 * fz * fga * fx * pb_yz

                     - pa_xx * fz * fgb * pb_yz

                     + 5.0 * pa_xx * fz * fx * pb_yz

                     + 20.0 * pa_x * fz * fx * pb_xyz

                     - fz * fga * pb_xxyz

                     + 5.0 * fz * fx * pb_xxyz

                     + 12.0 * pa_xx * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxzz,
                                   double pb_xzz,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.25 * pa_xx * fx * fx

                     + pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pb_zz

                     + 0.25 * fx * fx * pb_xx

                     + 0.5 * pa_xx * pb_xx * fx

                     + 0.5 * pa_xx * fx * pb_zz

                     + 2.0 * pa_x * fx * pb_xzz

                     + 0.5 * fx * pb_xxzz

                     + pa_xx * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xxzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxzz,
                                   double pb_xzz,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- fx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     - 0.25 * fz * fga * fx * fx

                     - pa_xx * fx * fz * fgb

                     - 2.0 * pa_x * fx * pb_x * fz * fgb

                     + 2.0 * pa_xx * fz * fx * fx

                     + 8.0 * pa_x * fz * fx * fx * pb_x

                     + 6.0 * fx * fx * fz * pb_zz

                     - 0.5 * fx * pb_xx * fz * fgb

                     - 0.5 * fx * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pb_xx * fx

                     - 0.5 * fz * fga * fx * pb_zz

                     - pa_xx * pb_xx * fz * fgb

                     - pa_xx * fz * fgb * pb_zz

                     + 2.0 * fz * fx * fx * pb_xx

                     + 5.0 * pa_xx * fz * pb_xx * fx

                     + 5.0 * pa_xx * fz * fx * pb_zz

                     + 20.0 * pa_x * fz * fx * pb_xzz

                     - fz * fga * pb_xxzz

                     + 5.0 * fz * fx * pb_xxzz

                     + 12.0 * pa_xx * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xyyy_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_xy,
                                   double pb_xyyy,
                                   double pb_y,
                                   double pb_yyy,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * fx * fx * pb_y

                     + 0.75 * fx * fx * pb_xy

                     + 1.5 * pa_xx * pb_xy * fx

                     + pa_x * fx * pb_yyy

                     + 0.5 * fx * pb_xyyy

                     + pa_xx * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xyyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_xy,
                                   double pb_xyyy,
                                   double pb_y,
                                   double pb_yyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_x * fx * pb_y * fz * fgb

                     + 12.0 * pa_x * fz * fx * fx * pb_y

                     - 1.5 * fx * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pb_xy * fx

                     - 3.0 * pa_xx * pb_xy * fz * fgb

                     + 6.0 * fz * fx * fx * pb_xy

                     + 15.0 * pa_xx * fz * pb_xy * fx

                     + 10.0 * pa_x * fz * fx * pb_yyy

                     - fz * fga * pb_xyyy

                     + 5.0 * fz * fx * pb_xyyy

                     + 12.0 * pa_xx * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xyyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_xyyz,
                                   double pb_xz,
                                   double pb_yyz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * fx * pb_z

                     + 0.25 * fx * fx * pb_xz

                     + 0.5 * pa_xx * pb_xz * fx

                     + pa_x * fx * pb_yyz

                     + 0.5 * fx * pb_xyyz

                     + pa_xx * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xyyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_xyyz,
                                   double pb_xz,
                                   double pb_yyz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_x * fx * fz * fgb * pb_z

                     + 4.0 * pa_x * fz * fx * fx * pb_z

                     - 0.5 * fx * pb_xz * fz * fgb

                     - 0.5 * fz * fga * pb_xz * fx

                     - pa_xx * pb_xz * fz * fgb

                     + 2.0 * fz * fx * fx * pb_xz

                     + 5.0 * pa_xx * fz * pb_xz * fx

                     + 10.0 * pa_x * fz * fx * pb_yyz

                     - fz * fga * pb_xyyz

                     + 5.0 * fz * fx * pb_xyyz

                     + 12.0 * pa_xx * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xyzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_xy,
                                   double pb_xyzz,
                                   double pb_y,
                                   double pb_yzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * fx * pb_y

                     + 0.25 * fx * fx * pb_xy

                     + 0.5 * pa_xx * pb_xy * fx

                     + pa_x * fx * pb_yzz

                     + 0.5 * fx * pb_xyzz

                     + pa_xx * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xyzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_xy,
                                   double pb_xyzz,
                                   double pb_y,
                                   double pb_yzz,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_x * fx * pb_y * fz * fgb

                     + 4.0 * pa_x * fz * fx * fx * pb_y

                     - 0.5 * fx * pb_xy * fz * fgb

                     - 0.5 * fz * fga * pb_xy * fx

                     - pa_xx * pb_xy * fz * fgb

                     + 2.0 * fz * fx * fx * pb_xy

                     + 5.0 * pa_xx * fz * pb_xy * fx

                     + 10.0 * pa_x * fz * fx * pb_yzz

                     - fz * fga * pb_xyzz

                     + 5.0 * fz * fx * pb_xyzz

                     + 12.0 * pa_xx * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xzzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_xz,
                                   double pb_xzzz,
                                   double pb_z,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_x * fx * fx * pb_z

                     + 0.75 * fx * fx * pb_xz

                     + 1.5 * pa_xx * pb_xz * fx

                     + pa_x * fx * pb_zzz

                     + 0.5 * fx * pb_xzzz

                     + pa_xx * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_xzzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pb_xz,
                                   double pb_xzzz,
                                   double pb_z,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_x * fx * pb_z * fz * fgb

                     + 12.0 * pa_x * fz * fx * fx * pb_z

                     - 1.5 * fx * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pb_xz * fx

                     - 3.0 * pa_xx * pb_xz * fz * fgb

                     + 6.0 * fz * fx * fx * pb_xz

                     + 15.0 * pa_xx * fz * pb_xz * fx

                     + 10.0 * pa_x * fz * fx * pb_zzz

                     - fz * fga * pb_xzzz

                     + 5.0 * fz * fx * pb_xzzz

                     + 12.0 * pa_xx * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yyyy_s_0(double fx,
                                   double pa_xx,
                                   double pb_yy,
                                   double pb_yyyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_xx * fx * fx

                     + 1.5 * fx * fx * pb_yy

                     + 3.0 * pa_xx * pb_yy * fx

                     + 0.5 * fx * pb_yyyy

                     + pa_xx * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yyyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pb_yy,
                                   double pb_yyyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * fz * fgb

                     - 0.75 * fz * fga * fx * fx

                     - 3.0 * pa_xx * fx * fz * fgb

                     + 2.25 * fz * fx * fx * fx

                     + 6.0 * pa_xx * fz * fx * fx

                     - 3.0 * fx * pb_yy * fz * fgb

                     - 3.0 * fz * fga * pb_yy * fx

                     - 6.0 * pa_xx * pb_yy * fz * fgb

                     + 12.0 * fz * fx * fx * pb_yy

                     + 30.0 * pa_xx * fz * pb_yy * fx

                     - fz * fga * pb_yyyy

                     + 5.0 * fz * fx * pb_yyyy

                     + 12.0 * pa_xx * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yyyz_s_0(double fx,
                                   double pa_xx,
                                   double pb_yyyz,
                                   double pb_yz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_yz

                     + 1.5 * pa_xx * pb_yz * fx

                     + 0.5 * fx * pb_yyyz

                     + pa_xx * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yyyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pb_yyyz,
                                   double pb_yz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pb_yz * fx

                     - 3.0 * pa_xx * pb_yz * fz * fgb

                     + 6.0 * fz * fx * fx * pb_yz

                     + 15.0 * pa_xx * fz * pb_yz * fx

                     - fz * fga * pb_yyyz

                     + 5.0 * fz * fx * pb_yyyz

                     + 12.0 * pa_xx * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yyzz_s_0(double fx,
                                   double pa_xx,
                                   double pb_yy,
                                   double pb_yyzz,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * pa_xx * fx * fx

                     + 0.25 * fx * fx * pb_yy

                     + 0.25 * fx * fx * pb_zz

                     + 0.5 * pa_xx * pb_yy * fx

                     + 0.5 * pa_xx * fx * pb_zz

                     + 0.5 * fx * pb_yyzz

                     + pa_xx * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yyzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pb_yy,
                                   double pb_yyzz,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * fz * fgb

                     - 0.25 * fz * fga * fx * fx

                     - pa_xx * fx * fz * fgb

                     + 0.75 * fz * fx * fx * fx

                     + 2.0 * pa_xx * fz * fx * fx

                     - 0.5 * fx * pb_yy * fz * fgb

                     - 0.5 * fx * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pb_yy * fx

                     - 0.5 * fz * fga * fx * pb_zz

                     - pa_xx * pb_yy * fz * fgb

                     - pa_xx * fz * fgb * pb_zz

                     + 2.0 * fz * fx * fx * pb_yy

                     + 2.0 * fz * fx * fx * pb_zz

                     + 5.0 * pa_xx * fz * pb_yy * fx

                     + 5.0 * pa_xx * fz * fx * pb_zz

                     - fz * fga * pb_yyzz

                     + 5.0 * fz * fx * pb_yyzz

                     + 12.0 * pa_xx * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yzzz_s_0(double fx,
                                   double pa_xx,
                                   double pb_yz,
                                   double pb_yzzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_yz

                     + 1.5 * pa_xx * pb_yz * fx

                     + 0.5 * fx * pb_yzzz

                     + pa_xx * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_yzzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pb_yz,
                                   double pb_yzzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pb_yz * fx

                     - 3.0 * pa_xx * pb_yz * fz * fgb

                     + 6.0 * fz * fx * fx * pb_yz

                     + 15.0 * pa_xx * fz * pb_yz * fx

                     - fz * fga * pb_yzzz

                     + 5.0 * fz * fx * pb_yzzz

                     + 12.0 * pa_xx * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_zzzz_s_0(double fx,
                                   double pa_xx,
                                   double pb_zz,
                                   double pb_zzzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_xx * fx * fx

                     + 1.5 * fx * fx * pb_zz

                     + 3.0 * pa_xx * pb_zz * fx

                     + 0.5 * fx * pb_zzzz

                     + pa_xx * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xx_zzzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pb_zz,
                                   double pb_zzzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * fz * fgb

                     - 0.75 * fz * fga * fx * fx

                     - 3.0 * pa_xx * fx * fz * fgb

                     + 2.25 * fz * fx * fx * fx

                     + 6.0 * pa_xx * fz * fx * fx

                     - 3.0 * fx * pb_zz * fz * fgb

                     - 3.0 * fz * fga * pb_zz * fx

                     - 6.0 * pa_xx * pb_zz * fz * fgb

                     + 12.0 * fz * fx * fx * pb_zz

                     + 30.0 * pa_xx * fz * pb_zz * fx

                     - fz * fga * pb_zzzz

                     + 5.0 * fz * fx * pb_zzzz

                     + 12.0 * pa_xx * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxxx_s_0(double fx,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double pb_xxxx,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xy * fx * fx

                     + 3.0 * fx * fx * pa_y * pb_x

                     + 3.0 * pa_xy * pb_xx * fx

                     + 2.0 * fx * pa_y * pb_xxx

                     + pa_xy * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxxx_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double pb_xxxx,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fx * fz * fgb

                     - 6.0 * fx * pa_y * pb_x * fz * fgb

                     + 6.0 * pa_xy * fz * fx * fx

                     + 24.0 * fx * fx * fz * pa_y * pb_x

                     - 6.0 * pa_xy * pb_xx * fz * fgb

                     + 30.0 * pa_xy * fz * pb_xx * fx

                     + 20.0 * fx * fz * pa_y * pb_xxx

                     + 12.0 * pa_xy * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxxy_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double pb_xxxy,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_y * pb_y

                     + 0.75 * fx * fx * pb_xx

                     + 1.5 * pa_xy * pb_xy * fx

                     + 0.5 * pa_x * fx * pb_xxx

                     + 1.5 * fx * pa_y * pb_xxy

                     + pa_xy * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxxy_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double pb_xxxy,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     - 1.5 * pa_x * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_y * fz * fgb * pb_y

                     + 6.0 * pa_x * fz * fx * fx * pb_x

                     + 6.0 * fx * fx * fz * pa_y * pb_y

                     + 6.0 * fx * fx * fz * pb_xx

                     - 3.0 * pa_xy * pb_xy * fz * fgb

                     + 15.0 * pa_xy * fz * pb_xy * fx

                     + 5.0 * pa_x * fz * fx * pb_xxx

                     + 15.0 * fx * fz * pa_y * pb_xxy

                     + 12.0 * pa_xy * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxxz_s_0(double fx,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_xxxz,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_y * pb_z

                     + 1.5 * pa_xy * pb_xz * fx

                     + 1.5 * fx * pa_y * pb_xxz

                     + pa_xy * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxxz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_xxxz,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_y * fz * fgb * pb_z

                     + 6.0 * fx * fx * fz * pa_y * pb_z

                     - 3.0 * pa_xy * pb_xz * fz * fgb

                     + 15.0 * pa_xy * fz * pb_xz * fx

                     + 15.0 * fx * fz * pa_y * pb_xxz

                     + 12.0 * pa_xy * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxyy_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_xxyy,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xy * fx * fx

                     + 0.5 * pa_x * fx * fx * pb_y

                     + 0.5 * fx * fx * pa_y * pb_x

                     + fx * fx * pb_xy

                     + 0.5 * pa_xy * pb_xx * fx

                     + 0.5 * pa_xy * fx * pb_yy

                     + pa_x * fx * pb_xxy

                     + fx * pa_y * pb_xyy

                     + pa_xy * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxyy_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_xxyy,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_xy * fx * fz * fgb

                     - pa_x * fx * fz * fgb * pb_y

                     - fx * pa_y * pb_x * fz * fgb

                     + 2.0 * pa_xy * fz * fx * fx

                     + 4.0 * pa_x * fz * fx * fx * pb_y

                     + 4.0 * fx * fx * fz * pa_y * pb_x

                     + 8.0 * fx * fx * fz * pb_xy

                     - pa_xy * pb_xx * fz * fgb

                     - pa_xy * fz * fgb * pb_yy

                     + 5.0 * pa_xy * fz * pb_xx * fx

                     + 5.0 * pa_xy * fz * fx * pb_yy

                     + 10.0 * pa_x * fz * fx * pb_xxy

                     + 10.0 * fx * fz * pa_y * pb_xyy

                     + 12.0 * pa_xy * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_xxyz,
                                   double pb_xxz,
                                   double pb_xyz,
                                   double pb_xz,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx * pb_z

                     + 0.5 * fx * fx * pb_xz

                     + 0.5 * pa_xy * fx * pb_yz

                     + 0.5 * pa_x * fx * pb_xxz

                     + fx * pa_y * pb_xyz

                     + pa_xy * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxyz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_xxyz,
                                   double pb_xxz,
                                   double pb_xyz,
                                   double pb_xz,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fz * fgb * pb_z

                     + 2.0 * pa_x * fz * fx * fx * pb_z

                     + 4.0 * fx * fx * fz * pb_xz

                     - pa_xy * fz * fgb * pb_yz

                     + 5.0 * pa_xy * fz * fx * pb_yz

                     + 5.0 * pa_x * fz * fx * pb_xxz

                     + 10.0 * fx * fz * pa_y * pb_xyz

                     + 12.0 * pa_xy * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxzz_s_0(double fx,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxzz,
                                   double pb_xzz,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xy * fx * fx

                     + 0.5 * fx * fx * pa_y * pb_x

                     + 0.5 * pa_xy * pb_xx * fx

                     + 0.5 * pa_xy * fx * pb_zz

                     + fx * pa_y * pb_xzz

                     + pa_xy * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xxzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxzz,
                                   double pb_xzz,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_xy * fx * fz * fgb

                     - fx * pa_y * pb_x * fz * fgb

                     + 2.0 * pa_xy * fz * fx * fx

                     + 4.0 * fx * fx * fz * pa_y * pb_x

                     - pa_xy * pb_xx * fz * fgb

                     - pa_xy * fz * fgb * pb_zz

                     + 5.0 * pa_xy * fz * pb_xx * fx

                     + 5.0 * pa_xy * fz * fx * pb_zz

                     + 10.0 * fx * fz * pa_y * pb_xzz

                     + 12.0 * pa_xy * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xyyy_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_xyyy,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_y * pb_y

                     + 0.75 * fx * fx * pb_yy

                     + 1.5 * pa_xy * pb_xy * fx

                     + 1.5 * pa_x * fx * pb_xyy

                     + 0.5 * fx * pa_y * pb_yyy

                     + pa_xy * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xyyy_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_xyyy,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     - 1.5 * pa_x * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_y * pb_y * fz * fgb

                     + 6.0 * pa_x * fz * fx * fx * pb_x

                     + 6.0 * fx * fx * fz * pa_y * pb_y

                     + 6.0 * fx * fx * fz * pb_yy

                     - 3.0 * pa_xy * pb_xy * fz * fgb

                     + 15.0 * pa_xy * fz * pb_xy * fx

                     + 15.0 * pa_x * fz * fx * pb_xyy

                     + 5.0 * fx * fz * pa_y * pb_yyy

                     + 12.0 * pa_xy * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xyyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_xyyz,
                                   double pb_xyz,
                                   double pb_xz,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_y * pb_z

                     + 0.5 * fx * fx * pb_yz

                     + 0.5 * pa_xy * pb_xz * fx

                     + pa_x * fx * pb_xyz

                     + 0.5 * fx * pa_y * pb_yyz

                     + pa_xy * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xyyz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_xyyz,
                                   double pb_xyz,
                                   double pb_xz,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_y * fz * fgb * pb_z

                     + 2.0 * fx * fx * fz * pa_y * pb_z

                     + 4.0 * fx * fx * fz * pb_yz

                     - pa_xy * pb_xz * fz * fgb

                     + 5.0 * pa_xy * fz * pb_xz * fx

                     + 10.0 * pa_x * fz * fx * pb_xyz

                     + 5.0 * fx * fz * pa_y * pb_yyz

                     + 12.0 * pa_xy * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xyzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyzz,
                                   double pb_xzz,
                                   double pb_y,
                                   double pb_yzz,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * pa_x * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_y * pb_y

                     + 0.25 * fx * fx * pb_zz

                     + 0.5 * pa_xy * pb_xy * fx

                     + 0.5 * pa_x * fx * pb_xzz

                     + 0.5 * fx * pa_y * pb_yzz

                     + pa_xy * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xyzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyzz,
                                   double pb_xzz,
                                   double pb_y,
                                   double pb_yzz,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * fz * fgb

                     + 0.75 * fx * fx * fx * fz

                     - 0.5 * pa_x * fx * pb_x * fz * fgb

                     - 0.5 * fx * pa_y * pb_y * fz * fgb

                     + 2.0 * pa_x * fz * fx * fx * pb_x

                     + 2.0 * fx * fx * fz * pa_y * pb_y

                     + 2.0 * fx * fx * fz * pb_zz

                     - pa_xy * pb_xy * fz * fgb

                     + 5.0 * pa_xy * fz * pb_xy * fx

                     + 5.0 * pa_x * fz * fx * pb_xzz

                     + 5.0 * fx * fz * pa_y * pb_yzz

                     + 12.0 * pa_xy * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xzzz_s_0(double fx,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_xz,
                                   double pb_xzzz,
                                   double pb_z,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_y * pb_z

                     + 1.5 * pa_xy * pb_xz * fx

                     + 0.5 * fx * pa_y * pb_zzz

                     + pa_xy * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_xzzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_xz,
                                   double pb_xzzz,
                                   double pb_z,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_y * pb_z * fz * fgb

                     + 6.0 * fx * fx * fz * pa_y * pb_z

                     - 3.0 * pa_xy * pb_xz * fz * fgb

                     + 15.0 * pa_xy * fz * pb_xz * fx

                     + 5.0 * fx * fz * pa_y * pb_zzz

                     + 12.0 * pa_xy * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yyyy_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double pb_yyyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xy * fx * fx

                     + 3.0 * pa_x * fx * fx * pb_y

                     + 3.0 * pa_xy * pb_yy * fx

                     + 2.0 * pa_x * fx * pb_yyy

                     + pa_xy * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yyyy_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double pb_yyyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fx * fz * fgb

                     - 6.0 * pa_x * fx * pb_y * fz * fgb

                     + 6.0 * pa_xy * fz * fx * fx

                     + 24.0 * pa_x * fz * fx * fx * pb_y

                     - 6.0 * pa_xy * pb_yy * fz * fgb

                     + 30.0 * pa_xy * fz * pb_yy * fx

                     + 20.0 * pa_x * fz * fx * pb_yyy

                     + 12.0 * pa_xy * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yyyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pb_yyyz,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * pb_z

                     + 1.5 * pa_xy * pb_yz * fx

                     + 1.5 * pa_x * fx * pb_yyz

                     + pa_xy * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yyyz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pb_yyyz,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fz * fgb * pb_z

                     + 6.0 * pa_x * fz * fx * fx * pb_z

                     - 3.0 * pa_xy * pb_yz * fz * fgb

                     + 15.0 * pa_xy * fz * pb_yz * fx

                     + 15.0 * pa_x * fz * fx * pb_yyz

                     + 12.0 * pa_xy * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yyzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyzz,
                                   double pb_yzz,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xy * fx * fx

                     + 0.5 * pa_x * fx * fx * pb_y

                     + 0.5 * pa_xy * pb_yy * fx

                     + 0.5 * pa_xy * fx * pb_zz

                     + pa_x * fx * pb_yzz

                     + pa_xy * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yyzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyzz,
                                   double pb_yzz,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_xy * fx * fz * fgb

                     - pa_x * fx * pb_y * fz * fgb

                     + 2.0 * pa_xy * fz * fx * fx

                     + 4.0 * pa_x * fz * fx * fx * pb_y

                     - pa_xy * pb_yy * fz * fgb

                     - pa_xy * fz * fgb * pb_zz

                     + 5.0 * pa_xy * fz * pb_yy * fx

                     + 5.0 * pa_xy * fz * fx * pb_zz

                     + 10.0 * pa_x * fz * fx * pb_yzz

                     + 12.0 * pa_xy * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yzzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pb_yz,
                                   double pb_yzzz,
                                   double pb_z,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * pb_z

                     + 1.5 * pa_xy * pb_yz * fx

                     + 0.5 * pa_x * fx * pb_zzz

                     + pa_xy * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_yzzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pb_yz,
                                   double pb_yzzz,
                                   double pb_z,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * pb_z * fz * fgb

                     + 6.0 * pa_x * fz * fx * fx * pb_z

                     - 3.0 * pa_xy * pb_yz * fz * fgb

                     + 15.0 * pa_xy * fz * pb_yz * fx

                     + 5.0 * pa_x * fz * fx * pb_zzz

                     + 12.0 * pa_xy * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_zzzz_s_0(double fx,
                                   double pa_xy,
                                   double pb_zz,
                                   double pb_zzzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xy * fx * fx

                     + 3.0 * pa_xy * pb_zz * fx

                     + pa_xy * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xy_zzzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xy,
                                   double pb_zz,
                                   double pb_zzzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xy * fx * fz * fgb

                     + 6.0 * pa_xy * fz * fx * fx

                     - 6.0 * pa_xy * pb_zz * fz * fgb

                     + 30.0 * pa_xy * fz * pb_zz * fx

                     + 12.0 * pa_xy * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxxx_s_0(double fx,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double pb_xxxx,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xz * fx * fx

                     + 3.0 * fx * fx * pa_z * pb_x

                     + 3.0 * pa_xz * pb_xx * fx

                     + 2.0 * fx * pa_z * pb_xxx

                     + pa_xz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxxx_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double pb_xxxx,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fx * fz * fgb

                     - 6.0 * fx * pa_z * pb_x * fz * fgb

                     + 6.0 * pa_xz * fz * fx * fx

                     + 24.0 * fx * fx * fz * pa_z * pb_x

                     - 6.0 * pa_xz * pb_xx * fz * fgb

                     + 30.0 * pa_xz * fz * pb_xx * fx

                     + 20.0 * fx * fz * pa_z * pb_xxx

                     + 12.0 * pa_xz * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxxy_s_0(double fx,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_xxxy,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z * pb_y

                     + 1.5 * pa_xz * pb_xy * fx

                     + 1.5 * fx * pa_z * pb_xxy

                     + pa_xz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxxy_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_xxxy,
                                   double pb_xxy,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_z * fz * fgb * pb_y

                     + 6.0 * fx * fx * fz * pa_z * pb_y

                     - 3.0 * pa_xz * pb_xy * fz * fgb

                     + 15.0 * pa_xz * fz * pb_xy * fx

                     + 15.0 * fx * fz * pa_z * pb_xxy

                     + 12.0 * pa_xz * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxxz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double pb_xxxz,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_z * pb_z

                     + 0.75 * fx * fx * pb_xx

                     + 1.5 * pa_xz * pb_xz * fx

                     + 0.5 * pa_x * fx * pb_xxx

                     + 1.5 * fx * pa_z * pb_xxz

                     + pa_xz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxxz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxx,
                                   double pb_xxxz,
                                   double pb_xxz,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     - 1.5 * pa_x * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_z * fz * fgb * pb_z

                     + 6.0 * pa_x * fz * fx * fx * pb_x

                     + 6.0 * fx * fx * fz * pa_z * pb_z

                     + 6.0 * fx * fx * fz * pb_xx

                     - 3.0 * pa_xz * pb_xz * fz * fgb

                     + 15.0 * pa_xz * fz * pb_xz * fx

                     + 5.0 * pa_x * fz * fx * pb_xxx

                     + 15.0 * fx * fz * pa_z * pb_xxz

                     + 12.0 * pa_xz * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxyy_s_0(double fx,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxyy,
                                   double pb_xyy,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xz * fx * fx

                     + 0.5 * fx * fx * pa_z * pb_x

                     + 0.5 * pa_xz * pb_xx * fx

                     + 0.5 * pa_xz * fx * pb_yy

                     + fx * pa_z * pb_xyy

                     + pa_xz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxyy_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxyy,
                                   double pb_xyy,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_xz * fx * fz * fgb

                     - fx * pa_z * pb_x * fz * fgb

                     + 2.0 * pa_xz * fz * fx * fx

                     + 4.0 * fx * fx * fz * pa_z * pb_x

                     - pa_xz * pb_xx * fz * fgb

                     - pa_xz * fz * fgb * pb_yy

                     + 5.0 * pa_xz * fz * pb_xx * fx

                     + 5.0 * pa_xz * fz * fx * pb_yy

                     + 10.0 * fx * fz * pa_z * pb_xyy

                     + 12.0 * pa_xz * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_xxy,
                                   double pb_xxyz,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_y,
                                   double pb_yz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_x * fx * fx * pb_y

                     + 0.5 * fx * fx * pb_xy

                     + 0.5 * pa_xz * fx * pb_yz

                     + 0.5 * pa_x * fx * pb_xxy

                     + fx * pa_z * pb_xyz

                     + pa_xz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxyz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_xxy,
                                   double pb_xxyz,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_y,
                                   double pb_yz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_x * fx * fz * fgb * pb_y

                     + 2.0 * pa_x * fz * fx * fx * pb_y

                     + 4.0 * fx * fx * fz * pb_xy

                     - pa_xz * fz * fgb * pb_yz

                     + 5.0 * pa_xz * fz * fx * pb_yz

                     + 5.0 * pa_x * fz * fx * pb_xxy

                     + 10.0 * fx * fz * pa_z * pb_xyz

                     + 12.0 * pa_xz * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_xxzz,
                                   double pb_xz,
                                   double pb_xzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xz * fx * fx

                     + 0.5 * pa_x * fx * fx * pb_z

                     + 0.5 * fx * fx * pa_z * pb_x

                     + fx * fx * pb_xz

                     + 0.5 * pa_xz * pb_xx * fx

                     + 0.5 * pa_xz * fx * pb_zz

                     + pa_x * fx * pb_xxz

                     + fx * pa_z * pb_xzz

                     + pa_xz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xxzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_xxzz,
                                   double pb_xz,
                                   double pb_xzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_xz * fx * fz * fgb

                     - pa_x * fx * fz * fgb * pb_z

                     - fx * pa_z * pb_x * fz * fgb

                     + 2.0 * pa_xz * fz * fx * fx

                     + 4.0 * pa_x * fz * fx * fx * pb_z

                     + 4.0 * fx * fx * fz * pa_z * pb_x

                     + 8.0 * fx * fx * fz * pb_xz

                     - pa_xz * pb_xx * fz * fgb

                     - pa_xz * fz * fgb * pb_zz

                     + 5.0 * pa_xz * fz * pb_xx * fx

                     + 5.0 * pa_xz * fz * fx * pb_zz

                     + 10.0 * pa_x * fz * fx * pb_xxz

                     + 10.0 * fx * fz * pa_z * pb_xzz

                     + 12.0 * pa_xz * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xyyy_s_0(double fx,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_xy,
                                   double pb_xyyy,
                                   double pb_y,
                                   double pb_yyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z * pb_y

                     + 1.5 * pa_xz * pb_xy * fx

                     + 0.5 * fx * pa_z * pb_yyy

                     + pa_xz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xyyy_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_xy,
                                   double pb_xyyy,
                                   double pb_y,
                                   double pb_yyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_z * pb_y * fz * fgb

                     + 6.0 * fx * fx * fz * pa_z * pb_y

                     - 3.0 * pa_xz * pb_xy * fz * fgb

                     + 15.0 * pa_xz * fz * pb_xy * fx

                     + 5.0 * fx * fz * pa_z * pb_yyy

                     + 12.0 * pa_xz * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xyyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xyy,
                                   double pb_xyyz,
                                   double pb_xz,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * pa_x * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_z * pb_z

                     + 0.25 * fx * fx * pb_yy

                     + 0.5 * pa_xz * pb_xz * fx

                     + 0.5 * pa_x * fx * pb_xyy

                     + 0.5 * fx * pa_z * pb_yyz

                     + pa_xz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xyyz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xyy,
                                   double pb_xyyz,
                                   double pb_xz,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * fz * fgb

                     + 0.75 * fx * fx * fx * fz

                     - 0.5 * pa_x * fx * pb_x * fz * fgb

                     - 0.5 * fx * pa_z * fz * fgb * pb_z

                     + 2.0 * pa_x * fz * fx * fx * pb_x

                     + 2.0 * fx * fx * fz * pa_z * pb_z

                     + 2.0 * fx * fx * fz * pb_yy

                     - pa_xz * pb_xz * fz * fgb

                     + 5.0 * pa_xz * fz * pb_xz * fx

                     + 5.0 * pa_x * fz * fx * pb_xyy

                     + 5.0 * fx * fz * pa_z * pb_yyz

                     + 12.0 * pa_xz * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xyzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_xyzz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_z * pb_y

                     + 0.5 * fx * fx * pb_yz

                     + 0.5 * pa_xz * pb_xy * fx

                     + pa_x * fx * pb_xyz

                     + 0.5 * fx * pa_z * pb_yzz

                     + pa_xz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xyzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_xyzz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_z * pb_y * fz * fgb

                     + 2.0 * fx * fx * fz * pa_z * pb_y

                     + 4.0 * fx * fx * fz * pb_yz

                     - pa_xz * pb_xy * fz * fgb

                     + 5.0 * pa_xz * fz * pb_xy * fx

                     + 10.0 * pa_x * fz * fx * pb_xyz

                     + 5.0 * fx * fz * pa_z * pb_yzz

                     + 12.0 * pa_xz * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xzzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double pb_xzzz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_z * pb_z

                     + 0.75 * fx * fx * pb_zz

                     + 1.5 * pa_xz * pb_xz * fx

                     + 1.5 * pa_x * fx * pb_xzz

                     + 0.5 * fx * pa_z * pb_zzz

                     + pa_xz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_xzzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double pb_xzzz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     - 1.5 * pa_x * fx * pb_x * fz * fgb

                     - 1.5 * fx * pa_z * pb_z * fz * fgb

                     + 6.0 * pa_x * fz * fx * fx * pb_x

                     + 6.0 * fx * fx * fz * pa_z * pb_z

                     + 6.0 * fx * fx * fz * pb_zz

                     - 3.0 * pa_xz * pb_xz * fz * fgb

                     + 15.0 * pa_xz * fz * pb_xz * fx

                     + 15.0 * pa_x * fz * fx * pb_xzz

                     + 5.0 * fx * fz * pa_z * pb_zzz

                     + 12.0 * pa_xz * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yyyy_s_0(double fx,
                                   double pa_xz,
                                   double pb_yy,
                                   double pb_yyyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xz * fx * fx

                     + 3.0 * pa_xz * pb_yy * fx

                     + pa_xz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yyyy_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xz,
                                   double pb_yy,
                                   double pb_yyyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fx * fz * fgb

                     + 6.0 * pa_xz * fz * fx * fx

                     - 6.0 * pa_xz * pb_yy * fz * fgb

                     + 30.0 * pa_xz * fz * pb_yy * fx

                     + 12.0 * pa_xz * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yyyz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yyy,
                                   double pb_yyyz,
                                   double pb_yz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * pb_y

                     + 1.5 * pa_xz * pb_yz * fx

                     + 0.5 * pa_x * fx * pb_yyy

                     + pa_xz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yyyz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yyy,
                                   double pb_yyyz,
                                   double pb_yz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * pb_y * fz * fgb

                     + 6.0 * pa_x * fz * fx * fx * pb_y

                     - 3.0 * pa_xz * pb_yz * fz * fgb

                     + 15.0 * pa_xz * fz * pb_yz * fx

                     + 5.0 * pa_x * fz * fx * pb_yyy

                     + 12.0 * pa_xz * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yyzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_yyzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xz * fx * fx

                     + 0.5 * pa_x * fx * fx * pb_z

                     + 0.5 * pa_xz * pb_yy * fx

                     + 0.5 * pa_xz * fx * pb_zz

                     + pa_x * fx * pb_yyz

                     + pa_xz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yyzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_yyzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_xz * fx * fz * fgb

                     - pa_x * fx * fz * fgb * pb_z

                     + 2.0 * pa_xz * fz * fx * fx

                     + 4.0 * pa_x * fz * fx * fx * pb_z

                     - pa_xz * pb_yy * fz * fgb

                     - pa_xz * fz * fgb * pb_zz

                     + 5.0 * pa_xz * fz * pb_yy * fx

                     + 5.0 * pa_xz * fz * fx * pb_zz

                     + 10.0 * pa_x * fz * fx * pb_yyz

                     + 12.0 * pa_xz * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yzzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double pb_yzzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * pb_y

                     + 1.5 * pa_xz * pb_yz * fx

                     + 1.5 * pa_x * fx * pb_yzz

                     + pa_xz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_yzzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double pb_yzzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * pb_y * fz * fgb

                     + 6.0 * pa_x * fz * fx * fx * pb_y

                     - 3.0 * pa_xz * pb_yz * fz * fgb

                     + 15.0 * pa_xz * fz * pb_yz * fx

                     + 15.0 * pa_x * fz * fx * pb_yzz

                     + 12.0 * pa_xz * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_zzzz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double pb_zzzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xz * fx * fx

                     + 3.0 * pa_x * fx * fx * pb_z

                     + 3.0 * pa_xz * pb_zz * fx

                     + 2.0 * pa_x * fx * pb_zzz

                     + pa_xz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xz_zzzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double pb_zzzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_xz * fx * fz * fgb

                     - 6.0 * pa_x * fx * pb_z * fz * fgb

                     + 6.0 * pa_xz * fz * fx * fx

                     + 24.0 * pa_x * fz * fx * fx * pb_z

                     - 6.0 * pa_xz * pb_zz * fz * fgb

                     + 30.0 * pa_xz * fz * pb_zz * fx

                     + 20.0 * pa_x * fz * fx * pb_zzz

                     + 12.0 * pa_xz * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxxx_s_0(double fx,
                                   double pa_yy,
                                   double pb_xx,
                                   double pb_xxxx,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_yy * fx * fx

                     + 1.5 * fx * fx * pb_xx

                     + 3.0 * pa_yy * pb_xx * fx

                     + 0.5 * fx * pb_xxxx

                     + pa_yy * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxxx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yy,
                                   double pb_xx,
                                   double pb_xxxx,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * fz * fgb

                     - 0.75 * fz * fga * fx * fx

                     - 3.0 * pa_yy * fx * fz * fgb

                     + 2.25 * fz * fx * fx * fx

                     + 6.0 * pa_yy * fz * fx * fx

                     - 3.0 * fx * pb_xx * fz * fgb

                     - 3.0 * fz * fga * pb_xx * fx

                     - 6.0 * pa_yy * pb_xx * fz * fgb

                     + 12.0 * fz * fx * fx * pb_xx

                     + 30.0 * pa_yy * fz * pb_xx * fx

                     - fz * fga * pb_xxxx

                     + 5.0 * fz * fx * pb_xxxx

                     + 12.0 * pa_yy * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxxy_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xxx,
                                   double pb_xxxy,
                                   double pb_xy,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * fx * fx * pb_x

                     + 0.75 * fx * fx * pb_xy

                     + 1.5 * pa_yy * pb_xy * fx

                     + pa_y * fx * pb_xxx

                     + 0.5 * fx * pb_xxxy

                     + pa_yy * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxxy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xxx,
                                   double pb_xxxy,
                                   double pb_xy,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_y * fx * pb_x * fz * fgb

                     + 12.0 * pa_y * fz * fx * fx * pb_x

                     - 1.5 * fx * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pb_xy * fx

                     - 3.0 * pa_yy * pb_xy * fz * fgb

                     + 6.0 * fz * fx * fx * pb_xy

                     + 15.0 * pa_yy * fz * pb_xy * fx

                     + 10.0 * pa_y * fz * fx * pb_xxx

                     - fz * fga * pb_xxxy

                     + 5.0 * fz * fx * pb_xxxy

                     + 12.0 * pa_yy * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxxz_s_0(double fx,
                                   double pa_yy,
                                   double pb_xxxz,
                                   double pb_xz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_xz

                     + 1.5 * pa_yy * pb_xz * fx

                     + 0.5 * fx * pb_xxxz

                     + pa_yy * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxxz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yy,
                                   double pb_xxxz,
                                   double pb_xz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pb_xz * fx

                     - 3.0 * pa_yy * pb_xz * fz * fgb

                     + 6.0 * fz * fx * fx * pb_xz

                     + 15.0 * pa_yy * fz * pb_xz * fx

                     - fz * fga * pb_xxxz

                     + 5.0 * fz * fx * pb_xxxz

                     + 12.0 * pa_yy * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxyy_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_xxyy,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.25 * pa_yy * fx * fx

                     + pa_y * fx * fx * pb_y

                     + 0.75 * fx * fx * pb_xx

                     + 0.25 * fx * fx * pb_yy

                     + 0.5 * pa_yy * pb_xx * fx

                     + 0.5 * pa_yy * fx * pb_yy

                     + 2.0 * pa_y * fx * pb_xxy

                     + 0.5 * fx * pb_xxyy

                     + pa_yy * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_xxyy,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- fx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     - 0.25 * fz * fga * fx * fx

                     - pa_yy * fx * fz * fgb

                     - 2.0 * pa_y * fx * fz * fgb * pb_y

                     + 2.0 * pa_yy * fz * fx * fx

                     + 8.0 * pa_y * fz * fx * fx * pb_y

                     + 6.0 * fx * fx * fz * pb_xx

                     - 0.5 * fx * pb_xx * fz * fgb

                     - 0.5 * fx * fz * fgb * pb_yy

                     - 0.5 * fz * fga * pb_xx * fx

                     - 0.5 * fz * fga * fx * pb_yy

                     - pa_yy * pb_xx * fz * fgb

                     - pa_yy * fz * fgb * pb_yy

                     + 2.0 * fz * fx * fx * pb_yy

                     + 5.0 * pa_yy * fz * pb_xx * fx

                     + 5.0 * pa_yy * fz * fx * pb_yy

                     + 20.0 * pa_y * fz * fx * pb_xxy

                     - fz * fga * pb_xxyy

                     + 5.0 * fz * fx * pb_xxyy

                     + 12.0 * pa_yy * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxyz_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_xxyz,
                                   double pb_xxz,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx * fx * pb_z

                     + 0.25 * fx * fx * pb_yz

                     + 0.5 * pa_yy * fx * pb_yz

                     + pa_y * fx * pb_xxz

                     + 0.5 * fx * pb_xxyz

                     + pa_yy * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_xxyz,
                                   double pb_xxz,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_y * fx * fz * fgb * pb_z

                     + 4.0 * pa_y * fz * fx * fx * pb_z

                     - 0.5 * fx * fz * fgb * pb_yz

                     - 0.5 * fz * fga * fx * pb_yz

                     - pa_yy * fz * fgb * pb_yz

                     + 2.0 * fz * fx * fx * pb_yz

                     + 5.0 * pa_yy * fz * fx * pb_yz

                     + 10.0 * pa_y * fz * fx * pb_xxz

                     - fz * fga * pb_xxyz

                     + 5.0 * fz * fx * pb_xxyz

                     + 12.0 * pa_yy * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxzz_s_0(double fx,
                                   double pa_yy,
                                   double pb_xx,
                                   double pb_xxzz,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * pa_yy * fx * fx

                     + 0.25 * fx * fx * pb_xx

                     + 0.25 * fx * fx * pb_zz

                     + 0.5 * pa_yy * pb_xx * fx

                     + 0.5 * pa_yy * fx * pb_zz

                     + 0.5 * fx * pb_xxzz

                     + pa_yy * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xxzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yy,
                                   double pb_xx,
                                   double pb_xxzz,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * fz * fgb

                     - 0.25 * fz * fga * fx * fx

                     - pa_yy * fx * fz * fgb

                     + 0.75 * fz * fx * fx * fx

                     + 2.0 * pa_yy * fz * fx * fx

                     - 0.5 * fx * pb_xx * fz * fgb

                     - 0.5 * fx * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pb_xx * fx

                     - 0.5 * fz * fga * fx * pb_zz

                     - pa_yy * pb_xx * fz * fgb

                     - pa_yy * fz * fgb * pb_zz

                     + 2.0 * fz * fx * fx * pb_xx

                     + 2.0 * fz * fx * fx * pb_zz

                     + 5.0 * pa_yy * fz * pb_xx * fx

                     + 5.0 * pa_yy * fz * fx * pb_zz

                     - fz * fga * pb_xxzz

                     + 5.0 * fz * fx * pb_xxzz

                     + 12.0 * pa_yy * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xyyy_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_xyyy,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * fx * fx * pb_x

                     + 2.25 * fx * fx * pb_xy

                     + 1.5 * pa_yy * pb_xy * fx

                     + 3.0 * pa_y * fx * pb_xyy

                     + 0.5 * fx * pb_xyyy

                     + pa_yy * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xyyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_xyyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_y * fx * pb_x * fz * fgb

                     + 12.0 * pa_y * fz * fx * fx * pb_x

                     + 18.0 * fx * fx * fz * pb_xy

                     - 1.5 * fx * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pb_xy * fx

                     - 3.0 * pa_yy * pb_xy * fz * fgb

                     + 15.0 * pa_yy * fz * pb_xy * fx

                     + 30.0 * pa_y * fz * fx * pb_xyy

                     - fz * fga * pb_xyyy

                     + 5.0 * fz * fx * pb_xyyy

                     + 12.0 * pa_yy * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xyyz_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_xyyz,
                                   double pb_xyz,
                                   double pb_xz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_xz

                     + 0.5 * pa_yy * pb_xz * fx

                     + 2.0 * pa_y * fx * pb_xyz

                     + 0.5 * fx * pb_xyyz

                     + pa_yy * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xyyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_xyyz,
                                   double pb_xyz,
                                   double pb_xz,
                                   double r_0_0)
    {
        return r_0_0 * (6.0 * fx * fx * fz * pb_xz

                     - 0.5 * fx * pb_xz * fz * fgb

                     - 0.5 * fz * fga * pb_xz * fx

                     - pa_yy * pb_xz * fz * fgb

                     + 5.0 * pa_yy * fz * pb_xz * fx

                     + 20.0 * pa_y * fz * fx * pb_xyz

                     - fz * fga * pb_xyyz

                     + 5.0 * fz * fx * pb_xyyz

                     + 12.0 * pa_yy * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xyzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyzz,
                                   double pb_xzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx * fx * pb_x

                     + 0.25 * fx * fx * pb_xy

                     + 0.5 * pa_yy * pb_xy * fx

                     + pa_y * fx * pb_xzz

                     + 0.5 * fx * pb_xyzz

                     + pa_yy * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xyzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyzz,
                                   double pb_xzz,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_y * fx * pb_x * fz * fgb

                     + 4.0 * pa_y * fz * fx * fx * pb_x

                     - 0.5 * fx * pb_xy * fz * fgb

                     - 0.5 * fz * fga * pb_xy * fx

                     - pa_yy * pb_xy * fz * fgb

                     + 2.0 * fz * fx * fx * pb_xy

                     + 5.0 * pa_yy * fz * pb_xy * fx

                     + 10.0 * pa_y * fz * fx * pb_xzz

                     - fz * fga * pb_xyzz

                     + 5.0 * fz * fx * pb_xyzz

                     + 12.0 * pa_yy * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xzzz_s_0(double fx,
                                   double pa_yy,
                                   double pb_xz,
                                   double pb_xzzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_xz

                     + 1.5 * pa_yy * pb_xz * fx

                     + 0.5 * fx * pb_xzzz

                     + pa_yy * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_xzzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yy,
                                   double pb_xz,
                                   double pb_xzzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pb_xz * fx

                     - 3.0 * pa_yy * pb_xz * fz * fgb

                     + 6.0 * fz * fx * fx * pb_xz

                     + 15.0 * pa_yy * fz * pb_xz * fx

                     - fz * fga * pb_xzzz

                     + 5.0 * fz * fx * pb_xzzz

                     + 12.0 * pa_yy * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yyyy_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double pb_yyyy,
                                   double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx

                     + 0.75 * pa_yy * fx * fx

                     + 6.0 * pa_y * fx * fx * pb_y

                     + 4.5 * fx * fx * pb_yy

                     + 3.0 * pa_yy * pb_yy * fx

                     + 4.0 * pa_y * fx * pb_yyy

                     + 0.5 * fx * pb_yyyy

                     + pa_yy * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yyyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double pb_yyyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * fz * fgb

                     + 11.25 * fx * fx * fx * fz

                     - 0.75 * fz * fga * fx * fx

                     - 3.0 * pa_yy * fx * fz * fgb

                     - 12.0 * pa_y * fx * pb_y * fz * fgb

                     + 6.0 * pa_yy * fz * fx * fx

                     + 48.0 * pa_y * fz * fx * fx * pb_y

                     + 36.0 * fx * fx * fz * pb_yy

                     - 3.0 * fx * pb_yy * fz * fgb

                     - 3.0 * fz * fga * pb_yy * fx

                     - 6.0 * pa_yy * pb_yy * fz * fgb

                     + 30.0 * pa_yy * fz * pb_yy * fx

                     + 40.0 * pa_y * fz * fx * pb_yyy

                     - fz * fga * pb_yyyy

                     + 5.0 * fz * fx * pb_yyyy

                     + 12.0 * pa_yy * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yyyz_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_yyyz,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * fx * fx * pb_z

                     + 2.25 * fx * fx * pb_yz

                     + 1.5 * pa_yy * pb_yz * fx

                     + 3.0 * pa_y * fx * pb_yyz

                     + 0.5 * fx * pb_yyyz

                     + pa_yy * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yyyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_yyyz,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_y * fx * fz * fgb * pb_z

                     + 12.0 * pa_y * fz * fx * fx * pb_z

                     + 18.0 * fx * fx * fz * pb_yz

                     - 1.5 * fx * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pb_yz * fx

                     - 3.0 * pa_yy * pb_yz * fz * fgb

                     + 15.0 * pa_yy * fz * pb_yz * fx

                     + 30.0 * pa_y * fz * fx * pb_yyz

                     - fz * fga * pb_yyyz

                     + 5.0 * fz * fx * pb_yyyz

                     + 12.0 * pa_yy * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yyzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyzz,
                                   double pb_yzz,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.25 * pa_yy * fx * fx

                     + pa_y * fx * fx * pb_y

                     + 0.75 * fx * fx * pb_zz

                     + 0.25 * fx * fx * pb_yy

                     + 0.5 * pa_yy * pb_yy * fx

                     + 0.5 * pa_yy * fx * pb_zz

                     + 2.0 * pa_y * fx * pb_yzz

                     + 0.5 * fx * pb_yyzz

                     + pa_yy * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yyzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyzz,
                                   double pb_yzz,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- fx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     - 0.25 * fz * fga * fx * fx

                     - pa_yy * fx * fz * fgb

                     - 2.0 * pa_y * fx * pb_y * fz * fgb

                     + 2.0 * pa_yy * fz * fx * fx

                     + 8.0 * pa_y * fz * fx * fx * pb_y

                     + 6.0 * fx * fx * fz * pb_zz

                     - 0.5 * fx * pb_yy * fz * fgb

                     - 0.5 * fx * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pb_yy * fx

                     - 0.5 * fz * fga * fx * pb_zz

                     - pa_yy * pb_yy * fz * fgb

                     - pa_yy * fz * fgb * pb_zz

                     + 2.0 * fz * fx * fx * pb_yy

                     + 5.0 * pa_yy * fz * pb_yy * fx

                     + 5.0 * pa_yy * fz * fx * pb_zz

                     + 20.0 * pa_y * fz * fx * pb_yzz

                     - fz * fga * pb_yyzz

                     + 5.0 * fz * fx * pb_yyzz

                     + 12.0 * pa_yy * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yzzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_yz,
                                   double pb_yzzz,
                                   double pb_z,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_y * fx * fx * pb_z

                     + 0.75 * fx * fx * pb_yz

                     + 1.5 * pa_yy * pb_yz * fx

                     + pa_y * fx * pb_zzz

                     + 0.5 * fx * pb_yzzz

                     + pa_yy * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_yzzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_yz,
                                   double pb_yzzz,
                                   double pb_z,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_y * fx * pb_z * fz * fgb

                     + 12.0 * pa_y * fz * fx * fx * pb_z

                     - 1.5 * fx * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pb_yz * fx

                     - 3.0 * pa_yy * pb_yz * fz * fgb

                     + 6.0 * fz * fx * fx * pb_yz

                     + 15.0 * pa_yy * fz * pb_yz * fx

                     + 10.0 * pa_y * fz * fx * pb_zzz

                     - fz * fga * pb_yzzz

                     + 5.0 * fz * fx * pb_yzzz

                     + 12.0 * pa_yy * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_zzzz_s_0(double fx,
                                   double pa_yy,
                                   double pb_zz,
                                   double pb_zzzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_yy * fx * fx

                     + 1.5 * fx * fx * pb_zz

                     + 3.0 * pa_yy * pb_zz * fx

                     + 0.5 * fx * pb_zzzz

                     + pa_yy * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yy_zzzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yy,
                                   double pb_zz,
                                   double pb_zzzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * fz * fgb

                     - 0.75 * fz * fga * fx * fx

                     - 3.0 * pa_yy * fx * fz * fgb

                     + 2.25 * fz * fx * fx * fx

                     + 6.0 * pa_yy * fz * fx * fx

                     - 3.0 * fx * pb_zz * fz * fgb

                     - 3.0 * fz * fga * pb_zz * fx

                     - 6.0 * pa_yy * pb_zz * fz * fgb

                     + 12.0 * fz * fx * fx * pb_zz

                     + 30.0 * pa_yy * fz * pb_zz * fx

                     - fz * fga * pb_zzzz

                     + 5.0 * fz * fx * pb_zzzz

                     + 12.0 * pa_yy * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxxx_s_0(double fx,
                                   double pa_yz,
                                   double pb_xx,
                                   double pb_xxxx,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_yz * fx * fx

                     + 3.0 * pa_yz * pb_xx * fx

                     + pa_yz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxxx_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yz,
                                   double pb_xx,
                                   double pb_xxxx,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fx * fz * fgb

                     + 6.0 * pa_yz * fz * fx * fx

                     - 6.0 * pa_yz * pb_xx * fz * fgb

                     + 30.0 * pa_yz * fz * pb_xx * fx

                     + 12.0 * pa_yz * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxxy_s_0(double fx,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xxx,
                                   double pb_xxxy,
                                   double pb_xy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z * pb_x

                     + 1.5 * pa_yz * pb_xy * fx

                     + 0.5 * fx * pa_z * pb_xxx

                     + pa_yz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxxy_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xxx,
                                   double pb_xxxy,
                                   double pb_xy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_z * pb_x * fz * fgb

                     + 6.0 * fx * fx * fz * pa_z * pb_x

                     - 3.0 * pa_yz * pb_xy * fz * fgb

                     + 15.0 * pa_yz * fz * pb_xy * fx

                     + 5.0 * fx * fz * pa_z * pb_xxx

                     + 12.0 * pa_yz * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxxz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xxx,
                                   double pb_xxxz,
                                   double pb_xz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx * pb_x

                     + 1.5 * pa_yz * pb_xz * fx

                     + 0.5 * pa_y * fx * pb_xxx

                     + pa_yz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxxz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xxx,
                                   double pb_xxxz,
                                   double pb_xz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * pb_x * fz * fgb

                     + 6.0 * pa_y * fz * fx * fx * pb_x

                     - 3.0 * pa_yz * pb_xz * fz * fgb

                     + 15.0 * pa_yz * fz * pb_xz * fx

                     + 5.0 * pa_y * fz * fx * pb_xxx

                     + 12.0 * pa_yz * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxyy_s_0(double fx,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_xxyy,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_yz * fx * fx

                     + 0.5 * fx * fx * pa_z * pb_y

                     + 0.5 * pa_yz * pb_xx * fx

                     + 0.5 * pa_yz * fx * pb_yy

                     + fx * pa_z * pb_xxy

                     + pa_yz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxyy_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_xxyy,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_yz * fx * fz * fgb

                     - fx * pa_z * fz * fgb * pb_y

                     + 2.0 * pa_yz * fz * fx * fx

                     + 4.0 * fx * fx * fz * pa_z * pb_y

                     - pa_yz * pb_xx * fz * fgb

                     - pa_yz * fz * fgb * pb_yy

                     + 5.0 * pa_yz * fz * pb_xx * fx

                     + 5.0 * pa_yz * fz * fx * pb_yy

                     + 10.0 * fx * fz * pa_z * pb_xxy

                     + 12.0 * pa_yz * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxyz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_xxyz,
                                   double pb_xxz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * pa_y * fx * fx * pb_y

                     + 0.25 * fx * fx * pa_z * pb_z

                     + 0.25 * fx * fx * pb_xx

                     + 0.5 * pa_yz * fx * pb_yz

                     + 0.5 * pa_y * fx * pb_xxy

                     + 0.5 * fx * pa_z * pb_xxz

                     + pa_yz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxyz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_xx,
                                   double pb_xxy,
                                   double pb_xxyz,
                                   double pb_xxz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * fz * fgb

                     + 0.75 * fx * fx * fx * fz

                     - 0.5 * pa_y * fx * fz * fgb * pb_y

                     - 0.5 * fx * pa_z * fz * fgb * pb_z

                     + 2.0 * pa_y * fz * fx * fx * pb_y

                     + 2.0 * fx * fx * fz * pa_z * pb_z

                     + 2.0 * fx * fx * fz * pb_xx

                     - pa_yz * fz * fgb * pb_yz

                     + 5.0 * pa_yz * fz * fx * pb_yz

                     + 5.0 * pa_y * fz * fx * pb_xxy

                     + 5.0 * fx * fz * pa_z * pb_xxz

                     + 12.0 * pa_yz * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_xxzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_yz * fx * fx

                     + 0.5 * pa_y * fx * fx * pb_z

                     + 0.5 * pa_yz * pb_xx * fx

                     + 0.5 * pa_yz * fx * pb_zz

                     + pa_y * fx * pb_xxz

                     + pa_yz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xxzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_xxzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_yz * fx * fz * fgb

                     - pa_y * fx * fz * fgb * pb_z

                     + 2.0 * pa_yz * fz * fx * fx

                     + 4.0 * pa_y * fz * fx * fx * pb_z

                     - pa_yz * pb_xx * fz * fgb

                     - pa_yz * fz * fgb * pb_zz

                     + 5.0 * pa_yz * fz * pb_xx * fx

                     + 5.0 * pa_yz * fz * fx * pb_zz

                     + 10.0 * pa_y * fz * fx * pb_xxz

                     + 12.0 * pa_yz * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xyyy_s_0(double fx,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_xyyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z * pb_x

                     + 1.5 * pa_yz * pb_xy * fx

                     + 1.5 * fx * pa_z * pb_xyy

                     + pa_yz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xyyy_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_xyyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_z * pb_x * fz * fgb

                     + 6.0 * fx * fx * fz * pa_z * pb_x

                     - 3.0 * pa_yz * pb_xy * fz * fgb

                     + 15.0 * pa_yz * fz * pb_xy * fx

                     + 15.0 * fx * fz * pa_z * pb_xyy

                     + 12.0 * pa_yz * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xyyz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_xyyz,
                                   double pb_xyz,
                                   double pb_xz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_y * fx * fx * pb_x

                     + 0.5 * fx * fx * pb_xy

                     + 0.5 * pa_yz * pb_xz * fx

                     + 0.5 * pa_y * fx * pb_xyy

                     + fx * pa_z * pb_xyz

                     + pa_yz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xyyz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyy,
                                   double pb_xyyz,
                                   double pb_xyz,
                                   double pb_xz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_y * fx * pb_x * fz * fgb

                     + 2.0 * pa_y * fz * fx * fx * pb_x

                     + 4.0 * fx * fx * fz * pb_xy

                     - pa_yz * pb_xz * fz * fgb

                     + 5.0 * pa_yz * fz * pb_xz * fx

                     + 5.0 * pa_y * fz * fx * pb_xyy

                     + 10.0 * fx * fz * pa_z * pb_xyz

                     + 12.0 * pa_yz * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xyzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_xyzz,
                                   double pb_xz,
                                   double pb_xzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_z * pb_x

                     + 0.5 * fx * fx * pb_xz

                     + 0.5 * pa_yz * pb_xy * fx

                     + pa_y * fx * pb_xyz

                     + 0.5 * fx * pa_z * pb_xzz

                     + pa_yz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xyzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_xyzz,
                                   double pb_xz,
                                   double pb_xzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_z * pb_x * fz * fgb

                     + 2.0 * fx * fx * fz * pa_z * pb_x

                     + 4.0 * fx * fx * fz * pb_xz

                     - pa_yz * pb_xy * fz * fgb

                     + 5.0 * pa_yz * fz * pb_xy * fx

                     + 10.0 * pa_y * fz * fx * pb_xyz

                     + 5.0 * fx * fz * pa_z * pb_xzz

                     + 12.0 * pa_yz * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xzzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double pb_xzzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx * pb_x

                     + 1.5 * pa_yz * pb_xz * fx

                     + 1.5 * pa_y * fx * pb_xzz

                     + pa_yz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_xzzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double pb_xzzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * pb_x * fz * fgb

                     + 6.0 * pa_y * fz * fx * fx * pb_x

                     - 3.0 * pa_yz * pb_xz * fz * fgb

                     + 15.0 * pa_yz * fz * pb_xz * fx

                     + 15.0 * pa_y * fz * fx * pb_xzz

                     + 12.0 * pa_yz * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yyyy_s_0(double fx,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double pb_yyyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_yz * fx * fx

                     + 3.0 * fx * fx * pa_z * pb_y

                     + 3.0 * pa_yz * pb_yy * fx

                     + 2.0 * fx * pa_z * pb_yyy

                     + pa_yz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yyyy_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double pb_yyyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fx * fz * fgb

                     - 6.0 * fx * pa_z * pb_y * fz * fgb

                     + 6.0 * pa_yz * fz * fx * fx

                     + 24.0 * fx * fx * fz * pa_z * pb_y

                     - 6.0 * pa_yz * pb_yy * fz * fgb

                     + 30.0 * pa_yz * fz * pb_yy * fx

                     + 20.0 * fx * fz * pa_z * pb_yyy

                     + 12.0 * pa_yz * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yyyz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double pb_yyyz,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_y * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_z * pb_z

                     + 0.75 * fx * fx * pb_yy

                     + 1.5 * pa_yz * pb_yz * fx

                     + 0.5 * pa_y * fx * pb_yyy

                     + 1.5 * fx * pa_z * pb_yyz

                     + pa_yz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yyyz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyy,
                                   double pb_yyyz,
                                   double pb_yyz,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     - 1.5 * pa_y * fx * pb_y * fz * fgb

                     - 1.5 * fx * pa_z * fz * fgb * pb_z

                     + 6.0 * pa_y * fz * fx * fx * pb_y

                     + 6.0 * fx * fx * fz * pa_z * pb_z

                     + 6.0 * fx * fx * fz * pb_yy

                     - 3.0 * pa_yz * pb_yz * fz * fgb

                     + 15.0 * pa_yz * fz * pb_yz * fx

                     + 5.0 * pa_y * fz * fx * pb_yyy

                     + 15.0 * fx * fz * pa_z * pb_yyz

                     + 12.0 * pa_yz * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yyzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_yyzz,
                                   double pb_yz,
                                   double pb_yzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_yz * fx * fx

                     + 0.5 * pa_y * fx * fx * pb_z

                     + 0.5 * fx * fx * pa_z * pb_y

                     + fx * fx * pb_yz

                     + 0.5 * pa_yz * pb_yy * fx

                     + 0.5 * pa_yz * fx * pb_zz

                     + pa_y * fx * pb_yyz

                     + fx * pa_z * pb_yzz

                     + pa_yz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yyzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_yyzz,
                                   double pb_yz,
                                   double pb_yzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_yz * fx * fz * fgb

                     - pa_y * fx * fz * fgb * pb_z

                     - fx * pa_z * pb_y * fz * fgb

                     + 2.0 * pa_yz * fz * fx * fx

                     + 4.0 * pa_y * fz * fx * fx * pb_z

                     + 4.0 * fx * fx * fz * pa_z * pb_y

                     + 8.0 * fx * fx * fz * pb_yz

                     - pa_yz * pb_yy * fz * fgb

                     - pa_yz * fz * fgb * pb_zz

                     + 5.0 * pa_yz * fz * pb_yy * fx

                     + 5.0 * pa_yz * fz * fx * pb_zz

                     + 10.0 * pa_y * fz * fx * pb_yyz

                     + 10.0 * fx * fz * pa_z * pb_yzz

                     + 12.0 * pa_yz * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yzzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double pb_yzzz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_y * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_z * pb_z

                     + 0.75 * fx * fx * pb_zz

                     + 1.5 * pa_yz * pb_yz * fx

                     + 1.5 * pa_y * fx * pb_yzz

                     + 0.5 * fx * pa_z * pb_zzz

                     + pa_yz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_yzzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double pb_yzzz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     - 1.5 * pa_y * fx * pb_y * fz * fgb

                     - 1.5 * fx * pa_z * pb_z * fz * fgb

                     + 6.0 * pa_y * fz * fx * fx * pb_y

                     + 6.0 * fx * fx * fz * pa_z * pb_z

                     + 6.0 * fx * fx * fz * pb_zz

                     - 3.0 * pa_yz * pb_yz * fz * fgb

                     + 15.0 * pa_yz * fz * pb_yz * fx

                     + 15.0 * pa_y * fz * fx * pb_yzz

                     + 5.0 * fx * fz * pa_z * pb_zzz

                     + 12.0 * pa_yz * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_zzzz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double pb_zzzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_yz * fx * fx

                     + 3.0 * pa_y * fx * fx * pb_z

                     + 3.0 * pa_yz * pb_zz * fx

                     + 2.0 * pa_y * fx * pb_zzz

                     + pa_yz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yz_zzzz_r_0(double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double pb_zzzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_yz * fx * fz * fgb

                     - 6.0 * pa_y * fx * pb_z * fz * fgb

                     + 6.0 * pa_yz * fz * fx * fx

                     + 24.0 * pa_y * fz * fx * fx * pb_z

                     - 6.0 * pa_yz * pb_zz * fz * fgb

                     + 30.0 * pa_yz * fz * pb_zz * fx

                     + 20.0 * pa_y * fz * fx * pb_zzz

                     + 12.0 * pa_yz * fz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxxx_s_0(double fx,
                                   double pa_zz,
                                   double pb_xx,
                                   double pb_xxxx,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_zz * fx * fx

                     + 1.5 * fx * fx * pb_xx

                     + 3.0 * pa_zz * pb_xx * fx

                     + 0.5 * fx * pb_xxxx

                     + pa_zz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxxx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_zz,
                                   double pb_xx,
                                   double pb_xxxx,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * fz * fgb

                     - 0.75 * fz * fga * fx * fx

                     - 3.0 * pa_zz * fx * fz * fgb

                     + 2.25 * fz * fx * fx * fx

                     + 6.0 * pa_zz * fz * fx * fx

                     - 3.0 * fx * pb_xx * fz * fgb

                     - 3.0 * fz * fga * pb_xx * fx

                     - 6.0 * pa_zz * pb_xx * fz * fgb

                     + 12.0 * fz * fx * fx * pb_xx

                     + 30.0 * pa_zz * fz * pb_xx * fx

                     - fz * fga * pb_xxxx

                     + 5.0 * fz * fx * pb_xxxx

                     + 12.0 * pa_zz * fz * pb_xxxx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxxy_s_0(double fx,
                                   double pa_zz,
                                   double pb_xxxy,
                                   double pb_xy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_xy

                     + 1.5 * pa_zz * pb_xy * fx

                     + 0.5 * fx * pb_xxxy

                     + pa_zz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxxy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_zz,
                                   double pb_xxxy,
                                   double pb_xy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pb_xy * fx

                     - 3.0 * pa_zz * pb_xy * fz * fgb

                     + 6.0 * fz * fx * fx * pb_xy

                     + 15.0 * pa_zz * fz * pb_xy * fx

                     - fz * fga * pb_xxxy

                     + 5.0 * fz * fx * pb_xxxy

                     + 12.0 * pa_zz * fz * pb_xxxy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxxz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xxx,
                                   double pb_xxxz,
                                   double pb_xz,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * fx * fx * pb_x

                     + 0.75 * fx * fx * pb_xz

                     + 1.5 * pa_zz * pb_xz * fx

                     + pa_z * fx * pb_xxx

                     + 0.5 * fx * pb_xxxz

                     + pa_zz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxxz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xxx,
                                   double pb_xxxz,
                                   double pb_xz,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_z * fx * pb_x * fz * fgb

                     + 12.0 * pa_z * fz * fx * fx * pb_x

                     - 1.5 * fx * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pb_xz * fx

                     - 3.0 * pa_zz * pb_xz * fz * fgb

                     + 6.0 * fz * fx * fx * pb_xz

                     + 15.0 * pa_zz * fz * pb_xz * fx

                     + 10.0 * pa_z * fz * fx * pb_xxx

                     - fz * fga * pb_xxxz

                     + 5.0 * fz * fx * pb_xxxz

                     + 12.0 * pa_zz * fz * pb_xxxz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxyy_s_0(double fx,
                                   double pa_zz,
                                   double pb_xx,
                                   double pb_xxyy,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * pa_zz * fx * fx

                     + 0.25 * fx * fx * pb_xx

                     + 0.25 * fx * fx * pb_yy

                     + 0.5 * pa_zz * pb_xx * fx

                     + 0.5 * pa_zz * fx * pb_yy

                     + 0.5 * fx * pb_xxyy

                     + pa_zz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_zz,
                                   double pb_xx,
                                   double pb_xxyy,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * fx * fz * fgb

                     - 0.25 * fz * fga * fx * fx

                     - pa_zz * fx * fz * fgb

                     + 0.75 * fz * fx * fx * fx

                     + 2.0 * pa_zz * fz * fx * fx

                     - 0.5 * fx * pb_xx * fz * fgb

                     - 0.5 * fx * fz * fgb * pb_yy

                     - 0.5 * fz * fga * pb_xx * fx

                     - 0.5 * fz * fga * fx * pb_yy

                     - pa_zz * pb_xx * fz * fgb

                     - pa_zz * fz * fgb * pb_yy

                     + 2.0 * fz * fx * fx * pb_xx

                     + 2.0 * fz * fx * fx * pb_yy

                     + 5.0 * pa_zz * fz * pb_xx * fx

                     + 5.0 * pa_zz * fz * fx * pb_yy

                     - fz * fga * pb_xxyy

                     + 5.0 * fz * fx * pb_xxyy

                     + 12.0 * pa_zz * fz * pb_xxyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxyz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_xxy,
                                   double pb_xxyz,
                                   double pb_y,
                                   double pb_yz,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_z * fx * fx * pb_y

                     + 0.25 * fx * fx * pb_yz

                     + 0.5 * pa_zz * fx * pb_yz

                     + pa_z * fx * pb_xxy

                     + 0.5 * fx * pb_xxyz

                     + pa_zz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_xxy,
                                   double pb_xxyz,
                                   double pb_y,
                                   double pb_yz,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_z * fx * fz * fgb * pb_y

                     + 4.0 * pa_z * fz * fx * fx * pb_y

                     - 0.5 * fx * fz * fgb * pb_yz

                     - 0.5 * fz * fga * fx * pb_yz

                     - pa_zz * fz * fgb * pb_yz

                     + 2.0 * fz * fx * fx * pb_yz

                     + 5.0 * pa_zz * fz * fx * pb_yz

                     + 10.0 * pa_z * fz * fx * pb_xxy

                     - fz * fga * pb_xxyz

                     + 5.0 * fz * fx * pb_xxyz

                     + 12.0 * pa_zz * fz * pb_xxyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxzz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_xxzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.25 * pa_zz * fx * fx

                     + pa_z * fx * fx * pb_z

                     + 0.75 * fx * fx * pb_xx

                     + 0.25 * fx * fx * pb_zz

                     + 0.5 * pa_zz * pb_xx * fx

                     + 0.5 * pa_zz * fx * pb_zz

                     + 2.0 * pa_z * fx * pb_xxz

                     + 0.5 * fx * pb_xxzz

                     + pa_zz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xxzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_xx,
                                   double pb_xxz,
                                   double pb_xxzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- fx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     - 0.25 * fz * fga * fx * fx

                     - pa_zz * fx * fz * fgb

                     - 2.0 * pa_z * fx * fz * fgb * pb_z

                     + 2.0 * pa_zz * fz * fx * fx

                     + 8.0 * pa_z * fz * fx * fx * pb_z

                     + 6.0 * fx * fx * fz * pb_xx

                     - 0.5 * fx * pb_xx * fz * fgb

                     - 0.5 * fx * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pb_xx * fx

                     - 0.5 * fz * fga * fx * pb_zz

                     - pa_zz * pb_xx * fz * fgb

                     - pa_zz * fz * fgb * pb_zz

                     + 2.0 * fz * fx * fx * pb_zz

                     + 5.0 * pa_zz * fz * pb_xx * fx

                     + 5.0 * pa_zz * fz * fx * pb_zz

                     + 20.0 * pa_z * fz * fx * pb_xxz

                     - fz * fga * pb_xxzz

                     + 5.0 * fz * fx * pb_xxzz

                     + 12.0 * pa_zz * fz * pb_xxzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xyyy_s_0(double fx,
                                   double pa_zz,
                                   double pb_xy,
                                   double pb_xyyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_xy

                     + 1.5 * pa_zz * pb_xy * fx

                     + 0.5 * fx * pb_xyyy

                     + pa_zz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xyyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_zz,
                                   double pb_xy,
                                   double pb_xyyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pb_xy * fz * fgb

                     - 1.5 * fz * fga * pb_xy * fx

                     - 3.0 * pa_zz * pb_xy * fz * fgb

                     + 6.0 * fz * fx * fx * pb_xy

                     + 15.0 * pa_zz * fz * pb_xy * fx

                     - fz * fga * pb_xyyy

                     + 5.0 * fz * fx * pb_xyyy

                     + 12.0 * pa_zz * fz * pb_xyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xyyz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xyy,
                                   double pb_xyyz,
                                   double pb_xz,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_z * fx * fx * pb_x

                     + 0.25 * fx * fx * pb_xz

                     + 0.5 * pa_zz * pb_xz * fx

                     + pa_z * fx * pb_xyy

                     + 0.5 * fx * pb_xyyz

                     + pa_zz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xyyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xyy,
                                   double pb_xyyz,
                                   double pb_xz,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_z * fx * pb_x * fz * fgb

                     + 4.0 * pa_z * fz * fx * fx * pb_x

                     - 0.5 * fx * pb_xz * fz * fgb

                     - 0.5 * fz * fga * pb_xz * fx

                     - pa_zz * pb_xz * fz * fgb

                     + 2.0 * fz * fx * fx * pb_xz

                     + 5.0 * pa_zz * fz * pb_xz * fx

                     + 10.0 * pa_z * fz * fx * pb_xyy

                     - fz * fga * pb_xyyz

                     + 5.0 * fz * fx * pb_xyyz

                     + 12.0 * pa_zz * fz * pb_xyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xyzz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_xyzz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_xy

                     + 0.5 * pa_zz * pb_xy * fx

                     + 2.0 * pa_z * fx * pb_xyz

                     + 0.5 * fx * pb_xyzz

                     + pa_zz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xyzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_xy,
                                   double pb_xyz,
                                   double pb_xyzz,
                                   double r_0_0)
    {
        return r_0_0 * (6.0 * fx * fx * fz * pb_xy

                     - 0.5 * fx * pb_xy * fz * fgb

                     - 0.5 * fz * fga * pb_xy * fx

                     - pa_zz * pb_xy * fz * fgb

                     + 5.0 * pa_zz * fz * pb_xy * fx

                     + 20.0 * pa_z * fz * fx * pb_xyz

                     - fz * fga * pb_xyzz

                     + 5.0 * fz * fx * pb_xyzz

                     + 12.0 * pa_zz * fz * pb_xyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xzzz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double pb_xzzz,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * fx * fx * pb_x

                     + 2.25 * fx * fx * pb_xz

                     + 1.5 * pa_zz * pb_xz * fx

                     + 3.0 * pa_z * fx * pb_xzz

                     + 0.5 * fx * pb_xzzz

                     + pa_zz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_xzzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_xzz,
                                   double pb_xzzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_z * fx * pb_x * fz * fgb

                     + 12.0 * pa_z * fz * fx * fx * pb_x

                     + 18.0 * fx * fx * fz * pb_xz

                     - 1.5 * fx * pb_xz * fz * fgb

                     - 1.5 * fz * fga * pb_xz * fx

                     - 3.0 * pa_zz * pb_xz * fz * fgb

                     + 15.0 * pa_zz * fz * pb_xz * fx

                     + 30.0 * pa_z * fz * fx * pb_xzz

                     - fz * fga * pb_xzzz

                     + 5.0 * fz * fx * pb_xzzz

                     + 12.0 * pa_zz * fz * pb_xzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yyyy_s_0(double fx,
                                   double pa_zz,
                                   double pb_yy,
                                   double pb_yyyy,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_zz * fx * fx

                     + 1.5 * fx * fx * pb_yy

                     + 3.0 * pa_zz * pb_yy * fx

                     + 0.5 * fx * pb_yyyy

                     + pa_zz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yyyy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_zz,
                                   double pb_yy,
                                   double pb_yyyy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fx * fz * fgb

                     - 0.75 * fz * fga * fx * fx

                     - 3.0 * pa_zz * fx * fz * fgb

                     + 2.25 * fz * fx * fx * fx

                     + 6.0 * pa_zz * fz * fx * fx

                     - 3.0 * fx * pb_yy * fz * fgb

                     - 3.0 * fz * fga * pb_yy * fx

                     - 6.0 * pa_zz * pb_yy * fz * fgb

                     + 12.0 * fz * fx * fx * pb_yy

                     + 30.0 * pa_zz * fz * pb_yy * fx

                     - fz * fga * pb_yyyy

                     + 5.0 * fz * fx * pb_yyyy

                     + 12.0 * pa_zz * fz * pb_yyyy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yyyz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yyy,
                                   double pb_yyyz,
                                   double pb_yz,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * fx * fx * pb_y

                     + 0.75 * fx * fx * pb_yz

                     + 1.5 * pa_zz * pb_yz * fx

                     + pa_z * fx * pb_yyy

                     + 0.5 * fx * pb_yyyz

                     + pa_zz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yyyz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yyy,
                                   double pb_yyyz,
                                   double pb_yz,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_z * fx * pb_y * fz * fgb

                     + 12.0 * pa_z * fz * fx * fx * pb_y

                     - 1.5 * fx * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pb_yz * fx

                     - 3.0 * pa_zz * pb_yz * fz * fgb

                     + 6.0 * fz * fx * fx * pb_yz

                     + 15.0 * pa_zz * fz * pb_yz * fx

                     + 10.0 * pa_z * fz * fx * pb_yyy

                     - fz * fga * pb_yyyz

                     + 5.0 * fz * fx * pb_yyyz

                     + 12.0 * pa_zz * fz * pb_yyyz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yyzz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_yyzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.25 * pa_zz * fx * fx

                     + pa_z * fx * fx * pb_z

                     + 0.75 * fx * fx * pb_yy

                     + 0.25 * fx * fx * pb_zz

                     + 0.5 * pa_zz * pb_yy * fx

                     + 0.5 * pa_zz * fx * pb_zz

                     + 2.0 * pa_z * fx * pb_yyz

                     + 0.5 * fx * pb_yyzz

                     + pa_zz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yyzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_yy,
                                   double pb_yyz,
                                   double pb_yyzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- fx * fx * fz * fgb

                     + 2.25 * fx * fx * fx * fz

                     - 0.25 * fz * fga * fx * fx

                     - pa_zz * fx * fz * fgb

                     - 2.0 * pa_z * fx * fz * fgb * pb_z

                     + 2.0 * pa_zz * fz * fx * fx

                     + 8.0 * pa_z * fz * fx * fx * pb_z

                     + 6.0 * fx * fx * fz * pb_yy

                     - 0.5 * fx * pb_yy * fz * fgb

                     - 0.5 * fx * fz * fgb * pb_zz

                     - 0.5 * fz * fga * pb_yy * fx

                     - 0.5 * fz * fga * fx * pb_zz

                     - pa_zz * pb_yy * fz * fgb

                     - pa_zz * fz * fgb * pb_zz

                     + 2.0 * fz * fx * fx * pb_zz

                     + 5.0 * pa_zz * fz * pb_yy * fx

                     + 5.0 * pa_zz * fz * fx * pb_zz

                     + 20.0 * pa_z * fz * fx * pb_yyz

                     - fz * fga * pb_yyzz

                     + 5.0 * fz * fx * pb_yyzz

                     + 12.0 * pa_zz * fz * pb_yyzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yzzz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double pb_yzzz,
                                   double s_0_0)
    {
        return s_0_0 * (1.5 * pa_z * fx * fx * pb_y

                     + 2.25 * fx * fx * pb_yz

                     + 1.5 * pa_zz * pb_yz * fx

                     + 3.0 * pa_z * fx * pb_yzz

                     + 0.5 * fx * pb_yzzz

                     + pa_zz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_yzzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_yzz,
                                   double pb_yzzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * pa_z * fx * pb_y * fz * fgb

                     + 12.0 * pa_z * fz * fx * fx * pb_y

                     + 18.0 * fx * fx * fz * pb_yz

                     - 1.5 * fx * pb_yz * fz * fgb

                     - 1.5 * fz * fga * pb_yz * fx

                     - 3.0 * pa_zz * pb_yz * fz * fgb

                     + 15.0 * pa_zz * fz * pb_yz * fx

                     + 30.0 * pa_z * fz * fx * pb_yzz

                     - fz * fga * pb_yzzz

                     + 5.0 * fz * fx * pb_yzzz

                     + 12.0 * pa_zz * fz * pb_yzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_zzzz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double pb_zzzz,
                                   double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx

                     + 0.75 * pa_zz * fx * fx

                     + 6.0 * pa_z * fx * fx * pb_z

                     + 4.5 * fx * fx * pb_zz

                     + 3.0 * pa_zz * pb_zz * fx

                     + 4.0 * pa_z * fx * pb_zzz

                     + 0.5 * fx * pb_zzzz

                     + pa_zz * pb_zzzz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zz_zzzz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_z,
                                   double pb_zz,
                                   double pb_zzz,
                                   double pb_zzzz,
                                   double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * fz * fgb

                     + 11.25 * fx * fx * fx * fz

                     - 0.75 * fz * fga * fx * fx

                     - 3.0 * pa_zz * fx * fz * fgb

                     - 12.0 * pa_z * fx * pb_z * fz * fgb

                     + 6.0 * pa_zz * fz * fx * fx

                     + 48.0 * pa_z * fz * fx * fx * pb_z

                     + 36.0 * fx * fx * fz * pb_zz

                     - 3.0 * fx * pb_zz * fz * fgb

                     - 3.0 * fz * fga * pb_zz * fx

                     - 6.0 * pa_zz * pb_zz * fz * fgb

                     + 30.0 * pa_zz * fz * pb_zz * fx

                     + 40.0 * pa_z * fz * fx * pb_zzz

                     - fz * fga * pb_zzzz

                     + 5.0 * fz * fx * pb_zzzz

                     + 12.0 * pa_zz * fz * pb_zzzz);

    }

    // SIMD elementary functions for (G|T|D) integrals

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xx_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pa_xxxx,
                                   double pb_x,
                                   double pb_xx,
                                   double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx

                     + 4.5 * pa_xx * fx * fx

                     + 6.0 * pa_x * fx * fx * pb_x

                     + 0.5 * pa_xxxx * fx

                     + 4.0 * pa_xxx * fx * pb_x

                     + 0.75 * fx * fx * pb_xx

                     + 3.0 * pa_xx * fx * pb_xx

                     + pa_xxxx * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pa_xxxx,
                                   double pb_x,
                                   double pb_xx,
                                   double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * fz * fga

                     + 11.25 * fx * fx * fx * fz

                     + 36.0 * pa_xx * fx * fx * fz

                     - 0.75 * fx * fx * fz * fgb

                     - 3.0 * pa_xx * fx * fz * fgb

                     - 3.0 * pa_xx * fz * fga * fx

                     - 12.0 * pa_x * fx * fz * fga * pb_x

                     - pa_xxxx * fz * fgb

                     + 48.0 * pa_x * fx * fx * fz * pb_x

                     + 5.0 * pa_xxxx * fz * fx

                     + 40.0 * pa_xxx * fz * fx * pb_x

                     - 3.0 * fx * fz * fga * pb_xx

                     - 6.0 * pa_xx * fz * fga * pb_xx

                     + 6.0 * fx * fx * fz * pb_xx

                     + 30.0 * pa_xx * fz * fx * pb_xx

                     + 12.0 * pa_xxxx * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xy_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pa_xxxx,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (3.0 * pa_x * fx * fx * pb_y

                     + 2.0 * pa_xxx * fx * pb_y

                     + 0.75 * fx * fx * pb_xy

                     + 3.0 * pa_xx * fx * pb_xy

                     + pa_xxxx * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xy_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pa_xxxx,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- 6.0 * pa_x * fx * fz * fga * pb_y

                     + 24.0 * pa_x * fx * fx * fz * pb_y

                     + 20.0 * pa_xxx * fz * fx * pb_y

                     - 3.0 * fx * fz * fga * pb_xy

                     - 6.0 * pa_xx * fz * fga * pb_xy

                     + 6.0 * fx * fx * fz * pb_xy

                     + 30.0 * pa_xx * fz * fx * pb_xy

                     + 12.0 * pa_xxxx * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pa_xxxx,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (3.0 * pa_x * fx * fx * pb_z

                     + 2.0 * pa_xxx * fx * pb_z

                     + 0.75 * fx * fx * pb_xz

                     + 3.0 * pa_xx * fx * pb_xz

                     + pa_xxxx * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_xz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pa_xxxx,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 6.0 * pa_x * fx * fz * fga * pb_z

                     + 24.0 * pa_x * fx * fx * fz * pb_z

                     + 20.0 * pa_xxx * fz * fx * pb_z

                     - 3.0 * fx * fz * fga * pb_xz

                     - 6.0 * pa_xx * fz * fga * pb_xz

                     + 6.0 * fx * fx * fz * pb_xz

                     + 30.0 * pa_xx * fz * fx * pb_xz

                     + 12.0 * pa_xxxx * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yy_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxxx,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 1.5 * pa_xx * fx * fx

                     + 0.5 * pa_xxxx * fx

                     + 0.75 * fx * fx * pb_yy

                     + 3.0 * pa_xx * fx * pb_yy

                     + pa_xxxx * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxxx,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 1.5 * fx * fx * fz * fga

                     - 3.0 * pa_xx * fx * fz * fgb

                     - 3.0 * pa_xx * fz * fga * fx

                     + 2.25 * fx * fx * fx * fz

                     - pa_xxxx * fz * fgb

                     + 12.0 * pa_xx * fz * fx * fx

                     + 5.0 * pa_xxxx * fz * fx

                     - 3.0 * fx * fz * fga * pb_yy

                     - 6.0 * pa_xx * fz * fga * pb_yy

                     + 6.0 * fx * fx * fz * pb_yy

                     + 30.0 * pa_xx * fz * fx * pb_yy

                     + 12.0 * pa_xxxx * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yz_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxxx,
                                   double pb_yz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_yz

                     + 3.0 * pa_xx * fx * pb_yz

                     + pa_xxxx * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_yz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxxx,
                                   double pb_yz,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * fx * fz * fga * pb_yz

                     - 6.0 * pa_xx * fz * fga * pb_yz

                     + 6.0 * fx * fx * fz * pb_yz

                     + 30.0 * pa_xx * fz * fx * pb_yz

                     + 12.0 * pa_xxxx * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_zz_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxxx,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 1.5 * pa_xx * fx * fx

                     + 0.5 * pa_xxxx * fx

                     + 0.75 * fx * fx * pb_zz

                     + 3.0 * pa_xx * fx * pb_zz

                     + pa_xxxx * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxx_zz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxxx,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 1.5 * fx * fx * fz * fga

                     - 3.0 * pa_xx * fx * fz * fgb

                     - 3.0 * pa_xx * fz * fga * fx

                     + 2.25 * fx * fx * fx * fz

                     - pa_xxxx * fz * fgb

                     + 12.0 * pa_xx * fz * fx * fx

                     + 5.0 * pa_xxxx * fz * fx

                     - 3.0 * fx * fz * fga * pb_zz

                     - 6.0 * pa_xx * fz * fga * pb_zz

                     + 6.0 * fx * fx * fz * pb_zz

                     + 30.0 * pa_xx * fz * fx * pb_zz

                     + 12.0 * pa_xxxx * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xx_s_0(double fx,
                                   double pa_xxxy,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xx,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_xy * fx * fx

                     + 1.5 * fx * fx * pa_y * pb_x

                     + 0.5 * pa_xxxy * fx

                     + 3.0 * pa_xxy * fx * pb_x

                     + 1.5 * pa_xy * fx * pb_xx

                     + pa_xxxy * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xxxy,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xx,
                                   double r_0_0)
    {
        return r_0_0 * (18.0 * pa_xy * fx * fx * fz

                     - 1.5 * pa_xy * fx * fz * fgb

                     - 1.5 * pa_xy * fz * fga * fx

                     - 3.0 * fx * fz * fga * pa_y * pb_x

                     - pa_xxxy * fz * fgb

                     + 12.0 * fx * fx * fz * pa_y * pb_x

                     + 5.0 * pa_xxxy * fz * fx

                     + 30.0 * pa_xxy * fx * fz * pb_x

                     - 3.0 * pa_xy * fz * fga * pb_xx

                     + 15.0 * pa_xy * fx * fz * pb_xx

                     + 12.0 * pa_xxxy * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xy_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pa_xxxy,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_xx * fx * fx

                     + 0.75 * pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_y * pb_y

                     + 0.5 * pa_xxx * fx * pb_x

                     + 1.5 * pa_xxy * fx * pb_y

                     + 1.5 * pa_xy * fx * pb_xy

                     + pa_xxxy * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xy_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pa_xxxy,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fga

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * pa_xx * fx * fx * fz

                     - 1.5 * pa_x * fz * fga * fx * pb_x

                     - 1.5 * fx * fz * fga * pa_y * pb_y

                     + 6.0 * pa_x * fx * fx * fz * pb_x

                     + 6.0 * fx * fx * fz * pa_y * pb_y

                     + 5.0 * pa_xxx * fz * fx * pb_x

                     + 15.0 * pa_xxy * fx * fz * pb_y

                     - 3.0 * pa_xy * fz * fga * pb_xy

                     + 15.0 * pa_xy * fx * fz * pb_xy

                     + 12.0 * pa_xxxy * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xz_s_0(double fx,
                                   double pa_xxxy,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_y * pb_z

                     + 1.5 * pa_xxy * fx * pb_z

                     + 1.5 * pa_xy * fx * pb_xz

                     + pa_xxxy * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_xz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_xxxy,
                                   double pa_xxy,
                                   double pa_xy,
                                   double pa_y,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fz * fga * pa_y * pb_z

                     + 6.0 * fx * fx * fz * pa_y * pb_z

                     + 15.0 * pa_xxy * fx * fz * pb_z

                     - 3.0 * pa_xy * fz * fga * pb_xz

                     + 15.0 * pa_xy * fx * fz * pb_xz

                     + 12.0 * pa_xxxy * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yy_s_0(double fx,
                                   double pa_x,
                                   double pa_xxx,
                                   double pa_xxxy,
                                   double pa_xy,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xy * fx * fx

                     + 1.5 * pa_x * fx * fx * pb_y

                     + 0.5 * pa_xxxy * fx

                     + pa_xxx * fx * pb_y

                     + 1.5 * pa_xy * fx * pb_yy

                     + pa_xxxy * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xxx,
                                   double pa_xxxy,
                                   double pa_xy,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xy * fx * fz * fgb

                     - 1.5 * pa_xy * fz * fga * fx

                     - 3.0 * pa_x * fz * fga * fx * pb_y

                     - pa_xxxy * fz * fgb

                     + 6.0 * pa_xy * fx * fx * fz

                     + 12.0 * pa_x * fx * fx * fz * pb_y

                     + 5.0 * pa_xxxy * fz * fx

                     + 10.0 * pa_xxx * fz * fx * pb_y

                     - 3.0 * pa_xy * fz * fga * pb_yy

                     + 15.0 * pa_xy * fx * fz * pb_yy

                     + 12.0 * pa_xxxy * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yz_s_0(double fx,
                                   double pa_x,
                                   double pa_xxx,
                                   double pa_xxxy,
                                   double pa_xy,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * pb_z

                     + 0.5 * pa_xxx * fx * pb_z

                     + 1.5 * pa_xy * fx * pb_yz

                     + pa_xxxy * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_yz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xxx,
                                   double pa_xxxy,
                                   double pa_xy,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fz * fga * fx * pb_z

                     + 6.0 * pa_x * fx * fx * fz * pb_z

                     + 5.0 * pa_xxx * fz * fx * pb_z

                     - 3.0 * pa_xy * fz * fga * pb_yz

                     + 15.0 * pa_xy * fx * fz * pb_yz

                     + 12.0 * pa_xxxy * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_zz_s_0(double fx,
                                   double pa_xxxy,
                                   double pa_xy,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xy * fx * fx

                     + 0.5 * pa_xxxy * fx

                     + 1.5 * pa_xy * fx * pb_zz

                     + pa_xxxy * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxy_zz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xxxy,
                                   double pa_xy,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xy * fx * fz * fgb

                     - 1.5 * pa_xy * fz * fga * fx

                     - pa_xxxy * fz * fgb

                     + 6.0 * pa_xy * fx * fx * fz

                     + 5.0 * pa_xxxy * fz * fx

                     - 3.0 * pa_xy * fz * fga * pb_zz

                     + 15.0 * pa_xy * fx * fz * pb_zz

                     + 12.0 * pa_xxxy * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xx_s_0(double fx,
                                   double pa_xxxz,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_xz * fx * fx

                     + 1.5 * fx * fx * pa_z * pb_x

                     + 0.5 * pa_xxxz * fx

                     + 3.0 * pa_xxz * fx * pb_x

                     + 1.5 * pa_xz * fx * pb_xx

                     + pa_xxxz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xxxz,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double r_0_0)
    {
        return r_0_0 * (18.0 * pa_xz * fx * fx * fz

                     - 1.5 * pa_xz * fx * fz * fgb

                     - 1.5 * pa_xz * fz * fga * fx

                     - 3.0 * fx * fz * fga * pa_z * pb_x

                     - pa_xxxz * fz * fgb

                     + 12.0 * fx * fx * fz * pa_z * pb_x

                     + 5.0 * pa_xxxz * fz * fx

                     + 30.0 * pa_xxz * fx * fz * pb_x

                     - 3.0 * pa_xz * fz * fga * pb_xx

                     + 15.0 * pa_xz * fx * fz * pb_xx

                     + 12.0 * pa_xxxz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xy_s_0(double fx,
                                   double pa_xxxz,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z * pb_y

                     + 1.5 * pa_xxz * fx * pb_y

                     + 1.5 * pa_xz * fx * pb_xy

                     + pa_xxxz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xy_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_xxxz,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fz * fga * pa_z * pb_y

                     + 6.0 * fx * fx * fz * pa_z * pb_y

                     + 15.0 * pa_xxz * fx * fz * pb_y

                     - 3.0 * pa_xz * fz * fga * pb_xy

                     + 15.0 * pa_xz * fx * fz * pb_xy

                     + 12.0 * pa_xxxz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pa_xxxz,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_xx * fx * fx

                     + 0.75 * pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_z * pb_z

                     + 0.5 * pa_xxx * fx * pb_x

                     + 1.5 * pa_xxz * fx * pb_z

                     + 1.5 * pa_xz * fx * pb_xz

                     + pa_xxxz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_xz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxx,
                                   double pa_xxxz,
                                   double pa_xxz,
                                   double pa_xz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fga

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * pa_xx * fx * fx * fz

                     - 1.5 * pa_x * fz * fga * fx * pb_x

                     - 1.5 * fx * fz * fga * pa_z * pb_z

                     + 6.0 * pa_x * fx * fx * fz * pb_x

                     + 6.0 * fx * fx * fz * pa_z * pb_z

                     + 5.0 * pa_xxx * fz * fx * pb_x

                     + 15.0 * pa_xxz * fx * fz * pb_z

                     - 3.0 * pa_xz * fz * fga * pb_xz

                     + 15.0 * pa_xz * fx * fz * pb_xz

                     + 12.0 * pa_xxxz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yy_s_0(double fx,
                                   double pa_xxxz,
                                   double pa_xz,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xz * fx * fx

                     + 0.5 * pa_xxxz * fx

                     + 1.5 * pa_xz * fx * pb_yy

                     + pa_xxxz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xxxz,
                                   double pa_xz,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xz * fx * fz * fgb

                     - 1.5 * pa_xz * fz * fga * fx

                     - pa_xxxz * fz * fgb

                     + 6.0 * pa_xz * fx * fx * fz

                     + 5.0 * pa_xxxz * fz * fx

                     - 3.0 * pa_xz * fz * fga * pb_yy

                     + 15.0 * pa_xz * fx * fz * pb_yy

                     + 12.0 * pa_xxxz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yz_s_0(double fx,
                                   double pa_x,
                                   double pa_xxx,
                                   double pa_xxxz,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * pb_y

                     + 0.5 * pa_xxx * fx * pb_y

                     + 1.5 * pa_xz * fx * pb_yz

                     + pa_xxxz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_yz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xxx,
                                   double pa_xxxz,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fz * fga * fx * pb_y

                     + 6.0 * pa_x * fx * fx * fz * pb_y

                     + 5.0 * pa_xxx * fz * fx * pb_y

                     - 3.0 * pa_xz * fz * fga * pb_yz

                     + 15.0 * pa_xz * fx * fz * pb_yz

                     + 12.0 * pa_xxxz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_zz_s_0(double fx,
                                   double pa_x,
                                   double pa_xxx,
                                   double pa_xxxz,
                                   double pa_xz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xz * fx * fx

                     + 1.5 * pa_x * fx * fx * pb_z

                     + 0.5 * pa_xxxz * fx

                     + pa_xxx * fx * pb_z

                     + 1.5 * pa_xz * fx * pb_zz

                     + pa_xxxz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxxz_zz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xxx,
                                   double pa_xxxz,
                                   double pa_xz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xz * fx * fz * fgb

                     - 1.5 * pa_xz * fz * fga * fx

                     - 3.0 * pa_x * fz * fga * fx * pb_z

                     - pa_xxxz * fz * fgb

                     + 6.0 * pa_xz * fx * fx * fz

                     + 12.0 * pa_x * fx * fx * fz * pb_z

                     + 5.0 * pa_xxxz * fz * fx

                     + 10.0 * pa_xxx * fz * fx * pb_z

                     - 3.0 * pa_xz * fz * fga * pb_zz

                     + 15.0 * pa_xz * fx * fz * pb_zz

                     + 12.0 * pa_xxxz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xx_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxyy,
                                   double pa_xyy,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xx,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * fx * fx * pa_yy

                     + 0.25 * pa_xx * fx * fx

                     + pa_x * fx * fx * pb_x

                     + 0.5 * pa_xxyy * fx

                     + 2.0 * pa_xyy * fx * pb_x

                     + 0.25 * fx * fx * pb_xx

                     + 0.5 * pa_xx * fx * pb_xx

                     + 0.5 * fx * pa_yy * pb_xx

                     + pa_xxyy * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxyy,
                                   double pa_xyy,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xx,
                                   double r_0_0)
    {
        return r_0_0 * (- fx * fx * fz * fga

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * fx * fx * pa_yy * fz

                     - 0.25 * fx * fx * fz * fgb

                     - 0.5 * pa_xx * fx * fz * fgb

                     - 0.5 * pa_xx * fz * fga * fx

                     - 2.0 * pa_x * fx * fz * fga * pb_x

                     - 0.5 * fx * pa_yy * fz * fgb

                     - 0.5 * fz * fga * pa_yy * fx

                     - pa_xxyy * fz * fgb

                     + 2.0 * pa_xx * fz * fx * fx

                     + 8.0 * pa_x * fx * fx * fz * pb_x

                     + 5.0 * pa_xxyy * fz * fx

                     + 20.0 * pa_xyy * fx * fz * pb_x

                     - fx * fz * fga * pb_xx

                     - pa_xx * fz * fga * pb_xx

                     + 2.0 * fx * fx * fz * pb_xx

                     - fz * fga * pa_yy * pb_xx

                     + 5.0 * pa_xx * fz * fx * pb_xx

                     + 5.0 * fx * pa_yy * fz * pb_xx

                     + 12.0 * pa_xxyy * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xy_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_xxyy,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (pa_xy * fx * fx

                     + 0.5 * pa_x * fx * fx * pb_y

                     + 0.5 * fx * fx * pa_y * pb_x

                     + pa_xxy * fx * pb_x

                     + pa_xyy * fx * pb_y

                     + 0.25 * fx * fx * pb_xy

                     + 0.5 * pa_xx * fx * pb_xy

                     + 0.5 * fx * pa_yy * pb_xy

                     + pa_xxyy * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xy_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_xxyy,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (8.0 * pa_xy * fx * fx * fz

                     - pa_x * fx * fz * fga * pb_y

                     - fz * fga * pa_y * fx * pb_x

                     + 4.0 * pa_x * fx * fx * fz * pb_y

                     + 4.0 * fx * fx * pa_y * fz * pb_x

                     + 10.0 * pa_xxy * fz * fx * pb_x

                     + 10.0 * pa_xyy * fx * fz * pb_y

                     - fx * fz * fga * pb_xy

                     - pa_xx * fz * fga * pb_xy

                     + 2.0 * fx * fx * fz * pb_xy

                     - fz * fga * pa_yy * pb_xy

                     + 5.0 * pa_xx * fz * fx * pb_xy

                     + 5.0 * fx * pa_yy * fz * pb_xy

                     + 12.0 * pa_xxyy * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxyy,
                                   double pa_xyy,
                                   double pa_yy,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * fx * pb_z

                     + pa_xyy * fx * pb_z

                     + 0.25 * fx * fx * pb_xz

                     + 0.5 * pa_xx * fx * pb_xz

                     + 0.5 * fx * pa_yy * pb_xz

                     + pa_xxyy * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_xz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxyy,
                                   double pa_xyy,
                                   double pa_yy,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_x * fx * fz * fga * pb_z

                     + 4.0 * pa_x * fx * fx * fz * pb_z

                     + 10.0 * pa_xyy * fx * fz * pb_z

                     - fx * fz * fga * pb_xz

                     - pa_xx * fz * fga * pb_xz

                     + 2.0 * fx * fx * fz * pb_xz

                     - fz * fga * pa_yy * pb_xz

                     + 5.0 * pa_xx * fz * fx * pb_xz

                     + 5.0 * fx * pa_yy * fz * pb_xz

                     + 12.0 * pa_xxyy * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yy_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_xxyy,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_xx * fx * fx

                     + 0.25 * fx * fx * pa_yy

                     + fx * fx * pa_y * pb_y

                     + 0.5 * pa_xxyy * fx

                     + 2.0 * pa_xxy * fx * pb_y

                     + 0.25 * fx * fx * pb_yy

                     + 0.5 * pa_xx * fx * pb_yy

                     + 0.5 * fx * pa_yy * pb_yy

                     + pa_xxyy * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_xxyy,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- fz * fga * fx * fx

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * pa_xx * fx * fx * fz

                     - 0.25 * fx * fx * fz * fgb

                     - 0.5 * pa_xx * fx * fz * fgb

                     - 0.5 * pa_xx * fz * fga * fx

                     - 0.5 * fx * pa_yy * fz * fgb

                     - 0.5 * fz * fga * pa_yy * fx

                     - 2.0 * fz * fga * pa_y * fx * pb_y

                     - pa_xxyy * fz * fgb

                     + 2.0 * fx * fx * pa_yy * fz

                     + 8.0 * fx * fx * pa_y * fz * pb_y

                     + 5.0 * pa_xxyy * fz * fx

                     + 20.0 * pa_xxy * fz * fx * pb_y

                     - fx * fz * fga * pb_yy

                     - pa_xx * fz * fga * pb_yy

                     + 2.0 * fx * fx * fz * pb_yy

                     - fz * fga * pa_yy * pb_yy

                     + 5.0 * pa_xx * fz * fx * pb_yy

                     + 5.0 * fx * pa_yy * fz * pb_yy

                     + 12.0 * pa_xxyy * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yz_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_xxyy,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * fx * fx * pa_y * pb_z

                     + pa_xxy * fx * pb_z

                     + 0.25 * fx * fx * pb_yz

                     + 0.5 * pa_xx * fx * pb_yz

                     + 0.5 * fx * pa_yy * pb_yz

                     + pa_xxyy * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_yz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_xxyy,
                                   double pa_y,
                                   double pa_yy,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- fz * fga * pa_y * fx * pb_z

                     + 4.0 * fx * fx * pa_y * fz * pb_z

                     + 10.0 * pa_xxy * fz * fx * pb_z

                     - fx * fz * fga * pb_yz

                     - pa_xx * fz * fga * pb_yz

                     + 2.0 * fx * fx * fz * pb_yz

                     - fz * fga * pa_yy * pb_yz

                     + 5.0 * pa_xx * fz * fx * pb_yz

                     + 5.0 * fx * pa_yy * fz * pb_yz

                     + 12.0 * pa_xxyy * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_zz_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxyy,
                                   double pa_yy,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * pa_xx * fx * fx

                     + 0.25 * fx * fx * pa_yy

                     + 0.5 * pa_xxyy * fx

                     + 0.25 * fx * fx * pb_zz

                     + 0.5 * pa_xx * fx * pb_zz

                     + 0.5 * fx * pa_yy * pb_zz

                     + pa_xxyy * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyy_zz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxyy,
                                   double pa_yy,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * fz * fgb

                     - 0.5 * fx * fx * fz * fga

                     - 0.5 * pa_xx * fx * fz * fgb

                     - 0.5 * pa_xx * fz * fga * fx

                     - 0.5 * fx * pa_yy * fz * fgb

                     + 0.75 * fx * fx * fx * fz

                     - 0.5 * fz * fga * pa_yy * fx

                     - pa_xxyy * fz * fgb

                     + 2.0 * pa_xx * fz * fx * fx

                     + 2.0 * fx * fx * pa_yy * fz

                     + 5.0 * pa_xxyy * fz * fx

                     - fx * fz * fga * pb_zz

                     - pa_xx * fz * fga * pb_zz

                     + 2.0 * fx * fx * fz * pb_zz

                     - fz * fga * pa_yy * pb_zz

                     + 5.0 * pa_xx * fz * fx * pb_zz

                     + 5.0 * fx * pa_yy * fz * pb_zz

                     + 12.0 * pa_xxyy * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xx_s_0(double fx,
                                   double pa_xxyz,
                                   double pa_xyz,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xx,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_yz

                     + 0.5 * pa_xxyz * fx

                     + 2.0 * pa_xyz * fx * pb_x

                     + 0.5 * fx * pa_yz * pb_xx

                     + pa_xxyz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xxyz,
                                   double pa_xyz,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xx,
                                   double r_0_0)
    {
        return r_0_0 * (6.0 * fx * fx * pa_yz * fz

                     - 0.5 * fx * pa_yz * fz * fgb

                     - 0.5 * fz * fga * pa_yz * fx

                     - pa_xxyz * fz * fgb

                     + 5.0 * pa_xxyz * fz * fx

                     + 20.0 * pa_xyz * fx * fz * pb_x

                     - fz * fga * pa_yz * pb_xx

                     + 5.0 * fx * pa_yz * fz * pb_xx

                     + 12.0 * pa_xxyz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xy_s_0(double fx,
                                   double pa_xxyz,
                                   double pa_xxz,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xz * fx * fx

                     + 0.25 * fx * fx * pa_z * pb_x

                     + 0.5 * pa_xxz * fx * pb_x

                     + pa_xyz * fx * pb_y

                     + 0.5 * fx * pa_yz * pb_xy

                     + pa_xxyz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xy_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_xxyz,
                                   double pa_xxz,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (4.0 * pa_xz * fx * fx * fz

                     - 0.5 * fz * fga * fx * pa_z * pb_x

                     + 2.0 * fx * fx * fz * pa_z * pb_x

                     + 5.0 * pa_xxz * fx * fz * pb_x

                     + 10.0 * pa_xyz * fx * fz * pb_y

                     - fz * fga * pa_yz * pb_xy

                     + 5.0 * fx * pa_yz * fz * pb_xy

                     + 12.0 * pa_xxyz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xz_s_0(double fx,
                                   double pa_xxy,
                                   double pa_xxyz,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xy * fx * fx

                     + 0.25 * fx * fx * pa_y * pb_x

                     + 0.5 * pa_xxy * fx * pb_x

                     + pa_xyz * fx * pb_z

                     + 0.5 * fx * pa_yz * pb_xz

                     + pa_xxyz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_xz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_xxy,
                                   double pa_xxyz,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (4.0 * pa_xy * fx * fx * fz

                     - 0.5 * fz * fga * pa_y * fx * pb_x

                     + 2.0 * fx * fx * pa_y * fz * pb_x

                     + 5.0 * pa_xxy * fz * fx * pb_x

                     + 10.0 * pa_xyz * fx * fz * pb_z

                     - fz * fga * pa_yz * pb_xz

                     + 5.0 * fx * pa_yz * fz * pb_xz

                     + 12.0 * pa_xxyz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yy_s_0(double fx,
                                   double pa_xxyz,
                                   double pa_xxz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_yz

                     + 0.5 * fx * fx * pa_z * pb_y

                     + 0.5 * pa_xxyz * fx

                     + pa_xxz * fx * pb_y

                     + 0.5 * fx * pa_yz * pb_yy

                     + pa_xxyz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xxyz,
                                   double pa_xxz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_yz * fz * fgb

                     - 0.5 * fz * fga * pa_yz * fx

                     - fz * fga * fx * pa_z * pb_y

                     - pa_xxyz * fz * fgb

                     + 2.0 * fx * fx * pa_yz * fz

                     + 4.0 * fx * fx * fz * pa_z * pb_y

                     + 5.0 * pa_xxyz * fz * fx

                     + 10.0 * pa_xxz * fx * fz * pb_y

                     - fz * fga * pa_yz * pb_yy

                     + 5.0 * fx * pa_yz * fz * pb_yy

                     + 12.0 * pa_xxyz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yz_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_xxyz,
                                   double pa_xxz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * pa_xx * fx * fx

                     + 0.25 * fx * fx * pa_y * pb_y

                     + 0.25 * fx * fx * pa_z * pb_z

                     + 0.5 * pa_xxy * fx * pb_y

                     + 0.5 * pa_xxz * fx * pb_z

                     + 0.5 * fx * pa_yz * pb_yz

                     + pa_xxyz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_yz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxy,
                                   double pa_xxyz,
                                   double pa_xxz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.25 * fz * fga * fx * fx

                     + 0.75 * fx * fx * fx * fz

                     + 2.0 * pa_xx * fx * fx * fz

                     - 0.5 * fz * fga * pa_y * fx * pb_y

                     - 0.5 * fz * fga * fx * pa_z * pb_z

                     + 2.0 * fx * fx * pa_y * fz * pb_y

                     + 2.0 * fx * fx * fz * pa_z * pb_z

                     + 5.0 * pa_xxy * fz * fx * pb_y

                     + 5.0 * pa_xxz * fx * fz * pb_z

                     - fz * fga * pa_yz * pb_yz

                     + 5.0 * fx * pa_yz * fz * pb_yz

                     + 12.0 * pa_xxyz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_zz_s_0(double fx,
                                   double pa_xxy,
                                   double pa_xxyz,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * fx * fx * pa_yz

                     + 0.5 * fx * fx * pa_y * pb_z

                     + 0.5 * pa_xxyz * fx

                     + pa_xxy * fx * pb_z

                     + 0.5 * fx * pa_yz * pb_zz

                     + pa_xxyz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxyz_zz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xxy,
                                   double pa_xxyz,
                                   double pa_y,
                                   double pa_yz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * fx * pa_yz * fz * fgb

                     - 0.5 * fz * fga * pa_yz * fx

                     - fz * fga * pa_y * fx * pb_z

                     - pa_xxyz * fz * fgb

                     + 2.0 * fx * fx * pa_yz * fz

                     + 4.0 * fx * fx * pa_y * fz * pb_z

                     + 5.0 * pa_xxyz * fz * fx

                     + 10.0 * pa_xxy * fz * fx * pb_z

                     - fz * fga * pa_yz * pb_zz

                     + 5.0 * fx * pa_yz * fz * pb_zz

                     + 12.0 * pa_xxyz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xx_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxzz,
                                   double pa_xzz,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xx,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * fx * fx * pa_zz

                     + 0.25 * pa_xx * fx * fx

                     + pa_x * fx * fx * pb_x

                     + 0.5 * pa_xxzz * fx

                     + 2.0 * pa_xzz * fx * pb_x

                     + 0.25 * fx * fx * pb_xx

                     + 0.5 * pa_xx * fx * pb_xx

                     + 0.5 * fx * pa_zz * pb_xx

                     + pa_xxzz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxzz,
                                   double pa_xzz,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xx,
                                   double r_0_0)
    {
        return r_0_0 * (- fx * fx * fz * fga

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * fx * fx * pa_zz * fz

                     - 0.25 * fx * fx * fz * fgb

                     - 0.5 * pa_xx * fx * fz * fgb

                     - 0.5 * pa_xx * fz * fga * fx

                     - 2.0 * pa_x * fx * fz * fga * pb_x

                     - 0.5 * fx * pa_zz * fz * fgb

                     - 0.5 * fz * fga * pa_zz * fx

                     - pa_xxzz * fz * fgb

                     + 2.0 * pa_xx * fz * fx * fx

                     + 8.0 * pa_x * fx * fx * fz * pb_x

                     + 5.0 * pa_xxzz * fz * fx

                     + 20.0 * pa_xzz * fx * fz * pb_x

                     - fx * fz * fga * pb_xx

                     - pa_xx * fz * fga * pb_xx

                     + 2.0 * fx * fx * fz * pb_xx

                     - fz * fga * pa_zz * pb_xx

                     + 5.0 * pa_xx * fz * fx * pb_xx

                     + 5.0 * fx * pa_zz * fz * pb_xx

                     + 12.0 * pa_xxzz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xy_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxzz,
                                   double pa_xzz,
                                   double pa_zz,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_x * fx * fx * pb_y

                     + pa_xzz * fx * pb_y

                     + 0.25 * fx * fx * pb_xy

                     + 0.5 * pa_xx * fx * pb_xy

                     + 0.5 * fx * pa_zz * pb_xy

                     + pa_xxzz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xy_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxzz,
                                   double pa_xzz,
                                   double pa_zz,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_x * fx * fz * fga * pb_y

                     + 4.0 * pa_x * fx * fx * fz * pb_y

                     + 10.0 * pa_xzz * fx * fz * pb_y

                     - fx * fz * fga * pb_xy

                     - pa_xx * fz * fga * pb_xy

                     + 2.0 * fx * fx * fz * pb_xy

                     - fz * fga * pa_zz * pb_xy

                     + 5.0 * pa_xx * fz * fx * pb_xy

                     + 5.0 * fx * pa_zz * fz * pb_xy

                     + 12.0 * pa_xxzz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xz_s_0(double fx,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_xxzz,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (pa_xz * fx * fx

                     + 0.5 * pa_x * fx * fx * pb_z

                     + 0.5 * fx * fx * pa_z * pb_x

                     + pa_xxz * fx * pb_x

                     + pa_xzz * fx * pb_z

                     + 0.25 * fx * fx * pb_xz

                     + 0.5 * pa_xx * fx * pb_xz

                     + 0.5 * fx * pa_zz * pb_xz

                     + pa_xxzz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_xz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_xxzz,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (8.0 * pa_xz * fx * fx * fz

                     - pa_x * fx * fz * fga * pb_z

                     - fz * fga * pa_z * fx * pb_x

                     + 4.0 * pa_x * fx * fx * fz * pb_z

                     + 4.0 * fx * fx * pa_z * fz * pb_x

                     + 10.0 * pa_xxz * fz * fx * pb_x

                     + 10.0 * pa_xzz * fx * fz * pb_z

                     - fx * fz * fga * pb_xz

                     - pa_xx * fz * fga * pb_xz

                     + 2.0 * fx * fx * fz * pb_xz

                     - fz * fga * pa_zz * pb_xz

                     + 5.0 * pa_xx * fz * fx * pb_xz

                     + 5.0 * fx * pa_zz * fz * pb_xz

                     + 12.0 * pa_xxzz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yy_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxzz,
                                   double pa_zz,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * pa_xx * fx * fx

                     + 0.25 * fx * fx * pa_zz

                     + 0.5 * pa_xxzz * fx

                     + 0.25 * fx * fx * pb_yy

                     + 0.5 * pa_xx * fx * pb_yy

                     + 0.5 * fx * pa_zz * pb_yy

                     + pa_xxzz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxzz,
                                   double pa_zz,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * fz * fgb

                     - 0.5 * fx * fx * fz * fga

                     - 0.5 * pa_xx * fx * fz * fgb

                     - 0.5 * pa_xx * fz * fga * fx

                     - 0.5 * fx * pa_zz * fz * fgb

                     + 0.75 * fx * fx * fx * fz

                     - 0.5 * fz * fga * pa_zz * fx

                     - pa_xxzz * fz * fgb

                     + 2.0 * pa_xx * fz * fx * fx

                     + 2.0 * fx * fx * pa_zz * fz

                     + 5.0 * pa_xxzz * fz * fx

                     - fx * fz * fga * pb_yy

                     - pa_xx * fz * fga * pb_yy

                     + 2.0 * fx * fx * fz * pb_yy

                     - fz * fga * pa_zz * pb_yy

                     + 5.0 * pa_xx * fz * fx * pb_yy

                     + 5.0 * fx * pa_zz * fz * pb_yy

                     + 12.0 * pa_xxzz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yz_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_xxzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yz,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * fx * fx * pa_z * pb_y

                     + pa_xxz * fx * pb_y

                     + 0.25 * fx * fx * pb_yz

                     + 0.5 * pa_xx * fx * pb_yz

                     + 0.5 * fx * pa_zz * pb_yz

                     + pa_xxzz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_yz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_xxzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yz,
                                   double r_0_0)
    {
        return r_0_0 * (- fz * fga * pa_z * fx * pb_y

                     + 4.0 * fx * fx * pa_z * fz * pb_y

                     + 10.0 * pa_xxz * fz * fx * pb_y

                     - fx * fz * fga * pb_yz

                     - pa_xx * fz * fga * pb_yz

                     + 2.0 * fx * fx * fz * pb_yz

                     - fz * fga * pa_zz * pb_yz

                     + 5.0 * pa_xx * fz * fx * pb_yz

                     + 5.0 * fx * pa_zz * fz * pb_yz

                     + 12.0 * pa_xxzz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_zz_s_0(double fx,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_xxzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_xx * fx * fx

                     + 0.25 * fx * fx * pa_zz

                     + fx * fx * pa_z * pb_z

                     + 0.5 * pa_xxzz * fx

                     + 2.0 * pa_xxz * fx * pb_z

                     + 0.25 * fx * fx * pb_zz

                     + 0.5 * pa_xx * fx * pb_zz

                     + 0.5 * fx * pa_zz * pb_zz

                     + pa_xxzz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xxzz_zz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xx,
                                   double pa_xxz,
                                   double pa_xxzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- fz * fga * fx * fx

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * pa_xx * fx * fx * fz

                     - 0.25 * fx * fx * fz * fgb

                     - 0.5 * pa_xx * fx * fz * fgb

                     - 0.5 * pa_xx * fz * fga * fx

                     - 0.5 * fx * pa_zz * fz * fgb

                     - 0.5 * fz * fga * pa_zz * fx

                     - 2.0 * fz * fga * pa_z * fx * pb_z

                     - pa_xxzz * fz * fgb

                     + 2.0 * fx * fx * pa_zz * fz

                     + 8.0 * fx * fx * pa_z * fz * pb_z

                     + 5.0 * pa_xxzz * fz * fx

                     + 20.0 * pa_xxz * fz * fx * pb_z

                     - fx * fz * fga * pb_zz

                     - pa_xx * fz * fga * pb_zz

                     + 2.0 * fx * fx * fz * pb_zz

                     - fz * fga * pa_zz * pb_zz

                     + 5.0 * pa_xx * fz * fx * pb_zz

                     + 5.0 * fx * pa_zz * fz * pb_zz

                     + 12.0 * pa_xxzz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xx_s_0(double fx,
                                   double pa_xy,
                                   double pa_xyyy,
                                   double pa_y,
                                   double pa_yyy,
                                   double pb_x,
                                   double pb_xx,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xy * fx * fx

                     + 1.5 * fx * fx * pa_y * pb_x

                     + 0.5 * pa_xyyy * fx

                     + fx * pa_yyy * pb_x

                     + 1.5 * pa_xy * fx * pb_xx

                     + pa_xyyy * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xy,
                                   double pa_xyyy,
                                   double pa_y,
                                   double pa_yyy,
                                   double pb_x,
                                   double pb_xx,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xy * fx * fz * fgb

                     - 1.5 * pa_xy * fz * fga * fx

                     - 3.0 * fx * pa_y * fz * fga * pb_x

                     - pa_xyyy * fz * fgb

                     + 6.0 * pa_xy * fz * fx * fx

                     + 12.0 * fx * fx * pa_y * fz * pb_x

                     + 5.0 * pa_xyyy * fz * fx

                     + 10.0 * fx * pa_yyy * fz * pb_x

                     - 3.0 * pa_xy * fz * fga * pb_xx

                     + 15.0 * pa_xy * fz * fx * pb_xx

                     + 12.0 * pa_xyyy * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xy_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_xyyy,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * fx * fx * pa_yy

                     + 0.75 * pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_y * pb_y

                     + 1.5 * pa_xyy * fx * pb_x

                     + 0.5 * fx * pa_yyy * pb_y

                     + 1.5 * pa_xy * fx * pb_xy

                     + pa_xyyy * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xy_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_xyyy,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fga

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * fx * fx * pa_yy * fz

                     - 1.5 * pa_x * fx * fz * fga * pb_x

                     - 1.5 * fx * pa_y * fz * fga * pb_y

                     + 6.0 * pa_x * fx * fx * fz * pb_x

                     + 6.0 * fx * fx * pa_y * fz * pb_y

                     + 15.0 * pa_xyy * fz * fx * pb_x

                     + 5.0 * fx * pa_yyy * fz * pb_y

                     - 3.0 * pa_xy * fz * fga * pb_xy

                     + 15.0 * pa_xy * fz * fx * pb_xy

                     + 12.0 * pa_xyyy * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xz_s_0(double fx,
                                   double pa_xy,
                                   double pa_xyyy,
                                   double pa_y,
                                   double pa_yyy,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_y * pb_z

                     + 0.5 * fx * pa_yyy * pb_z

                     + 1.5 * pa_xy * fx * pb_xz

                     + pa_xyyy * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_xz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_xy,
                                   double pa_xyyy,
                                   double pa_y,
                                   double pa_yyy,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_y * fz * fga * pb_z

                     + 6.0 * fx * fx * pa_y * fz * pb_z

                     + 5.0 * fx * pa_yyy * fz * pb_z

                     - 3.0 * pa_xy * fz * fga * pb_xz

                     + 15.0 * pa_xy * fz * fx * pb_xz

                     + 12.0 * pa_xyyy * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yy_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_xyyy,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_xy * fx * fx

                     + 1.5 * pa_x * fx * fx * pb_y

                     + 0.5 * pa_xyyy * fx

                     + 3.0 * pa_xyy * fx * pb_y

                     + 1.5 * pa_xy * fx * pb_yy

                     + pa_xyyy * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_xyyy,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (18.0 * pa_xy * fx * fx * fz

                     - 1.5 * pa_xy * fx * fz * fgb

                     - 1.5 * pa_xy * fz * fga * fx

                     - 3.0 * pa_x * fx * fz * fga * pb_y

                     - pa_xyyy * fz * fgb

                     + 12.0 * pa_x * fx * fx * fz * pb_y

                     + 5.0 * pa_xyyy * fz * fx

                     + 30.0 * pa_xyy * fz * fx * pb_y

                     - 3.0 * pa_xy * fz * fga * pb_yy

                     + 15.0 * pa_xy * fz * fx * pb_yy

                     + 12.0 * pa_xyyy * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yz_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_xyyy,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * pb_z

                     + 1.5 * pa_xyy * fx * pb_z

                     + 1.5 * pa_xy * fx * pb_yz

                     + pa_xyyy * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_yz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_xyyy,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fz * fga * pb_z

                     + 6.0 * pa_x * fx * fx * fz * pb_z

                     + 15.0 * pa_xyy * fz * fx * pb_z

                     - 3.0 * pa_xy * fz * fga * pb_yz

                     + 15.0 * pa_xy * fz * fx * pb_yz

                     + 12.0 * pa_xyyy * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_zz_s_0(double fx,
                                   double pa_xy,
                                   double pa_xyyy,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xy * fx * fx

                     + 0.5 * pa_xyyy * fx

                     + 1.5 * pa_xy * fx * pb_zz

                     + pa_xyyy * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyy_zz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xy,
                                   double pa_xyyy,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xy * fx * fz * fgb

                     - 1.5 * pa_xy * fz * fga * fx

                     - pa_xyyy * fz * fgb

                     + 6.0 * pa_xy * fz * fx * fx

                     + 5.0 * pa_xyyy * fz * fx

                     - 3.0 * pa_xy * fz * fga * pb_zz

                     + 15.0 * pa_xy * fz * fx * pb_zz

                     + 12.0 * pa_xyyy * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xx_s_0(double fx,
                                   double pa_xyyz,
                                   double pa_xz,
                                   double pa_yyz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xz * fx * fx

                     + 0.5 * fx * fx * pa_z * pb_x

                     + 0.5 * pa_xyyz * fx

                     + fx * pa_yyz * pb_x

                     + 0.5 * pa_xz * fx * pb_xx

                     + pa_xyyz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xyyz,
                                   double pa_xz,
                                   double pa_yyz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xx,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_xz * fx * fz * fgb

                     - 0.5 * pa_xz * fz * fga * fx

                     - fx * fz * fga * pa_z * pb_x

                     - pa_xyyz * fz * fgb

                     + 2.0 * pa_xz * fx * fx * fz

                     + 4.0 * fx * fx * fz * pa_z * pb_x

                     + 5.0 * pa_xyyz * fz * fx

                     + 10.0 * fx * pa_yyz * fz * pb_x

                     - pa_xz * fz * fga * pb_xx

                     + 5.0 * pa_xz * fx * fz * pb_xx

                     + 12.0 * pa_xyyz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xy_s_0(double fx,
                                   double pa_xyyz,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * fx * fx * pa_yz

                     + 0.25 * fx * fx * pa_z * pb_y

                     + pa_xyz * fx * pb_x

                     + 0.5 * fx * pa_yyz * pb_y

                     + 0.5 * pa_xz * fx * pb_xy

                     + pa_xyyz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xy_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_xyyz,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (4.0 * fx * fx * pa_yz * fz

                     - 0.5 * fx * fz * fga * pa_z * pb_y

                     + 2.0 * fx * fx * fz * pa_z * pb_y

                     + 10.0 * pa_xyz * fx * fz * pb_x

                     + 5.0 * fx * pa_yyz * fz * pb_y

                     - pa_xz * fz * fga * pb_xy

                     + 5.0 * pa_xz * fx * fz * pb_xy

                     + 12.0 * pa_xyyz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xz_s_0(double fx,
                                   double pa_x,
                                   double pa_xyy,
                                   double pa_xyyz,
                                   double pa_xz,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * fx * fx * pa_yy

                     + 0.25 * pa_x * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_z * pb_z

                     + 0.5 * pa_xyy * fx * pb_x

                     + 0.5 * fx * pa_yyz * pb_z

                     + 0.5 * pa_xz * fx * pb_xz

                     + pa_xyyz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_xz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xyy,
                                   double pa_xyyz,
                                   double pa_xz,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * fz * fga

                     + 0.75 * fx * fx * fx * fz

                     + 2.0 * fx * fx * pa_yy * fz

                     - 0.5 * pa_x * fz * fga * fx * pb_x

                     - 0.5 * fx * fz * fga * pa_z * pb_z

                     + 2.0 * pa_x * fx * fx * fz * pb_x

                     + 2.0 * fx * fx * fz * pa_z * pb_z

                     + 5.0 * pa_xyy * fz * fx * pb_x

                     + 5.0 * fx * pa_yyz * fz * pb_z

                     - pa_xz * fz * fga * pb_xz

                     + 5.0 * pa_xz * fx * fz * pb_xz

                     + 12.0 * pa_xyyz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yy_s_0(double fx,
                                   double pa_xyyz,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xz * fx * fx

                     + 0.5 * pa_xyyz * fx

                     + 2.0 * pa_xyz * fx * pb_y

                     + 0.5 * pa_xz * fx * pb_yy

                     + pa_xyyz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xyyz,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (6.0 * pa_xz * fx * fx * fz

                     - 0.5 * pa_xz * fx * fz * fgb

                     - 0.5 * pa_xz * fz * fga * fx

                     - pa_xyyz * fz * fgb

                     + 5.0 * pa_xyyz * fz * fx

                     + 20.0 * pa_xyz * fx * fz * pb_y

                     - pa_xz * fz * fga * pb_yy

                     + 5.0 * pa_xz * fx * fz * pb_yy

                     + 12.0 * pa_xyyz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yz_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_xyyz,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xy * fx * fx

                     + 0.25 * pa_x * fx * fx * pb_y

                     + 0.5 * pa_xyy * fx * pb_y

                     + pa_xyz * fx * pb_z

                     + 0.5 * pa_xz * fx * pb_yz

                     + pa_xyyz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_yz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyy,
                                   double pa_xyyz,
                                   double pa_xyz,
                                   double pa_xz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (4.0 * pa_xy * fx * fx * fz

                     - 0.5 * pa_x * fz * fga * fx * pb_y

                     + 2.0 * pa_x * fx * fx * fz * pb_y

                     + 5.0 * pa_xyy * fz * fx * pb_y

                     + 10.0 * pa_xyz * fx * fz * pb_z

                     - pa_xz * fz * fga * pb_yz

                     + 5.0 * pa_xz * fx * fz * pb_yz

                     + 12.0 * pa_xyyz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_zz_s_0(double fx,
                                   double pa_x,
                                   double pa_xyy,
                                   double pa_xyyz,
                                   double pa_xz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xz * fx * fx

                     + 0.5 * pa_x * fx * fx * pb_z

                     + 0.5 * pa_xyyz * fx

                     + pa_xyy * fx * pb_z

                     + 0.5 * pa_xz * fx * pb_zz

                     + pa_xyyz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyyz_zz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xyy,
                                   double pa_xyyz,
                                   double pa_xz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_xz * fx * fz * fgb

                     - 0.5 * pa_xz * fz * fga * fx

                     - pa_x * fz * fga * fx * pb_z

                     - pa_xyyz * fz * fgb

                     + 2.0 * pa_xz * fx * fx * fz

                     + 4.0 * pa_x * fx * fx * fz * pb_z

                     + 5.0 * pa_xyyz * fz * fx

                     + 10.0 * pa_xyy * fz * fx * pb_z

                     - pa_xz * fz * fga * pb_zz

                     + 5.0 * pa_xz * fx * fz * pb_zz

                     + 12.0 * pa_xyyz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xx_s_0(double fx,
                                   double pa_xy,
                                   double pa_xyzz,
                                   double pa_y,
                                   double pa_yzz,
                                   double pb_x,
                                   double pb_xx,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xy * fx * fx

                     + 0.5 * fx * fx * pa_y * pb_x

                     + 0.5 * pa_xyzz * fx

                     + fx * pa_yzz * pb_x

                     + 0.5 * pa_xy * fx * pb_xx

                     + pa_xyzz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xy,
                                   double pa_xyzz,
                                   double pa_y,
                                   double pa_yzz,
                                   double pb_x,
                                   double pb_xx,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_xy * fx * fz * fgb

                     - 0.5 * pa_xy * fz * fga * fx

                     - fx * pa_y * fz * fga * pb_x

                     - pa_xyzz * fz * fgb

                     + 2.0 * pa_xy * fz * fx * fx

                     + 4.0 * fx * fx * pa_y * fz * pb_x

                     + 5.0 * pa_xyzz * fz * fx

                     + 10.0 * fx * pa_yzz * fz * pb_x

                     - pa_xy * fz * fga * pb_xx

                     + 5.0 * pa_xy * fz * fx * pb_xx

                     + 12.0 * pa_xyzz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xy_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyzz,
                                   double pa_xzz,
                                   double pa_y,
                                   double pa_yzz,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * fx * fx * pa_zz

                     + 0.25 * pa_x * fx * fx * pb_x

                     + 0.25 * fx * fx * pa_y * pb_y

                     + 0.5 * pa_xzz * fx * pb_x

                     + 0.5 * fx * pa_yzz * pb_y

                     + 0.5 * pa_xy * fx * pb_xy

                     + pa_xyzz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xy_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyzz,
                                   double pa_xzz,
                                   double pa_y,
                                   double pa_yzz,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * fz * fga

                     + 0.75 * fx * fx * fx * fz

                     + 2.0 * fx * fx * pa_zz * fz

                     - 0.5 * pa_x * fx * fz * fga * pb_x

                     - 0.5 * fx * pa_y * fz * fga * pb_y

                     + 2.0 * pa_x * fx * fx * fz * pb_x

                     + 2.0 * fx * fx * pa_y * fz * pb_y

                     + 5.0 * pa_xzz * fx * fz * pb_x

                     + 5.0 * fx * pa_yzz * fz * pb_y

                     - pa_xy * fz * fga * pb_xy

                     + 5.0 * pa_xy * fz * fx * pb_xy

                     + 12.0 * pa_xyzz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xz_s_0(double fx,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_xyzz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * fx * fx * pa_yz

                     + 0.25 * fx * fx * pa_y * pb_z

                     + pa_xyz * fx * pb_x

                     + 0.5 * fx * pa_yzz * pb_z

                     + 0.5 * pa_xy * fx * pb_xz

                     + pa_xyzz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_xz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_xyzz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (4.0 * fx * fx * pa_yz * fz

                     - 0.5 * fx * pa_y * fz * fga * pb_z

                     + 2.0 * fx * fx * pa_y * fz * pb_z

                     + 10.0 * pa_xyz * fz * fx * pb_x

                     + 5.0 * fx * pa_yzz * fz * pb_z

                     - pa_xy * fz * fga * pb_xz

                     + 5.0 * pa_xy * fz * fx * pb_xz

                     + 12.0 * pa_xyzz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yy_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyzz,
                                   double pa_xzz,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.25 * pa_xy * fx * fx

                     + 0.5 * pa_x * fx * fx * pb_y

                     + 0.5 * pa_xyzz * fx

                     + pa_xzz * fx * pb_y

                     + 0.5 * pa_xy * fx * pb_yy

                     + pa_xyzz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyzz,
                                   double pa_xzz,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.5 * pa_xy * fx * fz * fgb

                     - 0.5 * pa_xy * fz * fga * fx

                     - pa_x * fx * fz * fga * pb_y

                     - pa_xyzz * fz * fgb

                     + 2.0 * pa_xy * fz * fx * fx

                     + 4.0 * pa_x * fx * fx * fz * pb_y

                     + 5.0 * pa_xyzz * fz * fx

                     + 10.0 * pa_xzz * fx * fz * pb_y

                     - pa_xy * fz * fga * pb_yy

                     + 5.0 * pa_xy * fz * fx * pb_yy

                     + 12.0 * pa_xyzz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yz_s_0(double fx,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_xyzz,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_xz * fx * fx

                     + 0.25 * pa_x * fx * fx * pb_z

                     + pa_xyz * fx * pb_y

                     + 0.5 * pa_xzz * fx * pb_z

                     + 0.5 * pa_xy * fx * pb_yz

                     + pa_xyzz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_yz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_xyzz,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (4.0 * pa_xz * fx * fx * fz

                     - 0.5 * pa_x * fx * fz * fga * pb_z

                     + 2.0 * pa_x * fx * fx * fz * pb_z

                     + 10.0 * pa_xyz * fz * fx * pb_y

                     + 5.0 * pa_xzz * fx * fz * pb_z

                     - pa_xy * fz * fga * pb_yz

                     + 5.0 * pa_xy * fz * fx * pb_yz

                     + 12.0 * pa_xyzz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_zz_s_0(double fx,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_xyzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xy * fx * fx

                     + 0.5 * pa_xyzz * fx

                     + 2.0 * pa_xyz * fx * pb_z

                     + 0.5 * pa_xy * fx * pb_zz

                     + pa_xyzz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xyzz_zz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xy,
                                   double pa_xyz,
                                   double pa_xyzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (6.0 * pa_xy * fx * fx * fz

                     - 0.5 * pa_xy * fx * fz * fgb

                     - 0.5 * pa_xy * fz * fga * fx

                     - pa_xyzz * fz * fgb

                     + 5.0 * pa_xyzz * fz * fx

                     + 20.0 * pa_xyz * fz * fx * pb_z

                     - pa_xy * fz * fga * pb_zz

                     + 5.0 * pa_xy * fz * fx * pb_zz

                     + 12.0 * pa_xyzz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xx_s_0(double fx,
                                   double pa_xz,
                                   double pa_xzzz,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_x,
                                   double pb_xx,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xz * fx * fx

                     + 1.5 * fx * fx * pa_z * pb_x

                     + 0.5 * pa_xzzz * fx

                     + fx * pa_zzz * pb_x

                     + 1.5 * pa_xz * fx * pb_xx

                     + pa_xzzz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xz,
                                   double pa_xzzz,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_x,
                                   double pb_xx,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xz * fx * fz * fgb

                     - 1.5 * pa_xz * fz * fga * fx

                     - 3.0 * fx * pa_z * fz * fga * pb_x

                     - pa_xzzz * fz * fgb

                     + 6.0 * pa_xz * fz * fx * fx

                     + 12.0 * fx * fx * pa_z * fz * pb_x

                     + 5.0 * pa_xzzz * fz * fx

                     + 10.0 * fx * pa_zzz * fz * pb_x

                     - 3.0 * pa_xz * fz * fga * pb_xx

                     + 15.0 * pa_xz * fz * fx * pb_xx

                     + 12.0 * pa_xzzz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xy_s_0(double fx,
                                   double pa_xz,
                                   double pa_xzzz,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_xy,
                                   double pb_y,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z * pb_y

                     + 0.5 * fx * pa_zzz * pb_y

                     + 1.5 * pa_xz * fx * pb_xy

                     + pa_xzzz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xy_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_xz,
                                   double pa_xzzz,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_xy,
                                   double pb_y,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_z * fz * fga * pb_y

                     + 6.0 * fx * fx * pa_z * fz * pb_y

                     + 5.0 * fx * pa_zzz * fz * pb_y

                     - 3.0 * pa_xz * fz * fga * pb_xy

                     + 15.0 * pa_xz * fz * fx * pb_xy

                     + 12.0 * pa_xzzz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pa_xzzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * fx * fx * pa_zz

                     + 0.75 * pa_x * fx * fx * pb_x

                     + 0.75 * fx * fx * pa_z * pb_z

                     + 1.5 * pa_xzz * fx * pb_x

                     + 0.5 * fx * pa_zzz * pb_z

                     + 1.5 * pa_xz * fx * pb_xz

                     + pa_xzzz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_xz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pa_xzzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_x,
                                   double pb_xz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fga

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * fx * fx * pa_zz * fz

                     - 1.5 * pa_x * fx * fz * fga * pb_x

                     - 1.5 * fx * pa_z * fz * fga * pb_z

                     + 6.0 * pa_x * fx * fx * fz * pb_x

                     + 6.0 * fx * fx * pa_z * fz * pb_z

                     + 15.0 * pa_xzz * fz * fx * pb_x

                     + 5.0 * fx * pa_zzz * fz * pb_z

                     - 3.0 * pa_xz * fz * fga * pb_xz

                     + 15.0 * pa_xz * fz * fx * pb_xz

                     + 12.0 * pa_xzzz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yy_s_0(double fx,
                                   double pa_xz,
                                   double pa_xzzz,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_xz * fx * fx

                     + 0.5 * pa_xzzz * fx

                     + 1.5 * pa_xz * fx * pb_yy

                     + pa_xzzz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_xz,
                                   double pa_xzzz,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_xz * fx * fz * fgb

                     - 1.5 * pa_xz * fz * fga * fx

                     - pa_xzzz * fz * fgb

                     + 6.0 * pa_xz * fz * fx * fx

                     + 5.0 * pa_xzzz * fz * fx

                     - 3.0 * pa_xz * fz * fga * pb_yy

                     + 15.0 * pa_xz * fz * fx * pb_yy

                     + 12.0 * pa_xzzz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pa_xzzz,
                                   double pb_y,
                                   double pb_yz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_x * fx * fx * pb_y

                     + 1.5 * pa_xzz * fx * pb_y

                     + 1.5 * pa_xz * fx * pb_yz

                     + pa_xzzz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_yz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pa_xzzz,
                                   double pb_y,
                                   double pb_yz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_x * fx * fz * fga * pb_y

                     + 6.0 * pa_x * fx * fx * fz * pb_y

                     + 15.0 * pa_xzz * fz * fx * pb_y

                     - 3.0 * pa_xz * fz * fga * pb_yz

                     + 15.0 * pa_xz * fz * fx * pb_yz

                     + 12.0 * pa_xzzz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_zz_s_0(double fx,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pa_xzzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_xz * fx * fx

                     + 1.5 * pa_x * fx * fx * pb_z

                     + 0.5 * pa_xzzz * fx

                     + 3.0 * pa_xzz * fx * pb_z

                     + 1.5 * pa_xz * fx * pb_zz

                     + pa_xzzz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_xzzz_zz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_x,
                                   double pa_xz,
                                   double pa_xzz,
                                   double pa_xzzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (18.0 * pa_xz * fx * fx * fz

                     - 1.5 * pa_xz * fx * fz * fgb

                     - 1.5 * pa_xz * fz * fga * fx

                     - 3.0 * pa_x * fx * fz * fga * pb_z

                     - pa_xzzz * fz * fgb

                     + 12.0 * pa_x * fx * fx * fz * pb_z

                     + 5.0 * pa_xzzz * fz * fx

                     + 30.0 * pa_xzz * fz * fx * pb_z

                     - 3.0 * pa_xz * fz * fga * pb_zz

                     + 15.0 * pa_xz * fz * fx * pb_zz

                     + 12.0 * pa_xzzz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xx_s_0(double fx,
                                   double pa_yy,
                                   double pa_yyyy,
                                   double pb_xx,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 1.5 * pa_yy * fx * fx

                     + 0.5 * pa_yyyy * fx

                     + 0.75 * fx * fx * pb_xx

                     + 3.0 * pa_yy * fx * pb_xx

                     + pa_yyyy * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yy,
                                   double pa_yyyy,
                                   double pb_xx,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 1.5 * fx * fx * fz * fga

                     - 3.0 * pa_yy * fx * fz * fgb

                     - 3.0 * pa_yy * fz * fga * fx

                     + 2.25 * fx * fx * fx * fz

                     - pa_yyyy * fz * fgb

                     + 12.0 * pa_yy * fz * fx * fx

                     + 5.0 * pa_yyyy * fz * fx

                     - 3.0 * fx * fz * fga * pb_xx

                     - 6.0 * pa_yy * fz * fga * pb_xx

                     + 6.0 * fx * fx * fz * pb_xx

                     + 30.0 * pa_yy * fz * fx * pb_xx

                     + 12.0 * pa_yyyy * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xy_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pa_yyyy,
                                   double pb_x,
                                   double pb_xy,
                                   double s_0_0)
    {
        return s_0_0 * (3.0 * pa_y * fx * fx * pb_x

                     + 2.0 * pa_yyy * fx * pb_x

                     + 0.75 * fx * fx * pb_xy

                     + 3.0 * pa_yy * fx * pb_xy

                     + pa_yyyy * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xy_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pa_yyyy,
                                   double pb_x,
                                   double pb_xy,
                                   double r_0_0)
    {
        return r_0_0 * (- 6.0 * pa_y * fx * fz * fga * pb_x

                     + 24.0 * pa_y * fx * fx * fz * pb_x

                     + 20.0 * pa_yyy * fz * fx * pb_x

                     - 3.0 * fx * fz * fga * pb_xy

                     - 6.0 * pa_yy * fz * fga * pb_xy

                     + 6.0 * fx * fx * fz * pb_xy

                     + 30.0 * pa_yy * fz * fx * pb_xy

                     + 12.0 * pa_yyyy * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xz_s_0(double fx,
                                   double pa_yy,
                                   double pa_yyyy,
                                   double pb_xz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_xz

                     + 3.0 * pa_yy * fx * pb_xz

                     + pa_yyyy * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_xz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_yy,
                                   double pa_yyyy,
                                   double pb_xz,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * fx * fz * fga * pb_xz

                     - 6.0 * pa_yy * fz * fga * pb_xz

                     + 6.0 * fx * fx * fz * pb_xz

                     + 30.0 * pa_yy * fz * fx * pb_xz

                     + 12.0 * pa_yyyy * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yy_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pa_yyyy,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx

                     + 4.5 * pa_yy * fx * fx

                     + 6.0 * pa_y * fx * fx * pb_y

                     + 0.5 * pa_yyyy * fx

                     + 4.0 * pa_yyy * fx * pb_y

                     + 0.75 * fx * fx * pb_yy

                     + 3.0 * pa_yy * fx * pb_yy

                     + pa_yyyy * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pa_yyyy,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * fz * fga

                     + 11.25 * fx * fx * fx * fz

                     + 36.0 * pa_yy * fx * fx * fz

                     - 0.75 * fx * fx * fz * fgb

                     - 3.0 * pa_yy * fx * fz * fgb

                     - 3.0 * pa_yy * fz * fga * fx

                     - 12.0 * pa_y * fx * fz * fga * pb_y

                     - pa_yyyy * fz * fgb

                     + 48.0 * pa_y * fx * fx * fz * pb_y

                     + 5.0 * pa_yyyy * fz * fx

                     + 40.0 * pa_yyy * fz * fx * pb_y

                     - 3.0 * fx * fz * fga * pb_yy

                     - 6.0 * pa_yy * fz * fga * pb_yy

                     + 6.0 * fx * fx * fz * pb_yy

                     + 30.0 * pa_yy * fz * fx * pb_yy

                     + 12.0 * pa_yyyy * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yz_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pa_yyyy,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (3.0 * pa_y * fx * fx * pb_z

                     + 2.0 * pa_yyy * fx * pb_z

                     + 0.75 * fx * fx * pb_yz

                     + 3.0 * pa_yy * fx * pb_yz

                     + pa_yyyy * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_yz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pa_yyyy,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 6.0 * pa_y * fx * fz * fga * pb_z

                     + 24.0 * pa_y * fx * fx * fz * pb_z

                     + 20.0 * pa_yyy * fz * fx * pb_z

                     - 3.0 * fx * fz * fga * pb_yz

                     - 6.0 * pa_yy * fz * fga * pb_yz

                     + 6.0 * fx * fx * fz * pb_yz

                     + 30.0 * pa_yy * fz * fx * pb_yz

                     + 12.0 * pa_yyyy * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_zz_s_0(double fx,
                                   double pa_yy,
                                   double pa_yyyy,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 1.5 * pa_yy * fx * fx

                     + 0.5 * pa_yyyy * fx

                     + 0.75 * fx * fx * pb_zz

                     + 3.0 * pa_yy * fx * pb_zz

                     + pa_yyyy * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyy_zz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yy,
                                   double pa_yyyy,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 1.5 * fx * fx * fz * fga

                     - 3.0 * pa_yy * fx * fz * fgb

                     - 3.0 * pa_yy * fz * fga * fx

                     + 2.25 * fx * fx * fx * fz

                     - pa_yyyy * fz * fgb

                     + 12.0 * pa_yy * fz * fx * fx

                     + 5.0 * pa_yyyy * fz * fx

                     - 3.0 * fx * fz * fga * pb_zz

                     - 6.0 * pa_yy * fz * fga * pb_zz

                     + 6.0 * fx * fx * fz * pb_zz

                     + 30.0 * pa_yy * fz * fx * pb_zz

                     + 12.0 * pa_yyyy * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xx_s_0(double fx,
                                   double pa_yyyz,
                                   double pa_yz,
                                   double pb_xx,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_yz * fx * fx

                     + 0.5 * pa_yyyz * fx

                     + 1.5 * pa_yz * fx * pb_xx

                     + pa_yyyz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yyyz,
                                   double pa_yz,
                                   double pb_xx,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_yz * fx * fz * fgb

                     - 1.5 * pa_yz * fz * fga * fx

                     - pa_yyyz * fz * fgb

                     + 6.0 * pa_yz * fx * fx * fz

                     + 5.0 * pa_yyyz * fz * fx

                     - 3.0 * pa_yz * fz * fga * pb_xx

                     + 15.0 * pa_yz * fx * fz * pb_xx

                     + 12.0 * pa_yyyz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xy_s_0(double fx,
                                   double pa_yyyz,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z * pb_x

                     + 1.5 * pa_yyz * fx * pb_x

                     + 1.5 * pa_yz * fx * pb_xy

                     + pa_yyyz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xy_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_yyyz,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_x,
                                   double pb_xy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * fz * fga * pa_z * pb_x

                     + 6.0 * fx * fx * fz * pa_z * pb_x

                     + 15.0 * pa_yyz * fx * fz * pb_x

                     - 3.0 * pa_yz * fz * fga * pb_xy

                     + 15.0 * pa_yz * fx * fz * pb_xy

                     + 12.0 * pa_yyyz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xz_s_0(double fx,
                                   double pa_y,
                                   double pa_yyy,
                                   double pa_yyyz,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx * pb_x

                     + 0.5 * pa_yyy * fx * pb_x

                     + 1.5 * pa_yz * fx * pb_xz

                     + pa_yyyz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_xz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yyy,
                                   double pa_yyyz,
                                   double pa_yz,
                                   double pb_x,
                                   double pb_xz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fz * fga * fx * pb_x

                     + 6.0 * pa_y * fx * fx * fz * pb_x

                     + 5.0 * pa_yyy * fz * fx * pb_x

                     - 3.0 * pa_yz * fz * fga * pb_xz

                     + 15.0 * pa_yz * fx * fz * pb_xz

                     + 12.0 * pa_yyyz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yy_s_0(double fx,
                                   double pa_yyyz,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_yz * fx * fx

                     + 1.5 * fx * fx * pa_z * pb_y

                     + 0.5 * pa_yyyz * fx

                     + 3.0 * pa_yyz * fx * pb_y

                     + 1.5 * pa_yz * fx * pb_yy

                     + pa_yyyz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yyyz,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (18.0 * pa_yz * fx * fx * fz

                     - 1.5 * pa_yz * fx * fz * fgb

                     - 1.5 * pa_yz * fz * fga * fx

                     - 3.0 * fx * fz * fga * pa_z * pb_y

                     - pa_yyyz * fz * fgb

                     + 12.0 * fx * fx * fz * pa_z * pb_y

                     + 5.0 * pa_yyyz * fz * fx

                     + 30.0 * pa_yyz * fx * fz * pb_y

                     - 3.0 * pa_yz * fz * fga * pb_yy

                     + 15.0 * pa_yz * fx * fz * pb_yy

                     + 12.0 * pa_yyyz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yz_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pa_yyyz,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_yy * fx * fx

                     + 0.75 * pa_y * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_z * pb_z

                     + 0.5 * pa_yyy * fx * pb_y

                     + 1.5 * pa_yyz * fx * pb_z

                     + 1.5 * pa_yz * fx * pb_yz

                     + pa_yyyz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_yz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyy,
                                   double pa_yyyz,
                                   double pa_yyz,
                                   double pa_yz,
                                   double pa_z,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fga

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * pa_yy * fx * fx * fz

                     - 1.5 * pa_y * fz * fga * fx * pb_y

                     - 1.5 * fx * fz * fga * pa_z * pb_z

                     + 6.0 * pa_y * fx * fx * fz * pb_y

                     + 6.0 * fx * fx * fz * pa_z * pb_z

                     + 5.0 * pa_yyy * fz * fx * pb_y

                     + 15.0 * pa_yyz * fx * fz * pb_z

                     - 3.0 * pa_yz * fz * fga * pb_yz

                     + 15.0 * pa_yz * fx * fz * pb_yz

                     + 12.0 * pa_yyyz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_zz_s_0(double fx,
                                   double pa_y,
                                   double pa_yyy,
                                   double pa_yyyz,
                                   double pa_yz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_yz * fx * fx

                     + 1.5 * pa_y * fx * fx * pb_z

                     + 0.5 * pa_yyyz * fx

                     + pa_yyy * fx * pb_z

                     + 1.5 * pa_yz * fx * pb_zz

                     + pa_yyyz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyyz_zz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yyy,
                                   double pa_yyyz,
                                   double pa_yz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_yz * fx * fz * fgb

                     - 1.5 * pa_yz * fz * fga * fx

                     - 3.0 * pa_y * fz * fga * fx * pb_z

                     - pa_yyyz * fz * fgb

                     + 6.0 * pa_yz * fx * fx * fz

                     + 12.0 * pa_y * fx * fx * fz * pb_z

                     + 5.0 * pa_yyyz * fz * fx

                     + 10.0 * pa_yyy * fz * fx * pb_z

                     - 3.0 * pa_yz * fz * fga * pb_zz

                     + 15.0 * pa_yz * fx * fz * pb_zz

                     + 12.0 * pa_yyyz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xx_s_0(double fx,
                                   double pa_yy,
                                   double pa_yyzz,
                                   double pa_zz,
                                   double pb_xx,
                                   double s_0_0)
    {
        return s_0_0 * (0.125 * fx * fx * fx

                     + 0.25 * pa_yy * fx * fx

                     + 0.25 * fx * fx * pa_zz

                     + 0.5 * pa_yyzz * fx

                     + 0.25 * fx * fx * pb_xx

                     + 0.5 * pa_yy * fx * pb_xx

                     + 0.5 * fx * pa_zz * pb_xx

                     + pa_yyzz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yy,
                                   double pa_yyzz,
                                   double pa_zz,
                                   double pb_xx,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.25 * fx * fx * fz * fgb

                     - 0.5 * fx * fx * fz * fga

                     - 0.5 * pa_yy * fx * fz * fgb

                     - 0.5 * pa_yy * fz * fga * fx

                     - 0.5 * fx * pa_zz * fz * fgb

                     + 0.75 * fx * fx * fx * fz

                     - 0.5 * fz * fga * pa_zz * fx

                     - pa_yyzz * fz * fgb

                     + 2.0 * pa_yy * fz * fx * fx

                     + 2.0 * fx * fx * pa_zz * fz

                     + 5.0 * pa_yyzz * fz * fx

                     - fx * fz * fga * pb_xx

                     - pa_yy * fz * fga * pb_xx

                     + 2.0 * fx * fx * fz * pb_xx

                     - fz * fga * pa_zz * pb_xx

                     + 5.0 * pa_yy * fz * fx * pb_xx

                     + 5.0 * fx * pa_zz * fz * pb_xx

                     + 12.0 * pa_yyzz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xy_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyzz,
                                   double pa_yzz,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xy,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * pa_y * fx * fx * pb_x

                     + pa_yzz * fx * pb_x

                     + 0.25 * fx * fx * pb_xy

                     + 0.5 * pa_yy * fx * pb_xy

                     + 0.5 * fx * pa_zz * pb_xy

                     + pa_yyzz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xy_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyzz,
                                   double pa_yzz,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xy,
                                   double r_0_0)
    {
        return r_0_0 * (- pa_y * fx * fz * fga * pb_x

                     + 4.0 * pa_y * fx * fx * fz * pb_x

                     + 10.0 * pa_yzz * fx * fz * pb_x

                     - fx * fz * fga * pb_xy

                     - pa_yy * fz * fga * pb_xy

                     + 2.0 * fx * fx * fz * pb_xy

                     - fz * fga * pa_zz * pb_xy

                     + 5.0 * pa_yy * fz * fx * pb_xy

                     + 5.0 * fx * pa_zz * fz * pb_xy

                     + 12.0 * pa_yyzz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xz_s_0(double fx,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_yyzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xz,
                                   double s_0_0)
    {
        return s_0_0 * (0.5 * fx * fx * pa_z * pb_x

                     + pa_yyz * fx * pb_x

                     + 0.25 * fx * fx * pb_xz

                     + 0.5 * pa_yy * fx * pb_xz

                     + 0.5 * fx * pa_zz * pb_xz

                     + pa_yyzz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_xz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_yyzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_x,
                                   double pb_xz,
                                   double r_0_0)
    {
        return r_0_0 * (- fz * fga * pa_z * fx * pb_x

                     + 4.0 * fx * fx * pa_z * fz * pb_x

                     + 10.0 * pa_yyz * fz * fx * pb_x

                     - fx * fz * fga * pb_xz

                     - pa_yy * fz * fga * pb_xz

                     + 2.0 * fx * fx * fz * pb_xz

                     - fz * fga * pa_zz * pb_xz

                     + 5.0 * pa_yy * fz * fx * pb_xz

                     + 5.0 * fx * pa_zz * fz * pb_xz

                     + 12.0 * pa_yyzz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yy_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyzz,
                                   double pa_yzz,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * fx * fx * pa_zz

                     + 0.25 * pa_yy * fx * fx

                     + pa_y * fx * fx * pb_y

                     + 0.5 * pa_yyzz * fx

                     + 2.0 * pa_yzz * fx * pb_y

                     + 0.25 * fx * fx * pb_yy

                     + 0.5 * pa_yy * fx * pb_yy

                     + 0.5 * fx * pa_zz * pb_yy

                     + pa_yyzz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyzz,
                                   double pa_yzz,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- fx * fx * fz * fga

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * fx * fx * pa_zz * fz

                     - 0.25 * fx * fx * fz * fgb

                     - 0.5 * pa_yy * fx * fz * fgb

                     - 0.5 * pa_yy * fz * fga * fx

                     - 2.0 * pa_y * fx * fz * fga * pb_y

                     - 0.5 * fx * pa_zz * fz * fgb

                     - 0.5 * fz * fga * pa_zz * fx

                     - pa_yyzz * fz * fgb

                     + 2.0 * pa_yy * fz * fx * fx

                     + 8.0 * pa_y * fx * fx * fz * pb_y

                     + 5.0 * pa_yyzz * fz * fx

                     + 20.0 * pa_yzz * fx * fz * pb_y

                     - fx * fz * fga * pb_yy

                     - pa_yy * fz * fga * pb_yy

                     + 2.0 * fx * fx * fz * pb_yy

                     - fz * fga * pa_zz * pb_yy

                     + 5.0 * pa_yy * fz * fx * pb_yy

                     + 5.0 * fx * pa_zz * fz * pb_yy

                     + 12.0 * pa_yyzz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yz_s_0(double fx,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_yyzz,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (pa_yz * fx * fx

                     + 0.5 * pa_y * fx * fx * pb_z

                     + 0.5 * fx * fx * pa_z * pb_y

                     + pa_yyz * fx * pb_y

                     + pa_yzz * fx * pb_z

                     + 0.25 * fx * fx * pb_yz

                     + 0.5 * pa_yy * fx * pb_yz

                     + 0.5 * fx * pa_zz * pb_yz

                     + pa_yyzz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_yz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_yyzz,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (8.0 * pa_yz * fx * fx * fz

                     - pa_y * fx * fz * fga * pb_z

                     - fz * fga * pa_z * fx * pb_y

                     + 4.0 * pa_y * fx * fx * fz * pb_z

                     + 4.0 * fx * fx * pa_z * fz * pb_y

                     + 10.0 * pa_yyz * fz * fx * pb_y

                     + 10.0 * pa_yzz * fx * fz * pb_z

                     - fx * fz * fga * pb_yz

                     - pa_yy * fz * fga * pb_yz

                     + 2.0 * fx * fx * fz * pb_yz

                     - fz * fga * pa_zz * pb_yz

                     + 5.0 * pa_yy * fz * fx * pb_yz

                     + 5.0 * fx * pa_zz * fz * pb_yz

                     + 12.0 * pa_yyzz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_zz_s_0(double fx,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_yyzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * pa_yy * fx * fx

                     + 0.25 * fx * fx * pa_zz

                     + fx * fx * pa_z * pb_z

                     + 0.5 * pa_yyzz * fx

                     + 2.0 * pa_yyz * fx * pb_z

                     + 0.25 * fx * fx * pb_zz

                     + 0.5 * pa_yy * fx * pb_zz

                     + 0.5 * fx * pa_zz * pb_zz

                     + pa_yyzz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yyzz_zz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yy,
                                   double pa_yyz,
                                   double pa_yyzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- fz * fga * fx * fx

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * pa_yy * fx * fx * fz

                     - 0.25 * fx * fx * fz * fgb

                     - 0.5 * pa_yy * fx * fz * fgb

                     - 0.5 * pa_yy * fz * fga * fx

                     - 0.5 * fx * pa_zz * fz * fgb

                     - 0.5 * fz * fga * pa_zz * fx

                     - 2.0 * fz * fga * pa_z * fx * pb_z

                     - pa_yyzz * fz * fgb

                     + 2.0 * fx * fx * pa_zz * fz

                     + 8.0 * fx * fx * pa_z * fz * pb_z

                     + 5.0 * pa_yyzz * fz * fx

                     + 20.0 * pa_yyz * fz * fx * pb_z

                     - fx * fz * fga * pb_zz

                     - pa_yy * fz * fga * pb_zz

                     + 2.0 * fx * fx * fz * pb_zz

                     - fz * fga * pa_zz * pb_zz

                     + 5.0 * pa_yy * fz * fx * pb_zz

                     + 5.0 * fx * pa_zz * fz * pb_zz

                     + 12.0 * pa_yyzz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xx_s_0(double fx,
                                   double pa_yz,
                                   double pa_yzzz,
                                   double pb_xx,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_yz * fx * fx

                     + 0.5 * pa_yzzz * fx

                     + 1.5 * pa_yz * fx * pb_xx

                     + pa_yzzz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yz,
                                   double pa_yzzz,
                                   double pb_xx,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_yz * fx * fz * fgb

                     - 1.5 * pa_yz * fz * fga * fx

                     - pa_yzzz * fz * fgb

                     + 6.0 * pa_yz * fz * fx * fx

                     + 5.0 * pa_yzzz * fz * fx

                     - 3.0 * pa_yz * fz * fga * pb_xx

                     + 15.0 * pa_yz * fz * fx * pb_xx

                     + 12.0 * pa_yzzz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xy_s_0(double fx,
                                   double pa_yz,
                                   double pa_yzzz,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_x,
                                   double pb_xy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pa_z * pb_x

                     + 0.5 * fx * pa_zzz * pb_x

                     + 1.5 * pa_yz * fx * pb_xy

                     + pa_yzzz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xy_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_yz,
                                   double pa_yzzz,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_x,
                                   double pb_xy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * fx * pa_z * fz * fga * pb_x

                     + 6.0 * fx * fx * pa_z * fz * pb_x

                     + 5.0 * fx * pa_zzz * fz * pb_x

                     - 3.0 * pa_yz * fz * fga * pb_xy

                     + 15.0 * pa_yz * fz * fx * pb_xy

                     + 12.0 * pa_yzzz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pa_yzzz,
                                   double pb_x,
                                   double pb_xz,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_y * fx * fx * pb_x

                     + 1.5 * pa_yzz * fx * pb_x

                     + 1.5 * pa_yz * fx * pb_xz

                     + pa_yzzz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_xz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pa_yzzz,
                                   double pb_x,
                                   double pb_xz,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_y * fx * fz * fga * pb_x

                     + 6.0 * pa_y * fx * fx * fz * pb_x

                     + 15.0 * pa_yzz * fz * fx * pb_x

                     - 3.0 * pa_yz * fz * fga * pb_xz

                     + 15.0 * pa_yz * fz * fx * pb_xz

                     + 12.0 * pa_yzzz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yy_s_0(double fx,
                                   double pa_yz,
                                   double pa_yzzz,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_y,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * pa_yz * fx * fx

                     + 1.5 * fx * fx * pa_z * pb_y

                     + 0.5 * pa_yzzz * fx

                     + fx * pa_zzz * pb_y

                     + 1.5 * pa_yz * fx * pb_yy

                     + pa_yzzz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_yz,
                                   double pa_yzzz,
                                   double pa_z,
                                   double pa_zzz,
                                   double pb_y,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- 1.5 * pa_yz * fx * fz * fgb

                     - 1.5 * pa_yz * fz * fga * fx

                     - 3.0 * fx * pa_z * fz * fga * pb_y

                     - pa_yzzz * fz * fgb

                     + 6.0 * pa_yz * fz * fx * fx

                     + 12.0 * fx * fx * pa_z * fz * pb_y

                     + 5.0 * pa_yzzz * fz * fx

                     + 10.0 * fx * pa_zzz * fz * pb_y

                     - 3.0 * pa_yz * fz * fga * pb_yy

                     + 15.0 * pa_yz * fz * fx * pb_yy

                     + 12.0 * pa_yzzz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pa_yzzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 0.75 * fx * fx * pa_zz

                     + 0.75 * pa_y * fx * fx * pb_y

                     + 0.75 * fx * fx * pa_z * pb_z

                     + 1.5 * pa_yzz * fx * pb_y

                     + 0.5 * fx * pa_zzz * pb_z

                     + 1.5 * pa_yz * fx * pb_yz

                     + pa_yzzz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_yz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pa_yzzz,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pb_y,
                                   double pb_yz,
                                   double pb_z,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fga

                     + 2.25 * fx * fx * fx * fz

                     + 6.0 * fx * fx * pa_zz * fz

                     - 1.5 * pa_y * fx * fz * fga * pb_y

                     - 1.5 * fx * pa_z * fz * fga * pb_z

                     + 6.0 * pa_y * fx * fx * fz * pb_y

                     + 6.0 * fx * fx * pa_z * fz * pb_z

                     + 15.0 * pa_yzz * fz * fx * pb_y

                     + 5.0 * fx * pa_zzz * fz * pb_z

                     - 3.0 * pa_yz * fz * fga * pb_yz

                     + 15.0 * pa_yz * fz * fx * pb_yz

                     + 12.0 * pa_yzzz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_zz_s_0(double fx,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pa_yzzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (2.25 * pa_yz * fx * fx

                     + 1.5 * pa_y * fx * fx * pb_z

                     + 0.5 * pa_yzzz * fx

                     + 3.0 * pa_yzz * fx * pb_z

                     + 1.5 * pa_yz * fx * pb_zz

                     + pa_yzzz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_yzzz_zz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_y,
                                   double pa_yz,
                                   double pa_yzz,
                                   double pa_yzzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (18.0 * pa_yz * fx * fx * fz

                     - 1.5 * pa_yz * fx * fz * fgb

                     - 1.5 * pa_yz * fz * fga * fx

                     - 3.0 * pa_y * fx * fz * fga * pb_z

                     - pa_yzzz * fz * fgb

                     + 12.0 * pa_y * fx * fx * fz * pb_z

                     + 5.0 * pa_yzzz * fz * fx

                     + 30.0 * pa_yzz * fz * fx * pb_z

                     - 3.0 * pa_yz * fz * fga * pb_zz

                     + 15.0 * pa_yz * fz * fx * pb_zz

                     + 12.0 * pa_yzzz * fz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xx_s_0(double fx,
                                   double pa_zz,
                                   double pa_zzzz,
                                   double pb_xx,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 1.5 * pa_zz * fx * fx

                     + 0.5 * pa_zzzz * fx

                     + 0.75 * fx * fx * pb_xx

                     + 3.0 * pa_zz * fx * pb_xx

                     + pa_zzzz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xx_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_zz,
                                   double pa_zzzz,
                                   double pb_xx,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 1.5 * fx * fx * fz * fga

                     - 3.0 * pa_zz * fx * fz * fgb

                     - 3.0 * pa_zz * fz * fga * fx

                     + 2.25 * fx * fx * fx * fz

                     - pa_zzzz * fz * fgb

                     + 12.0 * pa_zz * fz * fx * fx

                     + 5.0 * pa_zzzz * fz * fx

                     - 3.0 * fx * fz * fga * pb_xx

                     - 6.0 * pa_zz * fz * fga * pb_xx

                     + 6.0 * fx * fx * fz * pb_xx

                     + 30.0 * pa_zz * fz * fx * pb_xx

                     + 12.0 * pa_zzzz * fz * pb_xx);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xy_s_0(double fx,
                                   double pa_zz,
                                   double pa_zzzz,
                                   double pb_xy,
                                   double s_0_0)
    {
        return s_0_0 * (0.75 * fx * fx * pb_xy

                     + 3.0 * pa_zz * fx * pb_xy

                     + pa_zzzz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xy_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_zz,
                                   double pa_zzzz,
                                   double pb_xy,
                                   double r_0_0)
    {
        return r_0_0 * (- 3.0 * fx * fz * fga * pb_xy

                     - 6.0 * pa_zz * fz * fga * pb_xy

                     + 6.0 * fx * fx * fz * pb_xy

                     + 30.0 * pa_zz * fz * fx * pb_xy

                     + 12.0 * pa_zzzz * fz * pb_xy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pa_zzzz,
                                   double pb_x,
                                   double pb_xz,
                                   double s_0_0)
    {
        return s_0_0 * (3.0 * pa_z * fx * fx * pb_x

                     + 2.0 * pa_zzz * fx * pb_x

                     + 0.75 * fx * fx * pb_xz

                     + 3.0 * pa_zz * fx * pb_xz

                     + pa_zzzz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_xz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pa_zzzz,
                                   double pb_x,
                                   double pb_xz,
                                   double r_0_0)
    {
        return r_0_0 * (- 6.0 * pa_z * fx * fz * fga * pb_x

                     + 24.0 * pa_z * fx * fx * fz * pb_x

                     + 20.0 * pa_zzz * fz * fx * pb_x

                     - 3.0 * fx * fz * fga * pb_xz

                     - 6.0 * pa_zz * fz * fga * pb_xz

                     + 6.0 * fx * fx * fz * pb_xz

                     + 30.0 * pa_zz * fz * fx * pb_xz

                     + 12.0 * pa_zzzz * fz * pb_xz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yy_s_0(double fx,
                                   double pa_zz,
                                   double pa_zzzz,
                                   double pb_yy,
                                   double s_0_0)
    {
        return s_0_0 * (0.375 * fx * fx * fx

                     + 1.5 * pa_zz * fx * fx

                     + 0.5 * pa_zzzz * fx

                     + 0.75 * fx * fx * pb_yy

                     + 3.0 * pa_zz * fx * pb_yy

                     + pa_zzzz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yy_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_zz,
                                   double pa_zzzz,
                                   double pb_yy,
                                   double r_0_0)
    {
        return r_0_0 * (- 0.75 * fx * fx * fz * fgb

                     - 1.5 * fx * fx * fz * fga

                     - 3.0 * pa_zz * fx * fz * fgb

                     - 3.0 * pa_zz * fz * fga * fx

                     + 2.25 * fx * fx * fx * fz

                     - pa_zzzz * fz * fgb

                     + 12.0 * pa_zz * fz * fx * fx

                     + 5.0 * pa_zzzz * fz * fx

                     - 3.0 * fx * fz * fga * pb_yy

                     - 6.0 * pa_zz * fz * fga * pb_yy

                     + 6.0 * fx * fx * fz * pb_yy

                     + 30.0 * pa_zz * fz * fx * pb_yy

                     + 12.0 * pa_zzzz * fz * pb_yy);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pa_zzzz,
                                   double pb_y,
                                   double pb_yz,
                                   double s_0_0)
    {
        return s_0_0 * (3.0 * pa_z * fx * fx * pb_y

                     + 2.0 * pa_zzz * fx * pb_y

                     + 0.75 * fx * fx * pb_yz

                     + 3.0 * pa_zz * fx * pb_yz

                     + pa_zzzz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_yz_r_0(double fga,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pa_zzzz,
                                   double pb_y,
                                   double pb_yz,
                                   double r_0_0)
    {
        return r_0_0 * (- 6.0 * pa_z * fx * fz * fga * pb_y

                     + 24.0 * pa_z * fx * fx * fz * pb_y

                     + 20.0 * pa_zzz * fz * fx * pb_y

                     - 3.0 * fx * fz * fga * pb_yz

                     - 6.0 * pa_zz * fz * fga * pb_yz

                     + 6.0 * fx * fx * fz * pb_yz

                     + 30.0 * pa_zz * fz * fx * pb_yz

                     + 12.0 * pa_zzzz * fz * pb_yz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_zz_s_0(double fx,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pa_zzzz,
                                   double pb_z,
                                   double pb_zz,
                                   double s_0_0)
    {
        return s_0_0 * (1.875 * fx * fx * fx

                     + 4.5 * pa_zz * fx * fx

                     + 6.0 * pa_z * fx * fx * pb_z

                     + 0.5 * pa_zzzz * fx

                     + 4.0 * pa_zzz * fx * pb_z

                     + 0.75 * fx * fx * pb_zz

                     + 3.0 * pa_zz * fx * pb_zz

                     + pa_zzzz * pb_zz);

    }

    #pragma omp declare simd notinbranch
    inline double fvec_zzzz_zz_r_0(double fga,
                                   double fgb,
                                   double fx,
                                   double fz,
                                   double pa_z,
                                   double pa_zz,
                                   double pa_zzz,
                                   double pa_zzzz,
                                   double pb_z,
                                   double pb_zz,
                                   double r_0_0)
    {
        return r_0_0 * (- 4.5 * fx * fx * fz * fga

                     + 11.25 * fx * fx * fx * fz

                     + 36.0 * pa_zz * fx * fx * fz

                     - 0.75 * fx * fx * fz * fgb

                     - 3.0 * pa_zz * fx * fz * fgb

                     - 3.0 * pa_zz * fz * fga * fx

                     - 12.0 * pa_z * fx * fz * fga * pb_z

                     - pa_zzzz * fz * fgb

                     + 48.0 * pa_z * fx * fx * fz * pb_z

                     + 5.0 * pa_zzzz * fz * fx

                     + 40.0 * pa_zzz * fz * fx * pb_z

                     - 3.0 * fx * fz * fga * pb_zz

                     - 6.0 * pa_zz * fz * fga * pb_zz

                     + 6.0 * fx * fx * fz * pb_zz

                     + 30.0 * pa_zz * fz * fx * pb_zz

                     + 12.0 * pa_zzzz * fz * pb_zz);

    }


} // kinrecfunc namespace

