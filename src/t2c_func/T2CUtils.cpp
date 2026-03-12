//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "T2CUtils.hpp"

#include <algorithm>
#include <ranges>
#include <cmath>

#include <boost/math/special_functions/bessel.hpp>

#include "CustomViews.hpp"
#include "TensorComponents.hpp"
#include "MathConst.hpp"

namespace t2cfunc {  // t2cfunc namespace

inline auto
bessel_il_scaled_miller_positive(int l, double x) -> double
{
    // for l>=0, x>0

    const double s0 = (1.0 - std::exp(-2.0 * x)) / (2.0 * x);

    if (l == 0) return s0;

    // Empirical choice of starting order
    const int M = l + 80 + 10 * static_cast<int>(x / 100.0);

    // Backward recurrence
    // t_{n-1} = ((2n+1)/x) t_n + t_{n+1}

    double t_np1 = 0.0; // t_{M+1}
    double t_n   = 1.0; // t_M

    double t_l = 0.0;
    double t_0 = 0.0;

    for (int n = M; n >= 1; --n)
    {
        const double t_nm1 = ((2.0 * n + 1.0) / x) * t_n + t_np1;
        t_np1 = t_n;
        t_n   = t_nm1;

        if (n - 1 == l) t_l = t_nm1;
        if (n - 1 == 0) t_0 = t_nm1;
    }

    const double scale = s0 / t_0;

    return t_l * scale;
}

inline auto
bessel_il_scaled(int l, double x) -> double
{
    // for l>=0, x>0

    const double ax = std::abs(x);

    // parity: i_l(-x) = (-1)^l i_l(x)
    const double sgn = ((x < 0.0) && (l & 1)) ? -1.0 : 1.0;

    if (ax == 0.0) return (l == 0) ? 1.0 : 0.0;

    // Threshold for using Miller recurrence
    constexpr double MILLER_SWITCH = 30.0;

    if (ax >= MILLER_SWITCH) {
        return sgn * bessel_il_scaled_miller_positive(l, ax);
    }

    // For smaller x
    const double nu = l + 0.5;
    const double I = boost::math::cyl_bessel_i(nu, ax);
    const double i_l = std::sqrt(M_PI / (2.0 * ax)) * I;

    return sgn * std::exp(-ax) * i_l;
}

auto
comp_distances_ab(CSimdArray<double>& buffer, const size_t index_ab, const size_t index_b, const TPoint<double>& r_a) -> void
{
    // set up R(AB) distances

    auto ab_x = buffer.data(index_ab);

    auto ab_y = buffer.data(index_ab + 1);

    auto ab_z = buffer.data(index_ab + 2);

    // set up Cartesian B coordinates

    auto b_x = buffer.data(index_b);

    auto b_y = buffer.data(index_b + 1);

    auto b_z = buffer.data(index_b + 2);

    // set up Cartesian A coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // compute R(AB) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(ab_x, ab_y, ab_z, b_x, b_y, b_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ab_x[i] = a_x - b_x[i];

        ab_y[i] = a_y - b_y[i];

        ab_z[i] = a_z - b_z[i];
    }
}

auto
comp_coordinates_p(CSimdArray<double>& buffer, const size_t index_p, const size_t index_b, const TPoint<double>& r_a, const double a_exp) -> void
{
    // Set up exponents

    auto b_exps = buffer.data(0);

    // set up Cartesian P coordinates

    auto p_x = buffer.data(index_p);

    auto p_y = buffer.data(index_p + 1);

    auto p_z = buffer.data(index_p + 2);

    // set up Cartesian B coordinates

    auto b_x = buffer.data(index_b);

    auto b_y = buffer.data(index_b + 1);

    auto b_z = buffer.data(index_b + 2);

    // set up Cartesian A coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // compute R(AB) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(p_x, p_y, p_z, b_x, b_y, b_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double fact = 1.0 / (a_exp + b_exps[i]);

        p_x[i] = fact * (a_x * a_exp + b_x[i] * b_exps[i]);

        p_y[i] = fact * (a_y * a_exp + b_y[i] * b_exps[i]);

        p_z[i] = fact * (a_z * a_exp + b_z[i] * b_exps[i]);
    }
}

auto
comp_distances_pb(CSimdArray<double>& buffer, const size_t index_pb, const size_t index_ab, const double a_exp) -> void
{
    // Set up exponents

    auto b_exps = buffer.data(0);

    // set up R(PB) distances

    auto pb_x = buffer.data(index_pb);

    auto pb_y = buffer.data(index_pb + 1);

    auto pb_z = buffer.data(index_pb + 2);

    // set up R(AB) distances

    auto ab_x = buffer.data(index_ab);

    auto ab_y = buffer.data(index_ab + 1);

    auto ab_z = buffer.data(index_ab + 2);

    // compute R(PB) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(pb_x, pb_y, pb_z, ab_x, ab_y, ab_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double fact = a_exp / (a_exp + b_exps[i]);

        pb_x[i] = fact * ab_x[i];

        pb_y[i] = fact * ab_y[i];

        pb_z[i] = fact * ab_z[i];
    }
}

auto
comp_distances_pa(CSimdArray<double>& buffer, const size_t index_pa, const size_t index_ab, const double a_exp) -> void
{
    // Set up exponents

    auto b_exps = buffer.data(0);

    // set up R(PA) distances

    auto pa_x = buffer.data(index_pa);

    auto pa_y = buffer.data(index_pa + 1);

    auto pa_z = buffer.data(index_pa + 2);

    // set up R(AB) distances

    auto ab_x = buffer.data(index_ab);

    auto ab_y = buffer.data(index_ab + 1);

    auto ab_z = buffer.data(index_ab + 2);

    // compute R(PA) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(pa_x, pa_y, pa_z, ab_x, ab_y, ab_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double fact = -b_exps[i] / (a_exp + b_exps[i]);

        pa_x[i] = fact * ab_x[i];

        pa_y[i] = fact * ab_y[i];

        pa_z[i] = fact * ab_z[i];
    }
}

auto
comp_distances_pb_from_p(CSimdArray<double>& buffer, const size_t index_pb, const size_t index_p, const size_t index_b) -> void
{
    // set up R(PB) distances

    auto pb_x = buffer.data(index_pb);

    auto pb_y = buffer.data(index_pb + 1);

    auto pb_z = buffer.data(index_pb + 2);

    // set up Cartesian P coordinates

    auto p_x = buffer.data(index_p);

    auto p_y = buffer.data(index_p + 1);

    auto p_z = buffer.data(index_p + 2);

    // set up Cartesian B coordinates

    auto b_x = buffer.data(index_b);

    auto b_y = buffer.data(index_b + 1);

    auto b_z = buffer.data(index_b + 2);

    // compute R(PB) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(pb_x, pb_y, pb_z, p_x, p_y, p_z, b_x, b_y, b_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        pb_x[i] = p_x[i] - b_x[i];

        pb_y[i] = p_y[i] - b_y[i];

        pb_z[i] = p_z[i] - b_z[i];
    }
}

auto
comp_distances_pa_from_p(CSimdArray<double>& buffer, const size_t index_pa, const size_t index_p, const TPoint<double>& r_a) -> void
{
    // set up R(PA) distances

    auto pa_x = buffer.data(index_pa);

    auto pa_y = buffer.data(index_pa + 1);

    auto pa_z = buffer.data(index_pa + 2);

    // set up Cartesian P coordinates

    auto p_x = buffer.data(index_p);

    auto p_y = buffer.data(index_p + 1);

    auto p_z = buffer.data(index_p + 2);

    // set up Cartesian A coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // compute R(PA) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(pa_x, pa_y, pa_z, p_x, p_y, p_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        pa_x[i] = p_x[i] - a_x;

        pa_y[i] = p_y[i] - a_y;

        pa_z[i] = p_z[i] - a_z;
    }
}

auto
comp_distances_pc(CSimdArray<double>& buffer, const size_t index_pc, const size_t index_p, const TPoint<double>& r_c) -> void
{
    t2cfunc::comp_distances_pa_from_p(buffer, index_pc, index_p, r_c);
}

auto
comp_coordinates_r(CSimdArray<double>& buffer, const size_t index_r, const size_t index_b, const TPoint<double>& r_a, const double a_exp, const double c_exp)
    -> void
{
    // set up exponents

    auto b_exps = buffer.data(0);

    // set up Cartesian R coordinates

    auto r_x = buffer.data(index_r);

    auto r_y = buffer.data(index_r + 1);

    auto r_z = buffer.data(index_r + 2);

    // set up Cartesian B coordinates

    auto b_x = buffer.data(index_b);

    auto b_y = buffer.data(index_b + 1);

    auto b_z = buffer.data(index_b + 2);

    // set up Cartesian A coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // compute R coordinates

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(r_x, r_y, r_z, b_x, b_y, b_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double fact = 1.0 / (a_exp + b_exps[i] + c_exp);

        r_x[i] = fact * (a_x * a_exp + b_x[i] * b_exps[i]);

        r_y[i] = fact * (a_y * a_exp + b_y[i] * b_exps[i]);

        r_z[i] = fact * (a_z * a_exp + b_z[i] * b_exps[i]);
    }
}

auto
comp_distances_rb(CSimdArray<double>& buffer, const size_t index_rb, const size_t index_r, const size_t index_b) -> void
{
    // set up R(RB) distances

    auto rb_x = buffer.data(index_rb);

    auto rb_y = buffer.data(index_rb + 1);

    auto rb_z = buffer.data(index_rb + 2);

    // set up Cartesian R coordinates

    auto r_x = buffer.data(index_r);

    auto r_y = buffer.data(index_r + 1);

    auto r_z = buffer.data(index_r + 2);

    // set up Cartesian B coordinates

    auto b_x = buffer.data(index_b);

    auto b_y = buffer.data(index_b + 1);

    auto b_z = buffer.data(index_b + 2);

    // compute R(RB) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(rb_x, rb_y, rb_z, r_x, r_y, r_z, b_x, b_y, b_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        rb_x[i] = r_x[i] - b_x[i];

        rb_y[i] = r_y[i] - b_y[i];

        rb_z[i] = r_z[i] - b_z[i];
    }
}

auto
comp_distances_ra(CSimdArray<double>& buffer, const size_t index_ra, const size_t index_r, const TPoint<double>& r_a) -> void
{
    // set up R(RA) distances

    auto ra_x = buffer.data(index_ra);

    auto ra_y = buffer.data(index_ra + 1);

    auto ra_z = buffer.data(index_ra + 2);

    // set up Cartesian R coordinates

    auto r_x = buffer.data(index_r);

    auto r_y = buffer.data(index_r + 1);

    auto r_z = buffer.data(index_r + 2);

    // set up Cartesian A coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // compute R(PA) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(ra_x, ra_y, ra_z, r_x, r_y, r_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ra_x[i] = r_x[i] - a_x;

        ra_y[i] = r_y[i] - a_y;

        ra_z[i] = r_z[i] - a_z;
    }
}

auto
comp_inverted_zeta(CSimdArray<double>& buffer, const size_t index_fact, const double a_exp, const double c_exp) -> void
{
    // set up exponents

    auto b_exps = buffer.data(0);

    // set up zeta factors

    auto facts = buffer.data(index_fact);

    // compute R coordinates

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(facts, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        facts[i] = 0.5 / (a_exp + b_exps[i] + c_exp);
    }
}

auto
comp_coordinates_norm(CSimdArray<double>& buffer, const size_t index_mr, const size_t index_r) -> void
{
    // set up Cartesian R coordinates

    auto r_x = buffer.data(index_r);

    auto r_y = buffer.data(index_r + 1);

    auto r_z = buffer.data(index_r + 2);
    
    // set up |r| norms

    auto mr = buffer.data(index_mr);

    const auto nelems = buffer.number_of_active_elements();

    #pragma omp simd aligned(mr, r_x, r_y, r_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        mr[i] = std::sqrt(r_x[i] * r_x[i] + r_y[i] * r_y[i] + r_z[i] * r_z[i]);
    }
}

auto
comp_legendre_args(CSimdArray<double>& buffer, const size_t index_args, const size_t index_b, const size_t index_mb,  const TPoint<double>& r_a) -> void
{
    // set up Legendre arguments

    auto fargs = buffer.data(index_args);
    
    // set up Cartesian B coordinates

    auto b_x = buffer.data(index_b);

    auto b_y = buffer.data(index_b + 1);

    auto b_z = buffer.data(index_b + 2);
    
    // set up |B|

    auto mb = buffer.data(index_mb);

    // set up Cartesian A coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];
    
    // compute A.B

    const auto nelems = buffer.number_of_active_elements();

    #pragma omp simd aligned(fargs, b_x, b_y, b_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        fargs[i] = a_x * b_x[i] + a_y * b_y[i] + a_z * b_z[i];
    }
    
    // compute |A|
    
    const auto ma = std::sqrt(a_x * a_x + a_y * a_y + a_z * a_z);
    
    // conditionally scale A.B by |A||B|
    
    for (size_t i = 0; i < nelems; i++)
    {
        if (const auto fact = ma * mb[i]; fact > 1.0e-12)
        {
            fargs[i] /= fact;
        }
        else
        {
            fargs[i] = 1.0;
        }
    }
}

auto
comp_gamma_factors(CSimdArray<double>& buffer, const size_t index_gf, const size_t index_mb, const TPoint<double>& r_a, const double a_exp, const double c_exp) -> void
{
    // set up pi value
    
    const double fpi = mathconst::pi_value();
    
    // set up gamma factors

    auto fg = buffer.data(index_gf);
    
    // set up exponents

    auto b_exps = buffer.data(0);
    
    // set up B coordinates norms

    auto mb = buffer.data(index_mb);
    
    // set up Cartesian A coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];
    
    // compute R coordinates

    const auto nelems = buffer.number_of_active_elements();
    
    const double a2 = a_x * a_x + a_y * a_y + a_z * a_z;

    #pragma omp simd aligned(b_exps, mb, fg : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double fzi = 1.0 / (b_exps[i] + a_exp + c_exp);
        
        double fa = a_exp * (b_exps[i] + c_exp) * fzi;
        
        double fb = b_exps[i] * (a_exp + c_exp) * fzi;
        
        double fact = fpi * fzi * std::sqrt(fpi * fzi);
        
        fg[i] = fact * std::exp(-fa * a2 - fb * mb[i] * mb[i]);
    }
}

auto
comp_bessel_args(CSimdArray<double>& buffer, const size_t index_args, const size_t index_mb, const TPoint<double>& r_a, const double a_exp, const double c_exp) -> void
{
    // set up T arguments

    auto fargs = buffer.data(index_args);
    
    // set up exponents

    auto b_exps = buffer.data(0);
    
    // set up |B|

    auto mb = buffer.data(index_mb);
    
    // set up Cartesian A coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];
    
    // compute |A|
    
    const auto ma = std::sqrt(a_x * a_x + a_y * a_y + a_z * a_z);
    
    // compute Bessel arguments

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(b_exps, mb, fargs : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        fargs[i] = 2.0 * b_exps[i] * a_exp * ma * mb[i] / (b_exps[i] + a_exp + c_exp);
    }
}

auto
comp_i_vals(CSimdArray<double>& values, const int order, const CSimdArray<double>& buffer, const size_t index_args) -> void
{
    // set up T arguments

    auto fargs = buffer.data(index_args);
    
    // set up number of active elements
    
    const auto nelems = buffer.number_of_active_elements();
    
    // zero order
    
    auto f0vals = values.data(0);
    
    for (size_t i = 0; i < nelems; i++)
    {
//        if (const double fact = fargs[i]; fact > 150.0)
//        {
//            f0vals[i] = 0.0;
//        }
//        else
//        {
        const double fact = fargs[i];
        
        f0vals[i] = bessel_il_scaled(0, fact) * std::exp(fact);
        //}
    }
    
    if (order == 0) return;
    
    // first order
    
    auto f1vals = values.data(1);
    
    for (size_t i = 0; i < nelems; i++)
    {
//        if (const double fact = fargs[i]; fact > 150.0)
//        {
//            f1vals[i] = 0.0;
//        }
//        else
        if (const double fact = fargs[i]; fact <= 1.0e-12)
        {
            f1vals[i] = 1.0 / 3.0;
        }
        else
        {
            f1vals[i] = bessel_il_scaled(1, fact) * std::exp(fact) / fact;
        }
    }
    
    if (order == 1) return;
    
    // second order
    
    auto f2vals = values.data(2);
    
    for (size_t i = 0; i < nelems; i++)
    {
//        if (const double fact = fargs[i]; fact > 150.0)
//        {
//            f1vals[i] = 0.0;
//        }
//        else
        if (const double fact = fargs[i]; fact <= 1.0e-12)
        {
            f2vals[i] = 1.0/ 15.0;
        }
        else
        {
            f2vals[i] = bessel_il_scaled(2, fact) * std::exp(fact) / (fact * fact);
        }
    }
    
    if (order == 2) return;
    
    // higher orders
    
    double n2fact = 15.0;
    
    for (int k = 3; k <= order; k++)
    {
        auto fvals = values.data(size_t(k));
        
        n2fact *= 2 * (double)k + 1.0;
        
        for (size_t i = 0; i < nelems; i++)
        {
//            if (const double fact = fargs[i]; fact > 150.0)
//            {
//                f1vals[i] = 0.0;
//            }
//            else
            if (const double fact = fargs[i]; fact <= 1.0e-12)
            {
                fvals[i] = 1.0 / n2fact;
            }
            else
            {
                fvals[i] = bessel_il_scaled(k, fact) * std::exp(fact) * std::pow(fact, -(double)k);
            }
        }
    }
}

auto
comp_l_vals(CSimdArray<double>& values, const int order, const CSimdArray<double>& buffer, const size_t index_args, const size_t index_pab) -> void
{
    // set up T arguments

    auto fargs = buffer.data(index_args);
    
    // set up Legendre arguments

    auto pab = buffer.data(index_pab);
    
    // set up number of active elements
    
    const auto nelems = buffer.number_of_active_elements();
    
    // zero order
    
    auto f0vals = values.data(0);
    
    for (size_t i = 0; i < nelems; i++)
    {
        f0vals[i] = 1.0;
    }
    
    if (order == 0) return;
    
    // first order
    
    auto f1vals = values.data(1);
    
    for (size_t i = 0; i < nelems; i++)
    {
        f1vals[i] = 3.0 * fargs[i] * pab[i];
    }
    
    if (order == 1) return;
    
    // higher orders
    
    for (int k = 2; k <= order; k++)
    {
        // set up pointers
        
        auto frvals = values.data(size_t(k));
        
        auto fpvals = values.data(size_t(k - 1));
        
        auto fmvals = values.data(size_t(k - 2));
        
        const auto fk = double(k - 1);
        
        for (size_t i = 0; i < nelems; i++)
        {
            frvals[i] = (2.0 * fk + 3) * fargs[i] * (pab[i] * fpvals[i]
                                                     
                      - fk * fargs[i] * fmvals[i] / (2.0 * fk - 1.0)) / (fk + 1.0);
        }
    }
}

void
comp_boys_args(CSimdArray<double>& bf_data, const size_t index_args, CSimdArray<double>& buffer, const size_t index_pc, const double a_exp)
{
    // Set up exponents

    auto b_exps = buffer.data(0);

    // set up Boys function arguments

    auto bargs = bf_data.data(index_args);

    // set up R(PC) distances

    auto pc_x = buffer.data(index_pc);

    auto pc_y = buffer.data(index_pc + 1);

    auto pc_z = buffer.data(index_pc + 2);

    // compute Boys function arguments

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(bargs, pc_x, pc_y, pc_z, b_exps : 64)
    for (int i = 0; i < nelems; i++)
    {
        bargs[i] = (a_exp + b_exps[i]) * (pc_x[i] * pc_x[i] + pc_y[i] * pc_y[i] + pc_z[i] * pc_z[i]);
    }
}

void
comp_boys_args(CSimdArray<double>& bf_data,
               const size_t        index_args,
               CSimdArray<double>& buffer,
               const size_t        index_pc,
               const double        a_exp,
               const double        omega)
{
    // Set up exponents

    auto b_exps = buffer.data(0);

    // set up Boys function arguments

    auto bargs = bf_data.data(index_args);

    // set up R(PC) distances

    auto pc_x = buffer.data(index_pc);

    auto pc_y = buffer.data(index_pc + 1);

    auto pc_z = buffer.data(index_pc + 2);

    // compute Boys function arguments

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(bargs, pc_x, pc_y, pc_z, b_exps : 64)
    for (int i = 0; i < nelems; i++)
    {
        const double frho = a_exp + b_exps[i];

        bargs[i] = omega * omega / (omega * omega + frho);

        bargs[i] *= frho * (pc_x[i] * pc_x[i] + pc_y[i] * pc_y[i] + pc_z[i] * pc_z[i]);
    }
}

void
comp_boys_args_with_rho(CSimdArray<double>& bf_data, const size_t index_args, CSimdArray<double>& buffer, const size_t index_ab, const double a_exp)
{
    // Set up exponents

    auto b_exps = buffer.data(0);

    // set up Boys function arguments

    auto bargs = bf_data.data(index_args);

    // set up R(AB) distances

    auto ab_x = buffer.data(index_ab);

    auto ab_y = buffer.data(index_ab + 1);

    auto ab_z = buffer.data(index_ab + 2);

    // compute Boys function arguments

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(bargs, ab_x, ab_y, ab_z, b_exps : 64)
    for (int i = 0; i < nelems; i++)
    {
        bargs[i] = a_exp * b_exps[i] * (ab_x[i] * ab_x[i] + ab_y[i] * ab_y[i] + ab_z[i] * ab_z[i]) / (a_exp + b_exps[i]);
    }
}

auto
reduce(CSimdArray<double>& cbuffer, CSimdArray<double>& pbuffer, const size_t position, const size_t ndims, const size_t nblocks) -> void
{
    if (const auto nrows = cbuffer.number_of_rows(); nrows > 0)
    {
        std::ranges::for_each(views::rectangular(nrows, nblocks), [&](const auto& index) {
            const auto [i, j] = index;
            auto pdata        = pbuffer.data(position + i);
            auto cdata        = cbuffer.data(i);
            auto pvals        = &pdata[j * ndims];
#pragma omp simd
            for (size_t k = 0; k < ndims; k++)
            {
                cdata[k] += pvals[k];
            }
        });
    }
}

auto
reduce(CSimdArray<double>& cbuffer, const size_t cposition, CSimdArray<double>& pbuffer, const size_t pposition, const size_t nrows, const size_t ndims, const size_t nblocks) -> void
{
    if (nrows > 0)
    {
        std::ranges::for_each(views::rectangular(nrows, nblocks), [&](const auto& index) {
            const auto [i, j] = index;
            auto pdata        = pbuffer.data(pposition + i);
            auto cdata        = cbuffer.data(cposition + i);
            auto pvals        = &pdata[j * ndims];
#pragma omp simd
            for (size_t k = 0; k < ndims; k++)
            {
                cdata[k] += pvals[k];
            }
        });
    }
}

auto
reduce(CSimdArray<double>& cbuffer, CSimdArray<double>& pbuffer, const size_t position, const double factor, const size_t ndims, const size_t nblocks)
    -> void
{
    if (const auto nrows = cbuffer.number_of_rows(); nrows > 0)
    {
        std::ranges::for_each(views::rectangular(nrows, nblocks), [&](const auto& index) {
            const auto [i, j] = index;
            auto pdata        = pbuffer.data(position + i);
            auto cdata        = cbuffer.data(i);
            auto pvals        = &pdata[j * ndims];
#pragma omp simd
            for (size_t k = 0; k < ndims; k++)
            {
                cdata[k] += factor * pvals[k];
            }
        });
    }
}

auto
reduce(CSimdArray<double>&        cbuffer,
       CSimdArray<double>&        pbuffer,
       const size_t               position,
       const std::vector<double>& facts,
       const size_t               nfacts,
       const size_t               index,
       const size_t               ndims,
       const size_t               nblocks) -> void
{
    if (const auto nrows = cbuffer.number_of_rows(); nrows > 0)
    {
        const auto loc_rows = nrows / nfacts;

        const auto ncenters = facts.size() / nfacts;

        std::ranges::for_each(std::views::iota(size_t{0}, nfacts), [&](const size_t idx) {
            const auto factor = facts[idx * ncenters + index];
            std::ranges::for_each(views::rectangular(loc_rows, nblocks), [&](const auto& rindex) {
                const auto [i, j] = rindex;
                auto pdata        = pbuffer.data(idx * loc_rows + position + i);
                auto cdata        = cbuffer.data(idx * loc_rows + i);
                auto pvals        = &pdata[j * ndims];
#pragma omp simd
                for (size_t k = 0; k < ndims; k++)
                {
                    cdata[k] += factor * pvals[k];
                }
            });
        });
    }
}

auto
distribute(CSubMatrix*                      matrix,
           const double*                    buffer,
           const std::vector<size_t>&       indices,
           const int                        bra_comp,
           const int                        ket_comp,
           const size_t                     bra_igto,
           const std::pair<size_t, size_t>& ket_range) -> void
{
    const auto bra_idx = indices[0] * bra_comp + indices[bra_igto + 1];

    std::ranges::for_each(std::views::iota(ket_range.first, ket_range.second),
                          [&](const auto i) { matrix->operator[]({bra_idx, indices[0] * ket_comp + indices[i + 1]}) = buffer[i - ket_range.first]; });
}

auto
distribute(CSubMatrix*                      matrix,
           const double*                    buffer,
           const std::vector<size_t>&       bra_indices,
           const std::vector<size_t>&       ket_indices,
           const int                        bra_comp,
           const int                        ket_comp,
           const size_t                     bra_igto,
           const std::pair<size_t, size_t>& ket_range,
           const mat_t                      mat_type) -> void
{
    const auto bra_idx = bra_indices[0] * bra_comp + bra_indices[bra_igto + 1];

    std::ranges::for_each(std::views::iota(ket_range.first, ket_range.second), [&](const auto i) {
        const auto ket_idx                     = ket_indices[0] * ket_comp + ket_indices[i + 1];
        const auto fval                        = buffer[i - ket_range.first];
        matrix->operator[]({bra_idx, ket_idx}) = fval;
        if (mat_type == mat_t::symmetric)
        {
            matrix->operator[]({ket_idx, bra_idx}) = fval;
        }
        if (mat_type == mat_t::antisymmetric)
        {
            matrix->operator[]({ket_idx, bra_idx}) = -fval;
        }
    });
}

auto
distribute(CSubMatrix*                      matrix,
           const double*                    buffer,
           const std::vector<size_t>&       bra_indices,
           const std::vector<size_t>&       ket_indices,
           const int                        bra_comp,
           const int                        ket_comp,
           const size_t                     bra_igto,
           const std::pair<size_t, size_t>& ket_range,
           const bool                       ang_order) -> void
{
    const auto bra_idx = bra_indices[0] * bra_comp + bra_indices[bra_igto + 1];

    std::ranges::for_each(std::views::iota(ket_range.first, ket_range.second), [&](const auto i) {
        const auto ket_idx = ket_indices[0] * ket_comp + ket_indices[i + 1];
        if (ang_order)
        {
            matrix->operator[]({bra_idx, ket_idx}) = buffer[i - ket_range.first];
        }
        else
        {
            matrix->operator[]({ket_idx, bra_idx}) = buffer[i - ket_range.first];
        }
    });
}

auto
distribute(CSubMatrix*                      matrix,
           const CSimdArray<double>&        buffer,
           const size_t                     offset,
           const std::vector<size_t>&       indices,
           const int                        bra_angmom,
           const size_t                     bra_igto,
           const std::pair<size_t, size_t>& ket_range) -> void
{
    const auto bra_comps = tensor::number_of_spherical_components(std::array<int, 1>{bra_angmom});

    std::ranges::for_each(views::rectangular(bra_comps, bra_comps), [&](const auto& index) {
        const auto [i, j] = index;
        t2cfunc::distribute(matrix, buffer.data(offset + i * bra_comps + j), indices, i, j, bra_igto, ket_range);
    });
}

auto
distribute(CSubMatrix*                      matrix,
           const CSimdArray<double>&        buffer,
           const size_t                     offset,
           const std::vector<size_t>&       bra_indices,
           const std::vector<size_t>&       ket_indices,
           const int                        bra_angmom,
           const int                        ket_angmom,
           const size_t                     bra_igto,
           const std::pair<size_t, size_t>& ket_range,
           const mat_t                      mat_type) -> void
{
    const auto bra_comps = tensor::number_of_spherical_components(std::array<int, 1>{bra_angmom});

    const auto ket_comps = tensor::number_of_spherical_components(std::array<int, 1>{ket_angmom});

    std::ranges::for_each(views::rectangular(bra_comps, ket_comps), [&](const auto& index) {
        const auto [i, j] = index;
        t2cfunc::distribute(matrix, buffer.data(offset + i * ket_comps + j), bra_indices, ket_indices, i, j, bra_igto, ket_range, mat_type);
    });
}

auto
distribute(CSubMatrix*                      matrix,
           const CSimdArray<double>&        buffer,
           const size_t                     offset,
           const std::vector<size_t>&       bra_indices,
           const std::vector<size_t>&       ket_indices,
           const int                        bra_angmom,
           const int                        ket_angmom,
           const size_t                     bra_igto,
           const std::pair<size_t, size_t>& ket_range,
           const bool                       ang_order) -> void
{
    const auto bra_comps = tensor::number_of_spherical_components(std::array<int, 1>{bra_angmom});

    const auto ket_comps = tensor::number_of_spherical_components(std::array<int, 1>{ket_angmom});

    std::ranges::for_each(views::rectangular(bra_comps, ket_comps), [&](const auto& index) {
        const auto [i, j] = index;
        t2cfunc::distribute(matrix, buffer.data(offset + i * ket_comps + j), bra_indices, ket_indices, i, j, bra_igto, ket_range, ang_order);
    });
}

}  // namespace t2cfunc
