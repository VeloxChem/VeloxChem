#include "T3CUtils.hpp"

#include "MathConst.hpp"

namespace t3cfunc {  // t3cfunc namespace

auto
comp_distances_aq(CSimdArray<double>&   buffer,
                  const size_t          index_aq,
                  const size_t          index_q,
                  const TPoint<double>& r_a) -> void
{
    // set up R(AQ) distances

    auto aq_x = buffer.data(index_aq);

    auto aq_y = buffer.data(index_aq + 1);

    auto aq_z = buffer.data(index_aq + 2);

    // set up Cartesian Q coordinates

    auto q_x = buffer.data(index_q);

    auto q_y = buffer.data(index_q + 1);

    auto q_z = buffer.data(index_q + 2);

    // set up Cartesian P coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // compute R(PQ) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(aq_x, aq_y, aq_z, q_x, q_y, q_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        aq_x[i] = a_x - q_x[i];

        aq_y[i] = a_y - q_y[i];

        aq_z[i] = a_z - q_z[i];
    }
}

auto
comp_distances_wa(CSimdArray<double>&   buffer,
                  const size_t          index_wa,
                  const size_t          index_w,
                  const TPoint<double>& r_a) -> void
{
    // set up R(WA) distances

    auto wa_x = buffer.data(index_wa);

    auto wa_y = buffer.data(index_wa + 1);

    auto wa_z = buffer.data(index_wa + 2);

    // set up Cartesian W coordinates

    auto w_x = buffer.data(index_w);

    auto w_y = buffer.data(index_w + 1);

    auto w_z = buffer.data(index_w + 2);

    // set up Cartesian P coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // compute R(WQ) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(wa_x, wa_y, wa_z, w_x, w_y, w_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        wa_x[i] = w_x[i] - a_x;

        wa_y[i] = w_y[i] - a_y;

        wa_z[i] = w_z[i] - a_z;
    }
}

auto
comp_boys_args(CSimdArray<double>&       bf_data,
               const size_t              index_args,
               const CSimdArray<double>& buffer,
               const size_t              index_aq,
               const double              a_exp) -> void
{
    // Set up exponents

    auto c_exps = buffer.data(0);

    auto d_exps = buffer.data(1);

    // set up R(PQ) distances

    auto aq_x = buffer.data(index_aq);

    auto aq_y = buffer.data(index_aq + 1);

    auto aq_z = buffer.data(index_aq + 2);

    // set up Boys function arguments

    auto bargs = bf_data.data(index_args);

    // compute Boys function arguments

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(bargs, aq_x, aq_y, aq_z, c_exps, d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double cd_exp = c_exps[i] + d_exps[i];

        bargs[i] = a_exp * cd_exp * (aq_x[i] * aq_x[i] + aq_y[i] * aq_y[i] + aq_z[i] * aq_z[i]) / (a_exp + cd_exp);
    }
}

auto
comp_ovl_factors(CSimdArray<double>& buffer,
                 const size_t        index_ovl,
                 const size_t        index_ket_ovl,
                 const size_t        index_ket_norm,
                 const double        a_norm,
                 const double        a_exp) -> void
{
    // set up exponents

    auto c_exps = buffer.data(0);

    auto d_exps = buffer.data(1);

    // set up combined overlap

    auto fss = buffer.data(index_ovl);

    // set up ket data

    auto ket_ovls = buffer.data(index_ket_ovl);

    auto ket_norms = buffer.data(index_ket_norm);

    // set up pi constant

    const auto fpi = mathconst::pi_value();

    // compute combined overlap factors

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(fss, ket_ovls, ket_norms, c_exps, d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double cd_exp = c_exps[i] + d_exps[i];

        fss[i] = 2.0 * fpi * a_norm * ket_norms[i] * ket_ovls[i] * std::sqrt(cd_exp / (a_exp + cd_exp)) / a_exp;
    }
}


}  // namespace t3cfunc
