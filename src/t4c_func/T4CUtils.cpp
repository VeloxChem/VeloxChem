#include "T4CUtils.hpp"

#include <cmath>
#include <iostream>
#include <set>

#include "MathConst.hpp"
#include "T4CLocalDistributor.hpp"
#include "T4CLocalGeomDistributor.hpp"

namespace t4cfunc {  // t2cfunc namespace

auto
comp_coordinates_q(CSimdArray<double>& buffer, const size_t index_q, const size_t index_c, const size_t index_d) -> void
{
    // Set up exponents

    auto c_exps = buffer.data(0);

    auto d_exps = buffer.data(1);

    // set up Cartesian Q coordinates

    auto q_x = buffer.data(index_q);

    auto q_y = buffer.data(index_q + 1);

    auto q_z = buffer.data(index_q + 2);

    // set up Cartesian C coordinates

    auto c_x = buffer.data(index_c);

    auto c_y = buffer.data(index_c + 1);

    auto c_z = buffer.data(index_c + 2);

    // set up Cartesian D coordinates

    auto d_x = buffer.data(index_d);

    auto d_y = buffer.data(index_d + 1);

    auto d_z = buffer.data(index_d + 2);

    // compute Cartesian Q coordinates

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(q_x, q_y, q_z, c_x, c_y, c_z, d_x, d_y, d_z, c_exps, d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fact = 1.0 / (c_exps[i] + d_exps[i]);

        q_x[i] = fact * (c_x[i] * c_exps[i] + d_x[i] * d_exps[i]);

        q_y[i] = fact * (c_y[i] * c_exps[i] + d_y[i] * d_exps[i]);

        q_z[i] = fact * (c_z[i] * c_exps[i] + d_z[i] * d_exps[i]);
    }
}

auto
comp_coordinates_w(CSimdArray<double>&   buffer,
                   const size_t          index_w,
                   const size_t          index_q,
                   const TPoint<double>& r_p,
                   const double          a_exp,
                   const double          b_exp) -> void
{
    // Set up exponents

    auto c_exps = buffer.data(0);

    auto d_exps = buffer.data(1);

    // set up Cartesian W coordinates

    auto w_x = buffer.data(index_w);

    auto w_y = buffer.data(index_w + 1);

    auto w_z = buffer.data(index_w + 2);

    // set up Cartesian Q coordinates

    auto q_x = buffer.data(index_q);

    auto q_y = buffer.data(index_q + 1);

    auto q_z = buffer.data(index_q + 2);

    // set up Cartesian P coordinates

    const auto xyz = r_p.coordinates();

    const auto p_x = xyz[0];

    const auto p_y = xyz[1];

    const auto p_z = xyz[2];

    // compute Cartesian W center coordinates

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(w_x, w_y, w_z, q_x, q_y, q_z, c_exps, d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double ab_exp = a_exp + b_exp;

        const double cd_exp = c_exps[i] + d_exps[i];

        const double fact = 1.0 / (ab_exp + cd_exp);

        w_x[i] = fact * (p_x * ab_exp + q_x[i] * cd_exp);

        w_y[i] = fact * (p_y * ab_exp + q_y[i] * cd_exp);

        w_z[i] = fact * (p_z * ab_exp + q_z[i] * cd_exp);
    }
}

auto
comp_distances_pq(CSimdArray<double>& buffer, const size_t index_pq, const size_t index_q, const TPoint<double>& r_p) -> void
{
    // set up R(PQ) distances

    auto pq_x = buffer.data(index_pq);

    auto pq_y = buffer.data(index_pq + 1);

    auto pq_z = buffer.data(index_pq + 2);

    // set up Cartesian Q coordinates

    auto q_x = buffer.data(index_q);

    auto q_y = buffer.data(index_q + 1);

    auto q_z = buffer.data(index_q + 2);

    // set up Cartesian P coordinates

    const auto xyz = r_p.coordinates();

    const auto p_x = xyz[0];

    const auto p_y = xyz[1];

    const auto p_z = xyz[2];

    // compute R(PQ) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(pq_x, pq_y, pq_z, q_x, q_y, q_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        pq_x[i] = p_x - q_x[i];

        pq_y[i] = p_y - q_y[i];

        pq_z[i] = p_z - q_z[i];
    }
}

auto
comp_distances_wq(CSimdArray<double>& buffer, const size_t index_wq, const size_t index_w, const size_t index_q) -> void
{
    // set up R(WQ) distances

    auto wq_x = buffer.data(index_wq);

    auto wq_y = buffer.data(index_wq + 1);

    auto wq_z = buffer.data(index_wq + 2);

    // set up Cartesian W coordinates

    auto w_x = buffer.data(index_w);

    auto w_y = buffer.data(index_w + 1);

    auto w_z = buffer.data(index_w + 2);

    // set up Cartesian Q coordinates

    auto q_x = buffer.data(index_q);

    auto q_y = buffer.data(index_q + 1);

    auto q_z = buffer.data(index_q + 2);

    // compute R(WQ) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(wq_x, wq_y, wq_z, w_x, w_y, w_z, q_x, q_y, q_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        wq_x[i] = w_x[i] - q_x[i];

        wq_y[i] = w_y[i] - q_y[i];

        wq_z[i] = w_z[i] - q_z[i];
    }
}

auto
comp_distances_wp(CSimdArray<double>& buffer, const size_t index_wp, const size_t index_w, const TPoint<double>& r_p) -> void
{
    // set up R(WP) distances

    auto wp_x = buffer.data(index_wp);

    auto wp_y = buffer.data(index_wp + 1);

    auto wp_z = buffer.data(index_wp + 2);

    // set up Cartesian W coordinates

    auto w_x = buffer.data(index_w);

    auto w_y = buffer.data(index_w + 1);

    auto w_z = buffer.data(index_w + 2);

    // set up Cartesian P coordinates

    const auto xyz = r_p.coordinates();

    const auto p_x = xyz[0];

    const auto p_y = xyz[1];

    const auto p_z = xyz[2];

    // compute R(WQ) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(wp_x, wp_y, wp_z, w_x, w_y, w_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        wp_x[i] = w_x[i] - p_x;

        wp_y[i] = w_y[i] - p_y;

        wp_z[i] = w_z[i] - p_z;
    }
}

auto
comp_distances_cd(CSimdArray<double>& buffer, const size_t index_cd, const size_t index_c, const size_t index_d) -> void
{
    // set up R(CD) distances

    auto cd_x = buffer.data(index_cd);

    auto cd_y = buffer.data(index_cd + 1);

    auto cd_z = buffer.data(index_cd + 2);

    // set up Cartesian C coordinates

    auto c_x = buffer.data(index_c);

    auto c_y = buffer.data(index_c + 1);

    auto c_z = buffer.data(index_c + 2);

    // set up Cartesian D coordinates

    auto d_x = buffer.data(index_d);

    auto d_y = buffer.data(index_d + 1);

    auto d_z = buffer.data(index_d + 2);

    // compute R(WQ) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(cd_x, cd_y, cd_z, c_x, c_y, c_z, d_x, d_y, d_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        cd_x[i] = c_x[i] - d_x[i];

        cd_y[i] = c_y[i] - d_y[i];

        cd_z[i] = c_z[i] - d_z[i];
    }
}

auto
comp_distances_qd(CSimdArray<double>& buffer, const size_t index_qd, const size_t index_q, const size_t index_d) -> void
{
    // set up R(QD) distances

    auto qd_x = buffer.data(index_qd);

    auto qd_y = buffer.data(index_qd + 1);

    auto qd_z = buffer.data(index_qd + 2);

    // set up Cartesian Q coordinates

    auto q_x = buffer.data(index_q);

    auto q_y = buffer.data(index_q + 1);

    auto q_z = buffer.data(index_q + 2);

    // set up Cartesian D coordinates

    auto d_x = buffer.data(index_d);

    auto d_y = buffer.data(index_d + 1);

    auto d_z = buffer.data(index_d + 2);

    // compute R(WQ) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(qd_x, qd_y, qd_z, q_x, q_y, q_z, d_x, d_y, d_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        qd_x[i] = q_x[i] - d_x[i];

        qd_y[i] = q_y[i] - d_y[i];

        qd_z[i] = q_z[i] - d_z[i];
    }
}

auto
comp_distances_qc(CSimdArray<double>& buffer, const size_t index_qc, const size_t index_q, const size_t index_c) -> void
{
    // set up R(QC) distances

    auto qc_x = buffer.data(index_qc);

    auto qc_y = buffer.data(index_qc + 1);

    auto qc_z = buffer.data(index_qc + 2);

    // set up Cartesian Q coordinates

    auto q_x = buffer.data(index_q);

    auto q_y = buffer.data(index_q + 1);

    auto q_z = buffer.data(index_q + 2);

    // set up Cartesian C coordinates

    auto c_x = buffer.data(index_c);

    auto c_y = buffer.data(index_c + 1);

    auto c_z = buffer.data(index_c + 2);

    // compute R(WQ) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(qc_x, qc_y, qc_z, q_x, q_y, q_z, c_x, c_y, c_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        qc_x[i] = q_x[i] - c_x[i];

        qc_y[i] = q_y[i] - c_y[i];

        qc_z[i] = q_z[i] - c_z[i];
    }
}

auto
comp_boys_args(CSimdArray<double>&       bf_data,
               const size_t              index_args,
               const CSimdArray<double>& buffer,
               const size_t              index_pq,
               const double              a_exp,
               const double              b_exp) -> void
{
    // Set up exponents

    auto c_exps = buffer.data(0);

    auto d_exps = buffer.data(1);

    // set up R(PQ) distances

    auto pq_x = buffer.data(index_pq);

    auto pq_y = buffer.data(index_pq + 1);

    auto pq_z = buffer.data(index_pq + 2);

    // set up Boys function arguments

    auto bargs = bf_data.data(index_args);

    // compute Boys function arguments

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(bargs, pq_x, pq_y, pq_z, c_exps, d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double ab_exp = a_exp + b_exp;

        double cd_exp = c_exps[i] + d_exps[i];

        bargs[i] = ab_exp * cd_exp * (pq_x[i] * pq_x[i] + pq_y[i] * pq_y[i] + pq_z[i] * pq_z[i]) / (ab_exp + cd_exp);
    }
}

auto
comp_boys_args(CSimdArray<double>&       bf_data,
               const size_t              index_args,
               const CSimdArray<double>& buffer,
               const size_t              index_pq,
               const double              a_exp,
               const double              b_exp,
               const double              omega) -> void
{
    // Set up exponents

    auto c_exps = buffer.data(0);

    auto d_exps = buffer.data(1);

    // set up R(PQ) distances

    auto pq_x = buffer.data(index_pq);

    auto pq_y = buffer.data(index_pq + 1);

    auto pq_z = buffer.data(index_pq + 2);

    // set up Boys function arguments

    auto bargs = bf_data.data(index_args);

    // compute Boys function arguments

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(bargs, pq_x, pq_y, pq_z, c_exps, d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double ab_exp = a_exp + b_exp;

        double cd_exp = c_exps[i] + d_exps[i];

        double frho = ab_exp * cd_exp / (ab_exp + cd_exp);

        bargs[i] = omega * omega / (omega * omega + frho);

        bargs[i] *= frho * (pq_x[i] * pq_x[i] + pq_y[i] * pq_y[i] + pq_z[i] * pq_z[i]);
    }
}

auto
comp_ovl_factors(CSimdArray<double>& buffer,
                 const size_t        index_ovl,
                 const size_t        index_ket_ovl,
                 const size_t        index_ket_norm,
                 const double        bra_ovl,
                 const double        bra_norm,
                 const double        a_exp,
                 const double        b_exp) -> void
{
    // set up exponents

    auto c_exps = buffer.data(0);

    auto d_exps = buffer.data(1);

    // set up combined overlap

    auto fss = buffer.data(index_ovl);

    // set up ket data

    auto ket_ovls = buffer.data(index_ket_ovl);

    auto ket_norms = buffer.data(index_ket_norm);

    // set up inverted pi constant

    const auto invfpi = 1.0 / mathconst::pi_value();

    // compute combined overlap factors

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(fss, ket_ovls, ket_norms, c_exps, d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double ab_exp = a_exp + b_exp;

        double cd_exp = c_exps[i] + d_exps[i];

        fss[i] = 2.0 * bra_norm * ket_norms[i] * bra_ovl * ket_ovls[i] * std::sqrt(invfpi * ab_exp * cd_exp / (ab_exp + cd_exp));
    }
}

auto
update_max_values(std::vector<double>& max_values, const CSimdArray<double>& buffer, const size_t index) -> void
{
    for (int i = 0; i < buffer.number_of_rows(); i++)
    {
        if (const auto val = buffer.data(i)[0]; max_values[index] < val) max_values[index] = val;
    }
}

auto
store_values(CSimdArray<double>& buffer, const CSimdArray<double>& values, const int offset) -> void
{
    const auto ndims = buffer.number_of_active_elements();

    for (size_t i = 0; i < values.number_of_rows(); i++)
    {
        auto sdata = values.data(i);

        auto ddata = buffer.data(offset + i);

#pragma omp simd aligned(sdata, ddata : 64)
        for (size_t j = 0; j < ndims; j++)
        {
            ddata[j] = sdata[j];
        }
    }
}

auto
accumulate(CSubMatrix*                glob_matrix,
           const CSubMatrix*          loc_matrix,
           const std::vector<size_t>& bra_loc_indices,
           const std::vector<size_t>& ket_loc_indices,
           const std::vector<size_t>& bra_glob_indices,
           const std::vector<size_t>& ket_glob_indices,
           const int                  bra_comps,
           const int                  ket_comps,
           const bool                 ang_order) -> void
{
    // std::cout << "*** Accumulate " << bra_comps << "," << ket_comps << " *** " << std::endl;

    for (int i = 0; i < bra_comps; i++)
    {
        const auto bra_loff = i * bra_loc_indices[0];

        const auto bra_goff = i * bra_glob_indices[0];

        for (int j = 0; j < ket_comps; j++)
        {
            const auto ket_loff = j * ket_loc_indices[0];

            const auto ket_goff = j * ket_glob_indices[0];  // this need fix for mixed basis

            for (size_t k = 0; k < bra_loc_indices[0]; k++)
            {
                const auto kg = bra_goff + bra_glob_indices[k + 1];

                for (size_t l = 0; l < ket_loc_indices[0]; l++)
                {
                    const auto lg = ket_goff + ket_glob_indices[l + 1];

                    //                    std::cout << "(" << kg << "," << lg << ") = (" << bra_loff + k << "," << ket_loff + l << ") = ";
                    //
                    //                    std::cout << loc_matrix->at({bra_loff + k, ket_loff + l}) << std::endl;

                    if (ang_order)
                    {
                        glob_matrix->operator[]({kg, lg}) += loc_matrix->at({bra_loff + k, ket_loff + l});
                    }
                    else
                    {
                        glob_matrix->operator[]({lg, kg}) += loc_matrix->at({bra_loff + k, ket_loff + l});
                    }
                }
            }
        }
    }
}

auto
add_local_matrices(CMatrices&         matrices,
                   const std::string& label,
                   const mat_t        mtype,
                   const std::string& suffix,
                   const size_t       adims,
                   const size_t       bdims,
                   const size_t       cdims,
                   const size_t       ddims) -> void
{
    // Coulomb matrices

    if ((label == "2jk") || (label == "2jkx") || (label == "j") || (label == "j_rs"))
    {
        matrices.add(CMatrix(
                         {
                             {0, 0},
                         },
                         {
                             CSubMatrix({0, 0, adims, bdims}),
                         },
                         mtype),
                     "PQ_" + suffix);

        matrices.add(CMatrix(
                         {
                             {0, 0},
                         },
                         {
                             CSubMatrix({0, 0, cdims, ddims}),
                         },
                         mtype),
                     "RS_" + suffix);

        if (mtype == mat_t::general)
        {
            matrices.add(CMatrix(
                             {
                                 {0, 0},
                             },
                             {
                                 CSubMatrix({0, 0, bdims, adims}),
                             },
                             mtype),
                         "QP_" + suffix);

            matrices.add(CMatrix(
                             {
                                 {0, 0},
                             },
                             {
                                 CSubMatrix({0, 0, ddims, cdims}),
                             },
                             mtype),
                         "SR_" + suffix);
        }
    }

    // Exchange matrices

    if ((label == "2jk") || (label == "2jkx") || (label == "k") || (label == "kx") || (label == "k_rs") || (label == "kx_rs"))
    {
        matrices.add(CMatrix(
                         {
                             {0, 0},
                         },
                         {
                             CSubMatrix({0, 0, adims, cdims}),
                         },
                         mtype),
                     "PR_" + suffix);

        matrices.add(CMatrix(
                         {
                             {0, 0},
                         },
                         {
                             CSubMatrix({0, 0, adims, ddims}),
                         },
                         mtype),
                     "PS_" + suffix);

        matrices.add(CMatrix(
                         {
                             {0, 0},
                         },
                         {
                             CSubMatrix({0, 0, bdims, cdims}),
                         },
                         mtype),
                     "QR_" + suffix);

        matrices.add(CMatrix(
                         {
                             {0, 0},
                         },
                         {
                             CSubMatrix({0, 0, bdims, ddims}),
                         },
                         mtype),
                     "QS_" + suffix);

        if (mtype == mat_t::general)
        {
            matrices.add(CMatrix(
                             {
                                 {0, 0},
                             },
                             {
                                 CSubMatrix({0, 0, cdims, adims}),
                             },
                             mtype),
                         "RP_" + suffix);

            matrices.add(CMatrix(
                             {
                                 {0, 0},
                             },
                             {
                                 CSubMatrix({0, 0, ddims, adims}),
                             },
                             mtype),
                         "SP_" + suffix);

            matrices.add(CMatrix(
                             {
                                 {0, 0},
                             },
                             {
                                 CSubMatrix({0, 0, cdims, bdims}),
                             },
                             mtype),
                         "RQ_" + suffix);

            matrices.add(CMatrix(
                             {
                                 {0, 0},
                             },
                             {
                                 CSubMatrix({0, 0, ddims, bdims}),
                             },
                             mtype),
                         "SQ_" + suffix);
        }
    }
}

auto
local_distribute(CMatrices&                       focks,
                 const std::string&               suffix,
                 const CMatrix*                   density,
                 const std::string&               label,
                 const double                     exchange_factor,
                 const CSimdArray<double>&        buffer,
                 const size_t                     offset,
                 const std::vector<size_t>&       a_indices,
                 const std::vector<size_t>&       b_indices,
                 const std::vector<size_t>&       c_indices,
                 const std::vector<size_t>&       d_indices,
                 const std::vector<size_t>&       a_loc_indices,
                 const std::vector<size_t>&       b_loc_indices,
                 const std::vector<size_t>&       c_loc_indices,
                 const std::vector<size_t>&       d_loc_indices,
                 const int                        a_angmom,
                 const int                        b_angmom,
                 const int                        c_angmom,
                 const int                        d_angmom,
                 const size_t                     bra_igto,
                 const std::pair<size_t, size_t>& ket_range,
                 const bool                       diagonal) -> void
{
    if (density->get_type() == mat_t::symmetric)
    {
        if (label == "2jk")
        {
            t4cfunc::local_distribute_rest_jk(focks,
                                              suffix,
                                              density,
                                              buffer,
                                              offset,
                                              a_indices,
                                              b_indices,
                                              c_indices,
                                              d_indices,
                                              a_loc_indices,
                                              b_loc_indices,
                                              c_loc_indices,
                                              d_loc_indices,
                                              a_angmom,
                                              b_angmom,
                                              c_angmom,
                                              d_angmom,
                                              bra_igto,
                                              ket_range,
                                              diagonal);
        }

        if (label == "2jkx")
        {
            t4cfunc::local_distribute_rest_jkx(focks,
                                               suffix,
                                               density,
                                               buffer,
                                               offset,
                                               exchange_factor,
                                               a_indices,
                                               b_indices,
                                               c_indices,
                                               d_indices,
                                               a_loc_indices,
                                               b_loc_indices,
                                               c_loc_indices,
                                               d_loc_indices,
                                               a_angmom,
                                               b_angmom,
                                               c_angmom,
                                               d_angmom,
                                               bra_igto,
                                               ket_range,
                                               diagonal);
        }

        if ((label == "j") || (label == "j_rs"))
        {
            t4cfunc::local_distribute_rest_j(focks,
                                             suffix,
                                             density,
                                             buffer,
                                             offset,
                                             a_indices,
                                             b_indices,
                                             c_indices,
                                             d_indices,
                                             a_loc_indices,
                                             b_loc_indices,
                                             c_loc_indices,
                                             d_loc_indices,
                                             a_angmom,
                                             b_angmom,
                                             c_angmom,
                                             d_angmom,
                                             bra_igto,
                                             ket_range,
                                             diagonal);
        }

        if ((label == "k") || (label == "k_rs"))
        {
            t4cfunc::local_distribute_rest_k(focks,
                                             suffix,
                                             density,
                                             buffer,
                                             offset,
                                             a_indices,
                                             b_indices,
                                             c_indices,
                                             d_indices,
                                             a_loc_indices,
                                             b_loc_indices,
                                             c_loc_indices,
                                             d_loc_indices,
                                             a_angmom,
                                             b_angmom,
                                             c_angmom,
                                             d_angmom,
                                             bra_igto,
                                             ket_range,
                                             diagonal);
        }

        if ((label == "kx") || (label == "kx_rs"))
        {
            t4cfunc::local_distribute_rest_kx(focks,
                                              suffix,
                                              density,
                                              buffer,
                                              offset,
                                              exchange_factor,
                                              a_indices,
                                              b_indices,
                                              c_indices,
                                              d_indices,
                                              a_loc_indices,
                                              b_loc_indices,
                                              c_loc_indices,
                                              d_loc_indices,
                                              a_angmom,
                                              b_angmom,
                                              c_angmom,
                                              d_angmom,
                                              bra_igto,
                                              ket_range,
                                              diagonal);
        }
    }

    if (density->get_type() == mat_t::general)
    {
        if (label == "2jk")
        {
            t4cfunc::local_distribute_gen_jk(focks,
                                             suffix,
                                             density,
                                             buffer,
                                             offset,
                                             a_indices,
                                             b_indices,
                                             c_indices,
                                             d_indices,
                                             a_loc_indices,
                                             b_loc_indices,
                                             c_loc_indices,
                                             d_loc_indices,
                                             a_angmom,
                                             b_angmom,
                                             c_angmom,
                                             d_angmom,
                                             bra_igto,
                                             ket_range,
                                             diagonal);
        }

        if (label == "2jkx")
        {
            t4cfunc::local_distribute_gen_jkx(focks,
                                              suffix,
                                              density,
                                              buffer,
                                              offset,
                                              exchange_factor,
                                              a_indices,
                                              b_indices,
                                              c_indices,
                                              d_indices,
                                              a_loc_indices,
                                              b_loc_indices,
                                              c_loc_indices,
                                              d_loc_indices,
                                              a_angmom,
                                              b_angmom,
                                              c_angmom,
                                              d_angmom,
                                              bra_igto,
                                              ket_range,
                                              diagonal);
        }

        if ((label == "j") || (label == "j_rs"))
        {
            t4cfunc::local_distribute_gen_j(focks,
                                            suffix,
                                            density,
                                            buffer,
                                            offset,
                                            a_indices,
                                            b_indices,
                                            c_indices,
                                            d_indices,
                                            a_loc_indices,
                                            b_loc_indices,
                                            c_loc_indices,
                                            d_loc_indices,
                                            a_angmom,
                                            b_angmom,
                                            c_angmom,
                                            d_angmom,
                                            bra_igto,
                                            ket_range,
                                            diagonal);
        }

        if ((label == "k") || (label == "k_rs"))
        {
            t4cfunc::local_distribute_gen_k(focks,
                                            suffix,
                                            density,
                                            buffer,
                                            offset,
                                            a_indices,
                                            b_indices,
                                            c_indices,
                                            d_indices,
                                            a_loc_indices,
                                            b_loc_indices,
                                            c_loc_indices,
                                            d_loc_indices,
                                            a_angmom,
                                            b_angmom,
                                            c_angmom,
                                            d_angmom,
                                            bra_igto,
                                            ket_range,
                                            diagonal);
        }

        if ((label == "kx") || (label == "kx_rs"))
        {
            t4cfunc::local_distribute_gen_kx(focks,
                                             suffix,
                                             density,
                                             buffer,
                                             offset,
                                             exchange_factor,
                                             a_indices,
                                             b_indices,
                                             c_indices,
                                             d_indices,
                                             a_loc_indices,
                                             b_loc_indices,
                                             c_loc_indices,
                                             d_loc_indices,
                                             a_angmom,
                                             b_angmom,
                                             c_angmom,
                                             d_angmom,
                                             bra_igto,
                                             ket_range,
                                             diagonal);
        }
    }
}

auto
local_distribute_geom_ket_symm(CMatrices&                       focks,
                               const std::string&               suffix,
                               const CMatrix*                   density,
                               const std::string&               label,
                               const double                     exchange_factor,
                               const CSimdArray<double>&        buffer,
                               const size_t                     offset,
                               const std::vector<size_t>&       a_indices,
                               const std::vector<size_t>&       b_indices,
                               const std::vector<size_t>&       c_indices,
                               const std::vector<size_t>&       d_indices,
                               const std::vector<size_t>&       a_loc_indices,
                               const std::vector<size_t>&       b_loc_indices,
                               const std::vector<size_t>&       c_loc_indices,
                               const std::vector<size_t>&       d_loc_indices,
                               const int                        a_angmom,
                               const int                        b_angmom,
                               const int                        c_angmom,
                               const int                        d_angmom,
                               const size_t                     bra_igto,
                               const std::pair<size_t, size_t>& ket_range) -> void
{
    if (label == "2jk")
    {
        if (density->get_type() == mat_t::symmetric)
        {
            t4cfunc::local_distribute_jk_geom_ket_symm_den_symm(focks,
                                                                suffix,
                                                                density,
                                                                buffer,
                                                                offset,
                                                                a_indices,
                                                                b_indices,
                                                                c_indices,
                                                                d_indices,
                                                                a_loc_indices,
                                                                b_loc_indices,
                                                                c_loc_indices,
                                                                d_loc_indices,
                                                                a_angmom,
                                                                b_angmom,
                                                                c_angmom,
                                                                d_angmom,
                                                                bra_igto,
                                                                ket_range);
        }
    }

    if (label == "2jkx")
    {
        if (density->get_type() == mat_t::symmetric)
        {
            t4cfunc::local_distribute_jkx_geom_ket_symm_den_symm(focks,
                                                                 suffix,
                                                                 density,
                                                                 buffer,
                                                                 offset,
                                                                 exchange_factor,
                                                                 a_indices,
                                                                 b_indices,
                                                                 c_indices,
                                                                 d_indices,
                                                                 a_loc_indices,
                                                                 b_loc_indices,
                                                                 c_loc_indices,
                                                                 d_loc_indices,
                                                                 a_angmom,
                                                                 b_angmom,
                                                                 c_angmom,
                                                                 d_angmom,
                                                                 bra_igto,
                                                                 ket_range);
        }
    }

    if ((label == "j") || (label == "j_rs"))
    {
        if (density->get_type() == mat_t::symmetric)
        {
            t4cfunc::local_distribute_j_geom_ket_symm_den_symm(focks,
                                                               suffix,
                                                               density,
                                                               buffer,
                                                               offset,
                                                               a_indices,
                                                               b_indices,
                                                               c_indices,
                                                               d_indices,
                                                               a_loc_indices,
                                                               b_loc_indices,
                                                               c_loc_indices,
                                                               d_loc_indices,
                                                               a_angmom,
                                                               b_angmom,
                                                               c_angmom,
                                                               d_angmom,
                                                               bra_igto,
                                                               ket_range);
        }
    }

    if ((label == "k") || (label == "k_rs"))
    {
        if (density->get_type() == mat_t::symmetric)
        {
            t4cfunc::local_distribute_k_geom_ket_symm_den_symm(focks,
                                                               suffix,
                                                               density,
                                                               buffer,
                                                               offset,
                                                               a_indices,
                                                               b_indices,
                                                               c_indices,
                                                               d_indices,
                                                               a_loc_indices,
                                                               b_loc_indices,
                                                               c_loc_indices,
                                                               d_loc_indices,
                                                               a_angmom,
                                                               b_angmom,
                                                               c_angmom,
                                                               d_angmom,
                                                               bra_igto,
                                                               ket_range);
        }
    }

    if ((label == "kx") || (label == "kx_rs"))
    {
        if (density->get_type() == mat_t::symmetric)
        {
            t4cfunc::local_distribute_kx_geom_ket_symm_den_symm(focks,
                                                                suffix,
                                                                density,
                                                                buffer,
                                                                offset,
                                                                exchange_factor,
                                                                a_indices,
                                                                b_indices,
                                                                c_indices,
                                                                d_indices,
                                                                a_loc_indices,
                                                                b_loc_indices,
                                                                c_loc_indices,
                                                                d_loc_indices,
                                                                a_angmom,
                                                                b_angmom,
                                                                c_angmom,
                                                                d_angmom,
                                                                bra_igto,
                                                                ket_range);
        }
    }
}

auto
local_distribute_geom_ket_no_symm(CMatrices&                       focks,
                                  const std::string&               suffix,
                                  const CMatrix*                   density,
                                  const std::string&               label,
                                  const double                     exchange_factor,
                                  const CSimdArray<double>&        buffer,
                                  const size_t                     offset,
                                  const std::vector<size_t>&       a_indices,
                                  const std::vector<size_t>&       b_indices,
                                  const std::vector<size_t>&       c_indices,
                                  const std::vector<size_t>&       d_indices,
                                  const std::vector<size_t>&       a_loc_indices,
                                  const std::vector<size_t>&       b_loc_indices,
                                  const std::vector<size_t>&       c_loc_indices,
                                  const std::vector<size_t>&       d_loc_indices,
                                  const int                        a_angmom,
                                  const int                        b_angmom,
                                  const int                        c_angmom,
                                  const int                        d_angmom,
                                  const size_t                     bra_igto,
                                  const std::pair<size_t, size_t>& ket_range) -> void
{
    if (label == "2jk")
    {
        if (density->get_type() == mat_t::symmetric)
        {
            t4cfunc::local_distribute_jk_geom_ket_no_symm_den_symm(focks,
                                                                   suffix,
                                                                   density,
                                                                   buffer,
                                                                   offset,
                                                                   a_indices,
                                                                   b_indices,
                                                                   c_indices,
                                                                   d_indices,
                                                                   a_loc_indices,
                                                                   b_loc_indices,
                                                                   c_loc_indices,
                                                                   d_loc_indices,
                                                                   a_angmom,
                                                                   b_angmom,
                                                                   c_angmom,
                                                                   d_angmom,
                                                                   bra_igto,
                                                                   ket_range);
        }
    }

    if (label == "2jkx")
    {
        if (density->get_type() == mat_t::symmetric)
        {
            t4cfunc::local_distribute_jkx_geom_ket_no_symm_den_symm(focks,
                                                                    suffix,
                                                                    density,
                                                                    buffer,
                                                                    offset,
                                                                    exchange_factor,
                                                                    a_indices,
                                                                    b_indices,
                                                                    c_indices,
                                                                    d_indices,
                                                                    a_loc_indices,
                                                                    b_loc_indices,
                                                                    c_loc_indices,
                                                                    d_loc_indices,
                                                                    a_angmom,
                                                                    b_angmom,
                                                                    c_angmom,
                                                                    d_angmom,
                                                                    bra_igto,
                                                                    ket_range);
        }
    }

    if ((label == "j") || (label == "j_rs"))
    {
        if (density->get_type() == mat_t::symmetric)
        {
            t4cfunc::local_distribute_j_geom_ket_no_symm_den_symm(focks,
                                                                  suffix,
                                                                  density,
                                                                  buffer,
                                                                  offset,
                                                                  a_indices,
                                                                  b_indices,
                                                                  c_indices,
                                                                  d_indices,
                                                                  a_loc_indices,
                                                                  b_loc_indices,
                                                                  c_loc_indices,
                                                                  d_loc_indices,
                                                                  a_angmom,
                                                                  b_angmom,
                                                                  c_angmom,
                                                                  d_angmom,
                                                                  bra_igto,
                                                                  ket_range);
        }
    }

    if ((label == "k") || (label == "k_rs"))
    {
        if (density->get_type() == mat_t::symmetric)
        {
            t4cfunc::local_distribute_k_geom_ket_no_symm_den_symm(focks,
                                                                  suffix,
                                                                  density,
                                                                  buffer,
                                                                  offset,
                                                                  a_indices,
                                                                  b_indices,
                                                                  c_indices,
                                                                  d_indices,
                                                                  a_loc_indices,
                                                                  b_loc_indices,
                                                                  c_loc_indices,
                                                                  d_loc_indices,
                                                                  a_angmom,
                                                                  b_angmom,
                                                                  c_angmom,
                                                                  d_angmom,
                                                                  bra_igto,
                                                                  ket_range);
        }
    }

    if ((label == "kx") || (label == "kx_rs"))
    {
        if (density->get_type() == mat_t::symmetric)
        {
            t4cfunc::local_distribute_kx_geom_ket_no_symm_den_symm(focks,
                                                                   suffix,
                                                                   density,
                                                                   buffer,
                                                                   offset,
                                                                   exchange_factor,
                                                                   a_indices,
                                                                   b_indices,
                                                                   c_indices,
                                                                   d_indices,
                                                                   a_loc_indices,
                                                                   b_loc_indices,
                                                                   c_loc_indices,
                                                                   d_loc_indices,
                                                                   a_angmom,
                                                                   b_angmom,
                                                                   c_angmom,
                                                                   d_angmom,
                                                                   bra_igto,
                                                                   ket_range);
        }
    }
}

}  // namespace t4cfunc
