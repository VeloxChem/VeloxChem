#include "T4CLocalGeomDistributor.hpp"

#include <iomanip>
#include <iostream>

#include "TensorComponents.hpp"

namespace t4cfunc {  // t4cfunc namespace

auto
local_distribute_jk_geom_ket_symm_den_symm(CMatrices&                       focks,
                                           const std::string&               suffix,
                                           const CMatrix*                   density,
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
    // set up angular pairs

    const auto pq_pair = std::pair<int, int>({a_angmom, b_angmom});

    const auto rs_pair = std::pair<int, int>({c_angmom, d_angmom});

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    // set up Fock submatrices

    auto submat_pq = focks.matrix("PQ_" + suffix)->sub_matrix({0, 0});

    auto submat_qp = focks.matrix("QP_" + suffix)->sub_matrix({0, 0});

    auto submat_rs = focks.matrix("RS_" + suffix)->sub_matrix({0, 0});

    auto submat_sr = focks.matrix("SR_" + suffix)->sub_matrix({0, 0});

    auto submat_pr = focks.matrix("PR_" + suffix)->sub_matrix({0, 0});

    auto submat_rp = focks.matrix("RP_" + suffix)->sub_matrix({0, 0});

    auto submat_ps = focks.matrix("PS_" + suffix)->sub_matrix({0, 0});

    auto submat_sp = focks.matrix("SP_" + suffix)->sub_matrix({0, 0});

    auto submat_qr = focks.matrix("QR_" + suffix)->sub_matrix({0, 0});

    auto submat_rq = focks.matrix("RQ_" + suffix)->sub_matrix({0, 0});

    auto submat_qs = focks.matrix("QS_" + suffix)->sub_matrix({0, 0});

    auto submat_sq = focks.matrix("SQ_" + suffix)->sub_matrix({0, 0});

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    // set up angular orders

    const auto angord_pq = density->is_angular_order(pq_pair);

    const auto angord_rs = density->is_angular_order(rs_pair);

    const auto angord_pr = density->is_angular_order(pr_pair);

    const auto angord_ps = density->is_angular_order(ps_pair);

    const auto angord_qr = density->is_angular_order(qr_pair);

    const auto angord_qs = density->is_angular_order(qs_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // reference local indexes on bra side

    const auto loc_refp = a_loc_indices[bra_igto + 1];

    const auto loc_refq = b_loc_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // dimensions of local bra and ket orbital indexes

    const auto alocdim = a_loc_indices[0];

    const auto blocdim = b_loc_indices[0];

    const auto clocdim = c_loc_indices[0];

    const auto dlocdim = d_loc_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        const auto loc_p = i * alocdim + loc_refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            const auto loc_q = j * blocdim + loc_refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // reference local indexes on ket side

                        const auto loc_refr = c_loc_indices[m + 1];

                        const auto loc_refs = d_loc_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto loc_r = k * clocdim + loc_refr;

                        const auto s = l * ddim + refs;

                        const auto loc_s = l * dlocdim + loc_refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        if (angord_rs)
                        {
                            submat_pq->at({loc_p, loc_q}) += 4.0 * fval * denmat_rs->operator[]({r, s});

                            submat_qp->at({loc_q, loc_p}) += 4.0 * fval * denmat_rs->operator[]({r, s});
                        }
                        else
                        {
                            submat_pq->at({loc_p, loc_q}) += 4.0 * fval * denmat_rs->operator[]({s, r});

                            submat_qp->at({loc_q, loc_p}) += 4.0 * fval * denmat_rs->operator[]({s, r});
                        }

                        if (angord_pq)
                        {
                            submat_rs->at({loc_r, loc_s}) += 4.0 * fval * denmat_pq->operator[]({p, q});

                            submat_sr->at({loc_s, loc_r}) += 4.0 * fval * denmat_pq->operator[]({p, q});
                        }
                        else
                        {
                            submat_rs->at({loc_r, loc_s}) += 4.0 * fval * denmat_pq->operator[]({q, p});

                            submat_sr->at({loc_s, loc_r}) += 4.0 * fval * denmat_pq->operator[]({q, p});
                        }

                        // Exchange contribution (F_pr)

                        if (angord_qs)
                        {
                            submat_pr->at({loc_p, loc_r}) -= fval * denmat_qs->operator[]({q, s});

                            submat_rp->at({loc_r, loc_p}) -= fval * denmat_qs->operator[]({q, s});
                        }
                        else
                        {
                            submat_pr->at({loc_p, loc_r}) -= fval * denmat_qs->operator[]({s, q});

                            submat_rp->at({loc_r, loc_p}) -= fval * denmat_qs->operator[]({s, q});
                        }

                        // Exchange contribution (F_ps)

                        if (angord_qr)
                        {
                            submat_ps->at({loc_p, loc_s}) -= fval * denmat_qr->operator[]({q, r});

                            submat_sp->at({loc_s, loc_p}) -= fval * denmat_qr->operator[]({q, r});
                        }
                        else
                        {
                            submat_ps->at({loc_p, loc_s}) -= fval * denmat_qr->operator[]({r, q});

                            submat_sp->at({loc_s, loc_p}) -= fval * denmat_qr->operator[]({r, q});
                        }

                        // Exchange contribution (F_qr)

                        if (angord_ps)
                        {
                            submat_qr->at({loc_q, loc_r}) -= fval * denmat_ps->operator[]({p, s});

                            submat_rq->at({loc_r, loc_q}) -= fval * denmat_ps->operator[]({p, s});
                        }
                        else
                        {
                            submat_qr->at({loc_q, loc_r}) -= fval * denmat_ps->operator[]({s, p});

                            submat_rq->at({loc_r, loc_q}) -= fval * denmat_ps->operator[]({s, p});
                        }

                        // Exchange contribution (F_qs)

                        if (angord_pr)
                        {
                            submat_qs->at({loc_q, loc_s}) -= fval * denmat_pr->operator[]({p, r});

                            submat_sq->at({loc_s, loc_q}) -= fval * denmat_pr->operator[]({p, r});
                        }
                        else
                        {
                            submat_qs->at({loc_q, loc_s}) -= fval * denmat_pr->operator[]({r, p});

                            submat_sq->at({loc_s, loc_q}) -= fval * denmat_pr->operator[]({r, p});
                        }
                    }
                }
            }
        }
    }
}

auto
local_distribute_jkx_geom_ket_symm_den_symm(CMatrices&                       focks,
                                            const std::string&               suffix,
                                            const CMatrix*                   density,
                                            const CSimdArray<double>&        buffer,
                                            const size_t                     offset,
                                            const double                     factor,
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
    // set up angular pairs

    const auto pq_pair = std::pair<int, int>({a_angmom, b_angmom});

    const auto rs_pair = std::pair<int, int>({c_angmom, d_angmom});

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    // set up Fock submatrices

    auto submat_pq = focks.matrix("PQ_" + suffix)->sub_matrix({0, 0});

    auto submat_qp = focks.matrix("QP_" + suffix)->sub_matrix({0, 0});

    auto submat_rs = focks.matrix("RS_" + suffix)->sub_matrix({0, 0});

    auto submat_sr = focks.matrix("SR_" + suffix)->sub_matrix({0, 0});

    auto submat_pr = focks.matrix("PR_" + suffix)->sub_matrix({0, 0});

    auto submat_rp = focks.matrix("RP_" + suffix)->sub_matrix({0, 0});

    auto submat_ps = focks.matrix("PS_" + suffix)->sub_matrix({0, 0});

    auto submat_sp = focks.matrix("SP_" + suffix)->sub_matrix({0, 0});

    auto submat_qr = focks.matrix("QR_" + suffix)->sub_matrix({0, 0});

    auto submat_rq = focks.matrix("RQ_" + suffix)->sub_matrix({0, 0});

    auto submat_qs = focks.matrix("QS_" + suffix)->sub_matrix({0, 0});

    auto submat_sq = focks.matrix("SQ_" + suffix)->sub_matrix({0, 0});

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    // set up angular orders

    const auto angord_pq = density->is_angular_order(pq_pair);

    const auto angord_rs = density->is_angular_order(rs_pair);

    const auto angord_pr = density->is_angular_order(pr_pair);

    const auto angord_ps = density->is_angular_order(ps_pair);

    const auto angord_qr = density->is_angular_order(qr_pair);

    const auto angord_qs = density->is_angular_order(qs_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // reference local indexes on bra side

    const auto loc_refp = a_loc_indices[bra_igto + 1];

    const auto loc_refq = b_loc_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // dimensions of local bra and ket orbital indexes

    const auto alocdim = a_loc_indices[0];

    const auto blocdim = b_loc_indices[0];

    const auto clocdim = c_loc_indices[0];

    const auto dlocdim = d_loc_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        const auto loc_p = i * alocdim + loc_refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            const auto loc_q = j * blocdim + loc_refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // reference local indexes on ket side

                        const auto loc_refr = c_loc_indices[m + 1];

                        const auto loc_refs = d_loc_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto loc_r = k * clocdim + loc_refr;

                        const auto s = l * ddim + refs;

                        const auto loc_s = l * dlocdim + loc_refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        if (angord_rs)
                        {
                            submat_pq->at({loc_p, loc_q}) += 4.0 * fval * denmat_rs->operator[]({r, s});

                            submat_qp->at({loc_q, loc_p}) += 4.0 * fval * denmat_rs->operator[]({r, s});
                        }
                        else
                        {
                            submat_pq->at({loc_p, loc_q}) += 4.0 * fval * denmat_rs->operator[]({s, r});

                            submat_qp->at({loc_q, loc_p}) += 4.0 * fval * denmat_rs->operator[]({s, r});
                        }

                        if (angord_pq)
                        {
                            submat_rs->at({loc_r, loc_s}) += 4.0 * fval * denmat_pq->operator[]({p, q});

                            submat_sr->at({loc_s, loc_r}) += 4.0 * fval * denmat_pq->operator[]({p, q});
                        }
                        else
                        {
                            submat_rs->at({loc_r, loc_s}) += 4.0 * fval * denmat_pq->operator[]({q, p});

                            submat_sr->at({loc_s, loc_r}) += 4.0 * fval * denmat_pq->operator[]({q, p});
                        }

                        // Exchange contribution (F_pr)

                        fval *= factor;

                        if (angord_qs)
                        {
                            submat_pr->at({loc_p, loc_r}) -= fval * denmat_qs->operator[]({q, s});

                            submat_rp->at({loc_r, loc_p}) -= fval * denmat_qs->operator[]({q, s});
                        }
                        else
                        {
                            submat_pr->at({loc_p, loc_r}) -= fval * denmat_qs->operator[]({s, q});

                            submat_rp->at({loc_r, loc_p}) -= fval * denmat_qs->operator[]({s, q});
                        }

                        // Exchange contribution (F_ps)

                        if (angord_qr)
                        {
                            submat_ps->at({loc_p, loc_s}) -= fval * denmat_qr->operator[]({q, r});

                            submat_sp->at({loc_s, loc_p}) -= fval * denmat_qr->operator[]({q, r});
                        }
                        else
                        {
                            submat_ps->at({loc_p, loc_s}) -= fval * denmat_qr->operator[]({r, q});

                            submat_sp->at({loc_s, loc_p}) -= fval * denmat_qr->operator[]({r, q});
                        }

                        // Exchange contribution (F_qr)

                        if (angord_ps)
                        {
                            submat_qr->at({loc_q, loc_r}) -= fval * denmat_ps->operator[]({p, s});

                            submat_rq->at({loc_r, loc_q}) -= fval * denmat_ps->operator[]({p, s});
                        }
                        else
                        {
                            submat_qr->at({loc_q, loc_r}) -= fval * denmat_ps->operator[]({s, p});

                            submat_rq->at({loc_r, loc_q}) -= fval * denmat_ps->operator[]({s, p});
                        }

                        // Exchange contribution (F_qs)

                        if (angord_pr)
                        {
                            submat_qs->at({loc_q, loc_s}) -= fval * denmat_pr->operator[]({p, r});

                            submat_sq->at({loc_s, loc_q}) -= fval * denmat_pr->operator[]({p, r});
                        }
                        else
                        {
                            submat_qs->at({loc_q, loc_s}) -= fval * denmat_pr->operator[]({r, p});

                            submat_sq->at({loc_s, loc_q}) -= fval * denmat_pr->operator[]({r, p});
                        }
                    }
                }
            }
        }
    }
}

auto
local_distribute_j_geom_ket_symm_den_symm(CMatrices&                       focks,
                                          const std::string&               suffix,
                                          const CMatrix*                   density,
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
    // set up angular pairs

    const auto pq_pair = std::pair<int, int>({a_angmom, b_angmom});

    const auto rs_pair = std::pair<int, int>({c_angmom, d_angmom});

    // set up Fock submatrices

    auto submat_pq = focks.matrix("PQ_" + suffix)->sub_matrix({0, 0});

    auto submat_qp = focks.matrix("QP_" + suffix)->sub_matrix({0, 0});

    auto submat_rs = focks.matrix("RS_" + suffix)->sub_matrix({0, 0});

    auto submat_sr = focks.matrix("SR_" + suffix)->sub_matrix({0, 0});

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    // set up angular orders

    const auto angord_pq = density->is_angular_order(pq_pair);

    const auto angord_rs = density->is_angular_order(rs_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // reference local indexes on bra side

    const auto loc_refp = a_loc_indices[bra_igto + 1];

    const auto loc_refq = b_loc_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // dimensions of local bra and ket orbital indexes

    const auto alocdim = a_loc_indices[0];

    const auto blocdim = b_loc_indices[0];

    const auto clocdim = c_loc_indices[0];

    const auto dlocdim = d_loc_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        const auto loc_p = i * alocdim + loc_refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            const auto loc_q = j * blocdim + loc_refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // reference local indexes on ket side

                        const auto loc_refr = c_loc_indices[m + 1];

                        const auto loc_refs = d_loc_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto loc_r = k * clocdim + loc_refr;

                        const auto s = l * ddim + refs;

                        const auto loc_s = l * dlocdim + loc_refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        if (angord_rs)
                        {
                            submat_pq->at({loc_p, loc_q}) += 2.0 * fval * denmat_rs->operator[]({r, s});

                            submat_qp->at({loc_q, loc_p}) += 2.0 * fval * denmat_rs->operator[]({r, s});
                        }
                        else
                        {
                            submat_pq->at({loc_p, loc_q}) += 2.0 * fval * denmat_rs->operator[]({s, r});

                            submat_qp->at({loc_q, loc_p}) += 2.0 * fval * denmat_rs->operator[]({s, r});
                        }

                        if (angord_pq)
                        {
                            submat_rs->at({loc_r, loc_s}) += 2.0 * fval * denmat_pq->operator[]({p, q});

                            submat_sr->at({loc_s, loc_r}) += 2.0 * fval * denmat_pq->operator[]({p, q});
                        }
                        else
                        {
                            submat_rs->at({loc_r, loc_s}) += 2.0 * fval * denmat_pq->operator[]({q, p});

                            submat_sr->at({loc_s, loc_r}) += 2.0 * fval * denmat_pq->operator[]({q, p});
                        }
                    }
                }
            }
        }
    }
}

auto
local_distribute_k_geom_ket_symm_den_symm(CMatrices&                       focks,
                                          const std::string&               suffix,
                                          const CMatrix*                   density,
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
    // set up angular pairs

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    // set up Fock submatrices

    auto submat_pr = focks.matrix("PR_" + suffix)->sub_matrix({0, 0});

    auto submat_rp = focks.matrix("RP_" + suffix)->sub_matrix({0, 0});

    auto submat_ps = focks.matrix("PS_" + suffix)->sub_matrix({0, 0});

    auto submat_sp = focks.matrix("SP_" + suffix)->sub_matrix({0, 0});

    auto submat_qr = focks.matrix("QR_" + suffix)->sub_matrix({0, 0});

    auto submat_rq = focks.matrix("RQ_" + suffix)->sub_matrix({0, 0});

    auto submat_qs = focks.matrix("QS_" + suffix)->sub_matrix({0, 0});

    auto submat_sq = focks.matrix("SQ_" + suffix)->sub_matrix({0, 0});

    // set up AO density submatrix

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    // set up angular orders

    const auto angord_pr = density->is_angular_order(pr_pair);

    const auto angord_ps = density->is_angular_order(ps_pair);

    const auto angord_qr = density->is_angular_order(qr_pair);

    const auto angord_qs = density->is_angular_order(qs_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // reference local indexes on bra side

    const auto loc_refp = a_loc_indices[bra_igto + 1];

    const auto loc_refq = b_loc_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // dimensions of local bra and ket orbital indexes

    const auto alocdim = a_loc_indices[0];

    const auto blocdim = b_loc_indices[0];

    const auto clocdim = c_loc_indices[0];

    const auto dlocdim = d_loc_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        const auto loc_p = i * alocdim + loc_refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            const auto loc_q = j * blocdim + loc_refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // reference local indexes on ket side

                        const auto loc_refr = c_loc_indices[m + 1];

                        const auto loc_refs = d_loc_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto loc_r = k * clocdim + loc_refr;

                        const auto s = l * ddim + refs;

                        const auto loc_s = l * dlocdim + loc_refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Exchange contribution (F_pr)

                        if (angord_qs)
                        {
                            submat_pr->at({loc_p, loc_r}) += fval * denmat_qs->operator[]({q, s});

                            submat_rp->at({loc_r, loc_p}) += fval * denmat_qs->operator[]({q, s});
                        }
                        else
                        {
                            submat_pr->at({loc_p, loc_r}) += fval * denmat_qs->operator[]({s, q});

                            submat_rp->at({loc_r, loc_p}) += fval * denmat_qs->operator[]({s, q});
                        }

                        // Exchange contribution (F_ps)

                        if (angord_qr)
                        {
                            submat_ps->at({loc_p, loc_s}) += fval * denmat_qr->operator[]({q, r});

                            submat_sp->at({loc_s, loc_p}) += fval * denmat_qr->operator[]({q, r});
                        }
                        else
                        {
                            submat_ps->at({loc_p, loc_s}) += fval * denmat_qr->operator[]({r, q});

                            submat_sp->at({loc_s, loc_p}) += fval * denmat_qr->operator[]({r, q});
                        }

                        // Exchange contribution (F_qr)

                        if (angord_ps)
                        {
                            submat_qr->at({loc_q, loc_r}) += fval * denmat_ps->operator[]({p, s});

                            submat_rq->at({loc_r, loc_q}) += fval * denmat_ps->operator[]({p, s});
                        }
                        else
                        {
                            submat_qr->at({loc_q, loc_r}) += fval * denmat_ps->operator[]({s, p});

                            submat_rq->at({loc_r, loc_q}) += fval * denmat_ps->operator[]({s, p});
                        }

                        // Exchange contribution (F_qs)

                        if (angord_pr)
                        {
                            submat_qs->at({loc_q, loc_s}) += fval * denmat_pr->operator[]({p, r});

                            submat_sq->at({loc_s, loc_q}) += fval * denmat_pr->operator[]({p, r});
                        }
                        else
                        {
                            submat_qs->at({loc_q, loc_s}) += fval * denmat_pr->operator[]({r, p});

                            submat_sq->at({loc_s, loc_q}) += fval * denmat_pr->operator[]({r, p});
                        }
                    }
                }
            }
        }
    }
}

auto
local_distribute_kx_geom_ket_symm_den_symm(CMatrices&                       focks,
                                           const std::string&               suffix,
                                           const CMatrix*                   density,
                                           const CSimdArray<double>&        buffer,
                                           const size_t                     offset,
                                           const double                     factor,
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
    // set up angular pairs

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    // set up Fock submatrices

    auto submat_pr = focks.matrix("PR_" + suffix)->sub_matrix({0, 0});

    auto submat_rp = focks.matrix("RP_" + suffix)->sub_matrix({0, 0});

    auto submat_ps = focks.matrix("PS_" + suffix)->sub_matrix({0, 0});

    auto submat_sp = focks.matrix("SP_" + suffix)->sub_matrix({0, 0});

    auto submat_qr = focks.matrix("QR_" + suffix)->sub_matrix({0, 0});

    auto submat_rq = focks.matrix("RQ_" + suffix)->sub_matrix({0, 0});

    auto submat_qs = focks.matrix("QS_" + suffix)->sub_matrix({0, 0});

    auto submat_sq = focks.matrix("SQ_" + suffix)->sub_matrix({0, 0});

    // set up AO density submatrix

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    // set up angular orders

    const auto angord_pr = density->is_angular_order(pr_pair);

    const auto angord_ps = density->is_angular_order(ps_pair);

    const auto angord_qr = density->is_angular_order(qr_pair);

    const auto angord_qs = density->is_angular_order(qs_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // reference local indexes on bra side

    const auto loc_refp = a_loc_indices[bra_igto + 1];

    const auto loc_refq = b_loc_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // dimensions of local bra and ket orbital indexes

    const auto alocdim = a_loc_indices[0];

    const auto blocdim = b_loc_indices[0];

    const auto clocdim = c_loc_indices[0];

    const auto dlocdim = d_loc_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        const auto loc_p = i * alocdim + loc_refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            const auto loc_q = j * blocdim + loc_refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // reference local indexes on ket side

                        const auto loc_refr = c_loc_indices[m + 1];

                        const auto loc_refs = d_loc_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto loc_r = k * clocdim + loc_refr;

                        const auto s = l * ddim + refs;

                        const auto loc_s = l * dlocdim + loc_refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = factor * curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Exchange contribution (F_pr)

                        if (angord_qs)
                        {
                            submat_pr->at({loc_p, loc_r}) += fval * denmat_qs->operator[]({q, s});

                            submat_rp->at({loc_r, loc_p}) += fval * denmat_qs->operator[]({q, s});
                        }
                        else
                        {
                            submat_pr->at({loc_p, loc_r}) += fval * denmat_qs->operator[]({s, q});

                            submat_rp->at({loc_r, loc_p}) += fval * denmat_qs->operator[]({s, q});
                        }

                        // Exchange contribution (F_ps)

                        if (angord_qr)
                        {
                            submat_ps->at({loc_p, loc_s}) += fval * denmat_qr->operator[]({q, r});

                            submat_sp->at({loc_s, loc_p}) += fval * denmat_qr->operator[]({q, r});
                        }
                        else
                        {
                            submat_ps->at({loc_p, loc_s}) += fval * denmat_qr->operator[]({r, q});

                            submat_sp->at({loc_s, loc_p}) += fval * denmat_qr->operator[]({r, q});
                        }

                        // Exchange contribution (F_qr)

                        if (angord_ps)
                        {
                            submat_qr->at({loc_q, loc_r}) += fval * denmat_ps->operator[]({p, s});

                            submat_rq->at({loc_r, loc_q}) += fval * denmat_ps->operator[]({p, s});
                        }
                        else
                        {
                            submat_qr->at({loc_q, loc_r}) += fval * denmat_ps->operator[]({s, p});

                            submat_rq->at({loc_r, loc_q}) += fval * denmat_ps->operator[]({s, p});
                        }

                        // Exchange contribution (F_qs)

                        if (angord_pr)
                        {
                            submat_qs->at({loc_q, loc_s}) += fval * denmat_pr->operator[]({p, r});

                            submat_sq->at({loc_s, loc_q}) += fval * denmat_pr->operator[]({p, r});
                        }
                        else
                        {
                            submat_qs->at({loc_q, loc_s}) += fval * denmat_pr->operator[]({r, p});

                            submat_sq->at({loc_s, loc_q}) += fval * denmat_pr->operator[]({r, p});
                        }
                    }
                }
            }
        }
    }
}

auto
local_distribute_jk_geom_ket_gen(CMatrices&                       focks,
                                 const std::string&               suffix,
                                 const CMatrix*                   density,
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
    // set up angular pairs

    const auto pq_pair = std::pair<int, int>({a_angmom, b_angmom});

    const auto qp_pair = std::pair<int, int>({b_angmom, a_angmom});

    const auto rs_pair = std::pair<int, int>({c_angmom, d_angmom});

    const auto sr_pair = std::pair<int, int>({d_angmom, c_angmom});

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto rp_pair = std::pair<int, int>({c_angmom, a_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto sp_pair = std::pair<int, int>({d_angmom, a_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto rq_pair = std::pair<int, int>({c_angmom, b_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    const auto sq_pair = std::pair<int, int>({d_angmom, b_angmom});

    // set up Fock submatrices

    auto submat_pq = focks.matrix("PQ_" + suffix)->sub_matrix({0, 0});

    auto submat_qp = focks.matrix("QP_" + suffix)->sub_matrix({0, 0});

    auto submat_rs = focks.matrix("RS_" + suffix)->sub_matrix({0, 0});

    auto submat_sr = focks.matrix("SR_" + suffix)->sub_matrix({0, 0});

    auto submat_pr = focks.matrix("PR_" + suffix)->sub_matrix({0, 0});

    auto submat_rp = focks.matrix("RP_" + suffix)->sub_matrix({0, 0});

    auto submat_ps = focks.matrix("PS_" + suffix)->sub_matrix({0, 0});

    auto submat_sp = focks.matrix("SP_" + suffix)->sub_matrix({0, 0});

    auto submat_qr = focks.matrix("QR_" + suffix)->sub_matrix({0, 0});

    auto submat_rq = focks.matrix("RQ_" + suffix)->sub_matrix({0, 0});

    auto submat_qs = focks.matrix("QS_" + suffix)->sub_matrix({0, 0});

    auto submat_sq = focks.matrix("SQ_" + suffix)->sub_matrix({0, 0});

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_qp = density->sub_matrix(qp_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    auto denmat_sr = density->sub_matrix(sr_pair);

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_rp = density->sub_matrix(rp_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_sp = density->sub_matrix(sp_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_rq = density->sub_matrix(rq_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    auto denmat_sq = density->sub_matrix(sq_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // reference local indexes on bra side

    const auto loc_refp = a_loc_indices[bra_igto + 1];

    const auto loc_refq = b_loc_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // dimensions of local bra and ket orbital indexes

    const auto alocdim = a_loc_indices[0];

    const auto blocdim = b_loc_indices[0];

    const auto clocdim = c_loc_indices[0];

    const auto dlocdim = d_loc_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        const auto loc_p = i * alocdim + loc_refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            const auto loc_q = j * blocdim + loc_refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // reference local indexes on ket side

                        const auto loc_refr = c_loc_indices[m + 1];

                        const auto loc_refs = d_loc_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto loc_r = k * clocdim + loc_refr;

                        const auto s = l * ddim + refs;

                        const auto loc_s = l * dlocdim + loc_refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        const auto f2rs = 2.0 * fval * (denmat_rs->operator[]({r, s}) + denmat_sr->operator[]({s, r}));

                        submat_pq->at({loc_p, loc_q}) += f2rs;
                        submat_qp->at({loc_q, loc_p}) += f2rs;

                        // Coulomb contribution (F_rs)

                        const auto f2pq = 2.0 * fval * (denmat_pq->operator[]({p, q}) + denmat_qp->operator[]({q, p}));

                        submat_rs->at({loc_r, loc_s}) += f2pq;
                        submat_sr->at({loc_s, loc_r}) += f2pq;

                        // Exchange contribution (F_pr)

                        submat_pr->at({loc_p, loc_r}) -= fval * denmat_qs->operator[]({q, s});
                        submat_rp->at({loc_r, loc_p}) -= fval * denmat_sq->operator[]({s, q});

                        // Exchange contribution (F_ps)

                        submat_ps->at({loc_p, loc_s}) -= fval * denmat_qr->operator[]({q, r});
                        submat_sp->at({loc_s, loc_p}) -= fval * denmat_rq->operator[]({r, q});

                        // Exchange contribution (F_qr)

                        submat_qr->at({loc_q, loc_r}) -= fval * denmat_ps->operator[]({p, s});
                        submat_rq->at({loc_r, loc_q}) -= fval * denmat_sp->operator[]({s, p});

                        // Exchange contribution (F_qs)

                        submat_qs->at({loc_q, loc_s}) -= fval * denmat_pr->operator[]({p, r});
                        submat_sq->at({loc_s, loc_q}) -= fval * denmat_rp->operator[]({r, p});
                    }
                }
            }
        }
    }
}

auto
local_distribute_jkx_geom_ket_gen(CMatrices&                       focks,
                                  const std::string&               suffix,
                                  const CMatrix*                   density,
                                  const CSimdArray<double>&        buffer,
                                  const size_t                     offset,
                                  const double                     factor,
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
    // set up angular pairs

    const auto pq_pair = std::pair<int, int>({a_angmom, b_angmom});

    const auto qp_pair = std::pair<int, int>({b_angmom, a_angmom});

    const auto rs_pair = std::pair<int, int>({c_angmom, d_angmom});

    const auto sr_pair = std::pair<int, int>({d_angmom, c_angmom});

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto rp_pair = std::pair<int, int>({c_angmom, a_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto sp_pair = std::pair<int, int>({d_angmom, a_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto rq_pair = std::pair<int, int>({c_angmom, b_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    const auto sq_pair = std::pair<int, int>({d_angmom, b_angmom});

    // set up Fock submatrices

    auto submat_pq = focks.matrix("PQ_" + suffix)->sub_matrix({0, 0});

    auto submat_qp = focks.matrix("QP_" + suffix)->sub_matrix({0, 0});

    auto submat_rs = focks.matrix("RS_" + suffix)->sub_matrix({0, 0});

    auto submat_sr = focks.matrix("SR_" + suffix)->sub_matrix({0, 0});

    auto submat_pr = focks.matrix("PR_" + suffix)->sub_matrix({0, 0});

    auto submat_rp = focks.matrix("RP_" + suffix)->sub_matrix({0, 0});

    auto submat_ps = focks.matrix("PS_" + suffix)->sub_matrix({0, 0});

    auto submat_sp = focks.matrix("SP_" + suffix)->sub_matrix({0, 0});

    auto submat_qr = focks.matrix("QR_" + suffix)->sub_matrix({0, 0});

    auto submat_rq = focks.matrix("RQ_" + suffix)->sub_matrix({0, 0});

    auto submat_qs = focks.matrix("QS_" + suffix)->sub_matrix({0, 0});

    auto submat_sq = focks.matrix("SQ_" + suffix)->sub_matrix({0, 0});

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_qp = density->sub_matrix(qp_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    auto denmat_sr = density->sub_matrix(sr_pair);

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_rp = density->sub_matrix(rp_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_sp = density->sub_matrix(sp_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_rq = density->sub_matrix(rq_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    auto denmat_sq = density->sub_matrix(sq_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // reference local indexes on bra side

    const auto loc_refp = a_loc_indices[bra_igto + 1];

    const auto loc_refq = b_loc_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // dimensions of local bra and ket orbital indexes

    const auto alocdim = a_loc_indices[0];

    const auto blocdim = b_loc_indices[0];

    const auto clocdim = c_loc_indices[0];

    const auto dlocdim = d_loc_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        const auto loc_p = i * alocdim + loc_refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            const auto loc_q = j * blocdim + loc_refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // reference local indexes on ket side

                        const auto loc_refr = c_loc_indices[m + 1];

                        const auto loc_refs = d_loc_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto loc_r = k * clocdim + loc_refr;

                        const auto s = l * ddim + refs;

                        const auto loc_s = l * dlocdim + loc_refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        const auto f2rs = 2.0 * fval * (denmat_rs->operator[]({r, s}) + denmat_sr->operator[]({s, r}));

                        submat_pq->at({loc_p, loc_q}) += f2rs;
                        submat_qp->at({loc_q, loc_p}) += f2rs;

                        // Coulomb contribution (F_rs)

                        const auto f2pq = 2.0 * fval * (denmat_pq->operator[]({p, q}) + denmat_qp->operator[]({q, p}));

                        submat_rs->at({loc_r, loc_s}) += f2pq;
                        submat_sr->at({loc_s, loc_r}) += f2pq;

                        // Exchange contribution

                        fval *= factor;

                        // Exchange contribution (F_pr)

                        submat_pr->at({loc_p, loc_r}) -= fval * denmat_qs->operator[]({q, s});
                        submat_rp->at({loc_r, loc_p}) -= fval * denmat_sq->operator[]({s, q});

                        // Exchange contribution (F_ps)

                        submat_ps->at({loc_p, loc_s}) -= fval * denmat_qr->operator[]({q, r});
                        submat_sp->at({loc_s, loc_p}) -= fval * denmat_rq->operator[]({r, q});

                        // Exchange contribution (F_qr)

                        submat_qr->at({loc_q, loc_r}) -= fval * denmat_ps->operator[]({p, s});
                        submat_rq->at({loc_r, loc_q}) -= fval * denmat_sp->operator[]({s, p});

                        // Exchange contribution (F_qs)

                        submat_qs->at({loc_q, loc_s}) -= fval * denmat_pr->operator[]({p, r});
                        submat_sq->at({loc_s, loc_q}) -= fval * denmat_rp->operator[]({r, p});
                    }
                }
            }
        }
    }
}

auto
local_distribute_j_geom_ket_gen(CMatrices&                       focks,
                                const std::string&               suffix,
                                const CMatrix*                   density,
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
    // set up angular pairs

    const auto pq_pair = std::pair<int, int>({a_angmom, b_angmom});

    const auto qp_pair = std::pair<int, int>({b_angmom, a_angmom});

    const auto rs_pair = std::pair<int, int>({c_angmom, d_angmom});

    const auto sr_pair = std::pair<int, int>({d_angmom, c_angmom});

    // set up Fock submatrices

    auto submat_pq = focks.matrix("PQ_" + suffix)->sub_matrix({0, 0});

    auto submat_qp = focks.matrix("QP_" + suffix)->sub_matrix({0, 0});

    auto submat_rs = focks.matrix("RS_" + suffix)->sub_matrix({0, 0});

    auto submat_sr = focks.matrix("SR_" + suffix)->sub_matrix({0, 0});

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_qp = density->sub_matrix(qp_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    auto denmat_sr = density->sub_matrix(sr_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // reference local indexes on bra side

    const auto loc_refp = a_loc_indices[bra_igto + 1];

    const auto loc_refq = b_loc_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // dimensions of local bra and ket orbital indexes

    const auto alocdim = a_loc_indices[0];

    const auto blocdim = b_loc_indices[0];

    const auto clocdim = c_loc_indices[0];

    const auto dlocdim = d_loc_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        const auto loc_p = i * alocdim + loc_refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            const auto loc_q = j * blocdim + loc_refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // reference local indexes on ket side

                        const auto loc_refr = c_loc_indices[m + 1];

                        const auto loc_refs = d_loc_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto loc_r = k * clocdim + loc_refr;

                        const auto s = l * ddim + refs;

                        const auto loc_s = l * dlocdim + loc_refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        const auto f_rs = fval * (denmat_rs->operator[]({r, s}) + denmat_sr->operator[]({s, r}));

                        submat_pq->at({loc_p, loc_q}) += f_rs;
                        submat_qp->at({loc_q, loc_p}) += f_rs;

                        // Coulomb contribution (F_rs)

                        const auto f_pq = fval * (denmat_pq->operator[]({p, q}) + denmat_qp->operator[]({q, p}));

                        submat_rs->at({loc_r, loc_s}) += f_pq;
                        submat_sr->at({loc_s, loc_r}) += f_pq;
                    }
                }
            }
        }
    }
}

auto
local_distribute_k_geom_ket_gen(CMatrices&                       focks,
                                const std::string&               suffix,
                                const CMatrix*                   density,
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
    // set up angular pairs

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto rp_pair = std::pair<int, int>({c_angmom, a_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto sp_pair = std::pair<int, int>({d_angmom, a_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto rq_pair = std::pair<int, int>({c_angmom, b_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    const auto sq_pair = std::pair<int, int>({d_angmom, b_angmom});

    // set up Fock submatrices

    auto submat_pr = focks.matrix("PR_" + suffix)->sub_matrix({0, 0});

    auto submat_rp = focks.matrix("RP_" + suffix)->sub_matrix({0, 0});

    auto submat_ps = focks.matrix("PS_" + suffix)->sub_matrix({0, 0});

    auto submat_sp = focks.matrix("SP_" + suffix)->sub_matrix({0, 0});

    auto submat_qr = focks.matrix("QR_" + suffix)->sub_matrix({0, 0});

    auto submat_rq = focks.matrix("RQ_" + suffix)->sub_matrix({0, 0});

    auto submat_qs = focks.matrix("QS_" + suffix)->sub_matrix({0, 0});

    auto submat_sq = focks.matrix("SQ_" + suffix)->sub_matrix({0, 0});

    // set up AO density submatrix

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_rp = density->sub_matrix(rp_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_sp = density->sub_matrix(sp_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_rq = density->sub_matrix(rq_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    auto denmat_sq = density->sub_matrix(sq_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // reference local indexes on bra side

    const auto loc_refp = a_loc_indices[bra_igto + 1];

    const auto loc_refq = b_loc_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // dimensions of local bra and ket orbital indexes

    const auto alocdim = a_loc_indices[0];

    const auto blocdim = b_loc_indices[0];

    const auto clocdim = c_loc_indices[0];

    const auto dlocdim = d_loc_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        const auto loc_p = i * alocdim + loc_refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            const auto loc_q = j * blocdim + loc_refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // reference local indexes on ket side

                        const auto loc_refr = c_loc_indices[m + 1];

                        const auto loc_refs = d_loc_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto loc_r = k * clocdim + loc_refr;

                        const auto s = l * ddim + refs;

                        const auto loc_s = l * dlocdim + loc_refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Exchange contribution (F_pr)

                        submat_pr->at({loc_p, loc_r}) += fval * denmat_qs->operator[]({q, s});
                        submat_rp->at({loc_r, loc_p}) += fval * denmat_sq->operator[]({s, q});

                        // Exchange contribution (F_ps)

                        submat_ps->at({loc_p, loc_s}) += fval * denmat_qr->operator[]({q, r});
                        submat_sp->at({loc_s, loc_p}) += fval * denmat_rq->operator[]({r, q});

                        // Exchange contribution (F_qr)

                        submat_qr->at({loc_q, loc_r}) += fval * denmat_ps->operator[]({p, s});
                        submat_rq->at({loc_r, loc_q}) += fval * denmat_sp->operator[]({s, p});

                        // Exchange contribution (F_qs)

                        submat_qs->at({loc_q, loc_s}) += fval * denmat_pr->operator[]({p, r});
                        submat_sq->at({loc_s, loc_q}) += fval * denmat_rp->operator[]({r, p});
                    }
                }
            }
        }
    }
}

auto
local_distribute_kx_geom_ket_gen(CMatrices&                       focks,
                                 const std::string&               suffix,
                                 const CMatrix*                   density,
                                 const CSimdArray<double>&        buffer,
                                 const size_t                     offset,
                                 const double                     factor,
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
    // set up angular pairs

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto rp_pair = std::pair<int, int>({c_angmom, a_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto sp_pair = std::pair<int, int>({d_angmom, a_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto rq_pair = std::pair<int, int>({c_angmom, b_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    const auto sq_pair = std::pair<int, int>({d_angmom, b_angmom});

    // set up Fock submatrices

    auto submat_pr = focks.matrix("PR_" + suffix)->sub_matrix({0, 0});

    auto submat_rp = focks.matrix("RP_" + suffix)->sub_matrix({0, 0});

    auto submat_ps = focks.matrix("PS_" + suffix)->sub_matrix({0, 0});

    auto submat_sp = focks.matrix("SP_" + suffix)->sub_matrix({0, 0});

    auto submat_qr = focks.matrix("QR_" + suffix)->sub_matrix({0, 0});

    auto submat_rq = focks.matrix("RQ_" + suffix)->sub_matrix({0, 0});

    auto submat_qs = focks.matrix("QS_" + suffix)->sub_matrix({0, 0});

    auto submat_sq = focks.matrix("SQ_" + suffix)->sub_matrix({0, 0});

    // set up AO density submatrix

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_rp = density->sub_matrix(rp_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_sp = density->sub_matrix(sp_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_rq = density->sub_matrix(rq_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    auto denmat_sq = density->sub_matrix(sq_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // reference local indexes on bra side

    const auto loc_refp = a_loc_indices[bra_igto + 1];

    const auto loc_refq = b_loc_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // dimensions of local bra and ket orbital indexes

    const auto alocdim = a_loc_indices[0];

    const auto blocdim = b_loc_indices[0];

    const auto clocdim = c_loc_indices[0];

    const auto dlocdim = d_loc_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        const auto loc_p = i * alocdim + loc_refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            const auto loc_q = j * blocdim + loc_refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // reference local indexes on ket side

                        const auto loc_refr = c_loc_indices[m + 1];

                        const auto loc_refs = d_loc_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto loc_r = k * clocdim + loc_refr;

                        const auto s = l * ddim + refs;

                        const auto loc_s = l * dlocdim + loc_refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = factor * curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Exchange contribution (F_pr)

                        submat_pr->at({loc_p, loc_r}) += fval * denmat_qs->operator[]({q, s});
                        submat_rp->at({loc_r, loc_p}) += fval * denmat_sq->operator[]({s, q});

                        // Exchange contribution (F_ps)

                        submat_ps->at({loc_p, loc_s}) += fval * denmat_qr->operator[]({q, r});
                        submat_sp->at({loc_s, loc_p}) += fval * denmat_rq->operator[]({r, q});

                        // Exchange contribution (F_qr)

                        submat_qr->at({loc_q, loc_r}) += fval * denmat_ps->operator[]({p, s});
                        submat_rq->at({loc_r, loc_q}) += fval * denmat_sp->operator[]({s, p});

                        // Exchange contribution (F_qs)

                        submat_qs->at({loc_q, loc_s}) += fval * denmat_pr->operator[]({p, r});
                        submat_sq->at({loc_s, loc_q}) += fval * denmat_rp->operator[]({r, p});
                    }
                }
            }
        }
    }
}

auto
local_distribute_jk_geom_sym(std::vector<double>&             values,            
                             const int                        cart_ind,
                             const CMatrix*                   density,
                             const CMatrix*                   density_2,
                             const CSimdArray<double>&        buffer,
                             const size_t                     offset,
                             const std::vector<size_t>&       a_indices,
                             const std::vector<size_t>&       b_indices,
                             const std::vector<size_t>&       c_indices,
                             const std::vector<size_t>&       d_indices,
                             const int                        a_angmom,
                             const int                        b_angmom,
                             const int                        c_angmom,
                             const int                        d_angmom,
                             const size_t                     bra_igto,
                             const std::pair<size_t, size_t>& ket_range) -> void
{
    // set up angular pairs

    const auto pq_pair = std::pair<int, int>({a_angmom, b_angmom});

    const auto qp_pair = std::pair<int, int>({b_angmom, a_angmom});

    const auto rs_pair = std::pair<int, int>({c_angmom, d_angmom});

    const auto sr_pair = std::pair<int, int>({d_angmom, c_angmom});

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto rp_pair = std::pair<int, int>({c_angmom, a_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto sp_pair = std::pair<int, int>({d_angmom, a_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto rq_pair = std::pair<int, int>({c_angmom, b_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    const auto sq_pair = std::pair<int, int>({d_angmom, b_angmom});

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    auto denmat_2_pq = density_2->sub_matrix(pq_pair);

    auto denmat_2_qp = density_2->sub_matrix(qp_pair);

    auto denmat_2_rs = density_2->sub_matrix(rs_pair);

    auto denmat_2_sr = density_2->sub_matrix(sr_pair);

    auto denmat_2_pr = density_2->sub_matrix(pr_pair);

    auto denmat_2_rp = density_2->sub_matrix(rp_pair);

    auto denmat_2_ps = density_2->sub_matrix(ps_pair);

    auto denmat_2_sp = density_2->sub_matrix(sp_pair);

    auto denmat_2_qr = density_2->sub_matrix(qr_pair);

    auto denmat_2_rq = density_2->sub_matrix(rq_pair);

    auto denmat_2_qs = density_2->sub_matrix(qs_pair);

    auto denmat_2_sq = density_2->sub_matrix(sq_pair);

    // set up angular orders

    const auto angord_pq = density->is_angular_order(pq_pair);

    const auto angord_rs = density->is_angular_order(rs_pair);

    const auto angord_pr = density->is_angular_order(pr_pair);

    const auto angord_ps = density->is_angular_order(ps_pair);

    const auto angord_qr = density->is_angular_order(qr_pair);

    const auto angord_qs = density->is_angular_order(qs_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        if (angord_rs)
                        {
                            values[cart_ind] += denmat_2_pq->operator[]({p, q}) * 4.0 * fval * denmat_rs->operator[]({r, s});

                            values[cart_ind] += denmat_2_qp->operator[]({q, p}) * 4.0 * fval * denmat_rs->operator[]({r, s});
                        }
                        else
                        {
                            values[cart_ind] += denmat_2_pq->operator[]({p, q}) * 4.0 * fval * denmat_rs->operator[]({s, r});

                            values[cart_ind] += denmat_2_qp->operator[]({q, p}) * 4.0 * fval * denmat_rs->operator[]({s, r});
                        }

                        if (angord_pq)
                        {
                            values[cart_ind] += denmat_2_rs->operator[]({r, s}) * 4.0 * fval * denmat_pq->operator[]({p, q});

                            values[cart_ind] += denmat_2_sr->operator[]({s, r}) * 4.0 * fval * denmat_pq->operator[]({p, q});
                        }
                        else
                        {
                            values[cart_ind] += denmat_2_rs->operator[]({r, s}) * 4.0 * fval * denmat_pq->operator[]({q, p});

                            values[cart_ind] += denmat_2_sr->operator[]({s, r}) * 4.0 * fval * denmat_pq->operator[]({q, p});
                        }

                        // Exchange contribution (F_pr)

                        if (angord_qs)
                        {
                            values[cart_ind] -= denmat_2_pr->operator[]({p, r}) * fval * denmat_qs->operator[]({q, s});

                            values[cart_ind] -= denmat_2_rp->operator[]({r, p}) * fval * denmat_qs->operator[]({q, s});
                        }
                        else
                        {
                            values[cart_ind] -= denmat_2_pr->operator[]({p, r}) * fval * denmat_qs->operator[]({s, q});

                            values[cart_ind] -= denmat_2_rp->operator[]({r, p}) * fval * denmat_qs->operator[]({s, q});
                        }

                        // Exchange contribution (F_ps)

                        if (angord_qr)
                        {
                            values[cart_ind] -= denmat_2_ps->operator[]({p, s}) * fval * denmat_qr->operator[]({q, r});

                            values[cart_ind] -= denmat_2_sp->operator[]({s, p}) * fval * denmat_qr->operator[]({q, r});
                        }
                        else
                        {
                            values[cart_ind] -= denmat_2_ps->operator[]({p, s}) * fval * denmat_qr->operator[]({r, q});

                            values[cart_ind] -= denmat_2_sp->operator[]({s, p}) * fval * denmat_qr->operator[]({r, q});
                        }

                        // Exchange contribution (F_qr)

                        if (angord_ps)
                        {
                            values[cart_ind] -= denmat_2_qr->operator[]({q, r}) * fval * denmat_ps->operator[]({p, s});

                            values[cart_ind] -= denmat_2_rq->operator[]({r, q}) * fval * denmat_ps->operator[]({p, s});
                        }
                        else
                        {
                            values[cart_ind] -= denmat_2_qr->operator[]({q, r}) * fval * denmat_ps->operator[]({s, p});

                            values[cart_ind] -= denmat_2_rq->operator[]({r, q}) * fval * denmat_ps->operator[]({s, p});
                        }

                        // Exchange contribution (F_qs)

                        if (angord_pr)
                        {
                            values[cart_ind] -= denmat_2_qs->operator[]({q, s}) * fval * denmat_pr->operator[]({p, r});

                            values[cart_ind] -= denmat_2_sq->operator[]({s, q}) * fval * denmat_pr->operator[]({p, r});
                        }
                        else
                        {
                            values[cart_ind] -= denmat_2_qs->operator[]({q, s}) * fval * denmat_pr->operator[]({r, p});

                            values[cart_ind] -= denmat_2_sq->operator[]({s, q}) * fval * denmat_pr->operator[]({r, p});
                        }
                    }
                }
            }
        }
    }
}

auto
local_distribute_jkx_geom_sym(std::vector<double>&             values,            
                              const int                        cart_ind,
                              const CMatrix*                   density,
                              const CMatrix*                   density_2,
                              const CSimdArray<double>&        buffer,
                              const size_t                     offset,
                              const double                     factor,
                              const std::vector<size_t>&       a_indices,
                              const std::vector<size_t>&       b_indices,
                              const std::vector<size_t>&       c_indices,
                              const std::vector<size_t>&       d_indices,
                              const int                        a_angmom,
                              const int                        b_angmom,
                              const int                        c_angmom,
                              const int                        d_angmom,
                              const size_t                     bra_igto,
                              const std::pair<size_t, size_t>& ket_range) -> void
{
    // set up angular pairs

    const auto pq_pair = std::pair<int, int>({a_angmom, b_angmom});

    const auto qp_pair = std::pair<int, int>({b_angmom, a_angmom});

    const auto rs_pair = std::pair<int, int>({c_angmom, d_angmom});

    const auto sr_pair = std::pair<int, int>({d_angmom, c_angmom});

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto rp_pair = std::pair<int, int>({c_angmom, a_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto sp_pair = std::pair<int, int>({d_angmom, a_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto rq_pair = std::pair<int, int>({c_angmom, b_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    const auto sq_pair = std::pair<int, int>({d_angmom, b_angmom});

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    auto denmat_2_pq = density_2->sub_matrix(pq_pair);

    auto denmat_2_qp = density_2->sub_matrix(qp_pair);

    auto denmat_2_rs = density_2->sub_matrix(rs_pair);

    auto denmat_2_sr = density_2->sub_matrix(sr_pair);

    auto denmat_2_pr = density_2->sub_matrix(pr_pair);

    auto denmat_2_rp = density_2->sub_matrix(rp_pair);

    auto denmat_2_ps = density_2->sub_matrix(ps_pair);

    auto denmat_2_sp = density_2->sub_matrix(sp_pair);

    auto denmat_2_qr = density_2->sub_matrix(qr_pair);

    auto denmat_2_rq = density_2->sub_matrix(rq_pair);

    auto denmat_2_qs = density_2->sub_matrix(qs_pair);

    auto denmat_2_sq = density_2->sub_matrix(sq_pair);

    // set up angular orders

    const auto angord_pq = density->is_angular_order(pq_pair);

    const auto angord_rs = density->is_angular_order(rs_pair);

    const auto angord_pr = density->is_angular_order(pr_pair);

    const auto angord_ps = density->is_angular_order(ps_pair);

    const auto angord_qr = density->is_angular_order(qr_pair);

    const auto angord_qs = density->is_angular_order(qs_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        if (angord_rs)
                        {
                            values[cart_ind] += denmat_2_pq->operator[]({p, q}) * 4.0 * fval * denmat_rs->operator[]({r, s});

                            values[cart_ind] += denmat_2_qp->operator[]({q, p}) * 4.0 * fval * denmat_rs->operator[]({r, s});
                        }
                        else
                        {
                            values[cart_ind] += denmat_2_pq->operator[]({p, q}) * 4.0 * fval * denmat_rs->operator[]({s, r});

                            values[cart_ind] += denmat_2_qp->operator[]({q, p}) * 4.0 * fval * denmat_rs->operator[]({s, r});
                        }

                        if (angord_pq)
                        {
                            values[cart_ind] += denmat_2_rs->operator[]({r, s}) * 4.0 * fval * denmat_pq->operator[]({p, q});

                            values[cart_ind] += denmat_2_sr->operator[]({s, r}) * 4.0 * fval * denmat_pq->operator[]({p, q});
                        }
                        else
                        {
                            values[cart_ind] += denmat_2_rs->operator[]({r, s}) * 4.0 * fval * denmat_pq->operator[]({q, p});

                            values[cart_ind] += denmat_2_sr->operator[]({s, r}) * 4.0 * fval * denmat_pq->operator[]({q, p});
                        }

                        // Exchange contribution

                        fval *= factor;

                        // Exchange contribution (F_pr)

                        if (angord_qs)
                        {
                            values[cart_ind] -= denmat_2_pr->operator[]({p, r}) * fval * denmat_qs->operator[]({q, s});

                            values[cart_ind] -= denmat_2_rp->operator[]({r, p}) * fval * denmat_qs->operator[]({q, s});
                        }
                        else
                        {
                            values[cart_ind] -= denmat_2_pr->operator[]({p, r}) * fval * denmat_qs->operator[]({s, q});

                            values[cart_ind] -= denmat_2_rp->operator[]({r, p}) * fval * denmat_qs->operator[]({s, q});
                        }

                        // Exchange contribution (F_ps)

                        if (angord_qr)
                        {
                            values[cart_ind] -= denmat_2_ps->operator[]({p, s}) * fval * denmat_qr->operator[]({q, r});

                            values[cart_ind] -= denmat_2_sp->operator[]({s, p}) * fval * denmat_qr->operator[]({q, r});
                        }
                        else
                        {
                            values[cart_ind] -= denmat_2_ps->operator[]({p, s}) * fval * denmat_qr->operator[]({r, q});

                            values[cart_ind] -= denmat_2_sp->operator[]({s, p}) * fval * denmat_qr->operator[]({r, q});
                        }

                        // Exchange contribution (F_qr)

                        if (angord_ps)
                        {
                            values[cart_ind] -= denmat_2_qr->operator[]({q, r}) * fval * denmat_ps->operator[]({p, s});

                            values[cart_ind] -= denmat_2_rq->operator[]({r, q}) * fval * denmat_ps->operator[]({p, s});
                        }
                        else
                        {
                            values[cart_ind] -= denmat_2_qr->operator[]({q, r}) * fval * denmat_ps->operator[]({s, p});

                            values[cart_ind] -= denmat_2_rq->operator[]({r, q}) * fval * denmat_ps->operator[]({s, p});
                        }

                        // Exchange contribution (F_qs)

                        if (angord_pr)
                        {
                            values[cart_ind] -= denmat_2_qs->operator[]({q, s}) * fval * denmat_pr->operator[]({p, r});

                            values[cart_ind] -= denmat_2_sq->operator[]({s, q}) * fval * denmat_pr->operator[]({p, r});
                        }
                        else
                        {
                            values[cart_ind] -= denmat_2_qs->operator[]({q, s}) * fval * denmat_pr->operator[]({r, p});

                            values[cart_ind] -= denmat_2_sq->operator[]({s, q}) * fval * denmat_pr->operator[]({r, p});
                        }
                    }
                }
            }
        }
    }
}

auto
local_distribute_j_geom_sym(std::vector<double>&             values,            
                            const int                        cart_ind,
                            const CMatrix*                   density,
                            const CMatrix*                   density_2,
                            const CSimdArray<double>&        buffer,
                            const size_t                     offset,
                            const std::vector<size_t>&       a_indices,
                            const std::vector<size_t>&       b_indices,
                            const std::vector<size_t>&       c_indices,
                            const std::vector<size_t>&       d_indices,
                            const int                        a_angmom,
                            const int                        b_angmom,
                            const int                        c_angmom,
                            const int                        d_angmom,
                            const size_t                     bra_igto,
                            const std::pair<size_t, size_t>& ket_range) -> void
{
    // set up angular pairs

    const auto pq_pair = std::pair<int, int>({a_angmom, b_angmom});

    const auto qp_pair = std::pair<int, int>({b_angmom, a_angmom});

    const auto rs_pair = std::pair<int, int>({c_angmom, d_angmom});

    const auto sr_pair = std::pair<int, int>({d_angmom, c_angmom});

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    auto denmat_2_pq = density_2->sub_matrix(pq_pair);

    auto denmat_2_qp = density_2->sub_matrix(qp_pair);

    auto denmat_2_rs = density_2->sub_matrix(rs_pair);

    auto denmat_2_sr = density_2->sub_matrix(sr_pair);

    // set up angular orders

    const auto angord_pq = density->is_angular_order(pq_pair);

    const auto angord_rs = density->is_angular_order(rs_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        if (angord_rs)
                        {
                            values[cart_ind] += denmat_2_pq->operator[]({p, q}) * 2.0 * fval * denmat_rs->operator[]({r, s});

                            values[cart_ind] += denmat_2_qp->operator[]({q, p}) * 2.0 * fval * denmat_rs->operator[]({r, s});
                        }
                        else
                        {
                            values[cart_ind] += denmat_2_pq->operator[]({p, q}) * 2.0 * fval * denmat_rs->operator[]({s, r});

                            values[cart_ind] += denmat_2_qp->operator[]({q, p}) * 2.0 * fval * denmat_rs->operator[]({s, r});
                        }

                        if (angord_pq)
                        {
                            values[cart_ind] += denmat_2_rs->operator[]({r, s}) * 2.0 * fval * denmat_pq->operator[]({p, q});

                            values[cart_ind] += denmat_2_sr->operator[]({s, r}) * 2.0 * fval * denmat_pq->operator[]({p, q});
                        }
                        else
                        {
                            values[cart_ind] += denmat_2_rs->operator[]({r, s}) * 2.0 * fval * denmat_pq->operator[]({q, p});

                            values[cart_ind] += denmat_2_sr->operator[]({s, r}) * 2.0 * fval * denmat_pq->operator[]({q, p});
                        }
                    }
                }
            }
        }
    }
}

auto
local_distribute_k_geom_sym(std::vector<double>&             values,            
                            const int                        cart_ind,
                            const CMatrix*                   density,
                            const CMatrix*                   density_2,
                            const CSimdArray<double>&        buffer,
                            const size_t                     offset,
                            const std::vector<size_t>&       a_indices,
                            const std::vector<size_t>&       b_indices,
                            const std::vector<size_t>&       c_indices,
                            const std::vector<size_t>&       d_indices,
                            const int                        a_angmom,
                            const int                        b_angmom,
                            const int                        c_angmom,
                            const int                        d_angmom,
                            const size_t                     bra_igto,
                            const std::pair<size_t, size_t>& ket_range) -> void
{
    // set up angular pairs

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto rp_pair = std::pair<int, int>({c_angmom, a_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto sp_pair = std::pair<int, int>({d_angmom, a_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto rq_pair = std::pair<int, int>({c_angmom, b_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    const auto sq_pair = std::pair<int, int>({d_angmom, b_angmom});

    // set up AO density submatrix

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    auto denmat_2_pr = density_2->sub_matrix(pr_pair);

    auto denmat_2_rp = density_2->sub_matrix(rp_pair);

    auto denmat_2_ps = density_2->sub_matrix(ps_pair);

    auto denmat_2_sp = density_2->sub_matrix(sp_pair);

    auto denmat_2_qr = density_2->sub_matrix(qr_pair);

    auto denmat_2_rq = density_2->sub_matrix(rq_pair);

    auto denmat_2_qs = density_2->sub_matrix(qs_pair);

    auto denmat_2_sq = density_2->sub_matrix(sq_pair);

    // set up angular orders

    const auto angord_pr = density->is_angular_order(pr_pair);

    const auto angord_ps = density->is_angular_order(ps_pair);

    const auto angord_qr = density->is_angular_order(qr_pair);

    const auto angord_qs = density->is_angular_order(qs_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Exchange contribution (F_pr)

                        if (angord_qs)
                        {
                            values[cart_ind] += denmat_2_pr->operator[]({p, r}) * fval * denmat_qs->operator[]({q, s});

                            values[cart_ind] += denmat_2_rp->operator[]({r, p}) * fval * denmat_qs->operator[]({q, s});
                        }
                        else
                        {
                            values[cart_ind] += denmat_2_pr->operator[]({p, r}) * fval * denmat_qs->operator[]({s, q});

                            values[cart_ind] += denmat_2_rp->operator[]({r, p}) * fval * denmat_qs->operator[]({s, q});
                        }

                        // Exchange contribution (F_ps)

                        if (angord_qr)
                        {
                            values[cart_ind] += denmat_2_ps->operator[]({p, s}) * fval * denmat_qr->operator[]({q, r});

                            values[cart_ind] += denmat_2_sp->operator[]({s, p}) * fval * denmat_qr->operator[]({q, r});
                        }
                        else
                        {
                            values[cart_ind] += denmat_2_ps->operator[]({p, s}) * fval * denmat_qr->operator[]({r, q});

                            values[cart_ind] += denmat_2_sp->operator[]({s, p}) * fval * denmat_qr->operator[]({r, q});
                        }

                        // Exchange contribution (F_qr)

                        if (angord_ps)
                        {
                            values[cart_ind] += denmat_2_qr->operator[]({q, r}) * fval * denmat_ps->operator[]({p, s});

                            values[cart_ind] += denmat_2_rq->operator[]({r, q}) * fval * denmat_ps->operator[]({p, s});
                        }
                        else
                        {
                            values[cart_ind] += denmat_2_qr->operator[]({q, r}) * fval * denmat_ps->operator[]({s, p});

                            values[cart_ind] += denmat_2_rq->operator[]({r, q}) * fval * denmat_ps->operator[]({s, p});
                        }

                        // Exchange contribution (F_qs)

                        if (angord_pr)
                        {
                            values[cart_ind] += denmat_2_qs->operator[]({q, s}) * fval * denmat_pr->operator[]({p, r});

                            values[cart_ind] += denmat_2_sq->operator[]({s, q}) * fval * denmat_pr->operator[]({p, r});
                        }
                        else
                        {
                            values[cart_ind] += denmat_2_qs->operator[]({q, s}) * fval * denmat_pr->operator[]({r, p});

                            values[cart_ind] += denmat_2_sq->operator[]({s, q}) * fval * denmat_pr->operator[]({r, p});
                        }
                    }
                }
            }
        }
    }
}

auto
local_distribute_kx_geom_sym(std::vector<double>&             values,            
                             const int                        cart_ind,
                             const CMatrix*                   density,
                             const CMatrix*                   density_2,
                             const CSimdArray<double>&        buffer,
                             const size_t                     offset,
                             const double                     factor,
                             const std::vector<size_t>&       a_indices,
                             const std::vector<size_t>&       b_indices,
                             const std::vector<size_t>&       c_indices,
                             const std::vector<size_t>&       d_indices,
                             const int                        a_angmom,
                             const int                        b_angmom,
                             const int                        c_angmom,
                             const int                        d_angmom,
                             const size_t                     bra_igto,
                             const std::pair<size_t, size_t>& ket_range) -> void
{
    // set up angular pairs

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto rp_pair = std::pair<int, int>({c_angmom, a_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto sp_pair = std::pair<int, int>({d_angmom, a_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto rq_pair = std::pair<int, int>({c_angmom, b_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    const auto sq_pair = std::pair<int, int>({d_angmom, b_angmom});

    // set up AO density submatrix

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    auto denmat_2_pr = density_2->sub_matrix(pr_pair);

    auto denmat_2_rp = density_2->sub_matrix(rp_pair);

    auto denmat_2_ps = density_2->sub_matrix(ps_pair);

    auto denmat_2_sp = density_2->sub_matrix(sp_pair);

    auto denmat_2_qr = density_2->sub_matrix(qr_pair);

    auto denmat_2_rq = density_2->sub_matrix(rq_pair);

    auto denmat_2_qs = density_2->sub_matrix(qs_pair);

    auto denmat_2_sq = density_2->sub_matrix(sq_pair);

    // set up angular orders

    const auto angord_pr = density->is_angular_order(pr_pair);

    const auto angord_ps = density->is_angular_order(ps_pair);

    const auto angord_qr = density->is_angular_order(qr_pair);

    const auto angord_qs = density->is_angular_order(qs_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = factor * curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Exchange contribution (F_pr)

                        if (angord_qs)
                        {
                            values[cart_ind] += denmat_2_pr->operator[]({p, r}) * fval * denmat_qs->operator[]({q, s});

                            values[cart_ind] += denmat_2_rp->operator[]({r, p}) * fval * denmat_qs->operator[]({q, s});
                        }
                        else
                        {
                            values[cart_ind] += denmat_2_pr->operator[]({p, r}) * fval * denmat_qs->operator[]({s, q});

                            values[cart_ind] += denmat_2_rp->operator[]({r, p}) * fval * denmat_qs->operator[]({s, q});
                        }

                        // Exchange contribution (F_ps)

                        if (angord_qr)
                        {
                            values[cart_ind] += denmat_2_ps->operator[]({p, s}) * fval * denmat_qr->operator[]({q, r});

                            values[cart_ind] += denmat_2_sp->operator[]({s, p}) * fval * denmat_qr->operator[]({q, r});
                        }
                        else
                        {
                            values[cart_ind] += denmat_2_ps->operator[]({p, s}) * fval * denmat_qr->operator[]({r, q});

                            values[cart_ind] += denmat_2_sp->operator[]({s, p}) * fval * denmat_qr->operator[]({r, q});
                        }

                        // Exchange contribution (F_qr)

                        if (angord_ps)
                        {
                            values[cart_ind] += denmat_2_qr->operator[]({q, r}) * fval * denmat_ps->operator[]({p, s});

                            values[cart_ind] += denmat_2_rq->operator[]({r, q}) * fval * denmat_ps->operator[]({p, s});
                        }
                        else
                        {
                            values[cart_ind] += denmat_2_qr->operator[]({q, r}) * fval * denmat_ps->operator[]({s, p});

                            values[cart_ind] += denmat_2_rq->operator[]({r, q}) * fval * denmat_ps->operator[]({s, p});
                        }

                        // Exchange contribution (F_qs)

                        if (angord_pr)
                        {
                            values[cart_ind] += denmat_2_qs->operator[]({q, s}) * fval * denmat_pr->operator[]({p, r});

                            values[cart_ind] += denmat_2_sq->operator[]({s, q}) * fval * denmat_pr->operator[]({p, r});
                        }
                        else
                        {
                            values[cart_ind] += denmat_2_qs->operator[]({q, s}) * fval * denmat_pr->operator[]({r, p});

                            values[cart_ind] += denmat_2_sq->operator[]({s, q}) * fval * denmat_pr->operator[]({r, p});
                        }
                    }
                }
            }
        }
    }
}

auto
local_distribute_jk_geom_gen(std::vector<double>&             values,            
                             const int                        cart_ind,
                             const CMatrix*                   density,
                             const CMatrix*                   density_2,
                             const CSimdArray<double>&        buffer,
                             const size_t                     offset,
                             const std::vector<size_t>&       a_indices,
                             const std::vector<size_t>&       b_indices,
                             const std::vector<size_t>&       c_indices,
                             const std::vector<size_t>&       d_indices,
                             const int                        a_angmom,
                             const int                        b_angmom,
                             const int                        c_angmom,
                             const int                        d_angmom,
                             const size_t                     bra_igto,
                             const std::pair<size_t, size_t>& ket_range) -> void
{
    // set up angular pairs

    const auto pq_pair = std::pair<int, int>({a_angmom, b_angmom});

    const auto qp_pair = std::pair<int, int>({b_angmom, a_angmom});

    const auto rs_pair = std::pair<int, int>({c_angmom, d_angmom});

    const auto sr_pair = std::pair<int, int>({d_angmom, c_angmom});

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto rp_pair = std::pair<int, int>({c_angmom, a_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto sp_pair = std::pair<int, int>({d_angmom, a_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto rq_pair = std::pair<int, int>({c_angmom, b_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    const auto sq_pair = std::pair<int, int>({d_angmom, b_angmom});

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_qp = density->sub_matrix(qp_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    auto denmat_sr = density->sub_matrix(sr_pair);

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_rp = density->sub_matrix(rp_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_sp = density->sub_matrix(sp_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_rq = density->sub_matrix(rq_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    auto denmat_sq = density->sub_matrix(sq_pair);

    auto denmat_2_pq = density_2->sub_matrix(pq_pair);

    auto denmat_2_qp = density_2->sub_matrix(qp_pair);

    auto denmat_2_rs = density_2->sub_matrix(rs_pair);

    auto denmat_2_sr = density_2->sub_matrix(sr_pair);

    auto denmat_2_pr = density_2->sub_matrix(pr_pair);

    auto denmat_2_rp = density_2->sub_matrix(rp_pair);

    auto denmat_2_ps = density_2->sub_matrix(ps_pair);

    auto denmat_2_sp = density_2->sub_matrix(sp_pair);

    auto denmat_2_qr = density_2->sub_matrix(qr_pair);

    auto denmat_2_rq = density_2->sub_matrix(rq_pair);

    auto denmat_2_qs = density_2->sub_matrix(qs_pair);

    auto denmat_2_sq = density_2->sub_matrix(sq_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        const auto f2rs = 2.0 * fval * (denmat_rs->operator[]({r, s}) + denmat_sr->operator[]({s, r}));

                        values[cart_ind] += denmat_2_pq->operator[]({p, q}) * f2rs;
                        values[cart_ind] += denmat_2_qp->operator[]({q, p}) * f2rs;

                        // Coulomb contribution (F_rs)

                        const auto f2pq = 2.0 * fval * (denmat_pq->operator[]({p, q}) + denmat_qp->operator[]({q, p}));

                        values[cart_ind] += denmat_2_rs->operator[]({r, s}) * f2pq;
                        values[cart_ind] += denmat_2_sr->operator[]({s, r}) * f2pq;

                        // Exchange contribution (F_pr)

                        values[cart_ind] -= denmat_2_pr->operator[]({p, r}) * fval * denmat_qs->operator[]({q, s});
                        values[cart_ind] -= denmat_2_rp->operator[]({r, p}) * fval * denmat_sq->operator[]({s, q});

                        // Exchange contribution (F_ps)

                        values[cart_ind] -= denmat_2_ps->operator[]({p, s}) * fval * denmat_qr->operator[]({q, r});
                        values[cart_ind] -= denmat_2_sp->operator[]({s, p}) * fval * denmat_rq->operator[]({r, q});

                        // Exchange contribution (F_qr)

                        values[cart_ind] -= denmat_2_qr->operator[]({q, r}) * fval * denmat_ps->operator[]({p, s});
                        values[cart_ind] -= denmat_2_rq->operator[]({r, q}) * fval * denmat_sp->operator[]({s, p});

                        // Exchange contribution (F_qs)

                        values[cart_ind] -= denmat_2_qs->operator[]({q, s}) * fval * denmat_pr->operator[]({p, r});
                        values[cart_ind] -= denmat_2_sq->operator[]({s, q}) * fval * denmat_rp->operator[]({r, p});
                    }
                }
            }
        }
    }
}

auto
local_distribute_jkx_geom_gen(std::vector<double>&             values,            
                              const int                        cart_ind,
                              const CMatrix*                   density,
                              const CMatrix*                   density_2,
                              const CSimdArray<double>&        buffer,
                              const size_t                     offset,
                              const double                     factor,
                              const std::vector<size_t>&       a_indices,
                              const std::vector<size_t>&       b_indices,
                              const std::vector<size_t>&       c_indices,
                              const std::vector<size_t>&       d_indices,
                              const int                        a_angmom,
                              const int                        b_angmom,
                              const int                        c_angmom,
                              const int                        d_angmom,
                              const size_t                     bra_igto,
                              const std::pair<size_t, size_t>& ket_range) -> void
{
    // set up angular pairs

    const auto pq_pair = std::pair<int, int>({a_angmom, b_angmom});

    const auto qp_pair = std::pair<int, int>({b_angmom, a_angmom});

    const auto rs_pair = std::pair<int, int>({c_angmom, d_angmom});

    const auto sr_pair = std::pair<int, int>({d_angmom, c_angmom});

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto rp_pair = std::pair<int, int>({c_angmom, a_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto sp_pair = std::pair<int, int>({d_angmom, a_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto rq_pair = std::pair<int, int>({c_angmom, b_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    const auto sq_pair = std::pair<int, int>({d_angmom, b_angmom});

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_qp = density->sub_matrix(qp_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    auto denmat_sr = density->sub_matrix(sr_pair);

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_rp = density->sub_matrix(rp_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_sp = density->sub_matrix(sp_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_rq = density->sub_matrix(rq_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    auto denmat_sq = density->sub_matrix(sq_pair);

    auto denmat_2_pq = density_2->sub_matrix(pq_pair);

    auto denmat_2_qp = density_2->sub_matrix(qp_pair);

    auto denmat_2_rs = density_2->sub_matrix(rs_pair);

    auto denmat_2_sr = density_2->sub_matrix(sr_pair);

    auto denmat_2_pr = density_2->sub_matrix(pr_pair);

    auto denmat_2_rp = density_2->sub_matrix(rp_pair);

    auto denmat_2_ps = density_2->sub_matrix(ps_pair);

    auto denmat_2_sp = density_2->sub_matrix(sp_pair);

    auto denmat_2_qr = density_2->sub_matrix(qr_pair);

    auto denmat_2_rq = density_2->sub_matrix(rq_pair);

    auto denmat_2_qs = density_2->sub_matrix(qs_pair);

    auto denmat_2_sq = density_2->sub_matrix(sq_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        const auto f2rs = 2.0 * fval * (denmat_rs->operator[]({r, s}) + denmat_sr->operator[]({s, r}));

                        values[cart_ind] += denmat_2_pq->operator[]({p, q}) * f2rs;
                        values[cart_ind] += denmat_2_qp->operator[]({q, p}) * f2rs;

                        // Coulomb contribution (F_rs)

                        const auto f2pq = 2.0 * fval * (denmat_pq->operator[]({p, q}) + denmat_qp->operator[]({q, p}));

                        values[cart_ind] += denmat_2_rs->operator[]({r, s}) * f2pq;
                        values[cart_ind] += denmat_2_sr->operator[]({s, r}) * f2pq;

                        // Exchange contribution

                        fval *= factor;

                        // Exchange contribution (F_pr)

                        values[cart_ind] -= denmat_2_pr->operator[]({p, r}) * fval * denmat_qs->operator[]({q, s});
                        values[cart_ind] -= denmat_2_rp->operator[]({r, p}) * fval * denmat_sq->operator[]({s, q});

                        // Exchange contribution (F_ps)

                        values[cart_ind] -= denmat_2_ps->operator[]({p, s}) * fval * denmat_qr->operator[]({q, r});
                        values[cart_ind] -= denmat_2_sp->operator[]({s, p}) * fval * denmat_rq->operator[]({r, q});

                        // Exchange contribution (F_qr)

                        values[cart_ind] -= denmat_2_qr->operator[]({q, r}) * fval * denmat_ps->operator[]({p, s});
                        values[cart_ind] -= denmat_2_rq->operator[]({r, q}) * fval * denmat_sp->operator[]({s, p});

                        // Exchange contribution (F_qs)

                        values[cart_ind] -= denmat_2_qs->operator[]({q, s}) * fval * denmat_pr->operator[]({p, r});
                        values[cart_ind] -= denmat_2_sq->operator[]({s, q}) * fval * denmat_rp->operator[]({r, p});
                    }
                }
            }
        }
    }
}

auto
local_distribute_j_geom_gen(std::vector<double>&             values,            
                            const int                        cart_ind,
                            const CMatrix*                   density,
                            const CMatrix*                   density_2,
                            const CSimdArray<double>&        buffer,
                            const size_t                     offset,
                            const std::vector<size_t>&       a_indices,
                            const std::vector<size_t>&       b_indices,
                            const std::vector<size_t>&       c_indices,
                            const std::vector<size_t>&       d_indices,
                            const int                        a_angmom,
                            const int                        b_angmom,
                            const int                        c_angmom,
                            const int                        d_angmom,
                            const size_t                     bra_igto,
                            const std::pair<size_t, size_t>& ket_range) -> void
{
    // set up angular pairs

    const auto pq_pair = std::pair<int, int>({a_angmom, b_angmom});

    const auto qp_pair = std::pair<int, int>({b_angmom, a_angmom});

    const auto rs_pair = std::pair<int, int>({c_angmom, d_angmom});

    const auto sr_pair = std::pair<int, int>({d_angmom, c_angmom});

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_qp = density->sub_matrix(qp_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    auto denmat_sr = density->sub_matrix(sr_pair);

    auto denmat_2_pq = density_2->sub_matrix(pq_pair);

    auto denmat_2_qp = density_2->sub_matrix(qp_pair);

    auto denmat_2_rs = density_2->sub_matrix(rs_pair);

    auto denmat_2_sr = density_2->sub_matrix(sr_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        const auto f_rs = fval * (denmat_rs->operator[]({r, s}) + denmat_sr->operator[]({s, r}));

                        values[cart_ind] += denmat_2_pq->operator[]({p, q}) * f_rs;
                        values[cart_ind] += denmat_2_qp->operator[]({q, p}) * f_rs;

                        // Coulomb contribution (F_rs)

                        const auto f_pq = fval * (denmat_pq->operator[]({p, q}) + denmat_qp->operator[]({q, p}));

                        values[cart_ind] += denmat_2_rs->operator[]({r, s}) * f_pq;
                        values[cart_ind] += denmat_2_sr->operator[]({s, r}) * f_pq;
                    }
                }
            }
        }
    }
}

auto
local_distribute_k_geom_gen(std::vector<double>&             values,            
                            const int                        cart_ind,
                            const CMatrix*                   density,
                            const CMatrix*                   density_2,
                            const CSimdArray<double>&        buffer,
                            const size_t                     offset,
                            const std::vector<size_t>&       a_indices,
                            const std::vector<size_t>&       b_indices,
                            const std::vector<size_t>&       c_indices,
                            const std::vector<size_t>&       d_indices,
                            const int                        a_angmom,
                            const int                        b_angmom,
                            const int                        c_angmom,
                            const int                        d_angmom,
                            const size_t                     bra_igto,
                            const std::pair<size_t, size_t>& ket_range) -> void
{
    // set up angular pairs

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto rp_pair = std::pair<int, int>({c_angmom, a_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto sp_pair = std::pair<int, int>({d_angmom, a_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto rq_pair = std::pair<int, int>({c_angmom, b_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    const auto sq_pair = std::pair<int, int>({d_angmom, b_angmom});

    // set up AO density submatrix

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_rp = density->sub_matrix(rp_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_sp = density->sub_matrix(sp_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_rq = density->sub_matrix(rq_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    auto denmat_sq = density->sub_matrix(sq_pair);

    auto denmat_2_pr = density_2->sub_matrix(pr_pair);

    auto denmat_2_rp = density_2->sub_matrix(rp_pair);

    auto denmat_2_ps = density_2->sub_matrix(ps_pair);

    auto denmat_2_sp = density_2->sub_matrix(sp_pair);

    auto denmat_2_qr = density_2->sub_matrix(qr_pair);

    auto denmat_2_rq = density_2->sub_matrix(rq_pair);

    auto denmat_2_qs = density_2->sub_matrix(qs_pair);

    auto denmat_2_sq = density_2->sub_matrix(sq_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Exchange contribution (F_pr)

                        values[cart_ind] += denmat_2_pr->operator[]({p, r}) * fval * denmat_qs->operator[]({q, s});
                        values[cart_ind] += denmat_2_rp->operator[]({r, p}) * fval * denmat_sq->operator[]({s, q});

                        // Exchange contribution (F_ps)

                        values[cart_ind] += denmat_2_ps->operator[]({p, s}) * fval * denmat_qr->operator[]({q, r});
                        values[cart_ind] += denmat_2_sp->operator[]({s, p}) * fval * denmat_rq->operator[]({r, q});

                        // Exchange contribution (F_qr)

                        values[cart_ind] += denmat_2_qr->operator[]({q, r}) * fval * denmat_ps->operator[]({p, s});
                        values[cart_ind] += denmat_2_rq->operator[]({r, q}) * fval * denmat_sp->operator[]({s, p});

                        // Exchange contribution (F_qs)

                        values[cart_ind] += denmat_2_qs->operator[]({q, s}) * fval * denmat_pr->operator[]({p, r});
                        values[cart_ind] += denmat_2_sq->operator[]({s, q}) * fval * denmat_rp->operator[]({r, p});
                    }
                }
            }
        }
    }
}

auto
local_distribute_kx_geom_gen(std::vector<double>&             values,            
                             const int                        cart_ind,
                             const CMatrix*                   density,
                             const CMatrix*                   density_2,
                             const CSimdArray<double>&        buffer,
                             const size_t                     offset,
                             const double                     factor,
                             const std::vector<size_t>&       a_indices,
                             const std::vector<size_t>&       b_indices,
                             const std::vector<size_t>&       c_indices,
                             const std::vector<size_t>&       d_indices,
                             const int                        a_angmom,
                             const int                        b_angmom,
                             const int                        c_angmom,
                             const int                        d_angmom,
                             const size_t                     bra_igto,
                             const std::pair<size_t, size_t>& ket_range) -> void
{
    // set up angular pairs

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto rp_pair = std::pair<int, int>({c_angmom, a_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto sp_pair = std::pair<int, int>({d_angmom, a_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto rq_pair = std::pair<int, int>({c_angmom, b_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    const auto sq_pair = std::pair<int, int>({d_angmom, b_angmom});

    // set up AO density submatrix

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_rp = density->sub_matrix(rp_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_sp = density->sub_matrix(sp_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_rq = density->sub_matrix(rq_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    auto denmat_sq = density->sub_matrix(sq_pair);

    auto denmat_2_pr = density_2->sub_matrix(pr_pair);

    auto denmat_2_rp = density_2->sub_matrix(rp_pair);

    auto denmat_2_ps = density_2->sub_matrix(ps_pair);

    auto denmat_2_sp = density_2->sub_matrix(sp_pair);

    auto denmat_2_qr = density_2->sub_matrix(qr_pair);

    auto denmat_2_rq = density_2->sub_matrix(rq_pair);

    auto denmat_2_qs = density_2->sub_matrix(qs_pair);

    auto denmat_2_sq = density_2->sub_matrix(sq_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{b_angmom});

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = factor * curr_buffer[m - ket_range.first];

                        if (r == s) fval *= 0.5;

                        // Exchange contribution (F_pr)

                        values[cart_ind] += denmat_2_pr->operator[]({p, r}) * fval * denmat_qs->operator[]({q, s});
                        values[cart_ind] += denmat_2_rp->operator[]({r, p}) * fval * denmat_sq->operator[]({s, q});

                        // Exchange contribution (F_ps)

                        values[cart_ind] += denmat_2_ps->operator[]({p, s}) * fval * denmat_qr->operator[]({q, r});
                        values[cart_ind] += denmat_2_sp->operator[]({s, p}) * fval * denmat_rq->operator[]({r, q});

                        // Exchange contribution (F_qr)

                        values[cart_ind] += denmat_2_qr->operator[]({q, r}) * fval * denmat_ps->operator[]({p, s});
                        values[cart_ind] += denmat_2_rq->operator[]({r, q}) * fval * denmat_sp->operator[]({s, p});

                        // Exchange contribution (F_qs)

                        values[cart_ind] += denmat_2_qs->operator[]({q, s}) * fval * denmat_pr->operator[]({p, r});
                        values[cart_ind] += denmat_2_sq->operator[]({s, q}) * fval * denmat_rp->operator[]({r, p});
                    }
                }
            }
        }
    }
}

}  // namespace t4cfunc
