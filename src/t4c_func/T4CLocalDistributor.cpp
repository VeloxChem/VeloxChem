#include "T4CLocalDistributor.hpp"

#include <iostream>
#include <iomanip>

#include "TensorComponents.hpp"

namespace t4cfunc {  // t2cfunc namespace

auto
local_distribute_rest_jk(CMatrices&                       focks,
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
                         const std::pair<size_t, size_t>& ket_range,
                         const bool                       diagonal) -> void
{
    // set up angular pairs

    const auto pq_pair = std::pair<int, int>({a_angmom, b_angmom});

    const auto rs_pair = std::pair<int, int>({c_angmom, d_angmom});

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    // set up Fock submatrices

    auto submat_pq = focks.matrix("PQ")->sub_matrix({0, 0});

    auto submat_rs = focks.matrix("RS")->sub_matrix({0, 0});

    auto submat_pr = focks.matrix("PR")->sub_matrix({0, 0});

    auto submat_ps = focks.matrix("PS")->sub_matrix({0, 0});

    auto submat_qr = focks.matrix("QR")->sub_matrix({0, 0});

    auto submat_qs = focks.matrix("QS")->sub_matrix({0, 0});

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
            // impose angular symmetry on bra side

            if (refp == refq)
            {
                if (j < i) continue;
            }

            const auto q = j * bdim + refq;

            const auto loc_q = j * blocdim + loc_refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // skip repeating integrals in diagonal block

                        if (diagonal)
                        {
                            if (m < bra_igto) continue;
                        }

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

                        // impose angular symmetry for itentical bra and ket sides

                        if ((refp == refr) && (refq == refs))
                        {
                            if (k < i) continue;

                            if (i == k)
                            {
                                if (l < j) continue;
                            }
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto loc_r = k * clocdim + loc_refr;

                        const auto s = l * ddim + refs;

                        const auto loc_s = l * dlocdim + loc_refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (p == q) fval *= 0.5;

                        if (r == s) fval *= 0.5;

                        if ((p == r) && (q == s)) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        if (angord_rs)
                        {
                            submat_pq->at({loc_p, loc_q}) += 4.0 * fval * denmat_rs->operator[]({r, s});
                        }
                        else
                        {
                            submat_pq->at({loc_p, loc_q}) += 4.0 * fval * denmat_rs->operator[]({s, r});
                        }

                        // Coulomb contribution (F_rs)

                        if (angord_pq)
                        {
                            submat_rs->at({loc_r, loc_s}) += 4.0 * fval * denmat_pq->operator[]({p, q});
                        }
                        else
                        {
                            submat_rs->at({loc_r, loc_s}) += 4.0 * fval * denmat_pq->operator[]({q, p});
                        }

                        // Exchange contribution (F_pr)

                        if (angord_qs)
                        {
                            submat_pr->at({loc_p, loc_r}) -= fval * denmat_qs->operator[]({q, s});
                        }
                        else
                        {
                            submat_pr->at({loc_p, loc_r}) -= fval * denmat_qs->operator[]({s, q});
                        }

                        // Exchange contribution (F_ps)

                        if (angord_qr)
                        {
                            submat_ps->at({loc_p, loc_s}) -= fval * denmat_qr->operator[]({q, r});
                        }
                        else
                        {
                            submat_ps->at({loc_p, loc_s}) -= fval * denmat_qr->operator[]({r, q});
                        }

                        // Exchange contribution (F_qr)

                        if (angord_ps)
                        {
                            submat_qr->at({loc_q, loc_r}) -= fval * denmat_ps->operator[]({p, s});
                        }
                        else
                        {
                            submat_qr->at({loc_q, loc_r}) -= fval * denmat_ps->operator[]({s, p});
                        }

                        // Exchange contribution (F_qs)

                        if (angord_pr)
                        {
                            submat_qs->at({loc_q, loc_s}) -= fval * denmat_pr->operator[]({p, r});
                        }
                        else
                        {
                            submat_qs->at({loc_q, loc_s}) -= fval * denmat_pr->operator[]({r, p});
                        }
                    }
                }
            }
        }
    }
}

auto
local_distribute_rest_jkx(CMatrices&                       focks,
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
                          const std::pair<size_t, size_t>& ket_range,
                          const bool                       diagonal) -> void
{
    // set up angular pairs

    const auto pq_pair = std::pair<int, int>({a_angmom, b_angmom});

    const auto rs_pair = std::pair<int, int>({c_angmom, d_angmom});

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    // set up Fock submatrices

    auto submat_pq = focks.matrix("PQ")->sub_matrix({0, 0});

    auto submat_rs = focks.matrix("RS")->sub_matrix({0, 0});

    auto submat_pr = focks.matrix("PR")->sub_matrix({0, 0});

    auto submat_ps = focks.matrix("PS")->sub_matrix({0, 0});

    auto submat_qr = focks.matrix("QR")->sub_matrix({0, 0});

    auto submat_qs = focks.matrix("QS")->sub_matrix({0, 0});

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
            // impose angular symmetry on bra side

            if (refp == refq)
            {
                if (j < i) continue;
            }

            const auto q = j * bdim + refq;

            const auto loc_q = j * blocdim + loc_refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // skip repeating integrals in diagonal block

                        if (diagonal)
                        {
                            if (m < bra_igto) continue;
                        }

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

                        // impose angular symmetry for itentical bra and ket sides

                        if ((refp == refr) && (refq == refs))
                        {
                            if (k < i) continue;

                            if (i == k)
                            {
                                if (l < j) continue;
                            }
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto loc_r = k * clocdim + loc_refr;

                        const auto s = l * ddim + refs;

                        const auto loc_s = l * dlocdim + loc_refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (p == q) fval *= 0.5;

                        if (r == s) fval *= 0.5;

                        if ((p == r) && (q == s)) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        if (angord_rs)
                        {
                            submat_pq->at({loc_p, loc_q}) += 4.0 * fval * denmat_rs->operator[]({r, s});
                        }
                        else
                        {
                            submat_pq->at({loc_p, loc_q}) += 4.0 * fval * denmat_rs->operator[]({s, r});
                        }

                        // Coulomb contribution (F_rs)

                        if (angord_pq)
                        {
                            submat_rs->at({loc_r, loc_s}) += 4.0 * fval * denmat_pq->operator[]({p, q});
                        }
                        else
                        {
                            submat_rs->at({loc_r, loc_s}) += 4.0 * fval * denmat_pq->operator[]({q, p});
                        }

                        // Exchange contribution (F_pr)

                        fval *= factor;

                        if (angord_qs)
                        {
                            submat_pr->at({loc_p, loc_r}) -= fval * denmat_qs->operator[]({q, s});
                        }
                        else
                        {
                            submat_pr->at({loc_p, loc_r}) -= fval * denmat_qs->operator[]({s, q});
                        }

                        // Exchange contribution (F_ps)

                        if (angord_qr)
                        {
                            submat_ps->at({loc_p, loc_s}) -= fval * denmat_qr->operator[]({q, r});
                        }
                        else
                        {
                            submat_ps->at({loc_p, loc_s}) -= fval * denmat_qr->operator[]({r, q});
                        }

                        // Exchange contribution (F_qr)

                        if (angord_ps)
                        {
                            submat_qr->at({loc_q, loc_r}) -= fval * denmat_ps->operator[]({p, s});
                        }
                        else
                        {
                            submat_qr->at({loc_q, loc_r}) -= fval * denmat_ps->operator[]({s, p});
                        }

                        // Exchange contribution (F_qs)

                        if (angord_pr)
                        {
                            submat_qs->at({loc_q, loc_s}) -= fval * denmat_pr->operator[]({p, r});
                        }
                        else
                        {
                            submat_qs->at({loc_q, loc_s}) -= fval * denmat_pr->operator[]({r, p});
                        }
                    }
                }
            }
        }
    }
}

auto
local_distribute_rest_j(CMatrices&                       focks,
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
                        const std::pair<size_t, size_t>& ket_range,
                        const bool                       diagonal) -> void
{
    // set up angular pairs

    const auto pq_pair = std::pair<int, int>({a_angmom, b_angmom});

    const auto rs_pair = std::pair<int, int>({c_angmom, d_angmom});

    // set up Fock submatrices

    auto submat_pq = focks.matrix("PQ")->sub_matrix({0, 0});

    auto submat_rs = focks.matrix("RS")->sub_matrix({0, 0});

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
            // impose angular symmetry on bra side

            if (refp == refq)
            {
                if (j < i) continue;
            }

            const auto q = j * bdim + refq;

            const auto loc_q = j * blocdim + loc_refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // skip repeating integrals in diagonal block

                        if (diagonal)
                        {
                            if (m < bra_igto) continue;
                        }

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

                        // impose angular symmetry for itentical bra and ket sides

                        if ((refp == refr) && (refq == refs))
                        {
                            if (k < i) continue;

                            if (i == k)
                            {
                                if (l < j) continue;
                            }
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto loc_r = k * clocdim + loc_refr;

                        const auto s = l * ddim + refs;

                        const auto loc_s = l * dlocdim + loc_refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

//                        if ((a_angmom == 3) && (b_angmom == 3) && (c_angmom == 3) && (d_angmom == 3))
//                        {
//                            std::cout << std::setprecision(15);
//                            
//                            std::cout  << p << " " << q << " " << r << " " << s << " " << fval << std::endl;
//                        }

                        if (p == q) fval *= 0.5;

                        if (r == s) fval *= 0.5;

                        if ((p == r) && (q == s)) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        if (angord_rs)
                        {
                            submat_pq->at({loc_p, loc_q}) += 2.0 * fval * denmat_rs->operator[]({r, s});
                        }
                        else
                        {
                            submat_pq->at({loc_p, loc_q}) += 2.0 * fval * denmat_rs->operator[]({s, r});
                        }

                        // Coulomb contribution (F_rs)

                        if (angord_pq)
                        {
                            submat_rs->at({loc_r, loc_s}) += 2.0 * fval * denmat_pq->operator[]({p, q});
                        }
                        else
                        {
                            submat_rs->at({loc_r, loc_s}) += 2.0 * fval * denmat_pq->operator[]({q, p});
                        }
                    }
                }
            }
        }
    }
}

auto
local_distribute_rest_k(CMatrices&                       focks,
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
                        const std::pair<size_t, size_t>& ket_range,
                        const bool                       diagonal) -> void
{
    // set up angular pairs

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    // set up Fock submatrices

    auto submat_pr = focks.matrix("PR")->sub_matrix({0, 0});

    auto submat_ps = focks.matrix("PS")->sub_matrix({0, 0});

    auto submat_qr = focks.matrix("QR")->sub_matrix({0, 0});

    auto submat_qs = focks.matrix("QS")->sub_matrix({0, 0});

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
            // impose angular symmetry on bra side

            if (refp == refq)
            {
                if (j < i) continue;
            }

            const auto q = j * bdim + refq;

            const auto loc_q = j * blocdim + loc_refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // skip repeating integrals in diagonal block

                        if (diagonal)
                        {
                            if (m < bra_igto) continue;
                        }

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

                        // impose angular symmetry for itentical bra and ket sides

                        if ((refp == refr) && (refq == refs))
                        {
                            if (k < i) continue;

                            if (i == k)
                            {
                                if (l < j) continue;
                            }
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto loc_r = k * clocdim + loc_refr;

                        const auto s = l * ddim + refs;

                        const auto loc_s = l * dlocdim + loc_refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range.first];

                        if (p == q) fval *= 0.5;

                        if (r == s) fval *= 0.5;

                        if ((p == r) && (q == s)) fval *= 0.5;

                        // Exchange contribution (F_pr)

                        if (angord_qs)
                        {
                            submat_pr->at({loc_p, loc_r}) += fval * denmat_qs->operator[]({q, s});
                        }
                        else
                        {
                            submat_pr->at({loc_p, loc_r}) += fval * denmat_qs->operator[]({s, q});
                        }

                        // Exchange contribution (F_ps)

                        if (angord_qr)
                        {
                            submat_ps->at({loc_p, loc_s}) += fval * denmat_qr->operator[]({q, r});
                        }
                        else
                        {
                            submat_ps->at({loc_p, loc_s}) += fval * denmat_qr->operator[]({r, q});
                        }

                        // Exchange contribution (F_qr)

                        if (angord_ps)
                        {
                            submat_qr->at({loc_q, loc_r}) += fval * denmat_ps->operator[]({p, s});
                        }
                        else
                        {
                            submat_qr->at({loc_q, loc_r}) += fval * denmat_ps->operator[]({s, p});
                        }

                        // Exchange contribution (F_qs)

                        if (angord_pr)
                        {
                            submat_qs->at({loc_q, loc_s}) += fval * denmat_pr->operator[]({p, r});
                        }
                        else
                        {
                            submat_qs->at({loc_q, loc_s}) += fval * denmat_pr->operator[]({r, p});
                        }
                    }
                }
            }
        }
    }
}

auto
local_distribute_rest_kx(CMatrices&                       focks,
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
                         const std::pair<size_t, size_t>& ket_range,
                         const bool                       diagonal) -> void
{
    // set up angular pairs

    const auto pr_pair = std::pair<int, int>({a_angmom, c_angmom});

    const auto ps_pair = std::pair<int, int>({a_angmom, d_angmom});

    const auto qr_pair = std::pair<int, int>({b_angmom, c_angmom});

    const auto qs_pair = std::pair<int, int>({b_angmom, d_angmom});

    // set up Fock submatrices

    auto submat_pr = focks.matrix("PR")->sub_matrix({0, 0});

    auto submat_ps = focks.matrix("PS")->sub_matrix({0, 0});

    auto submat_qr = focks.matrix("QR")->sub_matrix({0, 0});

    auto submat_qs = focks.matrix("QS")->sub_matrix({0, 0});

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
            // impose angular symmetry on bra side

            if (refp == refq)
            {
                if (j < i) continue;
            }

            const auto q = j * bdim + refq;

            const auto loc_q = j * blocdim + loc_refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l);

                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // skip repeating integrals in diagonal block

                        if (diagonal)
                        {
                            if (m < bra_igto) continue;
                        }

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

                        // impose angular symmetry for itentical bra and ket sides

                        if ((refp == refr) && (refq == refs))
                        {
                            if (k < i) continue;

                            if (i == k)
                            {
                                if (l < j) continue;
                            }
                        }

                        // compute r and s indexes

                        const auto r = k * cdim + refr;

                        const auto loc_r = k * clocdim + loc_refr;

                        const auto s = l * ddim + refs;

                        const auto loc_s = l * dlocdim + loc_refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = factor * curr_buffer[m - ket_range.first];

                        if (p == q) fval *= 0.5;

                        if (r == s) fval *= 0.5;

                        if ((p == r) && (q == s)) fval *= 0.5;

                        // Exchange contribution (F_pr)

                        if (angord_qs)
                        {
                            submat_pr->at({loc_p, loc_r}) += fval * denmat_qs->operator[]({q, s});
                        }
                        else
                        {
                            submat_pr->at({loc_p, loc_r}) += fval * denmat_qs->operator[]({s, q});
                        }

                        // Exchange contribution (F_ps)

                        if (angord_qr)
                        {
                            submat_ps->at({loc_p, loc_s}) += fval * denmat_qr->operator[]({q, r});
                        }
                        else
                        {
                            submat_ps->at({loc_p, loc_s}) += fval * denmat_qr->operator[]({r, q});
                        }

                        // Exchange contribution (F_qr)

                        if (angord_ps)
                        {
                            submat_qr->at({loc_q, loc_r}) += fval * denmat_ps->operator[]({p, s});
                        }
                        else
                        {
                            submat_qr->at({loc_q, loc_r}) += fval * denmat_ps->operator[]({s, p});
                        }

                        // Exchange contribution (F_qs)

                        if (angord_pr)
                        {
                            submat_qs->at({loc_q, loc_s}) += fval * denmat_pr->operator[]({p, r});
                        }
                        else
                        {
                            submat_qs->at({loc_q, loc_s}) += fval * denmat_pr->operator[]({r, p});
                        }
                    }
                }
            }
        }
    }
}

}  // namespace t4cfunc
