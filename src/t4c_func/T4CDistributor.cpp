#include "T4CDistributor.hpp"

#include <iostream>

#include "TensorComponents.hpp"

namespace t4cfunc {  // t2cfunc namespace

auto
distribute_rest_jk(CMatrix*                  fock,
                   const CMatrix*            density,
                   const CSimdArray<double>& buffer,
                   const int                 offset,
                   const std::vector<int>&   a_indices,
                   const std::vector<int>&   b_indices,
                   const std::vector<int>&   c_indices,
                   const std::vector<int>&   d_indices,
                   const int                 a_angmom,
                   const int                 b_angmom,
                   const int                 c_angmom,
                   const int                 d_angmom,
                   const int                 bra_igto,
                   const std::array<int, 2>& ket_range,
                   const bool                diagonal) -> void
{
    // set up angular pairs

    const auto pq_pair = std::array<int, 2>({a_angmom, b_angmom});

    const auto rs_pair = std::array<int, 2>({c_angmom, d_angmom});

    const auto pr_pair = std::array<int, 2>({a_angmom, c_angmom});

    const auto ps_pair = std::array<int, 2>({a_angmom, d_angmom});

    const auto qr_pair = std::array<int, 2>({b_angmom, c_angmom});

    const auto qs_pair = std::array<int, 2>({b_angmom, d_angmom});

    // set up Fock submatrices

    auto submat_pq = fock->sub_matrix(pq_pair);

    auto submat_rs = fock->sub_matrix(rs_pair);

    auto submat_pr = fock->sub_matrix(pr_pair);

    auto submat_ps = fock->sub_matrix(ps_pair);

    auto submat_qr = fock->sub_matrix(qr_pair);

    auto submat_qs = fock->sub_matrix(qs_pair);

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    // set up angular orders

    const auto angord_pq = fock->is_angular_order(pq_pair);

    const auto angord_rs = fock->is_angular_order(rs_pair);

    const auto angord_pr = fock->is_angular_order(pr_pair);

    const auto angord_ps = fock->is_angular_order(ps_pair);

    const auto angord_qr = fock->is_angular_order(qr_pair);

    const auto angord_qs = fock->is_angular_order(qs_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(a_angmom);

    const auto bcomps = tensor::number_of_spherical_components(b_angmom);

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            // impose angular symmetry on bra side

            if (refp == refq)
            {
                if (j < i) continue;
            }

            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer[offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l];

                    for (int m = ket_range[0]; m < ket_range[1]; m++)
                    {
                        // skip repeating integrals in diagonal block

                        if (diagonal)
                        {
                            if (m < bra_igto) continue;
                        }

                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

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

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range[0]];

                        if (p == q) fval *= 0.5;

                        if (r == s) fval *= 0.5;

                        if ((p == r) && (q == s)) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        if (angord_pq)
                        {
                            if (angord_rs)
                            {
                                submat_pq->operator[]({p, q}) += 4.0 * fval * denmat_rs->operator[]({r, s});
                            }
                            else
                            {
                                submat_pq->operator[]({p, q}) += 4.0 * fval * denmat_rs->operator[]({s, r});
                            }
                        }
                        else
                        {
                            if (angord_rs)
                            {
                                submat_pq->operator[]({q, p}) += 4.0 * fval * denmat_rs->operator[]({r, s});
                            }
                            else
                            {
                                submat_pq->operator[]({q, p}) += 4.0 * fval * denmat_rs->operator[]({s, r});
                            }
                        }

                        // Coulomb contribution (F_rs)

                        if (angord_rs)
                        {
                            if (angord_pq)
                            {
                                submat_rs->operator[]({r, s}) += 4.0 * fval * denmat_pq->operator[]({p, q});
                            }
                            else
                            {
                                submat_rs->operator[]({r, s}) += 4.0 * fval * denmat_pq->operator[]({q, p});
                            }
                        }
                        else
                        {
                            if (angord_pq)
                            {
                                submat_rs->operator[]({s, r}) += 4.0 * fval * denmat_pq->operator[]({p, q});
                            }
                            else
                            {
                                submat_rs->operator[]({s, r}) += 4.0 * fval * denmat_pq->operator[]({q, p});
                            }
                        }

                        // Exchange contribution (F_pr)

                        if (angord_pr)
                        {
                            if (angord_qs)
                            {
                                submat_pr->operator[]({p, r}) -= fval * denmat_qs->operator[]({q, s});
                            }
                            else
                            {
                                submat_pr->operator[]({p, r}) -= fval * denmat_qs->operator[]({s, q});
                            }
                        }
                        else
                        {
                            if (angord_qs)
                            {
                                submat_pr->operator[]({r, p}) -= fval * denmat_qs->operator[]({q, s});
                            }
                            else
                            {
                                submat_pr->operator[]({r, p}) -= fval * denmat_qs->operator[]({s, q});
                            }
                        }

                        // Exchange contribution (F_ps)

                        if (angord_ps)
                        {
                            if (angord_qr)
                            {
                                submat_ps->operator[]({p, s}) -= fval * denmat_qr->operator[]({q, r});
                            }
                            else
                            {
                                submat_ps->operator[]({p, s}) -= fval * denmat_qr->operator[]({r, q});
                            }
                        }
                        else
                        {
                            if (angord_qr)
                            {
                                submat_ps->operator[]({s, p}) -= fval * denmat_qr->operator[]({q, r});
                            }
                            else
                            {
                                submat_ps->operator[]({s, p}) -= fval * denmat_qr->operator[]({r, q});
                            }
                        }

                        // Exchange contribution (F_qr)

                        if (angord_qr)
                        {
                            if (angord_ps)
                            {
                                submat_qr->operator[]({q, r}) -= fval * denmat_ps->operator[]({p, s});
                            }
                            else
                            {
                                submat_qr->operator[]({q, r}) -= fval * denmat_ps->operator[]({s, p});
                            }
                        }
                        else
                        {
                            if (angord_ps)
                            {
                                submat_qr->operator[]({r, q}) -= fval * denmat_ps->operator[]({p, s});
                            }
                            else
                            {
                                submat_qr->operator[]({r, q}) -= fval * denmat_ps->operator[]({s, p});
                            }
                        }

                        // Exchange contribution (F_qs)

                        if (angord_qs)
                        {
                            if (angord_pr)
                            {
                                submat_qs->operator[]({q, s}) -= fval * denmat_pr->operator[]({p, r});
                            }
                            else
                            {
                                submat_qs->operator[]({q, s}) -= fval * denmat_pr->operator[]({r, p});
                            }
                        }
                        else
                        {
                            if (angord_pr)
                            {
                                submat_qs->operator[]({s, q}) -= fval * denmat_pr->operator[]({p, r});
                            }
                            else
                            {
                                submat_qs->operator[]({s, q}) -= fval * denmat_pr->operator[]({r, p});
                            }
                        }
                    }
                }
            }
        }
    }
}

auto
distribute_rest_jkx(CMatrix*                  fock,
                    const CMatrix*            density,
                    const CSimdArray<double>& buffer,
                    const int                 offset,
                    const double              factor,
                    const std::vector<int>&   a_indices,
                    const std::vector<int>&   b_indices,
                    const std::vector<int>&   c_indices,
                    const std::vector<int>&   d_indices,
                    const int                 a_angmom,
                    const int                 b_angmom,
                    const int                 c_angmom,
                    const int                 d_angmom,
                    const int                 bra_igto,
                    const std::array<int, 2>& ket_range,
                    const bool                diagonal) -> void
{
    // set up angular pairs

    const auto pq_pair = std::array<int, 2>({a_angmom, b_angmom});

    const auto rs_pair = std::array<int, 2>({c_angmom, d_angmom});

    const auto pr_pair = std::array<int, 2>({a_angmom, c_angmom});

    const auto ps_pair = std::array<int, 2>({a_angmom, d_angmom});

    const auto qr_pair = std::array<int, 2>({b_angmom, c_angmom});

    const auto qs_pair = std::array<int, 2>({b_angmom, d_angmom});

    // set up Fock submatrices

    auto submat_pq = fock->sub_matrix(pq_pair);

    auto submat_rs = fock->sub_matrix(rs_pair);

    auto submat_pr = fock->sub_matrix(pr_pair);

    auto submat_ps = fock->sub_matrix(ps_pair);

    auto submat_qr = fock->sub_matrix(qr_pair);

    auto submat_qs = fock->sub_matrix(qs_pair);

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    // set up angular orders

    const auto angord_pq = fock->is_angular_order(pq_pair);

    const auto angord_rs = fock->is_angular_order(rs_pair);

    const auto angord_pr = fock->is_angular_order(pr_pair);

    const auto angord_ps = fock->is_angular_order(ps_pair);

    const auto angord_qr = fock->is_angular_order(qr_pair);

    const auto angord_qs = fock->is_angular_order(qs_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(a_angmom);

    const auto bcomps = tensor::number_of_spherical_components(b_angmom);

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            // impose angular symmetry on bra side

            if (refp == refq)
            {
                if (j < i) continue;
            }

            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer[offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l];

                    for (int m = ket_range[0]; m < ket_range[1]; m++)
                    {
                        // skip repeating integrals in diagonal block

                        if (diagonal)
                        {
                            if (m < bra_igto) continue;
                        }

                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

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

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range[0]];

                        if (p == q) fval *= 0.5;

                        if (r == s) fval *= 0.5;

                        if ((p == r) && (q == s)) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        if (angord_pq)
                        {
                            if (angord_rs)
                            {
                                submat_pq->operator[]({p, q}) += 4.0 * fval * denmat_rs->operator[]({r, s});
                            }
                            else
                            {
                                submat_pq->operator[]({p, q}) += 4.0 * fval * denmat_rs->operator[]({s, r});
                            }
                        }
                        else
                        {
                            if (angord_rs)
                            {
                                submat_pq->operator[]({q, p}) += 4.0 * fval * denmat_rs->operator[]({r, s});
                            }
                            else
                            {
                                submat_pq->operator[]({q, p}) += 4.0 * fval * denmat_rs->operator[]({s, r});
                            }
                        }

                        // Coulomb contribution (F_rs)

                        if (angord_rs)
                        {
                            if (angord_pq)
                            {
                                submat_rs->operator[]({r, s}) += 4.0 * fval * denmat_pq->operator[]({p, q});
                            }
                            else
                            {
                                submat_rs->operator[]({r, s}) += 4.0 * fval * denmat_pq->operator[]({q, p});
                            }
                        }
                        else
                        {
                            if (angord_pq)
                            {
                                submat_rs->operator[]({s, r}) += 4.0 * fval * denmat_pq->operator[]({p, q});
                            }
                            else
                            {
                                submat_rs->operator[]({s, r}) += 4.0 * fval * denmat_pq->operator[]({q, p});
                            }
                        }

                        // rescale exchange contribution

                        fval *= factor;

                        // Exchange contribution (F_pr)

                        if (angord_pr)
                        {
                            if (angord_qs)
                            {
                                submat_pr->operator[]({p, r}) -= fval * denmat_qs->operator[]({q, s});
                            }
                            else
                            {
                                submat_pr->operator[]({p, r}) -= fval * denmat_qs->operator[]({s, q});
                            }
                        }
                        else
                        {
                            if (angord_qs)
                            {
                                submat_pr->operator[]({r, p}) -= fval * denmat_qs->operator[]({q, s});
                            }
                            else
                            {
                                submat_pr->operator[]({r, p}) -= fval * denmat_qs->operator[]({s, q});
                            }
                        }

                        // Exchange contribution (F_ps)

                        if (angord_ps)
                        {
                            if (angord_qr)
                            {
                                submat_ps->operator[]({p, s}) -= fval * denmat_qr->operator[]({q, r});
                            }
                            else
                            {
                                submat_ps->operator[]({p, s}) -= fval * denmat_qr->operator[]({r, q});
                            }
                        }
                        else
                        {
                            if (angord_qr)
                            {
                                submat_ps->operator[]({s, p}) -= fval * denmat_qr->operator[]({q, r});
                            }
                            else
                            {
                                submat_ps->operator[]({s, p}) -= fval * denmat_qr->operator[]({r, q});
                            }
                        }

                        // Exchange contribution (F_qr)

                        if (angord_qr)
                        {
                            if (angord_ps)
                            {
                                submat_qr->operator[]({q, r}) -= fval * denmat_ps->operator[]({p, s});
                            }
                            else
                            {
                                submat_qr->operator[]({q, r}) -= fval * denmat_ps->operator[]({s, p});
                            }
                        }
                        else
                        {
                            if (angord_ps)
                            {
                                submat_qr->operator[]({r, q}) -= fval * denmat_ps->operator[]({p, s});
                            }
                            else
                            {
                                submat_qr->operator[]({r, q}) -= fval * denmat_ps->operator[]({s, p});
                            }
                        }

                        // Exchange contribution (F_qs)

                        if (angord_qs)
                        {
                            if (angord_pr)
                            {
                                submat_qs->operator[]({q, s}) -= fval * denmat_pr->operator[]({p, r});
                            }
                            else
                            {
                                submat_qs->operator[]({q, s}) -= fval * denmat_pr->operator[]({r, p});
                            }
                        }
                        else
                        {
                            if (angord_pr)
                            {
                                submat_qs->operator[]({s, q}) -= fval * denmat_pr->operator[]({p, r});
                            }
                            else
                            {
                                submat_qs->operator[]({s, q}) -= fval * denmat_pr->operator[]({r, p});
                            }
                        }
                    }
                }
            }
        }
    }
}

auto
distribute_rest_j(CMatrix*                  fock,
                  const CMatrix*            density,
                  const CSimdArray<double>& buffer,
                  const int                 offset,
                  const std::vector<int>&   a_indices,
                  const std::vector<int>&   b_indices,
                  const std::vector<int>&   c_indices,
                  const std::vector<int>&   d_indices,
                  const int                 a_angmom,
                  const int                 b_angmom,
                  const int                 c_angmom,
                  const int                 d_angmom,
                  const int                 bra_igto,
                  const std::array<int, 2>& ket_range,
                  const bool                diagonal) -> void
{
    // set up angular pairs

    const auto pq_pair = std::array<int, 2>({a_angmom, b_angmom});

    const auto rs_pair = std::array<int, 2>({c_angmom, d_angmom});

    // set up Fock submatrices

    auto submat_pq = fock->sub_matrix(pq_pair);

    auto submat_rs = fock->sub_matrix(rs_pair);

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    // set up angular orders

    const auto angord_pq = fock->is_angular_order(pq_pair);

    const auto angord_rs = fock->is_angular_order(rs_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(a_angmom);

    const auto bcomps = tensor::number_of_spherical_components(b_angmom);

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            // impose angular symmetry on bra side

            if (refp == refq)
            {
                if (j < i) continue;
            }

            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer[offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l];

                    for (int m = ket_range[0]; m < ket_range[1]; m++)
                    {
                        // skip repeating integrals in diagonal block

                        if (diagonal)
                        {
                            if (m < bra_igto) continue;
                        }

                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

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

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range[0]];

                        if (p == q) fval *= 0.5;

                        if (r == s) fval *= 0.5;

                        if ((p == r) && (q == s)) fval *= 0.5;

                        // Coulomb contribution (F_pq)

                        if (angord_pq)
                        {
                            if (angord_rs)
                            {
                                submat_pq->operator[]({p, q}) += 2.0 * fval * denmat_rs->operator[]({r, s});
                            }
                            else
                            {
                                submat_pq->operator[]({p, q}) += 2.0 * fval * denmat_rs->operator[]({s, r});
                            }
                        }
                        else
                        {
                            if (angord_rs)
                            {
                                submat_pq->operator[]({q, p}) += 2.0 * fval * denmat_rs->operator[]({r, s});
                            }
                            else
                            {
                                submat_pq->operator[]({q, p}) += 2.0 * fval * denmat_rs->operator[]({s, r});
                            }
                        }

                        // Coulomb contribution (F_rs)

                        if (angord_rs)
                        {
                            if (angord_pq)
                            {
                                submat_rs->operator[]({r, s}) += 2.0 * fval * denmat_pq->operator[]({p, q});
                            }
                            else
                            {
                                submat_rs->operator[]({r, s}) += 2.0 * fval * denmat_pq->operator[]({q, p});
                            }
                        }
                        else
                        {
                            if (angord_pq)
                            {
                                submat_rs->operator[]({s, r}) += 2.0 * fval * denmat_pq->operator[]({p, q});
                            }
                            else
                            {
                                submat_rs->operator[]({s, r}) += 2.0 * fval * denmat_pq->operator[]({q, p});
                            }
                        }
                    }
                }
            }
        }
    }
}

auto
distribute_rest_k(CMatrix*                  fock,
                  const CMatrix*            density,
                  const CSimdArray<double>& buffer,
                  const int                 offset,
                  const std::vector<int>&   a_indices,
                  const std::vector<int>&   b_indices,
                  const std::vector<int>&   c_indices,
                  const std::vector<int>&   d_indices,
                  const int                 a_angmom,
                  const int                 b_angmom,
                  const int                 c_angmom,
                  const int                 d_angmom,
                  const int                 bra_igto,
                  const std::array<int, 2>& ket_range,
                  const bool                diagonal) -> void
{
    // set up angular pairs

    const auto pr_pair = std::array<int, 2>({a_angmom, c_angmom});

    const auto ps_pair = std::array<int, 2>({a_angmom, d_angmom});

    const auto qr_pair = std::array<int, 2>({b_angmom, c_angmom});

    const auto qs_pair = std::array<int, 2>({b_angmom, d_angmom});

    // set up Fock submatrices

    auto submat_pr = fock->sub_matrix(pr_pair);

    auto submat_ps = fock->sub_matrix(ps_pair);

    auto submat_qr = fock->sub_matrix(qr_pair);

    auto submat_qs = fock->sub_matrix(qs_pair);

    // set up AO density submatrix

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    // set up angular orders

    const auto angord_pr = fock->is_angular_order(pr_pair);

    const auto angord_ps = fock->is_angular_order(ps_pair);

    const auto angord_qr = fock->is_angular_order(qr_pair);

    const auto angord_qs = fock->is_angular_order(qs_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(a_angmom);

    const auto bcomps = tensor::number_of_spherical_components(b_angmom);

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            // impose angular symmetry on bra side

            if (refp == refq)
            {
                if (j < i) continue;
            }

            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer[offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l];

                    for (int m = ket_range[0]; m < ket_range[1]; m++)
                    {
                        // skip repeating integrals in diagonal block

                        if (diagonal)
                        {
                            if (m < bra_igto) continue;
                        }

                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

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

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range[0]];

                        if (p == q) fval *= 0.5;

                        if (r == s) fval *= 0.5;

                        if ((p == r) && (q == s)) fval *= 0.5;

                        // Exchange contribution (F_pr)

                        if (angord_pr)
                        {
                            if (angord_qs)
                            {
                                submat_pr->operator[]({p, r}) += fval * denmat_qs->operator[]({q, s});
                            }
                            else
                            {
                                submat_pr->operator[]({p, r}) += fval * denmat_qs->operator[]({s, q});
                            }
                        }
                        else
                        {
                            if (angord_qs)
                            {
                                submat_pr->operator[]({r, p}) += fval * denmat_qs->operator[]({q, s});
                            }
                            else
                            {
                                submat_pr->operator[]({r, p}) += fval * denmat_qs->operator[]({s, q});
                            }
                        }

                        // Exchange contribution (F_ps)

                        if (angord_ps)
                        {
                            if (angord_qr)
                            {
                                submat_ps->operator[]({p, s}) += fval * denmat_qr->operator[]({q, r});
                            }
                            else
                            {
                                submat_ps->operator[]({p, s}) += fval * denmat_qr->operator[]({r, q});
                            }
                        }
                        else
                        {
                            if (angord_qr)
                            {
                                submat_ps->operator[]({s, p}) += fval * denmat_qr->operator[]({q, r});
                            }
                            else
                            {
                                submat_ps->operator[]({s, p}) += fval * denmat_qr->operator[]({r, q});
                            }
                        }

                        // Exchange contribution (F_qr)

                        if (angord_qr)
                        {
                            if (angord_ps)
                            {
                                submat_qr->operator[]({q, r}) += fval * denmat_ps->operator[]({p, s});
                            }
                            else
                            {
                                submat_qr->operator[]({q, r}) += fval * denmat_ps->operator[]({s, p});
                            }
                        }
                        else
                        {
                            if (angord_ps)
                            {
                                submat_qr->operator[]({r, q}) += fval * denmat_ps->operator[]({p, s});
                            }
                            else
                            {
                                submat_qr->operator[]({r, q}) += fval * denmat_ps->operator[]({s, p});
                            }
                        }

                        // Exchange contribution (F_qs)

                        if (angord_qs)
                        {
                            if (angord_pr)
                            {
                                submat_qs->operator[]({q, s}) += fval * denmat_pr->operator[]({p, r});
                            }
                            else
                            {
                                submat_qs->operator[]({q, s}) += fval * denmat_pr->operator[]({r, p});
                            }
                        }
                        else
                        {
                            if (angord_pr)
                            {
                                submat_qs->operator[]({s, q}) += fval * denmat_pr->operator[]({p, r});
                            }
                            else
                            {
                                submat_qs->operator[]({s, q}) += fval * denmat_pr->operator[]({r, p});
                            }
                        }
                    }
                }
            }
        }
    }
}

auto
distribute_rest_kx(CMatrix*                  fock,
                   const CMatrix*            density,
                   const CSimdArray<double>& buffer,
                   const int                 offset,
                   const double              factor,
                   const std::vector<int>&   a_indices,
                   const std::vector<int>&   b_indices,
                   const std::vector<int>&   c_indices,
                   const std::vector<int>&   d_indices,
                   const int                 a_angmom,
                   const int                 b_angmom,
                   const int                 c_angmom,
                   const int                 d_angmom,
                   const int                 bra_igto,
                   const std::array<int, 2>& ket_range,
                   const bool                diagonal) -> void
{
    // set up angular pairs

    const auto pr_pair = std::array<int, 2>({a_angmom, c_angmom});

    const auto ps_pair = std::array<int, 2>({a_angmom, d_angmom});

    const auto qr_pair = std::array<int, 2>({b_angmom, c_angmom});

    const auto qs_pair = std::array<int, 2>({b_angmom, d_angmom});

    // set up Fock submatrices

    auto submat_pr = fock->sub_matrix(pr_pair);

    auto submat_ps = fock->sub_matrix(ps_pair);

    auto submat_qr = fock->sub_matrix(qr_pair);

    auto submat_qs = fock->sub_matrix(qs_pair);

    // set up AO density submatrix

    auto denmat_pr = density->sub_matrix(pr_pair);

    auto denmat_ps = density->sub_matrix(ps_pair);

    auto denmat_qr = density->sub_matrix(qr_pair);

    auto denmat_qs = density->sub_matrix(qs_pair);

    // set up angular orders

    const auto angord_pr = fock->is_angular_order(pr_pair);

    const auto angord_ps = fock->is_angular_order(ps_pair);

    const auto angord_qr = fock->is_angular_order(qr_pair);

    const auto angord_qs = fock->is_angular_order(qs_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(a_angmom);

    const auto bcomps = tensor::number_of_spherical_components(b_angmom);

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            // impose angular symmetry on bra side

            if (refp == refq)
            {
                if (j < i) continue;
            }

            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer[offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l];

                    for (int m = ket_range[0]; m < ket_range[1]; m++)
                    {
                        // skip repeating integrals in diagonal block

                        if (diagonal)
                        {
                            if (m < bra_igto) continue;
                        }

                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

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

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = factor * curr_buffer[m - ket_range[0]];

                        if (p == q) fval *= 0.5;

                        if (r == s) fval *= 0.5;

                        if ((p == r) && (q == s)) fval *= 0.5;

                        // Exchange contribution (F_pr)

                        if (angord_pr)
                        {
                            if (angord_qs)
                            {
                                submat_pr->operator[]({p, r}) += fval * denmat_qs->operator[]({q, s});
                            }
                            else
                            {
                                submat_pr->operator[]({p, r}) += fval * denmat_qs->operator[]({s, q});
                            }
                        }
                        else
                        {
                            if (angord_qs)
                            {
                                submat_pr->operator[]({r, p}) += fval * denmat_qs->operator[]({q, s});
                            }
                            else
                            {
                                submat_pr->operator[]({r, p}) += fval * denmat_qs->operator[]({s, q});
                            }
                        }

                        // Exchange contribution (F_ps)

                        if (angord_ps)
                        {
                            if (angord_qr)
                            {
                                submat_ps->operator[]({p, s}) += fval * denmat_qr->operator[]({q, r});
                            }
                            else
                            {
                                submat_ps->operator[]({p, s}) += fval * denmat_qr->operator[]({r, q});
                            }
                        }
                        else
                        {
                            if (angord_qr)
                            {
                                submat_ps->operator[]({s, p}) += fval * denmat_qr->operator[]({q, r});
                            }
                            else
                            {
                                submat_ps->operator[]({s, p}) += fval * denmat_qr->operator[]({r, q});
                            }
                        }

                        // Exchange contribution (F_qr)

                        if (angord_qr)
                        {
                            if (angord_ps)
                            {
                                submat_qr->operator[]({q, r}) += fval * denmat_ps->operator[]({p, s});
                            }
                            else
                            {
                                submat_qr->operator[]({q, r}) += fval * denmat_ps->operator[]({s, p});
                            }
                        }
                        else
                        {
                            if (angord_ps)
                            {
                                submat_qr->operator[]({r, q}) += fval * denmat_ps->operator[]({p, s});
                            }
                            else
                            {
                                submat_qr->operator[]({r, q}) += fval * denmat_ps->operator[]({s, p});
                            }
                        }

                        // Exchange contribution (F_qs)

                        if (angord_qs)
                        {
                            if (angord_pr)
                            {
                                submat_qs->operator[]({q, s}) += fval * denmat_pr->operator[]({p, r});
                            }
                            else
                            {
                                submat_qs->operator[]({q, s}) += fval * denmat_pr->operator[]({r, p});
                            }
                        }
                        else
                        {
                            if (angord_pr)
                            {
                                submat_qs->operator[]({s, q}) += fval * denmat_pr->operator[]({p, r});
                            }
                            else
                            {
                                submat_qs->operator[]({s, q}) += fval * denmat_pr->operator[]({r, p});
                            }
                        }
                    }
                }
            }
        }
    }
}

auto
distribute_rest_gen_j(CMatrix*                  fock,
                      const CMatrix*            density,
                      const CSimdArray<double>& buffer,
                      const int                 offset,
                      const std::vector<int>&   a_indices,
                      const std::vector<int>&   b_indices,
                      const std::vector<int>&   c_indices,
                      const std::vector<int>&   d_indices,
                      const int                 a_angmom,
                      const int                 b_angmom,
                      const int                 c_angmom,
                      const int                 d_angmom,
                      const int                 bra_igto,
                      const std::array<int, 2>& ket_range,
                      const bool                diagonal) -> void
{
    // set up angular pairs

    const auto pq_pair = std::array<int, 2>({a_angmom, b_angmom});

    const auto qp_pair = std::array<int, 2>({b_angmom, a_angmom});

    const auto rs_pair = std::array<int, 2>({c_angmom, d_angmom});

    const auto sr_pair = std::array<int, 2>({d_angmom, c_angmom});

    // set up Fock submatrices

    auto submat_pq = fock->sub_matrix(pq_pair);

    auto submat_qp = fock->sub_matrix(qp_pair);

    auto submat_rs = fock->sub_matrix(rs_pair);

    auto submat_sr = fock->sub_matrix(sr_pair);

    // set up AO density submatrix

    auto denmat_pq = density->sub_matrix(pq_pair);

    auto denmat_qp = density->sub_matrix(qp_pair);

    auto denmat_rs = density->sub_matrix(rs_pair);

    auto denmat_sr = density->sub_matrix(sr_pair);

    // reference indexes on bra side

    const auto refp = a_indices[bra_igto + 1];

    const auto refq = b_indices[bra_igto + 1];

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(a_angmom);

    const auto bcomps = tensor::number_of_spherical_components(b_angmom);

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            // impose angular symmetry on bra side

            if (refp == refq)
            {
                if (j < i) continue;
            }

            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer[offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l];

                    for (int m = ket_range[0]; m < ket_range[1]; m++)
                    {
                        // skip repeating integrals in diagonal block

                        if (diagonal)
                        {
                            if (m < bra_igto) continue;
                        }

                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

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

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range[0]];

                        if (p == q) fval *= 0.5;

                        if (r == s) fval *= 0.5;

                        if ((p == r) && (q == s)) fval *= 0.5;

                        // precomputed integrals

                        const auto f2rs = fval * (denmat_rs->operator[]({r, s}) + denmat_sr->operator[]({s, r}));

                        const auto f2pq = fval * (denmat_pq->operator[]({p, q}) + denmat_qp->operator[]({q, p}));

                        // Coulomb contributions

                        submat_pq->operator[]({p, q}) += f2rs;

                        submat_qp->operator[]({q, p}) += f2rs;

                        submat_rs->operator[]({r, s}) += f2pq;

                        submat_sr->operator[]({s, r}) += f2pq;
                    }
                }
            }
        }
    }
}

auto
distribute_rest_gen_k(CMatrix*                  fock,
                      const CMatrix*            density,
                      const CSimdArray<double>& buffer,
                      const int                 offset,
                      const std::vector<int>&   a_indices,
                      const std::vector<int>&   b_indices,
                      const std::vector<int>&   c_indices,
                      const std::vector<int>&   d_indices,
                      const int                 a_angmom,
                      const int                 b_angmom,
                      const int                 c_angmom,
                      const int                 d_angmom,
                      const int                 bra_igto,
                      const std::array<int, 2>& ket_range,
                      const bool                diagonal) -> void
{
    // set up angular pairs

    const auto pr_pair = std::array<int, 2>({a_angmom, c_angmom});

    const auto rp_pair = std::array<int, 2>({c_angmom, a_angmom});

    const auto ps_pair = std::array<int, 2>({a_angmom, d_angmom});

    const auto sp_pair = std::array<int, 2>({d_angmom, a_angmom});

    const auto qr_pair = std::array<int, 2>({b_angmom, c_angmom});

    const auto rq_pair = std::array<int, 2>({c_angmom, b_angmom});

    const auto qs_pair = std::array<int, 2>({b_angmom, d_angmom});

    const auto sq_pair = std::array<int, 2>({d_angmom, b_angmom});

    // set up Fock submatrices

    auto submat_pr = fock->sub_matrix(pr_pair);

    auto submat_rp = fock->sub_matrix(rp_pair);

    auto submat_ps = fock->sub_matrix(ps_pair);

    auto submat_sp = fock->sub_matrix(sp_pair);

    auto submat_qr = fock->sub_matrix(qr_pair);

    auto submat_rq = fock->sub_matrix(rq_pair);

    auto submat_qs = fock->sub_matrix(qs_pair);

    auto submat_sq = fock->sub_matrix(sq_pair);

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

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(a_angmom);

    const auto bcomps = tensor::number_of_spherical_components(b_angmom);

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            // impose angular symmetry on bra side

            if (refp == refq)
            {
                if (j < i) continue;
            }

            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer[offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l];

                    for (int m = ket_range[0]; m < ket_range[1]; m++)
                    {
                        // skip repeating integrals in diagonal block

                        if (diagonal)
                        {
                            if (m < bra_igto) continue;
                        }

                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

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

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range[0]];

                        if (p == q) fval *= 0.5;

                        if (r == s) fval *= 0.5;

                        if ((p == r) && (q == s)) fval *= 0.5;

                        // Exchange contributions

                        submat_pr->operator[]({p, r}) += fval * denmat_qs->operator[]({q, s});

                        submat_rp->operator[]({r, p}) += fval * denmat_sq->operator[]({s, q});

                        submat_ps->operator[]({p, s}) += fval * denmat_qr->operator[]({q, r});

                        submat_sp->operator[]({s, p}) += fval * denmat_rq->operator[]({r, q});

                        submat_qr->operator[]({q, r}) += fval * denmat_ps->operator[]({p, s});

                        submat_rq->operator[]({r, q}) += fval * denmat_sp->operator[]({s, p});

                        submat_qs->operator[]({q, s}) += fval * denmat_pr->operator[]({p, r});

                        submat_sq->operator[]({s, q}) += fval * denmat_rp->operator[]({r, p});
                    }
                }
            }
        }
    }
}

auto
distribute_rest_gen_kx(CMatrix*                  fock,
                       const CMatrix*            density,
                       const CSimdArray<double>& buffer,
                       const int                 offset,
                       const double              factor,
                       const std::vector<int>&   a_indices,
                       const std::vector<int>&   b_indices,
                       const std::vector<int>&   c_indices,
                       const std::vector<int>&   d_indices,
                       const int                 a_angmom,
                       const int                 b_angmom,
                       const int                 c_angmom,
                       const int                 d_angmom,
                       const int                 bra_igto,
                       const std::array<int, 2>& ket_range,
                       const bool                diagonal) -> void
{
    // set up angular pairs

    const auto pr_pair = std::array<int, 2>({a_angmom, c_angmom});

    const auto rp_pair = std::array<int, 2>({c_angmom, a_angmom});

    const auto ps_pair = std::array<int, 2>({a_angmom, d_angmom});

    const auto sp_pair = std::array<int, 2>({d_angmom, a_angmom});

    const auto qr_pair = std::array<int, 2>({b_angmom, c_angmom});

    const auto rq_pair = std::array<int, 2>({c_angmom, b_angmom});

    const auto qs_pair = std::array<int, 2>({b_angmom, d_angmom});

    const auto sq_pair = std::array<int, 2>({d_angmom, b_angmom});

    // set up Fock submatrices

    auto submat_pr = fock->sub_matrix(pr_pair);

    auto submat_rp = fock->sub_matrix(rp_pair);

    auto submat_ps = fock->sub_matrix(ps_pair);

    auto submat_sp = fock->sub_matrix(sp_pair);

    auto submat_qr = fock->sub_matrix(qr_pair);

    auto submat_rq = fock->sub_matrix(rq_pair);

    auto submat_qs = fock->sub_matrix(qs_pair);

    auto submat_sq = fock->sub_matrix(sq_pair);

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

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(a_angmom);

    const auto bcomps = tensor::number_of_spherical_components(b_angmom);

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            // impose angular symmetry on bra side

            if (refp == refq)
            {
                if (j < i) continue;
            }

            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer[offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l];

                    for (int m = ket_range[0]; m < ket_range[1]; m++)
                    {
                        // skip repeating integrals in diagonal block

                        if (diagonal)
                        {
                            if (m < bra_igto) continue;
                        }

                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

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

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = factor * curr_buffer[m - ket_range[0]];

                        if (p == q) fval *= 0.5;

                        if (r == s) fval *= 0.5;

                        if ((p == r) && (q == s)) fval *= 0.5;

                        // Exchange contributions

                        submat_pr->operator[]({p, r}) += fval * denmat_qs->operator[]({q, s});

                        submat_rp->operator[]({r, p}) += fval * denmat_sq->operator[]({s, q});

                        submat_ps->operator[]({p, s}) += fval * denmat_qr->operator[]({q, r});

                        submat_sp->operator[]({s, p}) += fval * denmat_rq->operator[]({r, q});

                        submat_qr->operator[]({q, r}) += fval * denmat_ps->operator[]({p, s});

                        submat_rq->operator[]({r, q}) += fval * denmat_sp->operator[]({s, p});

                        submat_qs->operator[]({q, s}) += fval * denmat_pr->operator[]({p, r});

                        submat_sq->operator[]({s, q}) += fval * denmat_rp->operator[]({r, p});
                    }
                }
            }
        }
    }
}

auto
distribute_rest_gen_jk(CMatrix*                  fock,
                       const CMatrix*            density,
                       const CSimdArray<double>& buffer,
                       const int                 offset,
                       const std::vector<int>&   a_indices,
                       const std::vector<int>&   b_indices,
                       const std::vector<int>&   c_indices,
                       const std::vector<int>&   d_indices,
                       const int                 a_angmom,
                       const int                 b_angmom,
                       const int                 c_angmom,
                       const int                 d_angmom,
                       const int                 bra_igto,
                       const std::array<int, 2>& ket_range,
                       const bool                diagonal) -> void
{
    // set up angular pairs

    const auto pq_pair = std::array<int, 2>({a_angmom, b_angmom});

    const auto qp_pair = std::array<int, 2>({b_angmom, a_angmom});

    const auto rs_pair = std::array<int, 2>({c_angmom, d_angmom});

    const auto sr_pair = std::array<int, 2>({d_angmom, c_angmom});

    const auto pr_pair = std::array<int, 2>({a_angmom, c_angmom});

    const auto rp_pair = std::array<int, 2>({c_angmom, a_angmom});

    const auto ps_pair = std::array<int, 2>({a_angmom, d_angmom});

    const auto sp_pair = std::array<int, 2>({d_angmom, a_angmom});

    const auto qr_pair = std::array<int, 2>({b_angmom, c_angmom});

    const auto rq_pair = std::array<int, 2>({c_angmom, b_angmom});

    const auto qs_pair = std::array<int, 2>({b_angmom, d_angmom});

    const auto sq_pair = std::array<int, 2>({d_angmom, b_angmom});

    // set up Fock submatrices

    auto submat_pq = fock->sub_matrix(pq_pair);

    auto submat_qp = fock->sub_matrix(qp_pair);

    auto submat_rs = fock->sub_matrix(rs_pair);

    auto submat_sr = fock->sub_matrix(sr_pair);

    auto submat_pr = fock->sub_matrix(pr_pair);

    auto submat_rp = fock->sub_matrix(rp_pair);

    auto submat_ps = fock->sub_matrix(ps_pair);

    auto submat_sp = fock->sub_matrix(sp_pair);

    auto submat_qr = fock->sub_matrix(qr_pair);

    auto submat_rq = fock->sub_matrix(rq_pair);

    auto submat_qs = fock->sub_matrix(qs_pair);

    auto submat_sq = fock->sub_matrix(sq_pair);

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

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(a_angmom);

    const auto bcomps = tensor::number_of_spherical_components(b_angmom);

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            // impose angular symmetry on bra side

            if (refp == refq)
            {
                if (j < i) continue;
            }

            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer[offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l];

                    for (int m = ket_range[0]; m < ket_range[1]; m++)
                    {
                        // skip repeating integrals in diagonal block

                        if (diagonal)
                        {
                            if (m < bra_igto) continue;
                        }

                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

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

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range[0]];

                        if (p == q) fval *= 0.5;

                        if (r == s) fval *= 0.5;

                        if ((p == r) && (q == s)) fval *= 0.5;

                        // precomputed integrals

                        const auto f2rs = 2.0 * fval * (denmat_rs->operator[]({r, s}) + denmat_sr->operator[]({s, r}));

                        const auto f2pq = 2.0 * fval * (denmat_pq->operator[]({p, q}) + denmat_qp->operator[]({q, p}));

                        // Coulomb contributions

                        submat_pq->operator[]({p, q}) += f2rs;

                        submat_qp->operator[]({q, p}) += f2rs;

                        submat_rs->operator[]({r, s}) += f2pq;

                        submat_sr->operator[]({s, r}) += f2pq;

                        // Exchange contributions

                        submat_pr->operator[]({p, r}) -= fval * denmat_qs->operator[]({q, s});

                        submat_rp->operator[]({r, p}) -= fval * denmat_sq->operator[]({s, q});

                        submat_ps->operator[]({p, s}) -= fval * denmat_qr->operator[]({q, r});

                        submat_sp->operator[]({s, p}) -= fval * denmat_rq->operator[]({r, q});

                        submat_qr->operator[]({q, r}) -= fval * denmat_ps->operator[]({p, s});

                        submat_rq->operator[]({r, q}) -= fval * denmat_sp->operator[]({s, p});

                        submat_qs->operator[]({q, s}) -= fval * denmat_pr->operator[]({p, r});

                        submat_sq->operator[]({s, q}) -= fval * denmat_rp->operator[]({r, p});
                    }
                }
            }
        }
    }
}

auto
distribute_rest_gen_jkx(CMatrix*                  fock,
                        const CMatrix*            density,
                        const CSimdArray<double>& buffer,
                        const int                 offset,
                        const double              factor,
                        const std::vector<int>&   a_indices,
                        const std::vector<int>&   b_indices,
                        const std::vector<int>&   c_indices,
                        const std::vector<int>&   d_indices,
                        const int                 a_angmom,
                        const int                 b_angmom,
                        const int                 c_angmom,
                        const int                 d_angmom,
                        const int                 bra_igto,
                        const std::array<int, 2>& ket_range,
                        const bool                diagonal) -> void
{
    // set up angular pairs

    const auto pq_pair = std::array<int, 2>({a_angmom, b_angmom});

    const auto qp_pair = std::array<int, 2>({b_angmom, a_angmom});

    const auto rs_pair = std::array<int, 2>({c_angmom, d_angmom});

    const auto sr_pair = std::array<int, 2>({d_angmom, c_angmom});

    const auto pr_pair = std::array<int, 2>({a_angmom, c_angmom});

    const auto rp_pair = std::array<int, 2>({c_angmom, a_angmom});

    const auto ps_pair = std::array<int, 2>({a_angmom, d_angmom});

    const auto sp_pair = std::array<int, 2>({d_angmom, a_angmom});

    const auto qr_pair = std::array<int, 2>({b_angmom, c_angmom});

    const auto rq_pair = std::array<int, 2>({c_angmom, b_angmom});

    const auto qs_pair = std::array<int, 2>({b_angmom, d_angmom});

    const auto sq_pair = std::array<int, 2>({d_angmom, b_angmom});

    // set up Fock submatrices

    auto submat_pq = fock->sub_matrix(pq_pair);

    auto submat_qp = fock->sub_matrix(qp_pair);

    auto submat_rs = fock->sub_matrix(rs_pair);

    auto submat_sr = fock->sub_matrix(sr_pair);

    auto submat_pr = fock->sub_matrix(pr_pair);

    auto submat_rp = fock->sub_matrix(rp_pair);

    auto submat_ps = fock->sub_matrix(ps_pair);

    auto submat_sp = fock->sub_matrix(sp_pair);

    auto submat_qr = fock->sub_matrix(qr_pair);

    auto submat_rq = fock->sub_matrix(rq_pair);

    auto submat_qs = fock->sub_matrix(qs_pair);

    auto submat_sq = fock->sub_matrix(sq_pair);

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

    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto bdim = b_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];

    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(a_angmom);

    const auto bcomps = tensor::number_of_spherical_components(b_angmom);

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;

        for (int j = 0; j < bcomps; j++)
        {
            // impose angular symmetry on bra side

            if (refp == refq)
            {
                if (j < i) continue;
            }

            const auto q = j * bdim + refq;

            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer[offset + i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l];

                    for (int m = ket_range[0]; m < ket_range[1]; m++)
                    {
                        // skip repeating integrals in diagonal block

                        if (diagonal)
                        {
                            if (m < bra_igto) continue;
                        }

                        // reference indexes on ket side

                        const auto refr = c_indices[m + 1];

                        const auto refs = d_indices[m + 1];

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

                        const auto s = l * ddim + refs;

                        // prescale integral for accumulation to Fock matrix

                        auto fval = curr_buffer[m - ket_range[0]];

                        if (p == q) fval *= 0.5;

                        if (r == s) fval *= 0.5;

                        if ((p == r) && (q == s)) fval *= 0.5;

                        // precomputed integrals

                        const auto f2rs = 2.0 * fval * (denmat_rs->operator[]({r, s}) + denmat_sr->operator[]({s, r}));

                        const auto f2pq = 2.0 * fval * (denmat_pq->operator[]({p, q}) + denmat_qp->operator[]({q, p}));

                        // Coulomb contributions

                        submat_pq->operator[]({p, q}) += f2rs;

                        submat_qp->operator[]({q, p}) += f2rs;

                        submat_rs->operator[]({r, s}) += f2pq;

                        submat_sr->operator[]({s, r}) += f2pq;

                        // Exchange contributions

                        fval *= factor;

                        submat_pr->operator[]({p, r}) -= fval * denmat_qs->operator[]({q, s});

                        submat_rp->operator[]({r, p}) -= fval * denmat_sq->operator[]({s, q});

                        submat_ps->operator[]({p, s}) -= fval * denmat_qr->operator[]({q, r});

                        submat_sp->operator[]({s, p}) -= fval * denmat_rq->operator[]({r, q});

                        submat_qr->operator[]({q, r}) -= fval * denmat_ps->operator[]({p, s});

                        submat_rq->operator[]({r, q}) -= fval * denmat_sp->operator[]({s, p});

                        submat_qs->operator[]({q, s}) -= fval * denmat_pr->operator[]({p, r});

                        submat_sq->operator[]({s, q}) -= fval * denmat_rp->operator[]({r, p});
                    }
                }
            }
        }
    }
}

}  // namespace t4cfunc
