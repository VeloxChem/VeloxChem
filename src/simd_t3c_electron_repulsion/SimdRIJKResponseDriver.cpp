//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#include "SimdRIJKResponseDriver.hpp"

#include <algorithm>
#include <string>

#include "ErrorHandler.hpp"
#include "ThreadedDenseLinearAlgebra.hpp"

auto
CSimdRIJKResponseDriver::get_threshold() const -> double
{
    return _threshold;
}

auto
CSimdRIJKResponseDriver::get_block_size() const -> size_t
{
    return _block_size;
}

auto
CSimdRIJKResponseDriver::get_memory_budget() const -> size_t
{
    return _budget;
}

auto
CSimdRIJKResponseDriver::_check_factors(const CMolecularBasis            &basis,
                                        const CPackedMatrix              &left,
                                        const std::vector<CPackedMatrix> &rights) const -> void
{
    const auto nao = basis.dimensions_of_basis();

    errors::assertMsgCritical(
        left.number_of_rows() == nao,
        std::string("SimdRIJKResponseDriver: The left factor is not of the molecular basis"));

    errors::assertMsgCritical(
        left.number_of_columns() > 0,
        std::string("SimdRIJKResponseDriver: The left factor has no columns"));

    for (const auto &right : rights)
    {
        errors::assertMsgCritical(
            right.number_of_rows() == nao,
            std::string("SimdRIJKResponseDriver: A right factor is not of the molecular basis"));

        // NOTE: the two factors close over the same index, so a density whose
        // factors disagree on it is not a density. Caught here rather than by the
        // multiply, which would read past one of them.
        errors::assertMsgCritical(
            right.number_of_columns() == left.number_of_columns(),
            std::string("SimdRIJKResponseDriver: A right factor is not of the rank of the left factor"));
    }
}

auto
CSimdRIJKResponseDriver::compute_exchange(const CSparseTensor              &bq_vectors,
                                          const CMolecularBasis            &basis,
                                          const CMolecularBasis            &aux_basis,
                                          const CPackedMatrix              &left,
                                          const std::vector<CPackedMatrix> &rights) const -> std::vector<CPackedMatrix>
{
    return _exchange({{&bq_vectors, 1.0}}, basis, aux_basis, left, rights);
}

auto
CSimdRIJKResponseDriver::compute_exchange(const CSparseTensor              &bq_vectors,
                                          const CSparseTensor              &bq_vectors_erf,
                                          const CMolecularBasis            &basis,
                                          const CMolecularBasis            &aux_basis,
                                          const CPackedMatrix              &left,
                                          const std::vector<CPackedMatrix> &rights,
                                          const double                      exchange_scaling_factor,
                                          const double                      erf_exchange_scaling_factor) const
    -> std::vector<CPackedMatrix>
{
    return _exchange({{&bq_vectors, exchange_scaling_factor}, {&bq_vectors_erf, erf_exchange_scaling_factor}},
                     basis,
                     aux_basis,
                     left,
                     rights);
}

auto
CSimdRIJKResponseDriver::compute_unrestricted(const CSparseTensor              &bq_vectors,
                                              const CMolecularBasis            &basis,
                                              const CMolecularBasis            &aux_basis,
                                              const CPackedMatrix              &left_alpha,
                                              const std::vector<CPackedMatrix> &rights_alpha,
                                              const std::vector<CPackedMatrix> &transposed_rights_alpha,
                                              const CPackedMatrix              &left_beta,
                                              const std::vector<CPackedMatrix> &rights_beta,
                                              const std::vector<CPackedMatrix> &transposed_rights_beta,
                                              const double                      exchange_scaling_factor) const
    -> std::pair<std::vector<CPackedMatrix>, std::vector<CPackedMatrix>>
{
    return _unrestricted({{&bq_vectors, exchange_scaling_factor}},
                         bq_vectors,
                         basis,
                         aux_basis,
                         left_alpha,
                         rights_alpha,
                         transposed_rights_alpha,
                         left_beta,
                         rights_beta,
                         transposed_rights_beta);
}

auto
CSimdRIJKResponseDriver::compute_unrestricted_rs(const CSparseTensor              &bq_vectors,
                                                 const CSparseTensor              &bq_vectors_erf,
                                                 const CMolecularBasis            &basis,
                                                 const CMolecularBasis            &aux_basis,
                                                 const CPackedMatrix              &left_alpha,
                                                 const std::vector<CPackedMatrix> &rights_alpha,
                                                 const std::vector<CPackedMatrix> &transposed_rights_alpha,
                                                 const CPackedMatrix              &left_beta,
                                                 const std::vector<CPackedMatrix> &rights_beta,
                                                 const std::vector<CPackedMatrix> &transposed_rights_beta,
                                                 const double                      exchange_scaling_factor,
                                                 const double erf_exchange_scaling_factor) const
    -> std::pair<std::vector<CPackedMatrix>, std::vector<CPackedMatrix>>
{
    // NOTE: the Coulomb term is formed from the plain B vectors whatever the
    // functional. Only the exchange is split between the two operators; the
    // Coulomb of a range separated functional is the whole of 1 / r.
    return _unrestricted({{&bq_vectors, exchange_scaling_factor}, {&bq_vectors_erf, erf_exchange_scaling_factor}},
                         bq_vectors,
                         basis,
                         aux_basis,
                         left_alpha,
                         rights_alpha,
                         transposed_rights_alpha,
                         left_beta,
                         rights_beta,
                         transposed_rights_beta);
}

auto
CSimdRIJKResponseDriver::_unrestricted(const std::vector<std::pair<const CSparseTensor *, double>> &operators,
                                       const CSparseTensor              &coulomb_vectors,
                                       const CMolecularBasis            &basis,
                                       const CMolecularBasis            &aux_basis,
                                       const CPackedMatrix              &left_alpha,
                                       const std::vector<CPackedMatrix> &rights_alpha,
                                       const std::vector<CPackedMatrix> &transposed_rights_alpha,
                                       const CPackedMatrix              &left_beta,
                                       const std::vector<CPackedMatrix> &rights_beta,
                                       const std::vector<CPackedMatrix> &transposed_rights_beta) const
    -> std::pair<std::vector<CPackedMatrix>, std::vector<CPackedMatrix>>
{
    auto focks_alpha = std::vector<CPackedMatrix>();

    auto focks_beta = std::vector<CPackedMatrix>();

    if (rights_alpha.empty()) return {std::move(focks_alpha), std::move(focks_beta)};

    errors::assertMsgCritical(rights_alpha.size() == rights_beta.size(),
                              std::string("SimdRIJKResponseDriver: The two spins do not carry the same number of "
                                          "densities"));

    _check_factors(basis, left_alpha, rights_alpha);

    _check_factors(basis, left_beta, rights_beta);

    const auto two_termed_alpha = !transposed_rights_alpha.empty();

    const auto two_termed_beta = !transposed_rights_beta.empty();

    if (two_termed_alpha)
    {
        errors::assertMsgCritical(transposed_rights_alpha.size() == rights_alpha.size(),
                                  std::string("SimdRIJKResponseDriver: The two terms of the alpha densities are not "
                                              "of one count"));

        _check_factors(basis, left_alpha, transposed_rights_alpha);
    }

    if (two_termed_beta)
    {
        errors::assertMsgCritical(transposed_rights_beta.size() == rights_beta.size(),
                                  std::string("SimdRIJKResponseDriver: The two terms of the beta densities are not of "
                                              "one count"));

        _check_factors(basis, left_beta, transposed_rights_beta);
    }

    const auto nao = basis.dimensions_of_basis();

    const auto ndens = rights_alpha.size();

    // NOTE: the two spins do not share a transformation. Their left factors are
    // the occupied orbitals of each of them, which differ both in what they are
    // and in how many there are, so the exchange of each spin is formed from its
    // own. That is what an unrestricted batch costs over a restricted one, and it
    // is the operators and not the spins which are swept together.

    auto exchange_of = [&](const CPackedMatrix              &left,
                           const std::vector<CPackedMatrix> &rights,
                           const std::vector<CPackedMatrix> &transposed_rights,
                           const bool                        two_termed) {
        auto wanted = rights;

        if (two_termed)
        {
            wanted.insert(wanted.end(), transposed_rights.begin(), transposed_rights.end());
        }

        return _exchange(operators, basis, aux_basis, left, wanted);
    };

    auto any_exchange = false;

    for (const auto &op : operators)
    {
        if (op.second != 0.0) any_exchange = true;
    }

    auto exchanges_alpha = std::vector<CPackedMatrix>();

    auto exchanges_beta = std::vector<CPackedMatrix>();

    if (any_exchange)
    {
        exchanges_alpha = exchange_of(left_alpha, rights_alpha, transposed_rights_alpha, two_termed_alpha);

        exchanges_beta = exchange_of(left_beta, rights_beta, transposed_rights_beta, two_termed_beta);
    }

    // NOTE: the right factor the Coulomb of a spin is taken with, which for a
    // density of two terms is the sum of theirs: the two terms have the same
    // symmetric part as the single term whose right factor is that sum, and the
    // Coulomb sees nothing else of a density.

    auto summed_of = [&](const std::vector<CPackedMatrix> &rights,
                         const std::vector<CPackedMatrix> &transposed_rights,
                         const bool                        two_termed) {
        auto summed = std::vector<CPackedMatrix>();

        if (!two_termed) return summed;

        const auto nvec = rights.front().number_of_columns();

        for (size_t idens = 0; idens < rights.size(); idens++)
        {
            auto total = CPackedMatrix(nao, nvec, mat_t::general);

            const auto *first = rights[idens].data();

            const auto *second = transposed_rights[idens].data();

            auto *values = total.data();

            for (size_t at = 0; at < nao * nvec; at++) values[at] = first[at] + second[at];

            summed.push_back(std::move(total));
        }

        return summed;
    };

    const auto summed_alpha = summed_of(rights_alpha, transposed_rights_alpha, two_termed_alpha);

    const auto summed_beta = summed_of(rights_beta, transposed_rights_beta, two_termed_beta);

    const auto &coulomb_rights_alpha = two_termed_alpha ? summed_alpha : rights_alpha;

    const auto &coulomb_rights_beta = two_termed_beta ? summed_beta : rights_beta;

    const auto nvec_alpha = left_alpha.number_of_columns();

    const auto nvec_beta = left_beta.number_of_columns();

    auto expanded = std::vector<double>(nao * nao, 0.0);

    for (size_t idens = 0; idens < ndens; idens++)
    {
        // NOTE: the Coulomb is of the density of both spins added, and is formed
        // once for the pair rather than once for each of them. Both spins see the
        // same Coulomb matrix, which is what makes an unrestricted build cheaper
        // than two restricted ones.

        auto density = CPackedMatrix(nao, nao, mat_t::general);

        density.zero();

        tdenblas::threadedMultABt(nao, nao, nvec_alpha, 1.0, left_alpha.data(), nvec_alpha,
                                  coulomb_rights_alpha[idens].data(), nvec_alpha, 1.0, density.data(), nao);

        tdenblas::threadedMultABt(nao, nao, nvec_beta, 1.0, left_beta.data(), nvec_beta,
                                  coulomb_rights_beta[idens].data(), nvec_beta, 1.0, density.data(), nao);

        const auto yvector = _drv.compute_y_vector(coulomb_vectors, basis, aux_basis, density);

        const auto coulomb = _drv.compute_fock_matrix(coulomb_vectors, basis, aux_basis, yvector);

        coulomb.to_dense(expanded.data());

        // NOTE: the Coulomb enters **once** and is not doubled, where the
        // restricted entry above is handed one spin's density and doubles it.

        auto assemble = [&](std::vector<CPackedMatrix>       &focks,
                            const std::vector<CPackedMatrix> &exchanges,
                            const bool                        two_termed) {
            focks.push_back(CPackedMatrix(nao, nao, mat_t::general));

            auto *values = focks.back().data();

            const auto *kvalues = exchanges.empty() ? nullptr : exchanges[idens].data();

            const auto *tvalues =
                (two_termed && (kvalues != nullptr)) ? exchanges[ndens + idens].data() : nullptr;

            for (size_t row = 0; row < nao; row++)
            {
                for (size_t col = 0; col < nao; col++)
                {
                    const auto at = row * nao + col;

                    values[at] = expanded[at];

                    if (kvalues != nullptr) values[at] -= kvalues[at];

                    if (tvalues != nullptr) values[at] -= tvalues[col * nao + row];
                }
            }
        };

        assemble(focks_alpha, exchanges_alpha, two_termed_alpha);

        assemble(focks_beta, exchanges_beta, two_termed_beta);
    }

    return {std::move(focks_alpha), std::move(focks_beta)};
}

auto
CSimdRIJKResponseDriver::_exchange(const std::vector<std::pair<const CSparseTensor *, double>> &operators,
                                   const CMolecularBasis                                      &basis,
                                   const CMolecularBasis                                      &aux_basis,
                                   const CPackedMatrix                                        &left,
                                   const std::vector<CPackedMatrix> &rights) const -> std::vector<CPackedMatrix>
{
    auto exchanges = std::vector<CPackedMatrix>();

    if (rights.empty()) return exchanges;

    _check_factors(basis, left, rights);

    const auto nao = basis.dimensions_of_basis();

    const auto nvec = left.number_of_columns();

    const auto ndens = rights.size();

    const auto naux = aux_basis.dimensions_of_basis();

    for (size_t idens = 0; idens < ndens; idens++)
    {
        exchanges.push_back(CPackedMatrix(nao, nao, mat_t::general));

        exchanges.back().zero();
    }

    // NOTE: the right factors side by side as one matrix of the basis by the rank
    // times the densities. The transformation reads the B vectors, which is the
    // dearest thing it does, and one wide matrix reads them once where a matrix
    // for each density reads them once for each density.

    const auto nwide = nvec * ndens;

    auto stacked = CPackedMatrix(nao, nwide, mat_t::general);

    stacked.zero();

    for (size_t idens = 0; idens < ndens; idens++)
    {
        const auto *values = rights[idens].data();

        auto *target = stacked.data();

        for (size_t mu = 0; mu < nao; mu++)
        {
            for (size_t k = 0; k < nvec; k++)
            {
                target[mu * nwide + idens * nvec + k] = values[mu * nvec + k];
            }
        }
    }

    // NOTE: the auxiliary basis in batches, so the transformed vectors are bounded
    // by the budget and not by the problem. One batch holds the left factor
    // transformed and the right factors transformed, which is the basis by the
    // rank and the basis by the rank times the densities.

    const auto per_function = nao * nvec * (1 + ndens) * sizeof(double);

    const auto by_memory = _budget / std::max(per_function, size_t{1});

    const auto nbatch = std::min(naux, std::max(_min_batch, by_memory));

    // NOTE: one operator or two, inside one pass over the auxiliary basis and into
    // one set of exchange matrices. Two calls of this would form two sets and add
    // them, which for a batch of twenty trial vectors of five hundred functions is
    // eighty megabytes allocated to be summed and thrown away. Each operator has
    // its own B vectors and so its own transformation of the factors -- that part
    // cannot be shared -- but the batching, the stacking and the output are.

    for (size_t first = 0; first < naux; first += nbatch)
    {
        const auto last = std::min(first + nbatch, naux);

        const auto count = static_cast<int>(last - first);

        for (const auto &op : operators)
        {
            // NOTE: named rather than taken by a structured binding, which a
            // parallel region below cannot capture.
            const auto *tensor = op.first;

            const auto scale = op.second;

            if (scale == 0.0) continue;

            const auto uvecs = _drv.compute_w_vectors(*tensor, basis, aux_basis, left, first, last);

            const auto pvecs = _drv.compute_w_vectors(*tensor, basis, aux_basis, stacked, first, last);

        // NOTE: the densities and not the auxiliary functions are what the threads
        // divide, so that each of them writes into an exchange matrix of its own
        // and no two of them accumulate into the same one. A batch of trial vectors
        // is what a response calculation always has.

            const auto nrange = static_cast<int>(ndens);

#pragma omp parallel for schedule(static) if (nrange > 1)
            for (int at = 0; at < nrange; at++)
            {
                const auto idens = static_cast<size_t>(at);

                for (int q = 0; q < count; q++)
                {
                    const auto iq = static_cast<size_t>(q);

                    tdenblas::threadedMultABt(nao, nao, nvec, scale, uvecs[iq].data(), nvec,
                                              pvecs[iq].data() + idens * nvec, nwide, 1.0,
                                              exchanges[idens].data(), nao);
                }
            }
        }
    }

    return exchanges;
}

auto
CSimdRIJKResponseDriver::compute(const CSparseTensor              &bq_vectors,
                                 const CMolecularBasis            &basis,
                                 const CMolecularBasis            &aux_basis,
                                 const CPackedMatrix              &left,
                                 const std::vector<CPackedMatrix> &rights,
                                 const double                      exchange_scaling_factor,
                                 const CSparseTensor              *bq_vectors_erf,
                                 const double                      erf_exchange_scaling_factor) const
    -> std::vector<CPackedMatrix>
{
    return compute(bq_vectors,
                   basis,
                   aux_basis,
                   left,
                   rights,
                   {},
                   exchange_scaling_factor,
                   bq_vectors_erf,
                   erf_exchange_scaling_factor);
}

auto
CSimdRIJKResponseDriver::compute(const CSparseTensor              &bq_vectors,
                                 const CMolecularBasis            &basis,
                                 const CMolecularBasis            &aux_basis,
                                 const CPackedMatrix              &left,
                                 const std::vector<CPackedMatrix> &rights,
                                 const std::vector<CPackedMatrix> &transposed_rights,
                                 const double                      exchange_scaling_factor,
                                 const CSparseTensor              *bq_vectors_erf,
                                 const double                      erf_exchange_scaling_factor) const
    -> std::vector<CPackedMatrix>
{
    // NOTE: a caller which asks for the exchange of the attenuated operator has to
    // hand over the B vectors of it. Leaving the long-range term out in silence is
    // what this refuses; a response calculation which did would converge to the
    // wrong excitation energies with nothing to say so.
    errors::assertMsgCritical((erf_exchange_scaling_factor == 0.0) || (bq_vectors_erf != nullptr),
                              std::string("SimdRIJKResponseDriver: The exchange of the attenuated operator was asked "
                                          "for and no attenuated B vectors were given"));

    auto focks = std::vector<CPackedMatrix>();

    if (rights.empty()) return focks;

    _check_factors(basis, left, rights);

    const auto two_termed = !transposed_rights.empty();

    if (two_termed)
    {
        errors::assertMsgCritical(
            transposed_rights.size() == rights.size(),
            std::string("SimdRIJKResponseDriver: The two terms of the densities are not of one count"));

        _check_factors(basis, left, transposed_rights);
    }

    const auto nao = basis.dimensions_of_basis();

    const auto nvec = left.number_of_columns();

    const auto ndens = rights.size();

    // NOTE: the exchange of the whole batch first, so the left factor is
    // transformed once for all of it. A pure functional asks for none of it.

    // NOTE: what comes back here is already scaled and, for a hybrid range
    // separated functional, already the sum of the two operators. The assembly
    // below therefore subtracts it as it stands. Scaling at the end instead would
    // need the two operators kept apart all the way down, which is a second set of
    // matrices of the basis squared for every density of the batch.

    auto exchanges = std::vector<CPackedMatrix>();

    const auto attenuated = (erf_exchange_scaling_factor != 0.0);

    if ((exchange_scaling_factor != 0.0) || attenuated)
    {
        // NOTE: the two terms side by side in one batch, so the shared factor is
        // transformed once for both of them and not once for each.

        auto wanted = rights;

        if (two_termed)
        {
            wanted.insert(wanted.end(), transposed_rights.begin(), transposed_rights.end());
        }

        exchanges = attenuated ? compute_exchange(bq_vectors,
                                                  *bq_vectors_erf,
                                                  basis,
                                                  aux_basis,
                                                  left,
                                                  wanted,
                                                  exchange_scaling_factor,
                                                  erf_exchange_scaling_factor)
                               : _exchange({{&bq_vectors, exchange_scaling_factor}},
                                           basis,
                                           aux_basis,
                                           left,
                                           wanted);
    }

    // NOTE: the right factor the Coulomb is taken with. The two terms of a
    // density have the same symmetric part as the single term whose right factor
    // is the sum of theirs, and the Coulomb sees nothing else of a density, so it
    // is one closure of the B vectors rather than two.

    auto summed = std::vector<CPackedMatrix>();

    if (two_termed)
    {
        for (size_t idens = 0; idens < ndens; idens++)
        {
            auto total = CPackedMatrix(nao, nvec, mat_t::general);

            const auto *first = rights[idens].data();

            const auto *second = transposed_rights[idens].data();

            auto *values = total.data();

            for (size_t at = 0; at < nao * nvec; at++) values[at] = first[at] + second[at];

            summed.push_back(std::move(total));
        }
    }

    const auto &coulomb_rights = two_termed ? summed : rights;

    auto expanded = std::vector<double>(nao * nao, 0.0);

    for (size_t idens = 0; idens < ndens; idens++)
    {
        // the density this pair of factors stands for, which the Coulomb closes
        // whole. compute_y_vector reads the type and contracts a general density
        // as a general one.

        auto density = CPackedMatrix(nao, nao, mat_t::general);

        density.zero();

        tdenblas::threadedMultABt(nao, nao, nvec, 1.0, left.data(), nvec, coulomb_rights[idens].data(),
                                  nvec, 1.0, density.data(), nao);

        const auto yvector = _drv.compute_y_vector(bq_vectors, basis, aux_basis, density);

        const auto coulomb = _drv.compute_fock_matrix(bq_vectors, basis, aux_basis, yvector);

        // NOTE: the Coulomb comes back symmetric and packed and the exchange is
        // general, so the sum is taken over the expanded Coulomb. The matrix which
        // is returned is general: the density was not symmetric and neither is
        // what it gives.

        coulomb.to_dense(expanded.data());

        focks.push_back(CPackedMatrix(nao, nao, mat_t::general));

        auto *values = focks.back().data();

        const auto *kvalues = exchanges.empty() ? nullptr : exchanges[idens].data();

        // NOTE: the exchange of the second term is the exchange of its transpose,
        // transposed, which is what the batch above holds after the first count of
        // them. It is read the other way round as it is subtracted.

        const auto *tvalues = (two_termed && (kvalues != nullptr))
                                  ? exchanges[ndens + idens].data()
                                  : nullptr;

        for (size_t row = 0; row < nao; row++)
        {
            for (size_t col = 0; col < nao; col++)
            {
                const auto at = row * nao + col;

                values[at] = 2.0 * expanded[at];

                if (kvalues != nullptr) values[at] -= kvalues[at];

                if (tvalues != nullptr)
                {
                    values[at] -= tvalues[col * nao + row];
                }
            }
        }
    }

    return focks;
}
