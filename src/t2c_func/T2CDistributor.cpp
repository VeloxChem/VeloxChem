//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "T2CDistributor.hpp"

namespace t2cfunc {  // t2cfunc namespace

auto
distribute(CSubMatrix*                 matrix,
           const TDoubleArray&         buffer,
           const std::vector<int64_t>& indexes,
           const int64_t               bra_comp,
           const int64_t               ket_comp,
           const int64_t               bra_igto,
           const int64_t               ket_first,
           const int64_t               ket_last) -> void
{
    const auto ncgtos = indexes[0];

    const auto bra_idx = ncgtos * bra_comp + indexes[bra_igto + 1];

    const auto bra_off = matrix->getOffsetOfRows();

    const auto ket_off = matrix->getOffsetOfColumns();

    const auto ket_dim = matrix->getNumberOfColumns();

    auto data = matrix->getData();

    const auto ij_off = (bra_idx - bra_off) * ket_dim - ket_off;

    for (int64_t i = ket_first; i < ket_last; i++)
    {
        const auto ket_idx = ncgtos * ket_comp + indexes[i + 1];

        data[ij_off + ket_idx] += buffer[i - ket_first];
    }
}

auto
distribute(CSubMatrix*                 matrix,
           const TDoubleArray&         buffer,
           const double                factor,
           const std::vector<int64_t>& indexes,
           const int64_t               bra_comp,
           const int64_t               ket_comp,
           const int64_t               bra_igto,
           const int64_t               ket_first,
           const int64_t               ket_last) -> void
{
    const auto ncgtos = indexes[0];

    const auto bra_idx = ncgtos * bra_comp + indexes[bra_igto + 1];

    const auto bra_off = matrix->getOffsetOfRows();

    const auto ket_off = matrix->getOffsetOfColumns();

    const auto ket_dim = matrix->getNumberOfColumns();

    auto data = matrix->getData();

    const auto ij_off = (bra_idx - bra_off) * ket_dim - ket_off;

    for (int64_t i = ket_first; i < ket_last; i++)
    {
        const auto ket_idx = ncgtos * ket_comp + indexes[i + 1];

        data[ij_off + ket_idx] += factor * buffer[i - ket_first];
    }
}

auto
distribute(CSubMatrix*                 matrix,
           const TDoubleArray&         buffer,
           const std::vector<int64_t>& bra_indexes,
           const std::vector<int64_t>& ket_indexes,
           const int64_t               bra_comp,
           const int64_t               ket_comp,
           const int64_t               bra_igto,
           const int64_t               ket_first,
           const int64_t               ket_last,
           const mat_t                 mat_type) -> void
{
    const auto bra_ncgtos = bra_indexes[0];

    const auto ket_ncgtos = ket_indexes[0];

    const auto bra_idx = bra_ncgtos * bra_comp + bra_indexes[bra_igto + 1];

    const auto bra_off = matrix->getOffsetOfRows();

    const auto ket_off = matrix->getOffsetOfColumns();

    const auto ket_dim = matrix->getNumberOfColumns();

    auto data = matrix->getData();

    const auto ij_off = (bra_idx - bra_off) * ket_dim - ket_off;

    const auto ji_off = bra_idx - bra_off * ket_dim - ket_off;

    for (int64_t i = ket_first; i < ket_last; i++)
    {
        const auto ket_idx = ket_ncgtos * ket_comp + ket_indexes[i + 1];

        const auto fval = buffer[i - ket_first];

        data[ij_off + ket_idx] += fval;

        if (mat_type == mat_t::symm)
        {
            data[ji_off + ket_idx * ket_dim] += fval;
        }

        if (mat_type == mat_t::antisymm)
        {
            data[ji_off + ket_idx * ket_dim] += -fval;
        }
    }
}

auto
distribute(CSubMatrix*                 matrix,
           const TDoubleArray&         buffer,
           const double                factor,
           const std::vector<int64_t>& bra_indexes,
           const std::vector<int64_t>& ket_indexes,
           const int64_t               bra_comp,
           const int64_t               ket_comp,
           const int64_t               bra_igto,
           const int64_t               ket_first,
           const int64_t               ket_last,
           const mat_t                 mat_type) -> void
{
    const auto bra_ncgtos = bra_indexes[0];

    const auto ket_ncgtos = ket_indexes[0];

    const auto bra_idx = bra_ncgtos * bra_comp + bra_indexes[bra_igto + 1];

    const auto bra_off = matrix->getOffsetOfRows();

    const auto ket_off = matrix->getOffsetOfColumns();

    const auto ket_dim = matrix->getNumberOfColumns();

    auto data = matrix->getData();

    const auto ij_off = (bra_idx - bra_off) * ket_dim - ket_off;

    const auto ji_off = bra_idx - bra_off * ket_dim - ket_off;

    for (int64_t i = ket_first; i < ket_last; i++)
    {
        const auto ket_idx = ket_ncgtos * ket_comp + ket_indexes[i + 1];

        const auto fval = factor * buffer[i - ket_first];

        data[ij_off + ket_idx] += fval;

        if (mat_type == mat_t::symm)
        {
            data[ji_off + ket_idx * ket_dim] += fval;
        }

        if (mat_type == mat_t::antisymm)
        {
            data[ji_off + ket_idx * ket_dim] += -fval;
        }
    }
}

auto
distribute(CSubMatrix*                 matrix,
           const TDoubleArray&         buffer,
           const std::vector<int64_t>& bra_indexes,
           const std::vector<int64_t>& ket_indexes,
           const int64_t               bra_comp,
           const int64_t               ket_comp,
           const int64_t               bra_igto,
           const int64_t               ket_first,
           const int64_t               ket_last,
           const bool                  ang_order) -> void
{
    const auto bra_ncgtos = bra_indexes[0];

    const auto ket_ncgtos = ket_indexes[0];

    const auto bra_idx = bra_ncgtos * bra_comp + bra_indexes[bra_igto + 1];

    const auto bra_off = matrix->getOffsetOfRows();

    const auto ket_off = matrix->getOffsetOfColumns();

    const auto ket_dim = matrix->getNumberOfColumns();

    auto data = matrix->getData();

    const auto ij_off = (bra_idx - bra_off) * ket_dim - ket_off;

    const auto ji_off = bra_idx - bra_off * ket_dim - ket_off;

    for (int64_t i = ket_first; i < ket_last; i++)
    {
        const auto ket_idx = ket_ncgtos * ket_comp + ket_indexes[i + 1];

        if (ang_order)
        {
            data[ij_off + ket_idx] += buffer[i - ket_first];
        }
        else
        {
            data[ji_off + ket_idx * ket_dim] += buffer[i - ket_first];
        }
    }
}

auto
distribute(CSubMatrix*                 matrix,
           const TDoubleArray&         buffer,
           const double                factor,
           const std::vector<int64_t>& bra_indexes,
           const std::vector<int64_t>& ket_indexes,
           const int64_t               bra_comp,
           const int64_t               ket_comp,
           const int64_t               bra_igto,
           const int64_t               ket_first,
           const int64_t               ket_last,
           const bool                  ang_order) -> void
{
    const auto bra_ncgtos = bra_indexes[0];

    const auto ket_ncgtos = ket_indexes[0];

    const auto bra_idx = bra_ncgtos * bra_comp + bra_indexes[bra_igto + 1];

    const auto bra_off = matrix->getOffsetOfRows();

    const auto ket_off = matrix->getOffsetOfColumns();

    const auto ket_dim = matrix->getNumberOfColumns();

    auto data = matrix->getData();

    const auto ij_off = (bra_idx - bra_off) * ket_dim - ket_off;

    const auto ji_off = bra_idx - bra_off * ket_dim - ket_off;

    for (int64_t i = ket_first; i < ket_last; i++)
    {
        const auto ket_idx = ket_ncgtos * ket_comp + ket_indexes[i + 1];

        if (ang_order)
        {
            data[ij_off + ket_idx] += factor * buffer[i - ket_first];
        }
        else
        {
            data[ji_off + ket_idx * ket_dim] += factor * buffer[i - ket_first];
        }
    }
}

}  // namespace t2cfunc
