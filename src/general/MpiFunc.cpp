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

#include "MpiFunc.hpp"

#include <mpi.h>

#include <iostream>
#include <sstream>

#include "MathFunc.hpp"

namespace mpi {  // mpi namespace

auto
init(int argc, char** argv) -> bool
{
    if constexpr (ENABLE_MPI)
    {
        if (initialized())
        {
            return true;
        }

        int32_t mlevel_int32 = 0;

        auto merror = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mlevel_int32);

        if (merror != MPI_SUCCESS)
        {
            mpi::abort(merror, "mpi::init()");

            return false;
        }
    }

    return true;
}

auto
initialized() -> bool
{
    if constexpr (ENABLE_MPI)
    {
        int32_t minit_int32 = 0;

        MPI_Initialized(&minit_int32);

        if (minit_int32 == 1) return true;
    }

    return false;
}

auto
finalize() -> bool
{
    if constexpr (ENABLE_MPI)
    {
        auto merror = MPI_Finalize();

        if (merror != MPI_SUCCESS)
        {
            mpi::abort(merror, "mpi::finalize()");

            return false;
        }
    }

    return true;
}

auto
abort(const int errorcode, const char* label) -> void
{
    if constexpr (ENABLE_MPI)
    {
        int32_t errclass_int32 = 0;

        MPI_Error_class(errorcode, &errclass_int32);

        int32_t errlen_int32 = 0;

        char errstr[MPI_MAX_ERROR_STRING];

        MPI_Error_string(errorcode, errstr, &errlen_int32);

        std::stringstream sst;

        sst << "**** Critical Error in " << label << " ****" << std::endl;

        sst << "MPI ERROR " << errclass_int32 << ": " << errstr << std::endl;

        std::cerr << sst.str();

        MPI_Abort(MPI_COMM_WORLD, errorcode);
    }
}

auto
abort(const int errorcode, const std::string& label) -> void
{
    mpi::abort(errorcode, label.c_str());
}

auto
rank(MPI_Comm comm) -> int64_t
{
    int32_t mrank_int32 = 0;

    if constexpr (ENABLE_MPI) MPI_Comm_rank(comm, &mrank_int32);

    return static_cast<int64_t>(mrank_int32);
}

auto
nodes(MPI_Comm comm) -> int64_t
{
    int32_t mnodes_int32 = 1;

    if constexpr (ENABLE_MPI) MPI_Comm_size(comm, &mnodes_int32);

    return static_cast<int64_t>(mnodes_int32);
}

auto
bcastScalar(const int64_t val, MPI_Comm comm) -> int64_t
{
    int64_t bcast_val = 0;

    if (mpi::rank(comm) == mpi::master()) bcast_val = val;

    auto merror = MPI_Bcast(&bcast_val, 1, MPI_INT64_T, static_cast<int32_t>(mpi::master()), comm);

    if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::bcastScalar");

    return bcast_val;
}

auto
bcastScalar(const double val, MPI_Comm comm) -> double
{
    double bcast_val = 0.0;

    if (mpi::rank(comm) == mpi::master()) bcast_val = val;

    auto merror = MPI_Bcast(&bcast_val, 1, MPI_DOUBLE, static_cast<int32_t>(mpi::master()), comm);

    if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::bcastScalar");

    return bcast_val;
}

auto
bcastStdVectorInt(const std::vector<int64_t>& vector, MPI_Comm comm) -> std::vector<int64_t>
{
    if constexpr (ENABLE_MPI)
    {
        // MPI info

        auto mpi_master = mpi::master();

        auto mpi_master_int32 = static_cast<int32_t>(mpi_master);

        auto rank = mpi::rank(comm);

        auto nodes = mpi::nodes(comm);

        if (nodes == 1) return vector;

        // broadcast size

        int64_t num_elems = 0;

        if (rank == mpi_master) num_elems = static_cast<int64_t>(vector.size());

        auto merror = MPI_Bcast(&num_elems, 1, MPI_INT64_T, mpi_master_int32, comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::bcastStdVectorInt");

        // broadcast data

        std::vector<int64_t> bcast_vector;

        if (rank == mpi_master)
        {
            bcast_vector = vector;
        }
        else
        {
            bcast_vector = std::vector<int64_t>(num_elems);
        }

        merror = MPI_Bcast(bcast_vector.data(), static_cast<int32_t>(num_elems), MPI_INT64_T, mpi_master_int32, comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::bcastStdVectorInt");

        return bcast_vector;
    }

    return vector;
}

auto
bcastDenseMatrix(const CDenseMatrix& matrix, MPI_Comm comm) -> CDenseMatrix
{
    if constexpr (ENABLE_MPI)
    {
        // MPI info

        auto mpi_master = mpi::master();

        auto mpi_master_int32 = static_cast<int32_t>(mpi_master);

        auto rank = mpi::rank(comm);

        auto nodes = mpi::nodes(comm);

        if (nodes == 1) return matrix;

        // broadcast size

        std::vector<int64_t> num_rows_cols(2);

        if (rank == mpi_master)
        {
            num_rows_cols[0] = matrix.getNumberOfRows();

            num_rows_cols[1] = matrix.getNumberOfColumns();
        }

        auto merror = MPI_Bcast(num_rows_cols.data(), 2, MPI_INT64_T, mpi_master_int32, comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::bcastDenseMatrix");

        // broadcast data

        auto nrows = num_rows_cols[0];

        auto ncols = num_rows_cols[1];

        CDenseMatrix bcast_matrix;

        if (rank == mpi_master)
        {
            bcast_matrix = matrix;
        }
        else
        {
            bcast_matrix = CDenseMatrix(nrows, ncols);
        }

        merror = MPI_Bcast(bcast_matrix.values(), static_cast<int32_t>(nrows * ncols), MPI_DOUBLE, mpi_master_int32, comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::bcastDenseMatrix");

        return bcast_matrix;
    }

    return matrix;
}

auto
scatterStdVector(const std::vector<int64_t>& vec, MPI_Comm comm) -> std::vector<int64_t>
{
    if constexpr (ENABLE_MPI)
    {
        // MPI info

        auto mpi_master = mpi::master();

        auto mpi_master_int32 = static_cast<int32_t>(mpi_master);

        auto rank = mpi::rank(comm);

        auto nodes = mpi::nodes(comm);

        if (nodes == 1) return vec;

        // vector size, counts and displacements

        int64_t nelems = 0;

        std::vector<int32_t> counts, displs;

        if (rank == mpi_master)
        {
            nelems = static_cast<int64_t>(vec.size());

            for (const auto count : mathfunc::batch_sizes(nelems, nodes))
            {
                counts.push_back(static_cast<int32_t>(count));
            }

            for (const auto displ : mathfunc::batch_offsets(nelems, nodes))
            {
                displs.push_back(static_cast<int32_t>(displ));
            }
        }

        auto merror = MPI_Bcast(&nelems, 1, MPI_INT64_T, mpi_master_int32, comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::scatterStdVector");

        // receiving vector

        auto recv_size = mathfunc::batch_size(nelems, rank, nodes);

        auto recv_size_int32 = static_cast<int32_t>(recv_size);

        std::vector<int64_t> recv_vec(recv_size);

        // scatter data

        merror = MPI_Scatterv(
            vec.data(), counts.data(), displs.data(), MPI_INT64_T, recv_vec.data(), recv_size_int32, MPI_INT64_T, mpi_master_int32, comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::scatterStdVector");

        return recv_vec;
    }

    return vec;
}

auto
gatherDenseMatricesByColumns(const CDenseMatrix& matrix, MPI_Comm comm) -> CDenseMatrix
{
    if constexpr (ENABLE_MPI)
    {
        // MPI info

        auto mpi_master = mpi::master();

        auto mpi_master_int32 = static_cast<int32_t>(mpi_master);

        auto rank = mpi::rank(comm);

        auto nodes = mpi::nodes(comm);

        if (nodes == 1) return matrix;

        // matrix size

        auto nrows = matrix.getNumberOfRows();

        auto ncols = matrix.getNumberOfColumns();

        // number of columns

        std::vector<int64_t> col_counts;

        if (rank == mpi_master) col_counts = std::vector<int64_t>(nodes);

        auto merror = MPI_Gather(&ncols, 1, MPI_INT64_T, col_counts.data(), 1, MPI_INT64_T, mpi_master_int32, comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::gatherDenseMatricesByColumns MPI_Gather");

        // counts and displacements of elements

        std::vector<int32_t> counts, displs;

        CDenseMatrix gathered_data;

        if (rank == mpi_master)
        {
            int64_t ncols_sum = 0;

            for (const auto nc : col_counts)
            {
                counts.push_back(static_cast<int32_t>(nc));

                displs.push_back(static_cast<int32_t>(ncols_sum));

                ncols_sum += nc;
            }

            gathered_data = CDenseMatrix(nrows, ncols_sum);
        }

        // gather data

        for (int64_t row = 0; row < nrows; row++)
        {
            merror = MPI_Gatherv(matrix.row(row),
                                 static_cast<int32_t>(ncols),
                                 MPI_DOUBLE,
                                 gathered_data.row(row),
                                 counts.data(),
                                 displs.data(),
                                 MPI_DOUBLE,
                                 mpi_master_int32,
                                 comm);

            if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::gatherDenseMatricesByColumns MPI_Gatherv");
        }

        return gathered_data;
    }

    return matrix;
}

}  // namespace mpi
