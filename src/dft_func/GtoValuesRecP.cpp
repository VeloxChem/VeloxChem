#include "GtoValuesRecP.hpp"

#include <cmath>

#include "DftFunc.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"

namespace gtoval {  // gtoval namespace

auto
getLdaValuesRecP(const CGtoBlock&            gto_block,
                 const std::vector<double>&  grid_coords_x,
                 const std::vector<double>&  grid_coords_y,
                 const std::vector<double>&  grid_coords_z,
                 const std::vector<int64_t>& gtos_mask) -> CMatrix
{
    // set up GTO values storage

    const auto nrows = mathfunc::countSignificantElements(gtos_mask);

    const auto ncols = static_cast<int64_t>(grid_coords_x.size());

    auto gto_values = matfunc::makeMatrix("LDA", 3 * nrows, ncols);

    auto submat = gto_values.getSubMatrix({0, 0});

    submat->zero();

    // set up GTOs data

    const auto gto_exps = gto_block.getExponents();

    const auto gto_norms = gto_block.getNormalizationFactors();

    const auto gto_coords = gto_block.getCoordinates();

    // set up grid data

    auto g_x = grid_coords_x.data();

    auto g_y = grid_coords_y.data();

    auto g_z = grid_coords_z.data();

    // set GTOs block dimensions

    const auto ncgtos = gto_block.getNumberOfBasisFunctions();

    const auto npgtos = gto_block.getNumberOfPrimitives();

    // compute GTO values for S type GTOs on grid

    std::vector<double> buffer_x(ncols);

    std::vector<double> buffer_y(ncols);

    std::vector<double> buffer_z(ncols);

    auto ptr_buffer_x = buffer_x.data();

    auto ptr_buffer_y = buffer_y.data();

    auto ptr_buffer_z = buffer_z.data();

    int64_t irow = 0;

    for (int64_t i = 0; i < ncgtos; i++)
    {
        if (gtos_mask[i] == 1)
        {
            // set up GTO coordinates

            const auto r_x = gto_coords[i][0];

            const auto r_y = gto_coords[i][1];

            const auto r_z = gto_coords[i][2];

            // compute GTO values on grid

            mathfunc::zero(buffer_x);

            mathfunc::zero(buffer_y);

            mathfunc::zero(buffer_z);

            for (int64_t j = 0; j < npgtos; j++)
            {
                const auto fexp = gto_exps[j * ncgtos + i];

                const auto fnorm = gto_norms[j * ncgtos + i];

#pragma omp simd
                for (int64_t k = 0; k < ncols; k++)
                {
                    const auto gr_x = g_x[k] - r_x;

                    const auto gr_y = g_y[k] - r_y;

                    const auto gr_z = g_z[k] - r_z;

                    const auto fss = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

                    ptr_buffer_x[k] += gr_x * fss;

                    ptr_buffer_y[k] += gr_y * fss;

                    ptr_buffer_z[k] += gr_z * fss;
                }
            }

            // distribute GTO values into submatrix

            gtoval::distribute(submat, buffer_x, 2 * nrows + irow);

            gtoval::distribute(submat, buffer_y, irow);

            gtoval::distribute(submat, buffer_z, nrows + irow);

            irow++;
        }
    }

    return gto_values;
}

auto
getGgaValuesRecP(const CGtoBlock&            gto_block,
                 const std::vector<double>&  grid_coords_x,
                 const std::vector<double>&  grid_coords_y,
                 const std::vector<double>&  grid_coords_z,
                 const std::vector<int64_t>& gtos_mask) -> CMatrix
{
    // set up GTO values storage

    const auto nrows = mathfunc::countSignificantElements(gtos_mask);

    const auto ncols = static_cast<int64_t>(grid_coords_x.size());

    auto gto_values = matfunc::makeMatrix("GGA", 3 * nrows, ncols);

    auto submat_0 = gto_values.getSubMatrix({0, 0});
    auto submat_x = gto_values.getSubMatrix({1, 0});
    auto submat_y = gto_values.getSubMatrix({1, 1});
    auto submat_z = gto_values.getSubMatrix({1, 2});

    submat_0->zero();
    submat_x->zero();
    submat_y->zero();
    submat_z->zero();

    // set up GTOs data

    const auto gto_exps = gto_block.getExponents();

    const auto gto_norms = gto_block.getNormalizationFactors();

    const auto gto_coords = gto_block.getCoordinates();

    // set up grid data

    auto g_x = grid_coords_x.data();

    auto g_y = grid_coords_y.data();

    auto g_z = grid_coords_z.data();

    // set GTOs block dimensions

    const auto ncgtos = gto_block.getNumberOfBasisFunctions();

    const auto npgtos = gto_block.getNumberOfPrimitives();

    // compute GTO values for S type GTOs on grid

    std::vector<double> buffer_0_x(ncols);
    std::vector<double> buffer_0_y(ncols);
    std::vector<double> buffer_0_z(ncols);

    std::vector<double> buffer_x_x(ncols);
    std::vector<double> buffer_x_y(ncols);
    std::vector<double> buffer_x_z(ncols);

    std::vector<double> buffer_y_x(ncols);
    std::vector<double> buffer_y_y(ncols);
    std::vector<double> buffer_y_z(ncols);

    std::vector<double> buffer_z_x(ncols);
    std::vector<double> buffer_z_y(ncols);
    std::vector<double> buffer_z_z(ncols);

    auto ptr_buffer_0_x = buffer_0_x.data();
    auto ptr_buffer_0_y = buffer_0_y.data();
    auto ptr_buffer_0_z = buffer_0_z.data();

    auto ptr_buffer_x_x = buffer_x_x.data();
    auto ptr_buffer_x_y = buffer_x_y.data();
    auto ptr_buffer_x_z = buffer_x_z.data();

    auto ptr_buffer_y_x = buffer_y_x.data();
    auto ptr_buffer_y_y = buffer_y_y.data();
    auto ptr_buffer_y_z = buffer_y_z.data();

    auto ptr_buffer_z_x = buffer_z_x.data();
    auto ptr_buffer_z_y = buffer_z_y.data();
    auto ptr_buffer_z_z = buffer_z_z.data();

    int64_t irow = 0;

    for (int64_t i = 0; i < ncgtos; i++)
    {
        if (gtos_mask[i] == 1)
        {
            // set up GTO coordinates

            const auto r_x = gto_coords[i][0];

            const auto r_y = gto_coords[i][1];

            const auto r_z = gto_coords[i][2];

            // compute GTO values on grid

            mathfunc::zero(buffer_0_x);
            mathfunc::zero(buffer_0_y);
            mathfunc::zero(buffer_0_z);

            mathfunc::zero(buffer_x_x);
            mathfunc::zero(buffer_x_y);
            mathfunc::zero(buffer_x_z);

            mathfunc::zero(buffer_y_x);
            mathfunc::zero(buffer_y_y);
            mathfunc::zero(buffer_y_z);

            mathfunc::zero(buffer_z_x);
            mathfunc::zero(buffer_z_y);
            mathfunc::zero(buffer_z_z);

            for (int64_t j = 0; j < npgtos; j++)
            {
                const auto fexp = gto_exps[j * ncgtos + i];

                const auto fnorm = gto_norms[j * ncgtos + i];

#pragma omp simd
                for (int64_t k = 0; k < ncols; k++)
                {
                    const auto gr_x = g_x[k] - r_x;

                    const auto gr_y = g_y[k] - r_y;

                    const auto gr_z = g_z[k] - r_z;

                    const auto f00 = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

                    const auto fg0 = -2.0 * fexp;

                    ptr_buffer_0_x[k] += f00 * gr_x;
                    ptr_buffer_x_x[k] += f00 * (1.0 + gr_x * gr_x * fg0);
                    ptr_buffer_y_x[k] += f00 * gr_x * gr_y * fg0;
                    ptr_buffer_z_x[k] += f00 * gr_x * gr_z * fg0;

                    ptr_buffer_0_y[k] += f00 * gr_y;
                    ptr_buffer_x_y[k] += f00 * gr_y * gr_x * fg0;
                    ptr_buffer_y_y[k] += f00 * (1.0 + gr_y * gr_y * fg0);
                    ptr_buffer_z_y[k] += f00 * gr_y * gr_z * fg0;

                    ptr_buffer_0_z[k] += f00 * gr_z;
                    ptr_buffer_x_z[k] += f00 * gr_z * gr_x * fg0;
                    ptr_buffer_y_z[k] += f00 * gr_z * gr_y * fg0;
                    ptr_buffer_z_z[k] += f00 * (1.0 + gr_z * gr_z * fg0);
                }
            }

            // distribute GTO values into submatrix

            gtoval::distribute(submat_0, buffer_0_x, 2 * nrows + irow);
            gtoval::distribute(submat_0, buffer_0_y, irow);
            gtoval::distribute(submat_0, buffer_0_z, nrows + irow);

            gtoval::distribute(submat_x, buffer_x_x, 2 * nrows + irow);
            gtoval::distribute(submat_x, buffer_x_y, irow);
            gtoval::distribute(submat_x, buffer_x_z, nrows + irow);

            gtoval::distribute(submat_y, buffer_y_x, 2 * nrows + irow);
            gtoval::distribute(submat_y, buffer_y_y, irow);
            gtoval::distribute(submat_y, buffer_y_z, nrows + irow);

            gtoval::distribute(submat_z, buffer_z_x, 2 * nrows + irow);
            gtoval::distribute(submat_z, buffer_z_y, irow);
            gtoval::distribute(submat_z, buffer_z_z, nrows + irow);

            irow++;
        }
    }

    return gto_values;
}

}  // namespace gtoval
