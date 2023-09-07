#include "GtoValuesRecD.hpp"

#include <cmath>

#include "DftFunc.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"

namespace gtoval {  // gtoval namespace

auto
getLdaValuesRecD(const CGtoBlock&            gto_block,
                 const int64_t               n_points,
                 const double*               grid_coords_x,
                 const double*               grid_coords_y,
                 const double*               grid_coords_z,
                 const std::vector<int64_t>& gtos_mask) -> CMatrix
{
    // spherical transformation factors

    const double f2_3 = 2.0 * std::sqrt(3.0);

    // set up GTO values storage

    const auto nrows = mathfunc::countSignificantElements(gtos_mask);

    const auto ncols = n_points;

    auto gto_values = matfunc::makeMatrix("LDA", 5 * nrows, ncols);

    auto submat = gto_values.getSubMatrix({0, 0});

    // set up GTOs data

    const auto gto_exps = gto_block.getExponents();

    const auto gto_norms = gto_block.getNormalizationFactors();

    const auto gto_coords = gto_block.getCoordinates();

    // set GTOs block dimensions

    const auto ncgtos = gto_block.getNumberOfBasisFunctions();

    const auto npgtos = gto_block.getNumberOfPrimitives();

    // compute GTO values for S type GTOs on grid

    std::vector<double> buffer_xx(ncols);

    std::vector<double> buffer_xy(ncols);

    std::vector<double> buffer_xz(ncols);

    std::vector<double> buffer_yy(ncols);

    std::vector<double> buffer_yz(ncols);

    std::vector<double> buffer_zz(ncols);

    auto ptr_buffer_xx = buffer_xx.data();

    auto ptr_buffer_xy = buffer_xy.data();

    auto ptr_buffer_xz = buffer_xz.data();

    auto ptr_buffer_yy = buffer_yy.data();

    auto ptr_buffer_yz = buffer_yz.data();

    auto ptr_buffer_zz = buffer_zz.data();

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

            mathfunc::zero(buffer_xx);

            mathfunc::zero(buffer_xy);

            mathfunc::zero(buffer_xz);

            mathfunc::zero(buffer_yy);

            mathfunc::zero(buffer_yz);

            mathfunc::zero(buffer_zz);

            for (int64_t j = 0; j < npgtos; j++)
            {
                const auto fexp = gto_exps[j * ncgtos + i];

                const auto fnorm = gto_norms[j * ncgtos + i];

#pragma omp simd
                for (int64_t k = 0; k < ncols; k++)
                {
                    const auto gr_x = grid_coords_x[k] - r_x;

                    const auto gr_y = grid_coords_y[k] - r_y;

                    const auto gr_z = grid_coords_z[k] - r_z;

                    const auto fss = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

                    ptr_buffer_xx[k] += gr_x * gr_x * fss;

                    ptr_buffer_xy[k] += gr_x * gr_y * fss;

                    ptr_buffer_xz[k] += gr_x * gr_z * fss;

                    ptr_buffer_yy[k] += gr_y * gr_y * fss;

                    ptr_buffer_yz[k] += gr_y * gr_z * fss;

                    ptr_buffer_zz[k] += gr_z * gr_z * fss;
                }
            }

            // distribute GTO values into submatrix

            gtoval::distribute(submat, buffer_xx, -1.0, 2 * nrows + irow);

            gtoval::distribute(submat, buffer_xx, 0.5 * f2_3, 4 * nrows + irow);

            gtoval::distribute(submat, buffer_xy, f2_3, irow);

            gtoval::distribute(submat, buffer_xz, f2_3, 3 * nrows + irow);

            gtoval::distribute(submat, buffer_yy, -1.0, 2 * nrows + irow);

            gtoval::distribute(submat, buffer_yy, -0.5 * f2_3, 4 * nrows + irow);

            gtoval::distribute(submat, buffer_yz, f2_3, nrows + irow);

            gtoval::distribute(submat, buffer_zz, 2.0, 2 * nrows + irow);

            irow++;
        }
    }

    return gto_values;
}

}  // namespace gtoval
