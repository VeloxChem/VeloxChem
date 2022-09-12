//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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

#include "DensityGridGenerator.hpp"

#include <cstring>

#include "omp.h"

#include "DenseLinearAlgebra.hpp"

namespace dengridgen {  // dengridgen namespace

CDensityGrid
generateDensityGridForLDA(const int32_t       npoints,
                          const CDenseMatrix& gtoValues,
                          const CDenseMatrix& densityMatrix,
                          const xcfun         xcFunType,
                          CMultiTimer&        timer)
{
    CDensityGrid dengrid(npoints, 1, xcFunType, dengrid::ab);

    dengrid.zero();

    auto rhoa = dengrid.alphaDensity(0);

    auto rhob = dengrid.betaDensity(0);

    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    auto mat_F = denblas::multAB(densityMatrix, gtoValues);

    timer.stop("Density grid matmul");

    // eq.(27), JCTC 2021, 17, 1512-1521

    timer.start("Density grid rho");

    for (int32_t nu = 0; nu < gtoValues.getNumberOfRows(); nu++)
    {
        auto F_nu = mat_F.row(nu);

        auto chi_nu = gtoValues.row(nu);

        for (int32_t g = 0; g < npoints; g++)
        {
            rhoa[g] += F_nu[g] * chi_nu[g];
        }
    }

    std::memcpy(rhob, rhoa, npoints * sizeof(double));

    timer.stop("Density grid rho");

    return dengrid;
}

CDensityGrid
generateDensityGridForGGA(const int32_t       npoints,
                          const CDenseMatrix& gtoValues,
                          const CDenseMatrix& gtoValuesX,
                          const CDenseMatrix& gtoValuesY,
                          const CDenseMatrix& gtoValuesZ,
                          const CDenseMatrix& densityMatrix,
                          const xcfun         xcFunType,
                          CMultiTimer&        timer)
{
    CDensityGrid dengrid(npoints, 1, xcFunType, dengrid::ab);

    dengrid.zero();

    auto rhoa = dengrid.alphaDensity(0);

    auto rhob = dengrid.betaDensity(0);

    auto grada = dengrid.alphaDensityGradient(0);

    auto gradb = dengrid.betaDensityGradient(0);

    auto gradab = dengrid.mixedDensityGradient(0);

    auto gradax = dengrid.alphaDensityGradientX(0);

    auto graday = dengrid.alphaDensityGradientY(0);

    auto gradaz = dengrid.alphaDensityGradientZ(0);

    auto gradbx = dengrid.betaDensityGradientX(0);

    auto gradby = dengrid.betaDensityGradientY(0);

    auto gradbz = dengrid.betaDensityGradientZ(0);

    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    CDenseMatrix symmetricDensityMatrix(densityMatrix); 

    symmetricDensityMatrix.symmetrizeAndScale(0.5);

    auto mat_F = denblas::multAB(symmetricDensityMatrix, gtoValues);

    timer.stop("Density grid matmul");

    // eq.(27), JCTC 2021, 17, 1512-1521

    timer.start("Density grid rho");

    auto naos = gtoValues.getNumberOfRows();

    auto nthreads = omp_get_max_threads();

    auto F_val = mat_F.values();

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();

    auto chi_y_val = gtoValuesY.values();

    auto chi_z_val = gtoValuesZ.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(rhoa, gradax, graday, gradaz, F_val, \
                                     chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rhoa[g] += F_val[nu_offset + g] * chi_val[nu_offset + g];

                gradax[g] += 2.0 * F_val[nu_offset + g] * chi_x_val[nu_offset + g];

                graday[g] += 2.0 * F_val[nu_offset + g] * chi_y_val[nu_offset + g];

                gradaz[g] += 2.0 * F_val[nu_offset + g] * chi_z_val[nu_offset + g];
            }
        }

        #pragma omp simd aligned(rhoa, rhob, gradax, graday, gradaz, gradbx, gradby, gradbz, \
                                 grada, gradb, gradab : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            rhob[g] = rhoa[g];

            gradbx[g] = gradax[g];

            gradby[g] = graday[g];

            gradbz[g] = gradaz[g];

            grada[g] = std::sqrt(gradax[g] * gradax[g] + graday[g] * graday[g] + gradaz[g] * gradaz[g]);

            gradb[g] = std::sqrt(gradbx[g] * gradbx[g] + gradby[g] * gradby[g] + gradbz[g] * gradbz[g]);

            gradab[g] = gradax[g] * gradbx[g] + graday[g] * gradby[g] + gradaz[g] * gradbz[g];
        }
    }

    timer.stop("Density grid rho");

    return dengrid;
}

}  // namespace dengridgen
