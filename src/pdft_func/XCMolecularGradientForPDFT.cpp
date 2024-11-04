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

#include "XCMolecularGradientForPDFT.hpp"

#include <omp.h>

#include "DenseLinearAlgebra.hpp"
#include "DftSubMatrix.hpp"
#include "GtoFunc.hpp"
#include "GtoValues.hpp"
#include "PairDensityGridGenerator.hpp"

namespace xcgradpdft {  // xcgradpdft namespace

CDenseMatrix
integrateVxcPDFTGradientForLDA(const CMolecule&                molecule,
                               const CMolecularBasis&          basis,
                               const double*                   densityMatrixPointer,
                               const CDenseMatrix&             twoBodyDensityMatrix,
                               const CDenseMatrix&             activeMOs,
                               const CMolecularGrid&           molecularGrid,
                               const double                    screeningThresholdForGTOValues,
                               const CXCPairDensityFunctional& xcFunctional,
                               const double                    rs_omega)
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    auto nthreads = omp_get_max_threads();

    std::vector<CMultiTimer> omptimers(nthreads);

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // AO-to-atom mapping

    std::vector<int> ao_to_atom_ids(naos);

    _computeAOtoAtomMapping(ao_to_atom_ids, molecule, basis);

    // molecular gradient

    auto natoms = molecule.number_of_atoms();

    CDenseMatrix molgrad_threads(nthreads, natoms * 3);

    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    // density and functional derivatives

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(2 * max_npoints_per_box);
    std::vector<double> exc_data(1 * max_npoints_per_box); //Not needed but always provided for now
    std::vector<double> vrho_data(2 * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho  = rho_data.data();
    auto exc  = exc_data.data();
    auto vrho = vrho_data.data();

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    timer.stop("Preparation");

    for (int box_id = 0; box_id < counts.size(); box_id++)
    {
        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // pre-screening

        timer.start("GTO pre-screening");

        std::vector<std::vector<int>> cgto_mask_blocks, pre_ao_inds_blocks;

        std::vector<int> aoinds;

        for (const auto& gto_block : gto_blocks)
        {
            // 1st order GTO derivative
            auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(gto_block, 1, screeningThresholdForGTOValues, boxdim);

            cgto_mask_blocks.push_back(cgto_mask);

            pre_ao_inds_blocks.push_back(pre_ao_inds);

            for (const auto nu : pre_ao_inds)
            {
                aoinds.push_back(nu);
            }
        }

        const auto aocount = static_cast<int>(aoinds.size());

        timer.stop("GTO pre-screening");

        if (aocount == 0) continue;

        // GTO values on grid points

        timer.start("OMP GTO evaluation");

        CDenseMatrix mat_chi(aocount, npoints);
        CDenseMatrix mat_chi_x(aocount, npoints);
        CDenseMatrix mat_chi_y(aocount, npoints);
        CDenseMatrix mat_chi_z(aocount, npoints);

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            omptimers[thread_id].start("gtoeval");

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            const auto grid_x_ptr = xcoords + gridblockpos + grid_batch_offset;
            const auto grid_y_ptr = ycoords + gridblockpos + grid_batch_offset;
            const auto grid_z_ptr = zcoords + gridblockpos + grid_batch_offset;

            std::vector<double> grid_x(grid_x_ptr, grid_x_ptr + grid_batch_size);
            std::vector<double> grid_y(grid_y_ptr, grid_y_ptr + grid_batch_size);
            std::vector<double> grid_z(grid_z_ptr, grid_z_ptr + grid_batch_size);

            // go through GTO blocks

            for (size_t i_block = 0, idx = 0; i_block < gto_blocks.size(); i_block++)
            {
                const auto& gto_block = gto_blocks[i_block];

                const auto& cgto_mask = cgto_mask_blocks[i_block];

                const auto& pre_ao_inds = pre_ao_inds_blocks[i_block];

                auto cmat = gtoval::get_gto_values_for_gga(gto_block, grid_x, grid_y, grid_z, cgto_mask);

                if (cmat.is_empty()) continue;

                auto submat_0_ptr = cmat.sub_matrix({0, 0});
                auto submat_x_ptr = cmat.sub_matrix({1, 0});
                auto submat_y_ptr = cmat.sub_matrix({1, 1});
                auto submat_z_ptr = cmat.sub_matrix({1, 2});

                auto submat_0_data = submat_0_ptr->data();
                auto submat_x_data = submat_x_ptr->data();
                auto submat_y_data = submat_y_ptr->data();
                auto submat_z_data = submat_z_ptr->data();

                for (int nu = 0; nu < static_cast<int>(pre_ao_inds.size()); nu++, idx++)
                {
                    std::memcpy(mat_chi.row(idx) + grid_batch_offset, submat_0_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_x.row(idx) + grid_batch_offset, submat_x_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_y.row(idx) + grid_batch_offset, submat_y_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_z.row(idx) + grid_batch_offset, submat_z_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");
        }

        timer.stop("OMP GTO evaluation");

        // generate sub density matrix

        timer.start("Density matrix slicing");

        auto sub_dens_mat_a = dftsubmat::getSubDensityMatrix(densityMatrixPointer, aoinds, naos);

        auto sub_active_mos = dftsubmat::getSubMatrixByColumnSlicing(activeMOs, aoinds, naos);

        timer.stop("Density matrix slicing");

        // generate density and on-top pair density on the grid

        pairdengridgen::generatePairDensityForLDA(rho, mat_chi, sub_dens_mat_a, sub_active_mos, twoBodyDensityMatrix, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, npoints);
        CDenseMatrix dengrady(natoms, npoints);
        CDenseMatrix dengradz(natoms, npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(sub_dens_mat_a, mat_chi);

        timer.stop("Density grad. grid matmul");

        // Pair-density parts (this recomputes a lot of things)

        timer.start("Density grad mo pair");

        auto n_active = activeMOs.getNumberOfRows();

        //1) \phi_t(r) = C_mu^t \phi_\mu(r)
        CDenseMatrix mos_on_grid;
        if (n_active > 0)
        {
            mos_on_grid = denblas::multAB(sub_active_mos, mat_chi);
        }

        auto n_active2 = n_active * n_active;

        //2) \phi_tu(r) = \phi_t(r) \phi_u(r)
        CDenseMatrix mo_pair(n_active2, npoints);

        auto mo_pair_val = mo_pair.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            for (int32_t t = 0; t < n_active; t++)
            {
                auto MOt = mos_on_grid.row(t);

                auto t_offset = t * n_active * npoints;

                for (int32_t u = 0; u < n_active; u++)
                {
                    auto MOu = mos_on_grid.row(u);

                    auto tu_offset = t_offset + u * npoints;

                    #pragma omp simd aligned(mo_pair_val, MOt, MOu: VLX_ALIGN)
                    for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                    {
                        mo_pair_val[tu_offset + g] += MOt[g] * MOu[g];
                    }
                }
            }
        }

        timer.stop("Density grad mo pair");

        timer.start("Density grad pi matmul");

        //3) d_tu(r) = d_tuvw \phi_v(r) \phi_w(r)
        auto mat_d = denblas::multAB(twoBodyDensityMatrix, mo_pair);

        timer.stop("Density grad pi matmul");

        //4) g_t(r) = d_tu(r) phi_u(r)
        CDenseMatrix mat_g(n_active, npoints);

        auto g_val = mat_g.values();

        auto d_val = mat_d.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            for (int32_t v = 0; v < n_active; v++)
            {
                auto v_offset = v * npoints;

                for (int32_t w = 0; w < n_active; w++)
                {
                    auto vw_offset = (v*n_active + w) * npoints;

                    auto MOw = mos_on_grid.row(w);

                    #pragma omp simd aligned(g_val, d_val, MOw : VLX_ALIGN)
                    for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                    {
                        g_val[v_offset + g] += d_val[vw_offset + g] * MOw[g];
                    }
                }

            }
        }

        //5) k_mu(r) = C_mu^t g_t(r)
        CDenseMatrix mat_k(aocount, npoints);
        if (n_active > 0)
        {
            mat_k = denblas::multAtB(sub_active_mos, mat_g);
        }
        else
        {
            mat_k.zero();
        }

        timer.start("Density grad. grid rho");

        auto F_val = mat_F.values();
        auto k_val = mat_k.values();

        auto chi_x_val = mat_chi_x.values();
        auto chi_y_val = mat_chi_y.values();
        auto chi_z_val = mat_chi_z.values();

        auto gdenx = dengradx.values();
        auto gdeny = dengrady.values();
        auto gdenz = dengradz.values();

        CDenseMatrix dengradpi_x(natoms, npoints);
        CDenseMatrix dengradpi_y(natoms, npoints);
        CDenseMatrix dengradpi_z(natoms, npoints);
        auto gdenpi_x = dengradpi_x.values();
        auto gdenpi_y = dengradpi_y.values();
        auto gdenpi_z = dengradpi_z.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < aocount; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * npoints;

                auto nu_offset = nu * npoints;

                #pragma omp simd aligned(gdenx, gdeny, gdenz, gdenpi_x, gdenpi_y, gdenpi_z, F_val, k_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];
                    gdeny[atom_g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];
                    gdenz[atom_g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];
                    gdenpi_x[atom_g] -= 4.0 * k_val[nu_g] * chi_x_val[nu_g];
                    gdenpi_y[atom_g] -= 4.0 * k_val[nu_g] * chi_y_val[nu_g];
                    gdenpi_z[atom_g] -= 4.0 * k_val[nu_g] * chi_z_val[nu_g];
                }
            }
        }

        timer.stop("Density grad. grid rho");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_exc_vxc_for_plda(npoints, rho, exc, vrho, rs_omega);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {   
            auto thread_id = omp_get_thread_num();
            
            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);
            
            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);
            
            auto gatm = molgrad_threads.row(thread_id);
            
            for (int32_t iatom = 0; iatom < natoms; iatom++)
            {   
                auto atom_offset = iatom * npoints;
                
                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;
                
                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) aligned(local_weights, vrho, gdenx, gdeny, gdenz, gdenpi_x, gdenpi_y, gdenpi_z : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {   
                    auto atom_g = atom_offset + g;
                    
                    double prefac = local_weights[g] * vrho[2 * g + 0];
                    double prefacpi = local_weights[g] * vrho[2 * g + 1];
                    
                    gatmx += prefac * gdenx[atom_g] + prefacpi * gdenpi_x[atom_g];
                    gatmy += prefac * gdeny[atom_g] + prefacpi * gdenpi_y[atom_g];
                    gatmz += prefac * gdenz[atom_g] + prefacpi * gdenpi_z[atom_g];
                }

                gatm[iatom * 3 + 0] += gatmx;
                gatm[iatom * 3 + 1] += gatmy;
                gatm[iatom * 3 + 2] += gatmz;
            }
        }

        timer.stop("Accumulate gradient");
    }

    timer.stop("Total timing");

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;

    CDenseMatrix molgrad(natoms, 3);

    for (int32_t iatom = 0; iatom < natoms; iatom++)
    {
        for (int32_t thread_id = 0; thread_id < nthreads; thread_id++)
        {
            molgrad.row(iatom)[0] += molgrad_threads.row(thread_id)[iatom * 3 + 0];
            molgrad.row(iatom)[1] += molgrad_threads.row(thread_id)[iatom * 3 + 1];
            molgrad.row(iatom)[2] += molgrad_threads.row(thread_id)[iatom * 3 + 2];
        }
    }

    return molgrad;
}

CDenseMatrix
integrateVxcPDFTGradientForGGA(const CMolecule&                molecule,
                               const CMolecularBasis&          basis,
                               const double*                   densityMatrixPointer,
                               const CDenseMatrix&             twoBodyDensityMatrix,
                               const CDenseMatrix&             activeMOs,        
                               const CMolecularGrid&           molecularGrid,    
                               const double                    screeningThresholdForGTOValues,
                               const CXCPairDensityFunctional& xcFunctional,
                               const double                    rs_omega)
{
    return CDenseMatrix();
}

void
_computeAOtoAtomMapping(std::vector<int>& ao_to_atom_ids, const CMolecule& molecule, const CMolecularBasis& basis)
{
    auto natoms = molecule.number_of_atoms();

    auto max_angl = basis.max_angular_momentum();

    // azimuthal quantum number: s,p,d,f,...

    for (int angl = 0, aoidx = 0; angl <= max_angl; angl++)
    {
        auto nsph = angl * 2 + 1;

        // magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...

        for (int isph = 0; isph < nsph; isph++)
        {
            // atoms

            for (int atomidx = 0; atomidx < natoms; atomidx++)
            {
                auto nao = basis.number_of_basis_functions(std::vector<int>({atomidx}), angl);

                // atomic orbitals

                for (int iao = 0; iao < nao; iao++, aoidx++)
                {
                    ao_to_atom_ids[aoidx] = atomidx;
                }
            }
        }
    }
}

}   // namespace xcgradpdft
