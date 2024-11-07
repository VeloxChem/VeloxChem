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

#include "XCMolecularGradient.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "AODensityMatrix.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DensityGridGenerator.hpp"
#include "DftSubMatrix.hpp"
#include "ErrorHandler.hpp"
#include "FunctionalParser.hpp"
#include "GridScreener.hpp"
#include "GtoFunc.hpp"
#include "GtoValues.hpp"
#include "Prescreener.hpp"
#include "XCFunctional.hpp"

CXCMolecularGradient::CXCMolecularGradient()

    : _screeningThresholdForGTOValues(1.0e-12)
{
}

CDenseMatrix
CXCMolecularGradient::integrateVxcGradient(const CMolecule&        molecule,
                                           const CMolecularBasis&  basis,
                                           const std::vector<const double*>& gsDensityPointers,
                                           const CMolecularGrid&   molecularGrid,
                                           const std::string&      xcFuncLabel) const
{
    return integrateVxcGradient(molecule, basis, gsDensityPointers, gsDensityPointers, molecularGrid, xcFuncLabel);
}

CDenseMatrix
CXCMolecularGradient::integrateVxcGradient(const CMolecule&        molecule,
                                           const CMolecularBasis&  basis,
                                           const std::vector<const double*>& rwDensityPointers,
                                           const std::vector<const double*>& gsDensityPointers,
                                           const CMolecularGrid&   molecularGrid,
                                           const std::string&      xcFuncLabel) const
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (gsDensityPointers.size() == 1)
    {
        if (xcfuntype == xcfun::lda)
        {
            return _integrateVxcGradientForLDA(molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, fvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
            return _integrateVxcGradientForGGA(molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, fvxc);
        }
        else
        {
            std::string errxcfuntype("XCMolecularGradient.integrateVxcGradient: Only implemented for LDA/GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        if (xcfuntype == xcfun::lda)
        {
            return _integrateVxcGradientForLDAOpenShell(molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, fvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
            return _integrateVxcGradientForGGAOpenShell(molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, fvxc);
        }
        else
        {
            std::string errxcfuntype("XCMolecularGradient.integrateVxcGradient: Only implemented for LDA/GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }

    return CDenseMatrix();
}

auto
CXCMolecularGradient::_integrateVxcGradientForLDA(const CMolecule&        molecule,
                                                  const CMolecularBasis&  basis,
                                                  const std::vector<const double*>& rwDensityPointers,
                                                  const std::vector<const double*>& gsDensityPointers,
                                                  const CMolecularGrid&   molecularGrid,
                                                  const CXCFunctional&    xcFunctional) const -> CDenseMatrix
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

    auto       ldafunc = xcFunctional.getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);

    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho = rho_data.data();

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
            auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(gto_block, 1, _screeningThresholdForGTOValues, boxdim);

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

        // generate sub density matrix and density grid

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);
        auto rw_sub_dens_mat = dftsubmat::getSubDensityMatrix(rwDensityPointers[0], aoinds, naos);

        timer.stop("Density matrix slicing");

        dengridgen::generateDensityForLDA(rho, mat_chi, gs_sub_dens_mat, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, npoints);
        CDenseMatrix dengrady(natoms, npoints);
        CDenseMatrix dengradz(natoms, npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(rw_sub_dens_mat, mat_chi);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto F_val = mat_F.values();

        auto chi_x_val = mat_chi_x.values();
        auto chi_y_val = mat_chi_y.values();
        auto chi_z_val = mat_chi_z.values();

        auto gdenx = dengradx.values();
        auto gdeny = dengrady.values();
        auto gdenz = dengradz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            for (int nu = 0; nu < aocount; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * npoints;

                auto nu_offset = nu * npoints;

                #pragma omp simd
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];
                    gdeny[atom_g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];
                    gdenz[atom_g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];
                }
            }
        }

        timer.stop("Density grad. grid rho");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_lda(npoints, rho, vrho);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.row(thread_id);

            for (int iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz)
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac = local_weights[g] * vrho[2 * g + 0];

                    gatmx += prefac * gdenx[atom_g];
                    gatmy += prefac * gdeny[atom_g];
                    gatmz += prefac * gdenz[atom_g];
                }

                // factor of 2 from sum of alpha and beta contributions

                gatm[iatom * 3 + 0] += 2.0 * gatmx;
                gatm[iatom * 3 + 1] += 2.0 * gatmy;
                gatm[iatom * 3 + 2] += 2.0 * gatmz;
            }
        }

        timer.stop("Accumulate gradient");
    }

    timer.stop("Total timing");

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;

    CDenseMatrix molgrad(natoms, 3);

    for (int iatom = 0; iatom < natoms; iatom++)
    {
        for (int thread_id = 0; thread_id < nthreads; thread_id++)
        {
            molgrad.row(iatom)[0] += molgrad_threads.row(thread_id)[iatom * 3 + 0];
            molgrad.row(iatom)[1] += molgrad_threads.row(thread_id)[iatom * 3 + 1];
            molgrad.row(iatom)[2] += molgrad_threads.row(thread_id)[iatom * 3 + 2];
        }
    }

    return molgrad;
}

auto
CXCMolecularGradient::_integrateVxcGradientForLDAOpenShell(const CMolecule&        molecule,
                                                           const CMolecularBasis&  basis,
                                                           const std::vector<const double*>& rwDensityPointers,
                                                           const std::vector<const double*>& gsDensityPointers,
                                                           const CMolecularGrid&   molecularGrid,
                                                           const CXCFunctional&    xcFunctional) const -> CDenseMatrix
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

    std::vector<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    std::vector<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    std::vector<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho  = rho_data.data();
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
            auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(gto_block, 1, _screeningThresholdForGTOValues, boxdim);

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

        auto gs_sub_dens_mat_a = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);
        auto gs_sub_dens_mat_b = dftsubmat::getSubDensityMatrix(gsDensityPointers[1], aoinds, naos);

        auto rw_sub_dens_mat_a = dftsubmat::getSubDensityMatrix(rwDensityPointers[0], aoinds, naos);
        auto rw_sub_dens_mat_b = dftsubmat::getSubDensityMatrix(rwDensityPointers[1], aoinds, naos);

        timer.stop("Density matrix slicing");

        // generate density grid

        dengridgen::generateDensityForLDA(rho, mat_chi, gs_sub_dens_mat_a, gs_sub_dens_mat_b, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengrad_a_x(natoms, npoints);
        CDenseMatrix dengrad_a_y(natoms, npoints);
        CDenseMatrix dengrad_a_z(natoms, npoints);

        CDenseMatrix dengrad_b_x(natoms, npoints);
        CDenseMatrix dengrad_b_y(natoms, npoints);
        CDenseMatrix dengrad_b_z(natoms, npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F_a = denblas::multAB(rw_sub_dens_mat_a, mat_chi);
        auto mat_F_b = denblas::multAB(rw_sub_dens_mat_b, mat_chi);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto F_a_val = mat_F_a.values();
        auto F_b_val = mat_F_b.values();

        auto chi_x_val = mat_chi_x.values();
        auto chi_y_val = mat_chi_y.values();
        auto chi_z_val = mat_chi_z.values();

        auto gden_a_x = dengrad_a_x.values();
        auto gden_a_y = dengrad_a_y.values();
        auto gden_a_z = dengrad_a_z.values();

        auto gden_b_x = dengrad_b_x.values();
        auto gden_b_y = dengrad_b_y.values();
        auto gden_b_z = dengrad_b_z.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            for (int nu = 0; nu < aocount; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * npoints;

                auto nu_offset = nu * npoints;

                #pragma omp simd 
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gden_a_x[atom_g] -= 2.0 * F_a_val[nu_g] * chi_x_val[nu_g];
                    gden_a_y[atom_g] -= 2.0 * F_a_val[nu_g] * chi_y_val[nu_g];
                    gden_a_z[atom_g] -= 2.0 * F_a_val[nu_g] * chi_z_val[nu_g];

                    gden_b_x[atom_g] -= 2.0 * F_b_val[nu_g] * chi_x_val[nu_g];
                    gden_b_y[atom_g] -= 2.0 * F_b_val[nu_g] * chi_y_val[nu_g];
                    gden_b_z[atom_g] -= 2.0 * F_b_val[nu_g] * chi_z_val[nu_g];
                }
            }
        }

        timer.stop("Density grad. grid rho");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_lda(npoints, rho, vrho);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.row(thread_id);

            for (int iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) 
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    double prefac_a = local_weights[g] * vrho[2 * g + 0];
                    double prefac_b = local_weights[g] * vrho[2 * g + 1];

                    gatmx += prefac_a * gden_a_x[atom_g];
                    gatmy += prefac_a * gden_a_y[atom_g];
                    gatmz += prefac_a * gden_a_z[atom_g];

                    gatmx += prefac_b * gden_b_x[atom_g];
                    gatmy += prefac_b * gden_b_y[atom_g];
                    gatmz += prefac_b * gden_b_z[atom_g];
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

    for (int iatom = 0; iatom < natoms; iatom++)
    {
        for (int thread_id = 0; thread_id < nthreads; thread_id++)
        {
            molgrad.row(iatom)[0] += molgrad_threads.row(thread_id)[iatom * 3 + 0];
            molgrad.row(iatom)[1] += molgrad_threads.row(thread_id)[iatom * 3 + 1];
            molgrad.row(iatom)[2] += molgrad_threads.row(thread_id)[iatom * 3 + 2];
        }
    }

    return molgrad;
}

CDenseMatrix
CXCMolecularGradient::_integrateVxcGradientForGGA(const CMolecule&        molecule,
                                                  const CMolecularBasis&  basis,
                                                  const std::vector<const double*>& rwDensityPointers,
                                                  const std::vector<const double*>& gsDensityPointers,
                                                  const CMolecularGrid&   molecularGrid,
                                                  const CXCFunctional&    xcFunctional) const
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

    auto       ggafunc = xcFunctional.getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    std::vector<double> local_weights_data(max_npoints_per_box);

    std::vector<double> rho_data(dim->rho * max_npoints_per_box);
    std::vector<double> rhograd_data(dim->rho * 3 * max_npoints_per_box);
    std::vector<double> sigma_data(dim->sigma * max_npoints_per_box);

    std::vector<double> vrho_data(dim->vrho * max_npoints_per_box);
    std::vector<double> vsigma_data(dim->vsigma * max_npoints_per_box);

    auto local_weights = local_weights_data.data();

    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();

    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();

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

        // pre-screening of GTOs

        timer.start("GTO pre-screening");

        std::vector<std::vector<int>> cgto_mask_blocks, pre_ao_inds_blocks;

        std::vector<int> aoinds;

        for (const auto& gto_block : gto_blocks)
        {
            // 2nd order GTO derivative
            auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(gto_block, 2, _screeningThresholdForGTOValues, boxdim);

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

        CDenseMatrix mat_chi_xx(aocount, npoints);
        CDenseMatrix mat_chi_xy(aocount, npoints);
        CDenseMatrix mat_chi_xz(aocount, npoints);
        CDenseMatrix mat_chi_yy(aocount, npoints);
        CDenseMatrix mat_chi_yz(aocount, npoints);
        CDenseMatrix mat_chi_zz(aocount, npoints);

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

                auto cmat = gtoval::get_gto_values_for_mgga(gto_block, grid_x, grid_y, grid_z, cgto_mask);

                if (cmat.is_empty()) continue;

                auto submat_0_ptr = cmat.sub_matrix({0, 0});

                auto submat_x_ptr = cmat.sub_matrix({1, 0});
                auto submat_y_ptr = cmat.sub_matrix({1, 1});
                auto submat_z_ptr = cmat.sub_matrix({1, 2});

                auto submat_xx_ptr = cmat.sub_matrix({2, 0});
                auto submat_xy_ptr = cmat.sub_matrix({2, 1});
                auto submat_xz_ptr = cmat.sub_matrix({2, 2});
                auto submat_yy_ptr = cmat.sub_matrix({2, 3});
                auto submat_yz_ptr = cmat.sub_matrix({2, 4});
                auto submat_zz_ptr = cmat.sub_matrix({2, 5});

                auto submat_0_data = submat_0_ptr->data();

                auto submat_x_data = submat_x_ptr->data();
                auto submat_y_data = submat_y_ptr->data();
                auto submat_z_data = submat_z_ptr->data();

                auto submat_xx_data = submat_xx_ptr->data();
                auto submat_xy_data = submat_xy_ptr->data();
                auto submat_xz_data = submat_xz_ptr->data();
                auto submat_yy_data = submat_yy_ptr->data();
                auto submat_yz_data = submat_yz_ptr->data();
                auto submat_zz_data = submat_zz_ptr->data();

                for (int nu = 0; nu < static_cast<int>(pre_ao_inds.size()); nu++, idx++)
                {
                    std::memcpy(mat_chi.row(idx) + grid_batch_offset, submat_0_data + nu * grid_batch_size, grid_batch_size * sizeof(double));

                    std::memcpy(mat_chi_x.row(idx) + grid_batch_offset, submat_x_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_y.row(idx) + grid_batch_offset, submat_y_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_z.row(idx) + grid_batch_offset, submat_z_data + nu * grid_batch_size, grid_batch_size * sizeof(double));

                    std::memcpy(mat_chi_xx.row(idx) + grid_batch_offset, submat_xx_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xy.row(idx) + grid_batch_offset, submat_xy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xz.row(idx) + grid_batch_offset, submat_xz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yy.row(idx) + grid_batch_offset, submat_yy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yz.row(idx) + grid_batch_offset, submat_yz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_zz.row(idx) + grid_batch_offset, submat_zz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");
        }

        timer.stop("OMP GTO evaluation");

        // generate sub density matrix and density grid

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);
        auto rw_sub_dens_mat = dftsubmat::getSubDensityMatrix(rwDensityPointers[0], aoinds, naos);

        timer.stop("Density matrix slicing");

        dengridgen::generateDensityForGGA(rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengradx(natoms, npoints);
        CDenseMatrix dengrady(natoms, npoints);
        CDenseMatrix dengradz(natoms, npoints);

        CDenseMatrix dengradxx(natoms, npoints);
        CDenseMatrix dengradxy(natoms, npoints);
        CDenseMatrix dengradxz(natoms, npoints);

        CDenseMatrix dengradyx(natoms, npoints);
        CDenseMatrix dengradyy(natoms, npoints);
        CDenseMatrix dengradyz(natoms, npoints);

        CDenseMatrix dengradzx(natoms, npoints);
        CDenseMatrix dengradzy(natoms, npoints);
        CDenseMatrix dengradzz(natoms, npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F = denblas::multAB(rw_sub_dens_mat, mat_chi);

        auto mat_F_x = denblas::multAB(rw_sub_dens_mat, mat_chi_x);
        auto mat_F_y = denblas::multAB(rw_sub_dens_mat, mat_chi_y);
        auto mat_F_z = denblas::multAB(rw_sub_dens_mat, mat_chi_z);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto F_val = mat_F.values();

        auto F_x_val = mat_F_x.values();
        auto F_y_val = mat_F_y.values();
        auto F_z_val = mat_F_z.values();

        auto chi_x_val = mat_chi_x.values();
        auto chi_y_val = mat_chi_y.values();
        auto chi_z_val = mat_chi_z.values();

        auto chi_xx_val = mat_chi_xx.values();
        auto chi_xy_val = mat_chi_xy.values();
        auto chi_xz_val = mat_chi_xz.values();

        auto chi_yy_val = mat_chi_yy.values();
        auto chi_yz_val = mat_chi_yz.values();
        auto chi_zz_val = mat_chi_zz.values();

        auto gdenx = dengradx.values();
        auto gdeny = dengrady.values();
        auto gdenz = dengradz.values();

        auto gdenxx = dengradxx.values();
        auto gdenxy = dengradxy.values();
        auto gdenxz = dengradxz.values();

        auto gdenyx = dengradyx.values();
        auto gdenyy = dengradyy.values();
        auto gdenyz = dengradyz.values();

        auto gdenzx = dengradzx.values();
        auto gdenzy = dengradzy.values();
        auto gdenzz = dengradzz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            for (int nu = 0; nu < aocount; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * npoints;

                auto nu_offset = nu * npoints;

                #pragma omp simd
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gdenx[atom_g] -= 2.0 * F_val[nu_g] * chi_x_val[nu_g];
                    gdeny[atom_g] -= 2.0 * F_val[nu_g] * chi_y_val[nu_g];
                    gdenz[atom_g] -= 2.0 * F_val[nu_g] * chi_z_val[nu_g];

                    gdenxx[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xx_val[nu_g]);
                    gdenxy[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_xy_val[nu_g]);
                    gdenxz[atom_g] -= 2.0 * (F_x_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_xz_val[nu_g]);

                    gdenyx[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xy_val[nu_g]);
                    gdenyy[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_yy_val[nu_g]);
                    gdenyz[atom_g] -= 2.0 * (F_y_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_yz_val[nu_g]);

                    gdenzx[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_x_val[nu_g] + F_val[nu_g] * chi_xz_val[nu_g]);
                    gdenzy[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_y_val[nu_g] + F_val[nu_g] * chi_yz_val[nu_g]);
                    gdenzz[atom_g] -= 2.0 * (F_z_val[nu_g] * chi_z_val[nu_g] + F_val[nu_g] * chi_zz_val[nu_g]);
                }
            }
        }

        timer.stop("Density grad. grid rho");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.row(thread_id);

            for (int iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz)
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto vx = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                    auto vy = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                    auto vz = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                    gatmx += local_weights[g] * vrho[2 * g + 0] * gdenx[atom_g];
                    gatmy += local_weights[g] * vrho[2 * g + 0] * gdeny[atom_g];
                    gatmz += local_weights[g] * vrho[2 * g + 0] * gdenz[atom_g];

                    gatmx += local_weights[g] * (vx * gdenxx[atom_g] + vy * gdenyx[atom_g] + vz * gdenzx[atom_g]);
                    gatmy += local_weights[g] * (vx * gdenxy[atom_g] + vy * gdenyy[atom_g] + vz * gdenzy[atom_g]);
                    gatmz += local_weights[g] * (vx * gdenxz[atom_g] + vy * gdenyz[atom_g] + vz * gdenzz[atom_g]);
                }

                // factor of 2 from sum of alpha and beta contributions

                gatm[iatom * 3 + 0] += 2.0 * gatmx;
                gatm[iatom * 3 + 1] += 2.0 * gatmy;
                gatm[iatom * 3 + 2] += 2.0 * gatmz;
            }
        }

        timer.stop("Accumulate gradient");
    }

    timer.stop("Total timing");

    // std::cout << "Timing of new integrator" << std::endl;
    // std::cout << "------------------------" << std::endl;
    // std::cout << timer.getSummary() << std::endl;

    CDenseMatrix molgrad(natoms, 3);

    for (int iatom = 0; iatom < natoms; iatom++)
    {
        for (int thread_id = 0; thread_id < nthreads; thread_id++)
        {
            molgrad.row(iatom)[0] += molgrad_threads.row(thread_id)[iatom * 3 + 0];
            molgrad.row(iatom)[1] += molgrad_threads.row(thread_id)[iatom * 3 + 1];
            molgrad.row(iatom)[2] += molgrad_threads.row(thread_id)[iatom * 3 + 2];
        }
    }

    return molgrad;
}

auto
CXCMolecularGradient::_integrateVxcGradientForGGAOpenShell(const CMolecule&        molecule,
                                                           const CMolecularBasis&  basis,
                                                           const std::vector<const double*>& rwDensityPointers,
                                                           const std::vector<const double*>& gsDensityPointers,
                                                           const CMolecularGrid&   molecularGrid,
                                                           const CXCFunctional&    xcFunctional) const -> CDenseMatrix
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

    std::vector<double> local_weights_data(molecularGrid.getMaxNumberOfGridPointsPerBox());

    std::vector<double> rho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    std::vector<double> rhograd_data(6 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    std::vector<double> sigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    std::vector<double> vrho_data(2 * molecularGrid.getMaxNumberOfGridPointsPerBox());
    std::vector<double> vsigma_data(3 * molecularGrid.getMaxNumberOfGridPointsPerBox());

    auto local_weights = local_weights_data.data();

    auto rho     = rho_data.data();
    auto rhograd = rhograd_data.data();
    auto sigma   = sigma_data.data();

    auto vrho   = vrho_data.data();
    auto vsigma = vsigma_data.data();

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

        // pre-screening of GTOs

        timer.start("GTO pre-screening");

        std::vector<std::vector<int>> cgto_mask_blocks, pre_ao_inds_blocks;

        std::vector<int> aoinds;

        for (const auto& gto_block : gto_blocks)
        {
            // 2nd order GTO derivative
            auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(gto_block, 2, _screeningThresholdForGTOValues, boxdim);

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

        CDenseMatrix mat_chi_xx(aocount, npoints);
        CDenseMatrix mat_chi_xy(aocount, npoints);
        CDenseMatrix mat_chi_xz(aocount, npoints);
        CDenseMatrix mat_chi_yy(aocount, npoints);
        CDenseMatrix mat_chi_yz(aocount, npoints);
        CDenseMatrix mat_chi_zz(aocount, npoints);

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

                auto cmat = gtoval::get_gto_values_for_mgga(gto_block, grid_x, grid_y, grid_z, cgto_mask);

                if (cmat.is_empty()) continue;

                auto submat_0_ptr = cmat.sub_matrix({0, 0});

                auto submat_x_ptr = cmat.sub_matrix({1, 0});
                auto submat_y_ptr = cmat.sub_matrix({1, 1});
                auto submat_z_ptr = cmat.sub_matrix({1, 2});

                auto submat_xx_ptr = cmat.sub_matrix({2, 0});
                auto submat_xy_ptr = cmat.sub_matrix({2, 1});
                auto submat_xz_ptr = cmat.sub_matrix({2, 2});
                auto submat_yy_ptr = cmat.sub_matrix({2, 3});
                auto submat_yz_ptr = cmat.sub_matrix({2, 4});
                auto submat_zz_ptr = cmat.sub_matrix({2, 5});

                auto submat_0_data = submat_0_ptr->data();

                auto submat_x_data = submat_x_ptr->data();
                auto submat_y_data = submat_y_ptr->data();
                auto submat_z_data = submat_z_ptr->data();

                auto submat_xx_data = submat_xx_ptr->data();
                auto submat_xy_data = submat_xy_ptr->data();
                auto submat_xz_data = submat_xz_ptr->data();
                auto submat_yy_data = submat_yy_ptr->data();
                auto submat_yz_data = submat_yz_ptr->data();
                auto submat_zz_data = submat_zz_ptr->data();

                for (int nu = 0; nu < static_cast<int>(pre_ao_inds.size()); nu++, idx++)
                {
                    std::memcpy(mat_chi.row(idx) + grid_batch_offset, submat_0_data + nu * grid_batch_size, grid_batch_size * sizeof(double));

                    std::memcpy(mat_chi_x.row(idx) + grid_batch_offset, submat_x_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_y.row(idx) + grid_batch_offset, submat_y_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_z.row(idx) + grid_batch_offset, submat_z_data + nu * grid_batch_size, grid_batch_size * sizeof(double));

                    std::memcpy(mat_chi_xx.row(idx) + grid_batch_offset, submat_xx_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xy.row(idx) + grid_batch_offset, submat_xy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_xz.row(idx) + grid_batch_offset, submat_xz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yy.row(idx) + grid_batch_offset, submat_yy_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_yz.row(idx) + grid_batch_offset, submat_yz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                    std::memcpy(mat_chi_zz.row(idx) + grid_batch_offset, submat_zz_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }

            omptimers[thread_id].stop("gtoeval");
        }

        timer.stop("OMP GTO evaluation");

        // generate sub density matrix and density grid

        timer.start("Density matrix slicing");

        auto gs_sub_dens_mat_a = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);
        auto gs_sub_dens_mat_b = dftsubmat::getSubDensityMatrix(gsDensityPointers[1], aoinds, naos);

        auto rw_sub_dens_mat_a = dftsubmat::getSubDensityMatrix(rwDensityPointers[0], aoinds, naos);
        auto rw_sub_dens_mat_b = dftsubmat::getSubDensityMatrix(rwDensityPointers[1], aoinds, naos);

        timer.stop("Density matrix slicing");

        dengridgen::generateDensityForGGA(
            rho, rhograd, sigma, mat_chi, mat_chi_x, mat_chi_y, mat_chi_z, gs_sub_dens_mat_a, gs_sub_dens_mat_b, timer);

        // generate density gradient grid

        timer.start("Density grad. grid prep.");

        CDenseMatrix dengrad_a_x(natoms, npoints);
        CDenseMatrix dengrad_a_y(natoms, npoints);
        CDenseMatrix dengrad_a_z(natoms, npoints);

        CDenseMatrix dengrad_b_x(natoms, npoints);
        CDenseMatrix dengrad_b_y(natoms, npoints);
        CDenseMatrix dengrad_b_z(natoms, npoints);

        CDenseMatrix dengrad_a_xx(natoms, npoints);
        CDenseMatrix dengrad_a_xy(natoms, npoints);
        CDenseMatrix dengrad_a_xz(natoms, npoints);

        CDenseMatrix dengrad_b_xx(natoms, npoints);
        CDenseMatrix dengrad_b_xy(natoms, npoints);
        CDenseMatrix dengrad_b_xz(natoms, npoints);

        CDenseMatrix dengrad_a_yx(natoms, npoints);
        CDenseMatrix dengrad_a_yy(natoms, npoints);
        CDenseMatrix dengrad_a_yz(natoms, npoints);

        CDenseMatrix dengrad_b_yx(natoms, npoints);
        CDenseMatrix dengrad_b_yy(natoms, npoints);
        CDenseMatrix dengrad_b_yz(natoms, npoints);

        CDenseMatrix dengrad_a_zx(natoms, npoints);
        CDenseMatrix dengrad_a_zy(natoms, npoints);
        CDenseMatrix dengrad_a_zz(natoms, npoints);

        CDenseMatrix dengrad_b_zx(natoms, npoints);
        CDenseMatrix dengrad_b_zy(natoms, npoints);
        CDenseMatrix dengrad_b_zz(natoms, npoints);

        timer.stop("Density grad. grid prep.");

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid matmul");

        auto mat_F_a = denblas::multAB(rw_sub_dens_mat_a, mat_chi);
        auto mat_F_b = denblas::multAB(rw_sub_dens_mat_b, mat_chi);

        auto mat_F_a_x = denblas::multAB(rw_sub_dens_mat_a, mat_chi_x);
        auto mat_F_a_y = denblas::multAB(rw_sub_dens_mat_a, mat_chi_y);
        auto mat_F_a_z = denblas::multAB(rw_sub_dens_mat_a, mat_chi_z);

        auto mat_F_b_x = denblas::multAB(rw_sub_dens_mat_b, mat_chi_x);
        auto mat_F_b_y = denblas::multAB(rw_sub_dens_mat_b, mat_chi_y);
        auto mat_F_b_z = denblas::multAB(rw_sub_dens_mat_b, mat_chi_z);

        timer.stop("Density grad. grid matmul");

        // eq.(34), JCTC 2021, 17, 1512-1521

        timer.start("Density grad. grid rho");

        auto F_a_val = mat_F_a.values();
        auto F_b_val = mat_F_b.values();

        auto F_a_x_val = mat_F_a_x.values();
        auto F_a_y_val = mat_F_a_y.values();
        auto F_a_z_val = mat_F_a_z.values();

        auto F_b_x_val = mat_F_b_x.values();
        auto F_b_y_val = mat_F_b_y.values();
        auto F_b_z_val = mat_F_b_z.values();

        auto chi_x_val = mat_chi_x.values();
        auto chi_y_val = mat_chi_y.values();
        auto chi_z_val = mat_chi_z.values();

        auto chi_xx_val = mat_chi_xx.values();
        auto chi_xy_val = mat_chi_xy.values();
        auto chi_xz_val = mat_chi_xz.values();

        auto chi_yy_val = mat_chi_yy.values();
        auto chi_yz_val = mat_chi_yz.values();
        auto chi_zz_val = mat_chi_zz.values();

        auto gden_a_x = dengrad_a_x.values();
        auto gden_a_y = dengrad_a_y.values();
        auto gden_a_z = dengrad_a_z.values();

        auto gden_b_x = dengrad_b_x.values();
        auto gden_b_y = dengrad_b_y.values();
        auto gden_b_z = dengrad_b_z.values();

        auto gden_a_xx = dengrad_a_xx.values();
        auto gden_a_xy = dengrad_a_xy.values();
        auto gden_a_xz = dengrad_a_xz.values();

        auto gden_b_xx = dengrad_b_xx.values();
        auto gden_b_xy = dengrad_b_xy.values();
        auto gden_b_xz = dengrad_b_xz.values();

        auto gden_a_yx = dengrad_a_yx.values();
        auto gden_a_yy = dengrad_a_yy.values();
        auto gden_a_yz = dengrad_a_yz.values();

        auto gden_b_yx = dengrad_b_yx.values();
        auto gden_b_yy = dengrad_b_yy.values();
        auto gden_b_yz = dengrad_b_yz.values();

        auto gden_a_zx = dengrad_a_zx.values();
        auto gden_a_zy = dengrad_a_zy.values();
        auto gden_a_zz = dengrad_a_zz.values();

        auto gden_b_zx = dengrad_b_zx.values();
        auto gden_b_zy = dengrad_b_zy.values();
        auto gden_b_zz = dengrad_b_zz.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            for (int nu = 0; nu < aocount; nu++)
            {
                auto atomidx = ao_to_atom_ids[aoinds[nu]];

                auto atom_offset = atomidx * npoints;

                auto nu_offset = nu * npoints;

                #pragma omp simd 
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto nu_g = nu_offset + g;

                    gden_a_x[atom_g] -= 2.0 * F_a_val[nu_g] * chi_x_val[nu_g];
                    gden_a_y[atom_g] -= 2.0 * F_a_val[nu_g] * chi_y_val[nu_g];
                    gden_a_z[atom_g] -= 2.0 * F_a_val[nu_g] * chi_z_val[nu_g];

                    gden_b_x[atom_g] -= 2.0 * F_b_val[nu_g] * chi_x_val[nu_g];
                    gden_b_y[atom_g] -= 2.0 * F_b_val[nu_g] * chi_y_val[nu_g];
                    gden_b_z[atom_g] -= 2.0 * F_b_val[nu_g] * chi_z_val[nu_g];

                    gden_a_xx[atom_g] -= 2.0 * (F_a_x_val[nu_g] * chi_x_val[nu_g] + F_a_val[nu_g] * chi_xx_val[nu_g]);
                    gden_a_xy[atom_g] -= 2.0 * (F_a_x_val[nu_g] * chi_y_val[nu_g] + F_a_val[nu_g] * chi_xy_val[nu_g]);
                    gden_a_xz[atom_g] -= 2.0 * (F_a_x_val[nu_g] * chi_z_val[nu_g] + F_a_val[nu_g] * chi_xz_val[nu_g]);

                    gden_b_xx[atom_g] -= 2.0 * (F_b_x_val[nu_g] * chi_x_val[nu_g] + F_b_val[nu_g] * chi_xx_val[nu_g]);
                    gden_b_xy[atom_g] -= 2.0 * (F_b_x_val[nu_g] * chi_y_val[nu_g] + F_b_val[nu_g] * chi_xy_val[nu_g]);
                    gden_b_xz[atom_g] -= 2.0 * (F_b_x_val[nu_g] * chi_z_val[nu_g] + F_b_val[nu_g] * chi_xz_val[nu_g]);

                    gden_a_yx[atom_g] -= 2.0 * (F_a_y_val[nu_g] * chi_x_val[nu_g] + F_a_val[nu_g] * chi_xy_val[nu_g]);
                    gden_a_yy[atom_g] -= 2.0 * (F_a_y_val[nu_g] * chi_y_val[nu_g] + F_a_val[nu_g] * chi_yy_val[nu_g]);
                    gden_a_yz[atom_g] -= 2.0 * (F_a_y_val[nu_g] * chi_z_val[nu_g] + F_a_val[nu_g] * chi_yz_val[nu_g]);

                    gden_b_yx[atom_g] -= 2.0 * (F_b_y_val[nu_g] * chi_x_val[nu_g] + F_b_val[nu_g] * chi_xy_val[nu_g]);
                    gden_b_yy[atom_g] -= 2.0 * (F_b_y_val[nu_g] * chi_y_val[nu_g] + F_b_val[nu_g] * chi_yy_val[nu_g]);
                    gden_b_yz[atom_g] -= 2.0 * (F_b_y_val[nu_g] * chi_z_val[nu_g] + F_b_val[nu_g] * chi_yz_val[nu_g]);

                    gden_a_zx[atom_g] -= 2.0 * (F_a_z_val[nu_g] * chi_x_val[nu_g] + F_a_val[nu_g] * chi_xz_val[nu_g]);
                    gden_a_zy[atom_g] -= 2.0 * (F_a_z_val[nu_g] * chi_y_val[nu_g] + F_a_val[nu_g] * chi_yz_val[nu_g]);
                    gden_a_zz[atom_g] -= 2.0 * (F_a_z_val[nu_g] * chi_z_val[nu_g] + F_a_val[nu_g] * chi_zz_val[nu_g]);

                    gden_b_zx[atom_g] -= 2.0 * (F_b_z_val[nu_g] * chi_x_val[nu_g] + F_b_val[nu_g] * chi_xz_val[nu_g]);
                    gden_b_zy[atom_g] -= 2.0 * (F_b_z_val[nu_g] * chi_y_val[nu_g] + F_b_val[nu_g] * chi_yz_val[nu_g]);
                    gden_b_zz[atom_g] -= 2.0 * (F_b_z_val[nu_g] * chi_z_val[nu_g] + F_b_val[nu_g] * chi_zz_val[nu_g]);
                }
            }
        }

        timer.stop("Density grad. grid rho");

        // compute exchange-correlation functional derivative

        timer.start("XC functional eval.");

        xcFunctional.compute_vxc_for_gga(npoints, rho, sigma, vrho, vsigma);

        std::memcpy(local_weights, weights + gridblockpos, npoints * sizeof(double));

        timer.stop("XC functional eval.");

        // eq.(32), JCTC 2021, 17, 1512-1521

        timer.start("Accumulate gradient");

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            auto gatm = molgrad_threads.row(thread_id);

            for (int iatom = 0; iatom < natoms; iatom++)
            {
                auto atom_offset = iatom * npoints;

                double gatmx = 0.0, gatmy = 0.0, gatmz = 0.0;

                #pragma omp simd reduction(+ : gatmx, gatmy, gatmz) 
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    auto atom_g = atom_offset + g;

                    auto vxa = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 0] + vsigma[3 * g + 1] * rhograd[6 * g + 3];
                    auto vya = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 1] + vsigma[3 * g + 1] * rhograd[6 * g + 4];
                    auto vza = 2.0 * vsigma[3 * g + 0] * rhograd[6 * g + 2] + vsigma[3 * g + 1] * rhograd[6 * g + 5];

                    auto vxb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 3] + vsigma[3 * g + 1] * rhograd[6 * g + 0];
                    auto vyb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 4] + vsigma[3 * g + 1] * rhograd[6 * g + 1];
                    auto vzb = 2.0 * vsigma[3 * g + 2] * rhograd[6 * g + 5] + vsigma[3 * g + 1] * rhograd[6 * g + 2];

                    // alpha
                    gatmx += local_weights[g] * vrho[2 * g + 0] * gden_a_x[atom_g];
                    gatmy += local_weights[g] * vrho[2 * g + 0] * gden_a_y[atom_g];
                    gatmz += local_weights[g] * vrho[2 * g + 0] * gden_a_z[atom_g];

                    // beta
                    gatmx += local_weights[g] * vrho[2 * g + 1] * gden_b_x[atom_g];
                    gatmy += local_weights[g] * vrho[2 * g + 1] * gden_b_y[atom_g];
                    gatmz += local_weights[g] * vrho[2 * g + 1] * gden_b_z[atom_g];

                    // alpha
                    gatmx += local_weights[g] * (vxa * gden_a_xx[atom_g] + vya * gden_a_yx[atom_g] + vza * gden_a_zx[atom_g]);
                    gatmy += local_weights[g] * (vxa * gden_a_xy[atom_g] + vya * gden_a_yy[atom_g] + vza * gden_a_zy[atom_g]);
                    gatmz += local_weights[g] * (vxa * gden_a_xz[atom_g] + vya * gden_a_yz[atom_g] + vza * gden_a_zz[atom_g]);

                    // beta
                    gatmx += local_weights[g] * (vxb * gden_b_xx[atom_g] + vyb * gden_b_yx[atom_g] + vzb * gden_b_zx[atom_g]);
                    gatmy += local_weights[g] * (vxb * gden_b_xy[atom_g] + vyb * gden_b_yy[atom_g] + vzb * gden_b_zy[atom_g]);
                    gatmz += local_weights[g] * (vxb * gden_b_xz[atom_g] + vyb * gden_b_yz[atom_g] + vzb * gden_b_zz[atom_g]);
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

    for (int iatom = 0; iatom < natoms; iatom++)
    {
        for (int thread_id = 0; thread_id < nthreads; thread_id++)
        {
            molgrad.row(iatom)[0] += molgrad_threads.row(thread_id)[iatom * 3 + 0];
            molgrad.row(iatom)[1] += molgrad_threads.row(thread_id)[iatom * 3 + 1];
            molgrad.row(iatom)[2] += molgrad_threads.row(thread_id)[iatom * 3 + 2];
        }
    }

    return molgrad;
}

void
CXCMolecularGradient::_computeAOtoAtomMapping(std::vector<int>& ao_to_atom_ids, const CMolecule& molecule, const CMolecularBasis& basis) const
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
