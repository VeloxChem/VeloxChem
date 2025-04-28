#include "XCIntegratorForKx.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <iostream>

#include "GtoFunc.hpp"
#include "OpenMPFunc.hpp"
#include "Prescreener.hpp"
#include "DftSubMatrix.hpp"
#include "MultiTimer.hpp"
#include "GtoValues.hpp"
#include "DenseLinearAlgebra.hpp"
#include "NuclearPotentialDriver.hpp"

//#include <iostream>

namespace xcintkx {  // xcintkx namespace

auto
integrateKxFockForClosedShell(const CMolecule&                  molecule,
                              const CMolecularBasis&            basis,
                              const std::vector<const double*>& gsDensityPointers,
                              const CMolecularGrid&             molecularGrid,
                              const double                      screeningThresholdForGTOValues,
                              const double                      factor) -> CAOKohnShamMatrix
{
    CMultiTimer timer;
    
    timer.start("Total timing");
    
    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // Kohn-Sham matrix

    CAOKohnShamMatrix mat_kx(naos, naos, std::string("closedshell"));

    mat_kx.zero();

    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    
    auto ycoords = molecularGrid.getCoordinatesY();
    
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();
    
    // set up number of grid blocks
    
    const auto nboxes = counts.size();
    
    const auto ngblocks = gto_blocks.size();
    
    // set up pointers to OMP data
    
    auto ptr_counts = counts.data();
    
    auto ptr_displacements = displacements.data();
    
    auto ptr_gto_blocks = gto_blocks.data();
    
    auto ptr_gsDensityPointers = gsDensityPointers.data();
    
#pragma omp parallel shared(ptr_counts, ptr_displacements, xcoords, ycoords, zcoords, ptr_gto_blocks, ptr_gsDensityPointers, nboxes, ngblocks, naos, mat_kx)
    {
#pragma omp single nowait
        {
            for (size_t box_id = 0; box_id < nboxes; box_id++)
            {
#pragma omp task firstprivate(box_id)
                {
                    // grid points in box
                    
                    auto npoints = ptr_counts[box_id];
                    
                    auto gridblockpos = ptr_displacements[box_id];
                    
                    // dimension of grid box
                    
                    auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);
                    
                    // prescreening
                    
                    std::vector<std::vector<int>> cgto_mask_blocks, pre_ao_inds_blocks;
                    
                    std::vector<int> aoinds;
                    
                    // reserve prescreening vectors
                    
                    cgto_mask_blocks.reserve(ngblocks);
                    
                    pre_ao_inds_blocks.reserve(ngblocks);
                    
                    aoinds.reserve(naos);
                    
                    // set up reduced GTO blocks
                    
                    std::vector<CGtoBlock> red_gto_blocks;
                    
                    red_gto_blocks.reserve(ngblocks);
                    
                    for (size_t i = 0; i < ngblocks; i++)
                    {
                        // 0th order GTO derivative
                        auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(ptr_gto_blocks[i], 0, screeningThresholdForGTOValues, boxdim);
                        
                        auto red_gto_block = ptr_gto_blocks[i].reduce(cgto_mask);
                        
                        if (red_gto_block.number_of_basis_functions() > 0)
                        {
                            red_gto_blocks.push_back(red_gto_block);
                        }
                        
                        cgto_mask_blocks.push_back(cgto_mask);
                        
                        pre_ao_inds_blocks.push_back(pre_ao_inds);
                        
                        for (const auto nu : pre_ao_inds)
                        {
                            aoinds.push_back(nu);
                        }
                    }
                    
                    const auto aocount = static_cast<int>(aoinds.size());
                    
                    // compute VXC contributions
                    
                    if (aocount > 0)
                    {
                        std::map<size_t, size_t> mask;
                    
                        for (size_t i = 0; i < aocount; i++)
                        {
                            mask.insert({(size_t)aoinds[i], i});
                        }
                        
                        // density matrix slicing
                        
                        auto sub_dens_mat = dftsubmat::getSubDensityMatrix(ptr_gsDensityPointers[0], aoinds, naos);
                        
                        // GTO values on grid points
                        
                        CDenseMatrix mat_chi(aocount, npoints);
                        
                        const auto grid_x_ptr = xcoords + gridblockpos;
                        
                        const auto grid_y_ptr = ycoords + gridblockpos;
                        
                        const auto grid_z_ptr = zcoords + gridblockpos;
                        
                        std::vector<double> grid_x(grid_x_ptr, grid_x_ptr + npoints);
                        
                        std::vector<double> grid_y(grid_y_ptr, grid_y_ptr + npoints);
                        
                        std::vector<double> grid_z(grid_z_ptr, grid_z_ptr + npoints);
                        
                        // go through GTO blocks

                        for (size_t i_block = 0, idx = 0; i_block < ngblocks; i_block++)
                        {
                            const auto& gto_block = ptr_gto_blocks[i_block];

                            const auto& cgto_mask = cgto_mask_blocks[i_block];

                            const auto& pre_ao_inds = pre_ao_inds_blocks[i_block];

                            auto cmat = gtoval::get_gto_values_for_lda(gto_block, grid_x, grid_y, grid_z, cgto_mask);

                            if (cmat.is_empty()) continue;

                            auto submat_ptr = cmat.sub_matrix({0, 0});

                            auto submat_data = submat_ptr->data();

                            for (int nu = 0; nu < static_cast<int>(pre_ao_inds.size()); nu++, idx++)
                            {
                                std::memcpy(mat_chi.row(idx), submat_data + nu * npoints, npoints * sizeof(double));
                            }
                        }
                        
                        // compute F_vg matrix
                        
                        auto mat_fvg = denblas::serialMultAB(sub_dens_mat, mat_chi);
                        
                        // copy weights to local vector
                        
                        std::vector<double> local_weights(weights + gridblockpos, weights + gridblockpos + npoints);
                        
                        // compute G_gv matrix
                        
//                        CDenseMatrix mat_ggv(npoints, aocount);
//                        
//                        mat_ggv.zero();
                        
                        const auto npot_drv = CNuclearPotentialDriver();
                        
                        auto mat_ggv = npot_drv.compute(red_gto_blocks, mat_fvg, mask, grid_x, grid_y, grid_z, local_weights);
                        
//                        for (int g = 0; g < npoints; g++)
//                        {
//                            npot_drv.compute(mat_ggv, red_gto_blocks, mat_fvg, mask, g, grid_x[g], grid_y[g], grid_z[g], local_weights[g]);
//                        }
                        
                        // compute local Kx matrix
                        
                        auto local_mat_kx = denblas::serialMultAB(mat_chi, mat_ggv);
                        
                        // distribute local Kx matrix
                        
                        #pragma omp critical
                        dftsubmat::distributeSubMatrixToKohnSham(mat_kx, local_mat_kx, aoinds);
                    }
                }
            }
        }
    }
    
    mat_kx.inPlaceSymmetrizeAndScale(factor);

    timer.stop("Total timing");
   
//    std::cout << "Timing of new integrator" << std::endl;
//    std::cout << "------------------------" << std::endl;
//    std::cout << timer.getSummary() << std::endl;
    
    return mat_kx;
}

auto
integrateKxFockForOpenShell(const CMolecule&                  molecule,
                            const CMolecularBasis&            basis,
                            const std::vector<const double*>& gsDensityPointers,
                            const CMolecularGrid&             molecularGrid,
                            const double                      screeningThresholdForGTOValues,
                            const double                      factor) -> CAOKohnShamMatrix
{
    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");
    
    auto nthreads = omp_get_max_threads();
    
    std::vector<CMultiTimer> omptimers(nthreads);
    
    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);
    
    // non-symmetric exchange matrix

    CAOKohnShamMatrix mat_kx(naos, naos, std::string("openshell"));

    mat_kx.zero();
    
    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    auto omp_max_npoints = max_npoints_per_box / nthreads;
    
    if (max_npoints_per_box % nthreads != 0) omp_max_npoints++;
    
    // coordinates and weights of grid points

    auto xcoords = molecularGrid.getCoordinatesX();
    auto ycoords = molecularGrid.getCoordinatesY();
    auto zcoords = molecularGrid.getCoordinatesZ();

    auto weights = molecularGrid.getWeights();

    // counts and displacements of grid points in boxes

    auto counts = molecularGrid.getGridPointCounts();

    auto displacements = molecularGrid.getGridPointDisplacements();

    timer.stop("Preparation");

    for (size_t box_id = 0; box_id < counts.size(); box_id++)
    {
        // grid points in box

        auto npoints = counts.data()[box_id];

        auto gridblockpos = displacements.data()[box_id];

        // dimension of grid box

        auto boxdim = prescr::getGridBoxDimension(gridblockpos, npoints, xcoords, ycoords, zcoords);

        // prescreening

        timer.start("GTO pre-screening");

        std::vector<std::vector<int>> cgto_mask_blocks, pre_ao_inds_blocks;

        std::vector<int> aoinds;

        for (const auto& gto_block : gto_blocks)
        {
            // 0th order GTO derivative
            auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(gto_block, 0, screeningThresholdForGTOValues, boxdim);

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

        timer.start("Density matrix slicing");

        auto sub_dens_mat_a = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);
        auto sub_dens_mat_b = dftsubmat::getSubDensityMatrix(gsDensityPointers[1], aoinds, naos);

        timer.stop("Density matrix slicing");
    
    }

    timer.stop("Total timing");
    
    return mat_kx;
}

}  // namespace xcintkx
