#include "XCIntegratorForKx.hpp"

#include "GtoFunc.hpp"
#include "OpenMPFunc.hpp"
#include "Prescreener.hpp"
#include "DftSubMatrix.hpp"
#include "MultiTimer.hpp"
#include "GtoValues.hpp"
#include "DenseLinearAlgebra.hpp"
#include "NuclearPotentialDriver.hpp"

#include <iostream>

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

    timer.start("Preparation");
    
    auto nthreads = omp_get_max_threads();
    
    std::vector<CMultiTimer> omptimers(nthreads);
    
    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);
    
    // non-symmetric exchange matrix

    CAOKohnShamMatrix mat_kx(naos, naos, std::string("closedshell"));

    mat_kx.zero();
    
    // GTOs on grid points

    auto max_npoints_per_box = molecularGrid.getMaxNumberOfGridPointsPerBox();

    auto omp_max_npoints = max_npoints_per_box / nthreads;
    
    if (max_npoints_per_box % nthreads != 0) omp_max_npoints++;
    
    // density and functional derivatives

    std::vector<std::vector<double>> omp_local_weights_data(nthreads, std::vector<double>(omp_max_npoints));
    
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

        std::vector<CGtoBlock> red_gto_blocks;
        
        for (const auto& gto_block : gto_blocks)
        {
            // 0th order GTO derivative
            auto [cgto_mask, pre_ao_inds] = prescr::preScreenGtoBlock(gto_block, 0, screeningThresholdForGTOValues, boxdim);

            auto red_gto_block = gto_block.reduce(cgto_mask);
            
            if (red_gto_block.number_of_basis_functions() > 0)
            {
                red_gto_blocks.push_back(red_gto_block);
            }
            
            std::cout << " * Mask : " << cgto_mask.size() << " : ";
            for (const auto mask : cgto_mask)
            {
                std::cout << mask << " ";
            }
            std::cout << std::endl;
            
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
        
        std::map<size_t, size_t> mask;
        
        for (size_t i = 0; i < aocount; i++)
        {
            mask.insert({(size_t)aoinds[i], i});
        }
        
        auto ptr_mask = &mask;
        
//        std::cout << "INDICES : ";
//        for (const auto ao_idx : aoinds)
//        {
//            std::cout << ao_idx << " ";
//        }
//        std::cout << std::endl;
//        
//        std::cout << "MASK : ";
//        for (const auto m : mask)
//        {
//            std::cout << "(" << m.first << "," << m.second << ")" ;
//        }
//        std::cout << std::endl;

        timer.start("Density matrix slicing");

        auto sub_dens_mat = dftsubmat::getSubDensityMatrix(gsDensityPointers[0], aoinds, naos);

        timer.stop("Density matrix slicing");
        
        timer.start("OMP Kx calc.");
        
        CDenseMatrix sum_partial_mat_kx(aocount, aocount);

        sum_partial_mat_kx.zero();
        
#pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();
            
            omptimers[thread_id].start("gtoeval");
            
            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);
            
            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);
            
            CDenseMatrix mat_chi(aocount, grid_batch_size);
            
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
                
                auto cmat = gtoval::get_gto_values_for_lda(gto_block, grid_x, grid_y, grid_z, cgto_mask);
                
                if (cmat.is_empty()) continue;
                
                auto submat_ptr = cmat.sub_matrix({0, 0});
                
                auto submat_data = submat_ptr->data();
                
                for (int nu = 0; nu < static_cast<int>(pre_ao_inds.size()); nu++, idx++)
                {
                    std::memcpy(mat_chi.row(idx), submat_data + nu * grid_batch_size, grid_batch_size * sizeof(double));
                }
            }
            
            omptimers[thread_id].stop("gtoeval");
            
            omptimers[thread_id].start("Generate F_vg matrix");

            auto mat_fvg = denblas::serialMultAB(sub_dens_mat, mat_chi);

            omptimers[thread_id].stop("Generate F_vg matrix ");
            
            omptimers[thread_id].start("Copy grid weights");
            
            auto local_weights = omp_local_weights_data[thread_id].data();

            std::memcpy(local_weights, weights + gridblockpos + grid_batch_offset, grid_batch_size * sizeof(double));

            omptimers[thread_id].stop("Copy grid weights");
            
            omptimers[thread_id].start("Generate G_gv matrix");
        
            CDenseMatrix mat_ggv(grid_batch_size, aocount);
            
            mat_ggv.zero();
            
            const auto npot_drv = CNuclearPotentialDriver();
            
            for (int g = 0; g < grid_batch_size; g++)
            {
                npot_drv.compute(mat_ggv, red_gto_blocks, mat_fvg, *ptr_mask, g, grid_x[g], grid_y[g], grid_z[g], local_weights[g]);
            }
            
            omptimers[thread_id].stop("Generate G_gv matrix");
            
            omptimers[thread_id].start("Kx matmul");

            auto partial_mat_kx = denblas::serialMultAB(mat_chi, mat_ggv);

            omptimers[thread_id].stop("Kx matmul");

            omptimers[thread_id].start("Kx local matrix dist.");

#pragma omp critical
            denblas::serialInPlaceAddAB(sum_partial_mat_kx, partial_mat_kx);

            omptimers[thread_id].stop("Kx local matrix dist.");
        }
        
        std::cout << "Grid block : " << box_id  << " Grid dim. : " << npoints << " Local AOs : " << aocount << std::endl;
        
        timer.stop("OMP Kx calc.");

        timer.start("Kx matrix dist.");

        dftsubmat::distributeSubMatrixToKohnSham(mat_kx, sum_partial_mat_kx, aoinds);

        timer.stop("Kx matrix dist.");
    
    }
    
    mat_kx.inPlaceSymmetrizeAndScale(factor); 
    
    timer.stop("Total timing");
    
    std::cout << "Timing of new integrator" << std::endl;
    std::cout << "------------------------" << std::endl;
    std::cout << timer.getSummary() << std::endl;
    std::cout << "OpenMP timing" << std::endl;
    for (int thread_id = 0; thread_id < nthreads; thread_id++)
    {
         std::cout << "Thread " << thread_id << std::endl;
         std::cout << omptimers[thread_id].getSummary() << std::endl;
    }
    
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
