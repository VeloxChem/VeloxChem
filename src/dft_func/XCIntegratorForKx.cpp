#include "XCIntegratorForKx.hpp"

#include "GtoFunc.hpp"
#include "OpenMPFunc.hpp"

namespace xcintkx {  // xcintkx namespace

auto
integrateKxFockForClosedShell(const CMolecule&                  molecule,
                              const CMolecularBasis&            basis,
                              const std::vector<const double*>& gsDensityPointers,
                              const CMolecularGrid&             molecularGrid,
                              const double                      screeningThresholdForGTOValues,
                              const double                      factor) -> CAOKohnShamMatrix
{
    auto nthreads = omp_get_max_threads();
    
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
    auto nthreads = omp_get_max_threads();
    
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
    
    return mat_kx;
}

}  // namespace xcintkx
