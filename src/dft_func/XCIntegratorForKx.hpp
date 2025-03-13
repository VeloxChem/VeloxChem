#ifndef XCIntegratorForKx_hpp
#define XCIntegratorForKx_hpp

#include <string>

#include "AODensityMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"
#include "DensityGrid.hpp"
#include "DensityGridCubic.hpp"
#include "DensityGridQuad.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MultiTimer.hpp"
#include "XCFunctional.hpp"

namespace xcintkx {  // xcintkx namespace

/**
 Integrates  exchange contribution to closed-shell AO Kohn-Sham matrix.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param gsDensityPointers the pointers to AO density matrices.
 @param molecularGrid the molecular grid.
 @param factor the exchange scaling factor.
 @return the AO Kohn-Sham matrix.
 */
auto integrateKxFockForClosedShell(const CMolecule&                  molecule,
                                   const CMolecularBasis&            basis,
                                   const std::vector<const double*>& gsDensityPointers,
                                   const CMolecularGrid&             molecularGrid,
                                   const double                      screeningThresholdForGTOValues,
                                   const double                      factor) -> CAOKohnShamMatrix;

/**
 Integrates exchange contribution to open-shell AO Kohn-Sham matrix.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param gsDensityPointers the pointers to AO density matrices.
 @param molecularGrid the molecular grid.
 @param factor the exchange scaling factor.
 @return the AO Kohn-Sham matrix.
 */
auto integrateKxFockForOpenShell(const CMolecule&                  molecule,
                                 const CMolecularBasis&            basis,
                                 const std::vector<const double*>& gsDensityPointers,
                                 const CMolecularGrid&             molecularGrid,
                                 const double                      screeningThresholdForGTOValues,
                                 const double                      factor) -> CAOKohnShamMatrix;

}  // namespace xcintkx


#endif /* XCIntegratorForKx_hpp */
