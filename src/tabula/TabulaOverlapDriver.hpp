//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center overlap driver.
//

#ifndef TabulaOverlapDriver_hpp
#define TabulaOverlapDriver_hpp

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "TabulaDenseMatrix.hpp"

namespace tabula {  // tabula namespace

/// @brief Driver for the two-center overlap integral.
///
/// Intakes a VeloxChem molecule and molecular basis, builds the basis-function
/// blocks with `gtofunc::make_gto_blocks`, and for each block pair constructs
/// a screened `CGtoPairBlock` (the overlap screening estimator drops the
/// negligible contracted-GTO pairs). The overlap matrix is then evaluated on a
/// late-contraction path — the primitive integrals of each surviving
/// contracted-GTO pair are summed, weighted by the contraction coefficients.
///
/// This is the scaffold; the per-shell-pair primitive recursion is a separate
/// seam (see `TabulaOverlapDriver.cpp`).
class OverlapDriver
{
   public:
    /// @brief Creates an overlap driver.
    OverlapDriver() = default;

    /// @brief Computes the overlap matrix.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param threshold The screening threshold — a contracted-GTO pair is
    /// kept when its overlap screening estimate is at or above this value.
    /// `0` (the default) keeps every pair.
    /// @return The overlap matrix.
    auto compute(const CMolecule& molecule, const CMolecularBasis& basis, const double threshold = 0.0) const -> DenseMatrix;
};

}  // namespace tabula

#endif /* TabulaOverlapDriver_hpp */
