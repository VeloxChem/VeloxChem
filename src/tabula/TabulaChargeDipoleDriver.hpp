//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center charge-dipole driver.
//

#ifndef TabulaChargeDipoleDriver_hpp
#define TabulaChargeDipoleDriver_hpp

#include <array>
#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "TabulaBlockSparseMatrix.hpp"
#include "TabulaDenseMatrix.hpp"
#include "TabulaExternalCenterDriver.hpp"  // KernelProfile, ThreadBalance
#include "TabulaThreadBalance.hpp"

namespace tabula {  // tabula namespace

/// @brief Driver for the two-center charge-dipole integral
/// `Σ_N d_N · (a|(r−N)/|r−N|³|c)` — the field of the bra·ket density at each
/// point `N`, contracted with the site's dipole `d_N`, summed over sites. A
/// thin wrapper over the shared external-center core: it supplies the
/// charge-dipole kernel (the three field axes summed with the site moments) and
/// builds the point-dipole set. The dipoles are always external — a QM/MM
/// polarising environment — so there is no "from molecule" overload.
///
/// Separate from `NuclearAttractionDriver` (the two grow distinct special
/// methods), but both ride the same source-agnostic external-center scaffolding.
class ChargeDipoleDriver
{
   public:
    /// @brief Creates a charge-dipole driver.
    ChargeDipoleDriver() = default;

    /// @brief Computes the charge-dipole matrix over a set of point dipoles.
    /// @param molecule The molecule (supplies the basis-function centres).
    /// @param basis The molecular basis.
    /// @param moments The dipole vector `d_N` of each site, in atomic units.
    /// @param coordinates The site positions, in atomic units; one per moment.
    /// @param threshold The screening threshold. Negative (the default) selects
    /// automatically — a conservative screen for large molecules (≥ 100 atoms),
    /// exact otherwise; `0` forces exact dense; `> 0` is an explicit screen.
    /// @param profile Optional — receives the per-phase wall-time breakdown.
    auto compute(const CMolecule&                          molecule,
                 const CMolecularBasis&                    basis,
                 const std::vector<std::array<double, 3>>& moments,
                 const std::vector<std::array<double, 3>>& coordinates,
                 const double                              threshold = -1.0,
                 KernelProfile                            *profile   = nullptr) const -> DenseMatrix;

    /// @brief Computes the matrix in block-sparse storage, over the dipoles.
    auto computeSparse(const CMolecule&                          molecule,
                       const CMolecularBasis&                    basis,
                       const std::vector<std::array<double, 3>>& moments,
                       const std::vector<std::array<double, 3>>& coordinates,
                       const double                              threshold = 0.0,
                       KernelProfile                            *profile   = nullptr) const -> BlockSparseMatrix;
};

/// @brief Gets the per-thread balance captured by the most recent
/// `ChargeDipoleDriver` compute run.
auto charge_dipole_thread_balance() -> ThreadBalance;

}  // namespace tabula

#endif /* TabulaChargeDipoleDriver_hpp */
