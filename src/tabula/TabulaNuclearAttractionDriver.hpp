//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center nuclear-attraction driver.
//

#ifndef TabulaNuclearAttractionDriver_hpp
#define TabulaNuclearAttractionDriver_hpp

#include <array>
#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "TabulaBlockSparseMatrix.hpp"
#include "TabulaDenseMatrix.hpp"
#include "TabulaExternalCenterDriver.hpp"  // KernelProfile, ThreadBalance
#include "TabulaThreadBalance.hpp"

namespace tabula {  // tabula namespace

/// @brief Driver for the two-center nuclear-attraction integral
/// `Σ_N Z_N (a|1/|r−N||c)` (VeloxChem's positive convention). A thin wrapper
/// over the shared external-center core — supplies the nuclear-attraction
/// kernel and screening estimate, and builds the point-charge set from the
/// molecule's nuclei or from explicit external charges (a QM/MM embedding).
class NuclearAttractionDriver
{
   public:
    /// @brief Creates a nuclear-attraction driver.
    NuclearAttractionDriver() = default;

    /// @brief Computes the nuclear-attraction matrix over the molecule's nuclei.
    auto compute(const CMolecule&       molecule,
                 const CMolecularBasis& basis,
                 const double           threshold = 0.0,
                 KernelProfile         *profile   = nullptr) const -> DenseMatrix;

    /// @brief Computes the matrix in block-sparse storage, over the nuclei.
    auto computeSparse(const CMolecule&       molecule,
                       const CMolecularBasis& basis,
                       const double           threshold = 0.0,
                       KernelProfile         *profile   = nullptr) const -> BlockSparseMatrix;

    /// @brief Computes the matrix over an explicit set of external point
    /// charges — the QM/MM use. `coordinates` are in atomic units.
    auto compute(const CMolecule&                       molecule,
                 const CMolecularBasis&                 basis,
                 const std::vector<double>&             magnitudes,
                 const std::vector<std::array<double, 3>>& coordinates,
                 const double                           threshold = 0.0,
                 KernelProfile                         *profile   = nullptr) const -> DenseMatrix;

    /// @brief Block-sparse variant over external point charges.
    auto computeSparse(const CMolecule&                       molecule,
                       const CMolecularBasis&                 basis,
                       const std::vector<double>&             magnitudes,
                       const std::vector<std::array<double, 3>>& coordinates,
                       const double                           threshold = 0.0,
                       KernelProfile                         *profile   = nullptr) const -> BlockSparseMatrix;
};

/// @brief Gets the per-thread balance captured by the most recent
/// `NuclearAttractionDriver` compute run.
auto nuclear_attraction_thread_balance() -> ThreadBalance;

}  // namespace tabula

#endif /* TabulaNuclearAttractionDriver_hpp */
