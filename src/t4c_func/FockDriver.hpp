#ifndef FockDriver_hpp
#define FockDriver_hpp

#include <string>

#include "Matrices.hpp"
#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "T4CScreener.hpp"

/// Class CFockDriver provides methods for computing Fock matrices
/// using four center electron repulsion integrals..
class CFockDriver
{
   public:
    /// Creates a Fock matrices  driver.
    CFockDriver() = default;

    /// Computes Fock matrix for given density, basis and molecule.
    /// - Parameter basis: the molecular basis.
    /// - Parameter molecule: the molecule.
    /// - Parameter density: the density matrix to construct Fock matrix.
    /// - Parameter label: the label of Fock matrix type.
    /// - Parameter exchange_factor: the exchange-correlation factors.
    auto compute(const CMolecularBasis& basis,
                 const CMolecule&       molecule,
                 const CMatrix&         density,
                 const std::string&     label,
                 const double           exchange_factor) const -> CMatrix;

    /// Computes Fock matrix for given density, basis and molecule.
    /// - Parameter screener: the screener with blocked GTOs pairs data.
    /// - Parameter density: the density matrix to construct Fock matrix.
    /// - Parameter label: the label of Fock matrix type.
    /// - Parameter exchange_factor : the exchange-correlation factors.
    /// - Parameter ithreshold : the threshold for computing two-electron integrals.
    auto compute(const CT4CScreener& screener,
                 const CMatrix&      density,
                 const std::string&  label,
                 const double        exchange_factor,
                 const int           ithreshold) const -> CMatrix;

    /// Computes Fock matrix for given density, basis and molecule.
    /// - Parameter screener: the screener with blocked GTOs pairs data.
    /// - Parameter density: the density matrix to construct Fock matrix.
    /// - Parameter label: the label of Fock matrix type.
    /// - Parameter exchange_factor : the exchange-correlation factors.
    /// - Parameter ithreshold : the threshold for computing two-electron integrals.
    auto compute(const CT4CScreener&             screener,
                 const CMatrices&                densities,
                 const std::vector<std::string>& labels,
                 const std::vector<double>&      exchange_factors,
                 const int                       ithreshold) const -> CMatrices;

    /// Computes Fock matrix for given density, basis and molecule using ordered tasks grouping.
    /// - Parameter screener: the screener with blocked GTOs pairs data.
    /// - Parameter density: the density matrix to construct Fock matrix.
    /// - Parameter label: the label of Fock matrix type.
    /// - Parameter exchange_factor : the exchange-correlation factors.
    /// - Parameter ithreshold : the threshold for computing two-electron integrals.
    auto ordered_compute(const CT4CScreener& screener,
                         const CMatrix&      density,
                         const std::string&  label,
                         const double        exchange_factor,
                         const int           ithreshold) const -> CMatrix;
};

#endif /* FockDriver_hpp */
