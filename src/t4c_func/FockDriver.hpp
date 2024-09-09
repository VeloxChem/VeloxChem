#ifndef FockDriver_hpp
#define FockDriver_hpp

#include <string>

#include "Matrices.hpp"
#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

/// Class CFockDriver provides methods for computing Fock matrices
/// using four center electron repulsion integrals..
class CFockDriver
{
   public:
    /// Creates a Fock matrices  driver.
    CFockDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The Fock matrices driver to be copied.
    CFockDriver(const CFockDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The Fock matrices driver  to be moved.
    CFockDriver(CFockDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CFockDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The Fock matrices driver to be copy assigned.
    /// @return The assigned Fock matrices driver.
    auto operator=(const CFockDriver &other) -> CFockDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The Fock matrices driver to be move assigned.
    /// @return The assigned Fock matrices driver .
    auto operator=(CFockDriver &&other) noexcept -> CFockDriver & = delete;

    /// @brief The equality operator.
    /// @param other The Fock matrices driver  to be compared.
    /// @return True if Fock matrices drivers  are equal, False otherwise.
    auto operator==(const CFockDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The Fock matrices driver to be compared.
    /// @return True if Fock matrices drivers  are not equal, False otherwise.
    auto operator!=(const CFockDriver &other) const -> bool = delete;

    /// @brief Computes Fock matrix for given density, basis and molecule (N^4 scaling).
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param density The density matrix to construct Fock matrix.
    /// @param label The label of Fock matrix type.
    /// @param exchange_factor The exchange-correlation factors.
    /// @param omega The range separation factor.
    /// @return The Fock matrix.
    auto compute(const CMolecularBasis &basis,
                 const CMolecule       &molecule,
                 const CMatrix         &density,
                 const std::string     &label,
                 const double           exchange_factor,
                 const double           omega) const -> CMatrix;
};

#endif /* FockDriver_hpp */
