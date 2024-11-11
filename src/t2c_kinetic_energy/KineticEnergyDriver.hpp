#ifndef KineticEnergyDriver_hpp
#define KineticEnergyDriver_hpp

#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

/// @brief Class CKineticEnergyDriver provides methods for computing two-center kinetic energy integrals.
class CKineticEnergyDriver
{
   public:
    /// @brief Creates an kinetic energy integrals driver.
    CKineticEnergyDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The kinetic energy integrals driver to be copied.
    CKineticEnergyDriver(const CKineticEnergyDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The kinetic energy integrals driver  to be moved.
    CKineticEnergyDriver(CKineticEnergyDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CKineticEnergyDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The kinetic energy integrals driver to be copy assigned.
    /// @return The assigned kinetic energy integrals driver.
    auto operator=(const CKineticEnergyDriver &other) -> CKineticEnergyDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The kinetic energy integrals driver to be move assigned.
    /// @return The assigned kinetic energy integrals driver .
    auto operator=(CKineticEnergyDriver &&other) noexcept -> CKineticEnergyDriver & = delete;

    /// @brief The equality operator.
    /// @param other The kinetic energy integrals driver  to be compared.
    /// @return True if kinetic energy integrals drivers  are equal, False otherwise.
    auto operator==(const CKineticEnergyDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The kinetic energy integrals driver to be compared.
    /// @return True if kinetic energy integrals drivers  are not equal, False otherwise.
    auto operator!=(const CKineticEnergyDriver &other) const -> bool = delete;

    /// @brief Computes kinetic energy matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @return The kinetic energy matrix.
    auto compute(const CMolecularBasis &basis, const CMolecule &molecule) const -> CMatrix;
};

#endif /* KineticEnergyDriver_hpp */
