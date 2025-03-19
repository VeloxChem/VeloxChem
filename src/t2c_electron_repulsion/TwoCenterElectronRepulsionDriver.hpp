#ifndef TwoCenterElectronRepulsionDriver_hpp
#define TwoCenterElectronRepulsionDriver_hpp

#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

/// @brief Class CTwoCenterElectronRepulsionDriver provides methods for computing two-center electron repulsion integrals.
class CTwoCenterElectronRepulsionDriver
{
   public:
    /// @brief Creates an electron repulsion integrals driver.
    CTwoCenterElectronRepulsionDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The electron repulsion integrals driver to be copied.
    CTwoCenterElectronRepulsionDriver(const CTwoCenterElectronRepulsionDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The electron repulsion integrals driver  to be moved.
    CTwoCenterElectronRepulsionDriver(CTwoCenterElectronRepulsionDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CTwoCenterElectronRepulsionDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The electron repulsion integrals driver to be copy assigned.
    /// @return The assigned electron repulsion integrals driver.
    auto operator=(const CTwoCenterElectronRepulsionDriver &other) -> CTwoCenterElectronRepulsionDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The electron repulsion integrals driver to be move assigned.
    /// @return The assigned electron repulsion integrals driver .
    auto operator=(CTwoCenterElectronRepulsionDriver &&other) noexcept -> CTwoCenterElectronRepulsionDriver & = delete;

    /// @brief The equality operator.
    /// @param other The electron repulsion integrals driver  to be compared.
    /// @return True if electron repulsion integrals drivers  are equal, False otherwise.
    auto operator==(const CTwoCenterElectronRepulsionDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The electron repulsion integrals driver to be compared.
    /// @return True if electron repulsion integrals drivers  are not equal, False otherwise.
    auto operator!=(const CTwoCenterElectronRepulsionDriver &other) const -> bool = delete;

    /// @brief Computes electron repulsion matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @return The electron repulsion matrix.
    auto compute(const CMolecularBasis &basis, const CMolecule &molecule) const -> CMatrix;
};

#endif /* TwoCenterElectronRepulsionDriver_hpp */
