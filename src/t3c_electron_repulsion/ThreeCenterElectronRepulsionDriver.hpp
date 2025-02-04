#ifndef ThreeCenterElectronRepulsionDriver_hpp
#define ThreeCenterElectronRepulsionDriver_hpp

#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "T3FlatBuffer.hpp"

/// @brief Class CThreeCenterElectronRepulsionDriver provides methods for computing three-center electron repulsion integrals.
class CThreeCenterElectronRepulsionDriver
{
   public:
    /// @brief Creates an electron repulsion integrals driver.
    CThreeCenterElectronRepulsionDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The electron repulsion integrals driver to be copied.
    CThreeCenterElectronRepulsionDriver(const CThreeCenterElectronRepulsionDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The electron repulsion integrals driver  to be moved.
    CThreeCenterElectronRepulsionDriver(CThreeCenterElectronRepulsionDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CThreeCenterElectronRepulsionDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The electron repulsion integrals driver to be copy assigned.
    /// @return The assigned electron repulsion integrals driver.
    auto operator=(const CThreeCenterElectronRepulsionDriver &other) -> CThreeCenterElectronRepulsionDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The electron repulsion integrals driver to be move assigned.
    /// @return The assigned electron repulsion integrals driver .
    auto operator=(CThreeCenterElectronRepulsionDriver &&other) noexcept -> CThreeCenterElectronRepulsionDriver & = delete;

    /// @brief The equality operator.
    /// @param other The electron repulsion integrals driver  to be compared.
    /// @return True if electron repulsion integrals drivers  are equal, False otherwise.
    auto operator==(const CThreeCenterElectronRepulsionDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The electron repulsion integrals driver to be compared.
    /// @return True if electron repulsion integrals drivers  are not equal, False otherwise.
    auto operator!=(const CThreeCenterElectronRepulsionDriver &other) const -> bool = delete;

    /// @brief Computes electron repulsion matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxilary molecular basis for fiting of four-center repulsion integrals.
    /// @param molecule The molecule.
    /// @return The electron repulsion matrix.
    auto compute(const CMolecularBasis &basis, const CMolecularBasis &aux_basis, const CMolecule &molecule) const -> CT3FlatBuffer<double>;
};


#endif /* ThreeCenterElectronRepulsionDriver_hpp */
