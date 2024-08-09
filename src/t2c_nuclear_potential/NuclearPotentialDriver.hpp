#ifndef NuclearPotentialDriver_hpp
#define NuclearPotentialDriver_hpp

#include <vector>

#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

/// @brief Class CNuclearPotentialDriver provides methods for computing two-center nuclear potential integrals.
class CNuclearPotentialDriver
{
   public:
    /// @brief Creates an nuclear potential integrals driver.
    CNuclearPotentialDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The nuclear potential integrals driver to be copied.
    CNuclearPotentialDriver(const CNuclearPotentialDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The nuclear potential integrals driver  to be moved.
    CNuclearPotentialDriver(CNuclearPotentialDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CNuclearPotentialDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The nuclear potential integrals driver to be copy assigned.
    /// @return The assigned nuclear potential integrals driver.
    auto operator=(const CNuclearPotentialDriver &other) -> CNuclearPotentialDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The nuclear potential integrals driver to be move assigned.
    /// @return The assigned nuclear potential integrals driver .
    auto operator=(CNuclearPotentialDriver &&other) noexcept -> CNuclearPotentialDriver & = delete;

    /// @brief The equality operator.
    /// @param other The nuclear potential integrals driver  to be compared.
    /// @return True if nuclear potential integrals drivers  are equal, False otherwise.
    auto operator==(const CNuclearPotentialDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The nuclear potential integrals driver to be compared.
    /// @return True if nuclear potential integrals drivers  are not equal, False otherwise.
    auto operator!=(const CNuclearPotentialDriver &other) const -> bool = delete;

    /// @brief Computes nuclear potential matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @return The nuclear potential matrix.
    auto compute(const CMolecularBasis &basis, const CMolecule &molecule) const -> CMatrix;

    /// @brief Computes nuclear potential matrix for given set of external charges,  molecule and molecular basis.
    /// @param charges The vector of external charges.
    /// @param coordinates The vector of coordinates of external charges.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @return The nuclear potential matrix.
    auto compute(const std::vector<double>         &charges,
                 const std::vector<TPoint<double>> &coordinates,
                 const CMolecularBasis             &basis,
                 const CMolecule                   &molecule) const -> CMatrix;
};

#endif /* NuclearPotentialDriver_hpp */
