#ifndef NuclearPotentialErfDriver_hpp
#define NuclearPotentialErfDriver_hpp

#include <vector>

#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

/// @brief Class CNuclearPotentialErfDriver provides methods for computing two-center range separated nuclear potential integrals.
class CNuclearPotentialErfDriver
{
   public:
    /// @brief Creates an range separated nuclear potential integrals driver.
    CNuclearPotentialErfDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The range separated nuclear potential integrals driver to be copied.
    CNuclearPotentialErfDriver(const CNuclearPotentialErfDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The range separated nuclear potential integrals driver  to be moved.
    CNuclearPotentialErfDriver(CNuclearPotentialErfDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CNuclearPotentialErfDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The range separated nuclear potential integrals driver to be copy assigned.
    /// @return The assigned range separated nuclear potential integrals driver.
    auto operator=(const CNuclearPotentialErfDriver &other) -> CNuclearPotentialErfDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The range separated nuclear potential integrals driver to be move assigned.
    /// @return The assigned range separated nuclear potential integrals driver .
    auto operator=(CNuclearPotentialErfDriver &&other) noexcept -> CNuclearPotentialErfDriver & = delete;

    /// @brief The equality operator.
    /// @param other The range separated nuclear potential integrals driver  to be compared.
    /// @return True if range separated nuclear potential integrals drivers  are equal, False otherwise.
    auto operator==(const CNuclearPotentialErfDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The range separated nuclear potential integrals driver to be compared.
    /// @return True if range separated nuclear potential integrals drivers  are not equal, False otherwise.
    auto operator!=(const CNuclearPotentialErfDriver &other) const -> bool = delete;

    /// @brief Computes range separated nuclear potential matrix for given set of external charges,  molecule and molecular basis.
    /// @param charges The vector of external charges.
    /// @param coordinates The vector of coordinates of external charges.
    /// @param omegas The vector of range-separation factors.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @return The range separated nuclear potential matrix.
    auto compute(const std::vector<double>         &charges,
                 const std::vector<TPoint<double>> &coordinates,
                 const std::vector<double>         &omegas,
                 const CMolecularBasis             &basis,
                 const CMolecule                   &molecule) const -> CMatrix;

    /// @brief Computes range separated nuclear potential matrix for given set of external charges,  molecule and molecular basis.
    /// @param charges The vector of external charges.
    /// @param coordinates The vector of coordinates of external charges.
    /// @param omega The uniform range-separation factor.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @return The range separated nuclear potential matrix.
    auto compute(const std::vector<double>         &charges,
                 const std::vector<TPoint<double>> &coordinates,
                 const double                       omega,
                 const CMolecularBasis             &basis,
                 const CMolecule                   &molecule) const -> CMatrix;
};

#endif /* NuclearPotentialErfDriver_hpp */
