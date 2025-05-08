#ifndef ThreeCenterR2Driver_hpp
#define ThreeCenterR2Driver_hpp

#include <vector>

#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

/// @brief Class CThreeCenterOverlapDriver provides methods for computing three-center overlap integrals.
class CThreeCenterR2Driver
{
   public:
    /// @brief Creates an overlap integrals driver.
    CThreeCenterR2Driver() = default;

    /// @brief The default copy constructor.
    /// @param other The overlap integrals driver to be copied.
    CThreeCenterR2Driver(const CThreeCenterR2Driver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The overlap integrals driver  to be moved.
    CThreeCenterR2Driver(CThreeCenterR2Driver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CThreeCenterR2Driver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The overlap integrals driver to be copy assigned.
    /// @return The assigned overlap integrals driver.
    auto operator=(const CThreeCenterR2Driver &other) -> CThreeCenterR2Driver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The overlap integrals driver to be move assigned.
    /// @return The assigned overlap integrals driver .
    auto operator=(CThreeCenterR2Driver &&other) noexcept -> CThreeCenterR2Driver & = delete;

    /// @brief The equality operator.
    /// @param other The overlap integrals driver  to be compared.
    /// @return True if overlap integrals drivers  are equal, False otherwise.
    auto operator==(const CThreeCenterR2Driver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The overlap integrals driver to be compared.
    /// @return True if overlap integrals drivers  are not equal, False otherwise.
    auto operator!=(const CThreeCenterR2Driver &other) const -> bool = delete;

    /// @brief Computes overlap matrix for given molecule and molecular basis.
    /// @param exponents The vector of Gausian exponents of external centers.
    /// @param factors The vector of scaling factors of Gaussian exponents of external centers.
    /// @param coordinates The vector of external center coordinates.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @return The overlap matrix.
    auto compute(const std::vector<double>         &exponents,
                 const std::vector<double>         &factors,
                 const std::vector<TPoint<double>> &coordinates,
                 const CMolecularBasis             &basis,
                 const CMolecule                   &molecule) const -> CMatrix;
};


#endif /* ThreeCenterOverlapDriver_hpp */
