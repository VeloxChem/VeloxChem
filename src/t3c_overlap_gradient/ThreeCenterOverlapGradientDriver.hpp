#ifndef ThreeCenterOverlapGradientDriver_hpp
#define ThreeCenterOverlapGradientDriver_hpp

#include <vector>

#include "Matrices.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

/// @brief Class CThreeCenterOverlapGradientDriver provides methods for computing three-center overlap integrals.
class CThreeCenterOverlapGradientDriver
{
   public:
    /// @brief Creates an overlap integrals driver.
    CThreeCenterOverlapGradientDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The overlap integrals driver to be copied.
    CThreeCenterOverlapGradientDriver(const CThreeCenterOverlapGradientDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The overlap integrals driver  to be moved.
    CThreeCenterOverlapGradientDriver(CThreeCenterOverlapGradientDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CThreeCenterOverlapGradientDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The overlap integrals driver to be copy assigned.
    /// @return The assigned overlap integrals driver.
    auto operator=(const CThreeCenterOverlapGradientDriver &other) -> CThreeCenterOverlapGradientDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The overlap integrals driver to be move assigned.
    /// @return The assigned overlap integrals driver .
    auto operator=(CThreeCenterOverlapGradientDriver &&other) noexcept -> CThreeCenterOverlapGradientDriver & = delete;

    /// @brief The equality operator.
    /// @param other The overlap integrals driver  to be compared.
    /// @return True if overlap integrals drivers  are equal, False otherwise.
    auto operator==(const CThreeCenterOverlapGradientDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The overlap integrals driver to be compared.
    /// @return True if overlap integrals drivers  are not equal, False otherwise.
    auto operator!=(const CThreeCenterOverlapGradientDriver &other) const -> bool = delete;

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
                 const CMolecule                   &molecule) const -> CMatrices;
};



#endif /* ThreeCenterOverlapGradientDriver_hpp */
