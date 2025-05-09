#ifndef ThreeCenterRR2Driver_hpp
#define ThreeCenterRR2Driver_hpp

#include <vector>

#include "Matrices.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

/// @brief Class CThreeCenterOverlapGradientDriver provides methods for computing three-center overlap integrals.
class CThreeCenterRR2Driver
{
   public:
    /// @brief Creates an overlap integrals driver.
    CThreeCenterRR2Driver() = default;

    /// @brief The default copy constructor.
    /// @param other The overlap integrals driver to be copied.
    CThreeCenterRR2Driver(const CThreeCenterRR2Driver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The overlap integrals driver  to be moved.
    CThreeCenterRR2Driver(CThreeCenterRR2Driver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CThreeCenterRR2Driver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The overlap integrals driver to be copy assigned.
    /// @return The assigned overlap integrals driver.
    auto operator=(const CThreeCenterRR2Driver &other) -> CThreeCenterRR2Driver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The overlap integrals driver to be move assigned.
    /// @return The assigned overlap integrals driver .
    auto operator=(CThreeCenterRR2Driver &&other) noexcept -> CThreeCenterRR2Driver & = delete;

    /// @brief The equality operator.
    /// @param other The overlap integrals driver  to be compared.
    /// @return True if overlap integrals drivers  are equal, False otherwise.
    auto operator==(const CThreeCenterRR2Driver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The overlap integrals driver to be compared.
    /// @return True if overlap integrals drivers  are not equal, False otherwise.
    auto operator!=(const CThreeCenterRR2Driver &other) const -> bool = delete;

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



#endif /* ThreeCenterRR2Driver_hpp */
