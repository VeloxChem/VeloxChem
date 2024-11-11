#ifndef ElectricDipoleMomentumDriver_hpp
#define ElectricDipoleMomentumDriver_hpp

#include "Matrices.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

/// @brief Class CElectricDipoleMomentumDriver provides methods for computing two-center electric dipole momentum integrals.
class CElectricDipoleMomentumDriver
{
   public:
    /// @brief Creates an electric dipole momentum integrals driver.
    CElectricDipoleMomentumDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The electric dipole momentum integrals driver to be copied.
    CElectricDipoleMomentumDriver(const CElectricDipoleMomentumDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The electric dipole momentum integrals driver  to be moved.
    CElectricDipoleMomentumDriver(CElectricDipoleMomentumDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CElectricDipoleMomentumDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The electric dipole momentum integrals driver to be copy assigned.
    /// @return The assigned electric dipole momentum integrals driver.
    auto operator=(const CElectricDipoleMomentumDriver &other) -> CElectricDipoleMomentumDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The electric dipole momentum integrals driver to be move assigned.
    /// @return The assigned electric dipole momentum integrals driver .
    auto operator=(CElectricDipoleMomentumDriver &&other) noexcept -> CElectricDipoleMomentumDriver & = delete;

    /// @brief The equality operator.
    /// @param other The electric dipole momentum integrals driver  to be compared.
    /// @return True if electric dipole momentum integrals drivers  are equal, False otherwise.
    auto operator==(const CElectricDipoleMomentumDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The electric dipole momentum integrals driver to be compared.
    /// @return True if electric dipole momentum integrals drivers  are not equal, False otherwise.
    auto operator!=(const CElectricDipoleMomentumDriver &other) const -> bool = delete;

    /// @brief Computes electric dipole momentum matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param origin The origin of electric dipole momentum.
    /// @return The electric dipole momentum matrix.
    auto compute(const CMolecularBasis &basis, const CMolecule &molecule, const TPoint<double> &origin) const -> CMatrices;
};

#endif /* ElectricDipoleMomentumDriver_hpp */
