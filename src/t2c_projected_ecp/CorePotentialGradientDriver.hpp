#ifndef CorePotentialGradientDriver_hpp
#define CorePotentialGradientDriver_hpp

#include <vector>

#include "Matrices.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "AtomCorePotential.hpp"

/// @brief Class CCorePotentialDriver provides methods for computing two-center projected ECP integrals.
class CCorePotentialGradientDriver
{
   public:
    /// @brief Creates a ECP integrals driver.
    CCorePotentialGradientDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The ECP integrals driver to be copied.
    CCorePotentialGradientDriver(const CCorePotentialGradientDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The ECP integrals driver  to be moved.
    CCorePotentialGradientDriver(CCorePotentialGradientDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CCorePotentialGradientDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The ECP integrals driver to be copy assigned.
    /// @return The assigned ECP integrals driver.
    auto operator=(const CCorePotentialGradientDriver &other) -> CCorePotentialGradientDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The ECP integrals driver to be move assigned.
    /// @return The assigned ECP integrals driver .
    auto operator=(CCorePotentialGradientDriver &&other) noexcept -> CCorePotentialGradientDriver & = delete;

    /// @brief The equality operator.
    /// @param other The ECP integrals driver  to be compared.
    /// @return True if ECP integrals drivers  are equal, False otherwise.
    auto operator==(const CCorePotentialGradientDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The ECP integrals driver to be compared.
    /// @return True if ECP integrals drivers  are not equal, False otherwise.
    auto operator!=(const CCorePotentialGradientDriver &other) const -> bool = delete;

    /// @brief Computes ECP matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param atom_potential The atom core potential.
    /// @return The ECP matrix.
    auto compute_bra_grad(const CMolecularBasis &basis, const CMolecule &molecule, const CAtomCorePotential& atom_potential, const int iatom) const -> CMatrices;
    
    /// @brief Computes ECP matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param atom_potential The atom core potential.
    /// @return The ECP matrix.
    auto compute_bra_grad(const CMolecularBasis &basis, const CMolecule &molecule, const std::vector<int>& atoms, const int iatom) const -> CMatrices;
    
    /// @brief Computes ECP matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param atom_potential The atom core potential.
    /// @return The ECP matrix.
    auto compute_pot_grad(const CMolecularBasis &basis, const CMolecule &molecule, const CAtomCorePotential& atom_potential) const -> CMatrices;
    
    /// @brief Computes ECP matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param atom_potential The atom core potential.
    /// @return The ECP matrix.
    auto compute_pot_grad(const CMolecularBasis &basis, const CMolecule &molecule, const int iatom) const -> CMatrices;
};

#endif /* CorePotentialGradientDriver_hpp */
