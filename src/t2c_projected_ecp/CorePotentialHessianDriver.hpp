#ifndef CorePotentialHessianDriver_hpp
#define CorePotentialHessianDriver_hpp

#include <vector>

#include "Matrices.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "AtomCorePotential.hpp"

/// @brief Class CCorePotentialDriver provides methods for computing two-center projected ECP integrals.
class CCorePotentialHessianDriver
{
   public:
    /// @brief Creates a ECP integrals driver.
    CCorePotentialHessianDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The ECP integrals driver to be copied.
    CCorePotentialHessianDriver(const CCorePotentialHessianDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The ECP integrals driver  to be moved.
    CCorePotentialHessianDriver(CCorePotentialHessianDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CCorePotentialHessianDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The ECP integrals driver to be copy assigned.
    /// @return The assigned ECP integrals driver.
    auto operator=(const CCorePotentialHessianDriver &other) -> CCorePotentialHessianDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The ECP integrals driver to be move assigned.
    /// @return The assigned ECP integrals driver .
    auto operator=(CCorePotentialHessianDriver &&other) noexcept -> CCorePotentialHessianDriver & = delete;

    /// @brief The equality operator.
    /// @param other The ECP integrals driver  to be compared.
    /// @return True if ECP integrals drivers  are equal, False otherwise.
    auto operator==(const CCorePotentialHessianDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The ECP integrals driver to be compared.
    /// @return True if ECP integrals drivers  are not equal, False otherwise.
    auto operator!=(const CCorePotentialHessianDriver &other) const -> bool = delete;

    /// @brief Computes ECP matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param atom_potential The atom core potential.
    /// @return The ECP matrix.
    auto compute_geom_200(const CMolecularBasis &basis, const CMolecule &molecule, const CAtomCorePotential& atom_potential, const int iatom) const -> CMatrices;
    
    /// @brief Computes ECP matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param atom_potential The atom core potential.
    /// @return The ECP matrix.
    auto compute_geom_200(const CMolecularBasis &basis, const CMolecule &molecule, const std::vector<int>& atoms, const int iatom) const -> CMatrices;
    
    /// @brief Computes ECP matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param atom_potential The atom core potential.
    /// @return The ECP matrix.
    auto compute_geom_101(const CMolecularBasis &basis, const CMolecule &molecule, const CAtomCorePotential& atom_potential, const int iatom, const int jatom) const -> CMatrices;
    
    /// @brief Computes ECP matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param atom_potential The atom core potential.
    /// @return The ECP matrix.
    auto compute_geom_101(const CMolecularBasis &basis, const CMolecule &molecule, const std::vector<int>& atoms, const int iatom, const int jatom) const -> CMatrices;
    
    /// @brief Computes ECP matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param atom_potential The atom core potential.
    /// @return The ECP matrix.
    auto compute_geom_110(const CMolecularBasis &basis, const CMolecule &molecule, const CAtomCorePotential& atom_potential, const int iatom) const -> CMatrices;
    
    /// @brief Computes ECP matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param atom_potential The atom core potential.
    /// @return The ECP matrix.
    auto compute_geom_110(const CMolecularBasis &basis, const CMolecule &molecule, const int iatom, const int jatom) const -> CMatrices;
    
    /// @brief Computes ECP matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param atom_potential The atom core potential.
    /// @return The ECP matrix.
    auto compute_geom_020(const CMolecularBasis &basis, const CMolecule &molecule, const CAtomCorePotential& atom_potential) const -> CMatrices;
    
    /// @brief Computes ECP matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param atom_potential The atom core potential.
    /// @return The ECP matrix.
    auto compute_geom_020(const CMolecularBasis &basis, const CMolecule &molecule, const int iatom) const -> CMatrices;
};

#endif /* CorePotentialHessianDriver_hpp */
