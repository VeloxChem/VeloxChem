#ifndef OverlapDriver_hpp
#define OverlapDriver_hpp

#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

/// @brief Class COverlapDriver provides methods for computing two-center overlap integrals.
class COverlapDriver
{
   public:
    /// @brief Creates an overlap integrals driver.
    COverlapDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The overlap integrals driver to be copied.
    COverlapDriver(const COverlapDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The overlap integrals driver  to be moved.
    COverlapDriver(COverlapDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~COverlapDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The overlap integrals driver to be copy assigned.
    /// @return The assigned overlap integrals driver.
    auto operator=(const COverlapDriver &other) -> COverlapDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The overlap integrals driver to be move assigned.
    /// @return The assigned overlap integrals driver .
    auto operator=(COverlapDriver &&other) noexcept -> COverlapDriver & = delete;

    /// @brief The equality operator.
    /// @param other The overlap integrals driver  to be compared.
    /// @return True if overlap integrals drivers  are equal, False otherwise.
    auto operator==(const COverlapDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The overlap integrals driver to be compared.
    /// @return True if overlap integrals drivers  are not equal, False otherwise.
    auto operator!=(const COverlapDriver &other) const -> bool = delete;

    /// @brief Computes overlap matrix for given molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @return The overlap matrix.
    auto compute(const CMolecularBasis &basis, const CMolecule &molecule) const -> CMatrix;

    /// @brief Computes overlap matrix for given molecule and pair of molecular bases.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis the molecular basis on ket side.
    /// @param molecule The molecule.
    /// @return The overlap matrix.
    auto compute(const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis, const CMolecule &molecule) const -> CMatrix;

    /// @brief Computes overlap matrix for given pair of molecules and pair of molecular bases.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis the molecular basis on ket side.
    /// @param bra_molecule The molecule on bra side.
    /// @param ket_molecule The molecule on ket side.
    /// @return The overlap matrix.
    auto compute(const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis, const CMolecule& bra_molecule, const CMolecule& ket_molecule) const -> CMatrix;
};

#endif /* OverlapDriver_hpp */
