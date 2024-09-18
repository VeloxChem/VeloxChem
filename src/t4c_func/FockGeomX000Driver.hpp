#ifndef FockGeomX000Driver_hpp
#define FockGeomX000Driver_hpp

#include <string>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Matrix.hpp"
#include "Matrices.hpp"

/// Class CFockGeomX000Driver provides methods for computing Fock matrices
/// using derivatives of four center electron repulsion integrals.
template<int N>
class CFockGeomX000Driver
{
   public:
    /// Creates a Fock matrices  driver.
    CFockGeomX000Driver() = default;

    /// @brief The default copy constructor.
    /// @param other The Fock matrices driver to be copied.
    CFockGeomX000Driver(const CFockDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The Fock matrices driver  to be moved.
    CFockGeomX000Driver(CFockDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CFockGeomX000Driver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The Fock matrices driver to be copy assigned.
    /// @return The assigned Fock matrices driver.
    auto operator=(const CFockGeomX000Driver &other) -> CFockDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The Fock matrices driver to be move assigned.
    /// @return The assigned Fock matrices driver .
    auto operator=(CFockGeomX000Driver &&other) noexcept -> CFockDriver & = delete;

    /// @brief The equality operator.
    /// @param other The Fock matrices driver  to be compared.
    /// @return True if Fock matrices drivers  are equal, False otherwise.
    auto operator==(const CFockGeomX000Driver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The Fock matrices driver to be compared.
    /// @return True if Fock matrices drivers  are not equal, False otherwise.
    auto operator!=(const CFockGeomX000Driver &other) const -> bool = delete;

    /// @brief Computes Fock matrix derivatives for given density, basis and molecule (N^4 scaling).
    /// @param iatom The index of atom to compute derivatives of Fock matrix.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param density The density matrix to construct Fock matrix.
    /// @param label The label of Fock matrix type.
    /// @param exchange_factor The exchange-correlation factors.
    /// @param omega The range separation factor.
    /// @return The Fock matrices.
    auto compute(const int             iatom,
                 const CMolecularBasis &basis,
                 const CMolecule       &molecule,
                 const CMatrix         &density,
                 const std::string     &label,
                 const double          exchange_factor,
                 const double          omega) const -> CMatrices;
};

template<int N>
auto
CFockGeomX000Driver::compute(const int             iatom,
                             const CMolecularBasis &basis,
                             const CMolecule       &molecule,
                             const CMatrix         &density,
                             const std::string     &label,
                             const double          exchange_factor,
                             const double          omega) const -> CMatrices
{
    return CMatrices(); 
}

#endif /* FockGeomX000Driver_hpp */
