#ifndef FockDriver_hpp
#define FockDriver_hpp

#include <string>

#include "Matrices.hpp"
#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "T4CScreener.hpp"

/// Class CFockDriver provides methods for computing Fock matrices
/// using four center electron repulsion integrals.
class CFockDriver
{
   public:
    /// Creates a Fock matrices  driver.
    CFockDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The Fock matrices driver to be copied.
    CFockDriver(const CFockDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The Fock matrices driver  to be moved.
    CFockDriver(CFockDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CFockDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The Fock matrices driver to be copy assigned.
    /// @return The assigned Fock matrices driver.
    auto operator=(const CFockDriver &other) -> CFockDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The Fock matrices driver to be move assigned.
    /// @return The assigned Fock matrices driver .
    auto operator=(CFockDriver &&other) noexcept -> CFockDriver & = delete;

    /// @brief The equality operator.
    /// @param other The Fock matrices driver  to be compared.
    /// @return True if Fock matrices drivers  are equal, False otherwise.
    auto operator==(const CFockDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The Fock matrices driver to be compared.
    /// @return True if Fock matrices drivers  are not equal, False otherwise.
    auto operator!=(const CFockDriver &other) const -> bool = delete;

    /// @brief Computes Fock matrix for given density, basis and molecule (N^4 scaling).
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param density The density matrix to construct Fock matrix.
    /// @param label The label of Fock matrix type.
    /// @param exchange_factor The exchange-correlation factors.
    /// @param omega The range separation factor.
    /// @return The Fock matrix.
    auto compute(const CMolecularBasis &basis,
                 const CMolecule       &molecule,
                 const CMatrix         &density,
                 const std::string     &label,
                 const double           exchange_factor,
                 const double           omega) const -> CMatrix;

    /// @brief Computes Fock matrix for given density..
    /// @param screener The screener with basis function pairs data.
    /// @param density The density matrix to construct Fock matrix.
    /// @param label The label of Fock matrix type.
    /// @param exchange_factor The exchange-correlation factors.
    /// @param omega The range separation factor.
    /// @return The Fock matrix.
    auto compute(const CT4CScreener &screener,
                 const CMatrix      &density,
                 const std::string  &label,
                 const double        exchange_factor,
                 const double        omega,
                 const int           ithreshold) const -> CMatrix;

    /// @brief Computes Fock matrices for given densities.
    /// @param screener The screener with basis function pairs data.
    /// @param densities The density matrices to construct Fock matrix.
    /// @param labels The vector of  Fock matrix type labels.
    /// @param exchange_factor The exchange-correlation factors.
    /// @param omega The range separation factor.
    /// @return The Fock matrix.
    auto compute(const CT4CScreener             &screener,
                 const CMatrices                &densities,
                 const std::vector<std::string> &labels,
                 const double                    exchange_factor,
                 const double                    omega,
                 const int                       ithreshold) const -> CMatrices;

    /// @brief Computes Fock matrix for given density..
    /// @param screener The screener with basis function pairs data.
    /// @param rank The rank of specific node.
    /// @param nodes The number of nodes.
    /// @param density The density matrix to construct Fock matrix.
    /// @param label The label of Fock matrix type.
    /// @param exchange_factor The exchange-correlation factors.
    /// @param omega The range separation factor.
    /// @return The Fock matrix.
    auto compute(const CT4CScreener &screener,
                 const int           rank,
                 const int           nodes,
                 const CMatrix      &density,
                 const std::string  &label,
                 const double        exchange_factor,
                 const double        omega,
                 const int           ithreshold) const -> CMatrix;

    /// @brief Computes Fock matrices for given densities.
    /// @param screener The screener with basis function pairs data.
    /// @param rank The rank of specific node.
    /// @param nodes The number of nodes.
    /// @param densities The density matrices to construct Fock matrix.
    /// @param labels The vector of  Fock matrix type labels.
    /// @param exchange_factor The exchange-correlation factors.
    /// @param omega The range separation factor.
    /// @return The Fock matrix.
    auto compute(const CT4CScreener             &screener,
                 const int                       rank,
                 const int                       nodes,
                 const CMatrices                &densities,
                 const std::vector<std::string> &labels,
                 const double                    exchange_factor,
                 const double                    omega,
                 const int                       ithreshold) const -> CMatrices;

    auto set_block_size_factor(const int factor) -> void;

   private:
    int _block_size_factor = 16;

    auto _determine_block_size_factor(const int nao) const -> int;

    auto _get_nao(const CMatrix& mat) const -> int;

    auto _get_nao(const CMatrices& mats) const -> int;
};

#endif /* FockDriver_hpp */
