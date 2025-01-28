#ifndef RIFockDriver_hpp
#define RIFockDriver_hpp

#include <vector>

#include "Matrix.hpp"
#include "T3FlatBuffer.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"

/// Class CRIFockDriver provides methods for computing Fock matrices
/// using three center electron repulsion integrals.
class CRIFockDriver
{
   public:
    /// Creates a Fock matrices  driver.
    CRIFockDriver();
    
    /// Creates a Fock matrices  driver.
    /// @param j_metric The metric matrix for J fitting.
    CRIFockDriver(const CMatrix& j_metric);
    
    /// Creates a Fock matrices  driver.
    /// @param j_metric The metric matrix for J fitting.
    /// @param k_metric The metric matrix for K fitting.
    CRIFockDriver(const CMatrix& j_metric, const CMatrix& k_metric);

    /// @brief The default copy constructor.
    /// @param other The Fock matrices driver to be copied.
    CRIFockDriver(const CRIFockDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The Fock matrices driver  to be moved.
    CRIFockDriver(CRIFockDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CRIFockDriver();

    /// @brief The default copy assignment operator.
    /// @param other The Fock matrices driver to be copy assigned.
    /// @return The assigned Fock matrices driver.
    auto operator=(const CRIFockDriver &other) -> CRIFockDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The Fock matrices driver to be move assigned.
    /// @return The assigned Fock matrices driver .
    auto operator=(CRIFockDriver &&other) noexcept -> CRIFockDriver & = delete;

    /// @brief The equality operator.
    /// @param other The Fock matrices driver  to be compared.
    /// @return True if Fock matrices drivers  are equal, False otherwise.
    auto operator==(const CRIFockDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The Fock matrices driver to be compared.
    /// @return True if Fock matrices drivers  are not equal, False otherwise.
    auto operator!=(const CRIFockDriver &other) const -> bool = delete;
    
    /// @brief Computes three center electron repulsion integral buffers.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxilary molecular  basis.
    auto prepare_buffers(const CMolecule       &molecule,
                         const CMolecularBasis &basis,
                         const CMolecularBasis &aux_basis) -> void;
    
    /// @brief Computes Fock matrix for given density.
    /// @param density The density matrix to construct Fock matrix.
    /// @param label The label of Fock matrix type.
    /// @param exchange_factor The exchange-correlation factors.
    /// @return The Fock matrix.
    auto compute(const CMatrix     &density,
                 const std::string &label,
                 const double      exchange_factor) const -> CMatrix;
    
    /// @brief Computes Gamma vector  for given density.
    /// @param density The density matrix to construct Fock matrix.
    /// @return The Gamma vector.
    auto comp_gamma_vector(const CMatrix &density) const -> std::vector<double>;
    
    private:
    /// @brief Pointer to metric matrix for J fitting.
    CMatrix* _j_metric;
    
    /// @brief Pointer to metric matrix for K fitting.
    CMatrix* _k_metric;
    
    /// @brief Three center electron repulsion integrals buffer.
    CT3FlatBuffer<double> _eri_buffer;
    
    /// @brief Three transformed center electron repulsion integrals buffer.
    CT3FlatBuffer<double> _eri_trafo_buffer;
};

#endif /* RIFockDriver_hpp */
