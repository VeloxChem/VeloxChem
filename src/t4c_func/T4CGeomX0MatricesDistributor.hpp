#ifndef T4CGeomX0MatricesDistributor_hpp
#define T4CGeomX0MatricesDistributor_hpp

#include <array>
#include <string>
#include <vector>

#include "GtoPairBlock.hpp"
#include "Matrices.hpp"
#include "Matrix.hpp"
#include "SimdArray.hpp"

/// @brief Class CT4CGeomX0MatricesDistributor provides methods for distributing single Fock matrices associated
/// with density matrix.
class CT4CGeomX0MatricesDistributor
{
   public:
    /// @brief The default destructor.
    CT4CGeomX0MatricesDistributor() = default;

    /// Creates a Fock matrices distributor.
    /// @param focks  The Fock matrices.
    /// @param density  The density matrix.
    /// @param label  The standard  Fock matrix type labels.
    /// @param exchange_factor  The scaling factor of exchange contribution.
    /// @param omega The range separation factor.
    CT4CGeomX0MatricesDistributor(CMatrices*         focks,
                                  const CMatrix*     density,
                                  const std::string& label,
                                  const double       exchange_factor,
                                  const double       omega);

    CT4CGeomX0MatricesDistributor(CMatrices*         focks,
                                  const CMatrix*     density,
                                  const CMatrix*     density_2,
                                  const std::string& label,
                                  const double       exchange_factor,
                                  const double       omega);

    /// @brief The default copy constructor.
    /// @param other The distributor to be copied.
    CT4CGeomX0MatricesDistributor(const CT4CGeomX0MatricesDistributor& other) = delete;

    /// @brief The default move constructor.
    /// @param other The distributor to be moved.
    CT4CGeomX0MatricesDistributor(CT4CGeomX0MatricesDistributor&& other) noexcept = delete;

    /// @brief The default destructor.
    ~CT4CGeomX0MatricesDistributor() = default;

    /// @brief The default copy assignment operator.
    /// @param other The distributor to be copy assigned.
    /// @return The assigned distributor.
    auto operator=(const CT4CGeomX0MatricesDistributor& other) -> CT4CGeomX0MatricesDistributor& = delete;

    /// @brief The default move assignment operator.
    /// @param other The distributor to be move assigned.
    /// @return The assigned distributor.
    auto operator=(CT4CGeomX0MatricesDistributor&& other) noexcept -> CT4CGeomX0MatricesDistributor& = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are equal, False otherwise.
    auto operator==(const CT4CGeomX0MatricesDistributor& other) const -> bool = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are not equal, False otherwise.
    auto operator!=(const CT4CGeomX0MatricesDistributor& other) const -> bool = delete;

    /// @brief The gets range separation factors.
    /// @return The range separation
    auto get_omega() const -> double;

    /// @brief Checks if range separation factor is needed.
    /// @return The range separation
    auto need_omega() const -> bool;

    auto get_values() const -> std::vector<double>;

    /// Sets local matrices and their local/global indices.
    /// @param bra_gto_pair_block The basis function pairs block on bra side.
    /// @param ket_gto_pair_block The basis function pairs block on ket side.
    auto set_indices(const CGtoPairBlock& bra_gto_pair_block, const CGtoPairBlock& ket_gto_pair_block) -> void;

    /// @brief Distributes buffer of integrals into local Fock matrix.
    /// @param buffer The integrals buffer.
    /// @param a_indices The compressed basis function indexes on center A.
    /// @param b_indices The compressed basis function indexes on center B.
    /// @param c_indices The compressed basis function indexes on center C.
    /// @param d_indices The compressed basis function indexes on center D.
    /// @param a_angmom The angular momentum of integrals buffer on center A.
    /// @param b_angmom The angular momentum of integrals buffer on center B.
    /// @param c_angmom The angular momentum of integrals buffer on center C.
    /// @param d_angmom Tthe angular momentum of integrals buffer on center D.
    /// @param ibra_gto The index of basis function on bra side.
    /// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
    auto distribute(const CSimdArray<double>&        buffer,
                    const size_t                     offset,
                    const std::vector<size_t>&       a_indices,
                    const std::vector<size_t>&       b_indices,
                    const std::vector<size_t>&       c_indices,
                    const std::vector<size_t>&       d_indices,
                    const int                        a_angmom,
                    const int                        b_angmom,
                    const int                        c_angmom,
                    const int                        d_angmom,
                    const size_t                     ibra_gto,
                    const std::pair<size_t, size_t>& ket_range) -> void;

    /// Accumulates local Fock matrices contributions to targeted Fock matrix.
    /// @param bra_gto_pair_block The basis function pairs block on bra side.
    /// @param ket_gto_pair_block The basis function pairs block on ket side.
    auto accumulate(const CGtoPairBlock& bra_gto_pair_block, const CGtoPairBlock& ket_gto_pair_block) -> void;

    auto accumulate_values(double* values, const CMatrix* density, const CGtoPairBlock& bra_gto_pair_block, const CGtoPairBlock& ket_gto_pair_block) -> void;

   private:
    /// @brief The Fock matrices associated with distributor.
    CMatrices* _focks;

    std::vector<double> _values;

    /// @brief The density matrix associated with distributor.
    const CMatrix* _density;

    /// @brief The 2nd density matrix associated with distributor.
    const CMatrix* _density_2;

    /// @brief The standard Fock matrix type label.
    std::string _label;

    /// @brief The scalling factor for scaling exchange contribution.
    double _exchange_factor;

    /// @brief The range separation factor.
    double _omega;

    /// @brief The local storage matrices.
    CMatrices _matrices;

    /// @brief The local indices for center A.
    std::vector<size_t> _a_loc_indices;

    /// @brief The local indices for center B.
    std::vector<size_t> _b_loc_indices;

    /// @brief The local indices for center C.
    std::vector<size_t> _c_loc_indices;

    /// @brief The local indices for center D.
    std::vector<size_t> _d_loc_indices;

    /// @brief The global indices for center A.
    std::vector<size_t> _a_glob_indices;

    /// @brief The global indices for center B.
    std::vector<size_t> _b_glob_indices;

    /// @brief The global indices for center C.
    std::vector<size_t> _c_glob_indices;

    /// @brief The global indices for center D.
    std::vector<size_t> _d_glob_indices;
};

#endif /* T4CGeomX0MatricesDistributor_hpp */
