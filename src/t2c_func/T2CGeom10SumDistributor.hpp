#ifndef T2CGeom10SumDistributor_hpp
#define T2CGeom10SumDistributor_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "SimdArray.hpp"

/// @brief Class CT3CGeom100SumDistributor provides methods for distributing vector of into flat buffer.
class CT2CGeom10SumDistributor
{
   public:
    /// @brief Creates an integral shells distributor.
    /// @param values  The pointer to gradient values.
    /// @param gamma The pointer to B_q vector data.
    CT2CGeom10SumDistributor(double* values, const double* gamma);

    /// @brief The default copy constructor.
    /// @param other The distributor to be copied.
    CT2CGeom10SumDistributor(const CT2CGeom10SumDistributor& other) = delete;

    /// @brief The default move constructor.
    /// @param other The distributor to be moved.
    CT2CGeom10SumDistributor(CT2CGeom10SumDistributor&& other) noexcept = delete;

    /// @brief The default destructor.
    ~CT2CGeom10SumDistributor() = default;

    /// @brief The default copy assignment operator.
    /// @param other The distributor to be copy assigned.
    /// @return The assigned distributor.
    auto operator=(const CT2CGeom10SumDistributor& other) -> CT2CGeom10SumDistributor& = delete;

    /// @brief The default move assignment operator.
    /// @param other The distributor to be move assigned.
    /// @return The assigned distributor.
    auto operator=(CT2CGeom10SumDistributor&& other) noexcept -> CT2CGeom10SumDistributor& = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are equal, False otherwise.
    auto operator==(const CT2CGeom10SumDistributor& other) const -> bool = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are not equal, False otherwise.
    auto operator!=(const CT2CGeom10SumDistributor& other) const -> bool = delete;
    
    /// @brief Distributes buffer of integrals into storage.
    /// @param buffer The integrals buffer.
    /// @param bra_indices The compressed contracted basis functions indexes on bra side.
    /// @param ket_indices The compressed contracted basis functions indexes on ket side.
    /// @param bra_angmom The angular momentum of integrals buffer on bra side.
    /// @param ket_angmom The angular momentum of integrals buffer on ket side.
    /// @param bra_igto The index of the basis function on bra side.
    /// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
    /// @param diagonal True if basis functions blocks on bra and ket are the same, False otherwise.
    auto distribute(const CSimdArray<double>&        buffer,
                    const std::vector<size_t>&       bra_indices,
                    const std::vector<size_t>&       ket_indices,
                    const int                        bra_angmom,
                    const int                        ket_angmom,
                    const size_t                     bra_igto,
                    const std::pair<size_t, size_t>& ket_range,
                    const bool                       diagonal) -> void;
    
   private:
    /// @brief The pointer to gradient values.
    double* _grad_values;
    
    /// @brief The pointer to B_q vector values.
    const double* _ptr_gamma;
};

#endif /* T2CGeom10SumDistributor_hpp */
