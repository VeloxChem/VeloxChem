#ifndef T3CGeom0X0Distributor_hpp
#define T3CGeom0X0Distributor_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "SimdArray.hpp"
#include "T3RectFlatBuffer.hpp"

/// @brief Class CT3CGeom0X0Distributor provides methods for distributing vector of into flat buffer.
class CT3CGeom0X0Distributor
{
   public:
    /// @brief Creates an integral shells distributor.
    /// @param values  The pointer to flat tensor.
    CT3CGeom0X0Distributor(CT3RectFlatBuffer<double>* values);

    /// @brief The default copy constructor.
    /// @param other The distributor to be copied.
    CT3CGeom0X0Distributor(const CT3CGeom0X0Distributor& other) = delete;

    /// @brief The default move constructor.
    /// @param other The distributor to be moved.
    CT3CGeom0X0Distributor(CT3CGeom0X0Distributor&& other) noexcept = delete;

    /// @brief The default destructor.
    ~CT3CGeom0X0Distributor() = default;

    /// @brief The default copy assignment operator.
    /// @param other The distributor to be copy assigned.
    /// @return The assigned distributor.
    auto operator=(const CT3CGeom0X0Distributor& other) -> CT3CGeom0X0Distributor& = delete;

    /// @brief The default move assignment operator.
    /// @param other The distributor to be move assigned.
    /// @return The assigned distributor.
    auto operator=(CT3CGeom0X0Distributor&& other) noexcept -> CT3CGeom0X0Distributor& = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are equal, False otherwise.
    auto operator==(const CT3CGeom0X0Distributor& other) const -> bool = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are not equal, False otherwise.
    auto operator!=(const CT3CGeom0X0Distributor& other) const -> bool = delete;

    /// @brief Distributes buffer of integrals into local Fock matrix.
    /// @param buffer The integrals buffer.
    /// @param a_indices The compressed basis function indexes on center A.
    /// @param c_indices The compressed basis function indexes on center C.
    /// @param d_indices The compressed basis function indexes on center D.
    /// @param a_angmom The angular momentum of integrals buffer on center A.
    /// @param c_angmom The angular momentum of integrals buffer on center C.
    /// @param d_angmom Tthe angular momentum of integrals buffer on center D.
    /// @param ibra_gto The index of basis function on bra side.
    /// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
    auto distribute(const CSimdArray<double>&        buffer,
                    const size_t                     offset,
                    const std::vector<size_t>&       a_indices,
                    const std::vector<size_t>&       c_indices,
                    const std::vector<size_t>&       d_indices,
                    const int                        a_angmom,
                    const int                        c_angmom,
                    const int                        d_angmom,
                    const size_t                     ibra_gto,
                    const std::pair<size_t, size_t>& ket_range) -> void;
    
   private:
    /// @brief The pointer to flat buffer.
    CT3RectFlatBuffer<double>* _t3_values;
};

#endif /* T3CGeom0X0Distributor_hpp */
