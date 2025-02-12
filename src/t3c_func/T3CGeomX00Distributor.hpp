#ifndef T3CGeomX00Distributor_hpp
#define T3CGeomX00Distributor_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "SimdArray.hpp"
#include "T3FlatBuffer.hpp"

/// @brief Class CT3CGeomX00Distributor provides methods for distributing vector of into flat buffer.
class CT3CGeomX00Distributor
{
   public:
    /// @brief Creates an integral shells distributor.
    /// @param values  The pointer to flat tensor.
    CT3CGeomX00Distributor(CT3FlatBuffer<double>* values);

    /// @brief The default copy constructor.
    /// @param other The distributor to be copied.
    CT3CGeomX00Distributor(const CT3CGeomX00Distributor& other) = delete;

    /// @brief The default move constructor.
    /// @param other The distributor to be moved.
    CT3CGeomX00Distributor(CT3CGeomX00Distributor&& other) noexcept = delete;

    /// @brief The default destructor.
    ~CT3CGeomX00Distributor() = default;

    /// @brief The default copy assignment operator.
    /// @param other The distributor to be copy assigned.
    /// @return The assigned distributor.
    auto operator=(const CT3CGeomX00Distributor& other) -> CT3CGeomX00Distributor& = delete;

    /// @brief The default move assignment operator.
    /// @param other The distributor to be move assigned.
    /// @return The assigned distributor.
    auto operator=(CT3CGeomX00Distributor&& other) noexcept -> CT3CGeomX00Distributor& = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are equal, False otherwise.
    auto operator==(const CT3CGeomX00Distributor& other) const -> bool = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are not equal, False otherwise.
    auto operator!=(const CT3CGeomX00Distributor& other) const -> bool = delete;

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
    
    /// @brief Sets local index to flat buffer.
    auto set_index(const size_t local_index) -> void;
    
   private:
    /// @brief The pointer to flat buffer.
    CT3FlatBuffer<double>* _t3_values;
    
    /// @brief The local index.
    size_t _local_index;
};

#endif /* T3CGeomX00Distributor_hpp */
