#ifndef T2CDistributor_hpp
#define T2CDistributor_hpp

#include <algorithm>
#include <array>
#include <string>
#include <vector>

#include "Matrices.hpp"
#include "Matrix.hpp"
#include "Point.hpp"
#include "SimdArray.hpp"
#include "T2CUtils.hpp"
#include "TensorComponents.hpp"

/// @brief Class CT2CDistributor provides methods for distributing two-center integrals data.
template <class T>
class CT2CDistributor
{
   public:
    /// @brief The default constructor.
    /// @param storage  The storage associated with distributor.
    CT2CDistributor(T* storage) : _storage{storage}, _coordinates{}, _data{} {};

    /// @brief The constructor with storage, external coordinates and external data.
    /// @param storage  The storage associated with distributor.
    /// @param coordinates  The Cartesian coordinates of external points.
    /// @param data  The data assocated with external points.
    CT2CDistributor(T* storage, const std::vector<TPoint<double>>& coordinates, const std::vector<double>& data)
        : _storage{storage}

        , _coordinates{coordinates}

        , _data{data} {};
    
    /// @brief The default copy constructor.
    /// @param other The distributor to be copied.
    CT2CDistributor(const CT2CDistributor &other) = delete;

    /// @brief The default move constructor.
    /// @param other The distributor to be moved.
    CT2CDistributor(CT2CDistributor &&other) noexcept = delete;
    
    /// @brief The default destructor.
    ~CT2CDistributor() = default;

    /// @brief The default copy assignment operator.
    /// @param other The distributor to be copy assigned.
    /// @return The assigned distributor.
    auto operator=(const CT2CDistributor &other) -> CT2CDistributor & = delete;

    /// @brief The default move assignment operator.
    /// @param other The distributor to be move assigned.
    /// @return The assigned distributor.
    auto operator=(CT2CDistributor &&other) noexcept -> CT2CDistributor & = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are equal, False otherwise.
    auto operator==(const CT2CDistributor &other) const -> bool = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are not equal, False otherwise.
    auto operator!=(const CT2CDistributor &other) const -> bool = delete;
    
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
                    const bool                       diagonal) -> void {};

    /// @brief Gets Cartesian  coordinates of external points.
    /// @return The vector of Cartesian coordinates of external points.
    auto
    coordinates() const -> std::vector<TPoint<double>>
    {
        return _coordinates;
    };

    /// @brief Gets flat data associated with external points.
    /// @return The flat vector with data associated with external points.
    auto
    data() const -> std::vector<double>
    {
        return _data;
    };

   private:
    ///@brief  The storage container associated with distributor.
    T* _storage;

    ///@brief The Cartesian coordinates of external points.
    std::vector<TPoint<double>> _coordinates;

    /// @brief The flat array of external parameters associated with external points.
    std::vector<double> _data;
};

template <>
inline auto
CT2CDistributor<CMatrix>::distribute(const CSimdArray<double>&        buffer,
                                     const std::vector<size_t>&       bra_indices,
                                     const std::vector<size_t>&       ket_indices,
                                     const int                        bra_angmom,
                                     const int                        ket_angmom,
                                     const size_t                     bra_igto,
                                     const std::pair<size_t, size_t>& ket_range,
                                     const bool                       diagonal) -> void
{
    if (bra_angmom == ket_angmom)
    {
        if (diagonal)
        {
            t2cfunc::distribute(_storage->sub_matrix({bra_angmom, bra_angmom}), buffer, size_t{0}, bra_indices, bra_angmom, bra_igto, ket_range);
        }
        else
        {
            t2cfunc::distribute(_storage->sub_matrix({bra_angmom, ket_angmom}),
                                buffer,
                                size_t{0},
                                bra_indices,
                                ket_indices,
                                bra_angmom,
                                ket_angmom,
                                bra_igto,
                                ket_range,
                                _storage->get_type());
        }
    }
    else
    {
        t2cfunc::distribute(_storage->sub_matrix({bra_angmom, ket_angmom}),
                            buffer,
                            size_t{0},
                            bra_indices,
                            ket_indices,
                            bra_angmom,
                            ket_angmom,
                            bra_igto,
                            ket_range,
                            _storage->is_angular_order({bra_angmom, ket_angmom}));
    }
}

template <>
inline auto
CT2CDistributor<CMatrices>::distribute(const CSimdArray<double>&        buffer,
                                       const std::vector<size_t>&       bra_indices,
                                       const std::vector<size_t>&       ket_indices,
                                       const int                        bra_angmom,
                                       const int                        ket_angmom,
                                       const size_t                     bra_igto,
                                       const std::pair<size_t, size_t>& ket_range,
                                       const bool                       diagonal) -> void
{
    const auto tcomps = tensor::number_of_spherical_components(std::array<int, 2>({bra_angmom, ket_angmom}));

    size_t buff_off = 0;

    std::ranges::for_each(_storage->keys(), [&](const auto& key) {
        auto matrix = _storage->matrix(key);
        if (bra_angmom == ket_angmom)
        {
            if (diagonal)
            {
                t2cfunc::distribute(matrix->sub_matrix({bra_angmom, bra_angmom}), buffer, buff_off, bra_indices, bra_angmom, bra_igto, ket_range);
            }
            else
            {
                t2cfunc::distribute(matrix->sub_matrix({bra_angmom, ket_angmom}),
                                    buffer,
                                    buff_off,
                                    bra_indices,
                                    ket_indices,
                                    bra_angmom,
                                    ket_angmom,
                                    bra_igto,
                                    ket_range,
                                    matrix->get_type());
            }
        }
        else
        {
            t2cfunc::distribute(matrix->sub_matrix({bra_angmom, ket_angmom}),
                                buffer,
                                buff_off,
                                bra_indices,
                                ket_indices,
                                bra_angmom,
                                ket_angmom,
                                bra_igto,
                                ket_range,
                                matrix->is_angular_order({bra_angmom, ket_angmom}));
        }
        buff_off += tcomps;
    });
}

#endif /* T2CMatrixDistributor_hpp */
