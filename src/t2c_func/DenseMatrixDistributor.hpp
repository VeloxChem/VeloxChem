#ifndef DenseMatrixDistributor_hpp
#define DenseMatrixDistributor_hpp

#include <vector>

#include "Point.hpp"
#include "SimdArray.hpp"
#include "DenseMatrix.hpp"

/// @brief Class CDenseMatrixDistributor provides methods for distributing two-center integrals data.
class CDenseMatrixDistributor
{
   public:
    /// @brief The default constructor.
    CDenseMatrixDistributor();

    /// @brief The constructor with storage, external coordinates and external data.
    /// @param coordinates  The Cartesian coordinates of external points.
    /// @param data  The data assocated with external points.
    CDenseMatrixDistributor(CDenseMatrix*                      g_matrix,
                            const std::vector<TPoint<double>>& coordinates,
                            const std::vector<double>&         data,
                            const CDenseMatrix*                f_matrix,
                            const double                       weight);
    
    /// @brief The default copy constructor.
    /// @param other The distributor to be copied.
    CDenseMatrixDistributor(const CDenseMatrixDistributor &other) = delete;

    /// @brief The default move constructor.
    /// @param other The distributor to be moved.
    CDenseMatrixDistributor(CDenseMatrixDistributor &&other) noexcept = delete;
    
    /// @brief The default destructor.
    ~CDenseMatrixDistributor() = default;

    /// @brief The default copy assignment operator.
    /// @param other The distributor to be copy assigned.
    /// @return The assigned distributor.
    auto operator=(const CDenseMatrixDistributor &other) -> CDenseMatrixDistributor & = delete;

    /// @brief The default move assignment operator.
    /// @param other The distributor to be move assigned.
    /// @return The assigned distributor.
    auto operator=(CDenseMatrixDistributor &&other) noexcept -> CDenseMatrixDistributor & = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are equal, False otherwise.
    auto operator==(const CDenseMatrixDistributor &other) const -> bool = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are not equal, False otherwise.
    auto operator!=(const CDenseMatrixDistributor &other) const -> bool = delete;
    
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
    ///@brief The pointer to G matrix.
    CDenseMatrix* _g_matrix;
    
    ///@brief The Cartesian coordinates of external points.
    std::vector<TPoint<double>> _coordinates;

    /// @brief The flat array of external parameters associated with external points.
    std::vector<double> _data;
    
    ///@brief The pointer to F matrix.
    const CDenseMatrix* _f_matrix;
    
    ///@brief The weight of grid point.
    const double _weight;
};

#endif /* DenseMatrixDistributor_hpp */
