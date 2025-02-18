#ifndef T3RectFlatBuffer_hpp
#define T3RectFlatBuffer_hpp

#include <cstddef>
#include <vector>
#include <ranges>

/// @brief Class CT3RectFlatBuffer stores general semi-symmetric rank-3 tensor as 2D flattened data structure.
template <typename T>
class CT3RectFlatBuffer
{
   public:
    
    /// @brief The default constructor.
    /// @param other The SIMD array to be copied.
    CT3RectFlatBuffer() = default;
    
    /// @brief The default constructor.
    /// @param indices The vector of indices along x axis of tensor.
    /// @param row_indices The vector of indices along y axis of tensor.
    /// @param col_indices The vector of indices along z axis of tensor.
    CT3RectFlatBuffer(const std::vector<size_t>& indices,
                      const std::vector<size_t>& row_indices,
                      const std::vector<size_t>& col_indices)
    {
        _indices = indices;
        
        _row_indices = row_indices;
        
        _col_indices = col_indices;
        
        _data.reserve(_indices.size());
        
        if (const auto nelems = _row_indices.size() * _col_indices.size();  nelems > 0)
        {
            std::ranges::for_each(_indices, [&](const auto& index) {
                _data.push_back(std::vector<T>(nelems, T{0.0}));
            });
        }
    }
    
    /// @brief The default constructor.
    /// @param indices The vector of indices along x axis of tensor.
    /// @param row_indices The vector of indices along y axis of tensor.
    /// @param col_indices The vector of indices along z axis of tensor.
    /// @param nbatches The number of batches.
    CT3RectFlatBuffer(const std::vector<size_t>& indices,
                      const std::vector<size_t>& row_indices,
                      const std::vector<size_t>& col_indices,
                      const size_t               nbatches)
    {
        _indices = indices;
        
        _row_indices = row_indices;
        
        _col_indices = col_indices;
        
        _data.reserve(_indices.size() * nbatches);
        
        if (const auto nelems = _row_indices.size() * _col_indices.size();  nelems > 0)
        {
            std::ranges::for_each(std::views::iota(size_t{0}, nbatches), [&] (const auto index) {
                std::ranges::for_each(_indices, [&](const auto& index) {
                    _data.push_back(std::vector<T>(nelems, T{0.0}));
                });
            });
        }
    }

    /// @brief The default copy constructor.
    /// @param other The SIMD array to be copied.
    CT3RectFlatBuffer(const CT3RectFlatBuffer &other) = default;

    /// @brief The default move constructor.
    /// @param other The SIMD array to be moved.
    CT3RectFlatBuffer(CT3RectFlatBuffer &&other) noexcept = default;

    /// @brief The custom destructor.
    ~CT3RectFlatBuffer() = default;
    
    /// @brief The default copy assignment operator.
    /// @param other The SIMD array to be copy assigned.
    /// @return The assigned SIMD array.
    auto operator=(const CT3RectFlatBuffer &other) -> CT3RectFlatBuffer & = default;

    /// @brief The default move assignment operator.
    /// @param other The SIMD array to be move assigned.
    /// @return The assigned SIMD array.
    auto operator=(CT3RectFlatBuffer &&other) noexcept -> CT3RectFlatBuffer & = default;

    /// @brief The equality operator.
    /// @param other The SIMD array to be compared.
    /// @return True if SIMD arrays are equal, False otherwise.
    auto operator==(const CT3RectFlatBuffer &other) const -> bool = default;

    /// @brief The non-equality operator.
    /// @param other The SIMD array to be compared.
    /// @return True if SIMD arrays are not equal, False otherwise.
    auto operator!=(const CT3RectFlatBuffer &other) const -> bool = default;
    
    /// @brief The gets vector of indices along x axis.
    /// @return The vector of indices.
    inline auto
    indices() const -> std::vector<size_t>
    {
        return _indices;
    }
    
    /// @brief The gets vector of indices along y axis.
    /// @return The vector of indices.
    inline auto
    row_indices() const -> std::vector<size_t>
    {
        return _row_indices;
    }
    
    /// @brief The gets vector of indices along z axis.
    /// @return The vector of indices.
    inline auto
    col_indices() const -> std::vector<size_t>
    {
        return _col_indices;
    }
    
    /// @brief Gets the pointer to slice of tensor data.
    /// @param index The index of tensor slice.
    /// @return The pointer to slice of tensor.
    inline auto
    data(const size_t index) -> T*
    {
        return _data[index].data();
    }
    
    /// @brief Gets the constant pointer to slice of tensor data.
    /// @param index The index of tensor slice.
    /// @return The constant pointer to slice of tensor.
    inline auto
    data(const size_t index) const -> const T *
    {
        return _data[index].data();
    }
    
    /// @brief Gets tensor width along y axes.
    /// @return The width of tensor along  y axes.
    inline auto
    number_of_rows() const -> size_t
    {
        return _row_indices.size();
    }
    
    /// @brief Gets tensor width along y axes.
    /// @return The width of tensor along  z axes.
    inline auto
    number_of_columns() const -> size_t
    {
        return _col_indices.size();
    }
    
    /// @brief Gets tensor width along x axis.
    /// @return The width of tensor along x axis.
    inline auto
    aux_width() const -> size_t
    {
        return _indices.size();
    }
    
    /// @brief Gets number of blocks in tensor  width along x axis.
    /// @return The number of blocks in tensor.
    inline auto
    aux_blocks() const -> size_t
    {
        return _data.size() / _indices.size();
    }
    
   private:
    
    /// @brief Memory block for data storage of tensor slices along flatenned symmetrizes y,z axis.
    std::vector<std::vector<T>> _data;

    /// @brief The indices of compound tensor along x axis.
    std::vector<size_t> _indices;
    
    /// @brief The indices of compound tensor along y axis.
    std::vector<size_t> _row_indices;
    
    /// @brief The indices of compound tensor along z axis.
    std::vector<size_t> _col_indices;
};


#endif /* T3RectFlatBuffer_hpp */
