#ifndef T3FlatBuffer_hpp
#define T3FlatBuffer_hpp

#include <cstddef>
#include <vector>
#include <ranges>

/// @brief Class CT3FlatBuffer stores general semi-symmetric rank-3 tensor as 2D flattened data structure.
template <typename T>
class CT3FlatBuffer
{
   public:
    
    /// @brief The default constructor.
    /// @param other The SIMD array to be copied.
    CT3FlatBuffer() = default;
    
    /// @brief The default constructor.
    /// @param indices The vector of indices along x axis of tensor.
    /// @param width The width of tensor along  y,z axes.
    CT3FlatBuffer(const std::vector<size_t>& indices, const size_t width)
    {
        _indices = indices;
        
        _width = width;
        
        _data.reserve(_indices.size());
        
        if (const auto nelems = _width * (_width + 1) / 2;  nelems > 0)
        {
            std::ranges::for_each(_indices, [&](const auto& index) {
                _data.push_back(std::vector<T>(nelems, T{0.0}));
            });
        }
    }
    
    /// @brief The default constructor.
    /// @param indices The vector of indices along x axis of tensor.
    /// @param width The width of tensor along  y,z axes.
    /// @param nbatches The number of batches.
    CT3FlatBuffer(const std::vector<size_t>& indices, const size_t width, const size_t nbatches)
    {
        _indices = indices;
        
        _width = width;
        
        _data.reserve(_indices.size() * nbatches);
        
        if (const auto nelems = _width * (_width + 1) / 2;  nelems > 0)
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
    CT3FlatBuffer(const CT3FlatBuffer &other) = default;

    /// @brief The default move constructor.
    /// @param other The SIMD array to be moved.
    CT3FlatBuffer(CT3FlatBuffer &&other) noexcept = default;

    /// @brief The custom destructor.
    ~CT3FlatBuffer() = default;
    
    /// @brief The default copy assignment operator.
    /// @param other The SIMD array to be copy assigned.
    /// @return The assigned SIMD array.
    auto operator=(const CT3FlatBuffer &other) -> CT3FlatBuffer & = default;

    /// @brief The default move assignment operator.
    /// @param other The SIMD array to be move assigned.
    /// @return The assigned SIMD array.
    auto operator=(CT3FlatBuffer &&other) noexcept -> CT3FlatBuffer & = default;

    /// @brief The equality operator.
    /// @param other The SIMD array to be compared.
    /// @return True if SIMD arrays are equal, False otherwise.
    auto operator==(const CT3FlatBuffer &other) const -> bool = default;

    /// @brief The non-equality operator.
    /// @param other The SIMD array to be compared.
    /// @return True if SIMD arrays are not equal, False otherwise.
    auto operator!=(const CT3FlatBuffer &other) const -> bool = default;
    
    /// @brief The gets vector of indices along x axis.
    /// @return The vector of indices.
    inline auto
    indices() const -> std::vector<size_t>
    {
        return _indices;
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
    
    /// @brief Gets tensor width along y,z axes.
    /// @return The width of tensor along  y,z axes.
    inline auto
    width() const -> size_t
    {
        return _width;
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
    
    /// @brief The width of tensor along y,z axes.
    size_t _width;
};

#endif /* T3FlatBuffer_hpp */
