//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef T3FlatBuffer_hpp
#define T3FlatBuffer_hpp

#include <algorithm>
#include <cstddef>
#include <vector>
#include <map>
#include <ranges>
#include <array>

#include "SubMatrix.hpp"
#include "MathFunc.hpp"

#include <iostream>

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
        
        _reduced = false;
        
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
    /// @param dimensions The dimensions of reduced tensor along  y,z axes.
    CT3FlatBuffer(const std::vector<size_t>& indices, const std::array<size_t, 4>& dimensions)
    {
        _indices = indices;
        
        _width = dimensions[2] + dimensions[3];
        
        _reduced = true;
        
        _data.reserve(_indices.size());
        
        if (_width > 0)
        {
            std::ranges::for_each(_indices, [&](const auto& index) {
                _data.push_back(std::vector<T>(_width, T{0.0}));
            });
        }
    }
    
    /// @brief The default constructor.
    /// @param mask_indices The map of indices along x axis of tensor.
    /// @param width The width of tensor along  y,z axes.
    CT3FlatBuffer(const std::map<size_t, size_t>& mask_indices, const size_t width)
    {
        _mask_indices = mask_indices;
        
        _width = width;
        
        _reduced = false;
        
        _data.reserve(_mask_indices.size());
        
        if (const auto nelems = _width * (_width + 1) / 2;  nelems > 0)
        {
            std::ranges::for_each(_mask_indices, [&](const auto& index) {
                _data.push_back(std::vector<T>(nelems, T{0.0}));
            });
        }
    }
    
    /// @brief The default constructor.
    /// @param mask_indices The map of indices along x axis of tensor.
    /// @param dimensions The dimensions of reduced tensor along  y,z axes.
    CT3FlatBuffer(const std::map<size_t, size_t>& mask_indices, const std::array<size_t, 4>& dimensions)
    {
        _mask_indices = mask_indices;
        
        _width = dimensions[2] + dimensions[3];
        
        _reduced = true;
        
        _data.reserve(_mask_indices.size());
        
        if (_width > 0)
        {
            std::ranges::for_each(_mask_indices, [&](const auto& index) {
                _data.push_back(std::vector<T>(_width, T{0.0}));
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
        
        _reduced = false;
        
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
    
    /// @brief The default constructor.
    /// @param mask_indices The mask of indices along x axis of tensor.
    /// @param width The width of tensor along  y,z axes.
    /// @param nbatches The number of batches.
    CT3FlatBuffer(const std::map<size_t, size_t>& mask_indices, const size_t width, const size_t nbatches)
    {
        _mask_indices = mask_indices;
        
        _width = width;
        
        _reduced = false; 
        
        _data.reserve(_mask_indices.size() * nbatches);
        
        if (const auto nelems = _width * (_width + 1) / 2;  nelems > 0)
        {
            std::ranges::for_each(std::views::iota(size_t{0}, nbatches), [&] (const auto index) {
                std::ranges::for_each(_mask_indices, [&](const auto& index) {
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
    
    /// @brief The gets map of masked indices along y axis.
    /// @return The map of masked indices.
    inline auto
    mask_indices() const -> std::map<size_t, size_t>
    {
        return _mask_indices;
    }
    
    /// @brief Gets the pointer to slice of tensor data.
    /// @param index The index of tensor slice.
    /// @return The pointer to slice of tensor.
    inline auto
    data(const size_t index) -> T*
    {
        return _data[index].data();
    }
    
    /// @brief Unpacks specific tensor slice into matrix.
    /// @param matrix The matrix to unpack tensor slice.
    /// @param index The index of tensor slice.
    auto
    unpack_data(CSubMatrix& matrix, const size_t index) const -> void
    {
        if (_width == matrix.number_of_rows())
        {
            auto ptr_data = _data[index].data();
         
            for (size_t i = 0; i < _width; i++)
            {
                for (size_t j = i; j < _width; j++)
                {
                    const auto fact = ptr_data[mathfunc::uplo_rm_index(i, j, _width)];
                    
                    matrix.at({i, j}) = fact;
                    
                    matrix.at({j, i}) = fact;
                }
            }
        }
    }
    
    /// @brief Unpacks specific tensor slice into matrix.
    /// @param matrix The matrix to unpack tensor slice.
    /// @param indices The vector of reduction indices.
    /// @param index The index of tensor slice.
    auto
    reduced_unpack_data(CSubMatrix& matrix, const std::vector<std::pair<size_t, size_t>>& indices, const size_t index) const -> void
    {
        matrix.zero();
        
        auto ptr_data = _data[index].data();
         
        size_t idx = 0;
            
        for (const auto& index : indices)
        {
            const auto fact = ptr_data[idx];
                
            matrix.at(index) = fact;
                
            matrix.at({index.second, index.first}) = fact;
                
            idx++;
        }
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
    
    /// @brief Gets number of elements in tensor slice along y,z axes.
    /// @return The width of tensor along  y,z axes.
    inline auto
    elements() const -> size_t
    {
        if (_reduced)
        {
            return _width;
        }
        else
        {
            return _width * (_width + 1) / 2;
        }
    }
    
    /// @brief Gets tensor width along x axis.
    /// @return The width of tensor along x axis.
    inline auto
    aux_width() const -> size_t
    {
        if (!_indices.empty())
        {
            return _indices.size();
        }
        else
        {
            return _mask_indices.size();
        }
        
    }
    
    /// @brief Gets number of blocks in tensor  width along x axis.
    /// @return The number of blocks in tensor.
    inline auto
    aux_blocks() const -> size_t
    {
        if (!_indices.empty())
        {
            return _data.size() / _indices.size();
        }
        else
        {
            return _data.size() / _mask_indices.size();
        }
    }
    
    /// @brief Checks if  tensor slice along y,z axes.
    /// @return The flag of tensor slice reduction along  y,z axes.
    inline auto
    is_reduced() const -> bool
    {
        return _reduced;
    }
    
   private:
    
    /// @brief Memory block for data storage of tensor slices along flatenned symmetrizes y,z axis.
    std::vector<std::vector<T>> _data;

    /// @brief The indices of compound tensor along x axis.
    std::vector<size_t> _indices;
    
    /// @brief The indices of compound tensor along z axis.
    std::map<size_t, size_t> _mask_indices;
    
    /// @brief The width of tensor along y,z axes.
    size_t _width;
    
    /// @brief The form of tensor alongf y,z axis.
    bool _reduced;
};

#endif /* T3FlatBuffer_hpp */
