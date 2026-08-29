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


#ifndef SparseTensor_hpp
#define SparseTensor_hpp

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "AtomBasisGroup.hpp"
#include "AtomBasisPairGroup.hpp"
#include "AtomBasisTripleSparsity.hpp"
#include "ErrorHandler.hpp"
#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "ValuesState.hpp"

/// @brief Class CSparseTensor stores the sparsity pattern of a three-dimensional
/// tensor in atomic orbital basis, as the sparsity patterns of its partially
/// sparse blocks. Each block is described by CAtomBasisTripleSparsity, which is
/// sparse along the a and b sides and dense along the c side. The atom pairs of
/// the diagonal atom pairs on the a and b sides are part of those blocks, so no
/// separate diagonal blocks are needed.
/// @note The type of tensor refers to the a and b sides only, as the c side is
/// never interchangeable with them. A symmetric or antisymmetric tensor has a
/// single molecular basis on the a and b sides, while a general tensor has a
/// molecular basis on each of them.
/// @note The blocks without atom pairs are dropped, so a block index is not an
/// index of the atom basis pair group and atom basis group it originates from.
/// @note The values blocks are not allocated by the constructors, which set up
/// the sparsity patterns only. The values blocks are addressed as one pointer per
/// block, and the state of the values blocks is reported by get_values_state().
class CSparseTensor
{
   public:
    /// @brief The default constructor.
    CSparseTensor()

        : _blocks{}

        , _values{}

        , _type(mat_t::general)

        , _values_state(valstat::empty)
    {
    }

    /// @brief The constructor with sparsity patterns of the blocks and tensor type.
    /// @param blocks The sparsity patterns of the blocks.
    /// @param mat_type The type of tensor.
    CSparseTensor(const std::vector<CAtomBasisTripleSparsity> &blocks, const mat_t mat_type)

        : _blocks(blocks)

        , _values(blocks.size(), nullptr)

        , _type(mat_type)

        , _values_state(valstat::empty)
    {
    }

    /// @brief The constructor with molecule, molecular basis on a and b sides,
    /// molecular basis on c side, integral screener and screening threshold.
    /// @param molecule The molecule to compute interatomic distances from.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The molecular basis on c side.
    /// @param screener The integral bound, evaluated as screener(a_function,
    /// b_function, c_function, distance).
    /// @param threshold The screening threshold.
    /// @param mat_type The type of tensor.
    template <typename B>
    CSparseTensor(const CMolecule       &molecule,
                  const CMolecularBasis &basis,
                  const CMolecularBasis &aux_basis,
                  const B               &screener,
                  const double           threshold,
                  const mat_t            mat_type)

        : _blocks{}

        , _values{}

        , _type(mat_type)

        , _values_state(valstat::empty)
    {
        // NOTE: a symmetric or antisymmetric tensor needs the upper triangle of
        // the atom basis pair groups only, while a general tensor needs their
        // full direct product.

        auto groups = (mat_type == mat_t::general) ? basis.basis_pair_groups(basis) : basis.basis_pair_groups();

        _add_blocks(molecule, groups, aux_basis.basis_groups(), screener, threshold);
    }

    /// @brief The constructor with molecule, molecular bases on a and b sides,
    /// molecular basis on c side, integral screener and screening threshold.
    /// @param molecule The molecule to compute interatomic distances from.
    /// @param bra_basis The molecular basis on a side.
    /// @param ket_basis The molecular basis on b side.
    /// @param aux_basis The molecular basis on c side.
    /// @param screener The integral bound.
    /// @param threshold The screening threshold.
    /// @param mat_type The type of tensor, which must be general.
    template <typename B>
    CSparseTensor(const CMolecule       &molecule,
                  const CMolecularBasis &bra_basis,
                  const CMolecularBasis &ket_basis,
                  const CMolecularBasis &aux_basis,
                  const B               &screener,
                  const double           threshold,
                  const mat_t            mat_type)

        : _blocks{}

        , _values{}

        , _type(mat_type)

        , _values_state(valstat::empty)
    {
        errors::assertMsgCritical(mat_type == mat_t::general,
                                  std::string("SparseTensor: Tensor with two molecular bases must be of general type"));

        auto groups = bra_basis.basis_pair_groups(ket_basis);

        _add_blocks(molecule, groups, aux_basis.basis_groups(), screener, threshold);
    }

    /// @brief The copy constructor.
    /// @param other The sparse tensor to be copied.
    CSparseTensor(const CSparseTensor &other)

        : _blocks(other._blocks)

        , _values(other._values.size(), nullptr)

        , _type(other._type)

        , _values_state(other._values_state)
    {
        if (_values_state == valstat::allocated) _copy_values(other);
    }

    /// @brief The move constructor.
    /// @param other The sparse tensor to be moved.
    CSparseTensor(CSparseTensor &&other) noexcept

        : _blocks(std::move(other._blocks))

        , _values(std::move(other._values))

        , _type(other._type)

        , _values_state(other._values_state)
    {
        other._values.clear();

        other._values_state = valstat::empty;
    }

    /// @brief The destructor.
    ~CSparseTensor()
    {
        _deallocate();
    }

    /// @brief The copy assignment operator.
    /// @param other The sparse tensor to be copy assigned.
    /// @return The assigned sparse tensor.
    auto
    operator=(const CSparseTensor &other) -> CSparseTensor &
    {
        if (this != &other)
        {
            _deallocate();

            _blocks = other._blocks;

            _values.assign(other._values.size(), nullptr);

            _type = other._type;

            _values_state = other._values_state;

            if (_values_state == valstat::allocated) _copy_values(other);
        }

        return *this;
    }

    /// @brief The move assignment operator.
    /// @param other The sparse tensor to be move assigned.
    /// @return The assigned sparse tensor.
    auto
    operator=(CSparseTensor &&other) noexcept -> CSparseTensor &
    {
        if (this != &other)
        {
            _deallocate();

            _blocks = std::move(other._blocks);

            _values = std::move(other._values);

            _type = other._type;

            _values_state = other._values_state;

            other._values.clear();

            other._values_state = valstat::empty;
        }

        return *this;
    }

    /// @brief Allocates the values blocks of tensor, leaving their content
    /// undefined. Each values block is allocated separately, so no single
    /// allocation spans the whole tensor.
    auto
    allocate() -> void
    {
        if (_values_state == valstat::allocated) return;

        std::ranges::for_each(std::views::iota(size_t{0}, _values.size()),
                              [&](const auto i) { _values[i] = new double[number_of_elements(i)]; });

        _values_state = valstat::allocated;
    }

    /// @brief Deallocates the values blocks of tensor, leaving its sparsity
    /// patterns intact.
    auto
    deallocate() -> void
    {
        _deallocate();
    }

    /// @brief Sets the values of all values blocks of tensor to zero.
    auto
    zero() -> void
    {
        errors::assertMsgCritical(_values_state == valstat::allocated,
                                  std::string("SparseTensor.zero: Values blocks of tensor are not allocated"));

        std::ranges::for_each(std::views::iota(size_t{0}, _values.size()),
                              [&](const auto i) { std::fill(_values[i], _values[i] + number_of_elements(i), 0.0); });
    }

    /// @brief Scales the values of all values blocks of tensor.
    /// @param factor The scaling factor.
    auto
    scale(const double factor) -> void
    {
        errors::assertMsgCritical(_values_state == valstat::allocated,
                                  std::string("SparseTensor.scale: Values blocks of tensor are not allocated"));

        std::ranges::for_each(std::views::iota(size_t{0}, _values.size()), [&](const auto i) {
            std::ranges::for_each(std::views::iota(size_t{0}, number_of_elements(i)), [&](const auto k) { _values[i][k] *= factor; });
        });
    }

    /// @brief Sets type of tensor.
    /// @param mat_type The type of tensor.
    auto
    set_type(const mat_t mat_type) -> void
    {
        _type = mat_type;
    }

    /// @brief Gets type of tensor.
    /// @return The type of tensor.
    auto
    get_type() const -> mat_t
    {
        return _type;
    }

    /// @brief Gets state of the values blocks of tensor.
    /// @return The state of the values blocks.
    auto
    get_values_state() const -> valstat
    {
        return _values_state;
    }

    /// @brief Gets sparsity patterns of the blocks.
    /// @return The constant reference to vector of sparsity patterns.
    auto
    blocks() const -> const std::vector<CAtomBasisTripleSparsity> &
    {
        return _blocks;
    }

    /// @brief Gets sparsity pattern of the block with specific index.
    /// @param index The index of block.
    /// @return The constant reference to sparsity pattern.
    auto
    block(const size_t index) const -> const CAtomBasisTripleSparsity &
    {
        errors::assertMsgCritical(index < _blocks.size(), std::string("SparseTensor.block: Index of block is out of range"));

        return _blocks[index];
    }

    /// @brief Gets number of blocks in tensor.
    /// @return The number of blocks.
    auto
    number_of_blocks() const -> size_t
    {
        return _blocks.size();
    }

    /// @brief Gets values of the block with specific index.
    /// @param index The index of block.
    /// @return The pointer to the values of the block.
    auto
    values(const size_t index) -> double *
    {
        _check_values(index, std::string("values"));

        return _values[index];
    }

    /// @brief Gets values of the block with specific index.
    /// @param index The index of block.
    /// @return The constant pointer to the values of the block.
    auto
    values(const size_t index) const -> const double *
    {
        _check_values(index, std::string("values"));

        return _values[index];
    }

    /// @brief Gets values of specific combination of basis functions in the
    /// block with specific index.
    /// @param index The index of block.
    /// @param a_angular_momentum The angular momentum of basis function on a side.
    /// @param a_index The index of basis function on a side.
    /// @param b_angular_momentum The angular momentum of basis function on b side.
    /// @param b_index The index of basis function on b side.
    /// @param c_angular_momentum The angular momentum of basis function on c side.
    /// @param c_index The index of basis function on c side.
    /// @return The pointer to the values of the combination of basis functions.
    auto
    values(const size_t index,
           const int    a_angular_momentum,
           const size_t a_index,
           const int    b_angular_momentum,
           const size_t b_index,
           const int    c_angular_momentum,
           const size_t c_index) -> double *
    {
        return values(index) + _blocks[index].element_offset(a_angular_momentum, a_index, b_angular_momentum, b_index,
                                                             c_angular_momentum, c_index);
    }

    /// @brief Gets values of specific combination of basis functions in the
    /// block with specific index.
    /// @param index The index of block.
    /// @param a_angular_momentum The angular momentum of basis function on a side.
    /// @param a_index The index of basis function on a side.
    /// @param b_angular_momentum The angular momentum of basis function on b side.
    /// @param b_index The index of basis function on b side.
    /// @param c_angular_momentum The angular momentum of basis function on c side.
    /// @param c_index The index of basis function on c side.
    /// @return The constant pointer to the values of the combination of basis functions.
    auto
    values(const size_t index,
           const int    a_angular_momentum,
           const size_t a_index,
           const int    b_angular_momentum,
           const size_t b_index,
           const int    c_angular_momentum,
           const size_t c_index) const -> const double *
    {
        return values(index) + _blocks[index].element_offset(a_angular_momentum, a_index, b_angular_momentum, b_index,
                                                             c_angular_momentum, c_index);
    }

    /// @brief Gets number of values required to store the integrals of specific
    /// block.
    /// @param index The index of block.
    /// @return The number of values.
    auto
    number_of_elements(const size_t index) const -> size_t
    {
        errors::assertMsgCritical(index < _blocks.size(), std::string("SparseTensor.number_of_elements: Index of block is out of range"));

        return _blocks[index].number_of_elements();
    }

    /// @brief Gets number of values required to store the integrals of all blocks.
    /// @return The number of values.
    auto
    number_of_elements() const -> size_t
    {
        size_t nvals = 0;

        std::ranges::for_each(_blocks, [&](const auto &block) { nvals += block.number_of_elements(); });

        return nvals;
    }

    /// @brief Gets memory required to store the values blocks of tensor.
    /// @return The memory in bytes.
    /// @note Only the values blocks are counted, as they are what allocate()
    /// reserves. The sparsity patterns describing the blocks are not included.
    auto
    memory_size() const -> size_t
    {
        return number_of_elements() * sizeof(double);
    }

   private:
    /// @brief Checks that the values blocks are allocated and that a block index
    /// is in range.
    /// @param index The index of block.
    /// @param label The name of the accessor requesting the check.
    auto
    _check_values(const size_t index, const std::string &label) const -> void
    {
        errors::assertMsgCritical(_values_state == valstat::allocated,
                                  std::string("SparseTensor.") + label + std::string(": Values blocks of tensor are not allocated"));

        errors::assertMsgCritical(index < _blocks.size(),
                                  std::string("SparseTensor.") + label + std::string(": Index of block is out of range"));
    }

    /// @brief Deallocates the values blocks of tensor.
    auto
    _deallocate() -> void
    {
        std::ranges::for_each(_values, [](double *block) { delete[] block; });

        std::ranges::fill(_values, nullptr);

        _values_state = valstat::empty;
    }

    /// @brief Allocates the values blocks of tensor and copies the values of
    /// another tensor into them.
    /// @param other The sparse tensor to copy the values of.
    auto
    _copy_values(const CSparseTensor &other) -> void
    {
        std::ranges::for_each(std::views::iota(size_t{0}, _values.size()), [&](const auto i) {
            const auto nvals = number_of_elements(i);

            _values[i] = new double[nvals];

            std::copy(other._values[i], other._values[i] + nvals, _values[i]);
        });
    }

    /// @brief Adds the sparsity patterns of the non-empty blocks of the atom
    /// basis pair groups and atom basis groups.
    /// @param molecule The molecule to compute interatomic distances from.
    /// @param groups The atom basis pair groups on a and b sides.
    /// @param aux_groups The atom basis groups on c side.
    /// @param screener The integral bound.
    /// @param threshold The screening threshold.
    template <typename B>
    auto
    _add_blocks(const CMolecule                     &molecule,
                std::vector<CAtomBasisPairGroup>    &groups,
                const std::vector<CAtomBasisGroup>  &aux_groups,
                const B                             &screener,
                const double                         threshold) -> void
    {
        // NOTE: the atom basis pair groups are independent, so their atom pairs
        // are ordered by interatomic distance in parallel. Dynamic scheduling is
        // used as the groups differ widely in the number of atom pairs.

        const auto ngroups = static_cast<int>(groups.size());

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < ngroups; i++)
        {
            groups[i].sort_by_distance(molecule);
        }

        // NOTE: the blocks are added in the order of the groups, so that the
        // layout of the values blocks does not depend on the scheduling above.

        std::ranges::for_each(groups, [&](const auto &group) {
            std::ranges::for_each(aux_groups, [&](const auto &aux_group) {
                const auto block = CAtomBasisTripleSparsity(group, aux_group, screener, threshold);

                if (block.number_of_pairs() > 0) _blocks.push_back(block);
            });
        });

        _values.assign(_blocks.size(), nullptr);
    }

    /// @brief The sparsity patterns of the blocks.
    std::vector<CAtomBasisTripleSparsity> _blocks;

    /// @brief The values blocks of tensor, one pointer per block. The pointers
    /// are null while the values blocks are not allocated.
    std::vector<double *> _values;

    /// @brief The type of tensor.
    mat_t _type;

    /// @brief The state of the values blocks of tensor.
    valstat _values_state;
};

#endif /* SparseTensor_hpp */
