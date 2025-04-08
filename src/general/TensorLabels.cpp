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

#include "TensorLabels.hpp"

#include <algorithm>
#include <ranges>
#include <string>
#include <vector>

#include "TensorComponents.hpp"

namespace tensor {

auto
cartesian_labels(const int order) -> std::vector<std::string>
{
    if (order > 0)
    {
        std::vector<std::string> labels = {"X", "Y", "Z"};

        std::ranges::for_each(std::views::iota(1, order), [&](const int i) {
            std::vector<std::string> new_labels;
            std::ranges::for_each(labels, [&](const auto& label) { new_labels.push_back("X" + label); });
            std::ranges::for_each(labels, [&](const auto& label) {
                if (label[0] != 'X') new_labels.push_back("Y" + label);
            });
            new_labels.push_back("Z" + labels.back());
            labels = new_labels;
        });

        return labels;
    }
    else
    {
        return std::vector<std::string>();
    }
}

auto
spherical_labels(const int order) -> std::vector<std::string>
{
    const auto label = std::string(1, tensor::label(order));

    if (order == 0)
    {
        return {
            label,
        };
    }
    else if (order > 0)
    {
        std::vector<std::string> labels;

        const auto tcomps = tensor::number_of_spherical_components(std::array<int, 1>{order});

        std::ranges::for_each(std::views::iota(0, tcomps), [&](const int i) {
            if (i > order)
            {
                labels.push_back(label + "+" + std::to_string(i - order));
            }
            else
            {
                labels.push_back(label + std::to_string(i - order));
            }
        });

        return labels;
    }
    else
    {
        return std::vector<std::string>();
    }
}

auto
cartesian_index(const std::string& label) -> int
{
    if (auto order = static_cast<int>(label.size()); order > 0)
    {
        auto labels = tensor::cartesian_labels(order);

        auto pos = std::ranges::find(labels, label);

        return (pos != labels.end()) ? static_cast<int>(std::distance(labels.begin(), pos)) : -1;
    }
    else
    {
        return -1;
    }
}

auto
label(const int order) -> char
{
    const auto labels = std::string("SPDFGHIKLMNOQRTUV");

    return labels.at(order);
}

auto
order(const char label) -> int
{
    const auto labels = std::string("SPDFGHIKLMNOQRTUV");

    auto pos = labels.find_first_of(label);

    return (pos == std::string::npos) ? -1 : static_cast<int>(pos);
}

}  // namespace tensor
