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
