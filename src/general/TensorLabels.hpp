#ifndef TensorLabels_hpp
#define TensorLabels_hpp

#include <string>
#include <vector>

namespace tensor {  // tensor

/// @brief Generates all Cartesian labels for tensor.
/// @param order The order of tensor.
/// @return The vector of Cartesian tensor labels.
auto cartesian_labels(const int order) -> std::vector<std::string>;

/// @brief Generates all spherical labels for tensor.
/// @param order The order of tensor.
/// @return The vector of spherical tensor labels.
auto spherical_labels(const int order) -> std::vector<std::string>;

/// @brief Gets Cartesian index of canonical tensor component.
/// @param label The label of Cartesian tensor component.
/// @return The index of Cartesian tensor component.
auto cartesian_index(const std::string& label) -> int;

/// @brief Gets label of canonical tensor.
/// @param order The order of canonical tensor.
/// @return The label of canonical tensor.
auto label(const int order) -> char;

/// @brief Gets order of canonical tensor.
/// @param label The uppercased label of canonical tensor.
/// @return The order of canonical tensor.
auto order(const char label) -> int;

}  // namespace tensor

#endif /* TensorLabels_hpp */
