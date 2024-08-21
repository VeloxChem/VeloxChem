#ifndef MolecularBasis_hpp
#define MolecularBasis_hpp

#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "AtomBasis.hpp"

/// @brief Class CMolecularBasis stores data about molecular basis and provides
/// set of methods for handling of molecular basis data.
class CMolecularBasis
{
   public:
    /// @brief The default constructor.
    CMolecularBasis();

    /// @brief The constructor witj vector of unique atom bases and vector of atom bases indices.
    /// @param basis_sets The vector of unique atom bases.
    /// @param indices The vector of atom bases indices.
    CMolecularBasis(const std::vector<CAtomBasis> &basis_sets, const std::vector<int> &indices);

    /// @brief The default copy constructor.
    /// @param other The molecular basis to be copied.
    CMolecularBasis(const CMolecularBasis &other);

    /// @brief The default move constructor.
    /// @param other The molecular basis to be moved.
    CMolecularBasis(CMolecularBasis &&other) noexcept;

    /// @brief The default destructor.
    ~CMolecularBasis() = default;

    /// @brief The default copy assignment operator.
    /// @param other The molecular basis to be copy assigned.
    /// @return The assigned molecular basis.
    auto operator=(const CMolecularBasis &other) -> CMolecularBasis &;

    /// @brief The default move assignment operator.
    /// @param other The molecular basis to be move assigned.
    /// @return The assigned molecular basis.
    auto operator=(CMolecularBasis &&other) noexcept -> CMolecularBasis &;

    /// @brief The equality operator.
    /// @param other The molecular basis to be compared.
    /// @return True if molecular bases are equal, False otherwise.
    auto operator==(const CMolecularBasis &other) const -> bool;

    /// @brief The equality operator.
    /// @param other The molecular basis to be compared.
    /// @return True if molecular bases are not equal, False otherwise.
    auto operator!=(const CMolecularBasis &other) const -> bool;

    /// @brief Adds atom basis to molecular basis.
    /// @param basis The atom basis to be added.
    auto add(const CAtomBasis &basis) -> void;

    /// @brief Reduces molecular basis to valence molecular basis.
    /// @return The valence molecular basis.
    auto reduce_to_valence_basis() const -> CMolecularBasis;

    /// @brief Slice fraction of molecular basis for specific atoms.
    /// @param atoms The vector of atoms.
    /// @return The molecular basis for selected atoms.
    auto slice(const std::vector<int> &atoms) const -> CMolecularBasis;

    /// @brief Gets vector of unique atom bases.
    /// @return The vector of atom bases.
    auto basis_sets() const -> std::vector<CAtomBasis>;

    /// @brief Gets vector of indices for atom bases.
    /// @return The vector of indices for atom bases.
    auto basis_sets_indices() const -> std::vector<int>;

    /// @brief Gets maximum angular momentum of molecular basis.
    /// @return The maximum angular momentum.
    auto max_angular_momentum() const -> int;

    /// @brief Gets maximum angular momentum of molecular basis for selected
    /// atoms.
    /// @param atoms The vector of atoms to select atom bases.
    /// @return The maximum angular momentum.
    auto max_angular_momentum(const std::vector<int> &atoms) const -> int;

    /// Gets vector of contracted GTOs.

    /// @brief Gets vector of basis functions.
    /// @return The vector of basis functions.
    auto basis_functions() const -> std::vector<CBasisFunction>;

    /// @brief Gets vector of basis functions with specific angular momentum in
    /// molecular basis.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @return The vector of basis functions.
    auto basis_functions(const int angular_momentum) const -> std::vector<CBasisFunction>;

    /// @brief Gets vector of basis functions with specific angular momentum and
    /// number of primitives in molecular basis.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @param npgtos The number of primitive basis functions in basis function.
    /// @return The vector of basis functions.
    auto basis_functions(const int angular_momentum, const size_t npgtos) const -> std::vector<CBasisFunction>;

    /// @brief Gets vector of basis functions in molecular basis for vector of
    /// specific atoms.
    /// @param atoms The vector of selected atoms.
    /// @return The vector of basis functions.
    auto basis_functions(const std::vector<int> &atoms) const -> std::vector<CBasisFunction>;

    /// @brief Gets vector of basis functions with specific angular momentum in
    /// molecular basis for specific atoms.
    /// @param atoms The vector of selected atoms.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @return The vector of basis functions.
    auto basis_functions(const std::vector<int> &atoms, const int angular_momentum) const -> std::vector<CBasisFunction>;

    /// @brief Gets vector of basis functions with specific angular momentum and
    /// number of primitives in molecular basis for specific atoms.
    /// @param atoms The vector of selected atoms.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @param npgtos The number of primitive basis functions in basis function.
    /// @return The vector of basis functions.
    auto basis_functions(const std::vector<int> &atoms, const int angular_momentum, const size_t npgtos) const -> std::vector<CBasisFunction>;

    /// @brief Gets vector of atomic indices.
    /// @return The vector of atomic indices.
    auto atomic_indices() const -> std::vector<int>;

    /// @brief Gets vector of atomic indices for basis functions of specific
    /// angular momentum in molecular basis.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @return The vector of atomic indices.
    auto atomic_indices(const int angular_momentum) const -> std::vector<int>;

    /// @brief Gets vector of atomic indices for basis functions with specific
    /// angular momentum and number of primitives in molecular basis.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @param npgtos The number of primitive basis functions in basis function.
    /// @return The vector of atomic indices.
    auto atomic_indices(const int angular_momentum, const size_t npgtos) const -> std::vector<int>;

    /// @brief Gets vector of atomic indices for selected atoms.
    /// @param atoms The vector of selected atoms.
    /// @return The vector of atomic indices.
    auto atomic_indices(const std::vector<int> &atoms) const -> std::vector<int>;

    /// @brief Gets vector of atomic indices for basis functions of specific
    /// angular momentum in molecular basis for selected atoms.
    /// @param atoms The vector of selected atoms.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @return The vector of atomic indices.
    auto atomic_indices(const std::vector<int> &atoms, const int angular_momentum) const -> std::vector<int>;

    /// @brief Gets vector of atomic indices for basis functions of specific
    /// angular momentum and of number of primitives in molecular basis.
    /// @param atoms The vector of selected atoms.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @param npgtos The number of primitive basis functions in basis function.
    /// @return The vector of atomic indices.
    auto atomic_indices(const std::vector<int> &atoms, const int angular_momentum, const size_t npgtos) const -> std::vector<int>;

    /// @brief Determines number of basis functions with specific angular
    /// momentum in molecular basis.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @return The number of basis functions.
    auto number_of_basis_functions(const int angular_momentum) const -> size_t;

    /// @brief Determines number of basis functions with specific angular
    /// momentum and number of primitives in molecular basis.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @param npgtos The number of primitive basis functions in basis function.
    /// @return The number of basis functions.
    auto number_of_basis_functions(const int angular_momentum, const size_t npgtos) const -> size_t;

    /// @brief Determines number of basis functions with specific angular
    /// momentum in molecular basis for selected atoms.
    /// @param atoms The vector of selected atoms.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @return The number of basis functions.
    auto number_of_basis_functions(const std::vector<int> &atoms, const int angular_momentum) const -> size_t;

    /// @brief Determines number of basis functions with specific angular
    /// momentum and number of primitives in molecular basis for selected atoms.
    /// @param atoms The vector of selected atoms.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @param npgtos The number of primitive basis functions in basis function.
    /// @return The number of basis functions.
    auto number_of_basis_functions(const std::vector<int> &atoms, const int angular_momentum, const size_t npgtos) const -> size_t;

    /// @brief Determines number of primitive basis functions with specific
    /// angular momentum in molecular basis.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @return The number of primitive basis functions.
    auto number_of_primitive_functions(const int angular_momentum) const -> size_t;

    /// @brief Determines number of primitive basis functions with specific
    /// angular momentum in molecular basis for selected atoms.
    /// @param atoms The vector of selected atoms.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @return The number of primitive basis functions.
    auto number_of_primitive_functions(const std::vector<int> &atoms, const int angular_momentum) const -> size_t;

    /// @brief Gets set of unique contraction numbers of basis functions with
    /// given angular momentum in molecular basis.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @return The set of unique contraction numbers of basis functions.
    auto contraction_depths(const int angular_momentum) const -> std::set<size_t>;

    /// @brief Gets set of unique contraction numbers of basis functions with
    /// given angular momentum in molecular basis for selected atoms.
    /// @param atoms The vector of selected atoms.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @return The set of unique contraction numbers of basis functions.
    auto contraction_depths(const std::vector<int> &atoms, const int angular_momentum) const -> std::set<size_t>;

    /// @brief Determines size of atomic orbitals basis.
    /// @return The number of atomic orbitals.
    auto dimensions_of_basis() const -> size_t;

    /// @brief Determines partial size up to specific angular momentum of atomic
    /// orbitals basis.
    /// @param angular_momentum The angular momentum upper bound of atomic
    /// arbitals.
    /// @return The number of atomic orbitals.
    auto dimensions_of_basis(const int angular_momentum) const -> size_t;

    /// @brief Determines size of primitive atomic orbitals basis.
    /// @return The number of primitive atomic orbitals.
    auto dimensions_of_primitive_basis() const -> size_t;

    /// @brief Creates map indices for basis functions with specific angular
    /// momentum and number of primitives in molecular basis.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @param npgtos The number of primitive basis functions in basis function.
    /// @return The vector of basis function indices.
    auto index_map(const int angular_momentum, const size_t npgtos) const -> std::vector<size_t>;

    /// @brief Creates map indices for basis functions with specific angular
    /// momentum and number of primitives in molecular basis for selected atoms.
    /// @param atoms The vector of selected atoms.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @param npgtos The number of primitive basis functions in basis function.
    /// @return The vector of basis function indices.
    auto index_map(const std::vector<int> &atoms, const int angular_momentum, const size_t npgtos) const -> std::vector<size_t>;

    /// @brief Gets main atom basis label in molecular basis.
    /// @return The main atom basis label.
    auto main_basis_label() const -> std::string;

   private:
    /// @brief The vector of atom basis sets.
    std::vector<CAtomBasis> _basis_sets;

    /// @brief The vector of atom basis sets indices.
    std::vector<int> _indices;

    /// @brief Gets atom basis labels frequency map.
    /// @return The labels frequency map.
    auto _labels_frequency_map() const -> std::unordered_map<std::string, int>;
};

#endif /* MolecularBasis_hpp */
