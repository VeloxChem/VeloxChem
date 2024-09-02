#ifndef AtomBasis_hpp
#define AtomBasis_hpp

#include <set>
#include <string>
#include <vector>

#include "BasisFunction.hpp"

/// @brief Class CAtomBasis stores data about atom basis and provides set of
/// methods for handling of atom basis data.
class CAtomBasis
{
   public:
    /// @brief The default constructor.
    CAtomBasis();

    /// @brief The constructor with vector of basis functions, name of basis set,
    /// ECP label, and chemical element identifier.
    /// @param functions The vector of basis functions.
    /// @param name The name of atom basis.
    /// @param ecp_label The label of ECP in atom basis.
    /// @param identifier The chemical element identifier.
    CAtomBasis(const std::vector<CBasisFunction> &functions, const std::string &name, const std::string &ecp_label, const int identifier);

    /// @brief The default copy constructor.
    /// @param other The atom basis to be copied.
    CAtomBasis(const CAtomBasis &other);

    /// @brief The default move constructor.
    /// @param other The atom basis to be moved.
    CAtomBasis(CAtomBasis &&other) noexcept;

    /// @brief The default destructor.
    ~CAtomBasis() = default;

    /// @brief The default copy assignment operator.
    /// @param other The atom basis to be copy assigned.
    /// @return The assigned atom basis.
    auto operator=(const CAtomBasis &other) -> CAtomBasis &;

    /// @brief The default move assignment operator.
    /// @param other The atom basis to be move assigned.
    /// @return The assigned atom basis.
    auto operator=(CAtomBasis &&other) noexcept -> CAtomBasis &;

    /// @brief The equality operator.
    /// @param other The atom basis to be compared.
    /// @return True if atom bases are equal, False otherwise.
    auto operator==(const CAtomBasis &other) const -> bool;

    /// @brief The equality operator.
    /// @param other The atom basis to be compared.
    /// @return True if atom bases are not equal, False otherwise.
    auto operator!=(const CAtomBasis &other) const -> bool;

    /// @brief Sets identifier of chemical element in atom basis.
    /// @param identifier The identifier of chemical element.
    auto set_identifier(const int identifier) -> void;

    /// @brief Sets name of atom basis.
    /// @param name The name of atom basis.
    auto set_name(const std::string &name) -> void;

    /// @brief Sets effective core potential label of atom basis.
    /// @param label The label of effective core potential in atom basis.
    auto set_ecp_label(const std::string &label) -> void;

    /// @brief Adds basis function to atom basis.
    /// @param function The basis function.
    auto add(const CBasisFunction &function) -> void;

    /// @brief Reduces atom basis set to valence atom basis.
    /// @return The valence atom basis.
    auto reduce_to_valence_basis() const -> CAtomBasis;

    /// @brief Gets vector of basis functions.
    /// @return The vector of basis functions.
    auto basis_functions() const -> std::vector<CBasisFunction>;

    /// @brief Gets vector of basis functions with specific angular momentum.
    /// @param angular_momentum The angular momentum of requested basis
    /// functions.
    /// @return The vector of basis functions.
    auto basis_functions(const int angular_momentum) const -> std::vector<CBasisFunction>;

    /// @brief Gets vector of GTOs with specific angular momentum and number of
    /// primitive basis functions.
    /// @param angular_momentum The angular momentum of requested basis
    /// functions.
    /// @param npgtos The number primitive basis functions in requested basis
    /// functions.
    /// @return The vector of basis functions.
    auto basis_functions(const int angular_momentum, const size_t npgtos) const -> std::vector<CBasisFunction>;

    /// @brief Gets identifier of chemical element.
    /// @return The identifier of chemical element.
    auto get_identifier() const -> int;

    /// @brief Gets name of atom basis.
    /// @return The name of atom basis.
    auto get_name() const -> std::string;

    /// @brief Gets effective core potential label of atom basis.
    /// @return The label of effective core potential.
    auto get_ecp_label() const -> std::string;

    /// @brief Checks if atom basis requires effective core potential.
    /// @return Trrue if atom basis contains effective core potential, False
    /// otherwise.
    auto need_ecp() const -> bool;

    /// @brief Gets maximum angular momentum of basis functions in atom basis.
    /// @return The maximum angular momentum.
    auto max_angular_momentum() const -> int;

    /// @brief Gets number of basis functions with specific angular momentum in
    /// atom basis.
    /// @param angular_momentum The requested angular momentum of basis
    /// functions.
    /// @return The number of basis functions.
    auto number_of_basis_functions(const int angular_momentum) const -> size_t;

    /// @brief Gets number of basis functions with specific angular momentum and
    /// number of prmitive basis functions in atom basis.
    /// @param angular_momentum The requested angular momentum of basis
    /// functiions.
    /// @param npgtos The requested number of primitive basis functions in basis
    /// functions.
    /// @return The number of basis functions.
    auto number_of_basis_functions(const int angular_momentum, const size_t npgtos) const -> size_t;

    /// @brief Gets number of primitive basis functions with requested angular
    /// momentum.
    /// @param angular_momentum The requested angular momentum of primitive basis
    /// functiions.
    /// @return The number of primitive basis functions.
    auto number_of_primitive_functions(const int angular_momentum) const -> size_t;

    /// @brief Get set of unique contraction numbers of basis functions with
    /// given angular momentum in atom basis.
    /// @param angular_momentum The requested angular momentum of basis
    /// functions.
    /// @return The set of unique contraction numbers
    auto contraction_depths(const int angular_momentum) const -> std::set<size_t>;

    /// @brief Gets contraction string of atom basis.
    /// @return The contraction string.
    auto contraction_string() const -> std::string;

    /// @brief Gets primitive basis functions strig for atom basis.
    /// @return The primitive basis functions string.
    auto primitives_string() const -> std::string;

   private:
    /// @brief The vector of basis functions.
    std::vector<CBasisFunction> _functions;

    /// @brief The name of atomic basis.
    std::string _name;

    /// @brief The effective core potential label of atomic basis.
    std::string _ecp_label;

    /// @brief  identifier of chemical element.
    int _identifier;
};

#endif /* AtomBasis_hpp */
