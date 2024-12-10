#ifndef AtomCorePotential_hpp
#define AtomCorePotential_hpp

#include <cstddef>
#include <vector>

/// @brief Class CAtomCorePotential stores data about atom core potential and
/// provides set of methods for handling of base core potential data.
class CAtomCorePotential
{
   public:
    /// @brief The default constructor.
    CAtomCorePotential();

    /// @brief The constructor with exponents, expansion factors, and radial orders.
    /// @param exponents The vector of exponents of primitive local potentials.
    /// @param factors The vector of expansion factors of primitive local potential.
    /// @param radial_orders The vector of radial orders of primitive local potentials.
    CAtomCorePotential(const std::vector<double> &exponents,
                       const std::vector<double> &factors,
                       const std::vector<int>    &radial_orders);

    /// @brief The default copy constructor.
    /// @param other The base core potential to be copied.
    CAtomCorePotential(const CAtomCorePotential &other);

    /// @brief The default move constructor.
    /// @param other The base core potential to be moved.
    CAtomCorePotential(CAtomCorePotential &&other) noexcept;

    /// @brief The default destructor.
    ~CAtomCorePotential() = default;

    /// @brief The default copy assignment operator.
    /// @param other The base core potential to be copy assigned.
    /// @return The assigned base core potential.
    auto operator=(const CAtomCorePotential &other) -> CAtomCorePotential&;

    /// @brief The default move assignment operator.
    /// @param other The base core potential to be move assigned.
    /// @return The assigned base core potential.
    auto operator=(CAtomCorePotential &&other) noexcept -> CAtomCorePotential &;

    /// @brief The equality operator.
    /// @param other The base core potential to be compared.
    /// @return True if base core potentials are equal, False otherwise.
    auto operator==(const CAtomCorePotential &other) const -> bool;

    /// @brief The equality operator.
    /// @param other The base core potential to be compared.
    /// @return True if base core potentials are not equal, False otherwise.
    auto operator!=(const CAtomCorePotential &other) const -> bool;

    /// @brief Sets exponents of primittive base core potentials to specific vector
    /// of exponents.
    /// @param exponents The vector of exponents.
    auto set_exponents(const std::vector<double> &exponents) -> void;

    /// @brief Sets expansion factors of primitive base core potentials to
    /// specific vector of expansion factors.
    /// @param factors The vector of expansion factors.
    auto set_factors(const std::vector<double> &factors) -> void;
    
    /// @brief Sets radial orders of primitive base core potentials to
    /// specific vector of radial orders.
    /// @param radial_orders The vector of radial orders.
    auto set_radial_orders(const std::vector<int> &radial_orders) -> void;

    /// @brief Adds primittive base core potential to base core potential.
    /// @param exponent The exponent of primitive base core potential.
    /// @param factor The expansion factor of primitive base core potential.
    /// @param radial_order  The radial order of primitive base core potential.
    auto add(const double exponent, const double factor, const int radial_order) -> void;

    /// @brief Gets vector of exponents of primitive base core potentials.
    /// @return The vector of exponents of primitive base core potentials.
    auto get_exponents() const -> std::vector<double>;

    /// @brief Gets vector of expansion factors of primitive base core potentials.
    /// @return The vector of expansion factors of primitive base core potentials.
    auto get_factors() const -> std::vector<double>;

    /// @brief Gets vector of radial orders of primitive base core potentials.
    /// @return The vector of radial orders of primitive base core potentials.
    auto get_radial_orders() const -> std::vector<int>;

    /// @brief Gets number of primitive  base core potentials in base core potential.
    /// @return The number of primitive  base core potentials in base core potential.
    auto number_of_primitive_potentials() const -> size_t;

   private:
    /// @brief The vector of exponents of primitive local potentials.
    std::vector<double> _exponents;

    /// @brief The vector of expansion factors of primitive local potentials.
    std::vector<double> _factors;

    /// @brief The vector of radial orders of primitive local potentials.
    std::vector<int> _radial_orders;
};

#endif /* AtomCorePotential_hpp */
