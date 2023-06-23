#ifndef AtomBasis_hpp
#define AtomBasis_hpp

#include <cstdint>
#include <set>
#include <string>
#include <vector>

#include "BasisFunction.hpp"

/**
 Class CAtomBasis stores data about atomic basis set and provides set of methods
 for handling of atomic basis set data.

 @author Z. Rinkevicius
 */
class CAtomBasis
{
    /**
     The vector of GTOs.
     */
    std::vector<CBasisFunction> _functions;

    /**
     The identifier of chemical element.
     */
    int64_t _identifier;

    /**
     The name of atomic basis.
     */
    std::string _name;

    /**
     The effective core potential label of atomic basis.
     */
    std::string _ecp_label;

   public:
    /**
     Creates an empty atom basis.
     */
    CAtomBasis() = default;

    /**
     Creates an atom basis.

     @param functions the vector of GTOs.
     @param identifier the identifier of chemical element.
     @param name the name of atom basis.
     */
    CAtomBasis(const std::vector<CBasisFunction>& functions, const int64_t identifier, const std::string& name, const std::string& ecp_label);

    /**
     Sets identifier of chemical element in atom basis.

     @param identifier the identifier of chemical element.
     */
    auto setIdentifier(const int64_t identifier) -> void;

    /**
     Sets name of  atom basis.

     @param name the name of atom basis.
     */
    auto setName(const std::string& name) -> void;

    /**
     Sets effective core potential label of  atom basis.

     @param label the effective core potential label of atom basis.
     */
    auto setEffectiveCorePotentialLabel(const std::string& label) -> void;

    /**
     Adds basis function  to atom basis.

     @param function the basis function.
     */
    auto add(const CBasisFunction& function) -> void;

    /**
     Reduces atom basis set to valence atom basis.

     @return the valence atom basis.
     */
    auto reduceToValenceBasis() const -> CAtomBasis;

    /**
     Gets vector of GTOs.

     @return the vector of GTOs in atom basis.
     */
    auto getBasisFunctions() const -> std::vector<CBasisFunction>;

    /**
      Gets vector of GTOs with specific angular momentum in atom basis.

     @param angmom the angular momentum.
     @return the vector of GTOs.
     */
    auto getBasisFunctions(const int64_t angmom) const -> std::vector<CBasisFunction>;

    /**
      Gets vector of GTOs with specific angular momentum and number of primitives in atom basis.

     @param angmom the angular momentum.
     @param npgtos the number of primitive GTOs.
     @return the vector of GTOs.
     */
    auto getBasisFunctions(const int64_t angmom, const int64_t npgtos) const -> std::vector<CBasisFunction>;

    /**
     Gets identifier of chemical element.

     @return the identifier of chemical element.
     */
    auto getIdentifier() const -> int64_t;

    /**
     Gets name of atom basis.

     @return the name of atom basis.
     */
    auto getName() const -> std::string;

    /**
     Gets effective core potential label of atom basis.

     @return the label of effective core potential.
     */
    auto getEffectiveCorePotentialLabel() const -> std::string;

    /**
     Checks if atomic basis requires effective core potential.

     @return the label of effective core potential.
     */
    auto needEffectiveCorePotential() const -> bool;

    /**
     Gets maximum angular momentum.

     @return the maximum angular momentum.
     */
    auto getMaxAngularMomentum() const -> int64_t;

    /**
     Gets number of GTOs with specific angular momentum.

     @param angmom the angular momentum.
     @return the number of GTOs.
     */
    auto getNumberOfBasisFunctions(const int64_t angmom) const -> int64_t;

    /**
     Gets number of GTOs with specific angular momentum and number of primitve GTOs.

     @param angmom the angular momentum.
     @param npgtos the number of primitive GTOs.
     @return the number of GTOs.
     */
    auto getNumberOfBasisFunctions(const int64_t angmom, const int64_t npgtos) const -> int64_t;

    /**
     Gets number of primitive Gaussain functions with requested angular momentum.

     @param angmom the angular momentum.
     @return the number of primitive Gaussian functions.
     */
    auto getNumberOfPrimitiveFunctions(const int64_t angmom) const -> int64_t;

    /**
     Get vector of unique contraction numbers of GTOs with given angular momentum in atom basis.

     @param angmom the angular momentum.
     @return the vector of unique contraction numbers.
    */
    auto getContractionDepths(const int64_t angmom) const -> std::set<int64_t>;

    /**
     Gets contraction string for atom basis.

     @return the contraction string.
     */
    auto getContractionString() const -> std::string;

    /**
     Gets primitives strig for atom basis.

     @return the primitives string.
     */
    auto getPrimitivesString() const -> std::string;
};

#endif /* AtomBasis_hpp */
