#ifndef MolecularBasis_hpp
#define MolecularBasis_hpp

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "AtomBasis.hpp"

/**
Class CMolecularBasis stores data about molecular basis set and provides set of
methods for handling of molecular basis set data.

@author Z. Rinkevicius
*/
class CMolecularBasis
{
    /**
     The vector of atom basis sets.
     */
    std::vector<CAtomBasis> _basis_sets;

    /**
     The vector of atom basis sets indexes.
     */
    std::vector<int64_t> _indexes;

    /**
     Gets  basis set labels frequency map.

     @return the basis set labels frequency map of molecular basis.
     */
    auto _getLabelsFrequencyMap() const -> std::unordered_map<std::string, int64_t>;

   public:
    /**
     Creates an empty molecular basis.
     */
    CMolecularBasis() = default;

    /**
     Creates a molecular basis.

     @param basis_sets the vector of unique atom basis sets.
     @param indexes the vector of atom basis sets indexes.

     */
    CMolecularBasis(const std::vector<CAtomBasis>& basis_sets, const std::vector<int64_t>& indexes);

    /**
     Adds atom basis object to molecular basis.

     @param basis the atom basis object.
     */
    auto add(const CAtomBasis& basis) -> void;

    /**
     Reduces molecular basis to valence molecular basis.

     @return the valence molecular basis.
     */
    auto reduceToValenceBasis() const -> CMolecularBasis;
    
    /**
     Slice fraction of molecular basis for specific atoms.

     @param atoms the vector of atoms. 
     @return the fractional molecular basis.
     */
    auto slice(const std::vector<int64_t>& atoms) const -> CMolecularBasis;

    /**
     Gets vector of unique atomic basis sets.

     @return the vector of atomic basis sets.
     */
    auto getBasisSets() const -> std::vector<CAtomBasis>;

    /**
     Gets vector of indexes for atomic basis sets.

     @return the vector of indexes atomic basis sets.
     */
    auto getBasisSetsIndexes() const -> std::vector<int64_t>;

    /**
     Gets maximum angular momentum of molecular basis.

     @return the maximum angular momentum.
     */
    auto getMaxAngularMomentum() const -> int64_t;

    /**
     Gets maximum angular momentum of molecular basis.

     @param atoms the vector of atoms to select contracted GTOs.
     @return the maximum angular momentum.
     */
    auto getMaxAngularMomentum(const std::vector<int64_t>& atoms) const -> int64_t;

    /**
     Gets vector of contracted GTOs.

     @return the vector of contracted GTOs in molecular basis.
     */
    auto getBasisFunctions() const -> std::vector<CBasisFunction>;

    /**
      Gets vector of contracted GTOs with specific angular momentum in  molecular basis.

     @param angmom the angular momentum.
     @return the vector of contracted GTOs.
     */
    auto getBasisFunctions(const int64_t angmom) const -> std::vector<CBasisFunction>;

    /**
      Gets vector of contracted GTOs with specific angular momentum and number of primitives in molecular basis.

     @param angmom the angular momentum.
     @param npgtos the number of primitive GTOs.
     @return the vector of contracted GTOs.
     */
    auto getBasisFunctions(const int64_t angmom, const int64_t npgtos) const -> std::vector<CBasisFunction>;

    /**
     Gets vector of contracted GTOs in  molecular basis for vector of specific atoms.

     @param atoms  the vector of atoms to select contracted GTOs.
     @return the vector of contracted GTOs in molecular basis.
     */
    auto getBasisFunctions(const std::vector<int64_t>& atoms) const -> std::vector<CBasisFunction>;

    /**
      Gets vector of contracted GTOs with specific angular momentum in molecular basis for specific atoms.

     @param atoms the vector of atoms to select contracted GTOs.
     @param angmom the angular momentum.
     @return the vector of contracted GTOs.
     */
    auto getBasisFunctions(const std::vector<int64_t>& atoms, const int64_t angmom) const -> std::vector<CBasisFunction>;

    /**
      Gets vector of contracted GTOs with specific angular momentum and number of primitives in molecular basis for specific atoms.

     @param atoms the vector of atoms to select contracted GTOs.
     @param angmom the angular momentum.
     @param npgtos the number of primitive GTOs.
     @return the vector of contracted GTOs.
     */
    auto getBasisFunctions(const std::vector<int64_t>& atoms, const int64_t angmom, const int64_t npgtos) const -> std::vector<CBasisFunction>;

    /**
     Gets vector of atomic indexes.

     @return the vector of atomic indexes mapping contracted GTOs in molecular basis.
     */
    auto getAtomicIndexes() const -> std::vector<int64_t>;

    /**
     Gets vector of atomic indexes for contracted GTOs of specific angular momentum in molecular basis.

     @param angmom the angular momentum.
     @return the vector of atomic indexes mapping contracted GTOs in molecular basis.
     */
    auto getAtomicIndexes(const int64_t angmom) const -> std::vector<int64_t>;

    /**
     Gets vector of atomic indexes for contracted GTOs of specific angular momentum and of number of primitives in molecular basis.

     @param angmom the angular momentum.
     @param npgtos the number of primitive GTOs.
     @return the vector of atomic indexes mapping contracted GTOs in molecular basis.
     */
    auto getAtomicIndexes(const int64_t angmom, const int64_t npgtos) const -> std::vector<int64_t>;

    /**
     Gets vector of atomic indexes.

     @param atoms the vector of atoms to select contracted GTOs.
     @return the vector of atomic indexes mapping contracted GTOs in molecular basis.
     */
    auto getAtomicIndexes(const std::vector<int64_t>& atoms) const -> std::vector<int64_t>;

    /**
     Gets vector of atomic indexes for contracted GTOs of specific angular momentum in molecular basis.

     @param atoms the vector of atoms to select contracted GTOs.
     @param angmom the angular momentum.
     @return the vector of atomic indexes mapping contracted GTOs in molecular basis.
     */
    auto getAtomicIndexes(const std::vector<int64_t>& atoms, const int64_t angmom) const -> std::vector<int64_t>;

    /**
     Gets vector of atomic indexes for contracted GTOs of specific angular momentum and of number of primitives in molecular basis.

     @param atoms the vector of atoms to select contracted GTOs.
     @param angmom the angular momentum.
     @param npgtos the number of primitive GTOs.
     @return the vector of atomic indexes mapping contracted GTOs in molecular basis.
     */
    auto getAtomicIndexes(const std::vector<int64_t>& atoms, const int64_t angmom, const int64_t npgtos) const -> std::vector<int64_t>;

    /**
     Determines number of contracted GTOs with specific angular momentum in molecular basis.

     @param angmom the angular momentum.
     @return the number of contracted GTOs.
     */
    auto getNumberOfBasisFunctions(const int64_t angmom) const -> int64_t;

    /**
     Determines number of contracted GTOs with specific angular momentum in molecular basis.

     @param angmom the angular momentum.
     @param npgtos the number of primitive GTOs.
     @return the number of contracted GTOs.
     */
    auto getNumberOfBasisFunctions(const int64_t angmom, const int64_t npgtos) const -> int64_t;

    /**
     Determines number of contracted GTOs with specific angular momentum in molecular basis.

     @param atoms the vector of atoms to select contracted GTOs.
     @param angmom the angular momentum.
     @return the number of contracted GTOs.
     */
    auto getNumberOfBasisFunctions(const std::vector<int64_t>& atoms, const int64_t angmom) const -> int64_t;

    /**
     Determines number of contracted GTOs with specific angular momentum in molecular basis.

     @param atoms the vector of atoms to select contracted GTOs.
     @param angmom the angular momentum.
     @param npgtos the number of primitive GTOs.
     @return the number of contracted GTOs.
     */
    auto getNumberOfBasisFunctions(const std::vector<int64_t>& atoms, const int64_t angmom, const int64_t npgtos) const -> int64_t;

    /**
     Determines number of primitive Gaussian functions with specific angular momentum in molecular basis.

     @param angmom the angular momentum.
     @return the number of primitive Gaussian functions.
     */
    auto getNumberOfPrimitiveFunctions(const int64_t angmom) const -> int64_t;

    /**
     Determines number of primitive Gaussian functions with specific angular momentum in molecular basis.

     @param atoms the vector of atoms to select primitive Gaussian functions.
     @param angmom the angular momentum.
     @return the number of primitive Gaussian functions.
     */
    auto getNumberOfPrimitiveFunctions(const std::vector<int64_t>& atoms, const int64_t angmom) const -> int64_t;

    /**
     Get vector of unique contraction numbers of contracted GTOs with given angular momentum in molecular basis.

     @param angmom the angular momentum.
     @return the vector of unique contraction numbers.
    */
    auto getContractionDepths(const int64_t angmom) const -> std::set<int64_t>;

    /**
     Get vector of unique contraction numbers of contracted GTOs with given angular momentum in molecular basis.

     @param atoms the vector of atoms to select contracted GTOs.
     @param angmom the angular momentum.
     @return the vector of unique contraction numbers.
    */
    auto getContractionDepths(const std::vector<int64_t>& atoms, const int64_t angmom) const -> std::set<int64_t>;

    /**
     Determines size of contracted atomic orbitals basis.

     @return the size of contracted atomic orbitals basis.
     */
    auto getDimensionsOfBasis() const -> int64_t;

    /**
     Determines partial size up to specific angular momentum of contracted atomic orbitals basis.

     @param angmom the angular momentum.
     @return the partial size of contracted atomic orbitals basis.
     */
    auto getDimensionsOfBasis(const int64_t angmom) const -> int64_t;

    /**
     Determines size of primitive atomic orbitals basis.

     @return the size of primitive atomic orbitals basis.
     */
    auto getDimensionsOfPrimitiveBasis() const -> int64_t;

    /**
     Creates map indexes for contracted contracted GTOs with specific angular momentum and number
     of primitive Gaussian functions.

     @param angmom the angular momentum.
     @param npgtos the number of primitive GTOs.
     @return the vector of indexes (n_t, index_0,...,index_k).
     */
    auto getIndexMap(const int64_t angmom, const int64_t npgtos) const -> std::vector<int64_t>;

    /**
     Creates map indexes for contracted contracted GTOs with specific angular momentum and number
     of primitive Gaussian functions.

     @param atoms the vector of atoms to select contracted GTOs.
     @param angmom the angular momentum.
     @param npgtos the number of primitive GTOs.
     @return the vector of indexes (n_t, index_0,...,index_k ).
     */
    auto getIndexMap(const std::vector<int64_t>& atoms, const int64_t angmom, const int64_t npgtos) const -> std::vector<int64_t>;

    /**
     Prints atomic orbitals basis information to output stream for selected molecule.

     @param title the header line of atomic orbitals basis output.
     @return the basis set information string.
     */
    auto printBasis(const std::string& title) const -> std::string;

    /**
     Prints atomic orbitals basis information to output stream for selected molecule.

     @return the basis set information string.
     */
    auto printBasis() const -> std::string;
};

#endif /* MolecularBasis_hpp */
