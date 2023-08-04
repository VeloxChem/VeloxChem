#ifndef GtoBlock_hpp
#define GtoBlock_hpp

#include <cstdint>
#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

/**
 Class CGtoBlock stores data about contarcted GTOs block and provides set of methods
 for manipulating with contarcted  GTOs block.

 @author Z. Rinkevicius
 */
class CGtoBlock
{
    /**
     The vector of Cartesian coordinates of contracted GTOs.
     */
    std::vector<TPoint3D> _coordinates;

    /**
     The vector of exponents of primitive GTOs.
     */
    std::vector<double> _exponents;

    /**
     The vector of normalization factors of primitive GTOs.
     */
    std::vector<double> _norms;

    /**
     The vector of AO indexes of contracted GTOs.
     */
    std::vector<int64_t> _orb_indexes;

    /**
     The vector of atomic indexes of contracted GTOs.
     */
    std::vector<int64_t> _atm_indexes;

    /**
     The angular momentum of basis function.
     */
    int64_t _angmom;

    /**
     The number of primitive GTOs in contracted GTO.
     */
    int64_t _npgtos;

   public:
    /**
     Creates an empty contarcted GTOs block.
     */
    CGtoBlock() = default;

    /**
     Creates a contarcted GTOs block.

     @param coordinates the vector of basis function coordinates.
     @param exponents the vector of exponents of primitive GTOs.
     @param norms the vector of normalization factors of primitive GTOs.
     @param orb_indexes the vector of  AO indexes.
     @param atm_indexes the vector of  atomic indexes.
     @param angmom the angular momentum of GTOs.
     @param npgtos the number of primitive GTOs in contracted GTO.
     */
    CGtoBlock(const std::vector<TPoint3D>& coordinates,
              const std::vector<double>&   exponents,
              const std::vector<double>&   norms,
              const std::vector<int64_t>&  orb_indexes,
              const std::vector<int64_t>&  atm_indexes,
              const int64_t                angmom,
              const int64_t                npgtos);

    /**
     Creates a contarcted GTOs block.

     @param basis the molecular basis.
     @param molecule the molecule.
     @param angmom the angular momentum of basis functions.
     @param npgtos the number of primitive GTOs in contracted GTO.
     */
    CGtoBlock(const CMolecularBasis& basis, const CMolecule& molecule, const int64_t angmom, const int64_t npgtos);

    /**
     Creates  contarcted GTOs block.

     @param basis the molecular basis.
     @param molecule the molecule.
     @param atoms the vector of atoms to select basis functions.
     @param angmom the angular momentum of basis functions.
     @param npgtos the number of primitive GTOs in contracted GTO.
     */
    CGtoBlock(const CMolecularBasis& basis, const CMolecule& molecule, const std::vector<int64_t>& atoms, const int64_t angmom, const int64_t npgtos);

    /**
     Gets vector of GTOs coordinates.

     @return the vector of GTOs coordinates.
     */
    auto getCoordinates() const -> std::vector<TPoint3D>;

    /**
     Gets vector of GTOs exponents.

     @return the vector of GTOs exponents.
     */
    auto getExponents() const -> std::vector<double>;

    /**
     Gets vector of GTOs normalization factors.

     @return the vector of GTOs normalization factors.
     */
    auto getNormalizationFactors() const -> std::vector<double>;

    /**
     Gets vector of orbital indexes of contracted GTOs.

     @return the vector of orbital indexes of GTOs.
     */
    auto getOrbitalIndexes() const -> std::vector<int64_t>;

    /**
     Gets vector of atomic indexes of contracted GTOs.

     @return the vector of atomic indexes of GTOs.
     */
    auto getAtomicIndexes() const -> std::vector<int64_t>;

    /**
     Gets angular momentum.

     @return the angular momentum.
     */
    auto getAngularMomentum() const -> int64_t;

    /**
     Gets number of primitive GTOs in contracted GTO.

     @return the number of primtive GTOs in contracted GTO.
     */
    auto getNumberOfPrimitives() const -> int64_t;

    /**
     Gets number of GTOs in contracted GTOs block.

     @return the number of GTOs in contracted GTOs block.
     */
    auto getNumberOfBasisFunctions() const -> int64_t;
};

#endif /* GtoBlock_hpp */
