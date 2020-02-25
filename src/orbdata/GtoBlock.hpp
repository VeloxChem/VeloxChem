//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef GtoBlock_hpp
#define GtoBlock_hpp

#include <tuple>

#include "MemBlock.hpp"
#include "MemBlock2D.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

/**
 Class CGtoBlock stores data about basis functions and provides set of methods
 for manipulating with basis functions.

 @author Z. Rinkevicius
 */
class CGtoBlock
{
    /**
     The angular momentum.
     */
    int32_t _angularMomentum;

    /**
     The contraction scheme (contracted basis functions start/end position in
     primitive Gaussian functions vector, contracted basis functions indexes).
     */
    CMemBlock2D<int32_t> _contrPattern;

    /**
     The primitives Gaussian functions data (exponents, norm. factors and
     coordinates).
     */
    CMemBlock2D<double> _gtoPrimitives;

   public:
    /**
     Creates an empty GTOs block object.
     */
    CGtoBlock();

    /**
     Creates a GTOs block object.

     @param gtoPrimitives the primitives Gaussian functions data (exponents,
                          coordinates).
     @param contrPattern the contraction pattern (contracted basis functions
                         start/end position in primitive Gaussian functions
                         vector, basis functions indexes).
     @param angularMomentum the angular momentum of contracted basis functions.
     */
    CGtoBlock(const CMemBlock2D<double>& gtoPrimitives, const CMemBlock2D<int32_t>& contrPattern, const int32_t angularMomentum);

    /**
     Creates a GTOs block object.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param angularMomentum the angular momentum of contracted basis functions.
     */
    CGtoBlock(const CMolecule& molecule, const CMolecularBasis& basis, const int32_t angularMomentum);

    /**
     Creates a GTOs block object from list of atoms in molecule.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param iAtom the index of first atom in list of atoms.
     @param nAtoms the number of atoms in list of atoms.
     @param angularMomentum the angular momentum of contracted basis functions.
     */
    CGtoBlock(const CMolecule& molecule, const CMolecularBasis& basis, const int32_t iAtom, const int32_t nAtoms, const int32_t angularMomentum);

    /**
     Creates a GTOs block object by copying other GTOs block object.

     @param source the GTOs block object.
     */
    CGtoBlock(const CGtoBlock& source);

    /**
     Creates a GTOs block object by moving other GTOs block object.

     @param source the GTOs block object.
     */
    CGtoBlock(CGtoBlock&& source) noexcept;

    /**
     Destroys a GTOs block object.
     */
    ~CGtoBlock();

    /**
     Assigns a GTOs block object by copying other GTOs block object.

     @param source the GTOs block object.
     */
    CGtoBlock& operator=(const CGtoBlock& source);

    /**
     Assigns a GTOs block object by moving other GTOs block object.

     @param source the GTOs block object.
     */
    CGtoBlock& operator=(CGtoBlock&& source) noexcept;

    /**
     Compares GTOs block object with other GTOs block object.

     @param other the GTOs block object.
     @return true if GTOs block objects are equal, false otherwise.
     */
    bool operator==(const CGtoBlock& other) const;

    /**
     Compares GTOs block object with other GTOs block object.

     @param other the GTOs block object.
     @return true if GTOs block objects are not equal, false otherwise.
     */
    bool operator!=(const CGtoBlock& other) const;

    /**
     Set angular momentum of GTOs block.

     @param angularMomentum the angular momentum.
     */
    void setAngularMomentum(const int32_t angularMomentum);

    /**
     Compresses other GTOs block data into GTOs block object without changing
     dimensions of GTOs block object. Compression is performed using specified
     screening pattern.

     @param source the other GTOs block object.
     @param screeningFactors the vector of screening factors (absolute values).
     @param screeningThreshold the screening threshold.
     @return the number of primitive GTOs functions and number of contracted
             basis functions.
     */
    std::tuple<int32_t, int32_t> compress(const CGtoBlock& source, const CMemBlock<double>& screeningFactors, const double screeningThreshold);

    /**
     Gets angular momentum of GTOs block.

     @return the angular momentum.
     */
    int32_t getAngularMomentum() const;

    /**
     Checks if GTOs block is empty.

     @return true if GTOs block is empty, false otherwise.
     */
    bool empty() const;

    /**
     Gets number of primitive Gaussian functions in GTOs block.

     @return the number of primitive Gaussian functions.
     */
    int32_t getNumberOfPrimGtos() const;

    /**
     Gets number of contracted basis functions in GTOs block.

     @return the number of contracted basis functions.
     */
    int32_t getNumberOfContrGtos() const;

    /**
     Gets constant pointer to basis function start positions in primitive
     Gaussian functions vector.

     @return the start positions of basis fucntions.
     */
    const int32_t* getStartPositions() const;

    /**
     Gets pointer to basis function start positions in primitive
     Gaussian functions vector.

     @return the start positions of basis fucntions.
     */
    int32_t* getStartPositions();

    /**
     Gets constant pointer to basis function end positions in primitive
     Gaussian functions vector.

     @return the end positions of basis fucntions.
     */
    const int32_t* getEndPositions() const;

    /**
     Gets pointer to basis function end positions in primitive
     Gaussian functions vector.

     @return the end positions of basis fucntions.
     */
    int32_t* getEndPositions();

    /**
     Gets constant pointer to atomic identifiers vector.

     @return the atomic identifiers.
     */
    const int32_t* getAtomicIdentifiers() const;

    /**
     Gets pointer to atomic identifiers vector.

     @return the atomic identifiers.
     */
    int32_t* getAtomicIdentifiers();

    /**
     Gets constant pointer to basis function indexes in full AO basis for
     specific angular momentum component.

     @param iComponent the component of angular momentum.
     @return the indexes in full AO basis.
     */
    const int32_t* getIdentifiers(const int32_t iComponent) const;

    /**
     Gets pointer to basis function indexes in full AO basis for
     specific angular momentum component.

     @param iComponent the component of angular momentum.
     @return the indexes in full AO basis.
     */
    int32_t* getIdentifiers(const int32_t iComponent);

    /**
     Gets constant pointer to exponents of primitive Gaussian functions.

     @return the exponents of primitive Gaussian functions.
     */
    const double* getExponents() const;

    /**
     Gets pointer to exponents of primitive Gaussian functions.

     @return the exponents of primitive Gaussian functions.
     */
    double* getExponents();

    /**
     Gets constant pointer to normalization factors of primitive Gaussian
     functions.

     @return the normalization factors of primitive Gaussian functions.
     */
    const double* getNormFactors() const;

    /**
     Gets pointer to normalization factors of primitive Gaussian
     functions.

     @return the normalization factors of primitive Gaussian functions.
     */
    double* getNormFactors();

    /**
     Gets constant pointer to Cartesian X coordinates of primitive Gaussian
     functions.

     @return the exponents of primitive Gaussian functions.
     */
    const double* getCoordinatesX() const;

    /**
     Gets pointer to Cartesian X coordinates of primitive Gaussian functions.

     @return the exponents of primitive Gaussian functions.
     */
    double* getCoordinatesX();

    /**
     Gets constant pointer to Cartesian Y coordinates of primitive Gaussian
     functions.

     @return the exponents of primitive Gaussian functions.
     */
    const double* getCoordinatesY() const;

    /**
     Gets pointer to Cartesian Y coordinates of primitive Gaussian functions.

     @return the exponents of primitive Gaussian functions.
     */
    double* getCoordinatesY();

    /**
     Gets constant pointer to Cartesian Z coordinates of primitive Gaussian
     functions.

     @return the exponents of primitive Gaussian functions.
     */
    const double* getCoordinatesZ() const;

    /**
     Gets pointer to Cartesian Z coordinates of primitive Gaussian functions.

     @return the exponents of primitive Gaussian functions.
     */
    double* getCoordinatesZ();

    /**
     Gets maximum contraction depth of primitive Gaussian functions i.e. maximum
     number of primitive Gaussian functions in contracted basis fucntion.

     @return the contraction depth.
     */
    int32_t getMaxContractionDepth() const;

    /**
     Converts GTOs block object to text output and insert it into output
     text stream.

     @param output the output text stream.
     @param source the GTOs block object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CGtoBlock& source);
};

#endif /* GtoBlock_hpp */
