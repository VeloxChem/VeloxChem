//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#ifndef GtoContainer_hpp
#define GtoContainer_hpp

#include <vector>

#include "GtoBlock.hpp"
#include "MemBlock.hpp"
#include "SphericalMomentum.hpp"
#include "VecMemBlocks.hpp"

/**
 Class CGtoContainer stores vector of GTOs block objects and provides set of
 methods for manipulating with basis functions of various angular momentum.

 @author Z. Rinkevicius
 */
class CGtoContainer
{
    /**
     The maximum angular momentum.
     */
    int32_t _maxAngularMomentum;

    /**
     The vector of GTOs block objects.
     */
    std::vector<CGtoBlock> _gtoBlocks;

    /**
     Gets number of buffer components for specific angular momentum.

     @param angularMomentum the angular momentum.
     @return the number of buffer components.
     */
    int32_t _getPrimAngComponents(const int32_t angularMomentum) const;

   public:
    /**
     Creates an empty GTOs container object.
     */
    CGtoContainer();

    /**
     Creates a GTOs container object.

     @param gtoBlocks the vector GTOs block objects.
     */
    CGtoContainer(const std::vector<CGtoBlock>& gtoBlocks);

    /**
     Creates a GTOs container object.

     @param molecule the molecule.
     @param basis the molecular basis.
     */
    CGtoContainer(const CMolecule& molecule, const CMolecularBasis& basis);

    /**
     Creates a GTOs container object from list of atoms in molecule.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param iAtom the index of first atom in list of atoms.
     @param nAtoms the number of atoms in list of atoms.
     */
    CGtoContainer(const CMolecule& molecule, const CMolecularBasis& basis, const int32_t iAtom, const int32_t nAtoms);

    /**
     Creates a GTOs container object with atom and angular momentum specific
     GTO blocks.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param flag the flag.
     */
    CGtoContainer(const CMolecule&       molecule,
                  const CMolecularBasis& basis,
                  const std::string&     flag);

    /**
     Creates a GTOs container object according to given pattern of atomic batches
     in molecule.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param batches the distribution pattern for atomic batches.
     */
    CGtoContainer(const CMolecule& molecule, const CMolecularBasis& basis, const CMemBlock2D<int32_t>& batches);

    /**
     Creates a GTOs container object by copying other GTOs container object.

     @param source the GTOs container object.
     */
    CGtoContainer(const CGtoContainer& source);

    /**
     Creates a GTOs container object by moving other GTOs container object.

     @param source the GTOs container object.
     */
    CGtoContainer(CGtoContainer&& source) noexcept;

    /**
     Destroys a GTOs container object.
     */
    ~CGtoContainer();

    /**
     Assigns a GTOs container object by copying other GTOs container object.

     @param source the GTOs container object.
     */
    CGtoContainer& operator=(const CGtoContainer& source);

    /**
     Assigns a GTOs container object by moving other GTOs container object.

     @param source the GTOs container object.
     */
    CGtoContainer& operator=(CGtoContainer&& source) noexcept;

    /**
     Compares GTOs container object with other GTOs container object.

     @param other the GTOs container object.
     @return true if GTOs container objects are equal, false otherwise.
     */
    bool operator==(const CGtoContainer& other) const;

    /**
     Compares GTOs container object with other GTOs container object.

     @param other the GTOs container object.
     @return true if GTOs container objects are not equal, false otherwise.
     */
    bool operator!=(const CGtoContainer& other) const;

    /**
     Compresses the other GTOs container object into GTOs container object by
     applying specific screening pattern.

     @param source the other GTOs container.
     @param reducedDimensions the reduced dimensions of GTOs container.
     @param screeningFactors the vector of screening factors (absolute values).
     @param screeningThreshold the screening threshold.
     */
    void compress(const CGtoContainer&        source,
                  CMemBlock2D<int32_t>&       reducedDimensions,
                  const CVecMemBlock<double>& screeningFactors,
                  const double                screeningThreshold);

    /**
     Gets maximum angular momentum of GTOs block objects in GTOs container.

     @return the maxumum angular momentum.
     */
    int32_t getMaxAngularMomentum() const;

    /**
     Gets angular momentum of specific GTOs block object from GTOs container.

     @param iBlock the index of GTOs block object in GTOs container.
     @return the angular momentum.
     */
    int32_t getAngularMomentum(const int32_t iBlock) const;

    /**
     Gets number of GTOs blocks in GTOs container.

     @return the number of GTOs blocks.
     */
    int32_t getNumberOfGtoBlocks() const;

    /**
     Gets maximum number of primitive Gaussian functions within GTOs container.

     @return the number of primitive Gaussian functions.
     */
    int32_t getMaxNumberOfPrimGtos() const;

    /**
     Gets maximum number of contracted basis functions within GTOs container.

     @return the number of contracted GTOs.
     */
    int32_t getMaxNumberOfContrGtos() const;

    /**
     Gets number of primitive Gaussian functions in specific GTOs block object
     from GTOs container.

     @param iBlock the index of GTOs block object in GTOs container.
     @return the number of primitive Gaussian functions.
     */
    int32_t getNumberOfPrimGtos(const int32_t iBlock) const;

    /**
     Gets number of contracted basis functions in specific GTOs block object
     from GTOs container.

     @param iBlock the index of GTOs block object in GTOs container.
     @return the number of contracted basis functions.
     */
    int32_t getNumberOfContrGtos(const int32_t iBlock) const;

    /**
     Gets total number of atomic orbitals.

     @return the number of atomic orbitals.
     */
    int32_t getNumberOfAtomicOrbitals() const;

    /**
     Gets constant pointer to basis function start positions in primitive
     Gaussian functions vector from specific GTOs block object in GTOs
     container.

     @param iBlock the index of GTOs block object in GTOs container.
     @return the start positions of basis fucntions.
     */
    const int32_t* getStartPositions(const int32_t iBlock) const;

    /**
     Gets constant pointer to basis function end positions in primitive
     Gaussian functions vector from specific GTOs block object in GTOs
     container.

     @param iBlock the index of GTOs block object in GTOs container.
     @return the end positions of basis fucntions.
     */
    const int32_t* getEndPositions(const int32_t iBlock) const;

    /**
     Gets constant pointer to basis function indexes in full AO basis for
     specific angular momentum component from specific GTOs block object in
     GTOs container.

     @param iBlock the index of GTOs block object in GTOs container.
     @param iComponent the component of angular momentum.
     @return the indexes in full AO basis.
     */
    const int32_t* getIdentifiers(const int32_t iBlock, const int32_t iComponent) const;

    /**
     Gets constant pointer to exponents of primitive Gaussian functions from
     specific GTOs block object in GTOs container.

     @param iBlock the index of GTOs block object in GTOs container.
     @return the exponents of primitive Gaussian functions.
     */
    const double* getExponents(const int32_t iBlock) const;

    /**
     Gets constant pointer to normalization factors of primitive Gaussian
     functions from specific GTOs block object in GTOs container.

     @param iBlock the index of GTOs block object in GTOs container.
     @return the normalization factors of primitive Gaussian functions.
     */
    const double* getNormFactors(const int32_t iBlock) const;

    /**
     Gets constant pointer to Cartesian X coordinates of primitive Gaussian
     functions from specific GTOs block object in GTOs container.

     @param iBlock the index of GTOs block object in GTOs container.
     @return the exponents of primitive Gaussian functions.
     */
    const double* getCoordinatesX(const int32_t iBlock) const;

    /**
     Gets constant pointer to Cartesian Y coordinates of primitive Gaussian
     functions from specific GTOs block object in GTOs container.

     @param iBlock the index of GTOs block object in GTOs container.
     @return the exponents of primitive Gaussian functions.
     */
    const double* getCoordinatesY(const int32_t iBlock) const;

    /**
     Gets constant pointer to Cartesian Z coordinates of primitive Gaussian
     functions from specific GTOs block object in GTOs container.

     @param iBlock the index of GTOs block object in GTOs container.
     @return the exponents of primitive Gaussian functions.
     */
    const double* getCoordinatesZ(const int32_t iBlock) const;

    /**
     Creates vector of memory block objects according to dimensions of primitive
     Gaussian functions space of each GTOs block object in GTOs container.

     @return the vector of memory block objects.
     */
    CVecMemBlock<double> getPrimBuffer() const;

    /**
     Creates vector of 2D memory block objects according to dimensions of
     primitive Gaussian functions space and angular momentum of each GTOs block
     object in GTOs container.

     @param nComponents the number of buffer vectors for each angular momentum
             component.
     @return the vector of 2D memory block objects.
     */
    CVecMemBlock2D<double> getPrimAngBuffer(const int32_t nComponents) const;

    /**
     Creates vector of 2D memory block objects according to dimensions of
     contracted Gaussian functions space and Cartesian angular momentum of
     each GTOs block object in GTOs container.

     @param nComponents the number of buffer vectors for each angular momentum
            component.
     @return the vector of 2D memory block objects.
     */
    CVecMemBlock2D<double> getCartesianBuffer(const int32_t nComponents) const;

    /**
     Creates vector of 2D memory block objects according to dimensions of
     contracted Gaussian functions space and spherical angular momentum of
     each GTOs block object in GTOs container.

     @param nComponents the number of buffer vectors for each angular momentum
            component.
     @return the vector of 2D memory block objects.
     */
    CVecMemBlock2D<double> getSphericalBuffer(const int32_t nComponents) const;

    /**
     Creates vector of spherical momentum objects corresponding to angular
     momentum of each GTOs block object.

     @return the vector of spherical momentum objects.
     */
    std::vector<CSphericalMomentum> getSphericalMomentumVector() const;

    /**
     Gets specific GTOs block object from GTOs container.

     @param iBlock the index of GTOs block object in GTOs container.
     @return the GTOs block object.
     */
    CGtoBlock getGtoBlock(const int32_t iBlock) const;

    /**
     Converts GTOs container object to text output and insert it into output
     text stream.

     @param output the output text stream.
     @param source the GTOs container object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CGtoContainer& source);
};

#endif /* GtoContainer_hpp */
