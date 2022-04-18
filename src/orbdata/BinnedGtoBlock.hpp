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

#ifndef BinnedGtoBlock_hpp
#define BinnedGtoBlock_hpp

#include "Buffer.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "AngularMomentum.hpp"

/**
 Class CBinnedGtoBlock stores data about basis functions and provides set of methods
 for manipulating with basis functions.

 @author Z. Rinkevicius
 */
template <typename T, typename B = mem::Host>
class CBinnedGtoBlock
{
    /**
     The angular momentum.
     */
    int32_t _angularMomentum;
    
    /**
     The number of primitive GTOs in single basis function.
     */
    int32_t _nPrimitiveGtos;
    
    /**
     The contracted basis functions indexes.
     */
    BufferXY<int32_t, B> _gtoIndexes;
    
    /**
     The atomic indexes of contracted basis functions.
     */
    BufferX<int32_t, B> _atomicIndexes;

    /**
     The primitives Gaussian functions data (exponents, normalization
     factors and coordinates).
     */
    BufferMY<T, B, 5> _gtoPrimitives;
    
   public:
    /**
     Creates an empty binned GTOs block object.
     */
    CBinnedGtoBlock()
        : _angularMomentum(0)

        , _nPrimitiveGtos(0)

        , _gtoIndexes(BufferXY<int32_t, B>())

        , _atomicIndexes(BufferX<int32_t, B>())

        , _gtoPrimitives(BufferMY<T, B, 5>())

    {
        
    }
    
    /**
    Creates a binned GTOs block object.

    @param angularMomentum the angular momentum of contracted basis functions.
    @param nPrimitiveGtos the number of primitive GTOs in contracted basis function.
    @param gtoIndexes the indexes of contracted basis functions.
    @param atomicIndexes the atomic indexes of contracted basis functions.
    @param gtoPrimitives the primitives Gaussian functions data (exponents,
           coordinates).
    */
    CBinnedGtoBlock(const int32_t               angularMomentum,
                    const int32_t               nPrimitiveGtos,
                    const BufferXY<int32_t, B>& gtoIndexes,
                    const BufferX<int32_t, B>&  atomicIndexes,
                    const BufferMY<T, B, 5>&    gtoPrimitives)
    
        : _angularMomentum(angularMomentum)

        , _nPrimitiveGtos(nPrimitiveGtos)

        , _gtoIndexes(gtoIndexes)

        , _atomicIndexes(atomicIndexes)

        , _gtoPrimitives(gtoPrimitives)

    {
        
    }
    
    /**
     Creates a binned GTOs block object for molecule.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param angularMomentum the angular momentum of contracted basis functions.
     @param nPrimitiveGtos the number of primitive GTOs in contracted basis function.
     */
    CBinnedGtoBlock(const CMolecule&       molecule,
                    const CMolecularBasis& basis,
                    const int32_t          angularMomentum,
                    const int32_t          nPrimitiveGtos)
    
        : CBinnedGtoBlock<T, B>(molecule, basis, 0, molecule.getNumberOfAtoms(), angularMomentum, nPrimitiveGtos)
    {
        
    }
    
    /**
     Creates a binned GTOs block object from list of atoms in molecule.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param iAtom the index of first atom in list of atoms.
     @param nAtoms the number of atoms in list of atoms.
     @param angularMomentum the angular momentum of contracted basis functions.
     @param nPrimitiveGtos the number of primitive GTOs in contracted basis function.
     */
    CBinnedGtoBlock(const CMolecule&       molecule,
                    const CMolecularBasis& basis,
                    const int32_t          iAtom,
                    const int32_t          nAtoms,
                    const int32_t          angularMomentum,
                    const int32_t          nPrimitiveGtos)
    
        : CBinnedGtoBlock<T, B>()
    {
        if (const auto nfuncs = basis.getNumberOfBasisFunctions(molecule, iAtom, nAtoms, angularMomentum, nPrimitiveGtos); nfuncs > 0)
        {
            // intialize GTOs dimensions
                
            _angularMomentum = angularMomentum;
                
            _nPrimitiveGtos = nPrimitiveGtos;
                
            const auto angcomps = angmom::to_SphericalComponents(_angularMomentum);
                
            // intialize GTOs data
                
            _gtoIndexes = BufferXY<int32_t, B>(angcomps, nfuncs);
                
            _atomicIndexes = BufferX<int32_t, B>(nfuncs);
                
            _gtoPrimitives = BufferMY<T, B, 5>(nfuncs * _nPrimitiveGtos);
                
            // set up pointers to molecular data

            auto molrx = molecule.getCoordinatesX();

            auto molry = molecule.getCoordinatesY();

            auto molrz = molecule.getCoordinatesZ();

            auto idselem = molecule.getIdsElemental();
            
            // determine partial dimensions of AO basis
                
            const auto npartdim = basis.getPartialDimensionsOfBasis(molecule, _angularMomentum);
                
            // determine offset in contracted GTOs block
                
            const auto ncoff = basis.getNumberOfBasisFunctions(molecule, 0, iAtom, _angularMomentum);
                
            // determine total numner of contracted GTOs
                
            auto ntfuncs = basis.getNumberOfBasisFunctions(molecule, _angularMomentum);
                
            // loop over atoms in molecule

            int32_t icgto = 0;
                
            int32_t ncgto = 0;

            for (int32_t i = iAtom; i < (iAtom + nAtoms); i++)
            {
                // get atom coordinates
                    
                const auto rx = static_cast<T>(molrx[i]);
                   
                const auto ry = static_cast<T>(molry[i]);
                    
                const auto rz = static_cast<T>(molrz[i]);
                    
                // loop over basis functions of i-th atom
                    
                const auto gtos = basis.getBasisFunctions(idselem[i], _angularMomentum);
                    
                for (size_t j = 0; j < gtos.size(); j++)
                {
                    const auto nprim = gtos[j].getNumberOfPrimitiveFunctions();
                        
                    if (nprim == _nPrimitiveGtos)
                    {
                        // set up angular indexes of contracted GTOs
                            
                        for (int32_t k = 0; k < angcomps; k++)
                        {
                            _gtoIndexes(k, icgto) = npartdim + k * ntfuncs + ncoff + ncgto;
                        }
                            
                        // set up atomic indexes of contracted GTOs
                            
                        _atomicIndexes(icgto) = i;
                            
                        // retrieve primitve exponents, norm. factors
                            
                        auto pexp = gtos[j].getExponents();
                            
                        auto pnorm = gtos[j].getNormalizationFactors();
                            
                        // set up primitives data
                
                        const auto koff = icgto * _nPrimitiveGtos;
                            
                        for (int32_t k = 0; k < nprim; k++)
                        {
                            // assign exponent, norm. factor
                            
                            _gtoPrimitives(0, koff + k) = static_cast<T>(pexp[k]);
                            
                            _gtoPrimitives(1, koff + k)= static_cast<T>(pnorm[k]);
                            
                            // assign atom coordinates
                            
                            _gtoPrimitives(2, koff + k) = rx;
                            
                            _gtoPrimitives(3, koff + k) = ry;
                            
                            _gtoPrimitives(4, koff + k) = rz;
                        }
                            
                        icgto++;
                    }
                        
                    ncgto++;
                }
            }
        }
    }

    /**
     Creates a binned GTOs block object by copying other binned GTOs block object.

     @param source the binned GTOs block object.
     */
    CBinnedGtoBlock(const CBinnedGtoBlock<T, B>& source)
    
        : _angularMomentum(source._angularMomentum)

        , _nPrimitiveGtos(source._nPrimitiveGtos)

        , _gtoIndexes(source._gtoIndexes)

        , _atomicIndexes(source._atomicIndexes)

        , _gtoPrimitives(source._gtoPrimitives)
    {
        
    }

    /**
     Creates a binned GTOs block object by moving other binned GTOs block object.

     @param source the binned GTOs block object.
     */
    CBinnedGtoBlock(CBinnedGtoBlock<T, B>&& source) noexcept
    
        : _angularMomentum(std::move(source._angularMomentum))

        , _nPrimitiveGtos(std::move(source._nPrimitiveGtos))

        , _gtoIndexes(std::move(source._gtoIndexes))

        , _atomicIndexes(std::move(source._atomicIndexes))

        , _gtoPrimitives(std::move(source._gtoPrimitives))
    {
        
    }

    /**
     Destroys a binned GTOs block object.
     */
    ~CBinnedGtoBlock()
    {
        
    }

    /**
     Assigns a binned GTOs block object by copying other binned GTOs block object.

     @param source the binned GTOs block object.
     */
    auto
    operator=(const CBinnedGtoBlock<T, B>& source) -> CBinnedGtoBlock<T, B>&
    {
        if (this == &source) return *this;

        _angularMomentum = source._angularMomentum;

        _nPrimitiveGtos = source._nPrimitiveGtos;

        _gtoIndexes = source._gtoIndexes;
            
        _atomicIndexes = source._atomicIndexes;

        _gtoPrimitives = source._gtoPrimitives;

        return *this;
    }

    /**
     Assigns a binned GTOs block object by moving other binned GTOs block object.

     @param source the binned GTOs block object.
     */
    auto
    operator=(CBinnedGtoBlock<T, B>&& source) noexcept -> CBinnedGtoBlock<T, B>&
    {
        if (this == &source) return *this;

        _angularMomentum = std::move(source._angularMomentum);

        _nPrimitiveGtos = std::move(source._nPrimitiveGtos);

        _gtoIndexes = std::move(source._gtoIndexes);
            
        _atomicIndexes = std::move(source._atomicIndexes);

        _gtoPrimitives = std::move(source._gtoPrimitives);

        return *this;
    }

    /**
     Compares binned GTOs block object with other binned GTOs block object.

     @param other the binned GTOs block object.
     @return true if binned GTOs block objects are equal, false otherwise.
     */
    auto operator==(const CBinnedGtoBlock<T, B>& other) const -> bool
    {
        if (_angularMomentum != other._angularMomentum) return false;

        if (_nPrimitiveGtos != other._nPrimitiveGtos) return false;
            
        if (_gtoIndexes != other._gtoIndexes) return false;
            
        if (_atomicIndexes != other._atomicIndexes) return false;
            
        if (_gtoPrimitives != other._gtoPrimitives) return false;

        return true;
    }

    /**
     Compares binned GTOs block object with other binned GTOs block object.

     @param other the binned GTOs block object.
     @return true if binned GTOs block objects are not equal, false otherwise.
     */
    auto
    operator!=(const CBinnedGtoBlock<T, B>& other) const -> bool
    {
        return !(*this == other);
    }
    
    /**
     Gets angular momentum of binned GTOs block.

     @return the angular momentum.
     */
    auto
    getAngularMomentum() const -> int32_t
    {
        return _angularMomentum;
    }

    /**
     Checks if binned GTOs block is empty.

     @return true if binned GTOs block is empty, false otherwise.
     */
    auto empty() const -> bool
    {
        return _atomicIndexes.empty();
    }

    /**
     Gets number of primitive Gaussian functions in contracted basis function.

     @return the number of primitive Gaussian functions.
     */
    auto
    getNumberOfPrimGtos() const -> int32_t
    {
        return _nPrimitiveGtos;
    }
    
    /**
     Gets number of contracted basis functions in GTOs block.

     @return the number of contracted basis functions.
     */
    auto
    getNumberOfContrGtos() const -> int32_t
    {
        if (_atomicIndexes.empty())
        {
            return 0;
        }
        else
        {
            return static_cast<int32_t>((_atomicIndexes.getMDSpan()).extent(0));
        }
    }
    
    /**
     Gets constant pointer to atomic identifiers vector.

     @return the atomic identifiers.
     */
    auto
    getAtomicIdentifiers() const -> const int32_t*
    {
        return _atomicIndexes.data();
    }
    
    /**
     Gets constant pointer to basis function indexes in full AO basis for
     specific angular momentum component.

     @param iComponent the component of angular momentum.
     @return the indexes in full AO basis.
     */
    auto
    getIdentifiers(const int32_t iComponent) const -> const int32_t*
    {
        if (iComponent < angmom::to_SphericalComponents(_angularMomentum))
        {
            return _gtoIndexes.data(iComponent);
        }
        else
        {
            return nullptr;
        }
    }
    
    /**
     Gets constant pointer to exponents of primitive Gaussian functions.

     @return the exponents of primitive Gaussian functions.
     */
    auto
    getExponents() const -> const T*
    {
        return _gtoPrimitives.data(0);
    }
    
    /**
     Gets constant pointer to normalization factors of primitive Gaussian
     functions.

     @return the normalization factors of primitive Gaussian functions.
     */
    auto
    getNormFactors() const -> const T*
    {
        return _gtoPrimitives.data(1);
    }
    
    /**
     Gets constant pointer to Cartesian X coordinates of primitive Gaussian
     functions.

     @return the exponents of primitive Gaussian functions.
     */
    auto
    getCoordinatesX() const -> const T*
    {
        return _gtoPrimitives.data(2);
    }
    
    /**
     Gets constant pointer to Cartesian Y coordinates of primitive Gaussian
     functions.

     @return the exponents of primitive Gaussian functions.
     */
    auto
    getCoordinatesY() const -> const T*
    {
        return _gtoPrimitives.data(3);
    }
    
    /**
     Gets constant pointer to Cartesian Z coordinates of primitive Gaussian
     functions.

     @return the exponents of primitive Gaussian functions.
     */
    auto
    getCoordinatesZ() const -> const T*
    {
        return _gtoPrimitives.data(4);
    }
};

#endif /* BinnedGtoBlock_hpp */
