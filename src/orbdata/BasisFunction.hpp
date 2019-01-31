//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef BasisFunction_hpp
#define BasisFunction_hpp

#include <cstdint>
#include <ostream>
#include <vector>

#include "mpi.h"

/**
 Class CBasisFunction stores data about single contracted basis function and
 provides set of methods for handling of basis function data.
 
 @author Z. Rinkevicius
 */
class CBasisFunction
{
    /**
     The vector of exponents of primitive Gaussian functions.
     */
    std::vector<double> _exponents;

    /**
     The vector of normalization factors of primitive Gaussian functions.
     */
    std::vector<double> _normFactors;

    /**
     The angular momentum of basis function.
     */
    int32_t _angularMomentum;

    /**
     Rescales normalization factors to match normalization of spherical (l,0)
     component of basis function.
     */
    void _rescale();

    /**
     Computes overlap between two primitive Gaussian functions.

     @param iComponent the index of first primitve Gaussain function.
     @param jComponent the index of second primitve Gaussain function.
     @return the overlap between two primitive Gaussian functions.
     */
    double _overlap(const size_t iComponent,
                    const size_t jComponent) const;

public:

    /**
     Creates an empty basis function object.
     */
    CBasisFunction();

    /**
     Creates a basis function object.
     
     @param exponents the vector of exponents of primitive Gaussian functions.
     @param normFactors the vector of normalization factors of primitive
            Gaussian functions.
     @param angularMomentum the angular momentum of basis function.
     */
    CBasisFunction(const std::vector<double>& exponents,
                   const std::vector<double>& normFactors,
                   const int32_t              angularMomentum);

    /**
     Creates a basis function object by copying other basis function object.
     
     @param source the basis function object.
     */
    CBasisFunction(const CBasisFunction& source);

    /**
     Creates a basis function object by moving other basis function object.
     
     @param source the basis function object.
     */
    CBasisFunction(CBasisFunction&& source) noexcept;

    /**
     Destroys a basis function object.
     */
    ~CBasisFunction();

    /**
     Assigns a basis function object by copying other basis function object.
     
     @param source the basis function object.
     */
    CBasisFunction& operator=(const CBasisFunction& source);

    /**
     Assigns a basis function object by moving other basis function object.
     
     @param source the basis function object.
     */
    CBasisFunction& operator=(CBasisFunction&& source) noexcept;

    /**
     Compares basis function object with other basis function object.
     
     @param other the basis function object.
     @return true if basis function objects are equal, false otherwise.
     */
    bool operator==(const CBasisFunction& other) const;

    /**
     Compares basis function object with other basis function object.
     
     @param other the basis function object.
     @return true if basis function objects are not equal, false otherwise.
     */
    bool operator!=(const CBasisFunction& other) const;

    /**
     Sets exponents of primittive Gaussian functions with specific vector of
     exponents.

     @param exponents the vector of exponents.
     */
    void setExponents(const std::vector<double>& exponents);

    /**
     Sets normalization factors of primitive Gaussian functions with specific
     vector of normalization factors.

     @param normFactors the vector of normalization factors.
     */
    void setNormalizationFactors(const std::vector<double>& normFactors);

    /**
     Set angular momentum of basis function.

     @param angularMomentum the angular momentum.
     */
    void setAngularMomentum(const int32_t angularMomentum);

    /**
     Adds primittive Gaussian function to basis function.

     @param exponent the exponent of primitive Gaussian function.
     @param normFactor the normalization factor of primitive Gaussian function.
     */
    void add(const double exponent,
             const double normFactor);

    /**
     Normalizes basis function.
     */
    void normalize();

    /**
     Gets vector of exponents of primitive Gaussian functions.

     @return the vector of exponents.
     */
    std::vector<double> getExponents() const;

    /**
     Gets vector of normalization factors of primitive Gaussian functions.

     @return the vector of normalization factors.
     */
    std::vector<double> getNormalizationFactors() const;

    /**
     Gets angular momentum of basis function.

     @return the angular momentum.
     */
    int32_t getAngularMomentum() const;
    
    /**
     Gets number of primitive Gaussian functions in basis function.

     @return the number of primitive Gaussian functions.
     */
    int32_t getNumberOfPrimitiveFunctions() const;
    
    /**
     Broadcasts basis function object within domain of MPI communicator.

     @param rank the rank of MPI process.
     @param comm the MPI communicator.
     */
    void broadcast(int32_t  rank,
                   MPI_Comm comm);

    /**
     Converts basis function object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the basis function object.
     */
    friend std::ostream& operator<<(      std::ostream&   output,
                                    const CBasisFunction& source);
};

#endif /* BasisFunction_hpp */
