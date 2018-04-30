//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef BasisReader_hpp
#define BasisReader_hpp

#include <string>
#include <tuple>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "InputData.hpp"
#include "OutputStream.hpp"
#include "AtomBasis.hpp"

/**
 Class CBasisReader handles parsing of @basis control group and basis set
 reading from basis set library.
 
 @author Z. Rinkevicius
 */
class CBasisReader
{
    /**
     The state of basis reader object : true - no errors, false - otherwise.
     */
    bool _state;

    /**
     The name of basis set.
     */
    std::string _label;

    /**
     Reads atom basis from basis set library file for specific chemical element.
     Parsing errors are printed to output stream.

     @param idElemental the identifier of chemical element.
     @param fileName the name of basis set library file.
     @param oStream the output stream.
     @return the atom basis object.
     */
    CAtomBasis _readAtomBasis(const int32_t        idElemental,
                              const std::string&   fileName,
                                    COutputStream& oStream);

    /**
     Reads basis function data (angular momentum, number of primitive Gaussian
     functions) from input line object. Parsing errors are printed to output
     stream.

     @param inputLine the input line onject.
     @param oStream the output stream.
     @return the tuple (angular momentum, number of primitive Gaussian
             functions).
     */
    std::tuple<int32_t, int32_t> _readShellHeader(const CInputLine&    inputLine,
                                                        COutputStream& oStream);

    /**
     Reads primitive Gausian function data (exponent, normalization factor) from
     input line object. Parsing errors are printed to output stream.

     @param inputLine the input line object.
     @param oStream the output stream.
     @return the tuple (exponent, normalization factor).
     */
    std::tuple<double, double> _readPrimitveBasisFuction(const CInputLine& inputLine,
                                                               COutputStream& oStream);

    /**
     Prints multiple definitions error message for @basis control group to
     output stream and sets basis reader state to abnormal.

     @param oStream the output stream.
     */
    void _errorUniqueGroup(COutputStream& oStream);

    /**
     Prints unsupported basis set error message to output stream and sets basis
     reader object state to abnormal.

     @param oStream the output stream.
     */
    void _errorBasisSet(COutputStream& oStream);

    /**
     Prints corrupted basis set file error message to output stream and sets
     basis reader object state to abnormal.

     @param oStream the output stream.
     */
    void _errorCorruptedBasisSet(COutputStream& oStream);

    /**
     Prints syntax error massage for @basis control group to output stream and
     sets basis reader object state to abnormal.

     @param oStream the output stream.
     */
    void _syntaxBasisSet(COutputStream& oStream);

public:

    /**
     Creates a basis reader object.
     */
    CBasisReader();

    /**
     Destroys a basis reader object.
     */
    ~CBasisReader();

    /**
     Get a state of basis reader object.
     
     @return true if no errors, false otherwise.
     */
    bool getState() const;

    /**
     Parses @basis control group from input data and sets internal data of basis
     set reader object. Parsing errors are printed to output stream.

     @param inputData the input data object.
     @param oStream the output stream.
     */
    void parse(const CInputData&   inputData,
                     COutputStream& oStream);

    /**
     Creates molecular basis object by reading basis set data from basis library
     for selected molecule. Reading errors are printed to output stream.

     @param pathToBasisSets the path to basis set library.
     @param molecule the molecule.
     @param oStream the output stream.
     @return the molecular basis object.
     */
    CMolecularBasis getAOBasis(const std::string&   pathToBasisSets,
                               const CMolecule&     molecule,
                                     COutputStream& oStream);
};

#endif /* BasisReader_hpp */
