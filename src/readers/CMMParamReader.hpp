//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef CMMParamReader_hpp
#define CMMParamReader_hpp

#include <string>
#include <vector>
#include <tuple>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "InputData.hpp"
#include "OutputStream.hpp"
#include "CMMParameters.hpp"

/**
 Class CCMMParamReader handles parsing of @cmmparam control group and CMM
 parameters reading from force fields library.
 
 @author Z. Rinkevicius
 */
class CCMMParamReader
{
    /**
     The state of basis reader object : true - no errors, false - otherwise.
     */
    bool _state;

    /**
     The name of force field parameters file.
     */
    std::string _label;
    
    /**
     Reads atom CMM force field from force fields library file for specific
     chemical element. Parsing errors are printed to output stream.
     
     @param idElemental the identifier of chemical element.
     @param fileName the name of force field library file.
     @param oStream the output stream.
     @return the atom force field parameters.
     */
    CCMMParameters _readAtomForceField(const int32_t        idElemental,
                                       const std::string&   fileName,
                                             COutputStream& oStream);
    
    
    /**
     Prints syntax error massage for @cmmparam control group to output stream
     and sets CMM parameters reader object state to abnormal.
     
     @param oStream the output stream.
     */
    void _syntaxForceField(COutputStream& oStream);
    
    /**
     Prints multiple definitions error message for @cmmparam control group to
     output stream and sets CMM parameters reader state to abnormal.
     
     @param oStream the output stream.
     */
    void _errorUniqueGroup(COutputStream& oStream);
    
    /**
     Prints unsupported force field error message to output stream and sets CMM
     parameters reader object state to abnormal.
     
     @param oStream the output stream.
     */
    void _errorForceField(COutputStream& oStream);
    
    /**
     Prints corrupted force field file error message to output stream and sets
     CMM parameters reader object state to abnormal.
     
     @param oStream the output stream.
     */
    void _errorCorruptedForceField(COutputStream& oStream);
    
    /**
     Reads force field parameters (type, alpha, exponent, coordination number,
     coordination scaling factor) from input line object. Parsing errors are
     printed to output stream.

     @param ffParameters the CMM parameters object.
     @param inputLine the input line object.
     @param oStream the output stream.
     */
    void _readForceFieldHeader(      CCMMParameters& ffParameters,
                               const CInputLine&     inputLine,
                                     COutputStream&  oStream);
    
    /**
     Reads number of Lorentzians in dynamic atomic polarizability expansion.

     @param inputLine the input line object.
     @param oStream the output stream.
     @return the number of Lorentzians.
     */
    int32_t _readNumberOfLorentzians(const CInputLine&    inputLine,
                                           COutputStream& oStream);
    
    /**
     Reads Lorentzian parameters and adds it to vector of Lorentzians.

     @param ffParameters the CMM parameters object.
     @param inputLine the input line object.
     @param oStream the output stream.
     */
    void _readLorentzian(      CCMMParameters& ffParameters,
                         const CInputLine&     inputLine,
                               COutputStream&  oStream);
    
public:

    /**
     Creates a CMM parameters reader object.
     */
    CCMMParamReader();

    /**
     Destroys a CMM parameters reader object.
     */
    ~CCMMParamReader();

    /**
     Get a state of CMM parameters reader object.
     
     @return true if no errors, false otherwise.
     */
    bool getState() const;

    /**
     Parses @cmmparam control group from input data and sets internal data of
     CMM parameters reader object. Parsing errors are printed to output stream.

     @param inputData the input data object.
     @param oStream the output stream.
     */
    void parse(const CInputData&   inputData,
                     COutputStream& oStream);

    /**
     Creates vector of CMM parameters objects by reading force field data from
     force fields library for selected molecule. Reading errors are printed to
     output stream.

     @param pathToForceFields the path to force fields set library.
     @param molecule the molecule.
     @param oStream the output stream.
     @return the vector of CMM parameters objects.
     */
    std::vector<CCMMParameters> getCMMParameters(const std::string&   pathToForceFields,
                                                 const CMolecule&     molecule,
                                                       COutputStream& oStream);
};

#endif /* CMMParamReader_hpp */
