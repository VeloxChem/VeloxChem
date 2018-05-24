//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef EnvironmentReader_hpp
#define EnvironmentReader_hpp

#include <string>

#include "InputData.hpp"
#include "OutputStream.hpp"

/**
 Class CEnvironmentReader handles parsing of @progenv control group, which
 contains definitions of environmental variables.
 
 @author Z. Rinkevicius
 */
class CEnvironmentReader
{
    /**
     The state of environment reader object : true - no errors, false -
     otherwise.
     */
    bool _state;

    /**
     The path to basis set library.
     */
    std::string _pathToBasisSets;
    
    /**
     The path to force fields library.
     */
    std::string _pathToForceFields;

    /**
     Prints multiple definitions error message for @progenv control group to
     output stream and sets environment reader object state to abnormal.

     @param nGroups the number of @progenv control groups.
     @param oStream the output stream.
     */
    void _errorUniqueGroup(const size_t   nGroups,
                           COutputStream& oStream);

    /**
     Prints unknown environment variable error messeage to output stream and
     sets environment reader object state to abnormal.

     @param inputLine the input line object with unknown environmental variable.
     @param oStream the output stream.
     */
    void _errorUnknownVariable(const CInputLine&    inputLine,
                                     COutputStream& oStream);

    /**
     Reads path to basis set library from @progenv control group.

     @param inputLine the input line object with path to basis set library.
     @param oStream the output stream.
     @return true if parsing of path to basis set library is successful, false
             otherwise.
     */
    bool _addPathToBasisSets(const CInputLine&    inputLine,
                                   COutputStream& oStream);
    
    /**
     Reads path to force fields library from @progenv control group.
     
     @param inputLine the input line object with path to force fields library.
     @param oStream the output stream.
     @return true if parsing of path to force fields library is successful,
             false otherwise.
     */
    bool _addPathToForceFields(const CInputLine&    inputLine,
                                     COutputStream& oStream);

    /**
     Prints syntax error message for definition of path to basis set library to
     output stream and sets environment reader object state to abnormal.

     @param inputLine the input line object with syntax error.
     @param oStream the output stream.
     */
    void _syntaxBasisLibrary(const CInputLine&    inputLine,
                                   COutputStream& oStream);
    
    /**
     Prints syntax error message for definition of path to force fields library
     to output stream and sets environment reader object state to abnormal.
     
     @param inputLine the input line object with syntax error.
     @param oStream the output stream.
     */
    void _syntaxForceFieldsLibrary(const CInputLine&    inputLine,
                                         COutputStream& oStream);

public:

    /**
     Creates an environment reader object.
     */
    CEnvironmentReader();

    /**
     Destroys an environment reader object.
     */
    ~CEnvironmentReader();

    /**
     Get a state of environment reader object.

     @return true if no errors, false otherwise.
     */
    bool getState() const;

    /**
     Reads and parses @progenv control group from input data object. Parsing
     errors are printed to output stream.

     @param inputData the input data object.
     @param oStream the output stream.
     */
    void parse(const CInputData&    inputData,
                     COutputStream& oStream);

    /**
     Gets path to basis set library defined in @progenv group.

     @return the path to basis set library.
     */
    std::string getPathToBasisSets() const;
    
    /**
     Gets path to force fields library defined in @progenv group.
     
     @return the path to force fields library.
     */
    std::string getPathToForceFields() const;
};

#endif /* EnvironmentReader_hpp */
