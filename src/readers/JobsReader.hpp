//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef JobsReader_hpp
#define JobsReader_hpp

#include <vector>

#include "InputData.hpp"
#include "BaseJob.hpp"
#include "OutputStream.hpp"
#include "ExecMode.hpp"

/**
 Class CJobsReader handles parsing of @jobs control group, which
 contains list of jobs to be executed during calculations.
 
 @author Z. Rinkevicius
 */
class CJobsReader
{
    /**
     The state of jobs reader object : true - no errors, false - otherwise.
     */
    bool _state;
    
    /**
     The parsed execution mode of jobs.
     */
    execmode _runMode;
    
    /**
     Prints multiple definitions error message for @jobs control group to
     output stream and sets jobs reader object state to abnormal.

     @param nGroups the number of @jobs control groups.
     @param oStream the output stream.
     */
    void _errorUniqueGroup(const size_t         nGroups,
                                 COutputStream& oStream);

    /**
     Prints unknown job type error message to output stream and sets jobs
     reader object state to abnormal.

     @param inputLine the input line object with unknown job type.
     @param oStream the output stream.
     */
    void _errorUnknownJobType(const CInputLine&    inputLine,
                                    COutputStream& oStream);
    
    /**
     Reads a execution mode by parsing input line object. Errors are printed to
     output stream.
     
     @param inputLine the input line object.
     @param oStream the output stream.
     @return true if execution mode is read, false otherwise.
     */
    bool _addExecutionMode(const CInputLine&    inputLine,
                                 COutputStream& oStream);

    /**
     Adds a single point job identifier to vector of job identifiers by
     parsing input line object. Errors are printed to output stream.

     @param listOfJobIds the vector of job identifiers.
     @param inputLine the input line object.
     @param oStream the output stream.
     @return true if identifier is added to vector of job identifiers, false
             otherwise.
     */
    bool _addSinglePoint(      std::vector<int32_t>& listOfJobIds,
                         const CInputLine&           inputLine,
                               COutputStream&        oStream);
    
    /**
     Adds a optimization job identifier to vector of job identifiers by parsing
     input line object. Errors are printed to output stream.

     @param listOfJobIds the vector of job identifiers.
     @param inputLine the input line object.
     @param oStream the output stream.
     @return true if identifier is added to vector of job identifiers, false
             otherwise.
     */
    bool _addOptimization(      std::vector<int32_t>& listOfJobIds,
                          const CInputLine&           inputLine,
                                COutputStream&        oStream);

    /**
     Prints syntax error message for run mode selection to output stream and
     sets jobs reader object state to abnormal.
     
     @param inputLine the input line object with syntax error.
     @param oStream the output stream.
     */
    void _syntaxRunMode(const CInputLine&    inputLine,
                              COutputStream& oStream);
    
    /**
     Prints syntax error message for single point job to output stream and sets
     jobs reader object state to abnormal.

     @param inputLine the input line object with syntax error.
     @param oStream the output stream.
     */
    void _syntaxSinglePoint(const CInputLine&    inputLine,
                                  COutputStream& oStream);
    
    /**
     Prints syntax error message for definition of optimization job to output
     stream and sets jobs reader object state to abnormal.

     @param inputLine the input line object with syntax error.
     @param oStream the output stream.
     */
    void _syntaxOptimization(const CInputLine&    inputLine,
                                   COutputStream& oStream);

    /**
     Prints unknown calculation type error message to output stream and sets
     jobs reader object state to abnormal.

     @param calcType the unknown calculation type .
     @param inputLine the input line object with unknown calculation type.
     @param oStream the output stream.
     */
    void _errorUnknownCalculationType(const char*          calcType,
                                      const CInputLine&    inputLine,
                                            COutputStream& oStream);
public:

    /**
     Creates a jobs reader object.
     */
    CJobsReader();

    /**
     Destroys a jobs reader object.
     */
    ~CJobsReader();

    /**
     Gets a state of jobs reader object.

     @return true if no errors, false otherwise.
     */
    bool getState() const;
    
    /**
     Gets a parsed run mode of jobs.

     @return key value for run mode.
     */
    execmode getRunMode() const;

    /**
     Creates vector of job identifiers by parsing input data object. Parsing
     errors are printed to output stream.

     @param listOfJobIds the vector of jobs identifiers.
     @param inputData the input data object.
     @param oStream the output stream.
     */
    void parse(      std::vector<int32_t>& listOfJobIds,
               const CInputData&           inputData,
                     COutputStream&        oStream);
};

#endif /* JobsReader_hpp */
