//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef AppManager_hpp
#define AppManager_hpp

#include <cstdint>
#include <string>

#include "OutputStream.hpp"
#include "SystemClock.hpp"
#include "InputData.hpp"

/**
 Class CAppManager handles jobs manager, ordering and execution workflow in
 main() function of program.
 
 @author Z. Rinkevicius
  */
class CAppManager
{
    /**
     The execution state of application manager (true - no errors, false
     otherwise).
     */
    bool _state;
    
    /**
     The global rank of MPI process asssociate with application manager.
     */
    int32_t _globRank;
    
    /**
     The total number of MPI processes available to application manager.
     */
    int32_t _globNodes;
    
    /**
     The name of input file used to read input data.
     */
    std::string _iFilename;
    
    /**
     The name of output file used to write output data.
     */
    std::string _oFilename;
    
    /**
     The path to basis set library.
     */
    std::string _pathToBasLib;
    
    /**
     The application manager object's timer.
     */
    CSystemClock _sysClock;
    
    /**
     Prints Velox Chem MP start header to output stream.

     @param oStream the output stream.
     */
    void _printStartHeader(COutputStream& oStream);
    
    /**
     Prints Velox Chem MP finish header to output stream.

     @param oStream the output stream.
     */
    void _printFinishHeader(COutputStream& oStream);
    
    /**
     Sets environmental variables by reading @progenv control group from input
     data object. Errors are printed to output stream.

     @param inputData the input data object.
     @param oStream the output stream.
     */
    void _setEnvironment(const CInputData&    inputData,
                               COutputStream& oStream);
    
public:
    
    /**
     Creates an application manager object.
     
     @param argc the number of command line arguments.
     @param argv the array of command line arguments.
     */
    CAppManager(int    argc,
                char** argv);
    
    /**
     Destroys an application manager object.
     */
    ~CAppManager();
    
    /**
     Executes a list of jobs assigned to application manager object.
     */
    void execute();
    
    /**
     Updates execution status of application manager object.

     @param state the new execution status.
     */
    void updateState(const bool state);
    
    /**
     Gets execution status of application manager object.

     @return true - no errors, false - otherwise.
     */
    bool getState() const;
};

#endif /* AppManager_hpp */
