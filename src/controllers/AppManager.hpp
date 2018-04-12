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

/**
 Class CAppManager manages jobs list creation, ordering and execution workflow.
 
 @author Z. Rinkevicius
  */
class CAppManager
{
    /**
     The execution state of application manager (true - no
     errors, false otherwise)
     */
    bool _state;
    
    /**
     The global rank of MPI process asssociate with application
     manager.
     */
    int32_t _globRank;
    
    /**
     The total number of MPI processes available to application
     manager.
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
    
public:
    
    /**
     Creates an application manager object.
     
     @param argc the number of command line arguments.
     @param argv the array of command line arguments.
     */
    CAppManager(int argc, char** argv);
    
    /**
     Destroys an application manager objec.
     */
    ~CAppManager();
    
    /**
     Executes a list of jobs assigned to application manager.
     */
    void execute();
    
    /**
     Updates execution status of application manager.

     @param state the new execution status.
     */
    void updateState(const bool state);
    
    /**
     Gets execution status of application manager.

     @return true - no errors, false - otherwise.
     */
    bool getState() const;
};

#endif /* AppManager_hpp */
