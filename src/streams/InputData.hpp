//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef InputData_hpp
#define InputData_hpp

#include <ostream>
#include <vector>

#include "ControlGroup.hpp"

/**
 Class CInputData stores information about control groups obtained from input
 source and provides set of methods for handling control groups.
 
 @author Z. Rinkevicius
 */
class CInputData
{
    /**
     The vector of control group objects.
     */
    std::vector<CControlGroup> _controlGroups;

public:

    /**
     Creates an empty input data object.
     */
    CInputData();

    /**
     Destroys an input data object.
     */
    ~CInputData();

    /**
     Compares input data object with other input data object.
     
     @param other the input data object.
     @return true if input data objects are equal, false otherwise.
     */
    bool operator==(const CInputData& other) const;

    /**
     Compares input data object with other input data object.
     
     @param other the input data object.
     @return true if input data objects are not equal, false otherwise.
     */
    bool operator!=(const CInputData& other) const;

    /**
     Adds a control group object to vector of control group objects.

     @param controlGroup the control group object.
     */
    void addControlGroup(const CControlGroup& controlGroup);

    /**
     Determines number of control group objects in vector of control group
     objects.

     @return the number of control group objects.
     */
    int32_t getNumberOfControlGroups() const;

    /**
     Determines number of control group objects with requested name in vector
     of control group objects.

     @param nameOfControlGroup the name of control group.
     @return the number of control group objects.
     */
    int32_t getNumberOfControlGroups(const std::string& nameOfControlGroup) const;

    /**
     Determines number of control group objects with requested name in vector
     of control group objects.
     
     @param nameOfControlGroup the name of control group.
     @return the number of control group objects.
     */
    int32_t getNumberOfControlGroups(const char* nameOfControlGroup) const;

    // CControlGroup getControlGroup(const size_t indexOfControlGroup):
    //
    // Gets control group with requested index from vector of control groups.
    //
    // Input:
    // indexOfControlGroup (size_t) - the index of control group.
    //
    // Output:
    // (CControlGroup) - the control group.

    /**
     Gets requested control group object from vector of control group objects.

     @param indexOfControlGroup the index of control group object in vector of
            control group objects.
     @return the control group object.
     */
    CControlGroup getControlGroup(const size_t indexOfControlGroup) const;

    /**
     Gets requested control group with specified name from vector of control
     group objects..

     @param indexOfControlGroup the exclusive index of control group objects
            with specified name.
     @param nameOfControlGroup the name of control group object.
     @return the control group object.
     */
    CControlGroup getControlGroup(const size_t       indexOfControlGroup,
                                  const std::string& nameOfControlGroup) const;

    /**
     Gets requested control group with specified name from vector of control
     group objects..
     
     @param indexOfControlGroup the exclusive index of control group objects
            with specified name.
     @param nameOfControlGroup the name of control group object.
     @return the control group object.
     */
    CControlGroup getControlGroup(const size_t indexOfControlGroup,
                                  const char*  nameOfControlGroup) const;

    /**
     Converts input data object to formatted text and insert it into output text
     stream.
     
     @param output the output text stream.
     @param source the input data object.
     */
    friend std::ostream& operator<<(      std::ostream& output,
                                    const CInputData&   source);
};

#endif /* InputData_hpp */
