//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef CMMParameters_hpp
#define CMMParameters_hpp

#include <cstdint>
#include <ostream>
#include <vector>

#include "CMMType.hpp"

/**
 Class CMMParameters stores CMM model parameters for chemical element.
 
 @author Z. Rinkevicius
 */
class CCMMParameters
{
    /**
     The parameterization model.
     */
    cmmtyp _paramModel;
    
    /**
     The exponent of Gaussian charge function.
     */
    double _exponent;

    /**
     The reference coordination number of chemical element.
     */
    double _refCoordNumber;
    
    /**
     The scaling factor of coordination numbers difference.
     */
    double _factCoordNumber;
    
    /**
     The vector of Lorentzian frequencies.
     */
    std::vector<double> _refFrequencies;
    
    /**
     The vector of Lorentzian damping factors.
     */
    std::vector<double> _refGammas;
    
    /**
     The vector of Lorentzian scaling factors.
     */
    std::vector<double> _refScaleFactors;

public:

    /**
     Creates an empty CMM parameters object.
     */
    CCMMParameters();

    /**
     Creates a CMM parameters object.
     
     @param paramModel the parameterization of CMM model.
     @param exponent the exponent of Gaussian charge function.
     @param refCoordNumber the reference coordination number.
     @param factCoordNumber the coordination numbers difference scaling factor.
     @param refFrequencies the vector of Lorentzian frequencies.
     @param refGammas the vector of Lorentzian damping factors.
     @param refScaleFactors the vector of Lorentzian scaling factors.
     */
    CCMMParameters(const cmmtyp               paramModel,
                   const double               exponent,
                   const double               refCoordNumber,
                   const double               factCoordNumber,
                   const std::vector<double>& refFrequencies,
                   const std::vector<double>& refGammas,
                   const std::vector<double>& refScaleFactors); 
                   
    /**
     Creates a CMM parameters object by copying other CMM parameters object.
     
     @param source the CMM parameters object.
     */
    CCMMParameters(const CCMMParameters& source);

    /**
     Creates a CMM parameters object by moving other CMM parameters object.
     
     @param source the CMM parameters object.
     */
    CCMMParameters(CCMMParameters&& source) noexcept;

    /**
     Destroys a CMM parameters object.
     */
    ~CCMMParameters();

    /**
     Assigns a CMM parameters object by copying other CMM parameters object.
     
     @param source the CMM parameters object.
     */
    CCMMParameters& operator=(const CCMMParameters& source);

    /**
     Assigns a CMM parameters object by moving other CMM parameters object.
     
     @param source the CMM parameters object.
     */
    CCMMParameters& operator=(CCMMParameters&& source) noexcept;

    /**
     Compares CMM parameters object with other CMM parameters object.
     
     @param other the CMM parameters object.
     @return true if CMM parameters objects are equal, false otherwise.
     */
    bool operator==(const CCMMParameters& other) const;

    /**
     Compares CMM parameters object with other CMM parameters object.
     
     @param other the CMM parameters object.
     @return true if CMM parameters objects are not equal, false otherwise.
     */
    bool operator!=(const CCMMParameters& other) const;

    /**
     Converts CMM parameters object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the CMM parameters object.
     */
    friend std::ostream& operator<<(      std::ostream&   output,
                                    const CCMMParameters& source);
};

#endif /* CMMParameters_hpp */
