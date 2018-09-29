//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ScreeningContainer_hpp
#define ScreeningContainer_hpp

#include <vector>

#include "CauchySchwarzScreener.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "GtoPairsContainer.hpp"
#include "VecMemBlocks.hpp"

/**
Class CScreeningContainer class stores precomputed electron repulsion integrals
screener objects for bra and ket sides and provides set of methods for
manpiluting screener objects.

@author Z. Rinkevicius
*/
class CScreeningContainer
{
    /**
     The vector of screener object for bra/ket side.
     */
    std::vector<CCauchySchwarzScreener> _screeners;
    
public:
    
    /**
     Creates an empty screening container object.
     */
    CScreeningContainer();
    
    /**
     Creates a screening container object.
     
     @param braQValues the vector Q values on bra side.
     @param ketQValues the vector Q values on ket side.
     @param braGtoPairsContainer the GTOs pairs container on bra side.
     @param ketGtoPairsContainer the GTOs pairs container on ket side.
     @param screeningScheme the screening scheme.
     @param threshold the screening threshold.
     */
    CScreeningContainer(const CVecMemBlock<double>& braQValues,
                        const CVecMemBlock<double>& ketQValues,
                        const CGtoPairsContainer&   braGtoPairsContainer,
                        const CGtoPairsContainer&   ketGtoPairsContainer,
                        const ericut                screeningScheme,
                        const double                threshold);
    
    /**
     Creates a screening container object by copying other screening container
     object.
     
     @param source the screening container object.
     */
    CScreeningContainer(const CScreeningContainer& source);
    
    /**
     Creates a screening container object by moving other screening container
     object.
     
     @param source the screening container object.
     */
    CScreeningContainer(CScreeningContainer&& source) noexcept;
    
    /**
     Destroys a screening container object.
     */
    ~CScreeningContainer();
    
    /**
     Assigns a screening container object by copying other screening container
     object.
     
     @param source the screening container object.
     */
    CScreeningContainer& operator=(const CScreeningContainer& source);
    
    /**
     Assigns a screening container object by moving other screening container
     object.
     
     @param source the screening container object.
     */
    CScreeningContainer& operator=(CScreeningContainer&& source) noexcept;
    
    /**
     Compares screening container object with other screening container
     object.
     
     @param other the screening container object.
     @return true if screening container objects are equal, false otherwise.
     */
    bool operator==(const CScreeningContainer& other) const;
    
    /**
     Compares screening container object with other screening container
     object.
     
     @param other the screening container object.
     @return true if screening container objects are not equal, false otherwise.
     */
    bool operator!=(const CScreeningContainer& other) const;
    
    /**
     Converts screening container object to text output and insert it into
     output text stream.
     
     @param output the output text stream.
     @param source the screening container object.
     */
    friend std::ostream& operator<<(      std::ostream&        output,
                                    const CScreeningContainer& source);
};

#endif /* ScreeningContainer_hpp */
