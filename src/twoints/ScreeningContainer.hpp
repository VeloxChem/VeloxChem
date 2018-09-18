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
    
    /**
     Initializes vector of screener objects for bra and ket sides.

     @param braGtoPairsContainer the GTOs pairs container on bra side.
     @param ketGtoPairsContainer the GTOs pairs container on ket side.
     @param screeningScheme the screening scheme.
     @param threshold the screening threshold.
     */
    void _setScreeners(const CGtoPairsContainer& braGtoPairsContainer,
                       const CGtoPairsContainer& ketGtoPairsContainer,
                       const ericut              screeningScheme,
                       const double              threshold);
    
public:
    
    /**
     Creates an empty screening container object.
     */
    CScreeningContainer();
    
    /**
     Creates a screening container object.
     
     @param molecule the molecule.
     @param aoBasis the molecular AO basis.
     @param screeningScheme the screening scheme.
     @param threshold the screening threshold.
     */
    CScreeningContainer(const CMolecule&       molecule,
                        const CMolecularBasis& aoBasis,
                        const ericut           screeningScheme,
                        const double           threshold);
    
    /**
     Destroys a screening container object.
     */
    ~CScreeningContainer();
    
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
