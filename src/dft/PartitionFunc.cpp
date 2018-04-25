//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "PartitionFunc.hpp"

#include "MathFunc.hpp"

namespace partfunc { // partfunc namespace

    void
    ssf(      double* gridCoordsX,
              double* gridCoordsY,
              double* gridCoordsZ,
              double* gridWeights,
        const int32_t nGridPoints,
        const double* molCoordsX,
        const double* molCoordsY,
        const double* molCoordsZ,
        const int32_t nAtoms,
              double* partWeights,
        const double  minDistanceAB, 
        const int32_t idAtom)
    {
        // SSF parameter: 0.5 * (1 - a)
        
        const double parssf = 0.180;
        
        // loop over grid points
        
        for (int32_t i = 0; i < nGridPoints; i++)
        {
            // grid coordinates
            
            double rgx = gridCoordsX[i];
            
            double rgy = gridCoordsY[i];
            
            double rgz = gridCoordsZ[i];
            
            // weights screening
            
            auto rig = mathfunc::distance(molCoordsX[idAtom],
                                          molCoordsY[idAtom],
                                          molCoordsZ[idAtom],
                                          rgx, rgy, rgz);
            
            if (rig < parssf * minDistanceAB) continue;
    
            // initialize weights
            
            mathfunc::set_to(partWeights, 1.0, nAtoms);
            
            // outer loop over atoms
            
            for (int32_t j = 0; j < nAtoms; j++)
            {
                // molecular coodinates
                
                double rax = molCoordsX[j];
                
                double ray = molCoordsY[j];
                
                double raz = molCoordsZ[j];
                
                // distance from grid point to j-th atom
                
                double rag = mathfunc::distance(rax, ray, raz, rgx, rgy, rgz);
                
                // loop over atoms
                
                for (int32_t k = j + 1; k < nAtoms; k++)
                {
                    // molecular coodinates
                    
                    double rbx = molCoordsX[k];
                    
                    double rby = molCoordsY[k];
                    
                    double rbz = molCoordsZ[k];
                    
                    // distance from grid point to k-th atom
                    
                    double rbg = mathfunc::distance(rbx, rby, rbz, rgx, rgy, rgz);
                    
                    // distance from j-th atom to k-th atom
                    
                    double rab = mathfunc::distance(rax, ray, raz, rbx, rby, rbz);
                    
                    // eliptical coordinate
                    
                    double mab = (rag - rbg) / rab;
                    
                    // scale partial weight
                    
                    partWeights[j] *= 0.5 * (1.0 - partfunc::zeta(mab));
                    
                    partWeights[k] *= 0.5 * (1.0 + partfunc::zeta(mab));
                }
            }
            
            //  adjust weight of i-th grid point
            
            gridWeights[i] *= partWeights[idAtom] / mathfunc::sum(partWeights, nAtoms);
        }
    }
    
    double
    zeta(const double eRadius)
    {
        // efefctive SSF radius
        
        const double radssf = 0.64;
        
        // lower boundary
        
        if (eRadius <= -radssf) return -1.0;
        
        // upper boundary
        
        if (eRadius >= radssf) return 1.0;
        
        // middle interval
        
        auto mab  = eRadius / radssf;
        
        auto mab2 = mab * mab;
        
        auto gab  = 0.0625 * mab * (35.0 + mab2 * (-35.0 + mab2 * (21.0 - 5.0 * mab2)));
        
        return gab;
    }
    
} // partfunc namespace
