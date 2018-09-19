//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MolecularGrid_hpp
#define MolecularGrid_hpp

#include <cstdint>
#include <ostream>

#include "MemBlock2D.hpp"

/**
 Class CMolecularGrid class generates molecular grid.
 
 @author Z. Rinkevicius
 */
class CMolecularGrid
{
    /**
     The distribution status of molecular grid object.
     */
    bool _isDistributed;

    /**
     The grid points (coordinates, weights).
     */
    CMemBlock2D<double> _gridPoints;

public:
    
    /**
     Creates an empty molecular grid object.
     */
    CMolecularGrid();

    /**
     Creates a molecular grid object.
     
     @param gridPoints the 2D memory block object with grid points data.
     */
    CMolecularGrid(const CMemBlock2D<double>& gridPoints);

    /**
     Creates a molecular grid object by copying other molecular grid object.
     
     @param source the molecular grid object.
     */
    CMolecularGrid(const CMolecularGrid& source);

    /**
     Creates a molecular grid object by by moving other molecular grid object.
     
     @param source the molecular grid object.
     */
    CMolecularGrid(CMolecularGrid&& source) noexcept;

    /**
     Destroys a molecular grid object.
     */
    ~CMolecularGrid();

    /**
     Assigns a molecular grid object by copying other molecular grid object.
     
     @param source the molecular grid object.
     */
    CMolecularGrid& operator=(const CMolecularGrid& source);

    /**
     Assigns a molecular grid object by moving other molecular grid object.
     
     @param source the molecular grid object.
     */
    CMolecularGrid& operator=(CMolecularGrid&& source) noexcept;

    /**
     Compares molecular grid object with other molecular grid object.
     
     @param other the molecular grid object.
     @return true if molecular grid objects are equal, false otherwise.
     */
    bool operator==(const CMolecularGrid& other) const;

    /**
     Compares molecular grid object with other molecular grid object.
     
     @param other the molecular grid object.
     @return true if molecular grid objects are not equal, false otherwise.
     */
    bool operator!=(const CMolecularGrid& other) const;

    /**
     Gets number of grid points in molecular grid object.

     @return the number of grid points.
     */
    int32_t getNumberOfGridPoints() const;
    
    /**
     Gets Cartesian X coordinates of grid points in molecular grid object.

     @return he pointer to Cartesian X coordinates of grid points.
     */
    const double* getCoordinatesX() const;
    
    /**
     Gets Cartesian Y coordinates of grid points in molecular grid object.
     
     @return the pointer to Cartesian Y coordinates of grid points.
     */
    const double* getCoordinatesY() const;
    
    /**
     Gets Cartesian Z coordinates of grid points in molecular grid object.
     
     @return the pointer to Cartesian Z coordinates of grid points.
     */
    const double* getCoordinatesZ() const;
    
    /**
     Gets weights of grid points in molecular grid object.
     
     @return the pointer to weights of grid points.
     */
    const double* getWeights() const;

    /**
     Distributes grid points data across a molecular grid objects associated
     with MPI processes within domain of MPI communacator and sets distribution
     flag to true.

     @param rank the rank of MPI process.
     @param nodes the number of nodes in MPI domain.
     @param comm he MPI communicator.
     */
    void distribute(int32_t  rank,
                    int32_t  nodes,
                    MPI_Comm comm);

    /**
     Converts molecular grid object to text and insert it into output text
     stream.

     @param output the output text stream.
     @param source the molecular grid.
     */
    friend std::ostream& operator<<(      std::ostream&   output,
                                    const CMolecularGrid& source);
};

#endif /* MolecularGrid_hpp */
