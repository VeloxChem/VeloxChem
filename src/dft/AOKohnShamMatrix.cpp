//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "AOKohnShamMatrix.hpp"

CAOKohnShamMatrix::CAOKohnShamMatrix()
    : _xcMatrices(std::vector<CDenseMatrix>())

    , _xcRestricted(true)

    , _xcElectrons(0.0)

    , _xcEnergy(0.0)
{
    
}

CAOKohnShamMatrix::CAOKohnShamMatrix(const std::vector<CDenseMatrix>& xcMatrices,
                                     const bool                       xcRestricted,
                                     const double                     xcElectrons,
                                     const double                     xcEnergy)

    : _xcMatrices(xcMatrices)

    , _xcRestricted(xcRestricted)

    , _xcElectrons(xcElectrons)

    , _xcEnergy(xcEnergy)
{
    
}

CAOKohnShamMatrix::CAOKohnShamMatrix(const int32_t nRows,
                                     const int32_t nColumns,
                                     const bool    xcRestricted)
{
    _xcRestricted = xcRestricted;
    
    _xcElectrons = 0.0;
    
    _xcEnergy = 0.0;
    
    _xcMatrices.push_back(CDenseMatrix(nRows, nColumns));
    
    if (!_xcRestricted) _xcMatrices.push_back(CDenseMatrix(nRows, nColumns));
}

CAOKohnShamMatrix::CAOKohnShamMatrix(const CAOKohnShamMatrix& source)

    : _xcMatrices(source._xcMatrices)

    , _xcRestricted(source._xcRestricted)

    , _xcElectrons(source._xcElectrons)

    , _xcEnergy(source._xcEnergy)
{
    
}

CAOKohnShamMatrix::CAOKohnShamMatrix(CAOKohnShamMatrix&& source) noexcept

    : _xcMatrices(std::move(source._xcMatrices))

    , _xcRestricted(std::move(source._xcRestricted))

    , _xcElectrons(std::move(source._xcElectrons))

    , _xcEnergy(std::move(source._xcEnergy))
{
    
}

CAOKohnShamMatrix::~CAOKohnShamMatrix()
{
    
}

CAOKohnShamMatrix&
CAOKohnShamMatrix::operator=(const CAOKohnShamMatrix& source)
{
    if (this == &source) return *this;
    
    _xcMatrices = source._xcMatrices;
    
    _xcRestricted = source._xcRestricted;
    
    _xcElectrons = source._xcElectrons;
    
    _xcEnergy = source._xcEnergy;
    
    return *this;
}

CAOKohnShamMatrix&
CAOKohnShamMatrix::operator=(CAOKohnShamMatrix&& source) noexcept
{
    if (this == &source) return *this;
    
    _xcMatrices = std::move(source._xcMatrices);
    
    _xcRestricted = std::move(source._xcRestricted);
    
    _xcElectrons = std::move(source._xcElectrons);
    
    _xcEnergy = std::move(source._xcEnergy);
    
    return *this;
}

bool
CAOKohnShamMatrix::operator==(const CAOKohnShamMatrix& other) const
{
    if (_xcMatrices.size() != other._xcMatrices.size()) return false;
    
    for (size_t i = 0; i < _xcMatrices.size(); i++)
    {
        if (_xcMatrices[i] != other._xcMatrices[i]) return false;
    }
    
    if (_xcRestricted != other._xcRestricted) return false;
    
    if (std::fabs(_xcElectrons - other._xcElectrons) > 1.0e-13) return false;
    
    if (std::fabs(_xcEnergy - other._xcEnergy) > 1.0e-13) return false;
    
    return true;
}

bool
CAOKohnShamMatrix::operator!=(const CAOKohnShamMatrix& other) const
{
    return !(*this == other);
}

void
CAOKohnShamMatrix::zero()
{
    for (size_t i = 0; i < _xcMatrices.size(); i++)
    {
        _xcMatrices[i].zero();
    }
}

void
CAOKohnShamMatrix::symmetrize()
{
    for (size_t i = 0; i < _xcMatrices.size(); i++)
    {
        _xcMatrices[i].symmetrize();
    }
}

void
CAOKohnShamMatrix::setNumberOfElectrons(const double xcElectrons)
{
    _xcElectrons = xcElectrons;
}

void
CAOKohnShamMatrix::setExchangeCorrelationEnergy(const double xcEnergy)
{
    _xcEnergy = xcEnergy;
}

void
CAOKohnShamMatrix::reduce_sum(int32_t  rank,
                              int32_t  nodes,
                              MPI_Comm comm)
{
    for (size_t i = 0; i < _xcMatrices.size(); i++)
    {
        _xcMatrices[i].reduce_sum(rank, nodes, comm);
        
        MPI_Barrier(comm);
    }
    
    auto fsum = mpi::reduce_sum(_xcElectrons, comm);
    
    _xcElectrons = fsum;
    
    fsum = mpi::reduce_sum(_xcEnergy, comm);
    
    _xcEnergy = fsum; 
}

bool
CAOKohnShamMatrix::isRestricted() const
{
    return _xcRestricted;
}

double
CAOKohnShamMatrix::getNumberOfElectrons() const
{
    return _xcElectrons;
}

double
CAOKohnShamMatrix::getExchangeCorrelationEnergy() const
{
    return _xcEnergy; 
}

int32_t
CAOKohnShamMatrix::getNumberOfRows() const
{
    if (!_xcMatrices.empty())
    {
        return _xcMatrices[0].getNumberOfRows();
    }
    
    return 0;
}

int32_t
CAOKohnShamMatrix::getNumberOfColumns() const
{
    if (!_xcMatrices.empty())
    {
        return _xcMatrices[0].getNumberOfColumns();
    }
    
    return 0;
}

int32_t
CAOKohnShamMatrix::getNumberOfElements() const
{
    if (!_xcMatrices.empty())
    {
        return _xcMatrices[0].getNumberOfElements();
    }
    
    return 0;
}

const double*
CAOKohnShamMatrix::getKohnSham(const bool beta) const
{
    if (!_xcMatrices.empty())
    {
        if (_xcRestricted)
        {
            return _xcMatrices[0].values();
        }
        else if (!beta)
        {
            return _xcMatrices[0].values();
        }
        else
        {
            return _xcMatrices[1].values();
        }
    }
    
    return nullptr;
}

double*
CAOKohnShamMatrix::getKohnSham(const bool beta)
{
    if (!_xcMatrices.empty())
    {
        if (_xcRestricted)
        {
            return _xcMatrices[0].values();
        }
        else if (!beta)
        {
            return _xcMatrices[0].values();
        }
        else
        {
            return _xcMatrices[1].values();
        }
    }
    
    return nullptr;
}

const CDenseMatrix&
CAOKohnShamMatrix::getReferenceToKohnSham(const bool beta) const
{
    if (_xcRestricted)
    {
        return _xcMatrices[0];
    }
    else if (!beta)
    {
        return _xcMatrices[0];
    }
    else
    {
        return _xcMatrices[1];
    }
}

std::string
CAOKohnShamMatrix::getString() const
{
    if (_xcMatrices.empty()) return std::string();
 
    std::string ksmat_str;
    
    ksmat_str += "Is restricted: ";
    
    ksmat_str += (_xcRestricted) ? "Yes" : "No";

    ksmat_str += "Exchange-correlation energy: " + std::to_string(_xcEnergy) + "\n";
    
    if (_xcRestricted)
    {
        ksmat_str += "Total Kohn-Sham matrix: \n";
        
        ksmat_str += _xcMatrices[0].getString();
    }
    else
    {
        ksmat_str += "Alpha Kohn-Sham matrix: \n";
        
        ksmat_str += _xcMatrices[0].getString();
        
        ksmat_str += "Beta Kohn-Sham matrix: \n";
        
        ksmat_str += _xcMatrices[1].getString();
    }
    
    return ksmat_str;
}

std::ostream&
operator<<(      std::ostream&     output,
           const CAOKohnShamMatrix& source)
{
    output << std::endl;
    
    output << "[CAOKohnShamMatrix (Object):" << &source << "]" << std::endl;
    
    output << "_xcMatrices: " << std::endl;
    
    for (size_t i = 0; i < source._xcMatrices.size(); i++)
    {
        output << "_xcMatrices[" << i << "]: " << std::endl;
        
        output << source._xcMatrices[i] << std::endl;
    }
    
    output << "_xcRestricted: " << source._xcRestricted << std::endl;
    
    output << "_xcElectrons: " << source._xcElectrons << std::endl;
    
    output << "_xcEnergy" << source._xcEnergy;
    
    return output;
}
