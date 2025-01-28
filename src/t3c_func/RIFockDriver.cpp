#include "RIFockDriver.hpp"

#include "ThreeCenterElectronRepulsionDriver.hpp"

CRIFockDriver::CRIFockDriver()

    : _j_metric(nullptr)

    , _k_metric(nullptr)

    , _eri_buffer(CT3FlatBuffer<double>())

    , _eri_trafo_buffer(CT3FlatBuffer<double>())
{
    
}

CRIFockDriver::CRIFockDriver(const CMatrix& j_metric)

    : _j_metric(new CMatrix(j_metric))

    , _k_metric(nullptr)

    , _eri_buffer(CT3FlatBuffer<double>())

    , _eri_trafo_buffer(CT3FlatBuffer<double>())
{
    
}

CRIFockDriver::CRIFockDriver(const CMatrix& j_metric, const CMatrix& k_metric)

    : _j_metric(new CMatrix(j_metric))

    , _k_metric(new CMatrix(k_metric))

    , _eri_buffer(CT3FlatBuffer<double>())

    , _eri_trafo_buffer(CT3FlatBuffer<double>())
{
    
}

CRIFockDriver::~CRIFockDriver()
{
    if (_j_metric != nullptr)
    {
        delete _j_metric;
    }
    
    if (_k_metric != nullptr)
    {
        delete _k_metric;
    }
}

auto
CRIFockDriver::prepare_buffers(const CMolecule       &molecule,
                               const CMolecularBasis &basis,
                               const CMolecularBasis &aux_basis) -> void
{
    CThreeCenterElectronRepulsionDriver eri_drv;
    
    _eri_buffer = eri_drv.compute(basis, aux_basis, molecule);
    
    if (_k_metric != nullptr)
    {
        // FIX ME: Add transformation of integrals required for K fitting.
    }
}

auto
CRIFockDriver::compute(const CMatrix     &density,
                       const std::string &label,
                       const double      exchange_factor) const -> CMatrix
{
    CMatrix fmat(density);
    
    // compute Coulomb contribution to Fock matrix
    
    if ((label == "2jk") || (label == "2jkx") || (label == "j") || (label == "j_rs"))
    {
        
        
    }
    
    // compute exchange contribution to Fock matrix
    
    if ((label == "2jk") || (label == "2jkx") || (label == "k") || (label == "kx") || (label == "k_rs") || (label == "kx_rs"))
    {
        
    }
    
    return fmat;
}

auto
CRIFockDriver::comp_gamma_vector(const CMatrix &density) const -> std::vector<double>
{
    const auto ndim = _eri_buffer.aux_width();
    
    const auto nrows = _eri_buffer.width();
    
    const auto nelems = nrows * (nrows + 1) / 2;
    
    std::vector<double> gvec(ndim, 0.0);
    
    const auto dvec = density.flat_values();
    
    auto dvec_ptr = dvec.data();
    
    auto buff_ptr = &_eri_buffer;
    
    for (size_t i = 0; i < ndim; i++)
    {
        auto tint_ptr = buff_ptr->data(i);
         
        double fsum = 0.0;
        
#pragma omp simd reduction(+:fsum)
        for (size_t j = 0; j < nelems; j++)
        {
            fsum += dvec_ptr[j] * tint_ptr[j];
        }
        
        gvec[i] = fsum;
    }
    
    return gvec;
}
