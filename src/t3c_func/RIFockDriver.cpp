#include "RIFockDriver.hpp"

#include "ThreeCenterElectronRepulsionDriver.hpp"

CRIFockDriver::CRIFockDriver()

    : _j_metric(nullptr)

    , _eri_buffer(CT3FlatBuffer<double>())
{
    
}

CRIFockDriver::CRIFockDriver(const CSubMatrix& j_metric)

    : _j_metric(new CSubMatrix(j_metric))

    , _eri_buffer(CT3FlatBuffer<double>())
{
    
}

CRIFockDriver::~CRIFockDriver()
{
    if (_j_metric != nullptr)
    {
        delete _j_metric;
    }
}

auto
CRIFockDriver::prepare_buffers(const CMolecule       &molecule,
                               const CMolecularBasis &basis,
                               const CMolecularBasis &aux_basis) -> void
{
    CThreeCenterElectronRepulsionDriver eri_drv;
    
    _eri_buffer = eri_drv.compute(basis, aux_basis, molecule);
}

auto
CRIFockDriver::compute(const CMatrix     &density,
                       const std::string &label) const -> CMatrix
{
    CMatrix fmat(density);
    
    // compute Coulomb contribution to Fock matrix
    
    if ((label == "2jk") || (label == "2jkx") || (label == "j") || (label == "j_rs"))
    {
        auto t_vec = _trafo_gamma_vector(_comp_gamma_vector(density));
        
        fmat.assign_flat_values(_comp_j_vector(t_vec));
        
        if ((label == "2jk") || (label == "2jkx"))
        {
            fmat.scale(2.0); 
        }
    }
    
    return fmat;
}

auto
CRIFockDriver::compute_bq_vector(const CMatrix &density) const -> std::vector<double>
{
    return _trafo_gamma_vector(_comp_gamma_vector(density));
}

auto
CRIFockDriver::_comp_gamma_vector(const CMatrix &density) const -> std::vector<double>
{
    const auto ndim = _eri_buffer.aux_width();
    
    const auto nrows = _eri_buffer.width();
    
    const auto nelems = nrows * (nrows + 1) / 2;
    
    std::vector<double> gvec(ndim, 0.0);
    
    const auto dvec = density.flat_values();
    
    auto dvec_ptr = dvec.data();
    
    auto buff_ptr = &_eri_buffer;
    
    auto gvec_ptr = gvec.data();
    
    for (size_t i = 0; i < ndim; i++)
    {
        auto tint_ptr = buff_ptr->data(i);
         
        double fsum = 0.0;
        
#pragma omp simd reduction(+:fsum)
        for (size_t j = 0; j < nelems; j++)
        {
            fsum += dvec_ptr[j] * tint_ptr[j];
        }
        
        gvec_ptr[i] = fsum;
    }
    
    return gvec;
}

auto
CRIFockDriver::_trafo_gamma_vector(const std::vector<double>& gvector) const -> std::vector<double>
{
    const auto ndim = _eri_buffer.aux_width();
    
    std::vector<double> tvec(ndim, 0.0);
    
    auto gvec_ptr = gvector.data();
    
    auto tvec_ptr = tvec.data();
    
    auto tmat_ptr = _j_metric->data();
    
    for (size_t i = 0; i < ndim; i++)
    {
        auto row_ptr = &tmat_ptr[i * ndim];
        
        double fsum = 0.0;
        
#pragma omp simd reduction(+:fsum)
        for (size_t j = 0; j < ndim; j++)
        {
            fsum += row_ptr[j] * gvec_ptr[j];
        }
        
        tvec_ptr[i] = fsum;
    }
    
    return tvec;
}

auto
CRIFockDriver::_comp_j_vector(const std::vector<double>& gvector) const -> std::vector<double>
{
    const auto ndim = _eri_buffer.aux_width();
    
    const auto nrows = _eri_buffer.width();
    
    const auto nelems = nrows * (nrows + 1) / 2;
    
    std::vector<double> jvec(nelems, 0.0);
    
    auto jvec_ptr = jvec.data();
    
    auto gvec_ptr = gvector.data();
    
    auto buff_ptr = &_eri_buffer;
    
    for (size_t i = 0; i < ndim; i++)
    {
        auto tint_ptr = buff_ptr->data(i);
        
        const auto fact = gvec_ptr[i];
         
#pragma omp simd
        for (size_t j = 0; j < nelems; j++)
        {
            jvec_ptr[j] += tint_ptr[j] * fact;
        }
    }
        
    return jvec;
}

