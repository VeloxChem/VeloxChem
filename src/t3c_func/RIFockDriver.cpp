#include "RIFockDriver.hpp"

#include "ThreeCenterElectronRepulsionDriver.hpp"
#include "MatrixFunc.hpp"

#include <iostream>

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
CRIFockDriver::prepare_buffers(const CMolecule&        molecule,
                               const CMolecularBasis&  basis,
                               const CMolecularBasis&  aux_basis,
                               const std::vector<int>& atoms) -> void
{
    CThreeCenterElectronRepulsionDriver eri_drv;
    
    _eri_buffer = eri_drv.compute(basis, aux_basis, molecule, atoms);
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
CRIFockDriver::compute(const CMatrix     &density,
                       const std::vector<double>& gvector,
                       const std::string &label) const -> CMatrix
{
    CMatrix fmat(density);
    
    // compute Coulomb contribution to Fock matrix
    
    if ((label == "2jk") || (label == "2jkx") || (label == "j") || (label == "j_rs"))
    {
        fmat.assign_flat_values(_comp_j_vector(gvector));
        
        if ((label == "2jk") || (label == "2jkx"))
        {
            fmat.scale(2.0);
        }
    }
    
    return fmat;
}

auto
CRIFockDriver::local_compute(const CMatrix&             density,
                             const std::vector<double>& gvector,
                             const std::string&         label) const -> CMatrix
{
    CMatrix fmat(density);
    
    // compute Coulomb contribution to Fock matrix
    
    if ((label == "2jk") || (label == "2jkx") || (label == "j") || (label == "j_rs"))
    {
        fmat.assign_flat_values(_comp_local_j_vector(gvector));
        
        if ((label == "2jk") || (label == "2jkx"))
        {
            fmat.scale(2.0);
        }
    }
    
    return fmat;
}

auto
CRIFockDriver::compute(const CMolecule&        molecule,
                       const CMolecularBasis&  basis,
                       const CSubMatrix&       molorbs,
                       const std::string&      label) const -> CMatrix
{
    // set up exchange matrix

    auto exc_mat = matfunc::make_matrix(basis, mat_t::symmetric);
    
    // set up dimensions
    
    const auto nmos = molorbs.number_of_columns();
    
    const auto naos = molorbs.number_of_rows();
    
    const auto ndim = _eri_buffer.aux_width();
    
    const auto nelems = naos * (naos + 1) / 2;
    
    // compute B^Q_it matrices
    
    auto bmats = compute_bq_matrices(molorbs);
    
    // contract B^Q_it matrices to exchange matrix
    
    std::vector<double> exc_values(nelems, 0.0);
    
    for (size_t i = 0; i < ndim; i++)
    {
        auto bmat = bmats.data(i);
        
        for (size_t j = 0; j < naos; j++)
        {
            for (size_t k = j; k < naos; k++)
            {
                double fsum = 0;
                
                for (size_t l = 0; l < nmos; l++)
                {
                    fsum += bmat[l * naos + j] * bmat[l * naos + k];
                }
                
                exc_values[mathfunc::uplo_rm_index(j, k, naos)] = fsum;
            }
        }
    }
    
    exc_mat.assign_flat_values(exc_values);
    
    return exc_mat;
}

auto
CRIFockDriver::compute_bq_vector(const CMatrix &density) const -> std::vector<double>
{
    return _trafo_gamma_vector(_comp_gamma_vector(density));
}

auto
CRIFockDriver::compute_local_bq_vector(const CMatrix &density) const -> std::vector<double>
{
    return _trafo_local_gamma_vector(_comp_gamma_vector(density));
}

auto
CRIFockDriver::compute_bq_matrices(const CSubMatrix& molorbs) const -> CT3RectFlatBuffer<double>
{
    const auto nmos = molorbs.number_of_columns();
    
    const auto naos = molorbs.number_of_rows();
    
    CT3RectFlatBuffer<double> bmats(_eri_buffer.indices(), nmos, naos);
  
    const auto ndim = _eri_buffer.aux_width();
    
    auto buff_ptr = &_eri_buffer;
    
    auto bmats_ptr = &bmats;
    
    auto cmos_ptr = molorbs.data();
    
    for (size_t i = 0; i < ndim; i++)
    {
        auto tint_ptr = buff_ptr->data(i);
        
        auto bmat_ptr = bmats_ptr->data(i);
        
        for (size_t j = 0; j < nmos; j++)
        {
            for (size_t k = 0; k <  naos; k++)
            {
                double fsum = 0.0;
                
                for (size_t l = 0; l < naos; l++)
                {
                    if (l <= k)
                    {
                        fsum += cmos_ptr[l * nmos + j] * tint_ptr[mathfunc::uplo_rm_index(l, k, naos)];
                    }
                    else
                    {
                        fsum += cmos_ptr[l * nmos + j] * tint_ptr[mathfunc::uplo_rm_index(k, l, naos)];
                    }
                }
                
                bmat_ptr[j * naos + k] = fsum;
            }
        }
    }
    
    return bmats;
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
CRIFockDriver::_trafo_local_gamma_vector(const std::vector<double>& gvector) const -> std::vector<double>
{
    const auto ndim = _j_metric->number_of_rows();
    
    std::vector<double> tvec(ndim, 0.0);
    
    auto gvec_ptr = gvector.data();
    
    auto tvec_ptr = tvec.data();
    
    auto tmat_ptr = _j_metric->data();
    
    const auto mask_indices = _eri_buffer.mask_indices();
    
    for (size_t i = 0; i < ndim; i++)
    {
        auto row_ptr = &tmat_ptr[i * ndim];
        
        double fsum = 0.0;
        
        for (const auto [gidx, lidx] : mask_indices)
        {
            fsum += row_ptr[gidx] * gvec_ptr[lidx];
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

auto
CRIFockDriver::_comp_local_j_vector(const std::vector<double>& gvector) const -> std::vector<double>
{
    const auto nrows = _eri_buffer.width();
    
    const auto nelems = nrows * (nrows + 1) / 2;
    
    std::vector<double> jvec(nelems, 0.0);
    
    auto jvec_ptr = jvec.data();
    
    auto gvec_ptr = gvector.data();
    
    auto buff_ptr = &_eri_buffer;
    
    const auto mask_indices = _eri_buffer.mask_indices();
    
    for (const auto [gidx, lidx] : mask_indices)
    {
        auto tint_ptr = buff_ptr->data(lidx);
        
        const auto fact = gvec_ptr[gidx];
         
#pragma omp simd
        for (size_t j = 0; j < nelems; j++)
        {
            jvec_ptr[j] += tint_ptr[j] * fact;
        }
    }
        
    return jvec;
}
