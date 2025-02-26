#include "RIFockGradDriver.hpp"

#include <ranges>

#include "ThreeCenterElectronRepulsionGeomX00Driver.hpp"
#include "ThreeCenterElectronRepulsionGeom0X0Driver.hpp"
#include "TwoCenterElectronRepulsionGeomX00Driver.hpp"

auto
CRIFockGradDriver::compute(const CMolecularBasis&     basis,
                           const CMolecularBasis&     aux_basis,
                           const CMolecule&           molecule,
                           const std::vector<double>& gamma,
                           const CMatrix&             density,
                           const int                  iatom) const -> TPoint<double>
{
    auto geri = _comp_eri_grad(basis, aux_basis, molecule, gamma, density, iatom);
    
    return TPoint<double>(geri);
}

auto
CRIFockGradDriver::_comp_eri_grad(const CMolecularBasis&     basis,
                                  const CMolecularBasis&     aux_basis,
                                  const CMolecule&           molecule,
                                  const std::vector<double>& gamma,
                                  const CMatrix&             density,
                                  const int                  iatom) const -> std::array<double, 3>
{
    auto g_xyz = std::array<double, 3>({0.0, 0.0, 0.0});
    
   // set up flat density
    
    const auto dvec = density.flat_values();
    
    auto dvec_ptr = dvec.data();
    
    // set up accumulation for aux. basis derivatives contribution
    
    const auto t3cb_drv = CThreeCenterElectronRepulsionGeomX00Driver<1>();
    
    const auto gints = t3cb_drv.compute(basis, aux_basis, molecule, iatom);
    
    auto nelems = dvec.size();

    std::vector<double> gvec(nelems, 0.0);
    
    std::cout << "*** NELEMS " << nelems << std::endl;
    
    auto indices = gints.indices();
    
    std::cout << "*** NINDICES " << indices.size() << std::endl;
    
    for (const auto idx : indices)
    {
        std::cout << idx << " ";
    }
    
    std::cout << std::endl; 
    
    for (size_t i = 0; i < gints.aux_blocks(); i++)
    {
        std::ranges::fill(gvec, 0.0);
        
        std::cout << " INITIAL GVEC for component = " << i << std::endl;
        
        for (const auto g : gvec)
        {
            std::cout << g << std::endl;
        }
        
        auto gvec_ptr = gvec.data();
     
        for (size_t j = 0; j < indices.size(); j++)
        {
            const double fact = gamma[indices[j]];
            
            std::cout << "Gamma(" << j << " -> "<< indices[j] <<") = " << fact << std::endl; 
            
            auto eri_ptr = gints.data(i * indices.size() + j);
        
            #pragma omp simd
            for (size_t k = 0; k < nelems; k++)
            {
                gvec_ptr[k] += fact * eri_ptr[k];
            }
            
            double gsum = 0.0;
            
            for (size_t j = 0; j < nelems; j++)
            {
                gsum += eri_ptr[j];
            }
            
            std::cout << "SUM " << i * indices.size() + j << " = " << gsum << std::endl;
        }
        
        std::cout << "GVEC for component = " << i << std::endl;
        
        for (const auto g : gvec)
        {
            std::cout << g << std::endl;
        }
        
        double gsum = 0.0;
        
        for (size_t j = 0; j < nelems; j++)
        {
            gsum += dvec_ptr[j] * gvec_ptr[j];
        }
        
        g_xyz[i] = 2.0 * gsum;
        
        std::cout << " I was here..." << std::endl;
    }
    
    // set up accumulation of basis derivatives contribution
    
    const auto t3ck_drv = CThreeCenterElectronRepulsionGeom0X0Driver<1>();
    
    const auto tints = t3ck_drv.compute(basis, aux_basis, molecule, iatom);
    
    indices = tints.indices();
    
    auto width = tints.width();
    
    const auto mask_indices = tints.mask_indices();
    
    nelems = mask_indices.size() * width;

    gvec = std::vector<double>(nelems, 0.0);
    
    for (size_t i = 0; i < tints.aux_blocks(); i++)
    {
        std::ranges::fill(gvec, 0.0);
        
        auto gvec_ptr = gvec.data();
        
        for (size_t j = 0; j < indices.size(); j++)
        {
            const double fact = gamma[indices[j]];
            
            auto eri_ptr = tints.data(i * indices.size() + j);
        
            #pragma omp simd
            for (size_t k = 0; k < nelems; k++)
            {
                gvec_ptr[k] += fact * eri_ptr[k];
            }
        }
        
        double gsum = 0.0;
        
        for (const auto [gidx, lidx] : mask_indices)
        {
            for (size_t j = 0; j < width; j++)
            {
                if (gidx < j)
                {
                    gsum += gvec_ptr[lidx * width + j] * dvec_ptr[mathfunc::uplo_rm_index(gidx, j, width)];
                }
                else
                {
                    gsum += gvec_ptr[lidx * width + j] * dvec_ptr[mathfunc::uplo_rm_index(j, gidx, width)];
                }
            }
        }
        
        g_xyz[i] += 2.0 * gsum;
        
        std::cout << " I am here..." << std::endl;
    }
    
    // set up accumulation of basis derivatives contribution
    
    const auto labels = std::vector<std::string>({"X", "Y", "Z"});
    
    const auto t2c_drv = CTwoCenterElectronRepulsionGeomX00Driver<1>();
    
    const auto gmats = t2c_drv.compute(aux_basis, molecule, iatom);
    
    width = indices.size();
    
    for (size_t i = 0; i < 3; i++)
    {
        const auto gmat = gmats.matrix(labels[i])->full_matrix();
        
        double gsum = 0.0;
        
        for (size_t j = 0; j < width; j++)
        {
            for (size_t k = 0; k < width; k++)
            {
                gsum += gamma[i] * gamma[j] * (gmat.at({i, j}) + gmat.at({j, i}));
            }
        }
        
        g_xyz[i] += gsum;
    }
    
    return g_xyz;
}
