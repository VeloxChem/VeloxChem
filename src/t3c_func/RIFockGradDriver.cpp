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
    return TPoint<double>(_comp_eri_grad(basis, aux_basis, molecule, gamma, density, iatom));
}

auto
CRIFockGradDriver::compute(const CT4CScreener&        screener,
                           const CMolecularBasis&     basis,
                           const CMolecularBasis&     aux_basis,
                           const CMolecule&           molecule,
                           const std::vector<double>& gamma,
                           const CMatrix&             density,
                           const int                  iatom,
                           const int                  ithreshold) const -> TPoint<double>
{
    return TPoint<double>(_comp_eri_grad(screener, basis, aux_basis, molecule, gamma, density, iatom, ithreshold));
}

auto
CRIFockGradDriver::compute(const CMolecularBasis&     basis,
                           const CMolecularBasis&     aux_basis,
                           const CMolecule&           molecule,
                           const std::vector<double>& gamma,
                           const CMatrix&             density,
                           const std::vector<int>     atoms) const -> std::vector<TPoint<double>>
{
    std::vector<TPoint<double>> grads(atoms.size(), TPoint<double>({0.0, 0.0, 0.0}));
    
    // prepare pointers for OMP parallel region

    auto ptr_gamma = &gamma;

    auto ptr_density = &density;

    auto ptr_molecule = &molecule;
    
    auto ptr_basis = &basis;
    
    auto ptr_aux_basis = &aux_basis;
    
    auto ptr_atoms = atoms.data();
    
    auto ptr_grads = grads.data();

    // execute OMP tasks with static scheduling
    
    const auto natoms = atoms.size();

#pragma omp parallel shared(ptr_gamma, ptr_density, ptr_molecule, ptr_basis, ptr_aux_basis, ptr_grads, ptr_atoms, natoms)
    {
#pragma omp single nowait
        {
            std::ranges::for_each(std::views::iota(size_t{0}, natoms), [&] (const auto index) {
#pragma omp task firstprivate(index)
                {
                    ptr_grads[index] = TPoint<double>(_comp_eri_grad(*ptr_basis, *ptr_aux_basis, *ptr_molecule, *ptr_gamma, *ptr_density, ptr_atoms[index])); 
                }
            });
        }
    }
    
    return grads;
}

auto
CRIFockGradDriver::compute(const CT4CScreener&        screener,
                           const CMolecularBasis&     basis,
                           const CMolecularBasis&     aux_basis,
                           const CMolecule&           molecule,
                           const std::vector<double>& gamma,
                           const CMatrix&             density,
                           const std::vector<int>     atoms,
                           const int                  ithreshold) const -> std::vector<TPoint<double>>
{
    std::vector<TPoint<double>> grads(atoms.size(), TPoint<double>({0.0, 0.0, 0.0}));
    
    // prepare pointers for OMP parallel region

    auto ptr_gamma = &gamma;

    auto ptr_density = &density;

    auto ptr_molecule = &molecule;
    
    auto ptr_basis = &basis;
    
    auto ptr_aux_basis = &aux_basis;
    
    auto ptr_atoms = atoms.data();
    
    auto ptr_grads = grads.data();
    
    auto ptr_screener = &screener;

    // execute OMP tasks with static scheduling
    
    const auto natoms = atoms.size();

#pragma omp parallel shared(ptr_screener, ptr_gamma, ptr_density, ptr_molecule, ptr_basis, ptr_aux_basis, ptr_grads, ptr_atoms, natoms, ithreshold)
    {
#pragma omp single nowait
        {
            std::ranges::for_each(std::views::iota(size_t{0}, natoms), [&] (const auto index) {
#pragma omp task firstprivate(index)
                {
                    ptr_grads[index] = TPoint<double>(_comp_eri_grad(*ptr_screener, *ptr_basis, *ptr_aux_basis, *ptr_molecule, *ptr_gamma, *ptr_density, ptr_atoms[index], ithreshold));
                }
            });
        }
    }
    
    return grads;
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
    
    auto mask_indices = gints.mask_indices();

    for (size_t i = 0; i < gints.aux_blocks(); i++)
    {
        std::ranges::fill(gvec, 0.0);
        
        auto gvec_ptr = gvec.data();
     
        for (const auto [gidx, lidx] : mask_indices)
        {
            const double fact = gamma[gidx];

            auto eri_ptr = gints.data(i * mask_indices.size() + lidx);
        
            #pragma omp simd
            for (size_t k = 0; k < nelems; k++)
            {
                gvec_ptr[k] += fact * eri_ptr[k];
            }
        }
        
        double gsum = 0.0;
        
        for (size_t j = 0; j < nelems; j++)
        {
            gsum += dvec_ptr[j] * gvec_ptr[j];
        }
        
        g_xyz[i] = 4.0 * gsum;
    }
    
    // set up accumulation of basis derivatives contribution
    
    const auto t3ck_drv = CThreeCenterElectronRepulsionGeom0X0Driver<1>();
    
    const auto tints = t3ck_drv.compute(basis, aux_basis, molecule, iatom);
    
    const auto indices = tints.indices();
    
    auto width = tints.width();
    
    mask_indices = tints.mask_indices();
    
    nelems = mask_indices.size() * width;

    gvec = std::vector<double>(nelems, 0.0);
    
    for (size_t i = 0; i < tints.aux_blocks(); i++)
    {
        std::ranges::fill(gvec, 0.0);
        
        auto gvec_ptr = gvec.data();
        
        for (size_t j = 0; j < indices.size(); j++)
        {
            const double fact = gamma[j];
            
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
                else if (gidx == j)
                {
                    gsum += 2.0 * gvec_ptr[lidx * width + j] * dvec_ptr[mathfunc::uplo_rm_index(gidx, j, width)];
                }
                else
                {
                    gsum += gvec_ptr[lidx * width + j] * dvec_ptr[mathfunc::uplo_rm_index(j, gidx, width)];
                }
            }
        }
        
        g_xyz[i] += 4.0 * gsum;
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
                gsum += gamma[j] * gamma[k] * (gmat.at({j, k}) + gmat.at({k, j}));
            }
        }
        
        g_xyz[i] -=  2.0 * gsum;
    }
    
    return g_xyz;
}

auto
CRIFockGradDriver::_comp_eri_grad(const CT4CScreener&        screener,
                                  const CMolecularBasis&     basis,
                                  const CMolecularBasis&     aux_basis,
                                  const CMolecule&           molecule,
                                  const std::vector<double>& gamma,
                                  const CMatrix&             density,
                                  const int                  iatom,
                                  const int                  ithreshold) const -> std::array<double, 3>
{
    auto g_xyz = std::array<double, 3>({0.0, 0.0, 0.0});
    
    // set up flat density
    
    const auto dvec = density.flat_values();
    
    auto dvec_ptr = dvec.data();
    
    // set up accumulation for aux. basis derivatives contribution
    
    const auto t3cb_drv = CThreeCenterElectronRepulsionGeomX00Driver<1>();
    
    const auto gints = t3cb_drv.compute(screener, basis, aux_basis, molecule, iatom, ithreshold);
    
    auto nelems = dvec.size();

    std::vector<double> gvec(nelems, 0.0);
    
    auto mask_indices = gints.mask_indices();

    for (size_t i = 0; i < gints.aux_blocks(); i++)
    {
        std::ranges::fill(gvec, 0.0);
        
        auto gvec_ptr = gvec.data();
     
        for (const auto [gidx, lidx] : mask_indices)
        {
            const double fact = gamma[gidx];

            auto eri_ptr = gints.data(i * mask_indices.size() + lidx);
        
            #pragma omp simd
            for (size_t k = 0; k < nelems; k++)
            {
                gvec_ptr[k] += fact * eri_ptr[k];
            }
        }
        
        double gsum = 0.0;
        
        for (size_t j = 0; j < nelems; j++)
        {
            gsum += dvec_ptr[j] * gvec_ptr[j];
        }
        
        g_xyz[i] = 4.0 * gsum;
    }
    
    // set up accumulation of basis derivatives contribution
    
    const auto t3ck_drv = CThreeCenterElectronRepulsionGeom0X0Driver<1>();
    
    const auto tints = t3ck_drv.compute(basis, aux_basis, molecule, iatom);
    
    const auto indices = tints.indices();
    
    auto width = tints.width();
    
    mask_indices = tints.mask_indices();
    
    nelems = mask_indices.size() * width;

    gvec = std::vector<double>(nelems, 0.0);
    
    for (size_t i = 0; i < tints.aux_blocks(); i++)
    {
        std::ranges::fill(gvec, 0.0);
        
        auto gvec_ptr = gvec.data();
        
        for (size_t j = 0; j < indices.size(); j++)
        {
            const double fact = gamma[j];
            
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
                else if (gidx == j)
                {
                    gsum += 2.0 * gvec_ptr[lidx * width + j] * dvec_ptr[mathfunc::uplo_rm_index(gidx, j, width)];
                }
                else
                {
                    gsum += gvec_ptr[lidx * width + j] * dvec_ptr[mathfunc::uplo_rm_index(j, gidx, width)];
                }
            }
        }
        
        g_xyz[i] += 4.0 * gsum;
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
                gsum += gamma[j] * gamma[k] * (gmat.at({j, k}) + gmat.at({k, j}));
            }
        }
        
        g_xyz[i] -=  2.0 * gsum;
    }
    
    return g_xyz;
}
