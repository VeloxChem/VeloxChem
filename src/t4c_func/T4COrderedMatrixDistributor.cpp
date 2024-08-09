#include "T4COrderedMatrixDistributor.hpp"

#include <set>

#include "StringFormat.hpp"
#include "T4CUtils.hpp"
#include "TensorComponents.hpp"
#include "T4CLocalDistributor.hpp"

CT4COrderedMatrixDistributor::CT4COrderedMatrixDistributor(      CMatrix*       fock,
                                                           const CMatrix*       density,
                                                           const CGtoPairBlock& gto_pair_block,
                                                           const std::string&   label,
                                                           const double         exchange_factor)
{
    _fock = fock;
    
    _density = density;
    
    _label = format::lower_case(label);
    
    _exchange_factor = exchange_factor;
    
    _set_local_matrices(gto_pair_block, gto_pair_block);
}

CT4COrderedMatrixDistributor::CT4COrderedMatrixDistributor(      CMatrix*       fock,
                                                           const CMatrix*       density,
                                                           const CGtoPairBlock& bra_gto_pair_block,
                                                           const CGtoPairBlock& ket_gto_pair_block,
                                                           const std::string&   label,
                                                           const double         exchange_factor)
{
    _fock = fock;
    
    _density = density;
    
    _label = format::lower_case(label);
    
    _exchange_factor = exchange_factor;
    
    _set_local_matrices(bra_gto_pair_block, ket_gto_pair_block);
}

auto
CT4COrderedMatrixDistributor::distribute(const CSimdArray<double>& buffer,
                                         const std::vector<int>&   a_indices,
                                         const std::vector<int>&   b_indices,
                                         const int                 a_angmom,
                                         const int                 b_angmom,
                                         const std::array<int, 2>& bra_range,
                                         const std::array<int, 2>& ket_range) -> void
{
    const auto angpairs = std::array<int, 4>({a_angmom, b_angmom, a_angmom, b_angmom});
    
    const auto tcomps = tensor::number_of_spherical_components(angpairs);
    
    if (_density->get_type() == mat_t::symmetric)
    {
        if (_label == "2jk")
        {
            int offset = 0;
            
            for (int i = bra_range[0]; i < bra_range[1]; i++)
            {
                t4cfunc::local_distribute_rest_jk(_matrices, _density, buffer, offset,
                                                  a_indices, b_indices, a_indices, b_indices,
                                                  _a_loc_indices, _b_loc_indices, _a_loc_indices, _b_loc_indices,
                                                  a_angmom, b_angmom, a_angmom, b_angmom, i, ket_range, true);
                
                offset += tcomps;
            }
        }
    }
}

auto
CT4COrderedMatrixDistributor::distribute(const CSimdArray<double>& buffer,
                                         const std::vector<int>&   a_indices,
                                         const std::vector<int>&   b_indices,
                                         const std::vector<int>&   c_indices,
                                         const std::vector<int>&   d_indices,
                                         const int                 a_angmom,
                                         const int                 b_angmom,
                                         const int                 c_angmom,
                                         const int                 d_angmom,
                                         const std::array<int, 2>& bra_range,
                                         const std::array<int, 2>& ket_range) -> void
{
    const auto angpairs = std::array<int, 4>({a_angmom, b_angmom, c_angmom, d_angmom});
    
    const auto tcomps = tensor::number_of_spherical_components(angpairs);
    
    if (_density->get_type() == mat_t::symmetric)
    {
        if (_label == "2jk")
        {
            int offset = 0;
            
            for (int i = bra_range[0]; i < bra_range[1]; i++)
            {
                t4cfunc::local_distribute_rest_jk(_matrices, _density, buffer, offset,
                                                  a_indices, b_indices, c_indices, d_indices,
                                                  _a_loc_indices, _b_loc_indices, _c_loc_indices, _d_loc_indices,
                                                  a_angmom, b_angmom, c_angmom, d_angmom, i, ket_range, false);
                
                offset += tcomps;
            }
        }
    }
}

auto
CT4COrderedMatrixDistributor::accumulate(const CGtoPairBlock& gto_pair_block) -> void
{
    #pragma omp critical
    {
        const auto a_indices = _get_global_indices(gto_pair_block.bra_orbital_indices());
        
        const auto b_indices = _get_global_indices(gto_pair_block.ket_orbital_indices());
        
        const auto bra_ang_moms = gto_pair_block.angular_momentums();
        
        const auto a_comps = tensor::number_of_spherical_components(bra_ang_moms[0]);
        
        const auto b_comps = tensor::number_of_spherical_components(bra_ang_moms[1]);
        
        if (_density->get_type() == mat_t::symmetric)
        {
            if ((_label == "2jk") || (_label == "2jkx") || (_label == "j"))
            {
                // set up angular pairs
                
                const auto pq_pair = std::array<int, 2>({bra_ang_moms[0], bra_ang_moms[1]});
                
                const auto rs_pair = std::array<int, 2>({bra_ang_moms[0], bra_ang_moms[1]});
                
                // set up Fock submatrices
                
                auto submat_pq = _fock->sub_matrix(pq_pair);
                
                auto submat_rs = _fock->sub_matrix(rs_pair);
                
                // set up angular orders
                
                const auto angord_pq = _fock->is_angular_order(pq_pair);
                
                const auto angord_rs = _fock->is_angular_order(rs_pair);
                
                // acummulate contributions to Fock matrix
                
                t4cfunc::accumulate(submat_pq, _matrices.matrix("PQ")->sub_matrix({0, 0}),
                                    _a_loc_indices, _b_loc_indices, a_indices, b_indices,
                                    a_comps, b_comps, angord_pq);
                
                t4cfunc::accumulate(submat_rs, _matrices.matrix("RS")->sub_matrix({0, 0}),
                                    _a_loc_indices, _b_loc_indices, a_indices, b_indices,
                                    a_comps, b_comps, angord_rs);
            }
            
            if ((_label == "2jk") || (_label == "2jkx") || (_label == "k") || (_label == "kx"))
            {
                // set up angular pairs
                
                const auto pr_pair = std::array<int, 2>({bra_ang_moms[0], bra_ang_moms[0]});
                
                const auto ps_pair = std::array<int, 2>({bra_ang_moms[0], bra_ang_moms[1]});
                
                const auto qr_pair = std::array<int, 2>({bra_ang_moms[1], bra_ang_moms[0]});
                
                const auto qs_pair = std::array<int, 2>({bra_ang_moms[1], bra_ang_moms[1]});
                
                // set up Fock submatrices
                
                auto submat_pr = _fock->sub_matrix(pr_pair);
                
                auto submat_ps = _fock->sub_matrix(ps_pair);
                
                auto submat_qr = _fock->sub_matrix(qr_pair);
                
                auto submat_qs = _fock->sub_matrix(qs_pair);
                
                // set up angular orders
                
                const auto angord_pr = _fock->is_angular_order(pr_pair);
                
                const auto angord_ps = _fock->is_angular_order(ps_pair);
                
                const auto angord_qr = _fock->is_angular_order(qr_pair);
                
                const auto angord_qs = _fock->is_angular_order(qs_pair);
                
                // acummulate contributions to Fock matrix
                
                t4cfunc::accumulate(submat_pr, _matrices.matrix("PR")->sub_matrix({0, 0}),
                                    _a_loc_indices, _a_loc_indices, a_indices, a_indices,
                                    a_comps, a_comps, angord_pr);

                t4cfunc::accumulate(submat_ps, _matrices.matrix("PS")->sub_matrix({0, 0}),
                                    _a_loc_indices, _b_loc_indices, a_indices, b_indices,
                                    a_comps, b_comps, angord_ps);

                t4cfunc::accumulate(submat_qr, _matrices.matrix("QR")->sub_matrix({0, 0}),
                                    _b_loc_indices, _a_loc_indices, b_indices, a_indices,
                                    b_comps, a_comps, angord_qr);

                t4cfunc::accumulate(submat_qs, _matrices.matrix("QS")->sub_matrix({0, 0}),
                                    _b_loc_indices, _b_loc_indices, b_indices, b_indices,
                                    b_comps, b_comps, angord_qs);
            }
        }
    }
}

auto
CT4COrderedMatrixDistributor::accumulate(const CGtoPairBlock& bra_gto_pair_block,
                                         const CGtoPairBlock& ket_gto_pair_block) -> void
{
    #pragma omp critical
    {
        const auto a_indices = _get_global_indices(bra_gto_pair_block.bra_orbital_indices());
        
        const auto b_indices = _get_global_indices(bra_gto_pair_block.ket_orbital_indices());
        
        const auto c_indices = _get_global_indices(ket_gto_pair_block.bra_orbital_indices());
        
        const auto d_indices = _get_global_indices(ket_gto_pair_block.ket_orbital_indices());
        
        const auto bra_ang_moms = bra_gto_pair_block.angular_momentums();
        
        const auto ket_ang_moms = ket_gto_pair_block.angular_momentums();
        
        const auto a_comps = tensor::number_of_spherical_components(bra_ang_moms[0]);
        
        const auto b_comps = tensor::number_of_spherical_components(bra_ang_moms[1]);
        
        const auto c_comps = tensor::number_of_spherical_components(ket_ang_moms[0]);
        
        const auto d_comps = tensor::number_of_spherical_components(ket_ang_moms[1]);
        
        if (_density->get_type() == mat_t::symmetric)
        {
            if ((_label == "2jk") || (_label == "2jkx") || (_label == "j"))
            {
                // set up angular pairs
                
                const auto pq_pair = std::array<int, 2>({bra_ang_moms[0], bra_ang_moms[1]});
                
                const auto rs_pair = std::array<int, 2>({ket_ang_moms[0], ket_ang_moms[1]});
                
                // set up Fock submatrices
                
                auto submat_pq = _fock->sub_matrix(pq_pair);
                
                auto submat_rs = _fock->sub_matrix(rs_pair);
                
                // set up angular orders
                
                const auto angord_pq = _fock->is_angular_order(pq_pair);
                
                const auto angord_rs = _fock->is_angular_order(rs_pair);
                
                // acummulate contributions to Fock matrix
                
                t4cfunc::accumulate(submat_pq, _matrices.matrix("PQ")->sub_matrix({0, 0}),
                                    _a_loc_indices, _b_loc_indices, a_indices, b_indices,
                                    a_comps, b_comps, angord_pq);
                
                t4cfunc::accumulate(submat_rs, _matrices.matrix("RS")->sub_matrix({0, 0}),
                                    _c_loc_indices, _d_loc_indices, c_indices, d_indices,
                                    c_comps, d_comps, angord_rs);
            }
            
            if ((_label == "2jk") || (_label == "2jkx") || (_label == "k") || (_label == "kx"))
            {
                // set up angular pairs
                
                const auto pr_pair = std::array<int, 2>({bra_ang_moms[0], ket_ang_moms[0]});
                
                const auto ps_pair = std::array<int, 2>({bra_ang_moms[0], ket_ang_moms[1]});
                
                const auto qr_pair = std::array<int, 2>({bra_ang_moms[1], ket_ang_moms[0]});
                
                const auto qs_pair = std::array<int, 2>({bra_ang_moms[1], ket_ang_moms[1]});
                
                // set up Fock submatrices
                
                auto submat_pr = _fock->sub_matrix(pr_pair);
                
                auto submat_ps = _fock->sub_matrix(ps_pair);
                
                auto submat_qr = _fock->sub_matrix(qr_pair);
                
                auto submat_qs = _fock->sub_matrix(qs_pair);
                
                // set up angular orders
                
                const auto angord_pr = _fock->is_angular_order(pr_pair);
                
                const auto angord_ps = _fock->is_angular_order(ps_pair);
                
                const auto angord_qr = _fock->is_angular_order(qr_pair);
                
                const auto angord_qs = _fock->is_angular_order(qs_pair);
                
                // acummulate contributions to Fock matrix
                
                t4cfunc::accumulate(submat_pr, _matrices.matrix("PR")->sub_matrix({0, 0}),
                                    _a_loc_indices, _c_loc_indices, a_indices, c_indices,
                                    a_comps, c_comps, angord_pr);
                
                t4cfunc::accumulate(submat_ps, _matrices.matrix("PS")->sub_matrix({0, 0}),
                                    _a_loc_indices, _d_loc_indices, a_indices, d_indices,
                                    a_comps, d_comps, angord_ps);
                
                t4cfunc::accumulate(submat_qr, _matrices.matrix("QR")->sub_matrix({0, 0}),
                                    _b_loc_indices, _c_loc_indices, b_indices, c_indices,
                                    b_comps, c_comps, angord_qr);
                
                t4cfunc::accumulate(submat_qs, _matrices.matrix("QS")->sub_matrix({0, 0}),
                                    _b_loc_indices, _d_loc_indices, b_indices, d_indices,
                                    b_comps, d_comps, angord_qs);
            }
        }
    }
}

auto
CT4COrderedMatrixDistributor::_set_local_matrices(const CGtoPairBlock& bra_gto_pair_block,
                                                  const CGtoPairBlock& ket_gto_pair_block) -> void
{
    // set up local indices
    
    _a_loc_indices = t4cfunc::masked_indices(bra_gto_pair_block.bra_orbital_indices());
    
    _b_loc_indices = t4cfunc::masked_indices(bra_gto_pair_block.ket_orbital_indices());
    
    _c_loc_indices = t4cfunc::masked_indices(ket_gto_pair_block.bra_orbital_indices());
    
    _d_loc_indices = t4cfunc::masked_indices(ket_gto_pair_block.ket_orbital_indices());
    
    // set up local matrices
    
    const auto bra_ang_moms = bra_gto_pair_block.angular_momentums();
    
    const auto ket_ang_moms = ket_gto_pair_block.angular_momentums();
    
    const auto a_dims = _a_loc_indices[0] * tensor::number_of_spherical_components(bra_ang_moms[0]);
    
    const auto b_dims = _b_loc_indices[0] * tensor::number_of_spherical_components(bra_ang_moms[1]);
    
    const auto c_dims = _c_loc_indices[0] * tensor::number_of_spherical_components(ket_ang_moms[0]);
    
    const auto d_dims = _d_loc_indices[0] * tensor::number_of_spherical_components(ket_ang_moms[1]);
    
    if (_density->get_type() == mat_t::symmetric)
    {
        if ((_label == "2jk") || (_label == "2jkx") || (_label == "j"))
        {
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, a_dims, b_dims}), }, mat_t::symmetric), "PQ");
            
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, c_dims, d_dims}), }, mat_t::symmetric), "RS");
        }
        
        if ((_label == "2jk") || (_label == "2jkx") || (_label == "k") || (_label == "kx"))
        {
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, a_dims, c_dims}), }, mat_t::symmetric), "PR");
            
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, a_dims, d_dims}), }, mat_t::symmetric), "PS");
            
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, b_dims, c_dims}), }, mat_t::symmetric), "QR");
            
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, b_dims, d_dims}), }, mat_t::symmetric), "QS");
        }
    }

    if (_density->get_type() == mat_t::general)
    {
        if ((_label == "2jk") || (_label == "2jkx") || (_label == "j"))
        {
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, a_dims, b_dims}), }, mat_t::general), "PQ");
            
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, b_dims, a_dims}), }, mat_t::general), "QP");
            
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, c_dims, d_dims}), }, mat_t::general), "RS");
            
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, d_dims, c_dims}), }, mat_t::general), "SR");
        }
        
        if ((_label == "2jk") || (_label == "2jkx") || (_label == "k") || (_label == "kx"))
        {
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, a_dims, c_dims}), }, mat_t::symmetric), "PR");
            
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, c_dims, a_dims}), }, mat_t::symmetric), "RP");
            
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, a_dims, d_dims}), }, mat_t::symmetric), "PS");
            
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, d_dims, a_dims}), }, mat_t::symmetric), "SP");
            
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, b_dims, c_dims}), }, mat_t::symmetric), "QR");
            
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, c_dims, b_dims}), }, mat_t::symmetric), "RQ");
            
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, b_dims, d_dims}), }, mat_t::symmetric), "QS");
            
            _matrices.add(CMatrix({{0, 0}, }, {CSubMatrix({0, 0, d_dims, b_dims}), }, mat_t::symmetric), "SQ");
        }
    }
    
    _matrices.zero(); 
}

auto
CT4COrderedMatrixDistributor::_get_global_indices(const std::vector<int>& indices) -> std::vector<int>
{
    std::set<int> unique_indices(std::next(indices.cbegin()), indices.cend());
    
    std::vector<int> glob_indices;
    
    glob_indices.push_back(indices[0]);
    
    glob_indices.insert(glob_indices.end(), unique_indices.begin(), unique_indices.end());
    
    return glob_indices;
}
