function [ctr1,ctr2] = compute_trace_constant(model_data)

    H10 = model_data.df_infos{1}.l2_inner_product_matrix+...
        model_data.df_infos{1}.h10_inner_product_matrix;
    L2_GAMMA = model_data.gamma_inner_product_matrices{1};
    e = eigs(L2_GAMMA,H10);
    ctr2 = sqrt(max(e));
    
    H10 = model_data.df_infos{2}.l2_inner_product_matrix+...
        model_data.df_infos{2}.h10_inner_product_matrix;
    H10(model_data.df_infos{2}.dirichlet_gids,:) = [];
    H10(:,model_data.df_infos{2}.dirichlet_gids) = [];
    L2 = L2_GAMMA(model_data.gamma_dofs{1},model_data.gamma_dofs{1});
    L2_GAMMA = 0*H10;
    gamma_dofs2 = model_data.gamma_dofs{2};
    L2_GAMMA(gamma_dofs2,gamma_dofs2) = L2(2:end,2:end);
    e = eigs(L2_GAMMA,H10);
    ctr1 = sqrt(max(e));
    
end