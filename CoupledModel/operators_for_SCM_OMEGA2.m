function [AA,rr] = operators_for_SCM_OMEGA2(model,d)

if model.decomp_mode == 2
    model.model.decomp_mode = 2;
    mu = get_mu(model);
    
    AA(1,1)=1;
    AA(2,1) = -mu(4);
    rr=1;
    
end

if model.decomp_mode == 1
    
    S2 = d.AGG2{1} - d.AGI2{1}*(d.AII2{1}\d.AIG2{1});
    
    AA{1} = S2;
    
    Lgamma = d.gamma_inner_product_matrices{1};
    gamma_dofs2 = d.gamma_dofs{1};
    Lgamma2 = Lgamma(gamma_dofs2,gamma_dofs2);
    
    AA{2} = Lgamma2;
    
    r2 = d.rG2{1} - d.AGI2{1}*d.utilde2;
    rr{1} = r2;
    
end



if model.decomp_mode == 0
    mu = get_mu(model);
    S2 = d.AGG2{1} - d.AGI2{1}*(d.AII2{1}\d.AIG2{1});
    Lgamma = d.gamma_inner_product_matrices{1};
    gamma_dofs2 = d.gamma_dofs{1};
    Lgamma2 = Lgamma(gamma_dofs2,gamma_dofs2);
    AA = S2 - mu(4)*Lgamma2;
    rr = d.rG2{1} - d.AGI2{1}*d.utilde2;
end