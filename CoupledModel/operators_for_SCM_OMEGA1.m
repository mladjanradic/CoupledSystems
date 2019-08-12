function [AA,rr] = operators_for_SCM_OMEGA1(model,d)

if model.decomp_mode == 2
    
    model.model.decomp_mode = 2;
    mu = get_mu(model);
    
    AA(1,1) = 1;
    AA(2,1) = mu(3);
    rr=1;
    
end

if model.decomp_mode == 1
    
    S1 = d.AGG1{1} - d.AGI1{1}*(d.AII1{1}\d.AIG1{1});
    S1 = d.AGG1{1} - d.AGI1{1}*(d.AII1{1}\d.AIG1{1});
    
    AA{1} = S1;
    
    Lgamma = d.gamma_inner_product_matrices{1};
    gamma_dofs2 = d.gamma_dofs{1};
    Lgamma2 = Lgamma(gamma_dofs2,gamma_dofs2);
    Lgamma1 = Lgamma2;
    Lgamma1(1,:) = 0;
    
    Lgamma1 = 0.5*(Lgamma1 + Lgamma1');
    
    AA{2} = Lgamma1;
    
    r1 = d.rG1{1} - d.AGI1{1}*d.utilde1;
    rr{1} = r1;
    
    
end




if model.decomp_mode == 0
    mu = get_mu(model);
    S1 = d.AGG1{1} - d.AGI1{1}*(d.AII1{1}\d.AIG1{1});
    
    S1 = 0.5*(S1+S1');
    
    Lgamma = d.gamma_inner_product_matrices{1};
    gamma_dofs2 = d.gamma_dofs{1};
    Lgamma2 = Lgamma(gamma_dofs2,gamma_dofs2);
    Lgamma1 = Lgamma2;
    Lgamma1(1,:) = 0;
    
    Lgamma1 = 0.5*(Lgamma1 + Lgamma1');
    
    AA = S1 + mu(3)*Lgamma1;
    rr = d.rG1{1} - d.AGI1{1}*d.utilde1;
end