function [AA,rr] = operators_for_res_OMEGA2(model,d)

if model.decomp_mode == 2
    model.model.decomp_mode = 2;
    mu = get_mu(model.model);
    
    AA(1,1)=1;
    AA(2,1) = 0;
    AA(3,1) = 0;
    AA(4,1) = mu(4);
    AA(5,1) = -mu(4);
    rr=[0;1];
    
end

if model.decomp_mode == 1
    
S1 = d.AGG1{1} - d.AGI1{1}*(d.AII1{1}\d.AIG1{1});
S2 = d.AGG2{1} - d.AGI2{1}*(d.AII2{1}\d.AIG2{1});

AA{1} = blkdiag(0*S1,S2);

Lgamma = d.gamma_inner_product_matrices{1};
gamma_dofs2 = d.gamma_dofs{1};
Lgamma2 = Lgamma(gamma_dofs2,gamma_dofs2);
Lgamma1 = Lgamma2;
Lgamma1(1,:) = 0;

AA{2} = sparse([Lgamma1 0*Lgamma1 ; 0*Lgamma2 0*Lgamma2]);
AA{3} = sparse([0*Lgamma1 Lgamma1 ; 0*Lgamma2 0*Lgamma2]);
AA{4} = sparse([0*Lgamma1 0*Lgamma1 ; Lgamma2 0*Lgamma2]);
AA{5} = sparse([0*Lgamma1 0*Lgamma1 ; 0*Lgamma2 Lgamma2]);


r1 = d.rG1{1} - d.AGI1{1}*d.utilde1;
r2 = d.rG2{1} - d.AGI2{1}*d.utilde2;

rr{1} = [r1;0*r2];
rr{2} = [0*r1;r2];

end