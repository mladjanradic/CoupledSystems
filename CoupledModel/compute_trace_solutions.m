function [tru1,tru2,u1,u2] = compute_trace_solutions(model,model_data,utilde1,utilde2)

mu=model.get_mu(model);
model.decomp_mode = 1;
[A,r] = model.operators(model,model_data);

model.decomp_mode = 2;
[Acoeff,rcoeff] = model.operators(model,model_data);

A1_diff = Acoeff{2}(1)*A{2}{1};
A1_adv = Acoeff{2}(3)*A{2}{3};
A1_dir = Acoeff{2}(5)*A{2}{5};

A2_diff = Acoeff{1}(2)*A{1}{2};
A2_reac = Acoeff{1}(4)*A{1}{4};

r1 = lincomb_sequence(r{2},rcoeff{2});
r2 = lincomb_sequence(r{1},rcoeff{1});

A1 = A1_diff + A1_adv + A1_dir;
A2 = A2_diff + A2_reac;

gamma_dofs1 = [1;model_data.gamma_dofs{2}];
gamma_dofs2 = model_data.gamma_dofs{1};

Lgamma = model_data.gamma_inner_product_matrices{1};
LGG2 = Lgamma(gamma_dofs2,gamma_dofs2);
LGG1 = LGG2;
LGG1(1,:) = 0;
%LGG1(1,1) = 1;

AGI1 = A1(gamma_dofs1,:);
AGI1(:,gamma_dofs1) = [];
AIG1 = A1(:,gamma_dofs1);
AIG1(gamma_dofs1,:) = [];
AGG1 = A1(gamma_dofs1,gamma_dofs1);
rG1 = r1(gamma_dofs1) - AGI1*utilde1;

AGI2 = A2(gamma_dofs2,:);
AGI2(:,gamma_dofs2) = [];
AIG2 = A2(:,gamma_dofs2);
AIG2(gamma_dofs2,:) = [];
AGG2 = A2(gamma_dofs2,gamma_dofs2);
rG2 = r2(gamma_dofs2) - AGI2*utilde2;


AII1 = A1;
AII1(gamma_dofs1,:) = [];
AII1(:,gamma_dofs1) = [];
    
AII2 = A2;
AII2(gamma_dofs2,:) = [];
AII2(:,gamma_dofs2) = [];



AGG1 = AGG1 - AGI1*(AII1\AIG1);
AGG2 = AGG2 - AGI2*(AII2\AIG2);


AGG = [AGG1+mu(3)*LGG1 -mu(3)*LGG1 ; mu(4)*LGG2 AGG2-mu(4)*LGG2];
r = [rG1;rG2];
%uG can be seen as "[uG1;uG2]"
uG = AGG \ r;

tru1 = uG(gamma_dofs1);
uG(gamma_dofs1) = [];
tru2 = uG;



u1 = utilde1 - AII1\(AIG1*tru1);
u2 = utilde2 - AII2\(AIG2*tru2);

end
