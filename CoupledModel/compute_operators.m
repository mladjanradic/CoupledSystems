function [AA,r] = compute_operators(model,model_data)
    
    mu = get_mu(model);
    model.decomp_mode = 1;
    [A,r] = model.operators(model,model_data);

    model.decomp_mode = 2;
    [Acoeff,rcoeff] = model.operators(model,model_data);
    
    A1_diff = Acoeff{2}(1)*A{2}{1};
    A1_adv  = Acoeff{2}(3)*A{2}{3};
    A1_dir  = Acoeff{2}(5)*A{2}{5};
    
    A2_diff = Acoeff{1}(2)*A{1}{2};
    A2_reac = Acoeff{1}(4)*A{1}{4};
    
    r1 = lincomb_sequence(r{2},rcoeff{2});
    r2 = lincomb_sequence(r{1},rcoeff{1});
    
    A1 = A1_diff + A1_adv + A1_dir;
    A2 = A2_diff + A2_reac;
    
    
    Lgamma = model_data.gamma_inner_product_matrices{1};
   
    
    gamma_dofs1 = [1;model_data.gamma_dofs{2}];
    gamma_dofs2 = model_data.gamma_dofs{1};
    Lgamma2 = Lgamma(gamma_dofs2,gamma_dofs2);
    A1(gamma_dofs1(2:end),gamma_dofs1) = A1(gamma_dofs1(2:end),gamma_dofs1) + mu(3)*Lgamma2(2:end,:);
    A2(gamma_dofs2,gamma_dofs2) = A2(gamma_dofs2,gamma_dofs2) - mu(4)*Lgamma(gamma_dofs2,gamma_dofs2);
    
    global_index = size(A1,1);
    gamma_glob2 = global_index + gamma_dofs2;
    
    N = size(A1,1) + size(A2,1);
    
    AA=sparse(1,1,0,N,N);
    AA(1:size(A1,1),1:size(A1,1)) = A1;
    AA(size(A1,1)+1:end,size(A1,1)+1:end) = A2;
    
    AA(gamma_dofs1(2:end),gamma_glob2) = - mu(3)*Lgamma2(gamma_dofs1(2:end),:);
    AA(gamma_glob2,gamma_dofs1) = + mu(4)*Lgamma2(:,gamma_dofs1);
    
    r = [r1;r2];
    
end