function [u1 , u2] = compute_partial_solutions(model,model_data)

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
    
    A1(gamma_dofs1,:) = [];
    A1(:,gamma_dofs1) = [];
    r1(gamma_dofs1) = [];
    u1 = A1 \ r1;
    
    A2(gamma_dofs2,:) = [];
    A2(:,gamma_dofs2) = [];
    r2(gamma_dofs2) = [];
    u2 = A2 \ r2;
    
end