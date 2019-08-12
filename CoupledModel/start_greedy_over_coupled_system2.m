function d = start_greedy_over_coupled_system2(model,model_data,m,m1,m2,d)

% uNtilde1 = stabilized_lin_stat_rb_simulation(m1,d.red_data1);
% uNtilde2 = stabilized_lin_stat_rb_simulation(m2,d.red_data2);
% 
% if isnan(uNtilde2.uN)
%     uNtilde2.uN = 0;
% end
% 
% utilde1 = d.RB1*uNtilde1.uN;
% utilde2 = d.RB2*uNtilde2.uN;

W1 = model_data.df_infos{2}.l2_inner_product_matrix + ...
    model_data.df_infos{2}.h10_inner_product_matrix;
W2 = model_data.df_infos{1}.l2_inner_product_matrix + ...
    model_data.df_infos{1}.h10_inner_product_matrix;

mu = d.M_train(:,1);
m.model = m.model.set_mu(m.model,mu);
model = model.set_mu(model,mu);

[utilde1,utilde2] = compute_partial_solutions(model,model_data);
[tru1,tru2,u1,u2] = compute_trace_solutions(model,model_data,utilde1,utilde2);

d.RB1 = [tru1;u1];
d.RB2 = [u2;tru2];

d.RB = [d.RB1 0*d.RB1 ; 0*d.RB2 d.RB2];

red_data = stabilized_lin_stat_gen_reduced_data(m,d);
M_train = d.M_train;

flag=1;
err_seq = [];
ind_seq = [];
while flag
    for i=1:size(M_train,2)
        mu = d.M_train(:,i);
        m.model = m.model.set_mu(m.model,mu);
        model = model.set_mu(model,mu);
        red_sim = stabilized_lin_stat_rb_simulation(m,red_data);
        err(i) = red_sim.res_norm;
    end
    [val,ind] = max(err);
    val=val(1); ind=ind(1);
    err_seq = [err_seq val]
    ind_seq = [ind_seq ind]
    if val<1e-4
        flag=0;
    end
    mu = M_train(:,ind);
    m.model = m.model.set_mu(m.model,mu);
    model = model.set_mu(model,mu);
    [utilde1,utilde2] = compute_partial_solutions(model,model_data);
    [tru1,tru2,u1,u2] = compute_trace_solutions(model,model_data,utilde1,utilde2);
    d.RB1 = orthonormalize_gram_schmidt([d.RB1 [tru1;u1]],W1);
    d.RB2 = orthonormalize_gram_schmidt([d.RB2 [u2;tru2]],W2);
    d.RB = [d.RB1 0*d.RB1 ; 0*d.RB2 d.RB2];
    
    red_data = stabilized_lin_stat_gen_reduced_data(m,d);
    
end
end