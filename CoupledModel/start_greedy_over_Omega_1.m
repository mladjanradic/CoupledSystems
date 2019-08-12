function [detailed_data,m1] = start_greedy_over_Omega_1(model,detailed_data)

M_train = detailed_data.M_train;

m1.get_inner_product_matrix = @(d) d.W1;
m1.get_rb_size = @(m1,d) size(d.RB,2);
m1 = gen_m1(m1,detailed_data);
m1.varepsilon = model.base_model.varepsilon;
m1.mus = M_train(:,1);
uh = m1.detailed_simulation(m1,detailed_data);
detailed_data.RB = uh;

red_data1 = stabilized_lin_stat_gen_reduced_data(m1,detailed_data);

tol1 = 1e-4;
flag=1;
while flag
    for i=1:size(M_train,2)
        m1.mus = M_train(:,i);
        sim_red = stabilized_lin_stat_rb_simulation(m1,red_data1);
        err(i) = sim_red.res_norm;
    end
    [val,ind] = max(err);
    if max(err) < tol1
        flag=0;
    else
        u = m1.detailed_simulation(m1,detailed_data);
        detailed_data.RB = orthonormalize_gram_schmidt([detailed_data.RB,u],detailed_data.W1);
        red_data1 = stabilized_lin_stat_gen_reduced_data(m1,detailed_data);
    end
end

detailed_data.RB1 = detailed_data.RB;
detailed_data.red_data1 = red_data1;