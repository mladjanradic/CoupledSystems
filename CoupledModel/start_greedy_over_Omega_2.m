function [detailed_data,m2] = start_greedy_over_Omega_2(model,detailed_data)

M_train = detailed_data.M_train;

m2.get_inner_product_matrix = @(d) d.W2;
m2.get_rb_size = @(m2,d) size(d.RB,2);
m2 = gen_m2(m2,detailed_data);
m2.mus = M_train(:,1);
uh = m2.detailed_simulation(m2,detailed_data);
detailed_data.RB = uh;

red_data2 = stabilized_lin_stat_gen_reduced_data(m2,detailed_data);

tol1 = 1e-4;
flag=1;
while flag
    for i=1:size(M_train,2)
        m2.mus = M_train(:,i);
        sim_red = stabilized_lin_stat_rb_simulation(m2,red_data2);
        err(i) = sim_red.res_norm;
    end
    [val,ind] = max(err);
    if max(err) < tol1
        flag=0;
    else
        u = m2.detailed_simulation(m2,detailed_data);
        detailed_data.RB = orthonormalize_gram_schmidt([detailed_data.RB,u],detailed_data.W2);
        red_data2 = stabilized_lin_stat_gen_reduced_data(m2,detailed_data);
    end
end

detailed_data.RB2 = detailed_data.RB;
detailed_data.red_data2 = red_data2;