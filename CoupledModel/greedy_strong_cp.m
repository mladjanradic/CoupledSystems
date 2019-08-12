function d = greedy_strong_cp(m,model,model_data,d)


d.RB = m.detailed_simulation(m,model_data);

red_data = stabilized_lin_stat_gen_reduced_data(m,d);
M_train = d.M_train;

flag=1;
err_seq = [];
ind_seq = [];

K=d.W;

j=0;

while flag
    for i=1:size(M_train,2)
        mu = d.M_train(:,i);
        m.model = m.model.set_mu(m.model,mu);
        model = model.set_mu(model,mu);
        red_sim = stabilized_lin_stat_rb_simulation(m,red_data);
        load(['/Users/mladi/Desktop/CP_sol/sol' num2str(i) '.mat'])
        diff = uh - d.RB*red_sim.uN;
        err(i) = sqrt(diff'*K*diff);
    end
    j=j+1;
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
    uh = m.detailed_simulation(m,model_data);
    d.RB = orthonormalize_qr([d.RB uh],d.W);
    %d.RB = [d.RB data.uh];
    red_data = stabilized_lin_stat_gen_reduced_data(m,d);
    
end

d.RB_info.err_seq = err_seq;
d.RB_info.ind_seq = ind_seq;

end