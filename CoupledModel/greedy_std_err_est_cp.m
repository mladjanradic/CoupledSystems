function d = greedy_std_err_est_cp(m,model,model_data,d)

uh = m.detailed_simulation(m,model_data);
d.RB = uh;

red_data = stabilized_lin_stat_gen_reduced_data(m,d);
M_train = d.M_train;

flag=1;
err_seq = [];
ind_seq = [];

FN = '/Users/mladi/Desktop/HighDim/Std_Error';
mkdir(FN);

K = d.W;

flag_comp_sol = 0;
%load('beta.mat');
%BETA = BETA(end:-1:1);
BETA = 0.001;
j=0;

while flag
    for i=1:size(M_train,2)
        mu = d.M_train(:,i);
        m.model = m.model.set_mu(m.model,mu);
        model = model.set_mu(model,mu);
        
        if flag_comp_sol
            uh = m.detailed_simulation(m,model_data);
            save(['/Users/mladi/Desktop/Sol/sol' num2str(i)],'uh');
        else
            %load(['/Users/mladi/Desktop/Sol/sol' num2str(i) '.mat']);
            load(['/Users/mladi/Desktop/HighDim/Sol/sol' num2str(i) '.mat']);
        end
        red_sim = stabilized_lin_stat_rb_simulation(m,red_data);
        err(i) = (1/BETA)*red_sim.res_norm;
        
        UN = d.RB*red_sim.uN;
        diff = uh - UN;
        diff(model_data.df_infos{2}.dirichlet_gids) = 0;
        truth_err(i) = sqrt(diff'*K*diff);
    end
    j=j+1;
    save(['/Users/mladi/Desktop/HighDim/Std_Error/err' num2str(j)],'err');
    save(['/Users/mladi/Desktop/HighDim/Std_Error/truth_err' num2str(j)],'truth_err');
    
    flag_comp_sol = 0;
    
    [val,ind] = max(err);
    val=val(1); ind=ind(1);
    err_seq = [err_seq val]
    ind_seq = [ind_seq ind];
    if val<1e-4
        flag=0;
    end
    mu = M_train(:,ind);
    m.model = m.model.set_mu(m.model,mu);
    uh = m.detailed_simulation(m,model_data);
    d.RB = orthonormalize_gram_schmidt([d.RB uh],d.W);
    red_data = stabilized_lin_stat_gen_reduced_data(m,d);
    
end

d.RB_info.err_seq = err_seq;
d.RB_info.ind_seq = ind_seq;

end