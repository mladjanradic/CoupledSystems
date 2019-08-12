function d = start_std_greedy_over_GAMMA(m,d)

FN = '/Users/mladi/Desktop/GAMMA/';
FN_save = [FN 'Std_Err_Est_part/'];

K = blkdiag(d.LGG2,d.LGG2);
K1 = blkdiag(d.LGG2,0*d.LGG2);
K2 = blkdiag(0*d.LGG2,d.LGG2);
d.W = K;

data = load([FN 'PartSol/sol' num2str(1) '.mat']);
d.RB = [data.tru1;data.tru2];
d.RB_glob = [data.tru1;data.tru2];

m.operators = @operators;
red_data = stabilized_lin_stat_gen_reduced_data(m,d);
m.operators = @operators_for_res_OMEGA1;
red_data_OMEGA1 = stabilized_lin_stat_gen_reduced_data(m,d);
m.operators = @operators_for_res_OMEGA2;
red_data_OMEGA2 = stabilized_lin_stat_gen_reduced_data(m,d);

M_train = d.M_train;
val_OMEGA1 = 1;
val_OMEGA2 = 1;

flag=1;
err_seq = [];
ind_seq = [];
mu_seq = [];


err_seq_OMEGA1 = [];
ind_seq_OMEGA1 = [];
mu_seq_OMEGA1 = [];

err_seq_OMEGA2 = [];
ind_seq_OMEGA2 = [];
mu_seq_OMEGA2 = [];


mkdir(FN_save);
j=0;
N=1;



DELTA = 0.8;
alpha1 = 0.5;
alpha2 = 0.5;
gamma1 = 1;
gammab1 = 1;
gamma2 = 1;
gammab2 = 1;
beta1 = 0.5;
beta2 = 0.5;
gamma12 = 1;


while flag
    for i=1:size(M_train,2)
        mu = M_train(:,i);
        m.model = m.model.set_mu(m.model,mu);
        m.operators = @operators;
        red_sim = stabilized_lin_stat_rb_simulation(m,red_data);
        err(i) = red_sim.res_norm;
        uN = red_sim.uN;
        
        m.operators = @operators_for_res_OMEGA1;
        red_sim_OMEGA1 = stabilized_lin_stat_rb_simulation_OMEGA(m,red_data_OMEGA1,uN);
        m.operators = @operators_for_res_OMEGA2;
        red_sim_OMEGA2 = stabilized_lin_stat_rb_simulation_OMEGA(m,red_data_OMEGA2,uN);
        
        err_OMEGA1(i) = red_sim_OMEGA1.res_norm;
        err_OMEGA2(i) = red_sim_OMEGA2.res_norm;
        
        resnormf = err_OMEGA1(i);
        resnormg = err_OMEGA2(i);
        
        DELTAN1(i) =  (alpha1./(1-DELTA)) * ( resnormf + ...
            gammab2* ( ( (1/alpha1)*gamma2*resnormf + resnormg )/ ...
            (alpha2*(1/(gamma2^2))*(beta2^2)+alpha1  ) ) );
        
        DELTAN2(i) =  ( ((1/alpha1)*gamma2*resnormf + resnormg)/ ...
            (alpha2*(1/(gamma2^2))*beta2^2+alpha1 ) )...
            + (gamma12/(alpha2*(1/(gamma2^2))*beta2^2+alpha2)) * DELTAN1(i);
        
        load([FN 'PartSol/sol' num2str(i) '.mat'])
        uh = [tru1;tru2];
        U = d.RB*red_sim.uN;
        diff = uh - U;
        truth_err(i) = sqrt(abs(diff'*K*diff));
        
        truth_err_OMEGA1(i) = sqrt(abs(diff'*K1*diff));
        truth_err_OMEGA2(i) = sqrt(abs(diff'*K2*diff));
        
    end
    
    j=j+1;
    save([FN_save 'err' num2str(j)]);
    save([FN_save 'truth_err' num2str(j)]);
    
    
    [val,ind] = max(err);
    val=val(1); ind=ind(1);
    err_seq = [err_seq val]
    ind_seq = [ind_seq ind]
    mu = M_train(:,ind);
    mu_seq = [mu_seq mu];
    
    data = load([FN 'PartSol/sol' num2str(ind) '.mat']);
    u = [data.tru1;data.tru2];
    d.RB_glob = orthonormalize_qr([d.RB_glob u],K);
    
    
    
    %%Omega1:
    [val_OMEGA1,ind] = max(DELTAN1);
    val_OMEGA1=val_OMEGA1(1); ind=ind(1);
    err_seq_OMEGA1 = [err_seq_OMEGA1 val_OMEGA1]
    ind_seq_OMEGA1 = [ind_seq_OMEGA1 ind]
    mu_OMEGA1 = M_train(:,ind);
    mu_seq_OMEGA1 = [mu_seq_OMEGA1 mu_OMEGA1];
    
    data = load([FN 'PartSol/sol' num2str(ind) '.mat']);
    u1 = data.tru1;
    
    
    
    %%Omega2:
    [val_OMEGA2,ind] = max(DELTAN2);
    val_OMEGA2=val_OMEGA2(1); ind=ind(1);
    err_seq_OMEGA2 = [err_seq_OMEGA2 val_OMEGA2]
    ind_seq_OMEGA2 = [ind_seq_OMEGA2 ind]
    mu_OMEGA2 = M_train(:,ind);
    mu_seq_OMEGA2 = [mu_seq_OMEGA2 mu_OMEGA2];
    
    data = load([FN 'PartSol/sol' num2str(ind) '.mat']);
    u2 = data.tru2;
    
    d.RB = orthonormalize_qr([d.RB,[u1;u2]],K);
    
    
    m.operators = @operators;
    red_data = stabilized_lin_stat_gen_reduced_data(m,d);
    m.operators = @operators_for_res_OMEGA1;
    red_data_OMEGA1 = stabilized_lin_stat_gen_reduced_data(m,d);
    m.operators = @operators_for_res_OMEGA2;
    red_data_OMEGA2 = stabilized_lin_stat_gen_reduced_data(m,d);
    
    
    
    
    
    
    if N == 20
        flag=0;
    end
    N = N+1;
end

d.RB_info.err_seq = err_seq;
d.RB_info.ind_seq = ind_seq;

end