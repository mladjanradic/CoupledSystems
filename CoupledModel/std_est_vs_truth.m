clear, close, clc;

FN = '/Users/mladi/Desktop/GAMMA/Std_Err_Est_part/';

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


for i = 1:15
    data = load([FN 'truth_err' num2str(i) '.mat']);
    m=data.m;
    model_data = data.d;
    d = data.d;
    red_data = data.red_data;
    M = rand_uniform(500, m.model.mu_ranges);
    M(1,:) = 1;
    M(2,:) = 3;
    red_data_OMEGA1 = data.red_data_OMEGA1;
    red_data_OMEGA2 = data.red_data_OMEGA2;
    K1 = blkdiag(data.d.LGG2,0*data.d.LGG2);
    K2 = blkdiag(0*data.d.LGG2,data.d.LGG2);
    for k=1:size(M,2)
        mu = M(:,k);
        m.model = m.model.set_mu(m.model,mu);

        m.operators = @operators;
        red_sim = stabilized_lin_stat_rb_simulation(m,red_data);
        
        uN = red_sim.uN;
        m.operators = @operators_for_res_OMEGA1;
        red_sim_OMEGA1 = stabilized_lin_stat_rb_simulation_OMEGA(m,red_data_OMEGA1,uN);
        m.operators = @operators_for_res_OMEGA2;
        red_sim_OMEGA2 = stabilized_lin_stat_rb_simulation_OMEGA(m,red_data_OMEGA2,uN);
        
        err_OMEGA1(i,k) = red_sim_OMEGA1.res_norm;
        err_OMEGA2(i,k) = red_sim_OMEGA2.res_norm;
        
        resnormf = err_OMEGA1(i,k);
        resnormg = err_OMEGA2(i,k);
        
        DELTAN1(i,k) =  (alpha1./(1-DELTA)) * ( resnormf + ...
            gammab2* ( ( (1/alpha1)*gamma2*resnormf + resnormg )/ ...
            (alpha2*(1/(gamma2^2))*(beta2^2)+alpha1  ) ) );
        
        DELTAN2(i,k) =  ( ((1/alpha1)*gamma2*resnormf + resnormg)/ ...
            (alpha2*(1/(gamma2^2))*beta2^2+alpha1 ) )...
            + (gamma12/(alpha2*(1/(gamma2^2))*beta2^2+alpha2)) * DELTAN1(i);
        
        
        [tru1,tru2,u1,u2] = compute_trace_solutions(m.model,d,d.utilde1,d.utilde2);
        uh = [tru1;tru2];
        U = d.RB*red_sim.uN;
        diff = uh - U;
        
        truth_err_OMEGA1(i,k) = sqrt(abs(diff'*K1*diff));
        truth_err_OMEGA2(i,k) = sqrt(abs(diff'*K2*diff));
        
    end
    disp(i)
    clear data;
    
end


N=size(M,2);
for i=1:15
    DN1(i) = sum(DELTAN1(i,:))/N;
    ERR1(i) = sum(truth_err_OMEGA1(i,:))/N;
    DN2(i) = sum(DELTAN2(i,:))/N;
    ERR2(i) = sum(truth_err_OMEGA2(i,:))/N;
end

figure(1)
semilogy(DN1,'g');
hold on;
semilogy(ERR1,'r');

figure(2)
semilogy(DN2,'g');
hold on;
semilogy(ERR2,'r');


header = {'N','DELTAN1','DELTAN2','ErrOMEGA1','ErrOMEGA2'};
matrix(1:15,1) = 1:15;
matrix(1:15,2) = DN1;
matrix(1:15,3) = DN2;
matrix(1:15,4) = ERR1;
matrix(1:15,5) = ERR2;
matrix2txt(header, matrix, [FN 'std_vs_error.txt'])


