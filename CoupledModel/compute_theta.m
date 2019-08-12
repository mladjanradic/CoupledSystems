clear, close, clc;
%{
FN = '/Users/mladi/Desktop/HighDim/Std_Err_Est/';

i=1;
data = load([FN 'truth_err' num2str(i) '.mat']);
m=data.m;
model=data.model;
model_data = data.model_data;
M_train = data.M_train;
K = data.K;
clear data;

flag_dir=0;

for D = 2:5
    
    theta_part = zeros(2,15-D);
    theta = zeros(1,15-D);
    
    FN_save = ['/Users/mladi/Desktop/HighDim/PartTHETA_onGamma_D' num2str(D) '/'];
    mkdir(FN_save);
    
    for i=1:15-D
        dataN = load([FN 'truth_err' num2str(i) '.mat']);
        red_dataN = dataN.red_data;
        
        dataM = load([FN 'truth_err' num2str(i+D) '.mat']);
        red_dataM = dataM.red_data;
        N=i;
        M=i+D;
        if flag_dir
            for k = 1:size(M_train,2)
                mu = M_train(:,k);
                
                model = model.set_mu(model,mu);
                m.model = m.model.set_mu(m.model,mu);
                
                red_simN = stabilized_lin_stat_rb_simulation(m,red_dataN);
                red_simM = stabilized_lin_stat_rb_simulation(m,red_dataM);
                load(['/Users/mladi/Desktop/HighDim/PartSol/sol' num2str(k) '.mat'])
                UN1 = dataN.d.RB1*red_simN.uN(1:N);
                UN2 = dataN.d.RB2*red_simN.uN(N+1:end);
                UM1 = dataM.d.RB1*red_simM.uN(1:M);
                UM2 = dataM.d.RB2*red_simM.uN(M+1:end);
                
                UN1(model_data.df_infos{2}.dirichlet_gids) = model.base_model.uD;
                UM1(model_data.df_infos{2}.dirichlet_gids) = model.base_model.uD;
                
                G1 = length(tru1);
                G2 = length(tru2);
                
                diff = tru1 - UM1(1:G1); diff(1) = 0;
                zaehler(1,k) = sqrt(abs(diff'*dataN.d.LGG1*diff));
                diff = tru1 - UN1(1:G1); diff(1) = 0;
                nenner(1,k) = sqrt(abs(diff'*dataN.d.LGG1*diff));
                
                diff = tru2 - UM2(end-G2+1:end);
                zaehler(2,k) = sqrt(abs(diff'*dataN.d.LGG2*diff));
                diff = tru2 - UN2(end-G2+1:end);
                nenner(2,k) = sqrt(abs(diff'*dataN.d.LGG2*diff));
            end
        else
            FN_LOAD='/Users/mladi/Desktop/HighDim/PartTHETA_onGamma_D1/';
            DATA_N = load([FN_LOAD 'thetaNM' num2str(N) num2str(N+1)]);
            DATA_M = load([FN_LOAD 'thetaNM' num2str(M-1) num2str(M)]);
            zaehler(1,:) = DATA_M.zaehler(1,:);
            nenner(1,:) = DATA_N.nenner(1,:);
            
            zaehler(2,:) = DATA_M.zaehler(2,:);
            nenner(2,:) = DATA_N.nenner(2,:);
        end
        
        ind1 = find(nenner(1,:)>(norm(nenner(1,:))/100));
        %ind1 = find(nenner(1,:)>1e-6);
        MAX1 = max(zaehler(1,ind1)./nenner(1,ind1));
        theta_part(1,i) = MAX1(1);
        
        ind2 = find(nenner(2,:)>(norm(nenner(2,:))/100));
        %ind2 = find(nenner(2,:)>1e-6);
        MAX2 = max(zaehler(2,ind2)./nenner(2,ind2));
        theta_part(2,i) = MAX2(1);
        
%{
    nenner_glob = dataN.truth_err;
    zaehler_glob = dataM.truth_err;
    %ind = find(nenner_glob(1,:)>1e-6);
    ind = find(nenner_glob(1,:)>norm(nenner_glob(1,:))/100);
    MAX = max(zaehler_glob(1,ind)./nenner_glob(1,ind));
    theta(1,i) = MAX(1);
%}
        save([FN_save 'thetaNM' num2str(N) num2str(M)]);
        
        disp(i)
        
    end
    theta_part
    flag_dir = 0;
end
%}

%% Now control it on a test set:
%{

D=3;
FN_load = ['/Users/mladi/Desktop/HighDim/PartTHETA_D' num2str(D) '/'];

FN_save = ['/Users/mladi/Desktop/HighDim/HierErr_D' num2str(D) '/'];

mkdir(FN_save);
data = load([FN_load 'thetaNM1' num2str(1+D) '.mat']);
m=data.m;
model=data.model;
model_data = data.model_data;
Mtrain = rand_uniform(1000, model.mu_ranges);
FN = '/Users/mladi/Desktop/HighDim/Std_Err_Est/';

for i=1:15-D
    
    load([FN_load 'thetaNM' num2str(i) num2str(i+D) '.mat']);
    
    dataN = load([FN 'truth_err' num2str(i) '.mat']);
    red_dataN = dataN.red_data;
    
    dataM = load([FN 'truth_err' num2str(i+D) '.mat']);
    red_dataM = dataM.red_data;
    N=i;
    M=i+D;
    
    for k = 1:size(Mtrain,2)
        mu = Mtrain(:,k);
        m.model = m.model.set_mu(m.model,mu);
        model = model.set_mu(model,mu);
        
        red_sim_N = stabilized_lin_stat_rb_simulation(m,red_dataN);
        red_sim_M = stabilized_lin_stat_rb_simulation(m,red_dataM);

        UN = dataN.d.RB*red_sim_N.uN;
        UM = dataM.d.RB*red_sim_M.uN;
        
        
        UN1 = dataN.d.RB1*red_sim_N.uN(1:N);
        UN2 = dataN.d.RB2*red_sim_N.uN(N+1:end);
        UM1 = dataM.d.RB1*red_sim_M.uN(1:M);
        UM2 = dataM.d.RB2*red_sim_M.uN(M+1:end);
        
        
        diff1 = UN1 - UM1;
        diff1(model_data.df_infos{2}.dirichlet_gids) = 0;
        DELTA_part1(k) = sqrt(abs(diff1'*dataN.M1*diff1));
        
        diff2 = UN2 - UM2;
        diff2(model_data.df_infos{2}.dirichlet_gids) = 0;
        DELTA_part2(k) = sqrt(abs(diff2'*dataN.M2*diff2));
        
        
        diff = UN - UM;
        diff(model_data.df_infos{2}.dirichlet_gids) = 0;
        DELTA(k) = sqrt(abs(diff'*K*diff));
        
        
        [utilde1,utilde2] = compute_partial_solutions(model,model_data);
        [tru1,tru2,u1,u2] = compute_trace_solutions(model,model_data,utilde1,utilde2);
        
        diff_part1 = [tru1;u1] - UN1;
        diff_part2 = [u2;tru2] - UN2;
        diff_glob = [tru1;u1;u2;tru2] - UN;
        
        
        truth_err1(k) = sqrt(diff_part1'*dataN.M1*diff_part1);
        truth_err2(k) = sqrt(diff_part2'*dataN.M2*diff_part2);
        truth_err(k) = sqrt(diff_glob'*dataN.K*diff_glob);
        
    end
    
    hier_err_est = DELTA./(1-theta(i));
    hier_err_est1 = DELTA_part1./(1-theta_part(1,i));
    hier_err_est2 = DELTA_part2./(1-theta_part(2,i));
    FN_save = ['/Users/mladi/Desktop/HighDim/HierErr_D' num2str(D) '/'];
    save([FN_save 'data' num2str(i) '.mat']);
    disp(i)
end
%}




%% Now control it on a test set:

D=3;
FN_load = ['/Users/mladi/Desktop/HighDim/PartTHETA_onGamma_D' num2str(D) '/'];

FN_save = ['/Users/mladi/Desktop/HighDim/HierErr_onGamma_D' num2str(D) '/'];

mkdir(FN_save);
data = load([FN_load 'thetaNM1' num2str(1+D) '.mat']);
m=data.m;
model=data.model;
model_data = data.model_data;
Mtrain = rand_uniform(1000, model.mu_ranges);
FN = '/Users/mladi/Desktop/HighDim/Std_Err_Est/';

for i=1:15-D
    
    load([FN_load 'thetaNM' num2str(i) num2str(i+D) '.mat']);
    
    dataN = load([FN 'truth_err' num2str(i) '.mat']);
    red_dataN = dataN.red_data;
    
    dataM = load([FN 'truth_err' num2str(i+D) '.mat']);
    red_dataM = dataM.red_data;
    N=i;
    M=i+D;
    
    for k = 1:size(Mtrain,2)
        mu = Mtrain(:,k);
        m.model = m.model.set_mu(m.model,mu);
        model = model.set_mu(model,mu);
        
        red_sim_N = stabilized_lin_stat_rb_simulation(m,red_dataN);
        red_sim_M = stabilized_lin_stat_rb_simulation(m,red_dataM);
        
        UN = dataN.d.RB*red_sim_N.uN;
        UM = dataM.d.RB*red_sim_M.uN;
        
        
        UN1 = dataN.d.RB1*red_sim_N.uN(1:N);
        UN2 = dataN.d.RB2*red_sim_N.uN(N+1:end);
        UM1 = dataM.d.RB1*red_sim_M.uN(1:M);
        UM2 = dataM.d.RB2*red_sim_M.uN(M+1:end);
        
        
        [utilde1,utilde2] = compute_partial_solutions(model,model_data);
        [tru1,tru2,u1,u2] = compute_trace_solutions(model,model_data,utilde1,utilde2);
        G1 = length(tru1);
        G2 = length(tru2);
        UN1 = UN1(1:G1);
        UN2 = UN2(end-G2+1:end);
        UM1 = UM1(1:G1);
        UM2 = UM2(end-G2+1:end);
        
        diff1 = UN1 - UM1;
        diff1(1) = 0;
        DELTA_part1(k) = sqrt(abs(diff1'*dataN.d.LGG1*diff1));
        
        diff2 = UN2 - UM2;
        DELTA_part2(k) = sqrt(abs(diff2'*dataN.d.LGG2*diff2)); 
        
        diff = tru1 - UM1; diff(1) = 0;
        zaehler(1,k) = sqrt(abs(diff'*dataN.d.LGG1*diff));
        diff = tru1 - UN1; diff(1) = 0;
        nenner(1,k) = sqrt(abs(diff'*dataN.d.LGG1*diff));
        
        diff = tru2 - UM2;
        zaehler(2,k) = sqrt(abs(diff'*dataN.d.LGG2*diff));
        diff = tru2 - UN2;
        nenner(2,k) = sqrt(abs(diff'*dataN.d.LGG2*diff));
        
        
        truth_err1(k) = nenner(1,k);
        truth_err2(k) = nenner(2,k);
        
    end
    
    hier_err_est1 = DELTA_part1./(1-theta_part(1,i));
    hier_err_est2 = DELTA_part2./(1-theta_part(2,i));
    
    close all; 
    figure(i)
    semilogy(truth_err1,'r'); hold on; semilogy(hier_err_est1,'g'); hold off;
    figure(i+1)
    semilogy(truth_err2,'r'); hold on; semilogy(hier_err_est2,'g'); hold off;
    
    FN_save = ['/Users/mladi/Desktop/HighDim/HierErr_onGamma_D' num2str(D) '/'];
    save([FN_save 'data' num2str(i) '.mat']);
    disp(i)
end


%%
D=3;
for i=1:15-D
   FN_LOAD = ['/Users/mladi/Desktop/HighDim/HierErr_onGamma_D' num2str(D) '/'];
   data = load([FN_LOAD 'data' num2str(i) '.mat']); 
   
   HE1(i) = sum(data.hier_err_est1)/1000;
   HE2(i) = sum(data.hier_err_est2)/1000;
   
   TE1(i) = sum(data.truth_err1)/1000;
   TE2(i) = sum(data.truth_err2)/1000;
end

figure(1)
semilogy(HE1,'g'); hold on; semilogy(TE1,'r')
figure(2)
semilogy(HE2,'g'); hold on; semilogy(TE2,'r')