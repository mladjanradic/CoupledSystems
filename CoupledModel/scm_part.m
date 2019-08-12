% SCM:

for OMEGA = 1:2
    if OMEGA==1
        operator=@ (model) operators_for_SCM_OMEGA1(model,d);
        dir=1;
    else
        operator=@ (model) operators_for_SCM_OMEGA2(model,d);
    end
    %dir = model_data.df_infos{2}.dirichlet_gids;
    
    %dir(end) = [];
    %dir = [];
    K = d.LGG2;
    params.Ma = inf;
    params.Mp = 0;
    params.NumberOfParameters = 1000;
    params.eps=0.01;
    params.max_iterations=100;
    params.parallel=0;
    params.dirichlet_gids = dir;
    params.type='coercive';
    params.ParameterSet = M_train;
    
    SCM_data=SCM(model,params,operator,K);
    %load('SCM_data_OMEGA1.mat');
    
    %% Test:
    
    rmserror = 0;
    maderror = 0;
    
    Ntest = 20;
    Mtest = zeros(4,Ntest);
    ALPHAUB = zeros(1,Ntest);
    ALPHALB = zeros(1,Ntest);
    COER = zeros(1,Ntest);
    timeUB = zeros(1,Ntest);
    timeLB = zeros(1,Ntest);
    err_ub = zeros(1,Ntest);
    err_lb = zeros(1,Ntest);
    
    K=remove_dirichlet(K,0,dir);
    fprintf('Results:\ncoercive constant      upper bound       lower bound (scm coercive constant)\n');
    for i = 1:20
        mu = rand_uniform(1,model.mu_ranges);
        
        Mtest(:,i) = mu;
        
        model = set_mu(model,mu);
        model.decomp_mode=0;
        [A,~] = operator(model);
        A=remove_dirichlet(A,0,dir);
        A=0.5*(A+A');
        [v,e] = eig( full(A),full(K) );
        coer = min(max(e));
        COER(i) = coer;
        tic
        alpha_UB =  SCMonline(mu,model,operator,SCM_data,'ub');
        timeUB(i) = toc;
        ALPHAUB = alpha_UB;
        %ONLINE:
        tic
        alpha_LB = SCMonline(mu,model,operator,SCM_data,'lb');
        timeLB(i) = toc;
        ALPHALB(i) = alpha_UB;
        fprintf('%f               %f          %f\n',coer,alpha_UB,alpha_LB);
        
        rmserror = rmserror+(coer-alpha_LB)^2;
        maderror = maderror + abs(coer-alpha_LB);
        
        err_ub(i) = (coer-alpha_UB);
        err_lb(i) = (coer-alpha_LB);
    end
    rmserror = sqrt(rmserror/20);
    maderror = maderror/20;
    meantimeUB = sum(timeUB)/Ntest;
    meantimeLB = sum(timeLB)/Ntest;
    disp(meantimeUB)
    disp(meantimeLB)
    save(['SCM_data_OMEGA' num2str(OMEGA)]);
end