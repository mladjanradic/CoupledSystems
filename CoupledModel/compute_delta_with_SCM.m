function delta = compute_delta_with_SCM(model,d)


data = load('SCM_data_OMEGA1.mat');
operator1 = data.operator;
SCM_data1 = data.SCM_data;
data = load('SCM_data_OMEGA2.mat');
SCM_data2 = data.SCM_data;
operator2 = data.operator;
clear data;

delta = zeros(1,size(d.M_train,2));
TIME = zeros(1,size(d.M_train,2));

for i=1:size(d.M_train,2)
    
    mu = d.M_train(:,i);
    tic
    alpha1 = SCMonline(mu,model,operator1,SCM_data1,'lb');
    alpha2 = SCMonline(mu,model,operator2,SCM_data2,'lb');
    
    %model = set_mu(model,mu);
    %model.decomp_mode=0;
    %[A,~] = operator2(model);
    %gamma2 = eigs(A,d.LGG2,1,'lm');
    gamma2 = 0.5;
    beta2 = mu(4);
    gammab2b1 = abs(mu(3) - mu(4));
    
    zaehler = (1./alpha2)*beta2*gammab2b1;
    nenner = alpha2*(1/(gamma2^2))*(beta2^2)+alpha1;
    delta(i) = zaehler/nenner;
    TIME(i) = toc;
    if(delta(i)>1)
        disp(i)
        disp(gamma2)
    end
    
    %disp(i)
    
end
save('deltas.mat');
end
