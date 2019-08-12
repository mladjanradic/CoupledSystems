function d = greedy_hierarchical(model,model_data,m,d)

N=5;
M=7;

dataN = load(['/Users/mladi/Desktop/Error/err' num2str(N) '.mat']);
dataM = load(['/Users/mladi/Desktop/Error/err' num2str(M) '.mat']);

M_train = dataN.d.M_train;
dN = dataN.d;
dM = dataM.d;

rN = dataN.red_data;
rM = dataM.red_data;

clear dataN dataM;

P = size(M_train,2);
zaehler = zeros(1,P);
nenner = zeros(1,P);
hier_err = zeros(1,P);

for i=1:P
    mu = M_train(:,i);
    m.model = m.model.set_mu(m.model,mu);
    model = model.set_mu(model,mu);
    
    load(['/Users/mladi/Desktop/Sol/sol' num2str(i) '.mat']);
    red_sim = stabilized_lin_stat_rb_simulation(m,rN);
    uN = red_sim.uN;
    red_sim = stabilized_lin_stat_rb_simulation(m,rM);
    uM = red_sim.uN;
    
    UN = dN.RB*uN;
    UM = dM.RB*uM;
    
    diff = uh - UN;
    diff(model_data.df_infos{2}.dirichlet_gids) = 0;
    nenner(i) = abs(diff'*dN.W*diff);
    diff = uh - UM;
    diff(model_data.df_infos{2}.dirichlet_gids) = 0;
    zaehler(i) = abs(diff'*dN.W*diff);
    
    diff = UM - UN;
    hier_err(i) = sqrt(abs(diff'*dN.W*diff));
    
    disp(i);
end

ind = find(nenner > 1e-10);
d.theta = max(zaehler(ind)./nenner(ind));
d.truth_err = sqrt(nenner);
d.hier_err = hier_err;
plot(d.truth_err,'r'); hold on; plot(1/(1-d.theta)*d.hier_err,'g');
keyboard
