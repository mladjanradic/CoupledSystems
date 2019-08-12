% Uses the functions by Immanuel Martini (former Immanuel Maier)
% (domain decomposition), see "An Iterative Domain
% Decomposition Procedure for The Reduced Basis Method"
clear, close, clc;

%%

params.ynumintervals = 30;
params.xnumintervals = 50;
base_model=coupled_model(params);
params.dd_rect_corner1 = {[0,0]};
params.dd_rect_corner2 = {[1,1]};
model=dom_dec_model(base_model,params);

model.RB_numintervals=[100;100];
par.numintervals = model.RB_numintervals;
par.range = {model.mu_ranges{3:4}};
MMesh0 = cubegrid(par);
M_train = get(MMesh0,'vertex')';

M_train = [ones(size(M_train)) ; M_train];

M_train(2,:) = 3*M_train(2,:);

%%
model = model.set_mu(model, M_train(:,1));
model_data = model.gen_model_data(model);

model.detailed_simulation = @detailed_simulation_cp;

FileName1 = '/Users/mladi/Desktop/GAMMA/PartSol/';
FileName2 = '/Users/mladi/Desktop/GAMMA/Sol/';

mkdir(FileName1);
mkdir(FileName2);

for i=1:size(M_train,2)
    mu = M_train(:,i);
    model = model.set_mu(model,mu);
    
    sim_data = model.detailed_simulation(model,model_data);
    uh = sim_data.uh.dofs;
    FN = [FileName2 'sol' num2str(i) '.mat'];
    save(FN,'uh');
    
    [utilde1,utilde2] = compute_partial_solutions(model,model_data);
    [tru1,tru2,u1,u2] = compute_trace_solutions(model,model_data,utilde1,utilde2);
    
    FN = [FileName1 'sol' num2str(i) '.mat'];
    %u1 und u2 sind die auf dem Interface!
    save(FN,'utilde1','utilde2','tru1','tru2','u1','u2');
    disp(i);
end

FN = [FileName1 'Mtrain.mat'];
save(FN,'M_train');
FN = [FileName2 'Mtrain.mat'];
save(FN,'M_train');

