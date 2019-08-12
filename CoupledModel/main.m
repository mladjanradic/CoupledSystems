% Uses the functions by Immanuel Martini (former Immanuel Maier)
% (domain decomposition), see "An Iterative Domain
% Decomposition Procedure for The Reduced Basis Method"
clear, close, clc;

%% Preparation:
params.ynumintervals = 30;
params.xnumintervals = 50;
params.numintervals = 150;
base_model = coupled_model(params);
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

%% Generate model_data:

model = model.set_mu(model, M_train(:,1));
model_data = model.gen_model_data(model);
tic
model.detailed_simulation = @detailed_simulation_cp;
toc
%% next cell;
d = split_operators(model,model_data);
data = load('/Users/mladi/Desktop/GAMMA/PartSol/sol1.mat');
d.utilde1 = data.utilde1;
d.utilde2 = data.utilde2;

%% Greedy over Omega1 & Omega2
d.M_train = M_train;
[d,m1] = start_greedy_over_Omega_1(model,d);
[d,m2] = start_greedy_over_Omega_2(model,d);


%% Greedy over Gamma

model = model.set_mu(model, M_train(:,1334));
mcp = gen_mcp(model,d);
mcp.u = mcp.detailed_simulation(mcp,d);
d.gamma_inner_product_matrices = model_data.gamma_inner_product_matrices;
d.gamma_dofs = model_data.gamma_dofs;
W1 = model_data.df_infos{2}.l2_inner_product_matrix+ ...
    model_data.df_infos{2}.h10_inner_product_matrix;
W2 = model_data.df_infos{1}.l2_inner_product_matrix+ ...
    model_data.df_infos{1}.h10_inner_product_matrix;
d.W = blkdiag(W1,W2);
clear W1 W2;
d = start_std_greedy_over_GAMMA(mcp,d);

