edit main
% Uses the functions by Immanuel Martini (former Immanuel Maier)
% (domain decomposition), see "An Iterative Domain
% Decomposition Procedure for The Reduced Basis Method"
clear, close, clc;


params.ynumintervals = 30;
params.xnumintervals = 50;
params.numintervals = 15;
base_model=coupled_model(params);
params.dd_rect_corner1 = {[0,0]};
params.dd_rect_corner2 = {[1,1]};
model=dom_dec_model(base_model,params);
%M_train = rand_uniform(100, model.mu_ranges);

model.mu_ranges{2} = [-1,-2];
model.RB_numintervals=[10;10;10;10];
par.numintervals = model.RB_numintervals;
par.range = model.mu_ranges;
MMesh0 = cubegrid(par);
M_train = get(MMesh0,'vertex')';


model = model.set_mu(model, M_train(:,1));
model_data = model.gen_model_data(model);

model.detailed_simulation = @detailed_simulation_cp;


d = split_operators(model,model_data);


m.varepsilon = model.base_model.varepsilon;
m.mus = M_train(:,1);


%% Greedy over Omega1 & Omega2
d.M_train = M_train;


%% Greedy over Gamma

model = model.set_mu(model, M_train(:,1));
mcp = gen_mcp(model,model_data);

d.gamma_inner_product_matrices = model_data.gamma_inner_product_matrices;
d.gamma_dofs = model_data.gamma_dofs;
W1 = model_data.df_infos{2}.l2_inner_product_matrix+ ...
   model_data.df_infos{2}.h10_inner_product_matrix;
W2 = model_data.df_infos{1}.l2_inner_product_matrix+ ...
   model_data.df_infos{1}.h10_inner_product_matrix;
d.W = blkdiag(W1,W2);
clear W1 W2;

%d = greedy_std_err_est_cp(mcp,model,model_data,d);
%d = greedy_strong_cp(mcp,model,model_data,d);
d = start_greedy_over_coupled_system(mcp,model,model_data,d);
