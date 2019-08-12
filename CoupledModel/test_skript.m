% Uses the functions by Immanuel Martini (former Immanuel Maier)
% (domain decomposition), see "An Iterative Domain
% Decomposition Procedure for The Reduced Basis Method"
clear, close, clc;



params.numintervals = 15;
base_model=coupled_model(params);
params.dd_rect_corner1 = {[0,0]};
params.dd_rect_corner2 = {[1,1]};
model=dom_dec_model(base_model,params);
M_train = rand_uniform(100, model.mu_ranges);
M_train(:,1) = [2;-1;1;1];
model = model.set_mu(model, M_train(:,1));
model_data = model.gen_model_data(model);

model.detailed_simulation = @detailed_simulation_cp;

model.decomp_mode = 1;
[A,r] = model.operators(model,model_data);

model.decomp_mode = 2;
[Acoeff,rcoeff] = model.operators(model,model_data);

sim_data = model.detailed_simulation(model,model_data);

[utilde1,utilde2] = compute_partial_solutions(model,model_data);
[tru1,tru2,u1,u2] = compute_trace_solutions(model,model_data,utilde1,utilde2);
disp('Small comparison')
error = norm(sim_data.uh.dofs - [tru1;u1;u2;tru2]);
disp(error)





