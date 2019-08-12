function rb_sim_data = stabilized_lin_stat_rb_simulation_OMEGA(model,reduced_data,uN)


rb_sim_data = [];

tic
old_mode = model.decomp_mode;
model.decomp_mode = 2; % coefficients

[A_coeff,f_coeff] = ...
    model.operators(model,[]);

tic
% Computation of the Residual
%Q_r = size(reduced_data.G,1);
neg_auN_coeff = -A_coeff * uN.';
res_coeff = [f_coeff; neg_auN_coeff(:)];
res_norm_sqr = res_coeff' * reduced_data.G * res_coeff;

% prevent possibly negative numerical disturbances:
res_norm_sqr = max(res_norm_sqr,0);
res_norm = sqrt(abs(res_norm_sqr));
% if the SCM is being used perform an online-phase and use the resulting
% lower bound. Otherwise use the old code.

rb_sim_data.res_norm_sqr = abs(res_norm_sqr);
rb_sim_data.res_norm = res_norm;
rb_sim_data.res_norm_time = toc;

model.decomp_mode = old_mode;