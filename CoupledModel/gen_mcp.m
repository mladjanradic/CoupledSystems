function m = gen_mcp(model,model_data)
m.model = model;
m.model_data = model_data;
m.operators = @operators;
m.detailed_simulation = @detailed_simulation;
m.decomp_mode=1;
m.get_inner_product_matrix = @(d) d.W;
m.compute_output_functional = 0;
m.get_rb_size = @(md,d) size(d.RB,2);
m.use_scm = 0;
end


function u = detailed_simulation(m,d)

m.decomp_mode = 1;
[AA,rr] = m.operators(m,d);
m.decomp_mode = 2;
[Acoeff,rcoeff] = m.operators(m,d);

A = lincomb_sequence(AA,Acoeff);
r = lincomb_sequence(rr,rcoeff);

u=A\r;
end