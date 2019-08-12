function [A,r] = operators_for_SCM_global_2(model,d)

mu = model.get_mu(model);

if model.decomp_mode == 1
    A = d.AII2;
    r = d.r2;
end
if model.decomp_mode == 2
    A = mu([1,2]);
    r = ones(4,1);
end
if model.decomp_mode==0
    A = mu(1)*d.AII2{1} + mu(2)*d.AII2{2};
    r=0;
    for i=1:size(d.r2,2)
        r = r+d.r2{i};
    end
end
end