function [m2] = gen_m2(m2,d)

m2.decomp_mode=1;
m2.operators=@operators;
m2.compute_output_functional = 0;
m2.use_scm = 0;

m2.detailed_simulation = @detailed_simulation;
end

function [A,r] = operators(model,d)
    if model.decomp_mode == 1
        A = d.AII2;
        r = d.r2;
    end
    if model.decomp_mode == 2
        A = 1;
        r = 1;
    end
    if model.decomp_mode==0
        A = d.AII2{1};
        r=0;
        for i=1:size(d.r2,2)
            r = r+d.r2{i};
        end
    end
end



function u = detailed_simulation(model,detailed_data)
    model.decomp_mode=0;
    [A,r] = model.operators(model,detailed_data);
    u = A\r;
end
    