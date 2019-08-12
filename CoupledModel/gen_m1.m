function [m1] = gen_m1(m1,detailed_data)

m1.decomp_mode = 1;
m1.operators=@operators;
m1.compute_output_functional=0;
m1.use_scm = 0;
m1.detailed_simulation = @detailed_simulation;

end

function [A,r] = operators(model,detailed_data)
    if model.decomp_mode==1
        A = detailed_data.AII1;
        r = detailed_data.r1;
    end
    if model.decomp_mode==2
        A=1;
        r=1;
    end
    if model.decomp_mode==0
        A = detailed_data.AII1{1};
        r=0;
        for i=1:size(detailed_data.r1,2)
            r = r+detailed_data.r1{i};
        end
    end
    
end



function uh = detailed_simulation(model,detailed_data)

    model.decomp_mode = 0;
    [A,r] = model.operators(model,detailed_data);
    uh = A\r;

end