function [A,r] = operators_for_SCM_global_1(model,detailed_data)

if model.decomp_mode==1
    A = detailed_data.AII1;
    for i=1:length(A)
        A{i} = 0.5*(A{i}+A{i}');
    end
    r = detailed_data.r1;
end
if model.decomp_mode==2
    A=[model.base_model.varepsilon;1;1];
    r=ones(4,1);
end
if model.decomp_mode==0
    A = model.base_model.varepsilon*detailed_data.AII1{1} + ...
        detailed_data.AII1{2}+detailed_data.AII1{3};
    A=0.5*(A+A');
    r=0;
    for i=1:size(detailed_data.r1,2)
        r = r+detailed_data.r1{i};
    end
end

end