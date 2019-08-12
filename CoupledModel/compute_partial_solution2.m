function [u2] = compute_partial_solution2(model,model_data)

    model.decomp_mode = 0;
    [A,r] = model.operators(model,model_data);

    u2 = A\r;
    
end