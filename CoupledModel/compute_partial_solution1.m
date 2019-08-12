function [u1] = compute_partial_solution1(model,model_data)

    model.decomp_mode = 0;
    [A,r] = model.operators(model,model_data);

    u1 = A\r;
    
end