function [] = compute_infsup()

clear, close, clc;

FileName = '/Users/mladi/Desktop/HighDim/';

data = load('/Users/mladi/Desktop/HighDim/Sol/Mtrain.mat');
M_train = data.M_train;
clear data;

params.ynumintervals = 30;
params.xnumintervals = 50;
base_model=coupled_model(params);
params.dd_rect_corner1 = {[0,0]};
params.dd_rect_corner2 = {[1,1]};
model=dom_dec_model(base_model,params);

model.mu_ranges{2} = [-1 -2];
model_data = model.gen_model_data(model);

model.detailed_simulation = @detailed_simulation_cp;

W1 = model_data.df_infos{2}.l2_inner_product_matrix+ ...
   model_data.df_infos{2}.h10_inner_product_matrix;
W2 = model_data.df_infos{1}.l2_inner_product_matrix+ ...
   model_data.df_infos{1}.h10_inner_product_matrix;
W = blkdiag(W1,W2);
clear W1 W2;
D =model_data.df_infos{2}.dirichlet_gids;
W(D,:)=[];
W(:,D)=[];


BETA = zeros(1,size(M_train,2));
opts.maxit=1e4;
opts.isreal = false;
    
for i=1:size(M_train,2)
    mu = M_train(:,i);
    model = model.set_mu(model,mu);
    
    model.decomp_mode = 0;
    [A,~] = compute_operators(model,model_data);
    A(:,D)=[];
    A(D,:)=[];
    
    fun = @(X) Afun_inv(X,'notransp',A,W);
    [~,infsup] = eigs(fun, size(A,1), 1, 'sm', opts);
    
    infsup = abs(infsup);
    BETA(i) = infsup;
    disp(i);
end

save([FileName 'BETA.mat'],'BETA')

end


function x = Afun_inv(x,tflag,A,W)
    if strcmp(tflag,'transp')
        x = W*(A\(W*(A'\x)));
    else
        x = A\(W*(A'\(W*x)));
    end
end

