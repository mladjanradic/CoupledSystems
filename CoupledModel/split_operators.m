function d = split_operators(model,model_data)

d=model_data;

mu = model.get_mu(model);
epsilon = model.base_model.varepsilon;

model.decomp_mode = 1;
[A,r] = model.operators(model,model_data);

gamma_dofs1 = [1;model_data.gamma_dofs{2}];
gamma_dofs2 = model_data.gamma_dofs{1};

Lgamma = model_data.gamma_inner_product_matrices{1};
d.LGG2 = Lgamma(gamma_dofs2,gamma_dofs2);
d.LGG1 = d.LGG2;
%Here it is assumed, that the first node is a Dirichlet node.
%Therefore "eliminate it". This is not correct, if ther first node is not a
%Dirichlet node -> has to be adjusted otherwise.
d.LGG1(1,:) = 0;


for i=1:size(A{1},2)
    for j=1:2
        if j==1
            gamma_dofs = gamma_dofs2;
        else
            gamma_dofs = gamma_dofs1;
        end
        
        AA = A{j}{i};
        AGI_temp = AA(gamma_dofs,:);
        AGI_temp(:,gamma_dofs) = [];
        AIG_temp = AA(:,gamma_dofs);
        AIG_temp(gamma_dofs,:) = [];
        AGG_temp = AA(gamma_dofs,gamma_dofs);
        AA(gamma_dofs,:) = [];
        AA(:,gamma_dofs) = [];
        AGG{j}{i} = AGG_temp;
        AGI{j}{i} = AGI_temp;
        AIG{j}{i} = AIG_temp;
        AII{j}{i} = AA;
    end
end



for i=1:size(r{1},2)
    for j=1:2
        if j==1
            gamma_dofs = gamma_dofs2;
        else
            gamma_dofs = gamma_dofs1;
        end
        r_temp = r{j}{i};
        rG_temp = r_temp(gamma_dofs);
        r_temp(gamma_dofs) = [];
        rI_temp = r_temp;
        rI{j}{i} = rI_temp;
        rG{j}{i} = rG_temp;
    end
end

% Inner nodes "do not depend" on mu!
d.AII1 = {epsilon*AII{2}{1}+AII{2}{3}+AII{2}{5}};

d.AGI1 = {epsilon*AGI{2}{1}+AGI{2}{3}+AGI{2}{5}};
d.AIG1 = {epsilon*AIG{2}{1}+AIG{2}{3}+AIG{2}{5}};
d.AGG1 = {epsilon*AGG{2}{1}+AGG{2}{3}+AGG{2}{5}};

% Inner nodes "do not depend" on mu!
d.AII2 = {mu(1)*AII{1}{2}+mu(2)*AII{1}{4}};

d.AGI2 = {mu(1)*AGI{1}{2}+mu(2)*AGI{1}{4}};
d.AIG2 = {mu(1)*AIG{1}{2}+mu(2)*AIG{1}{4}};
d.AGG2 = {mu(1)*AGG{1}{2}+mu(2)*AGG{1}{4}};

d.r1 = 0*rI{2}{1};
d.rG1 = rG{2}{1}*0;
d.r2 = rI{1}{1}*0;
d.rG2 = rG{1}{1}*0;

for i=1:size(r{1},2)
    d.r1 = d.r1 + rI{2}{i};
    d.rG1 = d.rG1 + rG{2}{i};
    d.r2 = d.r2 + rI{1}{i};
    d.rG2 = d.rG2 + rG{1}{i};
end

d.r1 = {d.r1};
d.r2 = {d.r2};
d.rG1 = {d.rG1};
d.rG2 = {d.rG2};

d.W1 = model_data.df_infos{2}.h10_inner_product_matrix + ...
    model_data.df_infos{2}.l2_inner_product_matrix;
d.WG1 = d.W1(gamma_dofs1,gamma_dofs1);
d.W1(gamma_dofs1,:) = [];
d.W1(:,gamma_dofs1) = [];
d.W2 = model_data.df_infos{1}.h10_inner_product_matrix + ...
    model_data.df_infos{1}.l2_inner_product_matrix;
d.WG2 = d.W2(gamma_dofs2,gamma_dofs2);
d.W2(gamma_dofs2,:) = [];
d.W2(:,gamma_dofs2) = [];


end


