function model = coupled_model(params)

if nargin<1
    params = [];
end;

if ~isfield(params,'numintervals')
    params.numintervals = 16;
end;
if ~isfield(params,'pdeg')
    params.pdeg = 1;
end;
if ~isfield(params,'detailed_train_savepath')
    params.detailed_train_savepath = '/coupled_systems/tmp';
end;

model = lin_stat_model_default;
model.decomp_mode = 0; % default
model.dual_mode = 0;

% parameters:
model.varepsilon = 1e-3;
model.mu_names = {'d','r','e_1','e_2'};
model.mu_ranges = {[1,2],[3,4],[0.5,1.5],[0.5,1.5]};

model.RB_generation_mode = 'greedy_uniform_fixed';
model.RB_numintervals = 10;
model.detailed_train_savepath = params.detailed_train_savepath;
model.filecache_ignore_fields_in_model = {};


% default parameter values:
model = set_mu(model,[1,3,1,1]);
model.get_mu = @(model) model.mus;
model.mus = [1,3,1,1];

% data function specification:
model.has_diffusivity = 1;
model.has_advection = 1;
model.has_reaction = 1;
model.has_source = 1;
model.has_output_functional = 0;
model.has_dirichlet_values = 1;
model.has_neumann_values = 0;
model.has_robin_values = 0;
model.compute_output_functional = 0;

% source:
source_coefficients = @(grid,eindices,loc,params) 0;
source_components = @my_source_components;

model.source = @my_source_function;
model.source_function = @(grid,eindices,loc,params) ...
    eval_affine_decomp_general(source_components,source_coefficients,...
    grid,eindices,loc,params);

% diffusivity_tensor:
diffusivity_tensor_coefficients = @(grid,eincices,loc,params) ...
    [params.varepsilon;params.d];
diffusivity_tensor_components = @my_diffusivity_tensor_components;

model.diffusivity_tensor = @(grid,eindices,loc,params) ...
    eval_affine_decomp_general(diffusivity_tensor_components, ...
    diffusivity_tensor_coefficients, ...
    grid,eindices,loc,params);


% advection_tensor:
advection_tensor_coefficients = @(grid,eincices,loc,params) 1;
advection_tensor_components = @my_advection_tensor_components;

model.velocity = @(grid,eindices,loc,params) ...
    eval_affine_decomp_general(advection_tensor_components, ...
    advection_tensor_coefficients, ...
    grid,eindices,loc,params);


% reaction:
reaction_coefficients = @(grid,eincices,loc,params) params.r;
reaction_components = @my_reaction_components;

model.reaction = @(grid,eindices,loc,params) ...
    eval_affine_decomp_general(reaction_components, ...
    reaction_coefficients, ...
    grid,eindices,loc,params);



model.uD = 1;
% dirichlet_values:
dirichlet_coefficients = @(grid,eindices,face_index,lloc,params) 1;
dirichlet_components = @(grid,eindices,face_index,lloc,params) ...
    { params.uD*ones(size(eindices,1),1) };


model.dirichlet_values = @(grid,eindices,face_index,lloc,params) ...
    eval_affine_decomp_general(dirichlet_components, ...
    dirichlet_coefficients,grid, ...
    eindices,face_index,lloc,params);


model.boundary_type=@my_boundary_type;



% geometry settings:
model.gridtype = 'triagrid';
model.xnumintervals = params.xnumintervals;
model.ynumintervals = params.ynumintervals;
model.xrange = [0,1];
model.yrange = [0,1.5];

% discretization parameters:
model.pdeg = params.pdeg;
model.qdeg = 3;
model.dimrange = 1;

% for error estimation:
model.coercivity_alpha = @my_coercivity_alpha;
model.continuity_gamma = @my_continuity_gamma;

end

%%%%%%%%%%%% auxiliary functions:

function res = my_source_components(grid,eindices,loc,params)
% if isa(grid,'feminfo')
%     glob = local2global(grid.grid,eindices,loc,params);
% else
%     glob = local2global(grid,eindices,loc,params);
% end
res = {0};
end


function res = my_source_function(grid,eindices,loc,params)
if params.dual_mode
    res = params.output_function(grid,eindices,loc,params);
    if ~iscell(res)
        res = -res;
    end
else
    res = params.source_function(grid,eindices,loc,params);
end
end

% diffusion on Omega1 AND Omega2
function res = my_diffusivity_tensor_components(grid,eindices,loc,params)
CX = grid.CX(eindices);
CY = grid.CY(eindices);
res = {[CY>1, zeros(length(CX),2), CY>1], ...
    [CY<1, zeros(length(CX),2), CY<1]};
end

% advection only happens on Omega2
function res = my_advection_tensor_components(grid,eindices,loc,params)
CX = grid.CX(eindices);
CY = grid.CY(eindices);
res = {0.5*[(1-CX).*(CY>1), (CY-1.5).*(CY>1)]};
end


% reaction only happens on Omega2
function res = my_reaction_components(grid,eindices,loc,params)
if isa(grid,'feminfo')
grid = grid.grid;
end
CY = grid.CY(eindices);
res = {[CY<1]};
end




function res = my_coercivity_alpha(model)
mus = model.get_mu(model);
res = min([1;mus(1:end-1)]);
end


function res = my_continuity_gamma(model)
mus = model.get_mu(model);
res = max([1;mus(1:end-1)]);
end



function res = my_boundary_type(glob,params)
res = zeros(size(glob,1),1);
i  = find(glob(:,1)<=1e-10);
i  = [i; find(glob(:,1)>=1-1e-10)];
i  = [i; find(glob(:,2)<=1e-10)];
i  = [i; find(glob(:,2)>=1-1e-10)];
res(i) = -2;

res(glob(:,1)<=1e-10 & glob(:,2)>1-eps) = -1;
end