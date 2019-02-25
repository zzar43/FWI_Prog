include("program/model.jl");
include("program/make_diff_op.jl")
include("program/make_projection_op.jl")
include("program/make_source_vec.jl")
include("program/compute_gradient.jl")
include("program/lbfgsb_fwi.jl")
# include("draw_figure.jl")

# model 2
using JLD2
@load "/data/overthrust_small.jld2" vel_true vel_init
vel_true = copy(vel_true');
vel_init = copy(vel_init');

Nx, Ny = size(vel_true)
dx = 25; dy = 25
pml_len = 30;

# source and receiver
coor_source = ones(Int, 101, 2);
source = -1 .* ones(101);

for i = 1:101
    coor_source[i,1] = 2;
    coor_source[i,2] = (i-1)*4+1;
end

coor_receiver = ones(Int,101, 2);
for i = 1:101
    coor_receiver[i,1] = 2;
    coor_receiver[i,2] = (i-1)*4+1;
end

model = Model(Nx, Ny, dx, dy, coor_source, source, coor_receiver, pml_len);

vel_true = extend_vel(vel_true, model);
vel_init = extend_vel(vel_init, model);
Nx_pml = model.Nx + 2*pml_len;
Ny_pml = model.Ny + 2*pml_len;
vel_true = reshape(vel_true, Nx_pml*Ny_pml);
vel_init = reshape(vel_init, Nx_pml*Ny_pml);

# Minimization
lower = minimum(vel_true)
upper = maximum(vel_true)
fre_vec = [4 8 12 16]
result = lbfgsb_fwi(vel_true, vel_init, model, fre_vec, lower, upper; iterTime=10);

@save "example2.jld2" result