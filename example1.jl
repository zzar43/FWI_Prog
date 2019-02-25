include("program/model.jl");
include("program/make_diff_op.jl")
include("program/make_projection_op.jl")
include("program/make_source_vec.jl")
include("program/compute_gradient.jl")
include("program/lbfgsb_fwi.jl")
# include("program/draw_figure.jl")

# model 1
using ImageFiltering, JLD2

Nx = 51; Ny = 51;
dx = 20; dy = 20;
pml_len = 30;

# source and receiver
coor_source = ones(Int, 51, 2);
source = -1 .* ones(51);

for i = 1:51
    coor_source[i,1] = 2;
    coor_source[i,2] = i;
end

coor_receiver = ones(Int,102, 2);
for i = 1:102
    coor_receiver[i,1] = 50;
    coor_receiver[i,2] = i;
end
for i = 52:102
    coor_receiver[i,1] = 2;
    coor_receiver[i,2] = i-51;
end

model = Model(Nx, Ny, dx, dy, coor_source, source, coor_receiver, pml_len);

# Build the velocity model
pml_len = model.pml_len;
Nx_pml = model.Nx + 2*pml_len;
Ny_pml = model.Ny + 2*pml_len;

vel0 = 2000*ones(model.Nx, model.Ny);
vel = 2000*ones(model.Nx, model.Ny);
vel[21:36, :] .= 2100;
vel[37:51, :] .= 2200;
# vel[21:31, 21:31] .= 2100;
vel0_ex = extend_vel(vel0, model);
vel_ex = extend_vel(vel, model);
vel0_ex = imfilter(vel_ex, Kernel.gaussian(10));

vel_ex = reshape(vel_ex, Nx_pml*Ny_pml);
vel0_ex = reshape(vel0_ex, Nx_pml*Ny_pml);

fre_vec = [4 8 12 16 20]
lower = 2000
upper = 2200
result = lbfgsb_fwi(vel_ex, vel0_ex, model, fre_vec, lower, upper; iterTime=10);

@save "example1.jld2" result
