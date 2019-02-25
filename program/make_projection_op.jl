using SparseArrays;

# make projection operator
function make_projection_op(model)
    
    coor_receiver = model.coor_receiver;
    pml_len = model.pml_len;
    Nx_pml = model.Nx + 2*pml_len;
    Ny_pml = model.Ny + 2*pml_len;
    coor_receiver = coor_receiver .+ pml_len;
    N_r = size(coor_receiver,1);

    Ip = zeros(N_r);
    Jp = zeros(N_r);
    Vp = ones(N_r);

    for i = 1:N_r
        Ip[i] = i;
        Jp[i] = coor_receiver[i,1] + (coor_receiver[i,2]-1)*Nx_pml;
    end

    P = sparse(Ip,Jp,Vp,N_r,Nx_pml*Ny_pml);
    
    return P
end