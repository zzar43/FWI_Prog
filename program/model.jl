struct Model
    Nx::Int
    Ny::Int
    dx
    dy
    coor_source
    source
    coor_receiver
    pml_len
end

function extend_vel(vel, model)
    
    pml_len = model.pml_len;
    Nx_pml = model.Nx + 2*pml_len;
    Ny_pml = model.Ny + 2*pml_len;
    
    vel_ex = zeros(Nx_pml, Ny_pml);
    vel_ex[pml_len+1:end-pml_len, pml_len+1:end-pml_len] .= vel;
    for i = 1:pml_len
        vel_ex[i,:] = vel_ex[pml_len+1,:];
        vel_ex[end-i+1,:] = vel_ex[end-pml_len,:];
        vel_ex[:,i] = vel_ex[:,pml_len+1];
        vel_ex[:,end-i+1] = vel_ex[:,end-pml_len];
    end
    
    return vel_ex
end