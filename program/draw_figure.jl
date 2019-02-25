using PyPlot

function draw_model(vel, model; ex_model=false)

    Nx = model.Nx;
    Ny = model.Ny;
    Nx_pml = model.Nx + 2*model.pml_len;
    Ny_pml = model.Ny + 2*model.pml_len;
    coor_source = model.coor_source;
    coor_receiver = model.coor_receiver;
    pml_len = model.pml_len;

    if size(vel)[1] != Nx_pml*Ny_pml
        error("Check the model size. Need to input the extended model in vector form.")
    end
    vel = reshape(vel, Nx_pml, Ny_pml);

    if ex_model==true
        imshow(vel); colorbar(shrink=0.5)
        scatter(coor_source[:,2] .+ pml_len.-1, coor_source[:,1] .+ pml_len.-1, alpha=1)
        scatter(coor_receiver[:,2] .+ pml_len.-1, coor_receiver[:,1] .+ pml_len.-1, alpha=0.5)
    else
        vel = vel[(pml_len+1):(Nx+pml_len), (pml_len+1):(Ny+pml_len)];
        imshow(vel); colorbar(shrink=0.5)
        scatter(coor_source[:,2].-1, coor_source[:,1].-1)
        scatter(coor_receiver[:,2].-1, coor_receiver[:,1].-1, alpha=0.5)
    end
end

function draw_model2(vel1, vel2, model; vertial=true)
    Nx = model.Nx;
    Ny = model.Ny;
    Nx_pml = model.Nx + 2*model.pml_len;
    Ny_pml = model.Ny + 2*model.pml_len;
    coor_source = model.coor_source;
    coor_receiver = model.coor_receiver;
    pml_len = model.pml_len;

    if size(vel1)[1] != Nx_pml*Ny_pml
        error("Check the model size. Need to input the extended model in vector form.")
    end
    if size(vel2)[1] != Nx_pml*Ny_pml
        error("Check the model size. Need to input the extended model in vector form.")
    end
    vel1 = reshape(vel1, Nx_pml, Ny_pml);
    vel2 = reshape(vel2, Nx_pml, Ny_pml);
    vel1 = vel1[(pml_len+1):(Nx+pml_len), (pml_len+1):(Ny+pml_len)];
    vel2 = vel2[(pml_len+1):(Nx+pml_len), (pml_len+1):(Ny+pml_len)];

    if vertial==true
        subplot(211)
        imshow(vel1); colorbar(shrink=0.5)
        scatter(coor_source[:,2].-1, coor_source[:,1].-1)
        scatter(coor_receiver[:,2].-1, coor_receiver[:,1].-1, alpha=0.5)
        subplot(212)
        imshow(vel2); colorbar(shrink=0.5)
        scatter(coor_source[:,2].-1, coor_source[:,1].-1)
        scatter(coor_receiver[:,2].-1, coor_receiver[:,1].-1, alpha=0.5)
        tight_layout()
    else
        subplot(121)
        imshow(vel1); colorbar(shrink=0.5)
        scatter(coor_source[:,2].-1, coor_source[:,1].-1)
        scatter(coor_receiver[:,2].-1, coor_receiver[:,1].-1, alpha=0.5)
        subplot(122)
        imshow(vel2); colorbar(shrink=0.5)
        scatter(coor_source[:,2].-1, coor_source[:,1].-1)
        scatter(coor_receiver[:,2].-1, coor_receiver[:,1].-1, alpha=0.5)
        tight_layout()
    end

end



function draw_wavefield(u, model; index=1, ex_model=false)
    Nx = model.Nx;
    Ny = model.Ny;
    Nx_pml = model.Nx + 2*model.pml_len;
    Ny_pml = model.Ny + 2*model.pml_len;
    coor_source = model.coor_source;
    coor_receiver = model.coor_receiver;
    pml_len = model.pml_len;
    max_value = 0.8*maximum(real(u[:,index]))

    if index > size(u)[2] || index <= 0
        error("Check the index.")
    end
    if size(u)[1] != Nx_pml*Ny_pml
        error("Check the wavefield size. Need to input the extended wavefield in vector form.")
    end

    uu = u[:,index];
    uu = real(reshape(uu, Nx_pml, Ny_pml));

    if ex_model==true
        imshow(uu, cmap=:seismic, vmin=-max_value, vmax=max_value); colorbar(shrink=0.5)
        scatter(coor_source[:,2] .+ pml_len.-1, coor_source[:,1] .+ pml_len.-1)
        scatter(coor_receiver[:,2] .+ pml_len.-1, coor_receiver[:,1] .+ pml_len.-1, alpha=0.5)
        tight_layout()
    else
        uu = uu[(pml_len+1):(Nx+pml_len), (pml_len+1):(Ny+pml_len)];
        imshow(uu, cmap=:seismic, vmin=-max_value, vmax=max_value); colorbar(shrink=0.5)
        scatter(coor_source[:,2].-1, coor_source[:,1].-1)
        scatter(coor_receiver[:,2].-1, coor_receiver[:,1].-1, alpha=0.5)
        tight_layout()
    end
end