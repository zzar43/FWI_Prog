function make_source_vec(model)
    source = model.source;
    coor_source = model.coor_source;
    pml_len = model.pml_len;

    Nx_pml = model.Nx + 2*pml_len;
    Ny_pml = model.Ny + 2*pml_len;
    coor_source = coor_source .+ pml_len;

    f = zeros(Nx_pml * Ny_pml, length(source));

    for i = 1:length(source)
        f[coor_source[i,1] + (coor_source[i,2]-1)*Nx_pml,i] = source[i];    
    end
    
    return f
end