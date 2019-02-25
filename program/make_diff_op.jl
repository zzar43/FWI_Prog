using SparseArrays, LinearAlgebra;

diff1(M) = [[1.0 zeros(1,M-1)] ; diagm(1 => ones(M-1)) - Matrix(1.0I,M,M)];
sdiff1(M) = sparse(diff1(M));

function make_diff_op(vel_ex, model, fre)

    pml_coef = 1;
    pml_len = model.pml_len;
    Nx_pml = model.Nx + 2*pml_len;
    Ny_pml = model.Ny + 2*pml_len;

    if size(vel_ex)[1] != Nx_pml*Ny_pml
        error("Check the model size. Need to input the extended model in vector form.")
    end
    
    pml_value = range(0, stop=pml_coef, length=pml_len);
    
    beta = zeros(Nx_pml, Ny_pml);
    for i = 1:pml_len
        beta[pml_len+1-i,:] .= pml_value[i];
        beta[end-pml_len+i,:] .= pml_value[i];
        beta[:,pml_len+1-i] .= pml_value[i];
        beta[:,end-pml_len+i] .= pml_value[i];
    end
    
    coef = (1 .+ im*beta[:]) .* ((2*pi*fre).^2) ./ (vel_ex[:].^2);
    
    Dx = sdiff1(Nx_pml) / model.dx;
    Dy = sdiff1(Ny_pml) / model.dy;
    # D2x = Dx' * Dx;
    # D2y = Dy' * Dy;
    # eyex = spdiagm(0 => ones(Nx_pml))
    # eyey = spdiagm(0 => ones(Ny_pml))
    # L1 = kron(eyey, D2x); 
    # L2 = kron(D2y,eyex);
    # L = -1 .* (L1 + L2);
    L = -1 .* (kron(spdiagm(0 => ones(Ny_pml)), Dx' * Dx) + kron(Dy' * Dy,spdiagm(0 => ones(Nx_pml))))
    
    return L = spdiagm(0 => coef) + L;
end

function preconditioner(vel_ex, model, fre)

    pml_coef = 1;
    pml_len = model.pml_len;
    Nx_pml = model.Nx + 2*pml_len;
    Ny_pml = model.Ny + 2*pml_len;

    if size(vel_ex)[1] != Nx_pml*Ny_pml
        error("Check the model size. Need to input the extended model in vector form.")
    end
    
    # pml_value = range(0, stop=pml_coef, length=pml_len);
    
    # beta = zeros(Nx_pml, Ny_pml);
    # for i = 1:pml_len
    #     beta[pml_len+1-i,:] .= pml_value[i];
    #     beta[end-pml_len+i,:] .= pml_value[i];
    #     beta[:,pml_len+1-i] .= pml_value[i];
    #     beta[:,end-pml_len+i] .= pml_value[i];
    # end
    
    coef = (1 + im) .* ((2*pi*fre).^2) ./ (vel_ex[:].^2);
    
    Dx = sdiff1(Nx_pml) / model.dx;
    Dy = sdiff1(Ny_pml) / model.dy;
    # D2x = Dx' * Dx;
    # D2y = Dy' * Dy;
    # eyex = spdiagm(0 => ones(Nx_pml))
    # eyey = spdiagm(0 => ones(Ny_pml))
    # L1 = kron(eyey, D2x); 
    # L2 = kron(D2y,eyex);
    # L = -1 .* (L1 + L2);
    L = -1 .* (kron(spdiagm(0 => ones(Ny_pml)), Dx' * Dx) + kron(Dy' * Dy,spdiagm(0 => ones(Nx_pml))))
    
    return L = spdiagm(0 => coef) - L;
end