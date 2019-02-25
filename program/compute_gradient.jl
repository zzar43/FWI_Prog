function compute_gradient(vel_ex, model, fre, d, q, Pr)

    pml_len = model.pml_len;
    Nx_pml = model.Nx + 2*pml_len;
    Ny_pml = model.Ny + 2*pml_len;
    if size(vel_ex)[1] != Nx_pml*Ny_pml
        error("Check the model size. Need to input the extended model in vector form.")
    end

    A0 = make_diff_op(vel_ex, model, fre);
    u0 = A0 \ q;
    q_ad = -1 .* Pr' * (Pr*u0 - d);
    p = conj(A0) \ q_ad;
    g = real((2*pi*fre)^2 .* conj(u0) .* p);
    g = sum(g,dims=2);
    # g = -1 * g / maximum(abs.(g));
    # g = -1 * g / norm(g);
    g = -1 * g / norm(u0)^2;
    # f = 0.5 * norm(real(Pr*u0-d))^2;
    
    # return f, g
    return g
end

function compute_gradient_lbfgsb(vel_ex, model, fre, d, q, Pr)

    pml_len = model.pml_len;
    Nx_pml = model.Nx + 2*pml_len;
    Ny_pml = model.Ny + 2*pml_len;
    if size(vel_ex)[1] != Nx_pml*Ny_pml
        error("Check the model size. Need to input the extended model in vector form.")
    end

    A0 = make_diff_op(vel_ex, model, fre);
    u0 = A0 \ q;
    q_ad = -1 .* Pr' * (Pr*u0 - d);
    p = conj(A0) \ q_ad;
    g = real((2*pi*fre)^2 .* conj(u0) .* p);
    g = sum(g,dims=2);
    g = -1 * g / norm(u0)^2;
    # g = -1 * g / maximum(abs.(g));
    # g = -1 * g / norm(g);
    # f = 0.5 * norm(real(Pr*u0-d))^2;
    f = 0.5 * norm(real(Pr*real(u0)-d))^2;
    
    return f, g
end

function compute_misfit(vel_ex, model, fre, d, q, Pr)

    pml_len = model.pml_len;
    Nx_pml = model.Nx + 2*pml_len;
    Ny_pml = model.Ny + 2*pml_len;
    if size(vel_ex)[1] != Nx_pml*Ny_pml
        error("Check the model size. Need to input the extended model in vector form.")
    end

    A0 = make_diff_op(vel_ex, model, fre);
    u0 = A0 \ q;
    
    misfit = 0.5 * norm(real(Pr*u0-d))^2;
    return misfit
end


# 
function pivot_mat(model)
    coor_receiver = model.coor_receiver .+ model.pml_len;
    Nx_pml = model.Nx + 2*model.pml_len;
    Ny_pml = model.Ny + 2*model.pml_len;
    nr = size(coor_receiver,1);

    ii = 1:nr;
    jj = zeros(Int, 1,nr);
    for i = 1:nr
       jj[i] = coor_receiver[i,1]+(coor_receiver[i,2]-1)*Nx_pml;
    end
    kk = ones(1,nr);
    R = sparse(1I, Nx_pml*Ny_pml, Nx_pml*Ny_pml);
    for i = 1:nr
        temp = R[jj[i],:];
        R[jj[i],:] = R[i,:];
        R[i,:] = temp;
    end
    return nr, R
end

function compute_gradient_partial(vel, model, fre, d, q, Pr, Nr, R)
    
    lambda = 1;
    
    A0 = make_diff_op(vel0_ex, model, fre);
    A_temp = A0 * R';
    A1 = A_temp[:, 1:Nr];
    A2 = A_temp[:, Nr+1:end];
    q_temp = q - A1 * d;
#     u2 = A2 \ q_temp;
    u2 = (A2'*A2) \ (A2'*q_temp);
    u_temp = [d; u2];
    u0 = R' * u_temp;
    
    v = A0 * u0 - q;
    g = real(sum(conj(u0) .* v, dims=2));
    g = -1 * g / maximum(abs.(g));
    # f = 0.5 * norm(real(Pr*real(u0)-d))^2 + (lambda^2/2)*norm(real(v))^2;
    f = 0.5 * norm(real(Pr*real(u0)-d))^2;
    
    return f, g
end
# 
function compute_gradient_penalty(vel, model, fre, d, q, Pr, lambda)
    A = make_diff_op(vel, model, fre);
#     A_lambda = [lambda * A; Pr];
#     q_lamda = [lambda * q; d]; 
    A_lambda = lambda*(A'*A) + (Pr'*Pr);
    q_lambda = Pr'*d + lambda*A'*q;
    u_lambda = A_lambda \ q_lambda;
    v = lambda^2/2*(A * u_lambda - q);
    g = (2*pi*fre)^2 * real(sum(conj(u_lambda) .* v, dims=2));
    g = -1 * g / norm(u_lambda)^2;
    # f = 0.5 * norm(real(Pr*real(u_lambda)-d))^2 + (lambda^2/2)*norm(real(v))^2;
    f = 0.5 * norm(real(Pr*real(u_lambda)-d))^2 + lambda^2/2 * norm(v)^2;
    f = f * 10000;
    return f, g
end

