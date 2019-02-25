function linesearch(vel, fre, model, q, d, Pr, f0, g0, s0, lambda, max_value, min_value)
    lsiter = 0;
    c1 = 1e-15;
    c2 = 0.95;
    done = false;
    mu = 0;
    nu = Inf;
    lambda = .5*lambda;

    Nx_pml = model.Nx + 2*model.pml_len
    Ny_pml = model.Ny + 2*model.pml_len
    v = similar(vel)
    ft = 0;

    while done == false
        if nu < Inf
            lambda = (nu + mu)/2;
        else
            lambda = 2*lambda;
        end

        if lsiter < 10
            v = vel + lambda * reshape(s0, Nx_pml, Ny_pml);
            v[v .> max_value] .= max_value;
            v[v .< min_value] .= min_value;

            ft, gt = compute_gradient(v, model, fre, d, q, Pr);
            lsiter = lsiter + 1;
        else
            lambda = 0;
            println("Line search failed.")
            break;
        end

        if ft > (f0 .+ c1 .* lambda .* g0' * s0)[1]
            nu = lambda;
        elseif (gt' * s0)[1] < (c2 .* g0' * s0)[1]
            mu = lambda;
        else
            done = true;
        end
    end

    return lambda, ft
end

function sd_optimization(vel0, model, fre, d, q, Pr, lambda, max_value, min_value; itermax=10)

    Nx_pml = model.Nx + 2*model.pml_len
    Ny_pml = model.Ny + 2*model.pml_len

    vel = copy(vel0);
    iter_time = 0
    p = similar(q[:,1])
    while iter_time < itermax
        f0, g0 = compute_gradient(vel, model, fre, d, q, Pr)
        
        p = -1 .* g0 / maximum(abs.(g0));
        # draw_wavefield(p, model)

        step_size, ft = linesearch(vel, fre, model, q, d, Pr, f0, g0, p, lambda, max_value, min_value)
        
        # update
        vel = vel + step_size * reshape(p, Nx_pml, Ny_pml);
        vel[vel .> max_value] .= max_value;
        vel[vel .< min_value] .= min_value;
        iter_time += 1

        # disp
        println("iter time: ", iter_time, " | step size: ", step_size, " | misfit: ", ft)
    end

    return vel
end


function sd_opt(f, g, x0, lambda, lower, upper; maxIter=10)
    # initialization
    done = false
    iterTime = 0
    
    while done == false

        iterTime += 1
        if iterTime >= maxIter
            done = true
        end
    end

end
# function lbfgs_opt()

#     # initialization

#     n = length()



# end