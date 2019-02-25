using LBFGSB
# reference: https://github.com/Gnimuc/LBFGSB.jl

function lbfgsb_fwi(vel_ex, vel0_ex, model, fre_vec, lower, upper; iterTime=10)

    # Initialization
    Pr = make_projection_op(model);
    q = make_source_vec(model);
    x0 = copy(vel0_ex)

    Nx_pml = model.Nx + 2*model.pml_len;
    Ny_pml = model.Ny + 2*model.pml_len;
    n = Nx_pml*Ny_pml;
    optimizer = L_BFGS_B(Nx_pml*Ny_pml, 17);
    bounds = zeros(3, n)
    for i = 1:n
        bounds[1,i] = 2
        bounds[2,i] = lower
        bounds[3,i] = upper
    end

    println("Frequency content is: ", fre_vec')

    for iter_fre = 1:length(fre_vec)
        fre = fre_vec[iter_fre];
        println("==============================")
        println("Frequency is: ", fre, " Hz")

        # Solve data
        A = make_diff_op(vel_ex, model, fre);
        q = make_source_vec(model);
        u = A \ q;|
        d = Pr * u;

        # setup function
        # f(x) = compute_misfit(x, model, fre, d, q, Pr);
        # g(x) = compute_gradient(x, model, fre, d, q, Pr);
        # func(x) = f(x), g(x);
        func(x) = compute_gradient_lbfgsb(x, model, fre, d, q, Pr);

        # Optimization
        fout, xout = optimizer(func, x0, bounds, m=5, factr=1e12, pgtol=1e-10, iprint=1, maxfun=10000, maxiter=iterTime)

        # update
        x0 = copy(xout);
    end
    return x0
end


# using Optim
# # reference: http://julianlsolvers.github.io/Optim.jl/stable/#

# function lbfgs_fwi(vel_ex, vel0_ex, model, fre_vec; iterTime=10)

#     # Initialization
#     Pr = make_projection_op(model);
#     q = make_source_vec(model);
#     x0 = copy(vel0_ex)

#     println("Frequency content is: ", fre_vec')

#     for iter_fre = 1:length(fre_vec)
#         fre = fre_vec[iter_fre];
#         println("==============================")
#         println("Frequency is: ", fre, " Hz")

#         # Solve data
#         A = make_diff_op(vel_ex, model, fre);
#         q = make_source_vec(model);
#         u = A \ q;
#         d = Pr * u;

#         # setup function
#         f(x) = compute_misfit(x, model, fre, d, q, Pr)
#         g(x) = compute_gradient(x, model, fre, d, q, Pr)

#         # Optimization
#         results = optimize(f, g, x0, LBFGS(), Optim.Options(iterations=10, show_trace=true), inplace = false);

#         # update
#         x0 = Optim.minimizer(results);
#     end

#     return x0
# end