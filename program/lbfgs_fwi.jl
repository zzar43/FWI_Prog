using Optim
# reference: http://julianlsolvers.github.io/Optim.jl/stable/#

function lbfgs_fwi(vel_ex, vel0_ex, model, fre_vec; iterTime=10)

    # Initialization
    Pr = make_projection_op(model);
    q = make_source_vec(model);
    x0 = copy(vel0_ex)

    println("Frequency content is: ", fre_vec')

    for iter_fre = 1:length(fre_vec)
        fre = fre_vec[iter_fre];
        println("==============================")
        println("Frequency is: ", fre, " Hz")

        # Solve data
        A = make_diff_op(vel_ex, model, fre);
        q = make_source_vec(model);
        u = A \ q;
        d = Pr * u;

        # setup function
        f(x) = compute_misfit(x, model, fre, d, q, Pr)
        g(x) = compute_gradient(x, model, fre, d, q, Pr)

        # Optimization
        results = optimize(f, g, x0, LBFGS(), Optim.Options(iterations=10, show_trace=true), inplace = false);

        # update
        x0 = Optim.minimizer(results);
    end

    return x0
end