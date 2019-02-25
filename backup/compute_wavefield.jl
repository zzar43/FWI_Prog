function compute_wavefield(vel, fre, model, q, d, Pr)

    A0 = make_diff_op(vel_ex, model, fre);
    u0 = A0 \ q;
    q_ad = -1 .* Pr' * (Pr*u0 - d);
    p = conj(A0) \ q_ad;
    g = real((2*pi*fre)^2 .* conj(u0) .* p);
    g = sum(g,dims=2);

    phi = 0.5 * norm(Pr*u0-d)^2;
end