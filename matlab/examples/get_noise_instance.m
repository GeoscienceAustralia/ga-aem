function noise = get_noise_instance(d,mnp,an)
    % d   = [n x 1] response data [n x 1]
    % mnp = scalar relative noise percentage
    % an  = [n x 1] additive noise in units of d
    n=length(d);
    a_err = an .* randn(n,1);
    m_err = 0.01 * mnp * abs(d) .* randn(n,1);
    noise = a_err + m_err;
    % Add the noises
    % We are adding instances of noise not the population variances in
    % quadrature
end