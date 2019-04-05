function e = emgfxn_lu(x, a, x0, mu, sigma, B)
if isstruct(a)
    x0 = a.x_0;
    mu = a.mu_x;
    sigma = a.sigma_x;
    B = a.B;
    a = a.a;
end

e = a/(2 * x0) .* exp( mu / x0 + sigma.^2 / (2*x0.^2) - x/x0 )...
    .* erfc( -1/sqrt(2) * ((x-mu)/sigma - sigma/x0) ) + B;
e(isnan(e)) = Inf;
end