function f_S = S_weighted_t(mu,scale,nu,w,span)
% This function returns non-zero weighted sum of the univariate margins
% distribution from a p dimentional multivariate t distribution 
% Input
%      - mu is the p-by-1 location vector of bivariate t distribution
%      - scale is the p-by-p scale matrix
%      - nu is the degree of freedom
%      - w is the p-by-1 weight vector specifying each dimention weight
%      - span is the span you want to calculate the density(e.g.,-2:0.25:2)
% Output
%      - f_S returns the weighted sum of p univariate margins distribution

mu_S = w'*mu; kappa = sqrt(w'*scale*w);
n = length(span); f_S = zeros(n,1);

for i = 1:n
    %inversion formulae for density function
    f_S(i) = 1/(2*pi)*...
        quadgk(@(t)exp(-1i.*t.*span(i)).*... %exp(-it*s)
        exp(1i.*t.*mu_S).*... %exp(i*t*mu_S)
        besselk(nu/2,nu^0.5.*abs(kappa*t)).*... %K_v/2(v^0.5*|kappa*t|) 
        (nu^0.5.*abs(kappa*t)).^(nu/2)/... %(v^0.5*|kappa*t|)^(v/2)
        (gamma(nu/2)*2^(nu/2-1)),... %Gamma(v/2)*2^(v/2-1)
        -inf,inf);   
end
end
