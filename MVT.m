function f_T = MVT(X,mu,scale,df)
% Inputs  
%        - X is a n-by-2 matrix, where n is the points expected to be
%          calculated
%        - mu is the 2-by-1 location vector
%        - scale is the 2-by-2 scale matrix
%        - df is the degree of freedom
% Output
%        - returns the bivariate t distribution density
[n,p] = size(X);
G = gamma((df+p)/2)/(gamma(df/2)*(df*pi)^(p/2)*sqrt(det(scale)));
M = zeros(n);
for i = 1:n
    for j =1:n
        C = [X(i,1);X(j,2)];
        M(i,j) = (C-mu)'/scale*(C-mu);
    end
end
f_T = G*(1+M/df).^(-(df+p)/2);
end

