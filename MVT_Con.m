function [para,f_c]=MVT_Con(x1,x2,mu,scale,df)
% This function returns the Conditional distribution of X2 given X1 
% X1 and X2 are parts of Multivariate t distribution with 2 dimension 
% Input
%      - x1 is the vector that x2 conditioned on 
%      - x2 the span expected to calculate the conditional distribution
%      - mu is the location vector
%      - scale is the scale matrix
%      - df is degree of freedom
% Output
%      - para returns the conditional parameters of X2
%      - f_c returns the conditional distribution of X2
nc = length(x1);
nx = length(x2);
% location
mu_c = mu(2)+scale(2,1)/scale(1,1)*(x1-mu(1));
% scale
d1 = (x1-mu(1))./scale(1,1).*(x1-mu(1));
fac = (df+d1)/(df+1);
sigma_c = sqrt(fac*(scale(2,2)-scale(2,1)/scale(1,1)*scale(1,2)));
% degree of freedom
df_c = (df+1)*ones(nc,1);
para = [mu_c,sigma_c,df_c]; %conditional on x1 in row
f_c = zeros(nx,nc);
for i = 1:nc
    f_c(:,i) = (pdf('tLocationScale',x2,mu_c(i),sigma_c(i),df_c(i)))';
end
end
