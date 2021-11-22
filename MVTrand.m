function [r,fig] = MVTrand(n,mu,scale,df,seed)
% This function returns random numbers of Multivariate t distribution
% Input
%      - n is the number of random numbers you want to simulate
%      - mu is the p-by-1 location vector of bivariate t distribution
%      - scale is the p-by-p scale matrix
%      - df is the degree of freedom
%      - seed is the random seed
% Output
%      only specify one output variable will only return random numbers
%      - r returns n-by-p random numbers drawn from the bivariate t distribution
%      - fig returns the histogram of the simulated Bivariate t
if nargin < 5, seed='shuffle'; end
rng(seed);
p = length(mu);

v=df; T=n; Y=gamrnd(v/2,1,[T 1]) / (v/2); G=1./Y;
Z = mvnrnd(zeros(1,p),scale,n);
r = mu'+sqrt(G).*Z;

if p < 3
    if nargout > 1
        figure;
        if n > 1000
            nbins = round(sqrt(n)/2)*ones(1,2);
            fig = histogram2(r(:,1),r(:,2),nbins,'FaceColor','flat');
        else
            fig = histogram2(r(:,1),r(:,2),'FaceColor','flat');
        end
        title('Simulated Bivariate T distribution');
        if max(r)>15, xlim([-15 15]); end
        if min(r)<-15, ylim([-15 15]); end
        x1=xlabel('X1');
        x2=ylabel('X2');
        zlabel('Density');
        set(x1,'Rotation',30);
        set(x2,'Rotation',-30);
    end
end
end