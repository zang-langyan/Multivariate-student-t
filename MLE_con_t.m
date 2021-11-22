function [muhat,scalehat,dfhat,fitted_c,f_x2,fig] = MLE_con_t(X,x1,span,useM)
% This function returns MLE parameters of X2 given X1 conditional
% distribution from the bivariate t distribution - f
% Input
%      - X is a n-by-2 matrix of the simulated X1 and X2 
%      - x1 is the condition vector that X2 expect to condition on
%      - span is the span to fit the estimate on
%      - useM is the method to be used, default 'mle', 
%      if useM = 1, then use titer
% Output
%      - muhat returns mle location of X2 given X1
%      - chat returns the simulated variables in bivariate t distribution
%      - dfhat returns n random numbers drawn from the bivariate t distribution
%      - fitted_c returns fitted mle distribution density of X2 given X1
%      - fig returns fitted distribution density of X2

if nargin < 4, useM = 0; end

X1 = X(:,1); X2 = X(:,2);

n = length(span);
k = length(x1);

muhat = zeros(k,1);
scalehat = zeros(k,1);
dfhat = zeros(k,1);
fitted_c = zeros(n,k);
for i = 1:k
    eps = 1e-4;
    ind = find(X1<(x1(i)+eps) & X1>(x1(i)-eps));
    if length(ind) < 1000
        eps2 = 1e-3;
        ind = find(X1<(x1(i)+eps2) & X1>(x1(i)-eps2));
        if length(ind) < 1000
           eps3 = 1e-2;
           ind = find(X1<(x1(i)+eps3) & X1>(x1(i)-eps3));
        end
    end

    f_x2 = X2(ind);
    if useM == 0
    phat = mle(f_x2,"distribution","tlocationscale");
    muhat(i) = phat(1);
    scalehat(i) = phat(2);
    dfhat(i) = phat(3);
    else
    [d_hat,mu_hat,c_hat,~] = titer(f_x2);
    muhat(i) = mu_hat;
    scalehat(i) = c_hat;
    dfhat(i) = d_hat;
    end

    fitted = pdf("tLocationScale",span,muhat(i),scalehat(i),dfhat(i));
    fitted_c(:,i) = fitted';

    if nargout > 4, figure;
    fig = plot(span,fitted);
    title('X2 density conditional on X1','interpreter','latex','fontsize',15);
    xlabel(strcat('Condition $X_1$ ',num2str(x1(i))),'interpreter','latex','fontsize',12) 
    ylabel('Density','interpreter','latex','fontsize',12) 
    end
end
end


