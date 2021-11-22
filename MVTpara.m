function [muhat,Sigmahat,nuhat] = MVTpara(T,knownpara1,knownpara2)
% This function returns the approximate parameters from a n-by-p sample
% assuming to be a p dimensional multivariate t distribution, by method
% Batch Approximation Algorithm attributed to Aeschliman, Park & Kak(2010)
% Input
%      - T is a n-by-p sample matrix with n observations and p variables
%      - knownpara1 and knownpara2 (Optional) are the known parameters
%        (mu,sacle or degree of freedom) which force the estimate to be the
%        known true ones. See below for more information about Options
% Output
%      - muhat returns the approximate p-by-1 location vector
%      - Sigmahat returns the approximate p-by-p scale matrix
%      - nuhat returns the approximate degree of freedom
% Optional Inputs (knownpara1 & knownpara2)
%      - known mu must be p-by-1 vector
%      - known scale matrix must be p-by-p matrix, when you only know the
%        scale matrix is a scaled diagonal matrix, specify the 'knownpara'
%        Option to be zeros(p)
%      - known degree of freedom must be a 1-by-1 scaler
[n,p] = size(T);

%% No parameter is known
if nargin == 1
    c = trimmean(T,75); muhat = c'; %excluding 75% outliers
    
    z = log(vecnorm(T-c,2,2).^2);
    z_bar = mean(z);
    b = mean((z-z_bar).^2) - psi(1,p/2);
    nuhat = (1+sqrt(1+4*b))/b;

    alphahat = exp(z_bar-log(nuhat)+psi(nuhat/2)-psi(p/2));
    beta = 2*log2(p)/(nuhat^2+log2(p));
    S = zeros(p);
    for i = 1:n
        Si = (T(i,:)-c)'*(T(i,:)-c)/norm(T(i,:)-c)^beta;
        S = S + Si;
    end
    S_bar = S/n;
    Sigmahat = (alphahat*p/trace(S_bar))*S_bar;

    %% One parameter is known
elseif nargin == 2
    if isequal(size(knownpara1),[p,1]) % mu known
        c = knownpara1'; muhat = knownpara1;

        z = log(vecnorm(T-c,2,2).^2);
        z_bar = mean(z);
        b = mean((z-z_bar).^2) - psi(1,p/2);
        nuhat = (1+sqrt(1+4*b))/b;

        alphahat = exp(z_bar-log(nuhat)+psi(nuhat/2)-psi(p/2));
        beta = 2*log2(p)/(nuhat^2+log2(p));
        S = zeros(p);
        for i = 1:n
            Si = (T(i,:)-c)'*(T(i,:)-c)/norm(T(i,:)-c)^beta;
            S = S + Si;
        end
        S_bar = S/n;
        Sigmahat = (alphahat*p/trace(S_bar))*S_bar;

    elseif isequal(size(knownpara1),[p,p]) % scale known
        c = trimmean(T,75); muhat = c';

        z = log(vecnorm(T-c,2,2).^2);
        z_bar = mean(z);
        b = mean((z-z_bar).^2) - psi(1,p/2);
        nuhat = (1+sqrt(1+4*b))/b;

        if isequal(knownpara1,zeros(p)) % only know the scale matrix is a scaled I
            alphahat = exp(z_bar-log(nuhat)+psi(nuhat/2)-psi(p/2));
            Sigmahat = alphahat*eye(p);
        else % the whole scale matrix is known
            Sigmahat = knownpara1;
        end


    else % degree of freedom known
        c = trimmean(T,75); muhat = c';

        z = log(vecnorm(T-c,2,2).^2);
        z_bar = mean(z);
        nuhat = knownpara1;

        alphahat = exp(z_bar-log(nuhat)+psi(nuhat/2)-psi(p/2));
        beta = 2*log2(p)/(nuhat^2+log2(p));
        S = zeros(p);
        for i = 1:n
            Si = (T(i,:)-c)'*(T(i,:)-c)/norm(T(i,:)-c)^beta;
            S = S + Si;
        end
        S_bar = S/n;
        Sigmahat = (alphahat*p/trace(S_bar))*S_bar;
    end

    %% Two parameters are known
else
    if or(isequal(size(knownpara1),[p,1]) && isequal(size(knownpara2),[p,p]),...
            isequal(size(knownpara1),[p,p]) && isequal(size(knownpara2),[p,1]))
        % mu and scale matrix are known
        if isequal(size(knownpara1),[p,1]) %knownpara1 is mu
            c = knownpara1'; muhat = knownpara1;

            z = log(vecnorm(T-c,2,2).^2);
            z_bar = mean(z);
            b = mean((z-z_bar).^2) - psi(1,p/2);
            nuhat = (1+sqrt(1+4*b))/b;

            if isequal(knownpara2,zeros(p)) % only know the scale matrix is a scaled I
                alphahat = exp(z_bar-log(nuhat)+psi(nuhat/2)-psi(p/2));
                Sigmahat = alphahat*eye(p);
            else % the whole scale matrix is known
                Sigmahat = knownpara2;
            end

        else % knownpara1 is scale
            c = knownpara2'; muhat = knownpara2;

            z = log(vecnorm(T-c,2,2).^2);
            z_bar = mean(z);
            b = mean((z-z_bar).^2) - psi(1,p/2);
            nuhat = (1+sqrt(1+4*b))/b;

            if isequal(knownpara1,zeros(p)) % only know the scale matrix is a scaled I
                alphahat = exp(z_bar-log(nuhat)+psi(nuhat/2)-psi(p/2));
                Sigmahat = alphahat*eye(p);
            else % the whole scale matrix is known
                Sigmahat = knownpara1;
            end
        end

    elseif or(isequal(size(knownpara1),[p,1]) && isequal(size(knownpara2),[1,1]),...
            isequal(size(knownpara1),[1,1]) && isequal(size(knownpara2),[p,1]))
        % mu and degree of freedom are known
        if isequal(size(knownpara1),[p,1]) %knownpara1 is mu
            c = knownpara1'; muhat = knownpara1;

            z = log(vecnorm(T-c,2,2).^2);
            z_bar = mean(z);
            nuhat = knownpara2;

            alphahat = exp(z_bar-log(nuhat)+psi(nuhat/2)-psi(p/2));
            beta = 2*log2(p)/(nuhat^2+log2(p));
            S = zeros(p);
            for i = 1:n
                Si = (T(i,:)-c)'*(T(i,:)-c)/norm(T(i,:)-c)^beta;
                S = S + Si;
            end
            S_bar = S/n;
            Sigmahat = (alphahat*p/trace(S_bar))*S_bar;
        else %knownpara1 is degree of freedom
            c = knownpara2'; muhat = knownpara2;

            z = log(vecnorm(T-c,2,2).^2);
            z_bar = mean(z);
            nuhat = knownpara1;

            alphahat = exp(z_bar-log(nuhat)+psi(nuhat/2)-psi(p/2));
            beta = 2*log2(p)/(nuhat^2+log2(p));
            S = zeros(p);
            for i = 1:n
                Si = (T(i,:)-c)'*(T(i,:)-c)/norm(T(i,:)-c)^beta;
                S = S + Si;
            end
            S_bar = S/n;
            Sigmahat = (alphahat*p/trace(S_bar))*S_bar;
        end

    else
        % scale matrix and degree of freedom are known
        if isequal(size(knownpara1),[p,p]) %knownpara1 is scale
            c = trimmean(T,75); muhat = c';

            z = log(vecnorm(T-c,2,2).^2);
            z_bar = mean(z);
            nuhat = knownpara2;

            if isequal(knownpara1,zeros(p)) % only know the scale matrix is a scaled I
                alphahat = exp(z_bar-log(nuhat)+psi(nuhat/2)-psi(p/2));
                Sigmahat = alphahat*eye(p);
            else % the whole scale matrix is known
                Sigmahat = knownpara1;
            end

        else %knownpara1 is degree of freedom
            c = trimmean(T,75); muhat = c';

            z = log(vecnorm(T-c,2,2).^2);
            z_bar = mean(z);
            nuhat = knownpara1;

            if isequal(knownpara2,zeros(p)) % only know the scale matrix is a scaled I
                alphahat = exp(z_bar-log(nuhat)+psi(nuhat/2)-psi(p/2));
                Sigmahat = alphahat*eye(p);
            else % the whole scale matrix is known
                Sigmahat = knownpara2;
            end
        end
    end
end
end









