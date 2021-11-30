mu = [0.8;-0.5]; scale = [3,0.2;0.2,1]; df = 5; span = -2:0.01:3;
r = MVTrand(2e6,mu,scale,df,42); % use random numbers as assets returns 

w = [0.9;0.1]; % weights of each assest 
% sum weighted mu is 0.67
%% kernel density of portfolio returns
P = r*w; %weighted portfolio return
[kernel,~] = ksdensity(P,span);
figure;
plot(span,kernel,'--','LineWidth',2)
hold on

%% weighted sum of univariate margins density from assumed multivariate distribution(mu,scale,df)
f_S = S_weighted_t(mu,scale,df,w,span);
plot(span,f_S)
title({'Weighted sum of univariate margins density';'kernel versus true'},...
    'interpreter','latex','fontsize',15)
xlabel('Portfolio return','interpreter','latex','fontsize',12) 
ylabel('Density','interpreter','latex','fontsize',12) 
legend({'kernel',...
    'true'},...
    'Location','northwest',...
    'interpreter','latex',...
    'fontsize',12)
hold off