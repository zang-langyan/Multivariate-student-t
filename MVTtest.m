% X1=-2:0.1:2;X2=-1:0.1:3;
% X=[X1',X2'];
% T1=MVT(X,[0;0],[1,0;0,1],4);
% T2=MVT(X,[1;2],[1,0.3;0.3,1],4);
% figure;
% mesh(X2,X1,T1);
% figure; 
% mesh(X2,X1,T2);

%specify the conditions and Bivariate t parameters
con = (-0.5:0.25:0.5)';
span = (-2:0.05:2)';
mu = [0;0];
scale = [1,0.6;0.6,1];
df = 5;

%compute the conditional parameters
[para,f_t]=MVT_Con(con,span,mu,scale,df);

%simulate and estimate conditional parameters
n = 1e7;
[r,plot_sim] = MVTrand(n,mu,scale,df,42);
[muhat,scalehat,dfhat,fitted_c] = MLE_con_t(r,con,span,1);

figure;
tiledlayout(1,2);
%Therotical
axm1 = nexttile;
mesh(axm1,con,span,f_t)
title(axm1,'Theoretical X2 conditional on X1');
x1=xlabel(axm1,'X1');   
x2=ylabel(axm1,'X2');        
zlabel(axm1,'Density');        
set(x1,'Rotation',30);    
set(x2,'Rotation',-30);
%mle fitted
axm2 = nexttile;
mesh(axm2,con,span,fitted_c)
title(axm2,'MLE Fitted X2 conditional on X1');   
x1=xlabel(axm2,'X1');   
x2=ylabel(axm2,'X2');        
zlabel(axm2,'Density');        
set(x1,'Rotation',30);    
set(x2,'Rotation',-30);  

figure;
tiledlayout(1,3);
ax1 = nexttile;
plot(ax1,con,para(:,1),con,muhat,'-o')
title(ax1,'Theoretical versus Estimate $\mu$','interpreter','latex','fontsize',15);
xlabel(ax1,'Condition $X_1$','interpreter','latex','fontsize',12) 
ylabel(ax1,'$\mu$','interpreter','latex','fontsize',12) 

ax2 = nexttile;
plot(ax2,con,para(:,2),con,scalehat,'-o')
title(ax2,'Theoretical versus Estimate $\sigma$','interpreter','latex','fontsize',15);
xlabel(ax2,'Condition $X_1$','interpreter','latex','fontsize',12) 
ylabel(ax2,'$\sigma$','interpreter','latex','fontsize',12) 

ax3 = nexttile;
plot(ax3,con,para(:,3),con,dfhat,'-o')
title(ax3,'Theoretical versus Estimate $\nu$','interpreter','latex','fontsize',15);
xlabel(ax3,'Condition $X_1$','interpreter','latex','fontsize',12) 
ylabel(ax3,'$\nu$','interpreter','latex','fontsize',12) 

lgd = legend(ax3,'Theoretical','MLE','Orientation','horizontal');
lgd.Layout.Tile = 'south';