mu = [0;0];
scale = [1,0;0,1];
df = 4;
n = 1e6;
r = MVTrand(n,mu,scale,df,42);
[muhat,Sigmahat,nuhat] = MVTpara(r);

error = [norm(mu-muhat),norm(scale-Sigmahat),norm(df-nuhat)];
%              |                  |                |
%error--->     mu               scale              df

para = categorical({'mu','scale','degree of freedom'});
para = reordercats(para,{'mu','scale','degree of freedom'});
b = bar(para,error,'FaceColor','flat');
b.CData(1,:) = [0.4660 0.6740 0.1880];
b.CData(2,:) = [0.9290 0.6940 0.1250];
b.CData(3,:) = [0.8500 0.3250 0.0980];
ylabel('error')