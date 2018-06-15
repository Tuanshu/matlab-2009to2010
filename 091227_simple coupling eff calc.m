clear all;
k=10^-4;                    %coupling constant (micorn^-1)
lambda_center=1.55;     %(micron)
lambda=1.50:0.0001:1.6;
delta_n=0.02;
delta_b=2*pi*delta_n*(1./lambda-1/lambda_center);
kL=pi/3;
L=kL/k;                 %(micron)
eff=(abs(k)^2./(abs(k)^2+delta_b.^2)).*sin(abs(k)*L.*(1+(delta_b/2/abs(k)).^2).^0.5).^2;
env=(abs(k)^2./(abs(k)^2+delta_b.^2));
FWHM=lambda(find(eff>0.5*max(eff), 1, 'last'))-lambda(find(eff>0.5*max(eff), 1, 'first'));
plot(lambda,eff,lambda,env);