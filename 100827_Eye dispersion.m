clear all;

lambda=0.6:0.001:1.1;
D_n=0.0512-0.1455*lambda +0.0961*lambda.^2;
plot(lambda,D_n);
delta_n=max(D_n)-min(D_n);
dispersion= (delta_n/(max(lambda)-min(lambda)))*1E-3;            %/nm