clear all;

z=[0:0.01:2000]';       %micron
CCD_lambda=0.350/4048;   %micron
grating_lens_lambda=CCD_lambda;    %micron
center_lambda=0.76; %micron
R=log10((sinc(4*(CCD_lambda/center_lambda^2).*z).^2).*exp(-8/log(2)*((grating_lens_lambda/center_lambda^2).*z).^2));
plot(z,R);