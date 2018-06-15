clear all;
lambda=0.76;     %micron
k0=2*pi/lambda; %1/micron
a=20;           %core radius (micron)
dr=a/1000;
n1=1.77;        %of ti:sapphire
n2=1.46;        %of glass
V=k0*a*(n1^2-n2^2)^0.5;
ha=0:0.0001:100;
qa=(V^2-ha.^2).^0.5;
LHS=ha.*besselj(1,ha)./besselj(0,ha);
RHS=qa.*besselk(1,qa)./besselk(0,qa);
X=abs(LHS-RHS);
ha_0=ha(find(X<1,1,'first'));
%plot(ha,LHS,ha,RHS);
h_0=ha_0/a;         %h=(n1^2*k0^2-b^2)^0.5
b_0=(n1^2*k0^2-h_0.^2).^0.5;
neff_0=b_0/k0;

r=1:dr:a;
plot(r,besselj(0,h_0*r).^2);
FWHM=2*r(find(besselj(0,h_0*r).^2>0.5,1,'last'));
I=besselj(0,h_0*r).^2;
%h_0~=2.4/a;
%J0(hr=2.4)=0;
%J0(hr,r=a)=0;
total_power=sum(2*pi.*r.*I*dr);
integrated_power_ratio(1:length(r))=0;
for i=1:length(r)
integrated_power_ratio(i)=sum(2*pi.*r(1:i).*I(1:i)*dr)/total_power;
end
r_waist=r(find(integrated_power_ratio>(1-1/2.71828183),1,'first'));
plot(r,integrated_power_ratio);