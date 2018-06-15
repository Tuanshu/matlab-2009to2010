clear all;
close all;
clc;

N=100;
Range=5;







temp=importdata('D:\spectrum\2.txt');
%temp=temp.data;
lambda=temp(:,1)*1e-3;
S0=temp(:,2);
S0=10.^(S0/10); 
lambda=lambda(~isnan(lambda));
S=S0(~isnan(lambda));
S=S.^0.5;   %for field

%% Sellmeier

t=25;
		    fe=(t-24.5)*(t+570.82);
	        c1=5.35583      ;                          
		    c2=0.100473;
		    c3=0.20692;
		    c4=100;
		    c5=11.34927;
		    c6=-1.5334e-2;
		    d1=4.629e-7;
		    d2=3.826e-8;
		    d3=-8.9e-9;
		    d4=2.657e-5		 ;

			n=(c1+d1*fe+(c2+d2*fe)./(lambda.^2-(c3+d3*fe)^2)+(c4+d4*fe)./(lambda.^2-(c5)^2)+c6*lambda.^2).^0.5;
            n=2.1448;
%%
k=2*pi.*n./lambda;
delta_k=zeros(length(k),length(k));
S2=zeros(length(k),length(k));
for j=1:length(lambda)
    delta_k(j,:)=k-k(j);
    S2(j,:)=S*S(j);
end
X=-0.5*Range:Range/N:0.5*Range;
Q=zeros(1,length(X));
for j=1:length(X)
Q(j)=sum(sum(S2.*real(exp(i.*delta_k.*X(j)))));
end
Q=Q/max(Q);
%Q_en=abs(hilbert(Q));
%FWHM=X(find(Q_en>0.5, 1, 'last'))-X(find(Q_en>0.5, 1, 'first'));
plot(X,Q);
