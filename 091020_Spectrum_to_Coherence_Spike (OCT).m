clear all;
close all;
clc;

N_fx=4096*4;
N_x=N_fx*16;

temp=importdata('D:\spectrum\2.txt');
%temp=temp.data;
lambda=temp(:,1)*1e-3;
S0=temp(:,2);
S0=10.^(S0/10);       % convert to linear from dbm
lambda=lambda(~isnan(lambda));
S0=S0(~isnan(lambda));
S0=S0/max(S0);

ilambda=1./lambda;
dfx=max(ilambda)/(N_fx-1);
fx=0:dfx:max(ilambda);
S=interp1(ilambda,S0,fx);
S(isnan(S))=0;

CS=real(fftshift(fft(S,N_x)));
CS_normal=CS/max(abs(CS));
x=1/2*(-1/dfx*(1/2-1/N_x):1/dfx/N_x:1/dfx*(1/2));
CS_envelope=abs(hilbert(CS_normal));
plot(x,CS_normal,x,CS_envelope);

BW=lambda(find(S0>0.5,1,'last'))-lambda(find(S0>0.5,1,'first'));
x_Res=x(find(CS_envelope>0.5,1,'last'))-x(find(CS_envelope>0.5,1,'first'));