clear all;
close all;
clc;

N_f=4096*4;
N_t=N_f*16;

%temp=importdata('J:\100721_SDOCT with OO\100721_set2(1ms)_R+S(interfered)-3 (5times averaged).txt');
temp=importdata('D:\spectrum\2.txt');
%temp=temp.data;
lambda=temp(:,1)*1e-9;  %m
S0=temp(:,2);
S0=10.^(S0/10);       % convert to linear from dbm
lambda=lambda(~isnan(lambda));
S0=S0(~isnan(lambda));
S0=S0/max(S0);
%S0=(S0-min(S0))/(max(S0)-min(S0));

c=3E8;              %m/sec
freq=c./lambda;     %Hz
d_f=max(freq)/(N_f-1);
fx=0:d_f:max(freq);
S=interp1(freq,S0,fx);
S(isnan(S))=0;

zeros(1:N_f)=0;
S_padded=[S zeros];             %with minus frequency, 2*N_f
%CS=real(fft(S_padded,N_t));     %with minus time
CS=real(fft(S_padded));     %with minus time
CS_normal=CS/max(abs(CS));
%d_t=1/(d_f*N_t);
d_t=1/(d_f*2*N_f);
%time=[-0.5*(N_t-1)*d_t:d_t:0.5*N_t*d_t];
time=[-(N_f-1)*d_t:d_t:N_f*d_t];
space=c*time;
CS_envelope=abs(hilbert(CS_normal));
plot(space,CS_normal,space,CS_envelope);
FWHM=abs(2*(space(1)-space(find(CS_envelope<0.5, 1, 'first'))));

%BW=lambda(find(S0>0.5,1,'last'))-lambda(find(S0>0.5,1,'first'));
%x_Res=x(find(CS_envelope>0.5,1,'last'))-x(find(CS_envelope>0.5,1,'first'));