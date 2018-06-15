clear all;
close all;
clc;

N_fx=4096*4;
N_x=N_fx*16;

temp=importdata('D:\100723_TD data\100723_TDF3_fatch40micronsec_N5000_S25kHz_400mA_40dB.txt');
%temp=temp.data;
time=temp(:,1); %ms
S0=temp(:,2);
time=time(~isnan(time));
time=time-min(time);
S0=S0(~isnan(time));
S0=(S0-min(S0))/(max(S0)-min(S0));

N=length(time);
D_t=max(time);
d_t=D_t/N;
CS=real(fftshift(fft(S0)));
%CS=real(fftshift(fft(S0,N_x)));
%CS_normal=CS/max(abs(CS));
%x=1/2*(-1/dfx*(1/2-1/N_x):1/dfx/N_x:1/dfx*(1/2));
%CS_envelope=abs(hilbert(CS_normal));
%plot(x,CS_normal,x,CS_envelope);
plot(CS);

%BW=lambda(find(S0>0.5,1,'last'))-lambda(find(S0>0.5,1,'first'));
%x_Res=x(find(CS_envelope>0.5,1,'last'))-x(find(CS_envelope>0.5,1,'first'));