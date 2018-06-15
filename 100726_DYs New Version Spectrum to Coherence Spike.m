
clear all;
close all;
clc;

N_fx=4096;          % specify number of data points in spatial frequency domain
dx_max=0.001;       % specify upper limit of data resolution

%% Import spectrum data from text file
% text file format:
% first column: wavelength (nm)
% second column: spectral intensity (linear scale)

temp=importdata('D:\spectrum\2.txt');  % Enter the file path
if isfield(temp,'data')
    temp=temp.data;
end
lambda=temp(:,1)*1e-3;  % convert from nm to um
S0=temp(:,2);
S0=10.^(S0/10);         % convert from log scale to linear
S0=S0(~isnan(lambda));
lambda=lambda(~isnan(lambda));
S0=S0(lambda>0.550);          % discard missleading data
lambda=lambda(lambda>0.550);  % discard missleading data

%% Convert from wavelength domain to spatial frequency domain
%Increase number of data points for accuracy
lambda2=min(lambda):(max(lambda)-min(lambda))/(N_fx-1):max(lambda);
S0=interp1(lambda,S0,lambda2,'spline');
lambda=lambda2;
S0=S0/max(S0);          % Normalize
%plot(lambda,S0)    % confirm spectrum

ilambda=1./lambda;
dfx=max(ilambda)/(N_fx-1);
fx=0:dfx:max(ilambda);
S=interp1(ilambda,S0.*lambda.^2,fx);
S(isnan(S))=0;

%% FFT and Hilbert transform
% Calculate suitable N_x for specified data resolution
N_x=2^(ceil(log2(N_fx/dx_max/max(fx))));

CS=real(fftshift(fft(S,N_x)));
CS_normal=CS/max(abs(CS));
dx=1/dfx/N_x;
x=1/2*(-1/dfx/2:dx:1/dfx/2-dx);
CS_envelope=abs(hilbert(CS_normal));

%% Calculate 3-dB bandwidth and axial resolution
%BW=lambda(find(S0>0.5,1,'last'))-lambda(find(S0>0.5,1,'first'));
%x_Res=x(find(CS_envelope>0.5,1,'last'))-x(find(CS_envelope>0.5,1,'first'));

%% Crop interference data into 100-micron range
%CS_envelope=interp1(x,CS_envelope,-50:0.01:50);
%CS_normal=interp1(x,CS_normal,-50:0.01:50);
%x=-50:0.01:50;
plot(x,CS_normal,x,CS_envelope); 

%% Calculate crosstalk
%x_ct=[-3:-1,1:3]*x_Res;
%CRS_TLK=20*log10(interp1(x,CS_envelope,x_ct));
%CS_envelope_dB=20*log10(CS_envelope);
%figure
%plot(x,CS_envelope_dB,x_ct,CRS_TLK,'ro');
