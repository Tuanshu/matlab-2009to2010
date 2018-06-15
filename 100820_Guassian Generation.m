
clear all;
close all;
clc;

N_fx=4096;          % specify number of data points in spatial frequency domain
dx_max=0.001;       % specify upper limit of data resolution

%% Import spectrum data from text file
% text file format:
% first column: wavelength (nm)
% second column: spectral intensity (linear scale)
%temp=importdata('D:\100820_test output.txt');  % Enter the file path
temp=importdata('D:\1W 12milisecond Nd filter2+0.5 before.txt');  % Enter the file path
if isfield(temp,'data')
    temp=temp.data;
end
lambda=temp(:,1);  % convert from nm to um
c=3E8;
center=c/762;
fwhm=c/673-c/851;
frequency=c./lambda;
Sfalse=gaussmf(frequency,[fwhm/(2^0.5)/2/(log(2)^0.5) center]);

M=[lambda Sfalse];