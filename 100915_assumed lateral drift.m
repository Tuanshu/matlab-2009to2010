clear all;
close all;
clc;
% note:interpft
N_z=4*4096;        %for iFFT in the zero padding process
N_f=4*4096;

%%grating related (the center wavelength is the NOT the "blazed" wavelength)
grating_pitch=1000000/600;      %nm
center_wavelength=762;       %nm
incidence_angle=asin(center_wavelength/2/grating_pitch);     %���]�O�Ӧ����פJ�g, rad

assumed_BW=178;

short_coef=tan(asin((center_wavelength/2-assumed_BW/2)/grating_pitch)-incidence_angle);         %tan!
long_coef=tan(asin((center_wavelength/2+assumed_BW/2)/grating_pitch)-incidence_angle);          %tan!

pixel=[1:4096]';  %1~4096

%Q=178/(index_long-index_short);
%lambda=((-(pixel-index_peak)*Q+center_wavelength));   %�e�����t���ܭ��n! �v�Tfreq domain�O�_����guassian (���L������t�@�ӭt��dispersion��broaden�����n���]�|�ܤj? �]��carrier in lambda domain��chirp, ���bfreq domain���Ӥ��Ӧ�)

%Q=(long_coef-short_coef)/(index_long-index_short);      %roughly flens/pixel size
%lambda=grating_pitch*sin((asin((index_peak-pixel)*Q)+incidence_angle))+center_wavelength/2;

Q=(long_coef-short_coef)/(2000);      %roughly flens/pixel size
lambda=grating_pitch*sin((atan((2000-pixel)*Q)+incidence_angle))+center_wavelength/2;   %atan!
pixelr=(2048-pixel)/100;
plot(lambda,pixelr);
