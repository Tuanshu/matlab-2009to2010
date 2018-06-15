clear all;
close all;
clc;

N_f=4096*8;
N_t=N_f*16;


data_array=importdata('D:\101206_mirror under lens.txt');




%% FD data calculation

%ref=data_array(real_max_index,:);

%% Calibrate by asign 444 and 888 pixel
grating_pitch=1000000/600;      %nm
%center_wavelength=760;       %nm

wavelength_1=444;
wavelength_2=888;

index_1=3425;
index_2=1825;



spe_pixel=(wavelength_2-wavelength_1)/(index_2-index_1);   %spectral res of CCD (nm), can be minus, in fact its minus in current system
pixel=[1:4096]';  %1~4096
lambda_o=(pixel-index_1)*spe_pixel+wavelength_1;  %only used in peak finding

[dc_max_value dc_max_index]=max(data_array);
found_peak_wavelength=lambda_o(dc_max_index);
incidence_angle=asin(found_peak_wavelength/2/grating_pitch);     %���]�O�Ӧ����פJ�g,rad

coef_1=tan(asin((wavelength_1)/grating_pitch)-incidence_angle);         %tan!
coef_2=tan(asin((wavelength_2)/grating_pitch)-incidence_angle);

Q=-(coef_1-coef_2)/(index_1-index_2);      %roughly flens/pixel size
%%tan!
lambda=grating_pitch*sin((atan((dc_max_index-pixel)*Q)+incidence_angle))+found_peak_wavelength/2;   %atan!

%assumed_BW=137.28;

%short_coef=tan(asin((center_wavelength/2-assumed_BW/2)/grating_pitch)-incidence_angle);         %tan!
%long_coef=tan(asin((center_wavelength/2+assumed_BW/2)/grating_pitch)-incidence_angle);          %tan!



%[peak index_peak]=max(ref);
%index_long=find(ref>0.5*max(ref),1,'last');
%index_short=find(ref>0.5*max(ref),1,'first');

%Q=178/(index_long-index_short);
%lambda=((-(pixel-index_peak)*Q+found_peak_wavelength));   %�e�����t���ܭ��n! �v�Tfreq domain�O�_����guassian (���L������t�@�ӭt��dispersion��broaden�����n���]�|�ܤj? �]��carrier in lambda domain��chirp, ���bfreq domain���Ӥ��Ӧ�)

%Q=(long_coef-short_coef)/(index_long-index_short);      %roughly flens/pixel size
%lambda=grating_pitch*sin((asin((index_peak-pixel)*Q)+incidence_angle))+center_wavelength/2;

%Q=(long_coef-short_coef)/(index_long-index_short);      %roughly flens/pixel size
%lambda=grating_pitch*sin((atan((index_peak-pixel)*Q)+incidence_angle))+center_wavelength/2;   %atan!

%�ϥηs���t��k��, �̵uwavelength�q166micron -> 10micron, �i��y���I�Ƥ���, �ҥH���h�����n��data

%ref=ref(lambda>400);
%lambda2=lambda2(lambda>300);
lambda_2=lambda(lambda>300);

%% x-axis
c=3E8;                     %m/sec
freq=c./(lambda_2*1E-9);     %Hz orignal freq array
d_f=max(freq)/(N_f-1);
fx=0:d_f:max(freq);        %freq after interpolation
d_t=1/(d_f*N_t);
%d_t=1/(d_f*2*N_f);
time=[-0.5*(N_t-1)*d_t:d_t:0.5*N_t*d_t]'/2;%/2�O�]���@�Ӥ@�^
%time=[-(N_f-1)*d_t:d_t:N_f*d_t]/2;  %/2�O�]���@�Ӥ@�^
length_space_FD=round(length(time)/25);
space_FD=c*time(1:length_space_FD); % �ȩw�u��array�e1/100 (�j���]��100 micron���k�a)
space_FD=(space_FD-space_FD(1))*1E6;                     % shift to zero & m to micron
%plot(lambda,inter,lambda2,inter);

%% data related
S0=data_array-mean(data_array:100);
%S0=S0-mean(S0(lambda>300)); %background light substration
S0=S0(lambda>300);
%S0=inter;
%plot(lambda,S0/max(S0),lambda,inter/max(inter),lambda,ref/max(ref));
%S0=S0/max(S0);        %����normal���Ӥ]�S�t?
S=interp1(freq,S0,fx);
S(isnan(S))=0;
CS=fft(S,N_t)';     %with minus time
%CS=real(fft(S_padded));     %with minus time
CS_normal=CS/max(abs(CS));
CS_envelope=abs(CS(1:length_space_FD));                                        %data 1 to record (need decimate? maybe cut is better)
%% 

%plot(space,CS_normal,space,CS_envelope);

%% After the acquire of envelope


%%  DC PSF FWHM
value_DC_FD=CS_envelope;
FWHM_DC_FD=2*(space_FD(find(CS_envelope<0.5*value_DC_FD, 1, 'first'))-space_FD(1));

%%  Interfered PSF FWHM
space_min_for_inter_peak_FD=10;  %10fixed
space_min_for_inter_peak_index_FD=find(space_FD>space_min_for_inter_peak_FD, 1, 'first');
[inter_peakvalue_FD inter_peakindex_FD]=max(CS_envelope(space_min_for_inter_peak_index_FD:end)); %�`�N! ���B��inter_peakindex_FD�|���, �]��CS_envelope�qspace_min_for_inter_peak_index_FD�}�l�ӫD0�}�l
inter_position_FD=space_FD(inter_peakindex_FD+space_min_for_inter_peak_index_FD-1);
FWHM_FD_right=space_FD(find(CS_envelope(space_min_for_inter_peak_index_FD:end)>0.5*inter_peakvalue_FD, 1, 'last'));
FWHM_FD_left=space_FD(find(CS_envelope(space_min_for_inter_peak_index_FD:end)>0.5*inter_peakvalue_FD, 1, 'first'));
FWHM_inter_FD=FWHM_FD_right-FWHM_FD_left;
FWHM_FD=FWHM_inter_FD;
interference_efficiency=inter_peakvalue_FD/value_DC_FD*2;
%plot(space_FD,CS_normal,space_FD,CS_envelope);
%dlmwrite('Interferogram.txt',M,'delimiter','\t','newline','pc');
%BW=lambda(find(S0>0.5,1,'last'))-lambda(find(S0>0.5,1,'first'));
%x_Res=x(find(CS_envelope>0.5,1,'last'))-x(find(CS_envelope>0.5,1,'first'))
%;
%plot(space_TD,inter_position_FD,space_TD,interference_efficiency,space_TD,FWHM_FD);
%plot(space_TD,inter_position_FD,space_TD,FWHM_FD);
interference_efficiency=interference_efficiency';
inter_position_FD=inter_position_FD';

%imagesc(10*log10(CS_envelope./maxs)); figure(gcf)
%imagesc(10*log10(CS_envelope(2500:3500,16:116)./max(max(CS_envelope(2500:3500,16:116),[],1))),'xdata',[0:5:500],'ydata',[0:0.009375850226843:0.009375850226843*1000]); figure(gcf);
%dlmwrite('power_envelope.txt',power_envelope,'delimiter','\t','newline','pc');
%dlmwrite('power.txt',power,'delimiter','\t','newline','pc');
%dlmwrite('space_TD.txt',space_TD,'delimiter','\t','newline','pc');
%dlmwrite('CS_envelope.txt',CS_envelope,'delimiter','\t','newline','pc');
%dlmwrite('interference_efficiency.txt',interference_efficiency','delimiter','\t','newline','pc');
%dlmwrite('dispersion_expansion_ratio.txt',dispersion_expansion_ratio','delimiter','\t','newline','pc');
%dlmwrite('inter_position_FD.txt',inter_position_FD','delimiter','\t','newline','pc');