clear all;
close all;
clc;



data=importdata('D:\100913_2 arm but no interference (2048amp)_1kHz.txt');
data_array=data(1:500,2:end);
time_record=data(1:500,1);


N_f=4096*4;
%N_t=N_f;

%% ref&grating related (the center wavelength is the NOT the "blazed" wavelength)
grating_pitch=1000000/600;      %nm
center_wavelength=762;       %nm
incidence_angle=asin(center_wavelength/2/grating_pitch);     %假設是照此角度入射, rad

assumed_BW=178;

short_coef=tan(asin((center_wavelength/2-assumed_BW/2)/grating_pitch)-incidence_angle);         %tan!
long_coef=tan(asin((center_wavelength/2+assumed_BW/2)/grating_pitch)-incidence_angle);          %tan!

min_wavelength=550;     %nm
max_wavelength=1050;    %nm

pixel=[1:4096]';  %1~4096
[peak index_peak]=max(data_array(1,:));
index_long=find(data_array(1,:)>0.5*max(data_array(1,:)),1,'last');
index_short=find(data_array(1,:)>0.5*max(data_array(1,:)),1,'first');
%Q=178/(index_long-index_short);
%lambda=((-(pixel-index_peak)*Q+center_wavelength));   %前面的負號很重要! 影響freq domain是否接近guassian (不過為什麼差一個負號dispersion的broaden效應好像也會變大? 因為carrier in lambda domain有chirp, 但在freq domain應該不太有)

%Q=(long_coef-short_coef)/(index_long-index_short);      %roughly flens/pixel size
%lambda=grating_pitch*sin((asin((index_peak-pixel)*Q)+incidence_angle))+center_wavelength/2;

Q=(long_coef-short_coef)/(index_long-index_short);      %roughly flens/pixel size
lambda=grating_pitch*sin((atan((index_peak-pixel)*Q)+incidence_angle))+center_wavelength/2;   %atan!

%使用新的演算法後, 最短wavelength從166micron -> 10micron, 可能造成點數不夠, 所以先去掉不要的data

max_wavelength_index=find(lambda<max_wavelength,1,'first');
min_wavelength_index=find(lambda>min_wavelength,1,'last');
lambda=lambda(max_wavelength_index:min_wavelength_index);
data_array=data_array(:,max_wavelength_index:min_wavelength_index);
%lambda2=lambda2(lambda>300);

%% x-axis
c=3E8;                     %m/sec
freq=c./(lambda*1E-9);     %Hz orignal freq array
d_f=max(freq)/(N_f-1);
fx=0:d_f:max(freq);        %freq after interpolation


%% data related
[M N]=size(data_array);
S(1:length(fx),1:M)=0;
inter_position_FD(1:M)=0;
dispersion_expansion_ratio(1:M)=0;
interference_efficiency(1:M)=0;
%% SD summed in x-domain
for j=1:M
S0=interp1(freq,data_array(j,:),fx);
S0(isnan(S0))=0;
S(:,j)=S0;
end
S_mean0=mean(S,2);
S_mean(1:length(fx),1:M)=0;
peak_CS(1:M)=0;
for j=1:M
S_mean(:,j)=S_mean0;
end

RMS_S=sqrt(sum((S-S_mean).^2,2)/M);

RMS_total=sqrt(sum(RMS_S.^2));

plot(RMS_S);
%value_DC_FD=CS_envelope(1,j);
%FWHM_DC_FD=2*(space_FD(find(CS_envelope(:,j)<0.5*value_DC_FD, 1, 'first'))-space_FD(1));

%%  Interfered PSF FWHM
%space_min_for_inter_peak_FD=space_FD(1)+FWHM_DC_FD;  %要求well resolved
%space_min_for_inter_peak_index_FD=find(space_FD>space_min_for_inter_peak_FD, 1, 'first');
%[inter_peakvalue_FD inter_peakindex_FD]=max(CS_envelope(space_min_for_inter_peak_index_FD:round(length(CS_envelope(:,j))/2),j)); %注意! 此處的inter_peakindex_FD會抓錯, 因為CS_envelope從space_min_for_inter_peak_index_FD開始而非0開始
%inter_position_FD(j)=space_FD(inter_peakindex_FD+space_min_for_inter_peak_index_FD-1);
%FWHM_FD_right=space_FD(find(CS_envelope(space_min_for_inter_peak_index_FD:round(length(CS_envelope(:,j))/2),j)>0.5*inter_peakvalue_FD, 1, 'last'));
%FWHM_FD_left=space_FD(find(CS_envelope(space_min_for_inter_peak_index_FD:round(length(CS_envelope(:,j))/2),j)>0.5*inter_peakvalue_FD, 1, 'first'));
%FWHM_inter_FD=FWHM_FD_right-FWHM_FD_left;
%dispersion_expansion_ratio(j)=FWHM_inter_FD/FWHM_DC_FD;
%interference_efficiency(j)=inter_peakvalue_FD*2;
%plot(space_FD,CS_normal,space_FD,CS_envelope);
%dlmwrite('Interferogram.txt',M,'delimiter','\t','newline','pc');
%BW=lambda(find(S0>0.5,1,'last'))-lambda(find(S0>0.5,1,'first'));
%x_Res=x(find(CS_envelope>0.5,1,'last'))-x(find(CS_envelope>0.5,1,'first'))
%;
%plot(space_TD,inter_position_FD,space_TD,interference_efficiency,space_TD,dispersion_expansion_ratio);

%dlmwrite('CS_envelope.txt',CS_envelope,'delimiter','\t','newline','pc');
%dlmwrite('interference_efficiency.txt',interference_efficiency','delimiter','\t','newline','pc');
%dlmwrite('dispersion_expansion_ratio.txt',dispersion_expansion_ratio','delimiter','\t','newline','pc');
%dlmwrite('inter_position_FD.txt',inter_position_FD','delimiter','\t','newline','pc');
%dlmwrite('RMS_CS.txt',RMS_CS','delimiter','\t','newline','pc');
%dlmwrite('SNR.txt',SNR','delimiter','\t','newline','pc');