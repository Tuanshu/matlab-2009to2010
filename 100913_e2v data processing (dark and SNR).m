clear all;
close all;
clc;

N_f=4096;
N_t=N_f;


c=3E8;              %m/sec

min_wavelength=550;     %nm
max_wavelength=1050;    %nm

%inter=4096*gaussmf(1:4096,[100 2458])'+wgn(4096,1,10);
%temp=importdata('J:\100721_SDOCT with OO\100721_set2(1ms)_R+S(interfered)-3 (5times averaged).txt');
%inter=importdata('I:\100906_about SNR\100906_inter for SNR_325mA_1kHz.txt');
%inter=importdata('I:\100823\100823_inter9_600mA.txt');
%ref=0;
%ref=importdata('J:\100906_about SNR\100906_ref for SNR_325mA_30(26)kHz.txt');
%sample=importdata('J:\100906_about SNR\100906_sam for SNR_325mA_30(26)kHz.txt');
%AA=wgn(4096,1,1,'linear');
AA=gaussmf(1:4096,[10 2458])';
RMS_AA=sqrt(sum((AA-mean(AA)).^2)/length(AA));
sum_AA=sum(AA);
BB=abs(fft(AA));
RMS_BB=sqrt(sum((BB-mean(BB)).^2)/length(BB));
sum_BB=sum(BB);
QQQ=RMS_BB/RMS_AA;
CC=ifft(fft(AA));
RMS_CC=sqrt(sum((CC-mean(CC)).^2)/length(CC));
sum_CC=sum(CC);
RATIO=sum_BB/sum_AA;
dark1=importdata('I:\100908\100908_dark spectrum_1kHz_1.txt');
dark2=importdata('I:\100908\100908_dark spectrum_1kHz_2.txt');
dark3=importdata('I:\100908\100908_dark spectrum_1kHz_3.txt');
dark4=importdata('I:\100908\100908_dark spectrum_1kHz_4.txt');
dark5=importdata('I:\100908\100908_dark spectrum_1kHz_5.txt');
dark6=importdata('I:\100908\100908_dark spectrum_1kHz_6.txt');
dark7=importdata('I:\100908\100908_dark spectrum_1kHz_7.txt');
dark8=importdata('I:\100908\100908_dark spectrum_1kHz_8.txt');
dark9=importdata('I:\100908\100908_dark spectrum_1kHz_9.txt');
dark10=importdata('I:\100908\100908_dark spectrum_1kHz_10.txt');
dark11=importdata('I:\100908\100908_dark spectrum_1kHz_11_(blocked).txt');
dark12=importdata('I:\100908\100908_dark spectrum_1kHz_12_(blocked).txt');
dark13=importdata('I:\100908\100908_dark spectrum_1kHz_13_(blocked).txt');
dark14=importdata('I:\100908\100908_dark spectrum_1kHz_14_(blocked).txt');
dark15=importdata('I:\100908\100908_dark spectrum_1kHz_15_(blocked).txt');
dark16=importdata('I:\100908\100908_dark spectrum_20kHz_1_(blocked).txt');
%inter=(dark1+dark2+dark3+dark4+dark5)/5;
%inter=dark1;
%ref=0;
ref=importdata('I:\100906_about SNR\100906_ref for SNR_325mA_1kHz.txt');
sample=importdata('I:\100906_about SNR\100906_sam for SNR_325mA_1kHz.txt');
%temp=temp.data;



%%grating related (the center wavelength is the NOT the "blazed" wavelength)
grating_pitch=1000000/600;      %nm
center_wavelength=762;       %nm
incidence_angle=asin(center_wavelength/2/grating_pitch);     %假設是照此角度入射, rad

assumed_BW=178;

short_coef=tan(asin((center_wavelength/2-assumed_BW/2)/grating_pitch)-incidence_angle);         %tan!
long_coef=tan(asin((center_wavelength/2+assumed_BW/2)/grating_pitch)-incidence_angle);          %tan!

pixel=[1:4096]';  %1~4096
[peak index_peak]=max(ref);
index_long=find(ref>0.5*max(ref),1,'last');
index_short=find(ref>0.5*max(ref),1,'first');
%Q=178/(index_long-index_short);
%lambda=((-(pixel-index_peak)*Q+center_wavelength));   %前面的負號很重要! 影響freq domain是否接近guassian (不過為什麼差一個負號dispersion的broaden效應好像也會變大? 因為carrier in lambda domain有chirp, 但在freq domain應該不太有)

%Q=(long_coef-short_coef)/(index_long-index_short);      %roughly flens/pixel size
%lambda=grating_pitch*sin((asin((index_peak-pixel)*Q)+incidence_angle))+center_wavelength/2;

Q=(long_coef-short_coef)/(index_long-index_short);      %roughly flens/pixel size
lambda=grating_pitch*sin((atan((index_peak-pixel)*Q)+incidence_angle))+center_wavelength/2;   %atan!

%使用新的演算法後, 最短wavelength從166micron -> 10micron, 可能造成點數不夠, 所以先去掉不要的data
max_wavelength_index=find(lambda<max_wavelength,1,'first');
min_wavelength_index=find(lambda>min_wavelength,1,'last');
inter=inter(max_wavelength_index:min_wavelength_index);
ref=ref(max_wavelength_index:min_wavelength_index)-ref(min_wavelength_index);
sample=sample(max_wavelength_index:min_wavelength_index)-sample(min_wavelength_index);
lambda=lambda(max_wavelength_index:min_wavelength_index);

total_delta_frequency=c/(min_wavelength*1E-9)-c/(max_wavelength*1E-9);
d_lambda=(max_wavelength-min_wavelength)/(min_wavelength_index-max_wavelength_index);
d_f_max=c/(min_wavelength*1E-9)-c/((min_wavelength+d_lambda)*1E-9);
d_f_min=c/((max_wavelength-d_lambda)*1E-9)-c/(max_wavelength*1E-9);     %take the for interpolation
%lambda2=lambda2(lambda>550);
%ref=0;
%sample=0;
%index_i=(1:length(inter)).^-1;



%% z space
index_i=(501:(length(inter)+500)).^-1;
index_i=fliplr(index_i);
index_in=0:max(index_i)/length(index_i):max(index_i);
interi=interp1(index_i,inter',index_in);
CS_zspace=real(fft(inter,N_t))';     %with minus time
interi(isnan(interi))=0;
CS_zspace_i=real(fft(interi,N_t))';
inter_part(1:length(inter))=0;
inter_part(round(length(inter_part)/2):length(inter))=inter(round(length(inter_part)/2):length(inter));
CS_part=real(fft(inter_part,N_t))';
RMS_z_domain1=sqrt(sum((CS_zspace(101:2048)-mean(CS_zspace(101:2048))).^2)/length(CS_zspace(101:2048)));
RMS_z_domaini1=sqrt(sum((CS_zspace_i(101:1100)-mean(CS_zspace_i(101:1100))).^2)/length(CS_zspace_i(101:1100)));
%RMS_z_domain2=sqrt(sum((CS_zspace(60001:70000)-mean(CS_zspace(60001:70000))).^2)/length(CS_zspace(60001:70000)));
%RMS_z_domaini2=sqrt(sum((CS_zspace_i(60001:70000)-mean(CS_zspace_i(60001:70000))).^2)/length(CS_zspace_i(60001:70000)));
%plot(lambda,inter,lambda2,inter);

%S0=inter-ref-sample;
S0=inter;
%S0=inter;
%plot(lambda,S0/max(S0),lambda,inter/max(inter),lambda,ref/max(ref));
%S0=S0/max(S0);
%% Background subtraction (lambda)

%S0=S0-mean(S0);


freq=c./(lambda*1E-9);     %Hz

N_f_justok=max(freq)/d_f_min;

d_f=max(freq)/(N_f-1);
fx=0:d_f:max(freq);
S=interp1(freq,S0,fx);
S(isnan(S))=0;
%% Background subtraction (freq)

%S=S-mean(S);


freq_index_long=find(S>0.5*max(S),1,'last');
freq_index_short=find(S>0.5*max(S),1,'first');

CS=real(fft(S,N_t))';     %with minus time

%% Background subtraction (time)
%CS=CS-mean(CS);
%CS=real(fft(S_padded));     %with minus time
CS_normal=CS;
%CS_normal=CS/max(abs(CS));
d_t=1/(d_f*N_t);
%d_t=1/(d_f*2*N_f);
time=[-0.5*(N_t-1)*d_t:d_t:0.5*N_t*d_t]'/2;%/2是因為一來一回
%time=[-(N_f-1)*d_t:d_t:N_f*d_t]/2;  %/2是因為一來一回
space=c*time;
CS_envelope=abs(hilbert(CS_normal));



%%  DC PSF FWHM
value_DC_peak=CS_envelope(1);
FWHM_DC_peak=2*(space(find(CS_envelope<0.5*value_DC_peak, 1, 'first'))-space(1));

%%  Interfered PSF FWHM
space_min_for_inter_peak=space(1)+FWHM_DC_peak;  %要求well resolved
space_min_for_inter_peak_index=find(space>space_min_for_inter_peak, 1, 'first');
[inter_peakvalue inter_peakindex]=max(CS_envelope(space_min_for_inter_peak_index:round(length(CS_envelope)/2)));
FWHM_right=space(find(CS_envelope(space_min_for_inter_peak_index:round(length(CS_envelope)/2))>0.5*inter_peakvalue, 1, 'last'));
FWHM_left=space(find(CS_envelope(space_min_for_inter_peak_index:round(length(CS_envelope)/2))>0.5*inter_peakvalue, 1, 'first'));
FWHM_inter=FWHM_right-FWHM_left;
dispersion_expansion_ratio=FWHM_inter/FWHM_DC_peak;
interference_efficiency=inter_peakvalue*2;
plot(space,20*log10(CS_normal),space,20*log10(CS_envelope));
%dlmwrite('Interferogram.txt',M,'delimiter','\t','newline','pc');
%BW=lambda(find(S0>0.5,1,'last'))-lambda(find(S0>0.5,1,'first'));
%x_Res=x(find(CS_envelope>0.5,1,'last'))-x(find(CS_envelope>0.5,1,'first'));
%dlmwrite('CS.txt',CS,'delimiter','\t','newline','pc');
%dlmwrite('space.txt',space,'delimiter','\t','newline','pc','precision',16);
RMS_lambda_domain=sqrt(sum((S0-mean(S0)).^2)/length(S0));
RMS_frequency_domain=sqrt(sum((S-mean(S)).^2)/length(S));
%RMS_frequency_domain=sqrt(sum((S(9001:9100)-mean(S(9001:9100))).^2)/length(S(9001:9100)));
%RMS_frequency_domain=sqrt(sum((S(17501:17600)-mean(S(9001:9100))).^2)/length(S(9001:9100)));
%RMS_time_domain1=sqrt(sum((CS(round(1/100*length(CS)):round(11/100*length(CS)))-mean(CS(round(1/100*length(CS)):round(11/100*length(CS))))).^2)/length(CS(round(1/100*length(CS)):round(11/100*length(CS)))));
%RMS_time_domain1=sqrt(sum((CS(round(46/100*length(CS)):round(54/100*length(CS)))-mean(CS(round(46/100*length(CS)):round(54/100*length(CS))))).^2)/length(CS(round(46/100*length(CS)):round(54/100*length(CS)))));
%RMS_time_domain2=sqrt(sum((CS(round(30/100*length(CS)):round(38/100*length(CS)))-mean(CS(round(30/100*length(CS)):round(38/100*length(CS))))).^2)/length(CS(round(30/100*length(CS)):round(38/100*length(CS)))));
%RMS_time_domain1=sqrt(sum((CS(inter_peakindex*10:N_f_justok/2)-mean(CS(inter_peakindex*10:N_f_justok/2))).^2)/length(CS(inter_peakindex*10:N_f_justok/2)));
RMS_time_domain1=sqrt(sum((CS(201:N_f_justok/2)-mean(CS(201:N_f_justok/2))).^2)/length(CS(201:N_f_justok/2)));

%RMS_time_domain2=sqrt(sum((CS(50001:60000)-mean(CS(50001:60000))).^2)/length(CS(50001:60000)));
%RMS_time_domain3=sqrt(sum((CS(101:65536)-mean(CS(101:65536))).^2)/length(CS(101:65536)));
SNR_max1=20*log10(max(CS)/RMS_time_domain1);
SNR=20*log10(inter_peakvalue/RMS_time_domain1);
%RMS_time_domain_series(1:length(CS_normal/2)-(inter_peakindex*10))=0;
%for p=(inter_peakindex*10)+1:length(CS_normal/2)
%RMS_time_domain_series(p)=sqrt(sum((CS(inter_peakindex*10:p)-mean(CS(inter_peakindex*10:p))).^2)/length(CS(inter_peakindex*10:p)));
%end

RMS_time_domain_series(1:round(length(CS_normal)/2-(200)))=0;
for p=(200)+1:round(length(CS_normal)/2)
RMS_time_domain_series(p)=sqrt(sum((CS(inter_peakindex*10:p)-mean(CS(inter_peakindex*10:p))).^2)/length(CS(inter_peakindex*10:p)));
end

%SNR_max2=20*log10(max(CS)/RMS_time_domain2);
%SNR_max3=20*log10(max(CS)/RMS_time_domain3);
%plot(1:length(CS_normal),20*log10(CS_normal),1:length(CS_normal),20*log10(CS_zspace),1:length(CS_normal),20*log10(CS_zspace_i));
plot(1:length(CS_normal),20*log10(CS_normal));

%WIDTH_lambda_domain=index_long-index_short;
%WIDTH_inter_frequency_domain=freq_index_long-freq_index_short;
%WIDTH_DC_time_domain=2*(find(CS_envelope<0.5*value_DC_peak, 1, 'first')-1);
%WIDTH_inter_time_domain=find(CS_envelope(space_min_for_inter_peak_index:round(length(CS_envelope)/2))>0.5*inter_peakvalue, 1, 'last')-find(CS_envelope(space_min_for_inter_peak_index:round(length(CS_envelope)/2))>0.5*inter_peakvalue, 1, 'first');

