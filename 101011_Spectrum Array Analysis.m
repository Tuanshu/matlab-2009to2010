clear all;
close all;
clc;

N_f=4096*16;
N_t=N_f*16;

data=importdata('D:\101101_test(-0.715 to -0.735)(0.01).txt');
%data=data(750:1750,:);
space_TD=data(:,1)*1000;      %micron 
data_array=data(:,2:end);

%% TD data calculation
power=sum(data_array,2);
%power_DC=mean(power(1:500));
power_DC=mean(power(1:10));
power_AC=power-power_DC;
power_envelope=abs(hilbert(power_AC));              %做hilbert前一定要減掉DC
[max_envelope max_index]=max(power_envelope);
[real_max real_max_index]=max(power_AC);            %for reference spectrum
FWHM_right_index=find(power_envelope>0.5*max_envelope, 1, 'last');
FWHM_right=space_TD(FWHM_right_index);
FWHM_left_index=find(power_envelope>0.5*max_envelope, 1, 'first');
FWHM_left=space_TD(FWHM_left_index);
FWHM_TD=(FWHM_right-FWHM_left);       %micron          
space_TD=space_TD-space_TD(max_index);


%% FD data calculation

ref=data_array(real_max_index,:);

%% ref&grating related (the center wavelength is the NOT the "blazed" wavelength)
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

ref=ref(lambda>300);
%lambda2=lambda2(lambda>300);
lambda=lambda(lambda>300);

%% x-axis
c=3E8;                     %m/sec
freq=c./(lambda*1E-9);     %Hz orignal freq array
d_f=max(freq)/(N_f-1);
fx=0:d_f:max(freq);        %freq after interpolation
d_t=1/(d_f*N_t);
%d_t=1/(d_f*2*N_f);
time=[-0.5*(N_t-1)*d_t:d_t:0.5*N_t*d_t]'/2;%/2是因為一來一回
%time=[-(N_f-1)*d_t:d_t:N_f*d_t]/2;  %/2是因為一來一回
length_space_FD=round(length(time)/100);
space_FD=c*time(1:length_space_FD); % 暫定只用array前1/100 (大約也有100 micron左右吧)
space_FD=(space_FD-space_FD(1))*1E6;                     % shift to zero & m to micron
%plot(lambda,inter,lambda2,inter);

%% data related
[M N]=size(data_array);
CS_envelope(1:length_space_FD,1:M)=0;
inter_position_FD(1:M)=0;
dispersion_expansion_ratio(1:M)=0;
interference_efficiency(1:M)=0;
for j=1:M
S0=data_array(j,:);
S0=S0-mean(S0(lambda>300)); %background light substration
S0=S0(lambda>300);
%S0=inter;
%plot(lambda,S0/max(S0),lambda,inter/max(inter),lambda,ref/max(ref));
%S0=S0/max(S0);        %不先normal應該也沒差?
S=interp1(freq,S0,fx);
S(isnan(S))=0;
CS=real(fft(S,N_t))';     %with minus time
%CS=real(fft(S_padded));     %with minus time
CS_normal=CS/max(abs(CS));
CS_envelope(:,j)=abs(hilbert(CS_normal(1:round(length(CS_normal)/100))));                                        %data 1 to record (need decimate? maybe cut is better)
%% 暫定只用array前1/100 (大約也有100 micron左右吧)

%plot(space,CS_normal,space,CS_envelope);

%% After the acquire of envelope


%%  DC PSF FWHM
value_DC_FD=CS_envelope(1,j);
FWHM_DC_FD=2*(space_FD(find(CS_envelope(:,j)<0.5*value_DC_FD, 1, 'first'))-space_FD(1));

%%  Interfered PSF FWHM
space_min_for_inter_peak_FD=space_FD(1)+FWHM_DC_FD;  %要求well resolved
space_min_for_inter_peak_index_FD=find(space_FD>space_min_for_inter_peak_FD, 1, 'first');
[inter_peakvalue_FD inter_peakindex_FD]=max(CS_envelope(space_min_for_inter_peak_index_FD:round(length(CS_envelope(:,j))/2),j)); %注意! 此處的inter_peakindex_FD會抓錯, 因為CS_envelope從space_min_for_inter_peak_index_FD開始而非0開始
inter_position_FD(j)=space_FD(inter_peakindex_FD+space_min_for_inter_peak_index_FD-1);
FWHM_FD_right=space_FD(find(CS_envelope(space_min_for_inter_peak_index_FD:round(length(CS_envelope(:,j))/2),j)>0.5*inter_peakvalue_FD, 1, 'last'));
FWHM_FD_left=space_FD(find(CS_envelope(space_min_for_inter_peak_index_FD:round(length(CS_envelope(:,j))/2),j)>0.5*inter_peakvalue_FD, 1, 'first'));
FWHM_inter_FD=FWHM_FD_right-FWHM_FD_left;
dispersion_expansion_ratio(j)=FWHM_inter_FD/FWHM_DC_FD;
interference_efficiency(j)=inter_peakvalue_FD*2;
%plot(space_FD,CS_normal,space_FD,CS_envelope);
%dlmwrite('Interferogram.txt',M,'delimiter','\t','newline','pc');
%BW=lambda(find(S0>0.5,1,'last'))-lambda(find(S0>0.5,1,'first'));
%x_Res=x(find(CS_envelope>0.5,1,'last'))-x(find(CS_envelope>0.5,1,'first'))
%;
end
plot(space_TD,inter_position_FD,space_TD,interference_efficiency,space_TD,dispersion_expansion_ratio);

%dlmwrite('power_envelope.txt',power_envelope,'delimiter','\t','newline','pc');
%dlmwrite('power.txt',power,'delimiter','\t','newline','pc');
%dlmwrite('space_TD.txt',space_TD,'delimiter','\t','newline','pc');
%dlmwrite('CS_envelope.txt',CS_envelope,'delimiter','\t','newline','pc');
%dlmwrite('interference_efficiency.txt',interference_efficiency','delimiter','\t','newline','pc');
%dlmwrite('dispersion_expansion_ratio.txt',dispersion_expansion_ratio','delimiter','\t','newline','pc');
%dlmwrite('inter_position_FD.txt',inter_position_FD','delimiter','\t','newline','pc');