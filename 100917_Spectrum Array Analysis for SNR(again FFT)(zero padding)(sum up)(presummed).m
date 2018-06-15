clear all;
close all;
clc;

sum_up_factor=1000;

data=importdata('J:\101011\101011_500mA_0.02to-0.2.txt');
data_summed(1:round(size(data,1)/sum_up_factor),1:size(data,2))=0;
for j=1:size(data,1)/sum_up_factor
    data_summed(j,:)=sum(data((1+(j-1)*sum_up_factor):j*sum_up_factor,:),1);
end
%data_array_o=data_summed(1:250,2:end);
%time_record=data_summed(1:250,1);

data_array_o=data_summed(1:5,2:end);
time_record=data_summed(1:5,1);

clear data data_summed;

N_lambda=4096;    %zero padding
N_f=4*4096;         %interpolation
N_t=N_f;            %fft

%% ref&grating related (the center wavelength is the NOT the "blazed" wavelength)
grating_pitch=1000000/1200;      %nm
center_wavelength=762;       %nm
incidence_angle=asin(center_wavelength/2/grating_pitch);     %假設是照此角度入射, rad

assumed_BW=178;

short_coef=tan(asin((center_wavelength/2-assumed_BW/2)/grating_pitch)-incidence_angle);         %tan!
long_coef=tan(asin((center_wavelength/2+assumed_BW/2)/grating_pitch)-incidence_angle);          %tan!

min_wavelength=550;     %nm
max_wavelength=1050;    %nm

pixel_o=[1:4096]';  %1~4096
[peak_o index_peak_o]=max(data_array_o(1,:));
index_long_o=find(data_array_o(1,:)>0.5*max(data_array_o(1,:)),1,'last');
index_short_o=find(data_array_o(1,:)>0.5*max(data_array_o(1,:)),1,'first');
%Q=178/(index_long-index_short);
%lambda=((-(pixel-index_peak)*Q+center_wavelength));   %前面的負號很重要! 影響freq domain是否接近guassian (不過為什麼差一個負號dispersion的broaden效應好像也會變大? 因為carrier in lambda domain有chirp, 但在freq domain應該不太有)

%Q=(long_coef-short_coef)/(index_long-index_short);      %roughly flens/pixel size
%lambda=grating_pitch*sin((asin((index_peak-pixel)*Q)+incidence_angle))+center_wavelength/2;

Q_o=(long_coef-short_coef)/(index_long_o-index_short_o);      %roughly flens/pixel size
lambda_o=grating_pitch*sin((atan((index_peak_o-pixel_o)*Q_o)+incidence_angle))+center_wavelength/2;   %atan!

%使用新的演算法後, 最短wavelength從166micron -> 10micron, 可能造成點數不夠, 所以先去掉不要的data

max_wavelength_index_o=find(lambda_o<max_wavelength,1,'first');
min_wavelength_index_o=find(lambda_o>min_wavelength,1,'last');
lambda_o=lambda_o(max_wavelength_index_o:min_wavelength_index_o);
data_array_o=data_array_o(:,max_wavelength_index_o:min_wavelength_index_o);
%lambda2=lambda2(lambda>300);

%% Zero padding
data_array=interpft(data_array_o,N_lambda,2); %不知道為什麼這裡用real而非abs SNR會比較高

pixel=[1:N_lambda]';
%pixel_n=[1:N_lambda]';
[peak index_peak]=max(data_array(1,:));
index_long=find(data_array(1,:)>0.5*max(data_array(1,:)),1,'last');
index_short=find(data_array(1,:)>0.5*max(data_array(1,:)),1,'first');
Q=(long_coef-short_coef)/(index_long-index_short);      %roughly flens/pixel size
lambda=grating_pitch*sin((atan((index_peak-pixel)*Q)+incidence_angle))+center_wavelength/2;   %atan!
plot(abs(data_array(1,:)));

%% x-axis
c=3E8;                     %m/sec
freq=c./(lambda*1E-9);     %Hz orignal freq array
d_f=max(freq)/(N_f-1);
fx=0:d_f:max(freq);        %freq after interpolation
d_t=1/(d_f*N_t);
%d_t=1/(d_f*2*N_f);
time=[-0.5*(N_t-1)*d_t:d_t:0.5*N_t*d_t]'/2;%/2是因為一來一回
%time=[-(N_f-1)*d_t:d_t:N_f*d_t]/2;  %/2是因為一來一回
length_space_FD=round(length(time));
space_FD=c*time(1:length_space_FD); % 暫定只用array前1/100 (大約也有100 micron左右吧)
space_FD=(space_FD-space_FD(1))*1E6;                     % shift to zero & m to micron
%plot(lambda,inter,lambda2,inter);

%% data related
[M N]=size(data_array);
CS(1:length_space_FD,1:M)=0;
S(1:length(fx),1:M)=0;
inter_position_FD(1:M)=0;
%% SD summed in x-domain
for j=1:M
S0=interp1(freq,data_array(j,:),fx,'linear');
S0(isnan(S0))=0;
S(:,j)=S0-min(S0);
CS0=fft(S0,N_t)';     %with minus time
CS(:,j)=CS0(1:length_space_FD);
end

clear data_array S0 CS0;


CSr=real(CS);
CS_envelope=abs(CS);
CS_mean0=mean(real(CS),2);
CS_mean(1:length_space_FD,1:M)=0;
S_mean0=mean(S,2);
S_mean(1:length(fx),1:M)=0;
peak_CS(1:M)=0;

clear CS fx;

for j=1:M
peak_CS(j)=max(CS_envelope(:,j));
CS_mean(:,j)=CS_mean0;
S_mean(:,j)=S_mean0;
end

clear CS_mean0 S_mean0 CS_envelope;
RMS_CS=sqrt(sum((CSr-CS_mean).^2,2)/M);

plot(space_FD,20*log10(mean(peak_CS)/2./RMS_CS));

RMS_S=sqrt(sum((S-S_mean).^2,2)/M);

RMS_total=sqrt(sum(RMS_S.^2));

SNR=20*log10(mean(peak_CS)/2/RMS_total);
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