clear all;
close all;
clc;

N_f=4096*8;
N_t=N_f*8;

%temp=importdata('J:\100721_SDOCT with OO\100721_set2(1ms)_R+S(interfered)-3 (5times averaged).txt');
inter=importdata('I:\100823\100823_inter9_600mA.txt');
%ref=0;
ref=importdata('I:\100823\100823_inter9(ref)_600mA.txt');
sample=importdata('I:\100823\100823_inter9(sam)_600mA.txt');
%temp=temp.data;

%%grating related (the center wavelength is the NOT the "blazed" wavelength)
grating_pitch=1000000/600;      %nm
center_wavelength=762;       %nm
incidence_angle=asin(center_wavelength/2/grating_pitch);     %���]�O�Ӧ����פJ�g, rad

assumed_BW=178;

short_coef=tan(asin((center_wavelength/2-assumed_BW/2)/grating_pitch)-incidence_angle);         %tan!
long_coef=tan(asin((center_wavelength/2+assumed_BW/2)/grating_pitch)-incidence_angle);          %tan!

pixel=[1:4096]';  %1~4096
[peak index_peak]=max(ref);
index_long=find(ref>0.5*max(ref),1,'last');
index_short=find(ref>0.5*max(ref),1,'first');
%Q=178/(index_long-index_short);
%lambda=((-(pixel-index_peak)*Q+center_wavelength));   %�e�����t���ܭ��n! �v�Tfreq domain�O�_����guassian (���L������t�@�ӭt��dispersion��broaden�����n���]�|�ܤj? �]��carrier in lambda domain��chirp, ���bfreq domain���Ӥ��Ӧ�)

%Q=(long_coef-short_coef)/(index_long-index_short);      %roughly flens/pixel size
%lambda=grating_pitch*sin((asin((index_peak-pixel)*Q)+incidence_angle))+center_wavelength/2;

Q=(long_coef-short_coef)/(index_long-index_short);      %roughly flens/pixel size
lambda=grating_pitch*sin((atan((index_peak-pixel)*Q)+incidence_angle))+center_wavelength/2;   %atan!

%�ϥηs���t��k��, �̵uwavelength�q166micron -> 10micron, �i��y���I�Ƥ���, �ҥH���h�����n��data

inter=inter(lambda>300);
ref=ref(lambda>300);
sample=sample(lambda>300);
%lambda2=lambda2(lambda>300);
lambda=lambda(lambda>300);

ref=0;
sample=0;

%plot(lambda,inter,lambda2,inter);

S0=inter-ref-sample;
%S0=inter;
%plot(lambda,S0/max(S0),lambda,inter/max(inter),lambda,ref/max(ref));
S0=S0/max(S0);


c=3E8;              %m/sec
freq=c./(lambda*1E-9);     %Hz
d_f=max(freq)/(N_f-1);
fx=0:d_f:max(freq);
S=interp1(freq,S0,fx);
S(isnan(S))=0;

CS=real(fft(S,N_t))';     %with minus time
%CS=real(fft(S_padded));     %with minus time
CS_normal=CS/max(abs(CS));
d_t=1/(d_f*N_t);
%d_t=1/(d_f*2*N_f);
time=[-0.5*(N_t-1)*d_t:d_t:0.5*N_t*d_t]'/2;%/2�O�]���@�Ӥ@�^
%time=[-(N_f-1)*d_t:d_t:N_f*d_t]/2;  %/2�O�]���@�Ӥ@�^
space=c*time;
CS_envelope=abs(hilbert(CS_normal));
plot(space,CS_normal,space,CS_envelope);


%%  DC PSF FWHM
value_DC_peak=CS_envelope(1);
FWHM_DC_peak=2*(space(find(CS_envelope<0.5*value_DC_peak, 1, 'first'))-space(1));

%%  Interfered PSF FWHM
space_min_for_inter_peak=space(1)+FWHM_DC_peak;  %�n�Dwell resolved
space_min_for_inter_peak_index=find(space>space_min_for_inter_peak, 1, 'first');
[inter_peakvalue inter_peakindex]=max(CS_envelope(space_min_for_inter_peak_index:round(length(CS_envelope)/2)));
FWHM_right=space(find(CS_envelope(space_min_for_inter_peak_index:round(length(CS_envelope)/2))>0.5*inter_peakvalue, 1, 'last'));
FWHM_left=space(find(CS_envelope(space_min_for_inter_peak_index:round(length(CS_envelope)/2))>0.5*inter_peakvalue, 1, 'first'));
FWHM_inter=FWHM_right-FWHM_left;
dispersion_expansion_ratio=FWHM_inter/FWHM_DC_peak;
interference_efficiency=inter_peakvalue*2;
plot(space,CS_normal,space,CS_envelope);
%dlmwrite('Interferogram.txt',M,'delimiter','\t','newline','pc');
%BW=lambda(find(S0>0.5,1,'last'))-lambda(find(S0>0.5,1,'first'));
%x_Res=x(find(CS_envelope>0.5,1,'last'))-x(find(CS_envelope>0.5,1,'first'));
%dlmwrite('100823_ref.txt',ref,'delimiter','\t','newline','pc','precision',10);
%dlmwrite('100823_inter.txt',inter,'delimiter','\t','newline','pc','precision',10);
%dlmwrite('100823_inter.txt',lambda,'delimiter','\t','newline','pc','precision',10);
%dlmwrite('100823_sam.txt',sample,'delimiter','\t','newline','pc','precision',10);
%dlmwrite('100823_fx.txt',fx,'delimiter','\t','newline','pc','precision',10);
%dlmwrite('100823_S.txt',S,'delimiter','\t','newline','pc','precision',10);
%dlmwrite('100823_CS_envelope.txt',CS_envelope,'delimiter','\t','newline','pc','precision',10);
%dlmwrite('100823_CS_normal.txt',CS_normal,'delimiter','\t','newline','pc','precision',10);
%dlmwrite('100823_space.txt',space,'delimiter','\t','newline','pc','precision',10);