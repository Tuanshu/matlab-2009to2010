clear all;
close all;
clc;

N_f=4096*8;
N_t=N_f*8;

%temp=importdata('J:\100721_SDOCT with OO\100721_set2(1ms)_R+S(interfered)-3 (5times averaged).txt');
inter=importdata('J:\100823\100823_inter9_600mA.txt');
%ref=0;
ref=importdata('J:\100823\100823_inter9(ref)_600mA.txt');
sample=importdata('J:\100823\100823_inter9(sam)_600mA.txt');
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
CS_normal=CS;
d_t=1/(d_f*N_t);
%d_t=1/(d_f*2*N_f);
time=[-0.5*(N_t-1)*d_t:d_t:0.5*N_t*d_t]'/2;%/2是因為一來一回
%time=[-(N_f-1)*d_t:d_t:N_f*d_t]/2;  %/2是因為一來一回
space=c*time;
CS_envelope=abs(hilbert(CS_normal));
plot(space,CS_normal,space,CS_envelope);


%%  DC PSF FWHM
value_DC_peak=CS_envelope(1);
FWHM_DC_peak=2*(space(find(CS_envelope<0.5*value_DC_peak, 1, 'first'))-space(1));
End_of_DC=3.5*find(CS_envelope<0.5*value_DC_peak, 1, 'first');

space_min_for_inter_peak=space(1)+FWHM_DC_peak;  %要求well resolved
space_min_for_inter_peak_index=find(space>space_min_for_inter_peak, 1, 'first');
[inter_peakvalue inter_peakindex]=max(CS_envelope(space_min_for_inter_peak_index:round(length(CS_envelope)/2)));
CS_normaln=CS_normal;
CS_normaln(End_of_DC:length(CS_normal)-End_of_DC)=0;
CS_normaln=circshift(CS_normaln,space_min_for_inter_peak_index+inter_peakindex);
CS_envelopen=abs(hilbert(CS_normaln));

%%  Interfered PSF FWHM
%[peakvalue peakindex]=max(CS_envelope);
%FWHM_right=space(find(CS_envelopen>0.5*peakvalue, 1, 'last'));
%FWHM_left=space(find(CS_envelopen>0.5*peakvalue, 1, 'first'));
%FWHM_inter=FWHM_right-FWHM_left;
Sn=ifft(CS_normaln);
Sn=Sn(1:round(length(Sn)/2));
%dispersion_expansion_ratio=FWHM_inter/FWHM_DC_peak;
%interference_efficiency=peakvalue*2;
%plot(space,CS_normal,space,CS_normaln);
fxn=[0:d_f:d_f*(length(Sn)-1)]';
%plot(fxn,real(Sn));



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
%% add artifical dispersion (注意要取real的位置!!! 這裡假設corss talk為incoherent summation

left_f=2.5E14;
right_f=5.5E14;

left_index=find(fxn>left_f,1,'first');

Sn=Sn(1:find(fxn<right_f,1,'last'));
fxn=fxn(1:find(fxn<right_f,1,'last'));


Sn_real=real(Sn);
Sn_real_evelope=abs(hilbert(Sn_real));
Sn_real_total=Sn_real_evelope+Sn_real*0.8;

n_pixel=4096;
pixel_size=(length(fxn)-left_index)/4096;
o_size=pixel_size;
d_size=5*pixel_size/(length(fxn)-left_index);
%o_size=1;
%d_size=0;

% 也就是物理意義上的spot size
window_size(1:length(fxn)-left_index)=0;
for j=1:(length(fxn)-left_index)
    window_size(j)=round(o_size+d_size*j);
end


Sn_real_total_conv(1:(length(fxn)))=0;
for j=1:(length(fxn)-left_index)
    if j+window_size(j)<=(length(fxn)-left_index)
    Sn_real_total_conv(left_index+j)=sum(Sn_real_total(left_index+j:left_index+j+window_size(j)))/window_size(j).^2;
    %原本只要除一個window, 嘗試考慮collection eff效應, 多除一個R^2 (假設uniform disk)
    end
end

plot(fxn,real(Sn_real_total_conv));



%% again fft

CSnn=real(fft(Sn_real_total_conv,N_t))';     %with minus time
%CS=real(fft(S_padded));     %with minus time
CS_normalnn=CSnn;
%d_t=1/(d_f*2*N_f);
time=[-0.5*(N_t-1)*d_t:d_t:0.5*N_t*d_t]'/2;%/2是因為一來一回
%time=[-(N_f-1)*d_t:d_t:N_f*d_t]/2;  %/2是因為一來一回
CS_envelopenn=abs(hilbert(CS_normalnn));
plot(space,CS_normalnn,space,CS_envelopenn,space,CS_normal,space,CS_envelope);


%%  DC PSF FWHM
value_DC_peaknn=CS_envelopenn(1);
FWHM_DC_peaknn=2*(space(find(CS_envelopenn<0.5*value_DC_peaknn, 1, 'first'))-space(1));

space_min_for_inter_peak=space(1)+FWHM_DC_peak;  %要求well resolved
space_min_for_inter_peak_index=find(space>space_min_for_inter_peak, 1, 'first');
[inter_peakvalue inter_peakindex]=max(CS_envelope(space_min_for_inter_peak_index:round(length(CS_envelope)/2)));
FWHM_right=space(find(CS_envelope(space_min_for_inter_peak_index:round(length(CS_envelope)/2))>0.5*inter_peakvalue, 1, 'last'));
FWHM_left=space(find(CS_envelope(space_min_for_inter_peak_index:round(length(CS_envelope)/2))>0.5*inter_peakvalue, 1, 'first'));
FWHM_inter=FWHM_right-FWHM_left;

space_min_for_inter_peaknn=space(1)+FWHM_DC_peaknn;  %要求well resolved
space_min_for_inter_peak_indexnn=find(space>space_min_for_inter_peaknn, 1, 'first');
[inter_peakvaluenn inter_peakindexnn]=max(CS_envelopenn(space_min_for_inter_peak_indexnn:round(length(CS_envelopenn)/2)));
FWHM_rightnn=space(find(CS_envelopenn(space_min_for_inter_peak_indexnn:round(length(CS_envelopenn)/2))>0.5*inter_peakvaluenn, 1, 'last'));
FWHM_leftnn=space(find(CS_envelopenn(space_min_for_inter_peak_indexnn:round(length(CS_envelopenn)/2))>0.5*inter_peakvaluenn, 1, 'first'));
FWHM_internn=FWHM_rightnn-FWHM_leftnn;

FWHM_DC_peaknn=FWHM_DC_peaknn*1E6;
FWHM_DC_peak=FWHM_DC_peak*1E6;

FWHM_inter=FWHM_inter*1E6;
FWHM_internn=FWHM_internn*1E6;