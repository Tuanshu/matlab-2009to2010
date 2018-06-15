clear all;
close all;
clc;

power=5*10^-6;                                                                              %W, known power
N=4096;
cared_bandwidth=0.35;                                                                       %micron
delta_lambda=cared_bandwidth/N;
R=50;

temp=importdata('J:\100408\test-1W 12milisecond Nd filter 2+0.5 after.txt');
lambda=temp(:,1)*10^-3;                                                                     %(micron)
spectral_power_dB=temp(:,2);                                                                %dbm, a.u.
spectral_power=10.^(spectral_power_dB/10)-1;                                                %a.u., -1是把0dB減掉
spectral_power(length(spectral_power)+1:2*length(spectral_power))=0;                        %2N0
temporal_power=fft(spectral_power);                                                         %2N0
A=temporal_power(1:0.5*length(temporal_power));
B=temporal_power((0.5*length(temporal_power)+1):length(temporal_power));
Zeros(1:(2*R-2)*length(A))=0;
temporal_power=[A;Zeros';B];                                                                %zero padding, total 2RN0 -> R times
spectral_power=abs(ifft(temporal_power));
spectral_power=spectral_power(1:0.5*length(spectral_power));                                %RN0
lambda=interp(lambda,R,8,0.5);                                                              %interpolation參數我是try & error的
lambda=lambda(~isnan(lambda));                            
spectral_power=spectral_power(~isnan(lambda));
plot(lambda,spectral_power);
diff_lambda=diff(lambda);
spectral_power=spectral_power(1:(length(spectral_power)-1));
lambda=lambda(1:(length(lambda)-1));
integrated_spectral_power=sum(spectral_power.*diff_lambda);
normalized_spectral_power=spectral_power/integrated_spectral_power;                         %norm/micron
[peak_normalized_spectral_power,index_peak_lambda]=max(normalized_spectral_power);
peak_lambda=lambda(index_peak_lambda);
min_index=find(lambda>peak_lambda-0.5*cared_bandwidth,1,'first');
max_index=find(lambda<peak_lambda+0.5*cared_bandwidth,1,'last');
partial_integrated_spectral_power=sum(spectral_power(min_index:max_index).*diff_lambda(min_index:max_index))/integrated_spectral_power;
real_spectral_power=normalized_spectral_power*power;                                        %W/micron
pixel_power(1:N)=0;                                                                         %W!
pixel_lambda(1:N)=0;
initial_lambda=peak_lambda-0.5*cared_bandwidth;
j_max_index=-1;
for j=1:N
    j_min_index=max(min_index,j_max_index+1);
    j_max_index=find(lambda<(initial_lambda+j*delta_lambda),1,'last');
    pixel_power(j)=sum(real_spectral_power(j_min_index:j_max_index).*diff_lambda(j_min_index:j_max_index));
    %pixel_lambda(j)=lambda(round((j_min_index+j_max_index)/2));
    pixel_lambda(j)=lambda(j_min_index);
end
total_power=sum(pixel_power);

plot(pixel_lambda,pixel_power);                           

