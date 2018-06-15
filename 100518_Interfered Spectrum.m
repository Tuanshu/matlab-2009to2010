clear all;
close all;
clc;

power=5*10^-6;                                                                              %W, known power
N=4096;
cared_bandwidth=0.35;                                                                       %micron
delta_lambda=cared_bandwidth/N;
R=2;

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

%% Interfering (其實是假的 因為未轉成frequency base)

real_spectral_power((length(real_spectral_power)+1):(2*length(real_spectral_power)))=0;     %minus frequency
real_temporal_power=fft(real_spectral_power);                                               %complex
real_A=real_temporal_power(1:0.5*length(real_temporal_power));
                                                                %peak to center
real_B=real_temporal_power((0.5*length(real_temporal_power)+1):length(real_temporal_power));
real_temporal_power=[real_B;real_A];  
r(1:length(real_temporal_power))=0;                                                         %depth information
r(100)=1;                                                                                  %a peak somwhere
r(length(r)-100+1)=1;
real_temporal_power=convn(r,real_temporal_power,'same');

real_spectral_power=ifft(real_temporal_power)+real_spectral_power;
plot(real(real_spectral_power));



pixel_power(1:N)=0;
pixel_lambda(1:N)=0;
initial_lambda=peak_lambda-0.5*cared_bandwidth;
j_min_index(1:N)=0;
j_max_index(1:N)=0;
for j=1:N
    if j==1
        j_min_index(j)=min_index;
    else
        j_min_index(j)=j_max_index(j-1)+1;
    end
    j_max_index(j)=find(lambda<(initial_lambda+j*delta_lambda),1,'last');
    pixel_power(j)=sum(real_spectral_power(j_min_index(j):j_max_index(j)).*diff_lambda(j_min_index(j):j_max_index(j)));
    %pixel_lambda(j)=lambda(round((j_min_index+j_max_index)/2));
    pixel_lambda(j)=lambda(j_min_index(j));
end
total_power=sum(pixel_power);

%plot(pixel_lambda,pixel_power);                           

