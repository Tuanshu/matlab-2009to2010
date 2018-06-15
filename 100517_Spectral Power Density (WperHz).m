clear all;
close all;
clc;

power=5*10^-6;                                                                              %W, known power
temp=importdata('J:\100408\test-1W 12milisecond Nd filter 2+0.5 after.txt');
index_lambda=temp(:,1)*10^-3;                                                               %(micron)
spectral_power_dB=temp(:,2);                                                                %dbm, a.u.
spectral_power=10.^(spectral_power_dB/10);                                                  %a.u.
index_lambda=index_lambda(~isnan(index_lambda));                            
spectral_power=spectral_power(~isnan(index_lambda));
c=3*10^8;                                                                                   %m/sec
index_frequency=c./(index_lambda*10^-6);                                                    %Hz
diff_index_frequency=-diff(index_frequency);
spectral_power=spectral_power(1:(length(spectral_power)-1));
index_frequency=index_frequency(1:(length(index_frequency)-1));
integrated_spectral_power=sum(spectral_power.*diff_index_frequency);
normalized_spectral_power=spectral_power/integrated_spectral_power;
real_spectral_power=normalized_spectral_power*power;                                        %W/Hz
plot(index_frequency,real_spectral_power);                           

