clear all;
close all;
clc;

data=importdata('D:\100830\100830_set2_1kHz_0.01mms_340mA(1.55to1.58).txt');
space=data(:,1);
data_array=data(:,2:end);
power=sum(data_array,2);
power_DC=mean(power(1:500));
power_AC=power-power_DC;
[max_envelope max_index]=max(power_envelope);
power_envelope=abs(hilbert(power_AC));              %做hilbert前一定要減掉DC
FWHM_right=space(find(power_envelope>0.5*max_envelope, 1, 'last'));
FWHM_left=space(find(power_envelope>0.5*max_envelope, 1, 'first'));
FWHM=(FWHM_right-FWHM_left)*1000;       %micron          


plot(space,power_envelope,space,power_AC);