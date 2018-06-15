clear all;
close all;
clc;

N_f=4096*16;
N_t=N_f*16;

%% 假設速度uniform

data=importdata('D:\123.txt');
%data=data(750:1750,:);
space_TD=data(:,1)*1000;      %micron 
v=0.04; %mm/sec
%sr=4000;
%ns=2000;
%dspace=v/sr*2000*1000;
dspace=abs(space_TD(round(length(space_TD)/2),1)-space_TD(round(length(space_TD)/2)+1,1));
data_array=data(:,2:end);

data1D(1:size(data_array,1)*size(data_array,2))=0;
position1D(1:size(data_array,1)*size(data_array,2))=0;
for j=1:size(data_array,1)
    data1D(1+(j-1)*size(data_array,2):j*size(data_array,2))=data_array(j,:);
end
for j=1:size(data_array,1)*size(data_array,2)
    position1D(j)=dspace*j/size(data_array,2);
end


%% TD data calculation
power=data1D(1:end)';
position1D=position1D(1:end)';
%power_DC=mean(power(1:500));
power_DC=mean(power(1:10));
power_AC=power-power_DC;
power_envelope=abs(hilbert(power_AC));              %做hilbert前一定要減掉DC
[max_envelope max_index]=max(power_envelope);
[real_max real_max_index]=max(power_AC);            %for reference spectrum
FWHM_right_index=find(power_envelope>0.5*max_envelope, 1, 'last');
FWHM_right=position1D(FWHM_right_index);
FWHM_left_index=find(power_envelope>0.5*max_envelope, 1, 'first');
FWHM_left=position1D(FWHM_left_index);
FWHM_TD=(FWHM_right-FWHM_left);       %micron          
position1D=position1D-position1D(max_index);

plot(position1D,power_AC/max(power_AC),position1D,power_envelope/max(power_envelope));

