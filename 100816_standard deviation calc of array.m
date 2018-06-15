clear all;
close all;
clc;

%temp=importdata('J:\100721_SDOCT with OO\100721_set2(1ms)_R+S(interfered)-3 (5times averaged).txt');
data=importdata('J:\100824_dark noise_69deg_1line_1 of 50microS_1.txt');
dark=data(:,2)/16;
RMS=sqrt(sum((dark-mean(dark)).^2)/length(dark));
%ref=0;