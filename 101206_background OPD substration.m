clear all;
close all;
clc;


backgroundopd=importdata('D:\inter_index_FD.txt');
data=importdata('D:\CS_envelope.txt');
clear data_red data_new
data_red=data(2001:end,:);
data_new(1:size(data,1)-2000,1:size(data,2))=0;
backgroundopd=backgroundopd-min(backgroundopd)+2000;

for j=1:size(data,2)
    data_new(1:size(data_new,1)-backgroundopd(j),j)=data(((backgroundopd(j)+1)):end-2000,j)/max(data(((backgroundopd(j)+1)):end-2000,j));
end



clear maxs_2d;
maxs=max(data_red(1:8000,1:180),[],1);

clear data_norm
for j=1:180
data_red(1:8000,j)=data_red(1:8000,j)/maxs(j);
end

imagesc(10*log(data_red(3000:8000,1:180)),'xdata',[10:10:10*180],'ydata',[1:0.009375850226843:0.009375850226843*5000]); figure(gcf);
imagesc(10*log(data_new(2000:7000,1:180)),'xdata',[10:10:10*180],'ydata',[1:0.009375850226843:0.009375850226843*5000]); figure(gcf);