clear all
cd('D:\100224 3D\3D\ppln_ch1-z(10um)\');
material='ppln';
descript='ch1(10um)';
index=1:200;
sum(1:1600,1:200)=0;
sumpeakindex(1:200)=0;
for j=1:length(index)
filename=sprintf('%s_%s_%g',material,descript,index(j));
filename_txt=sprintf('%s.txt',filename);
%filename_png=sprintf('%s.png',filename);
data=dlmread(filename_txt,'');
%img=mat2gray(data,[min(min(data)) max(max(data))]);
%imwrite(img,filename_png,'Bitdepth',16);
[peak,peakindex]=max(data);
sumpeakindex=sumpeakindex+peakindex;
sum=sum+data;
end

ave=sum/200;
avepeakindex=sumpeakindex/200;
filename_png='ave.png';
dlmwrite('ave.txt',ave,'delimiter','\t','newline','pc');
dlmwrite('avepeakindex.txt',avepeakindex,'delimiter','\t','newline','pc');
img=mat2gray(ave,[min(min(ave)) max(max(ave))]);
imwrite(img,'ave.png','Bitdepth',16);
plot(avepeakindex);