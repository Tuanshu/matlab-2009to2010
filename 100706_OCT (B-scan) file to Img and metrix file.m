clear all
cd('D:\100224 3D\3D\ppln_ch1-z(10um)\');
material='ppln';
descript='ch1(10um)';
pixelsize=0.19/1600;    %(mm)
index=1:200;
peakindexarray(1:200,1:200)=0;
for j=1:length(index)
filename=sprintf('%s_%s_%g',material,descript,index(j));
filename_txt=sprintf('%s.txt',filename);
%filename_png=sprintf('%s.png',filename);
data=dlmread(filename_txt,'');
%img=mat2gray(data,[min(min(data)) max(max(data))]);
%imwrite(img,filename_png,'Bitdepth',16);
[peak,peakindex]=max(data);
peakindexarray(j,1:200)=peakindex;
end
peakarray=pixelsize*peakindexarray(1:j,1:200);
avepeakarray=mean(peakarray);
stdpeakarray=std(peakarray);
dlmwrite('avepeakarray.txt',avepeakarray,'delimiter','\t','newline','pc');
dlmwrite('stdpeakarray.txt',stdpeakarray,'delimiter','\t','newline','pc');
plot(1:200,avepeakarray);