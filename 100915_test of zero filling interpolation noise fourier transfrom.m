clear all;

f=wgn(100,1,1,'linear');

N_z=10000:10000:100000;       %for zero filling
N=100000;      %for interpolation


temp=1000000;
%plot(1:N_z,20*log10(Fz));
%% different N
for j=1:length(N_z)
fz=interpft(f,N_z(j));
fi_N=interp1(0.5+0.5./(100:-1:1),f,(1:N)/N);
fi_N(isnan(fi_N))=0;
fzi_N=interp1(0.5+0.5./(N_z(j):-1:1),fz,(1:N)/N);
fzi_N(isnan(fzi_N))=0;
Fi_N=abs(fft(fi_N));
Fzi_N=abs(fft(fzi_N));
SD_Fi_N=std(Fi_N);
SD_Fzi_N=std(Fzi_N);
if SD_Fzi_N<temp
    temp=SD_Fzi_N;
    temp_index=j;
end
end

%SD_fi_N=std(fi_N);



%SD_fzi_N=std(fzi_N);

F=abs(fft(f));

%Fz=abs(fft(fz));
N_z_best=N_z(temp_index);

plot(1:N,20*log10(Fi_N),1:N,20*log10(Fzi_N));
%plot(1:N,Fi_N,1:N,Fzi_N);
