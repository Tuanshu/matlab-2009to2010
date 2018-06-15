clear all;
temp=importdata('J:\100419\ch1.txt');
sig_X=temp(60,:);
X=[1:length(sig_X)]';
windowsize=2000;
[p,index_p]=max(sig_X);
sig_X_n(1:length(sig_X))=0;
sig_X_n(length(sig_X_n)/2-windowsize:length(sig_X_n)/2+windowsize)=sig_X(index_p-windowsize:index_p+windowsize);
A=sig_X_n(1:length(sig_X_n)/2);                     %負空間
B=sig_X_n(length(sig_X_n)/2+1:length(sig_X_n));     %正空間
sig_X_n_n=[B A]';
sig_X_n_n=sig_X_n_n/max(sig_X_n_n);
sig_F=fft(sig_X_n_n);       %0~length(sig_X_n)是正頻
zeros(1:length(sig_X_n)/2)=0;
P=sig_F(1:length(sig_F)/2);                     %正頻
Q=sig_F(length(sig_F)/2+1:length(sig_F));       %負頻
sig_F=[P;Q];
sig_X_r=ifft(sig_F,40960);
sig_X_r_r=real(sig_X_r)/max(real(sig_X_r));
X_2=[1:length(sig_X_r_r)]';
plot(X_2,sig_X_r_r,X,sig_X_n_n);
X_1=[1:40960/10000:length(sig_X_r_r)]';
%plot(X,angle(sig_F),X,sig_F/max(sig_F(1000:3000)));
%sig_X=abs(hilbert(sig_X));
%X=1:1000000;            %Dawson function
%sig_X=gaussmf(X,[1000 500000]);
%sig_X_h=hilbert(sig_X);
%sig_F=fft(hilbert(sig_X));
%sig_XX=(ifft(fft(sig_X)));
%plot(X,abs(sig_XX),X,phase(sig_XX));
%plot(X,abs(sig_X));
%sig_X=temp(60,:);
%resolution_X=0.04;      %micron
%position_X=0:resolution_X:(length(sig_X)-1)*resolution_X;
%sig_X_enve=abs(hilbert(sig_X));
%sig_F=fft(sig_X_enve);
%plot(angle(sig_F));
%X=1:1000000;            %Dawson function
%sig_X=gaussmf(X,[5 500000]);
%sig_F=fft(sig_X);
%plot(abs(sig_F));

%plot(X,ifft(hilbert(real(fft(sig_X)))));


