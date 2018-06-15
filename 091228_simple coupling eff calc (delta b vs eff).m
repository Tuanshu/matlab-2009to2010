clear all;
format long;
k=1;                    %coupling constant (1, as unit)
delta_b=-120:0.001:120;
kL=[pi/100:pi/100:pi/100];
L=kL/k;                 %(micron)
eff=zeros(length(kL),length(delta_b));
FWHM=zeros(length(kL),1);
peak=zeros(length(kL),1);
for j=1:length(kL)
L=kL(j)/k; 
eff(j,:)=(abs(k)^2./(abs(k)^2+(delta_b/2).^2)).*sin(abs(k)*L.*(1+(delta_b/2/abs(k)).^2).^0.5).^2;
FWHM(j)=delta_b(find(eff(j,:)>0.5*max(eff(j,:)), 1, 'last'))-delta_b(find(eff(j,:)>0.5*max(eff(j,:)), 1, 'first'));
peak(j)=max(eff(j,:));
end
eff=eff';
T=10*log10(1-eff);
delta_b=delta_b';
dlmwrite('eff.txt',eff,'delimiter','\t','newline','pc','precision','%.7f');
dlmwrite('delta_b.txt',delta_b,'delimiter','\t','newline','pc','precision','%.7f');
dlmwrite('T.txt',T,'delimiter','\t','newline','pc','precision','%.7f');
dlmwrite('FWHM.txt',FWHM,'delimiter','\t','newline','pc','precision','%.7f');
dlmwrite('peak.txt',peak,'delimiter','\t','newline','pc','precision','%.7f');
dlmwrite('kL.txt',kL','delimiter','\t','newline','pc','precision','%.7f');