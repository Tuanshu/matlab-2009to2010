clear all
cd('C:\');
A=dlmread('DCF_Confocal_ref.txt','');
A=A-0.192;
offset=min(min(A));
dlmwrite('DCF_Confocal_ref_new.txt',A,'delimiter','\t','newline','pc','precision','%.7f');           %.........為什麼是顛倒的