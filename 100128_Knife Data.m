clear all
cd('J:\100128\');
filename='100128_scan 15 (250mi).txt';
displacement=0.25;  %(mm)
A=dlmread(filename,'');
Power=A(:,1);
Div=-A(:,2);
Position(1:length(Div))=0;
for j=1:length(Div)
    Position(j)=displacement*(j-1);
end
Div_n=Div/max(Div);
FWHM=(find(Div_n>0.5, 1, 'last')-find(Div_n>0.5, 1, 'first'))*250;
FWHM_1e2=(find(Div_n>0.135335283, 1, 'last')-find(Div_n>0.135335283, 1, 'first'))*250;
CENTER=(find(Div_n>0.5, 1, 'last')+find(Div_n>0.5, 1, 'first'))/2;
Position(1:length(Div))=0;
for j=1:length(Div)
    Position(j)=displacement*(j-CENTER);
end
Position=Position';
plot(Position,Div_n);