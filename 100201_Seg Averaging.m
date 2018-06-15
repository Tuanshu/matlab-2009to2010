clear all
cd('J:\');
filename='100210_2Pol test (20sec).txt';
A=dlmread(filename,'');
Doff=1520.292000;
D360=347;   %for degree calibration
time=A(:,1);
data=A(:,2);        % usually power (W)
X=A(:,4);       % 3:given voltage, 4:rotator degree (relative)
averaging_time=1;   %(sec), only the data BEFORE each "change of X", and the time difference to the event: "change of X" is less than this averaging_time; will be taken
p=1;
for j=1:(length(X)-2) %to find the events: "change of X"
if (X(j)==X(j+1))&&(X(j+1)~=X(j+2))                                 %注意這一行包含兩個條件
    end_index(p)=j+1;
    p=p+1;
end 
end
averaged_data(1:length(end_index))=0;
start_data(1:length(end_index))=0;
start_index(1:length(end_index))=0;
for p=1:length(end_index) 
    start_index(p)=find(time>(time(end_index(p))-averaging_time),1,'first');
    averaged_data(p)=mean(data(start_index(p):end_index(p)));
    start_data(p)=data(start_index(p));
end 
X=360/D360*(X-Doff);
plot(X(end_index),averaged_data/7.918939,X(end_index),cos(X(end_index)/180*pi).^4);
M=[X(end_index)'; averaged_data;averaged_data/7.918939;cos(X(end_index)/180*pi).^4']';