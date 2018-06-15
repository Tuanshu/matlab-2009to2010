clear all;

lambda=0.44;                %micron

angleX=11.7;                %half
angleY=37.8;                %half, degree

wX0=lambda/pi/(angleX*pi/180);   %micron
wY0=lambda/pi/(angleY*pi/180);   %micron

%wX0=1;       %micron
%wY0=16;      %micron

%angleX=(lambda/pi/wX0)*180/pi;
%angleY=(lambda/pi/wY0)*180/pi;




zX0=(pi*wX0^2./lambda)*10^-6;   %m
zY0=(pi*wY0^2./lambda)*10^-6;   %m
f1=4.5;                       %mm
f2=8;                         %mm
increase=0.2;               %mm
spotdiff=100000;            %mm
for p=1:10/increase
    for q=1:150/increase
        for k=1:150/increase
d1=4.5+increase*(p-1);                        %mm, distance between sample and lens1
Md1=[1 d1*10^-3;0 1];         %m
Mlens1=[1 0;-1/(f1*10^-3) 1]; %m
M1=Mlens1*Md1;
A1=M1(1,1);
B1=M1(1,2);
C1=M1(2,1);
D1=M1(2,2);
qX0=i*zX0;
qY0=i*zY0;
qX1=(A1*qX0+B1)/(C1*qX0+D1);
qY1=(A1*qY0+B1)/(C1*qY0+D1);
%LX1=-real(qX1)*10^3;           %mm, the waist is L1 after the lens1
%zX1=imag(qX1);
%wX1=((zX1*10^6)*lambda/pi)^0.5;   %micron, half spot size
%LY1=-real(qY1)*10^3;           %mm, the waist is L1 after the lens1
%zY1=imag(qY1);
%wY1=((zY1*10^6)*lambda/pi)^0.5;   %micron, half spot size

d2=20+increase*(q-1);                     %mm, should be roughly L1+f2 for L2->inf
Md2=[1 d2*10^-3;0 1];        %m
Mlens2=[1 0;-1/(f2*10^-3) 1]; %m
M2=Mlens2*Md2;               %
A2=M2(1,1);
B2=M2(1,2);
C2=M2(2,1);
D2=M2(2,2);
qX2=(A2*qX1+B2)/(C2*qX1+D2);
qY2=(A2*qY1+B2)/(C2*qY1+D2);

d3=increase*(k-1);
Md3=[1 d3*10^-3;0 1];        %m
M3=Md3;               %
A3=M3(1,1);
B3=M3(1,2);
C3=M3(2,1);
D3=M3(2,2);
qX3=(A3*qX2+B3)/(C3*qX2+D3);
qY3=(A3*qY2+B3)/(C3*qY2+D3);

LX3=-real(qX3);           %m, the waist is L1 after the lens2
zX3=imag(qX3);
wX3=((zX3*10^6)*lambda/pi)^0.5;   %micron, half spot size
LY3=-real(qY3);           %m, the waist is L1 after the lens2
zY3=imag(qY3);
wY3=((zY3*10^6)*lambda/pi)^0.5;   %micron, half spot size
spotX3=wX3*(1+(LY3/zX3)^2)^0.5; %micron
spotY3=wY3*(1+(LY3/zY3)^2)^0.5; %micron

temp=abs(spotX3-spotY3);               %distance between 2 focal

if temp<=spotdiff
    spotdiff=temp;
    d1goal=d1;
    d2goal=d2;
    d3goal=d3;
    LX3goal=LX3;
    LY3goal=LY3;
    spotX3goal=spotX3;
    spotY3goal=spotY3;
end

        end
    end
end
