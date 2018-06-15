clear all;
theta=2/3*pi;
C=cos(theta);
S=sin(theta);
E=exp(i*theta);
Ca=[-C S 0; S C 0; 0 0 -1];
Cb=[-C -S 0; -S C 0; 0 0 -1];
Cc=[-1 0 0; 0 1 0; 0 0 -1];
C3=[C -S 0; S C 0; 0 0 1];
Ca_2=[0 1; 1 0];
Cb_2=[0 E^2; E 0];
Cc_2=[0 E; E^2 0];
C3_2=[E 0; 0 E^2];

Ca_O=[C -S 0; -S C 0; 0 0 1];
Cb_O=[C S 0; S -C 0; 0 0 1];
Cc_O=[1 0 0; 0 -1 0; 0 0 1];
C3_O=[C S 0; -S C 0; 0 0 1];

CaC3=Ca*C3;
CaC3_2=Ca_2*C3_2;
CaC3_O=Ca_O*C3_O;
C3Ca=C3*Ca;
C3Ca_2=C3_2*Ca_2;
C3Ca_O=C3_O*Ca_O;

Cb;CaC3; Cc; C3Ca;
