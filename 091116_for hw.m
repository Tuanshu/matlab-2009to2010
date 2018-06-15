clear all;
M1=[2 0 0; 0 2 0; 1 1 1];
M=[-1 1 1;1 -1 1; 1 1 -1];
M3=M1/M;
M3i=M3^-1;
DET=det(M3);
DETi=det(M3i);
a=