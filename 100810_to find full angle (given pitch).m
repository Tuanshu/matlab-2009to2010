clear all;

pitch=1000/300;        %micron, grating pitch
blazed_wavelength=0.76;     %micron, center wavelength
s_wavelength=0.6;           %shortest wavelength cared
l_wavelength=1.1;           %longest wavelength cared

lens_size=5;             %mm


s_angle=asin(s_wavelength./pitch-blazed_wavelength/2./pitch)*180/pi;   %degree
l_angle=asin(l_wavelength./pitch-blazed_wavelength/2./pitch)*180/pi;   %degree
full_angle=l_angle-s_angle;
blazed_angle=asin(blazed_wavelength/2/pitch)*180/pi;   %degree

required_focus=lens_size/2/tan(full_angle*pi/180/2);