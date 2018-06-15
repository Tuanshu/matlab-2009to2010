clear all;

given_full_angle=30;         %degree
pitch=0.5:0.00001:1.5;        %micron, grating pitch
blazed_wavelength=0.76;     %micron, center wavelength
s_wavelength=0.6;           %shortest wavelength cared
l_wavelength=1.1;           %longest wavelength cared

s_angle=asin(s_wavelength./pitch-blazed_wavelength/2./pitch)*180/pi;   %degree
l_angle=asin(l_wavelength./pitch-blazed_wavelength/2./pitch)*180/pi;   %degree
full_angle=l_angle-s_angle;
goal_pitch_index=find(full_angle<given_full_angle,1,'first');
goal_pitch=pitch(goal_pitch_index);
goal_blazed_angle=asin(blazed_wavelength/2/goal_pitch)*180/pi;   %degree
goal_s_angle=asin(s_wavelength/goal_pitch-blazed_wavelength/2/goal_pitch)*180/pi;   %degree
goal_l_angle=asin(l_wavelength/goal_pitch-blazed_wavelength/2/goal_pitch)*180/pi;   %degree

goal_full_angle=goal_l_angle-goal_s_angle;         %degree
goal_linedensity=1/goal_pitch*10^3;                %(line/mm)
goal_difference_of_blazed_angle_s_angle=goal_blazed_angle-goal_s_angle;
