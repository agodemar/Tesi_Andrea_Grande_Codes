clear all;close all; clc;
b = 16 * u.m;
c_r = 2.5 * u.m;
c_t = 1* u.m;
taper_ratio = c_t/c_r;
wing_surface = (b/2)*c_r*(1+taper_ratio);
wing_surface_mt_squared = wing_surface/u.m2;
aspect_ratio = b^2 / wing_surface;
mean_chord = (2/3)*c_r*(1+taper_ratio+taper_ratio^2)/(1+taper_ratio);
v_cl = [0.2148,0.1935,0.1569,0.1156,0.0732,0.0324,-0.0054,-0.0392,-0.0683,-0.0926,-0.1119,-0.1262,-0.1355,-0.1398,-0.1389,-0.1326,-0.1203,-0.1013,-0.0747,-0.0404,0.0000]*u.m;
v_y = [0.00,0.63,1.25,1.87,2.47,3.06,3.63,4.18,4.70,5.20,5.66,6.08,6.47,6.82,7.13,7.39,7.61,7.78,7.90,7.98,8.00]*u.m;
v_x_b = [-0.1864,-0.1570,-0.1278,-0.0989,-0.0706,-0.0429,-0.0162,0.0095,0.0340,0.0571,0.0787,0.0987,0.1169,0.1333,0.1477,0.1600,0.1702,0.1782,0.1839,0.1874,0.1886]*u.m;
nPoints = 100;
v_y_f = linspace(v_y(1),v_y(end),nPoints);
v_cl_f = interp1(v_y,v_cl,v_y_f);
%plot(v_y,v_cl,'ko',v_y_f,v_cl_f,'b');
v_x_b_f = interp1(v_y,v_x_b,v_y_f);
vec = (v_cl_f/u.m).*(v_x_b_f/u.m);
plot(v_y_f,vec,'b')
cm_ac_b = (2/((mean_chord/u.m)*(wing_surface/u.m^2)))* trapz(v_y_f,vec); 

%% Write data file
[status, msg] = mkdir("./pitching_moment_two"); % create folder first
fid = fopen('./pitching_moment_two/data.tex', 'w');
if (fid == -1)
    fprintf("Cannot open file.")
else
    fprintf(fid, "\\def\\myChordRootWingMT{%f}\n", c_r);
    fprintf(fid, "\\def\\myChordTipWingMT{%f}\n", c_t);
    fprintf(fid, "\\def\\myCmZeroBWing{%f}\n", cm_ac_b);
% ...
    fclose(fid);
end

 
