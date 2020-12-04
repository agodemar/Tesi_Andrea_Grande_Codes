clear all;close all; clc;

b = 26.8 * u.m;
c_r = 5.20 * u.m;
c_t = 2.18 * u.m;
lambda_C4 = 0*u.deg;
lambda_C4_rad = lambda_C4/u.rad;
lambda_C4_deg = lambda_C4/u.deg;
lambda_LE = 27.5*u.deg;
lambda_LE_rad = lambda_LE/u.rad;
lambda_LE_deg = lambda_LE/u.deg;
taper_ratio = c_t/c_r;
Mach = 0.7;
Mach_cr = 0.79;
wing_surface = (b/2)*c_r*(1+taper_ratio);
wing_surface_mt_squared = wing_surface/u.m2;
aspect_ratio = (b^2)/wing_surface; 
lambda_C2 = atan( tan(lambda_LE_rad)-(2*(1-taper_ratio)/(aspect_ratio*(1+taper_ratio)))) *u.rad;
lambda_C2_rad = lambda_C2/u.rad;
lambda_C2_deg = lambda_C2/u.deg;
kp = 1+((8.2-2.3*lambda_LE_rad)-aspect_ratio*(0.22-0.153*lambda_LE_rad))/100;
cl_alpha = 2*pi*aspect_ratio/(2+(4+(aspect_ratio*aspect_ratio*(1-Mach*Mach)/(kp*kp))*(1+(tan(lambda_C2_rad)*tan(lambda_C2_rad))/(1-Mach*Mach)))^(1/2)); 
cl_alpha_rad = cl_alpha*u.rad
cl_alpha_deg = cl_alpha*u.deg

   
  
  
 
%% Write data file
[status, msg] = mkdir("./lift_gradient_polhamus_formula"); % create folder first
fid = fopen('./lift_gradient_polhamus_formula/data.tex', 'w');
if (fid == -1)
    fprintf("Cannot open file.")
else
    fprintf(fid, "\\def\\mySpanWingMT{%f}\n", b);    
    fprintf(fid, "\\def\\myMach{%f}\n", Mach);   
    fprintf(fid, "\\def\\myAspectRatioWing{%f}\n", aspect_ratio);   
    fprintf(fid, "\\def\\myChordRootWingMT{%f}\n", c_r);
    fprintf(fid, "\\def\\myChordTipWingMT{%f}\n", c_t);
    fprintf(fid, "\\def\\myTaperRatioWing{%f}\n", taper_ratio);
    fprintf(fid, "\\def\\myAreaWingMTsquared{%f}\n", wing_surface_mt_squared);

    fprintf(fid, "\\def\\myCriticalMachNumberMACWing{%f}\n", Mach_cr);
    fprintf(fid, "\\def\\mySweepLEWingDEG{%f}\n", lambda_LE_deg);
    fprintf(fid, "\\def\\mySweepLEWingRAD{%f}\n", lambda_LE_rad);
    fprintf(fid, "\\def\\mySweepHalfChordWingRAD{%f}\n", lambda_C2_rad);
    fprintf(fid, "\\def\\mySweepHalfChordWingDEG{%f}\n", lambda_C2_deg);
    
    fprintf(fid, "\\def\\myKPolhamus{%f}\n", kp);
    fprintf(fid, "\\def\\myCLAlphaWingRAD{%f}\n", cl_alpha_rad);
    fprintf(fid, "\\def\\myCLAlphaWingDEG{%f}\n", cl_alpha_deg);
    % ...
    fclose(fid);
end
 
 