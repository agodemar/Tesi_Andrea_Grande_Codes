clear all;close all; clc;

b = 26.8 * u.m;
c_r = 5.20 * u.m;
c_t = 2.34 * u.m;
lambda_C4 = 0*u.deg;
lambda_C4_rad = lambda_C4/u.rad;
lambda_C4_deg = lambda_C4/u.deg;
taper_ratio = c_t/c_r;
Mach = 0.7;
wing_surface = (b/2)*c_r*(1+taper_ratio);
wing_surface_mt_squared = wing_surface/u.m2;
aspect_ratio = (b^2)/wing_surface;
coeff_a_c = (c_t-c_r)/(b/2);
coeff_b_c = c_r;
cl_alpha_r = 6.15*(1/u.rad);
cl_alpha_r_deg = cl_alpha_r*u.deg;
cl_alpha_r_rad = cl_alpha_r*u.rad;
cl_alpha_t = 6.05*(1*u.rad);
cl_alpha_t_deg = cl_alpha_t*u.deg;
cl_alpha_t_rad = cl_alpha_t*u.rad;
coeff_a_cl = (cl_alpha_t_rad-cl_alpha_r_rad)/(b/2);
coeff_a_cl_radm = coeff_a_cl*u.rad;
coeff_a_cl_degm = coeff_a_cl*u.deg;
coeff_b_cl = cl_alpha_r_rad;
fun=@(y) (((coeff_a_cl_radm*u.m)*y+coeff_b_cl).*((coeff_a_c)*y+coeff_b_c/u.m));
cl_alpha_med =( 2/(wing_surface/u.m^2))* integral(fun,0,(b/2)/u.m);
cl_alpha_med_rad = cl_alpha_med*u.rad;
cl_alpha_med_deg = cl_alpha_med*u.deg;
oswald_factor = 0.90;
cl_alpha = cl_alpha_med_rad/(1+(cl_alpha_med_rad/(pi*oswald_factor*aspect_ratio)));
cl_alpha_rad = cl_alpha*u.rad;
cl_alpha_deg = cl_alpha*u.deg;
 

%% Write data file
[status, msg] = mkdir("./lift_gradient_of_a_finite_wing"); % create folder first
fid = fopen('./lift_gradient_of_a_finite_wing/data.tex', 'w');
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

    fprintf(fid, "\\def\\myCoeffAChordWing{%f}\n", coeff_a_c);
    fprintf(fid, "\\def\\myCoeffBChordWingMT{%f}\n", coeff_b_c);
    fprintf(fid, "\\def\\myCLAlphaRootWingRAD{%f}\n", cl_alpha_r_rad);
    fprintf(fid, "\\def\\myCLAlphaTipWingRAD{%f}\n", cl_alpha_t_rad);
    
    fprintf(fid, "\\def\\myCoeffAClalphaWingRADMT{%f}\n", coeff_a_cl_radm);
    fprintf(fid, "\\def\\myCoeffAClalphaWingDEGMT{%f}\n", coeff_a_cl_degm);
     fprintf(fid, "\\def\\myCoeffBClalphaWingRAD{%f}\n", coeff_b_cl);
      fprintf(fid, "\\def\\mySweepQuarterChordWingDEG{%f}\n", lambda_C4_deg);
    fprintf(fid, "\\def\\mySweepQuarterChordWingRAD{%f}\n", lambda_C4_rad);
    fprintf(fid, "\\def\\myCLAlphaMeanWingRAD{%f}\n", cl_alpha_med_rad);
    fprintf(fid, "\\def\\myCLAlphaMeanWingDEG{%f}\n", cl_alpha_med_deg);
    fprintf(fid, "\\def\\myInducedDragFactorWing{%f}\n", oswald_factor);
    fprintf(fid, "\\def\\myCLAlphaWingRAD{%f}\n", cl_alpha_rad);
    fprintf(fid, "\\def\\myCLAlphaWingDEG{%f}\n", cl_alpha_deg);
    % ...
    fclose(fid);
end
 
 