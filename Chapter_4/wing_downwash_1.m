clear all;close all; clc;

b = 16 * u.m;
c_r = 2.5 * u.m;
c_t = 1* u.m;
taper_ratio = c_t/c_r;
wing_surface = (b/2)*c_r*(1+taper_ratio);
wing_surface_mt_squared = wing_surface/u.m2;
aspect_ratio = b^2 / wing_surface;
lambda_LE = 0 *u.deg;
lambda_LE_rad = lambda_LE /u.rad;
lambda_LE_deg = lambda_LE /u.deg;
cl_alpha_r = 6.15*(1/u.rad);
cl_alpha_r_deg = cl_alpha_r*u.deg;
cl_alpha_r_rad = cl_alpha_r*u.rad;
cl_alpha_t = 6.05*(1*u.rad);
cl_alpha_t_deg = cl_alpha_t*u.deg;
cl_alpha_t_rad = cl_alpha_t*u.rad;
x_ac_2D_r = 0.25;
x_ac_2D_t = 0.25;
cm_ac_r = -0.080;
cm_ac_t = -0.100;
alpha_0L_r = -3*u.deg;
alpha_0L_r_deg = alpha_0L_r/u.deg;
alpha_0L_r_rad = alpha_0L_r/u.rad;
alpha_0L_t_ = -1.5*u.deg;
alpha_0L_t_deg = alpha_0L_t_/u.deg;
alpha_0L_t_rad = alpha_0L_t_/u.rad;
mach = 0.40;
coeff_a_c = (c_t-c_r)/(b/2);
coeff_b_c = c_r;
coeff_a_cl = (cl_alpha_t_rad-cl_alpha_r_rad)/(b/2);
coeff_a_cl_radm = coeff_a_cl*u.rad;
coeff_a_cl_degm = coeff_a_cl*u.deg;
coeff_b_cl = cl_alpha_r_rad;
fun=@(y) (((coeff_a_cl_radm*u.m)*y+coeff_b_cl).*((coeff_a_c)*y+coeff_b_c/u.m));
cl_alpha_med =( 2/(wing_surface/u.m^2))* integral(fun,0,(b/2)/u.m);
cl_alpha_med_rad = cl_alpha_med*u.rad;
cl_alpha_med_deg = cl_alpha_med*u.deg;
oswald_factor = 2/(2-aspect_ratio+ (4+aspect_ratio*aspect_ratio*(1+0.009*0.009))^(1/2));
cl_alpha = cl_alpha_med_rad/(1+(cl_alpha_med_rad/(pi*oswald_factor*aspect_ratio)));
cl_alpha_rad = cl_alpha*u.rad;
cl_alpha_deg = cl_alpha*u.deg;
lambda_t_max = -0.54*u.deg;
lambda_t_max_deg = lambda_t_max/u.deg;
lambda_t_max_rad = lambda_t_max/u.rad;
d_espilon_zero = 2*cl_alpha_rad/(pi*aspect_ratio*oswald_factor);
d_espilon = d_espilon_zero*(1-mach*mach)^(1/2);

%% Write data file
[status, msg] = mkdir("./wing_downwash_1"); % create folder first
fid = fopen('./wing_downwash_1/data.tex', 'w');
if (fid == -1)
    fprintf("Cannot open file.")
else
    fprintf(fid, "\\def\\mySpanWingMT{%f}\n", b);    
    fprintf(fid, "\\def\\myMach{%f}\n", mach);   
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
    fprintf(fid, "\\def\\myCLAlphaMeanWingRAD{%f}\n", cl_alpha_med_rad);
    fprintf(fid, "\\def\\myCLAlphaMeanWingDEG{%f}\n", cl_alpha_med_deg);
    fprintf(fid, "\\def\\myInducedDragFactorWing{%f}\n", oswald_factor);
    fprintf(fid, "\\def\\myCLAlphaWingRAD{%f}\n", cl_alpha_rad);
    fprintf(fid, "\\def\\myCLAlphaWingDEG{%f}\n", cl_alpha_deg);
    fprintf(fid, "\\def\\mySweepTmaxWingDEG{%f}\n", lambda_t_max_deg);
    fprintf(fid, "\\def\\mySweepTmaxWingRAD{%f}\n", lambda_t_max_rad);
    fprintf(fid, "\\def\\myDownwashGradientLLTAtMachZeroWing{%f}\n", d_espilon_zero);
    fprintf(fid, "\\def\\myDownwashGradientLLTWing{%f}\n", d_espilon);
    % ...
    fclose(fid);
end
