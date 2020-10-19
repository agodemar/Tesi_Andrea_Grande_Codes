clear all;close all; clc;

b = 26.8 * u.m;
c_r = 5.2 * u.m;
c_t = 1.6 * u.m;
taper_ratio = c_t/c_r;
coeff_a = (c_t-c_r)/(b/2);
coeff_b = c_r;
wing_surface = (b/2)*c_r*(1+taper_ratio);
wing_surface_mt_squared = wing_surface/u.m2;
aspect_ratio = b^2 / wing_surface;
lambda_LE = 27.5 *u.deg;
lambda_LE_rad = lambda_LE /u.rad;
lambda_LE_deg = lambda_LE /u.deg;
lambda_C4 = atan( tan(lambda_LE_rad)-((1-0.31)/(aspect_ratio*1.31))) *u.rad;
lambda_C4_rad = lambda_C4/u.rad;
lambda_C4_deg = lambda_C4/u.deg;
lambda_C2 = atan( tan(lambda_LE_rad)-(2*(1-0.31)/(aspect_ratio*1.31))) *u.rad;
lambda_C2_rad = lambda_C2/u.rad;
lambda_C2_deg = lambda_C2/u.deg;
lambda_TE = atan( tan(lambda_LE_rad)-(4*(1-0.31)/(aspect_ratio*1.31))) *u.rad;
lambda_TE_rad = lambda_TE/u.rad;
lambda_TE_deg = lambda_TE/u.deg;
mean_chord = (2/3)*c_r*(1+taper_ratio+taper_ratio^2)/(1+taper_ratio);
x_le_mean_chord = (b/6)*((1+2*taper_ratio)/(1+taper_ratio))*tan(lambda_LE_rad);
y_mean_chord =(b/6)*((1+2*taper_ratio)/(1+taper_ratio));
%% Write data file
[status, msg] = mkdir("./geometric_characteristics_of_a_wing_with_arrow_not_null"); % create folder first
fid = fopen('./geometric_characteristics_of_a_wing_with_arrow_not_null/data.tex', 'w');
if (fid == -1)
    fprintf("Cannot open file.")
else
    fprintf(fid, "\\def\\mySpanWingMT{%f}\n", b);
    fprintf(fid, "\\def\\myChordRootWingMT{%f}\n", c_r);
    fprintf(fid, "\\def\\myChordTipWingMT{%f}\n", c_t);
    fprintf(fid, "\\def\\myTaperRatioWing{%f}\n", taper_ratio);
    fprintf(fid, "\\def\\myAreaWingMTsquared{%f}\n", wing_surface_mt_squared);
    fprintf(fid, "\\def\\myMACWingMT{%f}\n", mean_chord);
    fprintf(fid, "\\def\\myYMACWingMT{%f}\n", y_mean_chord);
    fprintf(fid, "\\def\\myXLEMACWingMT{%f}\n", x_le_mean_chord);
    fprintf(fid, "\\def\\myLambdaCQuarterDeg{%f}\n", lambda_C4_deg);
    fprintf(fid, "\\def\\myLambdaCQuarterRad{%f}\n", lambda_C4_rad);
    fprintf(fid, "\\def\\myLambdaCHalfDeg{%f}\n", lambda_C2_deg);
    fprintf(fid, "\\def\\myLambdaCHalfRad{%f}\n", lambda_C2_rad);
    fprintf(fid, "\\def\\myLambdaTEDeg{%f}\n", lambda_TE_deg);
    fprintf(fid, "\\def\\myLambdaTERad{%f}\n", lambda_TE_rad);
    fprintf(fid, "\\def\\myLambdaLEDeg{%f}\n", lambda_LE_deg);
    fprintf(fid, "\\def\\myLambdaLERad{%f}\n", lambda_LE_rad);
    fprintf(fid, "\\def\\myAspectRatioWing{%f}\n", aspect_ratio);
    fprintf(fid, "\\def\\myCoefficientA{%f}\n", coeff_a);
    fprintf(fid, "\\def\\myCoefficientB{%f}\n", coeff_b);

    % ...
    fclose(fid);
end
