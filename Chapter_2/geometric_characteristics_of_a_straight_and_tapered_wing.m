clear all;close all; clc;

b = 26 * u.m;
lambda_C4 = 0*u.deg;
c_r = 2.5 * u.m;
c_t = 1.25 * u.m;
taper_ratio = c_t/c_r;
coeff_a = (c_t-c_r)/(b/2);
coeff_b = c_r;
wing_surface = (b/2)*c_r*(1+taper_ratio);
wing_surface_mt_squared = wing_surface/u.m2;
aspect_ratio = b^2 / wing_surface;
lambda_LE = atan((0.5)/(aspect_ratio*1.5)) *u.rad;
lambda_LE_rad = lambda_LE /u.rad;
lambda_LE_deg = lambda_LE /u.deg;
mean_chord = (2/3)*c_r*(1+taper_ratio+taper_ratio^2)/(1+taper_ratio);
x_le_mean_chord = (b/6)*((1+2*taper_ratio)/(1+taper_ratio))*tan(lambda_LE);
y_mean_chord =(b/6)*((1+2*taper_ratio)/(1+taper_ratio));
%% Write data file
[status, msg] = mkdir("./geometric_characteristics_of_a_straight_and_tapered_wing"); % create folder first
fid = fopen('./geometric_characteristics_of_a_straight_and_tapered_wing/data.tex', 'w');
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
    fprintf(fid, "\\def\\myLambdaCQuarterDeg{%f}\n", lambda_C4);
    fprintf(fid, "\\def\\myLambdaLEDeg{%f}\n", lambda_LE_deg);
    fprintf(fid, "\\def\\myLambdaLERad{%f}\n", lambda_LE_rad);
    fprintf(fid, "\\def\\myAspectRatioWing{%f}\n", aspect_ratio);
    fprintf(fid, "\\def\\myCoefficientA{%f}\n", coeff_a);
    fprintf(fid, "\\def\\myCoefficientB{%f}\n", coeff_b);

    % ...
    fclose(fid);
end
