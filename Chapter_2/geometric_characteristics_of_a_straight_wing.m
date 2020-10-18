clear all;close all; clc;

b = 26 * u.m;
lambda_C4 = 0*u.deg;
taper_ratio = 1;
c = 2.5 * u.m;
lambda_LE = 0*u.deg;
wing_surface = (b/2)*c*(1+taper_ratio);
wing_surface_mt_squared = wing_surface/u.m2;
aspect_ratio = b^2 / wing_surface;
mean_chord = (2/3)*c*(1+taper_ratio+taper_ratio^2)/(1+taper_ratio);
x_le_mean_chord = (b/6)*((1+2*taper_ratio)/(1+taper_ratio))*tan(lambda_LE);
y_mean_chord =(b/6)*((1+2*taper_ratio)/(1+taper_ratio));
%% Write data file
[status, msg] = mkdir("./geometric_characteristics_of_a_straight_wing"); % create folder first
fid = fopen('./geometric_characteristics_of_a_straight_wing/data.tex', 'w');
if (fid == -1)
    fprintf("Cannot open file.")
else
    fprintf(fid, "\\def\\mySpanWingMT{%f}\n", b);
    fprintf(fid, "\\def\\myChordRootWingMT{%f}\n", c);
    fprintf(fid, "\\def\\myChordTipWingMT{%f}\n", c);
    fprintf(fid, "\\def\\myTaperRatioWing{%f}\n", taper_ratio);
    fprintf(fid, "\\def\\myAreaWingMTsquared{%f}\n", wing_surface_mt_squared);
    fprintf(fid, "\\def\\myMACWingMT{%f}\n", mean_chord);
    fprintf(fid, "\\def\\myYMACWingMT{%f}\n", y_mean_chord);
    fprintf(fid, "\\def\\myXLEMACWingMT{%f}\n", x_le_mean_chord);
    fprintf(fid, "\\def\\myLambdaCQuarter{%f}\n", lambda_C4);
    fprintf(fid, "\\def\\myLambdaLEDeg{%f}\n", lambda_LE);
    fprintf(fid, "\\def\\myLambdaLERad{%f}\n", lambda_LE/u.rad);
    fprintf(fid, "\\def\\myAspectRatioWing{%f}\n", aspect_ratio);

    % ...
    fclose(fid);
end
