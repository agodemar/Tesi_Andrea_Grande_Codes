clear all;close all; clc;

b = 26.8 * u.m;
b_1 = 2*7.37 * u.m;
b_2 = 2* 6.03 *u.m;
c_r_1 = 5.2 * u.m;
c_r = c_r_1;
c_t_1 = 3 * u.m;
c_r_2 = c_t_1;
c_t_2 = 2.2 * u.m;
c_t=c_t_2;
lambda_LE_1 = 32 *u.deg;
lambda_LE_1_rad = lambda_LE_1 /u.rad;
lambda_LE_1_deg = lambda_LE_1 /u.deg;
lambda_LE_2 = 12 *u.deg;
lambda_LE_2_rad = lambda_LE_2 /u.rad;
lambda_LE_2_deg = lambda_LE_2 /u.deg;
coeff_a_1 = (c_t_1-c_r_1)/(b_1/2);
coeff_b_1 = c_r_1;
coeff_a_2 = (c_t_2-c_r_2)/(b_2/2);
coeff_b_2 = c_r_2;
taper_ratio_1 = c_t_1/c_r_1;
taper_ratio_2 = c_t_2/c_r_2;
wing_surface_1 = (b_1/2)*c_r_1*(1+taper_ratio_1);
wing_surface_1_mt_squared = wing_surface_1/u.m2;
wing_surface_2 = (b_2/2)*c_r_2*(1+taper_ratio_2);
wing_surface_2_mt_squared = wing_surface_2/u.m2;
aspect_ratio_1= b_1^2 / wing_surface_1;
aspect_ratio_2= b_2^2 / wing_surface_2;
mean_chord_1 = (2/3)*c_r_1*(1+taper_ratio_1+taper_ratio_1^2)/(1+taper_ratio_1); 
x_le_mean_chord_1 = (b_1/6)*((1+2*taper_ratio_1)/(1+taper_ratio_1))*tan(lambda_LE_1_rad);
y_mean_chord_1 =(b_1/6)*((1+2*taper_ratio_1)/(1+taper_ratio_1));
mean_chord_2 = (2/3)*c_r_2*(1+taper_ratio_2+taper_ratio_2^2)/(1+taper_ratio_2); 
x_le_mean_chord_2 = (b_2/6)*((1+2*taper_ratio_2)/(1+taper_ratio_2))*tan(lambda_LE_2_rad);
y_mean_chord_2 = (b_2/6)*((1+2*taper_ratio_2)/(1+taper_ratio_2));
wing_surface = wing_surface_1+wing_surface_2;
wing_surface_mt_squared = wing_surface*u.m2;
aspect_ratio = b^2 / wing_surface ;
mean_chord = (wing_surface_1*mean_chord_1+wing_surface_2*mean_chord_2)/(wing_surface);
y_mean_chord = (mean_chord - coeff_b_1)/coeff_a_1;
x_le_mean_chord = y_mean_chord*tan(lambda_LE_1_rad);

%% Write data file
[status, msg] = mkdir("./geometric_characteristics_of_a_cranked_wing"); % create folder first
fid = fopen('./geometric_characteristics_of_a_cranked_wing/data.tex', 'w');
if (fid == -1)
    fprintf("Cannot open file.")
else
    fprintf(fid, "\\def\\mySpanWingMT{%f}\n", b);
    fprintf(fid, "\\def\\mySpanWingIMT{%f}\n", b_1);
    fprintf(fid, "\\def\\mySpanWingIIMT{%f}\n",b_2);

    fprintf(fid, "\\def\\myChordRootWingMT{%f}\n", c_r);
    fprintf(fid, "\\def\\myChordRootWingIMT{%f}\n", c_r_1);
    fprintf(fid, "\\def\\myChordRootWingIIMT{%f}\n", c_r_2);

    fprintf(fid, "\\def\\myChordTipWingMT{%f}\n", c_t);
    fprintf(fid, "\\def\\myChordTipWingIMT{%f}\n", c_t_1);
    fprintf(fid, "\\def\\myChordTipWingIIMT{%f}\n", c_t_2);
    
    fprintf(fid, "\\def\\mySweepLEWingIDEG{%f}\n", lambda_LE_1_deg);
    fprintf(fid, "\\def\\mySweepLEWingIIDEG{%f}\n", lambda_LE_2_deg);
    fprintf(fid, "\\def\\mySweepLEWingIRAD{%f}\n", lambda_LE_1_rad);
    fprintf(fid, "\\def\\mySweepLEWingIIRAD{%f}\n", lambda_LE_2_rad);
    
    fprintf(fid, "\\def\\myCoeffAChordWingI{%f}\n", coeff_a_1);
    fprintf(fid, "\\def\\myCoeffBChordWingIMT{%f}\n", coeff_b_1);
    fprintf(fid, "\\def\\myCoeffAChordWingII{%f}\n", coeff_a_2);
    fprintf(fid, "\\def\\myCoeffBChordWingIIMT{%f}\n", coeff_b_2);
    
    fprintf(fid, "\\def\\myTaperRatioWingI{%f}\n", taper_ratio_1);
    fprintf(fid, "\\def\\myTaperRatioWingII{%f}\n", taper_ratio_2);

    fprintf(fid, "\\def\\myAreaWingIMTsquared{%f}\n", wing_surface_1_mt_squared);
    fprintf(fid, "\\def\\myAreaWingIIMTsquared{%f}\n", wing_surface_2_mt_squared);
    fprintf(fid, "\\def\\myAreaWingCrankedMTsquared{%f}\n", wing_surface_mt_squared);
    fprintf(fid, "\\def\\myAspectRatioWingI{%f}\n", aspect_ratio_1);
    fprintf(fid, "\\def\\myAspectRatioWingII{%f}\n", aspect_ratio_2);
    fprintf(fid, "\\def\\myAspectRatioWingCranked{%f}\n", aspect_ratio);
    fprintf(fid, "\\def\\myMACWingIMT{%f}\n", mean_chord_1);
    fprintf(fid, "\\def\\myYMACWingIMT{%f}\n", y_mean_chord_1);
    fprintf(fid, "\\def\\myXMACLEToApexWingIMT{%f}\n", x_le_mean_chord_1);    
    fprintf(fid, "\\def\\myMACWingIIMT{%f}\n", mean_chord_2);
    fprintf(fid, "\\def\\myYMACWingIIMT{%f}\n", y_mean_chord_2);
    fprintf(fid, "\\def\\myXMACLEToApexWingIIMT{%f}\n", x_le_mean_chord_2);
    fprintf(fid, "\\def\\myMACWingCrankedMT{%f}\n", mean_chord);
    fprintf(fid, "\\def\\myYYMACWingCrankedMT{%f}\n", y_mean_chord);
    fprintf(fid, "\\def\\myXXMACLEToApexWingCrankedMT{%f}\n", x_le_mean_chord);

    % ...
    fclose(fid);
end


  