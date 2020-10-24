clear all;close all; clc;

b = 26.8 * u.m;
c_r = 5.20 * u.m;
c_t = 2.34 * u.m;
lambda_LE = 27.5*u.deg;
lambda_LE_deg = lambda_LE/u.deg;
lambda_LE_rad = lambda_LE/u.rad;
Mach = 0.7;
alfa_0L_r = -3*u.deg;
alfa_0L_r_deg = alfa_0L_r/u.deg;
alfa_0L_r_rad = alfa_0L_r/u.rad;
alfa_0L_t = -1.5*u.deg;
alfa_0L_t_deg = alfa_0L_t/u.deg;
alfa_0L_t_rad = alfa_0L_t/u.rad;
epsilon_G_0L_r = 0*u.deg;
epsilon_G_0L_r_deg = epsilon_G_0L_r/u.deg;
epsilon_G_0L_r_rad = epsilon_G_0L_r/u.rad;
epsilon_G_0L_t = -1.5*u.deg;
epsilon_G_0L_t_deg = epsilon_G_0L_t/u.deg;
epsilon_G_0L_t_rad = epsilon_G_0L_t/u.rad;
coeff_a_c = (c_t-c_r)/(b/2);
coeff_b_c = c_r;
coeff_a_0L = (alfa_0L_t_rad-alfa_0L_r_rad)/(b/2);
coeff_b_0L = alfa_0L_r_rad;
coeff_a_eps = (epsilon_G_0L_t_rad-epsilon_G_0L_r_rad)/(b/2);
coeff_b_eps = epsilon_G_0L_r_rad;
taper_ratio = c_t/c_r;
wing_surface = (b/2)*c_r*(1+taper_ratio);
wing_surface_mt_squared = wing_surface/u.m2;
tc_tip = 0.09;
tc_root = 0.15;
coeff_a_tc = (tc_tip-tc_root)/(b/2);
coeff_b_tc = tc_root;
aspect_ratio = (b^2)/wing_surface;
lambda_C4 = atan( tan(lambda_LE_rad)-((1-0.31)/(aspect_ratio*1.31))) *u.rad;
lambda_C4_rad = lambda_C4/u.rad;
lambda_C4_deg = lambda_C4/u.deg;

fun=@(y) (((coeff_a_0L*u.m)*y+coeff_b_0L).*((coeff_a_c)*y+coeff_b_c/u.m));
alpha_0L_med =( 2/(wing_surface/u.m^2))* integral(fun,0,13.4);
alpha_0L_med_rad = alpha_0L_med*u.rad;
alpha_0L_med_deg = alpha_0L_med_rad/u.deg;

fun=@(y) (((coeff_a_tc*u.m)*y+coeff_b_tc).*((coeff_a_c)*y+coeff_b_c/u.m));
tc_med =( 2/(wing_surface/u.m^2))* integral(fun,0,13.4);

k_alpha_1 = -0.381;
k_alpha_2 = 0.85;
alpha_0L = (alpha_0L_med_rad+k_alpha_1*epsilon_G_0L_t_rad)*k_alpha_2;
alpha_0L_rad = alpha_0L/u.rad;
alpha_0L_deg = alpha_0L/u.deg;


  
 
   
 
%% Write data file
[status, msg] = mkdir("./zero_lift_angle_graphic_method"); % create folder first
fid = fopen('./zero_lift_angle_graphic_method/data.tex', 'w');
if (fid == -1)
    fprintf("Cannot open file.")
else
    fprintf(fid, "\\def\\mySpanWingMT{%f}\n", b);    
    fprintf(fid, "\\def\\myMach{%f}\n", Mach);   
    fprintf(fid, "\\def\\myAspectRatioWing{%f}\n", aspect_ratio);   
    fprintf(fid, "\\def\\myChordRootWingMT{%f}\n", c_r);
    fprintf(fid, "\\def\\myChordTipWingMT{%f}\n", c_t);
    fprintf(fid, "\\def\\mySweepLEWingDEG{%f}\n", lambda_LE_deg);
    fprintf(fid, "\\def\\mySweepLEWingRAD{%f}\n", lambda_LE_rad);
    fprintf(fid, "\\def\\mySweepQuarterChordWingDEG{%f}\n", lambda_C4_deg);
    fprintf(fid, "\\def\\mySweepQuarterChordWingRAD{%f}\n", lambda_C4_rad);
    fprintf(fid, "\\def\\myCoeffAChordWing{%f}\n", coeff_a_c);
    fprintf(fid, "\\def\\myCoeffBChordWingMT{%f}\n", coeff_b_c);
    fprintf(fid, "\\def\\myAlphaZeroLiftRootWingDEG{%f}\n", alfa_0L_r_deg);
    fprintf(fid, "\\def\\myAlphaZeroLiftTipWingDEG{%f}\n", alfa_0L_t_deg);
    fprintf(fid, "\\def\\myAlphaZeroLiftRootWingRAD{%f}\n", alfa_0L_r_rad);
    fprintf(fid, "\\def\\myAlphaZeroLiftTipWingRAD{%f}\n", alfa_0L_t_rad);
    fprintf(fid, "\\def\\myTaperRatioWing{%f}\n", taper_ratio);
    fprintf(fid, "\\def\\myTwistWingDEG{%f}\n", epsilon_G_0L_t_deg);
    fprintf(fid, "\\def\\myTwistWingRAD{%f}\n", epsilon_G_0L_t_rad);
    fprintf(fid, "\\def\\myTwistWingDEG{%f}\n", taper_ratio);
    fprintf(fid, "\\def\\myAreaWingMTsquared{%f}\n", wing_surface_mt_squared);
    fprintf(fid, "\\def\\myCoeffAAeroTwistWingRADMT{%f}\n", coeff_a_0L);
    fprintf(fid, "\\def\\myCoeffATwistWingRADMT{%f}\n", coeff_a_eps);
    fprintf(fid, "\\def\\myCoeffBAeroTwistWingRAD{%f}\n", coeff_b_0L);
    
    fprintf(fid, "\\def\\myThicknessOverChordTipWing{%f}\n", tc_tip);
    fprintf(fid, "\\def\\myThicknessOverChordRootWing{%f}\n", tc_root);
    fprintf(fid, "\\def\\myThicknessOverChordMeanWing{%f}\n", tc_med);

    fprintf(fid, "\\def\\myCoeffAPercThicknessWingMT{%f}\n", coeff_a_tc);
    fprintf(fid, "\\def\\myCoeffBPercThicknessWing{%f}\n", coeff_b_tc);
    fprintf(fid, "\\def\\myAlphaZeroLiftMeanWingRAD{%f}\n", alpha_0L_med_rad );
    fprintf(fid, "\\def\\myAlphaZeroLiftMeanWingDEG{%f}\n", alpha_0L_med_rad);    
    fprintf(fid, "\\def\\myAlphaZeroLiftDatcomWingRAD{%f}\n",alpha_0L_rad );
    fprintf(fid, "\\def\\myAlphaZeroLiftDatcomWingDEG{%f}\n", alpha_0L_deg);
    fprintf(fid, "\\def\\myDeltaAlphaZeroLiftOverTwistWing{%f}\n",k_alpha_1 );
    fprintf(fid, "\\def\\myPrandtlGlauertCorrectionAlphaZeroLiftWing{%f}\n", k_alpha_2);
    % ...
    fclose(fid);
end


 