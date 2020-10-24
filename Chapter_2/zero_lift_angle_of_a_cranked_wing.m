clear all;close all; clc;

b = 10.6 * u.m;
b_1 = 2* 3.18 * u.m;
b_2 = 2* 2.12 *u.m;
c_r_1 = 1.44 * u.m;
c_r = c_r_1;
c_t_1 = 1.44 * u.m;
c_r_2 = c_t_1;
c_t_2 = 0.86 * u.m;
c_t=c_t_2;
alfa_0L_r_1 = -2.5*u.deg;
alfa_0L_r_1_deg = alfa_0L_r_1/u.deg;
alfa_0L_r_1_rad = alfa_0L_r_1/u.rad;
alfa_0L_r_2  = alfa_0L_r_1;
alfa_0L_r_2_deg = alfa_0L_r_2/u.deg;
alfa_0L_r_2_rad =alfa_0L_r_2/u.rad;
alfa_0L_t_1 = -2.5*u.deg;
alfa_0L_t_1_deg = alfa_0L_t_1/u.deg;
alfa_0L_t_1_rad = alfa_0L_t_1/u.rad;
alfa_0L_t_2 = -1*u.deg;
alfa_0L_t_2_deg = alfa_0L_t_2/u.deg;
alfa_0L_t_2_rad = alfa_0L_t_2/u.rad;
epsilon_G_0L_r_1 = 0*u.deg;
epsilon_G_0L_r_1_deg = epsilon_G_0L_r_1/u.deg;
epsilon_G_0L_r_1_rad = epsilon_G_0L_r_1/u.rad;
epsilon_G_0L_r_2  = epsilon_G_0L_r_1;
epsilon_G_0L_r_2_deg = epsilon_G_0L_r_2/u.deg;
epsilon_G_0L_r_2_rad = epsilon_G_0L_r_2/u.rad;
epsilon_G_0L_t_1 = 0*u.deg;
epsilon_G_0L_t_1_deg = epsilon_G_0L_t_1/u.deg;
epsilon_G_0L_t_1_rad = epsilon_G_0L_t_1/u.rad;
epsilon_G_0L_t_2 = -3*u.deg;
epsilon_G_0L_t_2_deg = epsilon_G_0L_t_2/u.deg;
epsilon_G_0L_t_2_rad = epsilon_G_0L_t_2/u.rad;
coeff_a_c_1 = (c_t_1-c_r_1)/(b_1/2);
coeff_a_c_2 = (c_t_2-c_r_2)/(b_2/2);
coeff_b_c_1 = c_r_1;
coeff_b_c_2 = c_r_2;
coeff_a_0L_1 = (alfa_0L_t_1_rad-alfa_0L_r_1_rad)/(b_1/2);
coeff_a_0L_2 = (alfa_0L_t_2_rad-alfa_0L_r_2_rad)/(b_2/2);
coeff_b_0L_1 = alfa_0L_r_1_rad;
coeff_b_0L_2 = alfa_0L_r_2_rad;
coeff_a_eps_1 = (epsilon_G_0L_t_1_rad-epsilon_G_0L_r_1_rad)/(b_1/2);
coeff_a_eps_2 = (epsilon_G_0L_t_2_rad-epsilon_G_0L_r_2_rad)/(b_2/2);
coeff_b_eps_1 = epsilon_G_0L_r_1_rad;
coeff_b_eps_2 = epsilon_G_0L_r_2_rad;
taper_ratio_1 = c_t_1/c_r_1;
taper_ratio_2 = c_t_2/c_r_2;
wing_surface_1 = (b_1/2)*c_r_1*(1+taper_ratio_1);
wing_surface_1_mt_squared = wing_surface_1/u.m2
wing_surface_2 = (b_2/2)*c_r_2*(1+taper_ratio_2);
wing_surface_2_mt_squared = wing_surface_2/u.m2
wing_surface = wing_surface_1+wing_surface_2
wing_surface_mt_squared = wing_surface*u.m2;
lambda_LE_1 = 0 *u.deg;
lambda_LE_1_rad = lambda_LE_1 /u.rad;
lambda_LE_1_deg = lambda_LE_1 /u.deg;
lambda_LE_2 = 0 *u.deg;
lambda_LE_2_rad = lambda_LE_2 /u.rad;
lambda_LE_2_deg = lambda_LE_2 /u.deg;
coeff_a_1 = (c_t_1-c_r_1)/(b_1/2);
coeff_b_1 = c_r_1;
coeff_a_2 = (c_t_2-c_r_2)/(b_2/2);
coeff_b_2 = c_r_2;
fun=@(y) ((((coeff_a_0L_1*u.m)*y+coeff_b_0L_1)-(coeff_a_eps_1*u.m)*y).*((coeff_a_c_1)*y+coeff_b_c_1/u.m));
alpha_0L_1 =( 2/(wing_surface/u.m^2))* integral(fun,0,3.2);
alpha_0L_1_rad = alpha_0L_1*u.rad;
alpha_0L_1_deg = alpha_0L_1_rad/u.deg;
fun=@(y) ((((coeff_a_0L_2*u.m)*y+coeff_b_0L_2)-(coeff_a_eps_2*u.m)*y).*((coeff_a_c_2)*y+coeff_b_c_2/u.m));
alpha_0L_2 =( 2/(wing_surface/u.m^2))* integral(fun,3.2,5.3);
alpha_0L_2_rad = alpha_0L_2*u.rad;
alpha_0L_2_deg = alpha_0L_2_rad/u.deg;
alpha_0L_rad = alpha_0L_1_rad+alpha_0L_2_rad;
alpha_0L_deg = alpha_0L_1_deg+alpha_0L_2_deg;



%% Write data file
[status, msg] = mkdir("./zero_lift_angle_of_a_cranked_wing"); % create folder first
fid = fopen('./zero_lift_angle_of_a_cranked_wing/data.tex', 'w');
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
    fprintf(fid, "\\def\\mySweepLEWingIIDEG{%f}\n", lambda_LE_2_deg);
    fprintf(fid, "\\def\\myCoeffAChordWingI{%f}\n", coeff_a_1);
    fprintf(fid, "\\def\\myCoeffBChordWingIMT{%f}\n", coeff_b_1);
    fprintf(fid, "\\def\\myCoeffAChordWingII{%f}\n", coeff_a_2);
    fprintf(fid, "\\def\\myCoeffBChordWingIIMT{%f}\n", coeff_b_2);
    fprintf(fid, "\\def\\myAlphaZeroLiftRootWingIDEG{%f}\n", alfa_0L_r_1_deg);
    fprintf(fid, "\\def\\myAlphaZeroLiftTipWingIDEG{%f}\n", alfa_0L_t_1_deg);
    fprintf(fid, "\\def\\myAlphaZeroLiftRootWingIRAD{%f}\n", alfa_0L_r_1_rad);
    fprintf(fid, "\\def\\myAlphaZeroLiftTipWingIRAD{%f}\n", alfa_0L_t_1_rad);
    fprintf(fid, "\\def\\myAlphaZeroLiftRootWingIIDEG{%f}\n", alfa_0L_r_2_deg);
    fprintf(fid, "\\def\\myAlphaZeroLiftTipWingIIDEG{%f}\n", alfa_0L_t_2_deg);
    fprintf(fid, "\\def\\myAlphaZeroLiftRootWingIIRAD{%f}\n", alfa_0L_r_2_rad);
    fprintf(fid, "\\def\\myAlphaZeroLiftTipWingIIRAD{%f}\n", alfa_0L_t_2_rad);
    
    fprintf(fid, "\\def\\myTaperRatioWingI{%f}\n", taper_ratio_1);
    fprintf(fid, "\\def\\myTaperRatioWingII{%f}\n", taper_ratio_2);
    fprintf(fid, "\\def\\myTwistWingIDEG{%f}\n", epsilon_G_0L_t_1_deg);
    fprintf(fid, "\\def\\myTwistWingIRAD{%f}\n", epsilon_G_0L_t_1_rad);
    fprintf(fid, "\\def\\myTwistWingIIDEG{%f}\n", epsilon_G_0L_t_2_deg);
    fprintf(fid, "\\def\\myTwistWingIIRAD{%f}\n", epsilon_G_0L_t_2_rad);
    fprintf(fid, "\\def\\myAreaWingIMTsquared{%f}\n", wing_surface_1_mt_squared);
    fprintf(fid, "\\def\\myAreaWingIIMTsquared{%f}\n", wing_surface_2_mt_squared);
    fprintf(fid, "\\def\\myAreaWingMTsquared{%f}\n", wing_surface_mt_squared);
    fprintf(fid, "\\def\\myCoeffAAeroTwistWingIRADMT{%f}\n", coeff_a_0L_1);
    fprintf(fid, "\\def\\myCoeffBAeroTwistWingIRAD{%f}\n", coeff_b_0L_1);
    fprintf(fid, "\\def\\myCoeffAAeroTwistWingIIRADMT{%f}\n", coeff_a_0L_2);
    fprintf(fid, "\\def\\myCoeffBAeroTwistWingIIRAD{%f}\n", coeff_b_0L_2);
    fprintf(fid, "\\def\\myCoeffATwistWingIRADMT{%f}\n", coeff_a_eps_1);
    fprintf(fid, "\\def\\myCoeffATwistWingIIRADMT{%f}\n", coeff_a_eps_2);
    fprintf(fid, "\\def\\myCoeffBTwistWingIRAD{%f}\n", coeff_b_eps_1);
    fprintf(fid, "\\def\\myCoeffBTwistWingIIRAD{%f}\n", coeff_b_eps_2);

    
    fprintf(fid, "\\def\\myAlphaZeroLiftWingIRAD{%f}\n", alpha_0L_1_rad );
    fprintf(fid, "\\def\\myAlphaZeroLiftWingIDEG{%f}\n", alpha_0L_1_deg);
    fprintf(fid, "\\def\\myAlphaZeroLiftWingIIRAD{%f}\n", alpha_0L_2_rad );
    fprintf(fid, "\\def\\myAlphaZeroLiftWingIIDEG{%f}\n", alpha_0L_2_deg); 
    fprintf(fid, "\\def\\myAlphaZeroLiftWingRAD{%f}\n", alpha_0L_rad );
    fprintf(fid, "\\def\\myAlphaZeroLiftWingDEG{%f}\n", alpha_0L_deg); 
    % ...
    fclose(fid);
end