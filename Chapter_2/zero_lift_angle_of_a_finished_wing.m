clear all;close all; clc;

b = 26.8 * u.m;
c_r = 5.20 * u.m;
c_t = 2.34 * u.m;
lambda_C4 = 0*u.deg;
lambda_C4_deg = lambda_C4/u.deg;
lambda_C4_rad = lambda_C4/u.rad;

alfa_0L_r = -3*u.deg;
alfa_0L_r_deg = alfa_0L_r/u.deg;
alfa_0L_r_rad = alfa_0L_r/u.rad;
alfa_0L_t = -2*u.deg;
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
fun=@(y) ((((coeff_a_0L*u.m)*y+coeff_b_0L)-(coeff_a_eps*u.m)*y).*((coeff_a_c)*y+coeff_b_c/u.m));
alpha_0L =( 2/(wing_surface/u.m^2))* integral(fun,0,(b/2)/u.m);
alpha_0L_rad = alpha_0L*u.rad;
alpha_0L_deg = alpha_0L_rad/u.deg;


%% Write data file
[status, msg] = mkdir("./zero_lift_angle_of_a_finished_wing"); % create folder first
fid = fopen('./zero_lift_angle_of_a_finished_wing/data.tex', 'w');
if (fid == -1)
    fprintf("Cannot open file.")
else
    fprintf(fid, "\\def\\mySpanWingMT{%f}\n", b);
    fprintf(fid, "\\def\\myChordRootWingMT{%f}\n", c_r);
    fprintf(fid, "\\def\\myChordTipWingMT{%f}\n", c_t);
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
    fprintf(fid, "\\def\\myAlphaZeroLiftWingRAD{%f}\n", alpha_0L_rad );
    fprintf(fid, "\\def\\myAlphaZeroLiftWingDEG{%f}\n", alpha_0L_deg);    
    
    % ...
    fclose(fid);
end

