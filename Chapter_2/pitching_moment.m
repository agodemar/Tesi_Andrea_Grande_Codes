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
coeff_a_cl = (cl_alpha_t_rad-cl_alpha_r_rad)/(b/2);
coeff_b_cl = cl_alpha_r_rad;
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

lambda_C4 = atan( tan(lambda_LE_rad)-((1-0.31)/(aspect_ratio*1.31))) *u.rad;
lambda_C4_rad = lambda_C4/u.rad;
lambda_C4_deg = lambda_C4/u.deg;
mean_chord = (2/3)*c_r*(1+taper_ratio+taper_ratio^2)/(1+taper_ratio);
x_le_mean_chord = (b/6)*((1+2*taper_ratio)/(1+taper_ratio))*tan(lambda_LE_rad);
y_mean_chord =(b/6)*((1+2*taper_ratio)/(1+taper_ratio));
k1 = 1.342;
k2 = 0;
x_ac_cr = 0.176;
x_ac_mean = k1*(x_ac_cr-k2);
X_ac = x_le_mean_chord + x_ac_mean*mean_chord;
coeff_a_c = (c_t-c_r)/(b/2);
coeff_b_c = c_r;
coeff_a_cm = (cm_ac_t-cm_ac_r)/(b/2);
coeff_b_cm = cm_ac_r;


coeff_d =  coeff_a_cm*coeff_a_c*coeff_a_c;
coeff_b =  coeff_b_c*coeff_b_c*coeff_a_cm+2*coeff_a_c*coeff_b_c*coeff_b_cm;
coeff_c =  2*coeff_a_c*coeff_b_c*coeff_a_cm+coeff_b_cm*coeff_a_c*coeff_a_c;
coeff_a =  coeff_b_cm*coeff_b_c*coeff_b_c;
fun=@(y) (coeff_a/u.m^2+coeff_b/u.m*y+coeff_c*y.^2 + coeff_d*u.m* y.^3);
cm_ac_a = (2/((mean_chord/u.m)*(wing_surface/u.m^2)))* integral(fun,0,(b/2)/u.m);

coeff_a_0L = (alpha_0L_t_rad-alpha_0L_r_rad)/(b/2);
coeff_b_0L = alpha_0L_r_rad;
epsilon_G_0L_r = 0*u.deg;
epsilon_G_0L_r_deg = epsilon_G_0L_r/u.deg;
epsilon_G_0L_r_rad = epsilon_G_0L_r/u.rad;
epsilon_G_0L_t = -2.5*u.deg;
epsilon_G_0L_t_deg = epsilon_G_0L_t/u.deg;
epsilon_G_0L_t_rad = epsilon_G_0L_t/u.rad;
coeff_a_eps = (epsilon_G_0L_t_rad-epsilon_G_0L_r_rad)/(b/2);
coeff_b_eps = epsilon_G_0L_r_rad; 

fun=@(y) ((((coeff_a_0L*u.m)*y+coeff_b_0L)-(coeff_a_eps*u.m)*y).*((coeff_a_c)*y+coeff_b_c/u.m));
alpha_0L =( 2/(wing_surface/u.m^2))* integral(fun,0,(b/2)/u.m);
alpha_0L_rad = alpha_0L*u.rad;
alpha_0L_deg = alpha_0L_rad/u.deg;
 

coeff_e =  pi*(coeff_b_c*alpha_0L_rad-coeff_b_0L*coeff_b_c)*(X_ac-0.25*coeff_b_c) ;
coeff_f =  pi*(coeff_a_c*alpha_0L_rad-coeff_b_0L*coeff_a_c-coeff_a_0L*coeff_b_c+coeff_b_c*coeff_a_eps)*(X_ac-0.25*coeff_b_c)+...
    pi*(-tan(lambda_LE_rad)-0.25*coeff_a_c)*(coeff_b_c*alpha_0L_rad-coeff_b_0L*coeff_b_c);
coeff_g = pi*(-coeff_a_c*coeff_a_0L+coeff_a_c*coeff_a_eps)*(X_ac-0.25*coeff_b_c)+...
    pi*(coeff_a_c*alpha_0L_rad-coeff_b_0L*coeff_a_c-coeff_a_0L*coeff_b_c+coeff_b_c*coeff_a_eps)*(-tan(lambda_LE_rad)-0.25*coeff_a_c);
 
coeff_h = pi*(-coeff_a_c*coeff_a_0L+coeff_a_c*coeff_a_eps)*(-tan(lambda_LE_rad)-0.25*coeff_a_c);
fun=@(y) (coeff_e/u.m^2+coeff_f/u.m*y+coeff_g*y.^2 + coeff_h*u.m* y.^3);
cm_ac_b = (2/((mean_chord/u.m)*(wing_surface/u.m^2)))* integral(fun,0,(b/2)/u.m);
 
cm_ac = cm_ac_a + cm_ac_b;
fun=@(y) (((coeff_a_cl*u.m)*y+coeff_b_cl).*((coeff_a_c)*y+coeff_b_c/u.m));
cl_alpha_med =( 2/(wing_surface/u.m^2))* integral(fun,0,(b/2)/u.m);
cl_alpha_med_rad = cl_alpha_med*u.rad;
cl_alpha_med_deg = cl_alpha_med*u.deg;
coeff_i = 0.5*(coeff_b_c*coeff_b_cl)*(alpha_0L_rad-coeff_b_0L);
coeff_j = 0.5*(coeff_a_c*coeff_b_cl+coeff_a_cl*coeff_b_c)*(alpha_0L_rad-coeff_b_0L)+0.5*(coeff_b_c*coeff_b_cl)*(coeff_a_eps-coeff_a_0L);
coeff_k = 0.5*( coeff_a_cl*coeff_a_c)*(alpha_0L_rad-coeff_b_0L)+0.5*(coeff_a_c*coeff_b_cl+coeff_a_cl*coeff_b_c)*(coeff_a_eps-coeff_a_0L); 
coeff_l = 0.5*( coeff_a_cl*coeff_a_c)*(coeff_a_eps-coeff_a_0L);

 
 
%% Write data file
[status, msg] = mkdir("./pitching_moment"); % create folder first
fid = fopen('./pitching_moment/data.tex', 'w');
if (fid == -1)
    fprintf("Cannot open file.")
else
    fprintf(fid, "\\def\\mySweepLEWingDEG{%f}\n", lambda_LE_deg);
    fprintf(fid, "\\def\\mySweepLEWingRAD{%f}\n", lambda_LE_rad);
    fprintf(fid, "\\def\\myChordRootWingMT{%f}\n", c_r);
    fprintf(fid, "\\def\\myChordTipWingMT{%f}\n", c_t);
    fprintf(fid, "\\def\\mySpanWingMT{%f}\n", b);
    fprintf(fid, "\\def\\myTaperRatioWing{%f}\n", taper_ratio); 
    fprintf(fid, "\\def\\myAreaWingMTsquared{%f}\n", wing_surface_mt_squared);
    fprintf(fid, "\\def\\myAspectRatioWing{%f}\n", aspect_ratio);
    fprintf(fid, "\\def\\myMach{%f}\n", mach);
    fprintf(fid, "\\def\\myCLAlphaRootWingRAD{%f}\n", cl_alpha_r_rad);
    fprintf(fid, "\\def\\myCLAlphaRootWingDEG{%f}\n", cl_alpha_r_deg);     
    fprintf(fid, "\\def\\myCLAlphaTipWingRAD{%f}\n", cl_alpha_t_rad);
    fprintf(fid, "\\def\\myCLAlphaTipWingDEG{%f}\n", cl_alpha_t_deg);
    fprintf(fid, "\\def\\myAlphaZeroLiftRootWingDEG{%f}\n", alpha_0L_r_deg);
    fprintf(fid, "\\def\\myAlphaZeroLiftTipWingDEG{%f}\n", alpha_0L_t_deg);
    fprintf(fid, "\\def\\myAlphaZeroLiftRootWingRAD{%f}\n", alpha_0L_r_rad);
    fprintf(fid, "\\def\\myAlphaZeroLiftTipWingRAD{%f}\n", alpha_0L_t_rad);
    fprintf(fid, "\\def\\myXsiacRootWing{%f}\n", x_ac_2D_r);
    fprintf(fid, "\\def\\myXsiacTipWing{%f}\n", x_ac_2D_t);
    fprintf(fid, "\\def\\myCmZeroRootWing{%f}\n", cm_ac_r);
    fprintf(fid, "\\def\\myCmZeroTipWing{%f}\n", cm_ac_t);
    fprintf(fid, "\\def\\mySweepQuarterChordWingDEG{%f}\n", lambda_C4_deg);
    fprintf(fid, "\\def\\mySweepQuarterChordWingRAD{%f}\n", lambda_C4_rad);
    fprintf(fid, "\\def\\myMACWingMT{%f}\n", mean_chord);
    fprintf(fid, "\\def\\myYMACWingMT{%f}\n", y_mean_chord);
    fprintf(fid, "\\def\\myXMACLEToApexWingMT{%f}\n", x_le_mean_chord);
    fprintf(fid, "\\def\\myKOneACDatcomWing{%f}\n", k1);   
    fprintf(fid, "\\def\\myKTwoACDatcomWing{%f}\n", k2);
    fprintf(fid, "\\def\\myXACOverChordRootDatcomWing{%f}\n", x_ac_cr);
    fprintf(fid, "\\def\\myXsiACWing{%f}\n", x_ac_mean);
    fprintf(fid, "\\def\\myACWingToApexWingMT{%f}\n", X_ac);
    fprintf(fid, "\\def\\myACWingToMACLEWingMT{%f}\n", x_ac_mean*mean_chord);
    fprintf(fid, "\\def\\myCoeffAChordWing{%f}\n", coeff_a_c);
    fprintf(fid, "\\def\\myCoeffBChordWingMT{%f}\n", coeff_b_c);
    fprintf(fid, "\\def\\myCoeffAAeroTwistWingRADMT{%f}\n", coeff_a_0L);
    fprintf(fid, "\\def\\myCoeffACmZeroWingMT{%f}\n", coeff_a_cm);
    fprintf(fid, "\\def\\myCoeffBCmZeroWing{%f}\n", coeff_b_cm);
    
    fprintf(fid, "\\def\\myCmZeroAWing{%f}\n", cm_ac_a);
    fprintf(fid, "\\def\\myCmZeroBWing{%f}\n", cm_ac_b);
    fprintf(fid, "\\def\\myCmZeroWing{%f}\n", cm_ac);
    fprintf(fid, "\\def\\myCoeffAWing{%f}\n", coeff_a);
    fprintf(fid, "\\def\\myCoeffBWing{%f}\n", coeff_b);
    fprintf(fid, "\\def\\myCoeffCWing{%f}\n", coeff_c);
    fprintf(fid, "\\def\\myCoeffDWing{%f}\n", coeff_d);     
    fprintf(fid, "\\def\\myCoeffEWing{%f}\n", coeff_e);     
    fprintf(fid, "\\def\\myCoeffFWing{%f}\n", coeff_f);     
    fprintf(fid, "\\def\\myCoeffHWing{%f}\n", coeff_h);
    fprintf(fid, "\\def\\myCoeffIWing{%f}\n", coeff_i);
    fprintf(fid, "\\def\\myCoeffJWing{%f}\n", coeff_j);
    fprintf(fid, "\\def\\myCoeffKWing{%f}\n", coeff_k);
    fprintf(fid, "\\def\\myCoeffLWing{%f}\n", coeff_l);
    fprintf(fid, "\\def\\myCLAlphaMeanWingRAD{%f}\n", cl_alpha_med_rad);
    fprintf(fid, "\\def\\myCoeffAClalphaWingRADMT{%f}\n", coeff_a_cl);
    fprintf(fid, "\\def\\myCoeffBClalphaWingRAD{%f}\n", coeff_b_cl);
    fprintf(fid, "\\def\\myTwistWingDEG{%f}\n", epsilon_G_0L_t_deg);
    fprintf(fid, "\\def\\myTwistWingRAD{%f}\n", epsilon_G_0L_t_rad);
    fprintf(fid, "\\def\\myCoeffATwistWingRADMT{%f}\n", coeff_a_eps);
    fprintf(fid, "\\def\\myCoeffBAeroTwistWingRAD{%f}\n", coeff_b_0L);
    fprintf(fid, "\\def\\myAlphaZeroLiftWingRAD{%f}\n", alpha_0L_rad );
    fprintf(fid, "\\def\\myAlphaZeroLiftWingDEG{%f}\n", alpha_0L_deg);
    % ...
    fclose(fid);
end

 



 
   
 