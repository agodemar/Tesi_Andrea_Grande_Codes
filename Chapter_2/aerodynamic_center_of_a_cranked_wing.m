clear all;close all; clc;

b = 58.67 * u.m;
b_1 = 19.36 * u.m;
b_2 = 39.31 *u.m;
c_r_1 = 11.87 * u.m;
c_r = c_r_1;
c_t_1 = 6.42 * u.m;
c_r_2 = c_t_1;
c_t_2 = 1.69 * u.m;
c_t=c_t_2;
lambda_LE_1 = 32.2 *u.deg;
lambda_LE_1_rad = lambda_LE_1 /u.rad;
lambda_LE_1_deg = lambda_LE_1 /u.deg;
lambda_TE_1 = 3.8 *u.deg;
lambda_TE_1_rad = lambda_TE_1 /u.rad;
lambda_TE_1_deg = lambda_TE_1 /u.deg;
cl_alpha_r_1 = 6.15*(1/u.rad);
cl_alpha_r_1_deg = cl_alpha_r_1*u.deg;
cl_alpha_r_1_rad = cl_alpha_r_1*u.rad;
cl_alpha_t_1 = 6.05*(1*u.rad);
cl_alpha_t_1_deg = cl_alpha_t_1*u.deg;
cl_alpha_t_1_rad = cl_alpha_t_1*u.rad;
coeff_a_cl_1 = (cl_alpha_t_1_rad-cl_alpha_r_1_rad)/(b_1/2);
coeff_a_cl_1_radm = coeff_a_cl_1*u.rad;
coeff_a_cl_1_degm = coeff_a_cl_1*u.deg;
lambda_LE_2 = 32.2 *u.deg;
lambda_LE_2_rad = lambda_LE_2 /u.rad;
lambda_LE_2_deg = lambda_LE_2 /u.deg;
lambda_TE_2 = 21.3 *u.deg;
lambda_TE_2_rad = lambda_TE_2 /u.rad;
lambda_TE_2_deg = lambda_TE_2 /u.deg;
cl_alpha_r_2 = 6.05*(1/u.rad);
cl_alpha_r_2_deg = cl_alpha_r_2*u.deg;
cl_alpha_r_2_rad = cl_alpha_r_2*u.rad;
cl_alpha_t_2 = 6.01*(1*u.rad);
cl_alpha_t_2_deg = cl_alpha_t_2*u.deg;
cl_alpha_t_2_rad = cl_alpha_t_2*u.rad;
coeff_a_cl_2 = (cl_alpha_t_2_rad-cl_alpha_r_2_rad)/(b_2/2);
coeff_a_cl_2_radm = coeff_a_cl_2*u.rad;
coeff_a_cl_2_degm = coeff_a_cl_2*u.deg;
x_ac_2D_1_r = 0.25;
x_ac_2D_1_t = 0.25;
cm_ac_r_1 = -0.080;
cm_ac_t_1 = -0.080;
x_ac_2D_2_r = 0.25;
x_ac_2D_2_t = 0.25;
cm_ac_r_2 = -0.080;
cm_ac_t_2 = -0.040;
alpha_0L_r_1 = -2.5*u.deg;
alpha_0L_r_1_deg = alpha_0L_r_1/u.deg;
alpha_0L_r_1_rad = alpha_0L_r_1/u.rad;
alpha_0L_r_2  = alpha_0L_r_1;
alpha_0L_r_2_deg = alpha_0L_r_2/u.deg;
alpha_0L_r_2_rad =alpha_0L_r_2/u.rad;
alpha_0L_t_1 = -2.5*u.deg;
alpha_0L_t_1_deg = alpha_0L_t_1/u.deg;
alpha_0L_t_1_rad = alpha_0L_t_1/u.rad;
alpha_0L_t_2 = -1*u.deg;
alpha_0L_t_2_deg = alpha_0L_t_2/u.deg;
alpha_0L_t_2_rad = alpha_0L_t_2/u.rad;
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
mach = 0.65;
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
y_mean_chord_2 = (b_1/2)+(b_2/6)*((1+2*taper_ratio_2)/(1+taper_ratio_2));
wing_surface = wing_surface_1+wing_surface_2;
wing_surface_mt_squared = wing_surface*u.m2;
aspect_ratio = b^2 / wing_surface ;
mean_chord = (wing_surface_1*mean_chord_1+wing_surface_2*mean_chord_2)/(wing_surface); 
y_mean_chord = (mean_chord - coeff_b_1)/coeff_a_1;
x_le_mean_chord = y_mean_chord*tan(lambda_LE_1_rad);
k1_1 = 1.259;
k2_1 = 0.242;
x_ac_cr_1 = 0.411;
x_ac_c_1 = k1_1*(x_ac_cr_1-k2_1);
b_2_datcom = 2*(b_2/2 + b_1/4);
c_r_2_datcom = -coeff_a_2*b_1/4 + coeff_b_2;
taper_ratio_2_datcom = c_t_2/c_r_2_datcom;
wing_surface_2_datcom = (b_2_datcom/2)*c_r_2_datcom*(1+taper_ratio_2_datcom);
wing_surface_2_datcom_mt_squared = wing_surface_2_datcom/u.m2;
aspect_ratio_2_datcom = b_2_datcom^2 / wing_surface_2_datcom ;
k1_2_datcom = 1.435;
k2_2_datcom = 0.820;
x_ac_cr_2_datcom = 0.958;
a0 = 6.088/(1-mach*mach*cos(lambda_LE_1_rad)*cos(lambda_LE_1_rad))^(1/2);
cl_alpha_1 = a0*cos(lambda_LE_1_rad)/((1-mach*mach*cos(lambda_LE_1_rad)*cos(lambda_LE_1_rad)+ ...
(a0*cos(lambda_LE_1_rad)/(pi*aspect_ratio_1))^2)^(1/2) +(a0*cos(lambda_LE_1_rad)/(pi*aspect_ratio_1) ));
cl_alpha_1_deg = cl_alpha_1*u.deg; 
cl_alpha_1_rad = cl_alpha_1*u.rad;
 cl_alpha_2_datcom = a0*cos(lambda_LE_2_rad)/((1-mach*mach*cos(lambda_LE_2_rad)*cos(lambda_LE_2_rad)+ ...
(a0*cos(lambda_LE_2_rad)/(pi*aspect_ratio_2_datcom))^2)^(1/2) +(a0*cos(lambda_LE_2_rad)/(pi*aspect_ratio_2_datcom) ));
cl_alpha_2_datcom_deg = cl_alpha_2_datcom*u.deg;
cl_alpha_2_datcom_rad = cl_alpha_2_datcom*u.rad;
X_ac_cr = (x_ac_cr_1*wing_surface_1*cl_alpha_1_rad+x_ac_cr_2_datcom*wing_surface_2_datcom*cl_alpha_2_datcom_rad)/(wing_surface_1*cl_alpha_1_rad+wing_surface_2_datcom*cl_alpha_2_datcom_rad);
X_ac = X_ac_cr *c_r_1;
x_ac_c = (X_ac-  x_le_mean_chord) / mean_chord;

%% Write data file
[status, msg] = mkdir("./aerodynamic_center_of_a_cranked_wing"); % create folder first
fid = fopen('./aerodynamic_center_of_a_cranked_wing/data.tex', 'w');
if (fid == -1)
    fprintf("Cannot open file.")
else
    fprintf(fid, "\\def\\mySweepLEWingIDEG{%f}\n", lambda_LE_1_deg);
    fprintf(fid, "\\def\\mySweepLEWingIIDEG{%f}\n", lambda_LE_2_deg);
    fprintf(fid, "\\def\\mySweepLEWingIRAD{%f}\n", lambda_LE_1_rad);
    fprintf(fid, "\\def\\mySweepLEWingIIRAD{%f}\n", lambda_LE_2_rad);
    fprintf(fid, "\\def\\mySweepTEWingIDEG{%f}\n", lambda_TE_1_deg);
    fprintf(fid, "\\def\\mySweepTEWingIIDEG{%f}\n", lambda_TE_2_deg);
    fprintf(fid, "\\def\\mySweepTEWingIRAD{%f}\n", lambda_TE_1_rad);
    fprintf(fid, "\\def\\mySweepTEWingIIRAD{%f}\n", lambda_TE_2_rad);
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
    fprintf(fid, "\\def\\myAlphaZeroLiftRootWingIDEG{%f}\n", alpha_0L_r_1_deg);
    fprintf(fid, "\\def\\myAlphaZeroLiftTipWingIDEG{%f}\n", alpha_0L_t_1_deg);
    fprintf(fid, "\\def\\myAlphaZeroLiftRootWingIRAD{%f}\n", alpha_0L_r_1_rad);
    fprintf(fid, "\\def\\myAlphaZeroLiftTipWingIRAD{%f}\n", alpha_0L_t_1_rad);
    fprintf(fid, "\\def\\myAlphaZeroLiftRootWingIIDEG{%f}\n", alpha_0L_r_2_deg);
    fprintf(fid, "\\def\\myAlphaZeroLiftTipWingIIDEG{%f}\n", alpha_0L_t_2_deg);
    fprintf(fid, "\\def\\myAlphaZeroLiftRootWingIIRAD{%f}\n", alpha_0L_r_2_rad);
    fprintf(fid, "\\def\\myAlphaZeroLiftTipWingIIRAD{%f}\n", alpha_0L_t_2_rad);
    
    fprintf(fid, "\\def\\myTaperRatioWingI{%f}\n", taper_ratio_1);
    fprintf(fid, "\\def\\myTaperRatioWingII{%f}\n", taper_ratio_2);
    fprintf(fid, "\\def\\myTwistWingIDEG{%f}\n", epsilon_G_0L_t_1_deg);
    fprintf(fid, "\\def\\myTwistWingIRAD{%f}\n", epsilon_G_0L_t_1_rad);
    fprintf(fid, "\\def\\myTwistWingIIDEG{%f}\n", epsilon_G_0L_t_2_deg);
    fprintf(fid, "\\def\\myTwistWingIIRAD{%f}\n", epsilon_G_0L_t_2_rad);
    fprintf(fid, "\\def\\myAreaWingIMTsquared{%f}\n", wing_surface_1_mt_squared);
    fprintf(fid, "\\def\\myAreaWingIIMTsquared{%f}\n", wing_surface_2_mt_squared);
       
    fprintf(fid, "\\def\\myCLAlphaRootWingIRAD{%f}\n",  cl_alpha_r_1_rad );
    fprintf(fid, "\\def\\myCLAlphaRootWingIDEG{%f}\n",  cl_alpha_r_1_deg );
    fprintf(fid, "\\def\\myCLAlphaRootWingIIRAD{%f}\n", cl_alpha_r_2_rad    );
    fprintf(fid, "\\def\\myCLAlphaRootWingIIDEG{%f}\n", cl_alpha_r_2_deg ); 
    fprintf(fid, "\\def\\myCLAlphaTipWingIRAD{%f}\n", cl_alpha_t_1_rad   );
    fprintf(fid, "\\def\\myCLAlphaTipWingIDEG{%f}\n", cl_alpha_t_1_deg ); 
    fprintf(fid, "\\def\\myCLAlphaTipWingIIRAD{%f}\n",  cl_alpha_t_2_rad );
    fprintf(fid, "\\def\\myCLAlphaTipWingIIDEG{%f}\n", cl_alpha_t_2_deg ); 
    
    fprintf(fid, "\\def\\myXsiacRootWingI{%f}\n", x_ac_2D_1_r ); 
    fprintf(fid, "\\def\\myXsiacTipWingI{%f}\n",x_ac_2D_1_t ); 
  fprintf(fid, "\\def\\myXsiacRootWingII{%f}\n", x_ac_2D_2_r ); 
  fprintf(fid, "\\def\\myXsiacTipWingII{%f}\n",x_ac_2D_2_t ); 
  fprintf(fid, "\\def\\myCmZeroRootWingI{%f}\n",cm_ac_r_1 ); 
  fprintf(fid, "\\def\\myCmZeroTipWingI{%f}\n",cm_ac_t_1 ); 
  fprintf(fid, "\\def\\myCmZeroRootWingII{%f}\n",cm_ac_r_2 ); 
  fprintf(fid, "\\def\\myCmZeroTipWingII{%f}\n",cm_ac_t_2 ); 
  fprintf(fid, "\\def\\myMach{%f}\n",mach ); 
  fprintf(fid, "\\def\\myMACWingIMT{%f}\n", mean_chord_1 ); 
  fprintf(fid, "\\def\\myMACWingIIMT{%f}\n",mean_chord_2 ); 
  fprintf(fid, "\\def\\myMACWingCrankedMT{%f}\n",mean_chord ); 
  fprintf(fid, "\\def\\myAreaWingCrankedMTsquared{%f}\n", wing_surface ); 
  fprintf(fid, "\\def\\myAspectRatioWingCranked{%f}\n",aspect_ratio ); 
  fprintf(fid, "\\def\\myXMACLEToApexWingIMT{%f}\n",x_le_mean_chord_1 ); 
  fprintf(fid, "\\def\\myAspectRatioWingI{%f}\n",aspect_ratio_1 ); 
  fprintf(fid, "\\def\\myAspectRatioWingII{%f}\n",aspect_ratio_2 ); 

  fprintf(fid, "\\def\\myYMACWingIMT{%f}\n",y_mean_chord_1 ); 
  fprintf(fid, "\\def\\myXMACLEToApexWingIIMT{%f}\n",x_le_mean_chord_2 ); 
  fprintf(fid, "\\def\\myYMACWingIIMT{%f}\n",y_mean_chord_2 ); 

  fprintf(fid, "\\def\\myYYMACWingCrankedMT{%f}\n",y_mean_chord ); 
  fprintf(fid, "\\def\\myXXMACLEToApexWingCrankedMT{%f}\n",x_le_mean_chord ); 
  fprintf(fid, "\\def\\myKOneACDatcomWingI{%f}\n",k1_1 );
  fprintf(fid, "\\def\\myKTwoACDatcomWingI{%f}\n",k2_1 ); 
  fprintf(fid, "\\def\\myXACOverChordRootDatcomWingI{%f}\n",x_ac_cr_1 ); 
  fprintf(fid, "\\def\\myXACOverChordRootDatcomWingII{%f}\n",x_ac_cr_2_datcom );
   
    fprintf(fid, "\\def\\myXACOverChordRootDatcomWing{%f}\n",X_ac_cr );

   fprintf(fid, "\\def\\myChordRootWingIIPrimeMT{%f}\n",c_r_2_datcom ); 
   fprintf(fid, "\\def\\myTaperRatioWingIIPrime{%f}\n",taper_ratio_2_datcom );
  fprintf(fid, "\\def\\mySpanWingIIPrimeMT{%f}\n", b_2_datcom); 
  fprintf(fid, "\\def\\myAreaWingIIPrimeMTsquared{%f}\n",wing_surface_2_datcom); 
  fprintf(fid, "\\def\\myAspectRatioWingIIPrime{%f}\n",aspect_ratio_2_datcom );
  fprintf(fid, "\\def\\myCLAlphaAtMACWingRAD{%f}\n",a0 );
  fprintf(fid, "\\def\\myCLAlphaPolhamusWingIRAD{%f}\n",cl_alpha_1_rad );
  fprintf(fid, "\\def\\myCLAlphaPolhamusWingIDEG{%f}\n",cl_alpha_1_deg);
  fprintf(fid, "\\def\\myCLAlphaPolhamusWingIIRAD{%f}\n",cl_alpha_2_datcom_rad );
  fprintf(fid, "\\def\\myCLAlphaPolhamusWingIIDEG{%f}\n",cl_alpha_2_datcom_deg );
  fprintf(fid, "\\def\\myXsiACWing{%f}\n",x_ac_c );
  fprintf(fid, "\\def\\myACWingToApexWingMT{%f}\n",X_ac );
  fprintf(fid, "\\def\\myXsiACWingI{%f}\n",x_ac_c_1);
 
  fprintf(fid, "\\def\\myKOneACDatcomWingII{%f}\n",k1_2_datcom);
  fprintf(fid, "\\def\\myKTwoACDatcomWingII{%f}\n",k2_2_datcom);
 
    % ...
    fclose(fid);
end  


 
 
 









 




 


 

 

