clear all ; close all ; clc ;
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
X_P_1 = 1.37 * u.m;
X_P_2 = 6.02 * u.m;
c_r_eq = X_P_2-X_P_1;
taper_ratio_eq = c_t/c_r_eq;
lambda_LE_E = atan((5.89*u.m -X_P_1)/(b/2)); 
lambda_LE_E_deg  = lambda_LE_E/u.rad;
lambda_LE_E_rad =  lambda_LE_E/u.deg;


  
%% Write data file
[status, msg] = mkdir("./equivalent_plan_form_of_a_cranked_wing"); % create folder first
fid = fopen('./equivalent_plan_form_of_a_cranked_wing/data.tex', 'w');
if (fid == -1)
    fprintf("Cannot open file.")
else
    fprintf(fid, "\\def\\mySpanWingMT{%f}\n", b);
    fprintf(fid, "\\def\\mySpanWingIMT{%f}\n", b_1);
    fprintf(fid, "\\def\\mySpanWingIIMT{%f}\n",b_2);

    fprintf(fid, "\\def\\myRootChordWingLongitudinalPlaneMT{%f}\n", c_r_eq);
    fprintf(fid, "\\def\\myChordRootWingIMT{%f}\n", c_r_1);
    fprintf(fid, "\\def\\myChordRootWingIIMT{%f}\n", c_r_2);

    fprintf(fid, "\\def\\myChordTipWingMT{%f}\n", c_t);
    fprintf(fid, "\\def\\myChordTipWingIMT{%f}\n", c_t_1);
    fprintf(fid, "\\def\\myChordTipWingIIMT{%f}\n", c_t_2);
    
    fprintf(fid, "\\def\\mySweepLEWingIDEG{%f}\n", lambda_LE_1_deg);
    fprintf(fid, "\\def\\mySweepLEWingIIDEG{%f}\n", lambda_LE_2_deg);
    fprintf(fid, "\\def\\mySweepLEWingIRAD{%f}\n", lambda_LE_1_rad);
    fprintf(fid, "\\def\\mySweepLEWingIIRAD{%f}\n", lambda_LE_2_rad);
    fprintf(fid, "\\def\\myEquivalentSweepLEWingRAD{%f}\n", lambda_LE_E_rad);
    fprintf(fid, "\\def\\myEquivalentSweepLEWingDEG{%f}\n", lambda_LE_E_deg);

    fprintf(fid, "\\def\\myXEquivalentChordLEToApexWingMT{%f}\n", X_P_1);
    fprintf(fid, "\\def\\myXEquivalentChordTEToApexWingMT{%f}\n",X_P_2);
    
    fprintf(fid, "\\def\\myEquivalentTaperRatioWing{%f}\n", taper_ratio_eq);
    


    % ...
    fclose(fid);
end


  