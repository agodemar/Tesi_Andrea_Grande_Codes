clear all;close all; clc;
altitude_mt= 9000 * u.m;
altitude_ft=altitude_mt/u.ft;
mass_kg=79000 * u.kg;
mass_lb=mass_kg/u.lb;
wing_surface_mt_squared= 124.6 * u.m^2;
wing_surface_ft_squared=wing_surface_mt_squared/u.ft^3;
speed_kts=200 * u.knot;
speed_kmh=speed_kts/u.kmh;
speed_ms= speed_kts/u.meterPerSecond;
density_kg_mt_cubed=0.466 * u.kg/u.m3; 
mach_number = 0.70;
gravitational_acceleration=9.81*(u.meterPerSecond/u.second);
dynamic_pressure_pa=0.5*density_kg_mt_cubed*(speed_ms)^2;
dynamic_pressure_bar=dynamic_pressure_pa/u.bar;
lift_coefficient=(gravitational_acceleration*mass_kg)/(dynamic_pressure_pa*wing_surface_mt_squared) * u.s^2/u.m^2;
aspect_ratio_wing=9.46;
oswald_factor=0.85;
k_polar=0.0124;
drag_coefficient_zero=0.0190;
drag_coefficient=drag_coefficient_zero + k_polar*(lift_coefficient)^2;
drag_N=drag_coefficient*dynamic_pressure_pa*wing_surface_mt_squared;
drag_kgf=drag_N/u.kgf;
drag_lbf=drag_N/u.lbf;
%% Write data file
[status, msg] = mkdir("./aircraft_aerodynamic_coefficients"); % create folder first
fid = fopen('./aircraft_aerodynamic_coefficients/data.tex', 'w');
if (fid == -1)
    fprintf("Cannot open file.")
else
    fprintf(fid, "\\def\\myAltitudeM{%f}\n", altitude_mt);
    fprintf(fid, "\\def\\myAltitudeFt{%f}\n", altitude_ft);
    fprintf(fid, "\\def\\myMassKg{%f}\n", mass_kg);
    fprintf(fid, "\\def\\myMassLb{%f}\n", mass_lb);
    fprintf(fid, "\\def\\myWingSurfaceMTSD{%f}\n", wing_surface_mt_squared);
    fprintf(fid, "\\def\\myWingSurfaceFTSD{%f}\n", wing_surface_ft_squared);
    fprintf(fid, "\\def\\myMach{%f}\n", mach_number);
    fprintf(fid, "\\def\\mySpeedKts{%f}\n", speed_kts);
    fprintf(fid, "\\def\\mySpeedMS{%f}\n", speed_ms);
    fprintf(fid, "\\def\\mySpeedKmH{%f}\n", speed_kmh);
    fprintf(fid, "\\def\\myGravitationalAcceleration{%f}\n", gravitational_acceleration);
    fprintf(fid, "\\def\\myDensityKGMC{%f}\n", density_kg_mt_cubed);
    fprintf(fid, "\\def\\myDynamicPressurePa{%f}\n", dynamic_pressure_pa);
    fprintf(fid, "\\def\\myDynamicPressureBar{%f}\n", dynamic_pressure_bar);
    fprintf(fid, "\\def\\myLiftCoefficient{%f}\n", lift_coefficient);
    fprintf(fid, "\\def\\myAspectRatio{%f}\n", aspect_ratio_wing);
    fprintf(fid, "\\def\\myOswaldFactor{%f}\n", oswald_factor);
    fprintf(fid, "\\def\\myKpolar{%f}\n", k_polar);
    fprintf(fid, "\\def\\myDragCoefficientZero{%f}\n", drag_coefficient_zero);
    fprintf(fid, "\\def\\myDragN{%f}\n", drag_N);
    fprintf(fid, "\\def\\myDragKgf{%f}\n", drag_kgf);
    fprintf(fid, "\\def\\myDragLbf{%f}\n", drag_lbf);









    % ...
    fclose(fid);
end
