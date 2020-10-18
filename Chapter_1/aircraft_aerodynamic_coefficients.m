clear all;close all; clc;
altitude = 9000 * u.m;
altitude_mt = altitude/u.m;
altitude_ft = altitude/u.ft;
mass = 79000 * u.kg;
mass_kg = mass/u.kg;
mass_lb = mass/u.lb;
wing_surface = 124.6 * u.m^2;
wing_surface_mt_squared = wing_surface /u.m^2;
wing_surface_ft_squared = wing_surface_mt_squared/u.ft^3;
speed = 350 * u.knot;
speed_kts = 350 / u.knot;
speed_kmh = speed_kts/u.kmh;
speed_ms = speed_kts/u.meterPerSecond;
[T, a, P, rho] = atmosisa(altitude_mt);
sound_speed = a * u.meterPerSecond;
sound_speed_ms = a;
density = rho* u.kg/u.m3;
density_kg_mt_cubed = rho;

mach_number = speed_ms/a;
gravitational_acceleration = 9.81*(u.meterPerSecond/u.second);
dynamic_pressure = 0.5*density*(speed)^2;
dynamic_pressure_pa = dynamic_pressure/ u.pascal;
dynamic_pressure_bar = dynamic_pressure/u.bar;
lift_coefficient = (gravitational_acceleration.*mass)/(dynamic_pressure*wing_surface) ;
aspect_ratio_wing = 9.46;
oswald_factor = 0.85;
k_polar = 1/(0.85*aspect_ratio_wing*oswald_factor);
drag_coefficient_zero = 0.0190;
drag_coefficient = drag_coefficient_zero + k_polar*(lift_coefficient)^2;
drag = drag_coefficient*dynamic_pressure*wing_surface;
drag_N = drag/u.N;
drag_kgf = drag/u.kgf;
drag_lbf = drag/u.lbf;
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
    fprintf(fid, "\\def\\myDragCoefficient{%f}\n", drag_coefficient);
    fprintf(fid, "\\def\\myDragN{%f}\n", drag_N);
    fprintf(fid, "\\def\\myDragKgf{%f}\n", drag_kgf);
    fprintf(fid, "\\def\\myDragLbf{%f}\n", drag_lbf);









    % ...
    fclose(fid);
end
