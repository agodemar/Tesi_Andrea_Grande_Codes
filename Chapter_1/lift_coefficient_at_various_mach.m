clear all; close all; clc;

altitude = 4000 * u.m;
altitude_ft = altitude/u.ft;
altitude_mt = altitude/u.m;
mass = 79000 * u.kg;
mass_kg = mass/u.kg;
mass_lb = mass/u.lb;
wing_surface = 124.6 * u.m^2;
wing_surface_mt_squared = wing_surface/u.m^2;
wing_surface_ft_squared = wing_surface/u.ft^3;
v_mach_number = [0.3; 0.4; 0.5; 0.6; 0.7; 0.8];
[T, a, P, rho] = atmosisa(altitude_mt);
sound_speed = a * u.meterPerSecond;
sound_speed_ms = a;
v_speed = v_mach_number * sound_speed ;
v_speed_ms = v_speed/u.meterPerSecond;
v_speed_kmh = v_speed/u.kmh;
density = rho * (u.kg/u.m3);
density_kg_mt_cubed = rho;
gravitational_acceleration = 9.81*(u.meterPerSecond/u.second);
v_dynamic_pressure = 0.5.*density.*(v_speed).^2;
v_dynamic_pressure_bar = v_dynamic_pressure./u.pascal;
v_lift_coefficient = (gravitational_acceleration.*mass)./(v_dynamic_pressure.*wing_surface); 

%writing data file for pgfplot
data = [v_mach_number, v_lift_coefficient, v_speed_ms, v_speed_kmh, v_dynamic_pressure_bar];

%% Write data file
[status, msg] = mkdir("./lift_coefficient_at_various_mach"); % create folder first
save("./lift_coefficient_at_various_mach/dataset.txt", 'data', '-ascii');
fid = fopen('./lift_coefficient_at_various_mach/data.tex', 'w');
if (fid == -1)
    fprintf("Cannot open file.")
else
    fprintf(fid, "\\def\\myAltitudeM{%f}\n", altitude_mt);
    fprintf(fid, "\\def\\myAltitudeFt{%f}\n", altitude_ft);
    fprintf(fid, "\\def\\myMassKg{%f}\n", mass_kg);
    fprintf(fid, "\\def\\myMassLb{%f}\n", mass_lb);
    fprintf(fid, "\\def\\myWingSurfaceMTSD{%f}\n", wing_surface_mt_squared);
    fprintf(fid, "\\def\\myWingSurfaceFTSD{%f}\n", wing_surface_ft_squared);
    fprintf(fid, "\\def\\mySoundSpeed{%f}\n", sound_speed);
    fprintf(fid, "\\def\\myGravitationalAcceleration{%f}\n", gravitational_acceleration);
    fprintf(fid, "\\def\\myDensityKGMC{%f}\n", density_kg_mt_cubed);
    % ...
    fclose(fid);
end
