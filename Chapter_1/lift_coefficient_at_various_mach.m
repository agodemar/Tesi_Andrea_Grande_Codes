clear all;close all; clc;
altitude_mt= 4000 * u.m;
altitude_ft=altitude_mt/u.ft;
mass_kg=79000 * u.kg;
mass_lb=mass_kg/u.lb;
wing_surface_mt_squared= 124.6 * u.m^2;
wing_surface_ft_squared=wing_surface_mt_squared/u.ft^3;
mach_number = [0.3;0.4;0.5;0.6;0.7;0.8];
[T, a, P, rho] = atmosisa(altitude_mt/u.m);
sound_speed = a * u.meterPerSecond;
speed_ms=mach_number*a ;
density_kg_mt_cubed = rho * (u.kg/u.m3);
gravitational_acceleration=9.81*(u.meterPerSecond/u.second);
dynamic_pressure_pa=0.5.*density_kg_mt_cubed.*(speed_ms).^2;
dynamic_pressure_bar=dynamic_pressure_pa/u.bar;
lift_coefficient=(gravitational_acceleration*mass_kg)./(dynamic_pressure_pa*wing_surface_mt_squared)* u.s^2/u.m^2; 
x = mach_number';
y = lift_coefficient';

%writing data file for pgfplot
data_M_CL = [x', y'];
save("data_M_CL.txt", 'data_M_CL', '-ascii');

%% Write data file
[status, msg] = mkdir("./lift_coefficient_at_various_mach"); % create folder first
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
    fprintf(fid, "\\def\\myMach{%f}\n", mach_number);
    fprintf(fid, "\\def\\mySoundSpeed{%f}\n", sound_speed);
    fprintf(fid, "\\def\\mySpeedMS{%f}\n", speed_ms);
    fprintf(fid, "\\def\\myGravitationalAcceleration{%f}\n", gravitational_acceleration);
    fprintf(fid, "\\def\\myDensityKGMC{%f}\n", density_kg_mt_cubed);
    fprintf(fid, "\\def\\myDynamicPressurePa{%f}\n", dynamic_pressure_pa);
    fprintf(fid, "\\def\\myDynamicPressureBar{%f}\n", dynamic_pressure_bar);
    fprintf(fid, "\\def\\myLiftCoefficient{%f}\n", lift_coefficient);




    % ...
    fclose(fid);
end
