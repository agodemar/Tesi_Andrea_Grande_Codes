clearvars; close all; clc

% see https://it.mathworks.com/matlabcentral/fileexchange/38977-physical-units-toolbox

speedEquivalent = 350 * str2u('km/h');

fprintf("Equivalent airspeed in m/s: %f\n", speedEquivalent / u.meterPerSecond)

altitude = 5000 * u.ft;

[T, a, P, rho] = atmosisa(altitude/u.m);

temperature = T * u.kelvin;
soundSpeed = a * u.meterPerSecond;
pressure = P * u.newton/u.m2;
density = rho * u.kg/u.m3;

fprintf("Air density in kg/m3: %f\n", density / (u.kg/u.m3))
fprintf("Air density in slug/ft3: %f\n", density / (u.slug/u.ft3))

%% Write data file
[status, msg] = mkdir("./equivalent_airspeed_basic_2"); % create folder first
fid = fopen('./equivalent_airspeed_basic_2/data.tex', 'w');
if (fid == -1)
    fprintf("Cannot open file.")
else
    fprintf(fid, "\\def\\myAltitudeMT{%f}\n", altitude/u.m);
    fprintf(fid, "\\def\\myAltitudeFT{%f}\n", altitude/u.ft);
    % ...
    fclose(fid);
end
