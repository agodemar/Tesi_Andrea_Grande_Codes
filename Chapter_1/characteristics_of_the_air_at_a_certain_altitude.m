clear all;close all; clc;
altitude = 8000 * u.m;
temperatureSeaLevel=288.16*u.K;
temperatureGradient=-0.0065*(u.K/u.m);
adiabaticIndex=1.4;
gravitationalacceleration=9.81*(u.meterPerSecond/u.second);
densitysealevel=1.225*u.kg/u.m3;
airSpecificConstant= 287*(u.N*u.m)/(u.kg*u.K);
mach_number = 0.70;

[T, a, P, rho] = atmosisa(altitude/u.m);

temperature = T * u.kelvin;
soundSpeed = a * u.meterPerSecond;
pressure = P * (u.newton/u.m2);
density = rho * (u.kg/u.m3);
speed=soundSpeed*mach_number;
relativedensity= density/densitysealevel;
dynamicpressure=0.5*density*(speed^2);
c1=1.458*u.kg / (u.meterPerSecond* (u.K)^(1/2)); %coeff legge di sutherland
c2=110.4*u.K;
dynamicviscosity=c1*temperature^(3/2)/(temperature+c2);  
reynoldsnumeberperunitoflenght=density*speed/dynamicviscosity;

fprintf("Temperature at 8000m in kelvin : %f\n", temperature )
fprintf("soundspeed at 8000m in m/s : %f\n", soundSpeed )
fprintf("speed at 8000m in m/s  : %f\n", speed )
fprintf("density at 8000m in kg/m^3 : %f\n", density )
fprintf("dynamic pressure at 8000m in N/m^2 : %f\n", dynamicpressure )
fprintf("dynamic viscosity at 8000m in kg/ms : %f\n", dynamicviscosity )
fprintf("reynolds number per unit of lenght at 8000m : %f\n",reynoldsnumeberperunitoflenght )









%% Write data file
[status, msg] = mkdir("./characteristics_of_the_air_at_a_certain_altitude"); % create folder first
fid = fopen('./characteristics_of_the_air_at_a_certain_altitude/data.tex', 'w');
if (fid == -1)
    fprintf("Cannot open file.")
else
    fprintf(fid, "\\def\\myAdiabaticindex{%f}\n", adiabaticIndex);
    fprintf(fid, "\\def\\myAltitudeM{%f}\n", altitude);
    fprintf(fid, "\\def\\myTemperaturesealevelK{%f}\n", temperatureSeaLevel);
    fprintf(fid, "\\def\\myTemperaturegradient{%f}\n", temperatureGradient);
    fprintf(fid, "\\def\\myAirspecificconstant{%f}\n", airSpecificConstant);
    fprintf(fid, "\\def\\myMachnumber{%f}\n", mach_number);
    fprintf(fid, "\\def\\myTemperatureK{%f}\n", temperature);
    fprintf(fid, "\\def\\mySoundspeed{%f}\n", soundSpeed);
    fprintf(fid, "\\def\\mySpeedMs{%f}\n", speed);
    fprintf(fid, "\\def\\mySpeedkmh{%f}\n", speed/u.kmh);
    fprintf(fid, "\\def\\myGravitationalacceleration{%f}\n", gravitationalacceleration);
    fprintf(fid, "\\def\\myDensitysealevel{%f}\n", densitysealevel);
    fprintf(fid, "\\def\\myRelativedensity{%f}\n", relativedensity);
    fprintf(fid, "\\def\\myDynamicpressureNm{%f}\n", dynamicpressure);
    fprintf(fid, "\\def\\myDynamicpressureBar{%f}\n", dynamicpressure/u.bar);
    fprintf(fid, "\\def\\myDensity{%f}\n", density);
    fprintf(fid, "\\def\\myCone{%f}\n", c1);
    fprintf(fid, "\\def\\myCtwo{%f}\n", c2);
    fprintf(fid, "\\def\\myReynoldsnumeberperunitoflenght{%f}\n", reynoldsnumeberperunitoflenght);
    fprintf(fid, "\\def\\myDynamicviscosity{%f}\n", dynamicviscosity);







    % ...
    fclose(fid);
end

