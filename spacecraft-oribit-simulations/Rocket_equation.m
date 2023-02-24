% Ideal rocket equation
% http://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node103.html

% dv = -v_exhaust * dm/m(i)
% a = dv/dt = -v_exhaust * dm/dt /m(i)
clear;
dt = 1;
Time_max = 200;
Length = Time_max / dt;

Time = 0:dt:Time_max - dt;
Y = zeros(1, Length);
V = zeros(1, Length);
M = zeros(1, Length);

% data from the legendary Saturn V rocket
% http://www.mnealon.eosc.edu/RocketSciencePage5.htm
Empty_rocket_mass = 0.136 * 10 ^ 6; 
Fuel_mass = 2.04 * 10 ^ 6;
Fuel_burn_time = 150;
V_exhaust = -2456;

dm = Fuel_mass / Fuel_burn_time * dt;

Y(1) = 0;
V(1) = 0;
M(1) = Empty_rocket_mass + Fuel_mass;

for i = 1:Length - 1
    if (i <= Fuel_burn_time / dt)
        M(i + 1) = M(i) - dm;
        V(i + 1) = V(i) - V_exhaust * dm / M(i);
    else
        M(i + 1) = M(i);
        V(i + 1) = V(i);
    end
    
    Y(i + 1) = Y(i) + V(i + 1);
end

figure(1);
plot(Time, Y);
title('Position vs Time');
xlabel('Time (s)');
ylabel('Height (m)');

figure(2);
plot(Time, V);
title('Velocity vs Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');