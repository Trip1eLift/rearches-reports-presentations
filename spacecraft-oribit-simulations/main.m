% All units are SI

%% Section 1: simulate ISS station with real world data
% assume no drag
clear;
obj = Motion2D;
% set Earth postion at (0, 0) stablly
% height of iss = 408000 m
% radius of earth = 6371000 m
% earth mass = 5.972*10^24 kg
% gravtational constant = 6.67408*10^-11 m^3 kg^-1 s^-2

%{
Time_max = 5565; % 5565 second is about 90 minutes which is a period of ISS
dt = 1;          % time resolution is 1 sec
Length = Time_max/dt;

Time = 0:dt:Time_max-dt;
position(1, Length) = Motion2D;
velocity(1, Length) = Motion2D;
acceleration(1 , Length-1) = Motion2D;


position(1, 1) = obj.setRT(6779000, 0);
v_0 = 7667*dt;
velocity(1, 1) = obj.setXY(0, v_0);

for i = 1:Length - 1
    acceleration(1, i) = obj.setRT( -gravity_earth(position(1, i).getR()), position(1,i).getT() );
    velocity(1, i+1) = velocity(1, i) + acceleration(1, i);
    position(1, i+1) = position(1, i) + velocity(1, i+1);
end

figure(1);
Motion2D_Theta_vs_Time_plot(position, Time);
title('Position of ISS in Theta vs Time');
xlabel('Time (s)');
ylabel('Theta (rad)');

figure(2);
Motion2D_Radius_vs_Time_plot(position, Time);
title('Position of ISS in Radius vs Time');
xlabel('Time (s)');
ylabel('Radius (m)');

figure(3);
Motion2D_XY_plot(position);
title('Position of ISS X vs Y');
xlabel('X (m)');
ylabel('Y (m)');

figure(4);
Motion2D_XY_plot(velocity);
title('Velocity X vs Y');
xlabel('V_x (m/s)');
ylabel('V_y (m/s)');

figure(5);
Motion2D_XY_plot(acceleration);
title('Acceleration X vs Y');
xlabel('a_x (m/s^2)');
ylabel('a_y (m/s^2)');

figure(6);
% ISS weight 419700 kg
Motion2D_Energy_of_system_vs_Time_plot(position, velocity, 419700, Time);
title('Energy of ISS system vs Time');
xlabel('Time (s)');
ylabel('Energy (J)');
%}

%% Section 2: Let ISS travels for more periods

%{
Time_max = 100000;
dt = 1;          % time resolution is 1 sec
Length = Time_max/dt;
Time = 0:dt:Time_max-dt;

position(1, Length) = Motion2D;
velocity(1, Length) = Motion2D;
acceleration(1 , Length-1) = Motion2D;

position(1, 1) = obj.setRT(6779000, 0);
v_0 = 7667*dt;
velocity(1, 1) = obj.setXY(0, v_0);

for i = 1:Length - 1
    acceleration(1, i) = obj.setRT( -gravity_earth(position(1, i).getR()), position(1,i).getT() );
    velocity(1, i+1) = velocity(1, i) + acceleration(1, i);
    position(1, i+1) = position(1, i) + velocity(1, i+1);
end

earth(1, 300) = Motion2D;
for i = 1:300
    earth(i) = obj.setRT(6371000, 2*pi/299 * i);
end

figure(7);
Motion2D_XY_plot(position);
hold on
Motion2D_XY_plot(earth);
hold off
title('Position of ISS');
legend('ISS orbit', 'Earth');
xlabel('X (m)');
ylabel('Y (m)');

figure(8);
Motion2D_Energy_of_system_vs_Time_plot(position, velocity, 419700, Time);
title('Energy of ISS system vs Time');
xlabel('Time (s)');
ylabel('Energy (J)');
%}


%% Section 3: Apply drag on ISS while let it travels for 4 periods

%{
Time_max = 22260;
dt = 1;          % time resolution is 1 sec
Length = Time_max/dt;
Time = 0:dt:Time_max-dt;

position(1, Length) = Motion2D;
velocity(1, Length) = Motion2D;
acceleration(1 , Length-1) = Motion2D;

position(1, 1) = obj.setRT(6779000, 0);
v_0 = 7667*dt;
velocity(1, 1) = obj.setXY(0, v_0);

for i = 1:Length - 1
    % gravitational and drag acceleration:
    grav_acc = obj.setRT( -gravity_earth(position(1, i).getR()), position(1,i).getT() );
    % Known: ISS weight 419700kg, the solar panel has area of 2500m^2
    % Assume ISS is a long cylinder (has a drag coefficient of 0.82)
    drag_acc = drag_on_height_and_velocity(position(1, i).getR() - 6371000, velocity(1, i), 419700, 2500, 0.82);
    acceleration(1, i) = grav_acc + drag_acc;
    velocity(1, i+1) = velocity(1, i) + acceleration(1, i);
    position(1, i+1) = position(1, i) + velocity(1, i+1);
end

earth(1, 300) = Motion2D;
for i = 1:300
    earth(i) = obj.setRT(6371000, 2*pi/299 * i);
end

figure(9);
Motion2D_XY_plot(position);
hold on
Motion2D_XY_plot(earth);
hold off
title('Position of ISS');
legend('ISS', 'Earth');
xlabel('X (m)');
ylabel('Y (m)');

figure(91);
Motion2D_XY_plot(position(1:5565));
hold on
Motion2D_XY_plot(position(5566:11130));
Motion2D_XY_plot(position(11131:16695));
Motion2D_XY_plot(position(16696:22260));
Motion2D_XY_plot(earth);
hold off
title('Position of ISS');
legend('Period 1', 'Period 2', 'Period 3', 'Period 4', 'Earth');
xlabel('X (m)');
ylabel('Y (m)');

figure(10);
Motion2D_Energy_of_system_vs_Time_plot(position, velocity, 419700, Time);
title('Energy of ISS system with drag vs Time');
xlabel('Time (s)');
ylabel('Energy (J)');

figure(11);
Motion2D_Velocity_vs_Time_plot(velocity, Time);
title('Velocity of ISS with drag vs Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');

figure(12);
Motion2D_Radius_vs_Time_plot(position, Time);
title('Radius of ISS orbit with drag vs Time');
xlabel('Time (s)');
ylabel('Radius (m)');
%}


%% Section 4: Rocket propulsion

%{
% data from the legendary Saturn V rocket
% http://www.mnealon.eosc.edu/RocketSciencePage5.htm

% stage 1 data: 
% frame weight: 0.136 * 10 ^ 6 kg
% fuel weight: 2.04 * 10 ^ 6 kg
% fuel burn duration: 150 s
% exhaust velocity: 2456 m/s
% Assumptions:
% 1. The rocket is launched vertically
% 2. The propulsion engine of stage one can only push the stage upward
% 3. The rocket is lauched from equator and carries the orbit velocity of
%    the equator

Time_max_stage_1 = 1700;
dt = 1;
Length_stage_1 = Time_max_stage_1 / dt;
Time_stage_1 = 0:dt:Time_max_stage_1 - dt;

stage_one = SpaceCraft;
stage_one = stage_one.properties(Time_max_stage_1, dt, 0.136 * 10 ^ 6, 2.04 * 10 ^ 6, 150, -2456);
stage_one = stage_one.initial_positionRT(6371000, 0);
dTheta_0 = dt/(24 * 3600);
stage_one = stage_one.initial_velocityXY(0, 6371000 * dTheta_0 / cos(dTheta_0));

% Stage 2 data:
% Engine ignore: 180s
% Burn duration: 360s
% Frame mass: 0.0432 * 10^6 kg
% Fuel mass: 0.428 * 10^6 kg
% Exhaust velocity: 4220 m/s
% Assumptions:
% 1. Stage two rocket pushes itself forward
% 2. Stage two rocket carries the position and velocity of stage one when
%    sparates.
Time_start_stage_2 = 180;
Time_max_stage_2 = 2000;
Time_stage_2 = Time_start_stage_2:dt:Time_max_stage_2 - dt;
Length_stage_2 = length(Time_stage_2);
stage_two = SpaceCraft;
stage_two = stage_two.properties(Time_max_stage_2 - Time_start_stage_2, dt, 0.0432 * 10^6, 0.428 * 10^6, 360, -4220);

% Stage 1 & Stage 2 separation
stage_one.mass(Time_start_stage_2 / dt : Length_stage_1) = stage_one.mass(Time_start_stage_2 / dt : Length_stage_1) - 0.4712 * 10^6;

% stage one simulation
acc_temp = Motion2D;
for i = 1:Length_stage_1 - 1
    % gravitational, drag, and propulsion acceleration:
    grav_acc = obj.setRT( -gravity_earth(stage_one.position(i).getR()), stage_one.position(i).getT() );
    % Known: Diameter of Saturn V = 10.1 m
    %        => Cross section area = pi * (10.1/2)^2 = 80.12
    % Online bloger says the drag of Saturn V at max Q is about 0.51
    % https://forum.nasaspaceflight.com/index.php?topic=13543.880
    drag_acc = drag_on_height_and_velocity(stage_one.position(i).getR() - 6371000, stage_one.velocity(i), stage_one.mass(i), 80.12, 0.51);
    
    % The rocket only pushed itself when there is fuel left
    if (i <= stage_one.fuel_burn_duration / dt)
        prop_acc = stage_one.upward_aceleration_at(i);
        acc_temp = grav_acc  + drag_acc + prop_acc;
    else
        acc_temp = grav_acc  + drag_acc;
    end
    
    stage_one.velocity(i+1) = stage_one.velocity(i) + acc_temp;
    stage_one.position(i+1) = stage_one.position(i) + stage_one.velocity(i+1);
end

stage_two.position(1) = stage_one.position(Time_start_stage_2 / dt);
stage_two.velocity(1) = stage_one.velocity(Time_start_stage_2 / dt);

% stage one simulation
for i = 1:Length_stage_2 - 1
    % gravitational, drag, and propulsion acceleration:
    grav_acc = obj.setRT( -gravity_earth(stage_two.position(i).getR()), stage_two.position(i).getT() );
    % Known that stage 2 behaive simular to stage 1 on drag
    drag_acc = drag_on_height_and_velocity(stage_two.position(i).getR() - 6371000, stage_two.velocity(i), stage_two.mass(i), 80.12, 0.51);
    if (i <= stage_two.fuel_burn_duration / dt)
        prop_acc = stage_two.upward_aceleration_at(i);
        acc_temp = grav_acc  + drag_acc + prop_acc;
    else
        acc_temp = grav_acc  + drag_acc;
    end
    
    stage_two.velocity(i+1) = stage_two.velocity(i) + acc_temp;
    stage_two.position(i+1) = stage_two.position(i) + stage_two.velocity(i+1);
end

figure(13);
Motion2D_XY_plot(stage_one.position);
earth(1, 300) = Motion2D;
for i = 1:300
    earth(i) = obj.setRT(6371000, 2*pi/299 * i);
end
hold on
Motion2D_XY_plot(stage_two.position);
Motion2D_XY_plot(earth);
hold off
title('Position of Saturn V');
legend('Stage 1', 'Stage 2', 'Earth');
xlabel('X (m)');
ylabel('Y (m)');

figure(14);
Motion2D_Earth_Height_vs_Time_plot(stage_one.position, Time_stage_1);
title('Height of Saturn V stage 1');
xlabel('Time (s)');
ylabel('Height (m)');


figure(15);
Motion2D_Earth_Height_vs_Time_plot(stage_one.position(1:200), Time_stage_1(1:200));
title('Height of Saturn V stage 1 (first 200 sec)');
xlabel('Time (s)');
ylabel('Height (m)');

figure(16);
Motion2D_Earth_Height_vs_Time_plot(stage_two.position, Time_stage_2);
title('Height of Saturn V stage 2');
xlabel('Time (s)');
ylabel('Height (m)');

figure(17);
Motion2D_Earth_Height_vs_Time_plot(stage_one.position, Time_stage_1);
hold on
Motion2D_Earth_Height_vs_Time_plot(stage_two.position, Time_stage_2);
hold off
title('Height of Saturn V');
xlabel('Time (s)');
ylabel('Height (m)');
legend('Stage 1', 'Stage 2');

figure(18);
Motion2D_Velocity_vs_Time_plot(stage_one.velocity, Time_stage_1);
hold on
Motion2D_Velocity_vs_Time_plot(stage_two.velocity, Time_stage_2);
hold off
title('Speed of Saturn V');
xlabel('Time (s)');
ylabel('Speed (m/s)');
legend('Stage 1', 'Stage 2');

figure(19);
Motion2D_Theta_vs_Time_plot(stage_one.position, Time_stage_1);
earth_angle = Time_stage_1 / (24 * 3600);
hold on
plot(Time_stage_1, earth_angle);
hold off
title('Angle vs Time (Predict the drop point of Saturn V stage 1)');
legend('stage 1 angle', 'Earth angle');
xlabel('Time (s)');
ylabel('Angle (rad)');

figure(20);
SpaceCraft_Energy_vs_Time_plot(stage_one, Time_stage_1);
hold on
SpaceCraft_Energy_vs_Time_plot(stage_two, Time_stage_2);
hold off
title('Energy of Stage one vs Time');
xlabel('Time (s)');
ylabel('Energy (J)');
legend('Stage 1', 'Stage 2');
%}

%% Helper functions for plot
function Motion2D_XY_plot(motion2D_arr)
    X = zeros(1, length(motion2D_arr));
    Y = zeros(1, length(motion2D_arr));
    for i = 1:length(X)
        X(i) = motion2D_arr(i).getX();
        Y(i) = motion2D_arr(i).getY();
    end
    plot(X, Y);
end

function Motion2D_Theta_vs_Time_plot(motion2D_arr, time_arr)
    theta_arr = zeros(1, length(motion2D_arr));
    for i = 1:length(theta_arr)
        theta_arr(i) = motion2D_arr(i).getT();
    end
    plot(time_arr, theta_arr);
end

function Motion2D_Radius_vs_Time_plot(motion2D_arr, time_arr)
    radius_arr = zeros(1, length(motion2D_arr));
    for i = 1:length(radius_arr)
        radius_arr(i) = motion2D_arr(i).getR();
    end
    plot(time_arr, radius_arr);
end

function Motion2D_Earth_Height_vs_Time_plot(motion2D_arr, time_arr)
    radius_arr = zeros(1, length(motion2D_arr));
    for i = 1:length(radius_arr)
        radius_arr(i) = motion2D_arr(i).getR() - 6371000;
    end
    plot(time_arr, radius_arr);
end

function Motion2D_Velocity_vs_Time_plot(velocity, time)
    v = zeros(1, length(velocity));
    for i = 1:length(v)
        v(i) = sqrt(velocity(i).getX()^2 + velocity(i).getY()^2);
    end
    plot(time, v);
end

function Motion2D_Energy_of_system_vs_Time_plot(position, velocity, mass, time)
    % E = T + U
    % T = 1/2*m*(vx^2 + vy^2)
    % U = -GMm/R
    Len = length(position);
    E = zeros(1, Len);
    for i = 1:Len
        T = 1/2 * mass * (velocity(i).getX()^2 + velocity(i).getY()^2);
        U = -Gravity_constant * earth_mass * mass / position(i).getR();
        E(i) = T + U;
    end
    plot(time, E);
end

function SpaceCraft_Energy_vs_Time_plot(craft, time)
    Len = length(time);
    E = zeros(1, Len);
    for i = 1:Len
        T = 1/2 * craft.mass(i) * (craft.velocity.getX()^2 + craft.velocity.getY()^2);
        U = -Gravity_constant * earth_mass * craft.mass(i) / craft.position(i).getR();
        E(i) = T + U;
    end
    plot(time, E);
end

%% Helper functions for calculations

function acc = drag_on_height_and_velocity(height, velocity, mass, area, Cd)
    acc = Motion2D;
    v_square = velocity.getX()^2 + velocity.getY()^2;
    
    FD = 0.5 * air_density_at_height(height) * v_square * Cd * area;
    
    acc_mag = FD / mass;
    
    v_mag = sqrt(v_square);
    acc = acc.setXY( -velocity.getX()/v_mag * acc_mag, -velocity.getY()/v_mag * acc_mag );
end

function Rho = air_density_at_height(height)
    % This function is only valid at height range from -1,000 to 2,500,000
    Rho = exp(-35.49) .* (height + 2000) .^ 10.63 .* (height + 2000) .^ (-0.7334 * log(height + 2000));
end

function acc = gravity_earth(r)
    acc = gravity_any(earth_mass, r);
end

function acc = gravity_any(M, r)
    acc = Gravity_constant * M / (r^2);
end

function m = earth_mass()
    m = 5.972*10^24;
end

function g = Gravity_constant
    g = 6.67408*10^-11;
end