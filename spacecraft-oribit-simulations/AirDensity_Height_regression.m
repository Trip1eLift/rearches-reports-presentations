% Data points
clear;
% https://www.engineeringtoolbox.com/standard-atmosphere-d_604.html
H_standard = [-1000 0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 15000 20000 25000 30000 40000 50000 60000 70000 80000];
Rho_standard = [1.347 1.225 1.112 1.007 0.9093 0.8194 0.7364 0.6601 0.5900 0.5258 0.4671 0.4135 0.1948 0.08891 0.04008 0.01841 0.003996 0.001027 0.0003097 0.00008283 0.00001846];

%plot(H_standard, Rho_standard);

% Over 80 KM: Jacchia Reference Atmosphere
% https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19830012203.pdf
H_Jacchia = [90 100 110 125 140 180 420 500 700 1500 2500];
H_Jacchia = H_Jacchia * 1000;
Rho_Jacchia = [0.344*10^-5 0.524*10^-6 0.966*10^-7 0.134*10^-7 0.384*10^-8 0.559*10^-9 0.218*10^-11 0.574*10^-12 0.336*10^-13 0.577*10^-15 0.650*10^-16];

%plot(H_Jacchia, Rho_Jacchia);

H_mix = zeros(1, length(H_standard) + length(H_Jacchia));
Rho_mix = zeros(1, length(H_mix));
H_mix(1:length(H_standard)) = H_standard;
H_mix((length(H_standard)+1):length(H_mix)) = H_Jacchia;
Rho_mix(1:length(H_standard)) = Rho_standard;
Rho_mix((length(H_standard)+1):length(H_mix)) = Rho_Jacchia;

figure(1);
plot(H_mix, Rho_mix, 'o');

H_mix = transpose(H_mix);
H_mix_shift = H_mix + 2000;
Rho_mix = transpose(Rho_mix);
logH_shift = log(H_mix_shift);
logRho = log(Rho_mix);

% linear log fit
%{
f = fit(logH_shift, logRho, 'poly1');
hold on
plot(H_mix, exp(f.p2) .* (H_mix + 2000) .^ f.p1);
xlim([-1000, 2500000]);
ylim([0, 1.5]);
hold off

legend('data' ,'fit');

figure(2);
plot(log(H_mix + 2000), log(Rho_mix), 'o');
hold on
plot(log(H_mix + 2000), log(exp(f.p2) .* (H_mix + 2000) .^ f.p1));
hold off
%}

% quadratic log fit

f2 = fit(logH_shift, logRho, 'poly2')
hold on
plot(H_mix, exp(f2.p3) .* (H_mix + 2000) .^ f2.p2 .* (H_mix + 2000) .^ (f2.p1 * log(H_mix + 2000)) );
hold off

figure(2);
plot(log(H_mix + 2000), log(Rho_mix), 'o');
hold on
plot(log(H_mix + 2000), log(exp(f2.p3) .* (H_mix + 2000) .^ f2.p2 .* (H_mix + 2000) .^ (f2.p1 * log(H_mix + 2000))));
hold off


%{
f3 = fit(logH_shift, logRho, 'poly3')
hold on
plot(H_mix, exp(f3.p4) .* ((H_mix + 2000) .^ f3.p3) .* ((H_mix + 2000) .^ (f3.p2 * log(H_mix + 2000))) .* ((H_mix + 2000) .^ (f3.p1 * log(H_mix + 2000) .^ 2)));
hold off

figure(2);
plot(log(H_mix + 2000), log(Rho_mix), 'o');
hold on
plot(log(H_mix + 2000), log(exp(f3.p4) .* ((H_mix + 2000) .^ f3.p3) .* ((H_mix + 2000) .^ (f3.p2 * log(H_mix + 2000))) .* ((H_mix + 2000) .^ (f3.p1 * log(H_mix + 2000) .^ 2))) );
%}