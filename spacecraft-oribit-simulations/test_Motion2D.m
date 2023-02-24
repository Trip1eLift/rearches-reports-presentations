test = Motion2D;

theta = -pi:0.05:pi;
theta_after = zeros(1, length(theta));

for i = 1:length(theta)
    test = test.setRT(1, theta(i));
    theta_after(i) = test.getT();
end

plot(theta, theta_after);