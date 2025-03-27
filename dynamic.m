function [theta1_ddot, theta2_ddot] = dynamic(theta1, theta1_dot, theta2, theta2_dot)

g = 9.81;
l1 = 1.0;
l2 = 1.0;
m1 = 1.0;
m2 = 1.0;

theta1_ddot = (-g * (2 * m1 + m2) * sin(theta1) - m2 * g * sin(theta1 - 2 * theta2) - ...
    2 * sin(theta1 - theta2) * m2 * (theta2_dot^2 * l2 + theta1_dot^2 * l1 * cos(theta1 - theta2))) / ...
    (l1 * (2 * m1 + m2 - m2 * cos(2 * theta1 - 2 * theta2)));

theta2_ddot = (2 * sin(theta1 - theta2) * (theta1_dot^2 * l1 * (m1 + m2) + g * (m1 + m2) * cos(theta1) + ...
    theta2_dot^2 * l2 * m2 * cos(theta1 - theta2))) / ...
    (l2 * (2 * m1 + m2 - m2 * cos(2 * theta1 - 2 * theta2)));

end