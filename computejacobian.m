function Jacobian_matrix = computejacobian(x_hat)

g = 9.81;
l1 = 1.0;
l2 = 1.0;   
mu_M1 = 0.5;
mu_M2 = 0.5;

syms x1 x2 x3 x4 real

x1_val = x_hat(1);
x2_val = x_hat(2);
x3_val = x_hat(3);
x4_val = x_hat(4);

beta = x3 - x1;

dx1 = x2;
dx2 = (g * (mu_M2 * cos(beta) * sin(x3) - sin(x1)) + ...
       mu_M2 * l1 * sin(beta) * cos(beta) * x2^2 + ...
       mu_M2 * l2 * sin(beta) * x4^2) / ...
       (l1 * (mu_M1 + mu_M2 * sin(beta)^2));

dx3 = x4;
dx4 = -(l1 * sin(beta) * x2^2 + ...
        l2 * mu_M2 * sin(beta) * cos(beta) * x4^2 + ...
        g * cos(x1) * sin(beta)) / ...
        (l2 * (mu_M1 + mu_M2 * sin(beta)^2));

state_vector = [x1; x2; x3; x4];
equations = [dx1; dx2; dx3; dx4];

Jacobian = jacobian(equations, state_vector);
Jacobian_func = matlabFunction(Jacobian, 'Vars', {x1, x2, x3, x4});

Jacobian_matrix = Jacobian_func(x1_val, x2_val, x3_val, x4_val);


end