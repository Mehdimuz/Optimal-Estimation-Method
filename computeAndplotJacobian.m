function Jacobian_matrix = computeAndplotJacobian(state)

g = 9.81;
l1 = 1.0;
l2 = 1.0;   
mu_M1 = 0.5;
mu_M2 = 0.5;

syms x1 x2 x3 x4 real

x1_val = state(1);
x2_val = state(2);
x3_val = state(3);
x4_val = state(4);

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

disp('Jacobian Matrix at the given state:');
disp(Jacobian_matrix)

figure;
[X1, X2] = meshgrid(linspace(-1, 1, 50), linspace(-1, 1, 50));
for i = 1:4
    for j = 1:4
       J_element = matlabFunction(Jacobian(i, j), 'Vars', {x1, x2, x3, x4});
       Z = arrayfun(@(x1, x2) J_element(x1, x2, x3_val, x4_val), X1, X2);
       subplot(4, 4, (i-1) * 4 + j);
       surf(X1, X2, Z)
       xlabel('x1'); ylabel('x2'); zlabel('Value');
       title(sprintf('J(%d , %d)', i, j));
       shading interp
    end
end

sgtitle('Jacobian at Initial State');

end