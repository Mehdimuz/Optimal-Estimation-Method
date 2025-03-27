function R = adapt_R(innovation)
base_noise = [0.05, 0.05, 0.05, 0.05];
adaptive_noise = 0.01 * abs(innovation);
R = diag(base_noise + adaptive_noise);
end