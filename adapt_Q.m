function Q = adapt_Q(innovation)
base_noise = [0.01, 0.01, 0.01, 0.01];
adaptive_noise = 0.1 * abs(innovation);
Q = diag(base_noise + adaptive_noise);
end
