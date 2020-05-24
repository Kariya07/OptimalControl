function res = system1(x)
global k_first;
global k_second;
global t1;
global t2;
res(1) = psi_2(t1, 0, k_first, x(1), x(2)) - 1/2;
res(2) = x(1) + psi_2(t2, t1, k_second, x(1), psi_2(t1, 0, k_first, x(1), x(2)))*...
    x_2(t2, t1, k_second, x(1), psi_2(t1, 0, k_first, x(1), x(2)), x_2(t1, 0, k_first, x(1), x(2), 0));
end

