function res = system6(x)
global k_first;
global k_second;
global t;
global t1;
global gL;
global gS;
global geps;
res(1) = x_1(t, t1, k_second, x(1), psi_2(t1, 0, k_first, x(1), x(2)), x_1(t1, 0, k_first, x(1), x(2), 0, 0), x_2(t1, 0, k_first, x(1), x(2), 0)) - gL - geps;
res(2) = x_2(t, t1, k_second, x(1), psi_2(t1, 0, k_first, x(1), x(2)), x_2(t1, 0, k_first, x(1), x(2), 0)) - gS - geps;
end

