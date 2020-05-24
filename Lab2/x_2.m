function res = x_2(t, t0, k, psi1, psi20, x0)
res = (1/(2*k) + psi1/k^2 - exp(k*(t-t0))/(2*k)) + psi20*sinh(k*(t-t0))/k -... 
psi1*cosh(k*(t-t0))/k^2 + x0*exp(k*(t-t0));
end

