function res = u0_x_1(t, t0, k, x10, x20)
res = k*(t - t0) + exp(k * (t - t0))*x20/k + (x10 - x20/k);
end

