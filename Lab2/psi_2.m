function res = psi_2(t, t0, k, psi1, psi20)
res = -psi1/k + (psi20 + psi1/k)*exp(-k*(t-t0));
end
