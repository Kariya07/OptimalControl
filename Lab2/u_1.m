function res = u_1(t, t0, k, psi1, psi20, flag)
psi2 = -1/2 + psi_2(t, t0, k, psi1, psi20);
if (flag == 1)
    res = (psi2 >= 0);
    res = psi2.*res;
else
    res = psi2;
end
end

