function res = u_2(t, t0, k, k1, k2, psi1, psi20, x0)
if(xi(t, t0, k, psi1, psi20, x0) > 0)
    res = k2;
else
    if(xi(t, t0, k, psi1, psi20, x0) < 0)
        res = k1;
    else
        error ("xi == 0");
    end
end
end

