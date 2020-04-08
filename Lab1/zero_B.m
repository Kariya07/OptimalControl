function res = zero_B(B, B_eps)
    [V, J] = jordan(B);
    disp(J);
    cond = abs(diag(J)) < B_eps;
    J = J + diag(cond*eps_B);
    res = V*J/V;
end

