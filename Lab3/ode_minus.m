function f = ode_minus(~, x, alpha)
    f = zeros(4, 1);
    f(1) = x(2);
    f(2) = -alpha - x(2) - 2*x(1) - x(1).*sin(x(1).^2) + 2*x(1).^2.*cos(x(1));
    f(3) = x(4) .* (2 + sin(x(1).^2) + 2*x(1).^2.*cos(x(1).^2) - 4*x(1).*cos(x(1)) + 2*x(1).^2.*sin(x(1)));
    f(4) = x(4) - x(3);
end

