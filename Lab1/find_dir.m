function res = find_dir(x, rho, n)
    phi = linspace(0, 2 * pi, n + 1);
    phi(n + 1) = [];
    l = [cos(phi); sin(phi)];
    dif = zeros(1, n);
    for i = 1:n
        [val, ~] = rho(l(:, i));
        dif(i) = abs(val - l(:, i)' * x);
    end
    [~, idx] = min(dif);
    res = l(:, idx);
end

