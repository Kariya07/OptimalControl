function drawSet(rho, N, color)
    i = 1 : N;
    L = [cos(2 * pi * i / N); sin(2 * pi * i / N)];
    vals = zeros(N, 1);
    points = zeros(2, N + 1);
    for i = 1 : N
        [vals(i), points(:, i)] = rho(L(:, i));
    end
    points(:, N + 1) = points(:, 1);
    fill(points(1, :), points(2, :), color);
end

