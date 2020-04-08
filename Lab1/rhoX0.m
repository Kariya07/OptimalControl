function [val, point] = rhoX0(l, k, x0, eps)
    val = k/2 * (abs(l(1)) + abs(l(2))) + eps*norm(l) + l' * x0;
    points = [k/2 k/2 -k/2 -k/2; k/2 -k/2 k/2 -k/2];
    [~, ind] = max(l' * points);
    point = x0 + eps*l/norm(l) + points(:, ind);
end

