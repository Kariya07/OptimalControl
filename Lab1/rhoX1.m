function [val, point] = rhoX1(l, x1, eps)
    val = eps*norm(l)+ l' * x1;
    point = (eps/norm(l)) * l + x1;
    %val = l' * x1;
    %point = x1;
end

