function res = intersection(A1, A2, B1, B2)
    x = [A1(1) A2(1) B1(1) B2(1)];
    y = [A1(2) A2(2) B1(2) B2(2)];
    k = convhull(x, y);
    k(end) = [];
    if (numel(k) == 4 && k(3) == 2)
        disp(k)
        res = true;
    else
        res = false;
    end
end
