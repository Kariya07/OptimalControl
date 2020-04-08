function res = in_square(point, x0, k, flag)
    res = 0;
    point = point - x0;
    if flag
        if abs(point(1)) < k/2 && abs(point(2)) < k/2
            res = 1;
        end
    else
        if abs(point(1)) <= k/2 && abs(point(2)) <= k/2
            res = 1;
        end
    end
end
