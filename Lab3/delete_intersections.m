function [new_X, new_Y] = delete_intersections(X, Y)
    P = [X'; Y'];
    n = numel(X);
    ind = find(X - max(X) == 0, 1, 'first');
    if ind ~= 1
        P = [P(:, 1:ind-1), P(:, ind:n)];
    end
    diam = max(max(abs(X)), max(abs(Y)));
    k = 30;
    k2 = 5;
    for i = 2:n-1
        if i > size(P, 2)
            break;
        end
        for j = size(P, 2)-1:(-1):i+1
            if (norm(P(:, i) - P(:, j)) < (diam / k2)  && mod(j - i, n) > n / k && mod(i - j, n) > n /k)
                if intersection(P(:, i - 1), P(:, i + 1), P(:, j - 1), P(:, j + 1))
                    P(:,(i+1):j) = [];
                    break;
                end
            end
        end
    end
    new_X = P(1, :);
    new_Y = P(2, :);
end

