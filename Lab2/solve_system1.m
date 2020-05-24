function [res, ex] = solve_system1(system, begin, dim)
x = fsolve(system,[begin, 0], optimset('Display', 'off'));
res = zeros([dim, 2]);
res(1, :) = x;
num_of_answers = 1;
eps = 0.5;
iter = 0;
ex = 1;
h = 50;
while (num_of_answers < dim && iter < 4)
    a = x(1);
    b = a + h;
    n = 30;
    x = linspace(a, b, n);
    for x0 = x
        x1 = fsolve(system,[x0, 0], optimset('Display', 'off'));
        if(abs(x1(1) - x(1)) > eps)
            num_of_answers = num_of_answers + 1;
            x = x1;
            res(num_of_answers, :) = x1;
            if(num_of_answers == dim)
                ex = ex + 1;
                break;
            end
        end
    end
    iter = iter + 1;
    b = a - h;
    n = 30;
    x = linspace(b, a, n);
    for x0 = x
        x1 = fsolve(system,[x0, 0], optimset('Display', 'off'));
        if(abs(x1(1) - x(1)) > eps)
            num_of_answers = num_of_answers + 1;
            x = x1;
            res(num_of_answers, :) = x1;
            if(num_of_answers == dim)
                ex = ex + 1;
                break;
            end
        end
    end
    iter = iter + 1;
end
end

