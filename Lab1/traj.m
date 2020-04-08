function [X, U, Psi, opt_traj_ind, t1, t_res] = traj(psi_0, psi_0_n, A, B, f, h_t, T_max, p1, p2, alpha, beta, gamma, delta, k, x0, x1, eps)
    t = 0:h_t:1;
    T = numel(t);
    U = zeros(2, T, psi_0_n);
    X = U;
    Psi = U;
    
    t1_opt = 0;
    opt_traj_ind = 0;
    is_reached = 0;
    conj_syst = @(t, y) -A' * y;
    first_time = 1;
    while(~is_reached && t(T) <= T_max) 
        for i = 1:psi_0_n
            if first_time
                psi_0_init = psi_0(:, i);
            else
                psi_0_init = Psi(:, t1_opt, i);
            end
            [~, psi] = ode45(conj_syst, t, psi_0_init);
            psi = psi';
            u = zeros(2, T);
            for j = 1:T
                 [~, u(:, j)] = rhoP(B' * psi(:, j), p1, p2, alpha, beta, gamma, delta);
            end
            if first_time
                x_init = x1;
            else
                x_init = X(:, t1_opt, i);
            end
            u_func = @(t) u(:, min(1 + round(t / h_t), T));
            [~, x] = ode45(@(t, y) A * y + B * u_func(t) + f,  0:h_t:1, x_init);
            x = x';
            for j = 1:T
                if in_square(x(:, j), x0, k, 0)
                    is_reached = 1;
                    T = j;
                    opt_traj_ind = i;
                    X = X(:, 1 : t1_opt + T, :);
                    U = U(:, 1 : t1_opt + T, :);
                    Psi = Psi(:, 1 : t1_opt + T, :);
                    break
                end
            end

            x = x(:, 1:T);
            u = u(:, 1:T);
            psi = psi(:, 1:T);
            Psi(:, t1_opt + (1:T), i) = psi;
            X(:, t1_opt + (1:T), i) = x;
            U(:, t1_opt + (1:T), i) = u;

        end
        t1_opt = t1_opt + T;
        if ~is_reached 
            if (1 + t(T)) <= T_max
                Psi(:, t1_opt+1:t1_opt+T, :) = zeros(2, T, psi_0_n);
                X(:, t1_opt+1:t1_opt+T, :) = zeros(2, T, psi_0_n);
                U(:, t1_opt+1:t1_opt+T, :) = zeros(2, T, psi_0_n);
            end
            first_time = 0;
            t = t + 1;
        end
    end
    disp("Ñompleted")
    if is_reached
        disp("Solved");
        disp("t1 = ");
        disp(t(T));
    else
        disp("No solution");
    end
    t1 = t(T);
    t_res = linspace(0, t1, size(X, 2));
end

