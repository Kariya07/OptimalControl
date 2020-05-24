function [min_J, opt_x1, opt_x2, opt_u1, opt_u2, opt_psi1, opt_psi2, opt_t1, opt_t2] = task2(k1, k2, T, L, S, eps, t_h)
global k_first;
global k_second;
global t1;
global t2;
iter = 0;
t_n = T/t_h + 1;
opt_x1 = zeros([1, t_n]);
opt_x2 = zeros([1, t_n]);
opt_u1 = zeros([1, t_n]);
opt_u2 = zeros([1, t_n]);
opt_psi1 = zeros([1, t_n]);
opt_psi2 = zeros([1, t_n]);
min_J = 0;
opt_t1 = 0;
opt_t2 = 0;

%psi1 = 0
disp('Особый случай')
for c = [S - eps, S + eps]
    for k = [k1, k2]
        psi20 = (-1/2 + exp(k*T)/2 + c*k)/sinh(k*T);
        if abs(x_1(T, 0, k, 0, psi20, 0, 0) - L) <= eps
            u1 = u_1(linspace(0, T, t_n), 0, k, 0, psi20, 2);
            J = trapz(linspace(0, T, t_n), u1 + u1.^2);
            disp([J, min_J]);
            if (iter == 0) || ((iter > 0) && (J < min_J))
                iter = iter + 1;
                min_J = J;
                opt_x1 = x_1(linspace(0, T, t_n), 0, k, 0, psi20, 0, 0);
                opt_x2 = x_2(linspace(0, T, t_n), 0, k, 0, psi20, 0);
                opt_u1 = u1;
                opt_u2 = k*ones(1, t_n);
                opt_psi1 = zeros(1, t_n);
                opt_psi2 = psi_2(linspace(0, T, t_n), 0, k, 0, psi20);
            end
        end
    end
end

%%Одно переключение 
disp('Без переключений или одно переключение');
%одно переключение по u2
global t;
global gL;
global gS;
global geps;
gL = L;
gS = S;
geps = eps;
t = T;
for m = 1:2
    if m == 1
        k_first = k1;
        k_second = k2;
    else
        k_first = k2;
        k_second = k1;
    end
    for t1 = linspace(0, T, t_n)
        disp("----------------------- NEW ITERATION11 ------------------------");
        for j = 1:4
            if j == 1
                psi = fsolve(@system3,[0, 0], optimset('Display', 'off'));
            end
            if j == 2
                psi = fsolve(@system4,[0, 0], optimset('Display', 'off'));
            end
            if j == 3
                psi = fsolve(@system5,[0, 0], optimset('Display', 'off'));
            end
            if j == 4
                psi = fsolve(@system6,[0, 0], optimset('Display', 'off'));
            end
            psi1 = psi(1, 1);
            psi20 = psi(1, 2);
            x1_1 = x_1(t1, 0, k_first, psi1, psi20, 0, 0);
            x2_1 = x_2(t1, 0, k_first, psi1, psi20, 0);
            psi21 = psi_2(t1, 0, k_first, psi1, psi20);
            x1_2 = x_1(T, t1, k_second, psi1, psi21, x1_1, x2_1);
            x2_2 = x_2(T, t1, k_second, psi1, psi21, x2_1);
            if (abs(x1_2 - L) <= eps && abs(x2_2 - S) <= eps)
                t_t1 = linspace(0, t1, round(t1/t_h + 1));
                t1_T = linspace(t1+t_h, T, round((T - t1)/t_h));
                u1(1 : round(t1/t_h+1)) = u_1(t_t1, 0, k_first, psi1, psi20, 2);
                u1(round(t1/t_h+2) : t_n) = u_1(t1_T, t1, k_second, psi1, psi21, 2);
                if t1 == 0
                    u1 = u_1(linspace(0, T, t_n), 0, k_second, psi1, psi20, 2);
                end
                if t1 == T
                    u1 = u_1(linspace(0, T, t_n), 0, k_first, psi1, psi20, 2);
                end
                J = trapz(linspace(0, T, t_n), u1 + u1.^2);
                disp([J, min_J]);
                if (iter == 0) || ((iter > 0) && (J < min_J))
                    iter = iter + 1;
                    min_J = J;
                    opt_u1 = u1;
                    if t1 == 0
                        opt_x1 = x_1(linspace(0, T, t_n), 0, k_second, psi1, psi20, 0, 0);
                        opt_x2 = x_2(linspace(0, T, t_n), 0, k_second, psi1, psi20, 0);
                        opt_u2 = k_second*ones(1, t_n);
                        opt_t1 = t1;
                        opt_psi1 = psi1*ones(1, t_n);
                        opt_psi2 = psi_2(linspace(0, T, t_n), 0, k_second, psi1, psi20);
                        continue;
                    end
                    if t1 == T
                        opt_x1 = x_1(linspace(0, T, t_n), 0, k_first, psi1, psi20, 0, 0);
                        opt_x2 = x_2(linspace(0, T, t_n), 0, k_first, psi1, psi20, 0);
                        opt_u2 = k_first*ones(1, t_n);
                        opt_t1 = t1;
                        opt_psi1 = psi1*ones(1, t_n);
                        opt_psi2 = psi_2(linspace(0, T, t_n), 0, k_first, psi1, psi20);
                        continue;
                    end
                    opt_x1(1 : round(t1/t_h+1)) = x_1(t_t1, 0, k_first, psi1, psi20, 0, 0);
                    opt_x1(round(t1/t_h+2) : t_n) = x_1(t1_T, t1, k_second, psi1, psi21, x1_1, x2_1);
                    opt_x2(1 : round(t1/t_h+1)) = x_2(t_t1, 0, k_first, psi1, psi20, 0);
                    opt_x2(round(t1/t_h+2) : t_n) = x_2(t1_T, t1, k_second, psi1, psi21, x2_1);
                    opt_u2(1 : round(t1/t_h+1)) = k_first;
                    opt_u2(round(t1/t_h+2) : t_n) = k_second;
                    opt_t1 = t1;
                    opt_psi1 = psi1*ones(1, t_n);
                    opt_psi2(1 : round(t1/t_h+1)) = psi_2(t_t1, 0, k_first, psi1, psi20);
                    opt_psi2(round(t1/t_h+2) : t_n) = psi_2(t1_T, t1, k_second, psi1, psi21);
                end
            end
        end
    end
end


%%Два переключения
disp('2 переключения');
for t1 = linspace(t_h, T, t_n-1)
    disp('----------------------- NEW ITERATION2 ----------------------');
    for t2 = linspace(t1+t_h, T, round(t_n - t1/t_h - 1))
        k_first = k1;
        k_second = k2;
        [psi2_1, ex21] = solve_system(@system2, 0, 4, false);
        if (ex21 < 4)
            [psi2_2, ex22] = solve_system(@system2, 0, 4-ex21, true);
        else
            ex22 = 0;
        end
        k_first = k2;
        k_second = k1;
        [psi4_1, ex41] = solve_system(@system2, 0, 4, false);
        if (ex41 < 4)
            [psi4_2, ex42] = solve_system(@system2, 0, 4-ex41, true);
        else
            ex42 = 0;
        end
        for i = 1:(ex21+ex22)
            if i <= ex21
                psi1 = psi2_1(i, 1);
                psi20 = psi2_1(i, 2);
            else
                psi1 = psi2_2(i-ex21, 1);
                psi20 = psi2_2(i-ex21, 2);
            end
            
            x1_1 = x_1(t1, 0, k1, psi1, psi20, 0, 0);
            x2_1 = x_2(t1, 0, k1, psi1, psi20, 0);
            psi21 = psi_2(t1, 0, k1, psi1, psi20);
            x1_2 = x_1(t2, t1, k2, psi1, psi21, x1_1, x2_1);
            x2_2 = x_2(t2, t1, k2, psi1, psi21, x2_1);
            psi22 = psi_2(t2, t1, k2, psi1, psi21);
            x1_3 = x_1(T, t2, k1, psi1, psi22, x1_2, x2_2);
            x2_3 = x_2(T, t2, k1, psi1, psi22, x2_2);
            if (abs(x1_3 - L) <= eps && abs(x2_3 - S) <= eps)
                t_t1 = linspace(0, t1, round(t1/t_h + 1));
                t1_t2 = linspace(t1+t_h, t2, round((t2 - t1)/t_h));
                t2_T = linspace(t2+t_h, T, round((T - t2)/t_h));
                u1(1 : round(t1/t_h+1)) = u_1(t_t1, 0, k1, psi1, psi20, 2);
                u1(round(t1/t_h+2) : round(t2/t_h+1)) = u_1(t1_t2, t1, k2, psi1, psi21, 2);
                u1(round(t2/t_h+2) : t_n) = u_1(t2_T, t2, k1, psi1, psi22, 2);
                J = trapz(linspace(0, T, t_n), u1 + u1.^2);
                disp([J, min_J]);
                if (iter == 0) || ((iter > 0) && (J < min_J))
                    iter = iter + 1;
                    min_J = J;
                    opt_x1(1 : round(t1/t_h+1)) = x_1(t_t1, 0, k1, psi1, psi20, 0, 0);
                    opt_x1(round(t1/t_h+2) : round(t2/t_h+1)) = x_1(t1_t2, t1, k2, psi1, psi21, x1_1, x2_1);
                    opt_x1(round(t2/t_h+2) : t_n) = x_1(t2_T, t2, k1, psi1, psi22, x1_2, x2_2);
                    opt_x2(1 : round(t1/t_h+1)) = x_2(t_t1, 0, k1, psi1, psi20, 0);
                    opt_x2(round(t1/t_h+2) : round(t2/t_h+1)) = x_2(t1_t2, t1, k2, psi1, psi21, x2_1);
                    opt_x2(round(t2/t_h+2) : t_n) = x_2(t2_T, t2, k1, psi1, psi22, x2_2);
                    opt_u1 = u1;
                    opt_u2(1 : round(t1/t_h+1)) = k1;
                    opt_u2(round(t1/t_h+2) : round(t2/t_h+1)) = k2;
                    opt_u2(round(t2/t_h+2) : t_n) = k1;
                    opt_t1 = t1;
                    opt_t2 = t2;
                    opt_psi1 = psi1*ones(1, t_n);
                    opt_psi2(1 : round(t1/t_h+1)) = psi_2(t_t1, 0, k1, psi1, psi20);
                    opt_psi2(round(t1/t_h+2) : round(t2/t_h+1)) = psi_2(t1_t2, t1, k2, psi1, psi21);
                    opt_psi2(round(t2/t_h+2) : t_n) = psi_2(t2_T, t2, k1, psi1, psi22);
                end
            end
        end
        
        for i = 1:(ex41+ex42)
            if i <= ex41
                psi1 = psi4_1(i, 1);
                psi20 = psi4_1(i, 2);
            else
                psi1 = psi4_2(i-ex41, 1);
                psi20 = psi4_2(i-ex41, 2); 
            end
            
            x1_1 = x_1(t1, 0, k2, psi1, psi20, 0, 0);
            x2_1 = x_2(t1, 0, k2, psi1, psi20, 0);
            psi21 = psi_2(t1, 0, k2, psi1, psi20);
            x1_2 = x_1(t2, t1, k1, psi1, psi21, x1_1, x2_1);
            x2_2 = x_2(t2, t1, k1, psi1, psi21, x2_1);
            psi22 = psi_2(t2, t1, k1, psi1, psi21);
            x1_3 = x_1(T, t2, k2, psi1, psi22, x1_2, x2_2);
            x2_3 = x_2(T, t2, k2, psi1, psi22, x2_2);
            if (abs(x1_3 - L) <= eps && abs(x2_3 - S) <= eps)
                t_t1 = linspace(0, t1, round(t1/t_h + 1));
                t1_t2 = linspace(t1+t_h, t2, round((t2 - t1)/t_h));
                t2_T = linspace(t2+t_h, T, round((T - t2)/t_h));
                u1(1 : round(t1/t_h+1)) = u_1(t_t1, 0, k2, psi1, psi20, 2);
                u1(round(t1/t_h+2) : round(t2/t_h+1)) = u_1(t1_t2, t1, k1, psi1, psi21, 2);
                u1(round(t2/t_h+2) : t_n) = u_1(t2_T, t2, k2, psi1, psi22, 2);
                J = trapz(linspace(0, T, t_n), u1 + u1.^2);
                disp([J, min_J]);
                if (iter == 0) || ((iter > 0) && (J < min_J))
                    iter = iter + 1;
                    min_J = J;
                    opt_x1(1 : round(t1/t_h+1)) = x_1(t_t1, 0, k2, psi1, psi20, 0, 0);
                    opt_x1(round(t1/t_h+2) : round(t2/t_h+1)) = x_1(t1_t2, t1, k1, psi1, psi21, x1_1, x2_1);
                    opt_x1(round(t2/t_h+2) : t_n) = x_1(t2_T, t2, k2, psi1, psi22, x1_2, x2_2);
                    opt_x2(1 : round(t1/t_h+1)) = x_2(t_t1, 0, k2, psi1, psi20, 0);
                    opt_x2(round(t1/t_h+2) : round(t2/t_h+1)) = x_2(t1_t2, t1, k1, psi1, psi21, x2_1);
                    opt_x2(round(t2/t_h+2) : t_n) = x_2(t2_T, t2, k2, psi1, psi22, x2_2);
                    opt_u1 = u1;
                    opt_u2(1 : round(t1/t_h+1)) = k2;
                    opt_u2(round(t1/t_h+2) : round(t2/t_h+1)) = k1;
                    opt_u2(round(t2/t_h+2) : t_n) = k2;
                    opt_t1 = t1;
                    opt_t2 = t2;
                    opt_psi1 = psi1*ones(1, t_n);
                    opt_psi2(1 : round(t1/t_h+1)) = psi_2(t_t1, 0, k2, psi1, psi20);
                    opt_psi2(round(t1/t_h+2) : round(t2/t_h+1)) = psi_2(t1_t2, t1, k1, psi1, psi21);
                    opt_psi2(round(t2/t_h+2) : t_n) = psi_2(t2_T, t2, k2, psi1, psi22);
                end
            end
        end
    end 
end
if iter == 0
    disp("No solution");
else
    disp("Found solution");
end
end

