function [min_J, opt_x1, opt_x2, opt_u1, opt_u2, opt_psi1, opt_psi2, opt_t1, opt_t2] = task1(k1, k2, T, L, S, eps, t_h)
global k_first;
global k_second;
global t1;
global t2;
t_n = T/t_h + 1;
opt_x1 = zeros([1, t_n]);
opt_x2 = zeros([1, t_n]);
opt_u1 = zeros([1, t_n]);
opt_u2 = zeros([1, t_n]);
opt_psi1 = zeros([1, t_n]);
opt_psi2 = zeros([1, t_n]);
min_J = 0;
J_eps = 0.001;
opt_t1 = 0;
opt_t2 = 0;

%1
%без переключений
disp('Без переключений')
J = 0;
if ((abs(k2*T - L) <= eps) && (abs(S) <= eps))
    min_J = J;
    opt_x1 = k2 * linspace(0, T, t_n);
    opt_x2 = zeros([1, t_n]);
    opt_u1 = zeros([1, t_n]);
    opt_u2 = k2 * ones([1, t_n]);
    return
end

J = 0;
if ((abs(k1*T - L) <= eps) && (abs(S) <= eps))
    min_J = J;
    opt_x1 = k1 * linspace(0, T, t_n);
    opt_x2 = zeros([1, t_n]);
    opt_u1 = zeros([1, t_n]);
    opt_u2 = k1 * ones([1, t_n]);
    return
end

%psi1 = 0
disp('Особый случай')
J = 0;
for i = 1:50
    u2 = randi([k1, k2], 1, t_n);
    x1 = trapz(linspace(0, T, t_n), u2);
    disp(x1);
    if ((x1 <= k2*T) && (x1 >= k1*T)&&(abs(x1 - L) <= eps) && (abs(S) <= eps))
        min_J = J;
        time = linspace(0, T, t_n);
        for i = 1:t_n
            opt_x1(i) = trapz(time(1:i), u2(1:i));
        end
        opt_x2 = zeros([1, t_n]);
        opt_u1 = zeros([1, t_n]);
        opt_u2 = u2;
        return
    end
end

%2
disp('Одно переключение');
iter = 0;
%одно переключение по u2
global t;
global gL;
global gS;
global geps;
k_first = k1;
k_second = k2;
gL = L;
gS = S;
geps = eps;
t = T;
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
        if(psi1<0 && psi20>= 0.5)
            x1_1 = x_1(t1, 0, k1, psi1, psi20, 0, 0);
            x2_1 = x_2(t1, 0, k1, psi1, psi20, 0);
            psi21 = psi_2(t1, 0, k1, psi1, psi20);
            x1_2 = x_1(T, t1, k2, psi1, psi21, x1_1, x2_1);
            x2_2 = x_2(T, t1, k2, psi1, psi21, x2_1);
            disp([abs(x1_2 - L), abs(x2_2 - S)]);
            if (abs(x1_2 - L) <= eps && abs(x2_2 - S) <= eps)
                t_t1 = linspace(0, t1, round(t1/t_h + 1));
                t1_T = linspace(t1+t_h, T, round((T - t1)/t_h));
                u1(1 : round(t1/t_h+1)) = u_1(t_t1, 0, k1, psi1, psi20, 1);
                u1(round(t1/t_h+2) : t_n) = u_1(t1_T, t1, k2, psi1, psi21, 1);
                J = trapz(linspace(0, T, t_n), u1 + u1.^2);
                if (iter == 0) || ((iter > 0) && (J < min_J))
                    iter = iter + 1;
                    min_J = J;
                    opt_x1(1 : t1/t_h+1) = x_1(t_t1, 0, k1, psi1, psi20, 0, 0);
                    opt_x1(t1/t_h+2 : t_n) = x_1(t1_T, t1, k2, psi1, psi21, x1_1, x2_1);
                    opt_x2(1 : t1/t_h+1) = x_2(t_t1, 0, k1, psi1, psi20, 0);
                    opt_x2(t1/t_h+2 : t_n) = x_2(t1_T, t1, k2, psi1, psi21, x2_1);
                    opt_u1 = u1;
                    opt_u2(1 : round(t1/t_h+1)) = k1;
                    opt_u2(round(t1/t_h+2) : t_n) = k2;
                    opt_t1 = t1;
                    opt_psi1 = psi1*ones(1, t_n);
                    opt_psi2(1 : round(t1/t_h+1)) = psi_2(t_t1, 0, k1, psi1, psi20);
                    opt_psi2(round(t1/t_h+2) : t_n) = psi_2(t1_T, t1, k2, psi1, psi21);
                    if (J < J_eps)
                        return;
                    end
                end
            end
        end
    end
end

%3
% psi1 = 0
% одно переключение по u1
k_first = k2;
k_second = k2;
gL = L;
gS = S;
geps = eps;
t = T;
for t1 = linspace(0, T, t_n)
    disp("----------------------- NEW ITERATION1 ------------------------");
    psi20 = 0.5/(exp(-k2*t1));
    if psi20 > 0.5
        x1_1 = x_1(t1, 0, k2, 0, psi20, 0, 0);
        x2_1 = x_2(t1, 0, k2, 0, psi20, 0);
        x1_2 = u0_x_1(T, t1, k2, x1_1, x2_1);
        x2_2 = u0_x_2(T, t1, k2, x2_1);
        if(abs(x1_2 - L) <= eps && abs(x2_2 - S) <= eps)
            t_t1 = linspace(0, t1, round(t1/t_h + 1));
            t1_T = linspace(t1+t_h, T, round((T - t1)/t_h));
            u1 = u_1(linspace(0, T, t_n), 0, k2, 0, psi20, 1);
            J = trapz(linspace(0, T, t_n), u1 + u1.^2);
            if (iter == 0) || ((iter > 0) && (J < min_J))
                iter = iter + 1;
                min_J = J;
                opt_x1(1 : round(t1/t_h+1)) = x_1(t_t1, 0, k2, 0, psi20, 0, 0);
                opt_x1(round(t1/t_h+2) : t_n) = u0_x_1(t1_T, t1, k2, x1_1, x2_1);
                opt_x2(1 : round(t1/t_h+1)) = x_2(t_t1, 0, k2, 0, psi20, 0);
                opt_x2(round(t1/t_h+2) : t_n) = u0_x_2(t1_T, t1, k2, x2_1);
                opt_u1 = u1;
                opt_u2 = k2;
                opt_psi1 = psi1*ones(1, t_n);
                opt_psi2 = psi_2(linspace(0, T, t_n), 0, k2, psi1, psi20);
                opt_t1 = t1;
                if (J < J_eps)
                    return;
                end
            end
        end
    end
end

%4
%два переключения
disp('2 переключения');
for t1 = linspace(t_h, T, t_n-1)
    disp('----------------------- NEW ITERATION2 ----------------------');
    for t2 = linspace(t1+t_h, T, round(t_n - t1/t_h - 1))
        k_first = k2;
        k_second = k2;
        [psi_1, ex_1] = solve_system(@system1, 0, 2, true);
        k_first = k1;
        k_second = k2;
        [psi2, ex2] = solve_system(@system2, 0, 4, false);
        k_first = k1;
        k_second = k1;
        [psi_3, ex_3] = solve_system(@system1, 0, 2, false);
        for i = 1:ex_1
            psi1 = psi_1(i, 1);
            psi20 = psi_1(i, 2);
            if(psi1>0 && psi20>= 0.5)
               t_t1 = linspace(0, t1, round(t1/t_h + 1));
               t1_t2 = linspace(t1+t_h, t2, round((t2 - t1)/t_h));
               t2_T = linspace(t2+t_h, T, round((T - t2)/t_h));
               x1_1 = x_1(t1, 0, k2, psi1, psi20, 0, 0);
               x2_1 = x_2(t1, 0, k2, psi1, psi20, 0);
               psi21 = psi_2(t1, 0, k2, psi1, psi20);
               psi22 = psi_2(t2, t1, k2, psi1, psi21);
               x1_2 = u0_x_1(t2, t1, k2, x1_1, x2_1);
               x2_2 = u0_x_2(t2, t1, k2, x2_1);
               x1_3 = u0_x_1(T, t2, k1, x1_2, x2_2);
               x2_3 = u0_x_2(T, t2, k1, x2_2);
               disp([abs(x1_3 - L), abs(x2_3 - S)]);
               if (abs(x1_3 - L) <= eps && abs(x2_3 - S) <= eps)
                   u1(1 : round(t2/t_h + 1)) = u_1(linspace(0, t2, round(t2/t_h + 1)), 0, k2, psi1, psi20, 1);
                   u1(round(t2/t_h + 2) : t_n) = u_1(linspace(t2+t_h, T, round((T - t2)/t_h)), t2, k1, psi1, psi22, 1);
                   J = trapz(linspace(0, T, t_n), u1 + u1.^2);
                   if (iter == 0) || ((iter > 0) && (J < min_J))
                       iter = iter + 1;
                       min_J = J; 
                       opt_x1(1 : round(t1/t_h + 1)) = x_1(t_t1, 0, k2, psi1, psi20, 0, 0);
                       opt_x1(round(t1/t_h + 2) : round(t2/t_h + 1)) = u0_x_1(t1_t2, t1, k2, x1_1, x2_1);
                       opt_x1(round(t2/t_h + 2) : t_n) = u0_x_1(t2_T, t2, k1, x1_2, x2_2);
                       
                       opt_x2(1 : round(t1/t_h + 1)) = x_2(t_t1, 0, k2, psi1, psi20, 0);
                       opt_x2(round(t1/t_h + 2) : round(t2/t_h + 1)) = u0_x_2(t1_t2, t1, k2, x2_1);
                       opt_x2(round(t2/t_h + 2) : t_n) = u0_x_2(t2_T, t2, k1, x2_2);
                       opt_u1 = u1;
                       opt_u2(1:round(t2/t_h + 1)) = k2;
                       opt_u2(round(t2/t_h + 2):t_n) = k1;
                       opt_psi2(1:round(t2/t_h + 1)) = psi_2(linspace(0, t2, round(t2/t_h + 1)), 0, k2, psi1, psi20);
                       opt_psi2(round(t2/t_h + 2):t_n) = psi_2(linspace(round(t2+t_h), T, round((T - t2)/t_h)), t2, k1, psi1, psi22);
                       opt_psi1 = psi1*ones(1, t_n);
                       opt_t1 = t1;
                       opt_t2 = t2;
                       if (J < J_eps)
                           return;
                       end
                   end
               end
            end
        end

        for i = 1:ex2
            psi1 = psi2(i, 1);
            psi20 = psi2(i, 2);
            if(psi1<0 && psi20 >= 0.5)
               x1_1 = x_1(t1, 0, k1, psi1, psi20, 0, 0);
               x2_1 = x_2(t1, 0, k1, psi1, psi20, 0);
               psi21 = psi_2(t1, 0, k1, psi1, psi20);
               x1_2 = x_1(t2, t1, k2, psi1, psi21, x1_1, x2_1);
               x2_2 = x_2(t2, t1, k2, psi1, psi21, x2_1);
               psi22 = psi_2(t2, t1, k2, psi1, psi21);
               x1_3 = x_1(T, t2, k1, psi1, psi22, x1_2, x2_2);
               x2_3 = x_2(T, t2, k1, psi1, psi22, x2_2);
               disp([abs(x1_3 - L), abs(x2_3 - S)]);
               if (abs(x1_3 - L) <= eps && abs(x2_3 - S) <= eps)
                   t_t1 = linspace(0, t1, round(t1/t_h + 1));
                   t1_t2 = linspace(t1+t_h, t2, round((t2 - t1)/t_h));
                   t2_T = linspace(t2+t_h, T, round((T - t2)/t_h));                   
                   u1(1 : round(t1/t_h+1)) = u_1(t_t1, 0, k1, psi1, psi20, 1);
                   u1(round(t1/t_h+2) : round(t2/t_h+1)) = u_1(t1_t2, t1, k2, psi1, psi21, 1);
                   u1(round(t2/t_h+2) : t_n) = u_1(t2_T, t2, k1, psi1, psi22, 1);
                   J = trapz(linspace(0, T, t_n), u1 + u1.^2);
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
                       if (J < J_eps)
                           return;
                       end
                   end
               end
            end
        end
        
        for i = 1:ex_3
            psi1 = psi_3(i, 1);
            psi20 = psi_3(i, 2);
            if(psi1 < 0 && psi20 <= 0.5)
               t_t1 = linspace(0, t1, t1/t_h + 1);
               t1_t2 = linspace(t1+t_h, t2, round((t2 - t1)/t_h));
               t2_T = linspace(t2+t_h, T, round((T - t2)/t_h));
               x1_1 = u0_x_1(t1, 0, k1, 0, 0);
               x2_1 = u0_x_2(t1, 0, k1, 0);
               psi21 = psi_2(t1, 0, k1, psi1, psi20);
               x1_2 = x_1(t2, t1, k1, psi1, psi21, x1_1, x2_1);
               x2_2 = x_2(t2, t1, k1, psi1, psi21, x2_1);
               psi22 = psi_2(t2, t1, k1, psi1, psi21);
               x1_3 = x_1(T, t2, k2, psi1, psi22, x1_2, x2_2);
               x2_3 = x_2(T, t2, k2, psi1, psi22, x2_2);
               disp([abs(x1_3 - L), abs(x2_3 - S)]);
               if (abs(x1_3 - L) <= eps && abs(x2_3 - S) <= eps)
                   u1(1 : round(t2/t_h+1)) = u_1(linspace(0, t2, round(t2/t_h+1)), 0, k1, psi1, psi20, 1);
                   u1(round(t2/t_h+2) : t_n) = u_1(t2_T, t2, k2, psi1, psi22, 1);
                   J = trapz(linspace(0, T, t_n), u1 + u1.^2);
                   if (iter == 0) || ((iter > 0) && (J < min_J))
                       iter = iter + 1;
                       min_J = J;
                       opt_x1(1 : round(t1/t_h+1)) = u0_x_1(t_t1, 0, k1, 0, 0);
                       opt_x1(round(t1/t_h+2) : round(t2/t_h+1)) = x_1(t1_t2, t1, k1, psi1, psi21, x1_1, x2_1);
                       opt_x1(round(t2/t_h+2) : t_n) = x_1(t2_T, t2, k2, psi1, psi22, x1_2, x2_2);
                       
                       opt_x2(1 : round(t1/t_h+1)) = u0_x_2(t_t1, 0, k1, 0);
                       opt_x2(round(t1/t_h+2) : round(t2/t_h+1)) = x_2(t1_t2, t1, k1, psi1, psi21, x2_1);
                       opt_x2(round(t2/t_h+2) : t_n) = x_2(t2_T, t2, k1, psi1, psi22, x2_2);
                       opt_u1 = u1;
                       opt_u2(1 : round(t2/t_h+1)) = k1;
                       opt_u2(round(t2/t_h+2) : t_n) = k2;
                       opt_t1 = t1;
                       opt_t2 = t2;
                       opt_psi1 = psi1*ones(1, t_n);
                       opt_psi2(1:round(t2/t_h + 1)) = psi_2(linspace(0, t2, round(t2/t_h + 1)), 0, k1, psi1, psi20);
                       opt_psi2(round(t2/t_h + 2) : t_n) = psi_2(t2_T, t2, k2, psi1, psi22);
                       if (J < J_eps)
                           return;
                       end                       
                   end
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

