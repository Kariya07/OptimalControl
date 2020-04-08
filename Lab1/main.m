% число точек для отрисовки одного множества
drawSet_n = 200; 
%размер сетки для psi_0
psi_0_n = 50;
%шаг сетки по времени
h_t = 0.001;
n_l = 128;
T_max = 4;
%для квадрата и отрисовки точки
eps = 0.05;
%для B
B_eps = 0.000001;
eps_B = 0.1;
%% Входные параметры1
A = -[-1 0; 0 -1]/10;
B = -[10 0; 0 10];
f = -[0; 0];
p1 = 0;
p2 = 0;
alpha = 2;
beta = 1/2;
gamma = 1;
delta = 1;
k = 2;
x0 = [7 4]';
x1 = [6 4]';
while det(B) < B_eps
    B = zero_B(B, B_eps);
end
%% Входные параметры 2
A = -[-1 -3; 3 1];
B = -[3 0; 0 3];
f = -[2; -2];
p1 = 1;
p2 = 1;
alpha = 1;
beta = 3;
gamma = 1;
delta = 1;
k = 3;
x1 = [5 -1]';
x0 = [-2 2]';
while det(B) < B_eps
    B = zero_B(B);
end
%% Входные параметры 3
A = -[1 3; -3 -1];
B = -[-2 0; 0 -2]/5;
%f = -[-7; 7];
f = -[0; 0];
p1 = 1;
p2 = 1;
alpha = 2;
beta = 1/2;
gamma = 1;
delta = 1;
k = 2;
%x0 = [-1.39 3]';
x1 = [1 1]';
x0 = [-3.4 2]';
while det(B) < B_eps
    B = zero_B(B);
end
%%
clc;
phi = linspace(0, 2 * pi, psi_0_n + 1);
phi(psi_0_n + 1) = [];
psi_0 = [cos(phi); sin(phi)];
[X, U, Psi, opt_ind, t1, t] = traj(psi_0, psi_0_n, A, B, f, h_t, T_max, p1, p2, alpha, beta, gamma, delta, k, x0, x1, eps);

%% Условие трансверсальности
x_1 = X(:, end, opt_ind);
psi1 = Psi(:, end, opt_ind);
psi1 = psi1 / norm(psi1);
l1 = find_dir(x_1, @(l) rhoX0(l, k, x0, eps), n_l);
angle = acos(-l1' * psi1) * 180 / pi;
disp('Inaccuracy = ');
disp(angle);

%% Локализация
phi_step = phi(2) - phi(1);
phi = linspace(phi(opt_ind) - phi_step * psi_0_n / 10, phi(opt_ind) + phi_step * psi_0_n / 10, psi_0_n);
psi_0 = [cos(phi); sin(phi)];

[X, U, Psi, opt_ind, t1, t] = traj(psi_0, psi_0_n, A, B, f, h_t, T_max, p1, p2, alpha, beta, gamma, delta, k, x0, x1, eps);

%% Оптимальная траектория
eps1 = 0.1;
drawSet(@(l) rhoX0(l, k, x0, 0), drawSet_n, 'm');
hold on;
drawSet(@(l) rhoX1(l, x1, eps1), drawSet_n, 'g');
xlabel('x_1');
ylabel('x_2');
title('(x_1, x_2) - Оптимальная траектория')
axis equal;

if opt_ind
    plot(X(1, :, opt_ind), X(2, :, opt_ind), 'y', 'LineWidth', 2);
end
%Отрисовка вектора пси и опорной гиперплоскости для квадрата
x_0 = X(:, size(X, 2), opt_ind);
if (~in_square(x_0, x0, k, 1))
    psi1 = Psi(:, size(Psi, 2), opt_ind);
    psi1 = psi1 / norm(psi1);
    l1 = find_dir(x_0, @(l) rhoX0(l, k, x0, eps), n_l);
    tau = [-l1(2); l1(1)];
    plot(x_0(1) + linspace(-2, 2, drawSet_n) * tau(1), x_0(2) + linspace(-2, 2, drawSet_n) * tau(2), 'black', 'LineWidth', 1);
    hold on;
    plot(x_0(1) + linspace(0, 1, drawSet_n) * (-psi1(1)), x_0(2) + linspace(0, 1, drawSet_n) * (-psi1(2)), 'red', 'LineWidth', 1);
    hold on;
    plot(x_0(1) + linspace(0, 1, drawSet_n) * (l1(1)), x_0(2) + linspace(0, 1, drawSet_n) * (l1(2)), 'blue', 'LineWidth', 1);
    hold on;
    legend('X_0', 'X_1', 'опт. траектория', 'опорн. гиперпл.', '-\psi_1', 'нормаль');
else
    legend('X_0', 'X_1');
end
%% Все траектории
drawSet(@(l) rhoX1(l, x1, eps1), drawSet_n, 'g');
hold on;
drawSet(@(l) rhoX0(l, k, x0, 0), drawSet_n, 'm');
xlabel('x_1');
ylabel('x_2');
title('(x_1, x_2) - Все траектории')
axis equal;
for i = 1:psi_0_n
    plot(X(1, :, i), X(2, :, i), 'b');
    hold on;
end

if opt_ind
    plot(X(1, :, opt_ind), X(2, :, opt_ind), 'y', 'LineWidth', 2);
end

x_0 = X(:, size(X, 2), opt_ind);
if (~in_square(x_0, x0, k, 1))
    psi1 = Psi(:, size(Psi, 2), opt_ind);
    psi1 = psi1 / norm(psi1);
    l1 = find_dir(x_0, @(l) rhoX0(l, k, x0, eps), n_l);
    tau = [-l1(2); l1(1)];
    plot(x_0(1) + linspace(-2, 2, drawSet_n) * tau(1), x_0(2) + linspace(-2, 2, drawSet_n) * tau(2), 'black', 'LineWidth', 1);
    hold on;
    plot(x_0(1) + linspace(0, 1, drawSet_n) * (-psi1(1)), x_0(2) + linspace(0, 1, drawSet_n) * (-psi1(2)), 'red', 'LineWidth', 1);
    hold on;
    plot(x_0(1) + linspace(0, 1, drawSet_n) * (l1(1)), x_0(2) + linspace(0, 1, drawSet_n) * (l1(2)), 'blue', 'LineWidth', 1);
    hold on;
end
%% (t, x1)
axis equal;
for i = 1:psi_0_n
    plot(t, X(1, end:-1:1, i), 'b');
    hold on;
end

if opt_ind
    plot(t, X(1, end:-1:1, opt_ind), 'y', 'LineWidth', 2);
end
xlabel('t');
ylabel('x1');
title('График (t, x_1)')
%% (t, x2)
axis equal;

for i = 1:psi_0_n
    plot(t, X(2, end:-1:1, i), 'b');
    hold on;
end

if opt_ind
    plot(t, X(2, end:-1:1, opt_ind), 'y', 'LineWidth', 2);
end
xlabel('t');
ylabel('x2');
title('График (t, x_2)')
%% (u1, u2)
axis equal;
drawSet(@(l) rhoP(l, p1, p2, alpha, beta, gamma, delta), drawSet_n, 'c');
hold on;
plot(U(1, :, opt_ind), U(2, :, opt_ind), 'or', 'MarkerSize', 5);
hold on;
xlabel('u_1');
ylabel('u_2');
title ('(u_1, u_2) - Перебор управлений')
%% (t, u1)
axis equal;

for i = 1:psi_0_n
    plot(t, U(1, end:-1:1, i), 'b');
    hold on;
end

if opt_ind
    plot(t, U(1, end:-1:1, opt_ind), 'y', 'LineWidth', 2);
end
xlabel('t');
ylabel('u1');
title('График (t, u_1)')
%% (t, u2)
axis equal;
for i = 1:psi_0_n
    plot(t, U(2, end:-1:1, i), 'b');
    hold on;
end

if opt_ind
    plot(t, U(2, end:-1:1, opt_ind), 'y', 'LineWidth', 2);
end
xlabel('t');
ylabel('u2');
title('График (t, u_2)');
%% (psi1, psi2)
for i = 1:psi_0_n
    plot(Psi(1, :, i), Psi(2, :, i), 'b');
    hold on;
end

if opt_ind
    plot(Psi(1, :, opt_ind), Psi(2, :, opt_ind), 'y', 'LineWidth', 2);
end
xlabel('\psi_1');
ylabel('\psi_2');
title('График (\psi_1, \psi_2)')
%% (t, psi1)
axis equal;

for i = 1:psi_0_n
    plot(t, Psi(1, end:-1:1, i), 'b');
    hold on;
end

if opt_ind
    plot(t, Psi(1, end:-1:1, opt_ind), 'y', 'LineWidth', 2);
end
xlabel('t');
ylabel('\psi_1');
title('График (t, \psi_1)')
%% (t, psi2)
axis equal;

for i = 1:psi_0_n
    plot(t, Psi(2, end:-1:1, i), 'b');
    hold on;
end

if opt_ind
    plot(t, Psi(2, end:-1:1, opt_ind), 'y', 'LineWidth', 2);
end
xlabel('t');
ylabel('\psi_2');
title('График (t, \psi_2)')