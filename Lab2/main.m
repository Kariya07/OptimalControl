%% Два переключения
k1 = -1;
k2 = 3;
L = 8;
S = 2.5;
T = 3;
eps = 1.5;
t_h = 0.1;
%J = 0.0103
%%
[J, opt_x1, opt_x2, opt_u1, opt_u2, opt_psi1, opt_psi2, opt_t1, opt_t2] = task1(k1, k2, T, L, S, eps, t_h);
%%
t_n = T/t_h + 1;
ind1 = round(t_n - (T - opt_t1)/t_h);
ind2 = round(t_n - (T - opt_t2)/t_h);
%%
t = linspace(0, T, T/t_h + 1);
plot(opt_x1, opt_x2, 0, 0, 'og', opt_x1(end), opt_x2(end), 'or', opt_x1(ind1), opt_x2(ind1), '*r', opt_x1(ind2), opt_x2(ind2), '*r');
xlabel('x_1');
ylabel('x_2');
legend('trajectory', 'begin', 'end');
%%
plot(t, opt_x1, t, opt_x2, opt_t1, opt_x1(ind1), '*r', opt_t1, opt_x2(ind1), '*r', ...
opt_t2, opt_x1(ind2), '*r', opt_t2, opt_x2(ind2), '*r');
xlabel('t');
legend('x_1', 'x_2', 'Location', 'west');
%%
plot(t, opt_psi1, t, opt_psi2, opt_t1, opt_psi1(ind1), '*r', opt_t1, opt_psi2(ind1), '*r', ...
opt_t2, opt_psi1(ind2), '*r', opt_t2, opt_psi2(ind2), '*r');
xlabel('t');
legend('\psi_1', '\psi_2', 'Location', 'west');
%%
%opt_u2 = opt_u2.*ones(1, t_n);
plot(t, opt_u1, t, opt_u2, opt_t1, opt_u1(ind1), '*r', opt_t1, opt_u2(ind1), '*r', ...
opt_t2, opt_u1(ind2), '*r', opt_t2, opt_u2(ind2), '*r');
xlabel('t');
legend('u_1', 'u_2', 'Location', 'northwest');
xlim([0, T]);
ylim([-1.3, 3.3])
%%
plot(0, opt_psi2(1), 'og', opt_x2(end), opt_psi2(end), 'or', opt_x2, opt_psi2, opt_x2(ind1), opt_psi2(ind1), '*r', opt_x2(ind2),...
    opt_psi2(ind2), '*r');
ylabel('\psi_2');
xlabel('x_2');
legend('begin', 'end');

%% Случайные переключения по u2
k1 = -1;
k2 = 3;
L = 6;
S = 0;
T = 4;
eps = 0.5;
t_h = 0.1;
%J = 0;
[J, opt_x1, opt_x2, opt_u1, opt_u2, opt_psi1, opt_psi2, opt_t1, opt_t2] = task1(k1, k2, T, L, S, eps, t_h);
%%
t = linspace(0, T, T/t_h + 1);
plot(opt_x1, opt_x2, 0, 0, 'og', opt_x1(end), opt_x2(end), 'or');
xlabel('x_1');
ylabel('x_2');
legend('trajectory', 'begin', 'end');
%%
plot(t, opt_x1, t, opt_x2);
xlabel('t');
legend('x_1', 'x_2', 'Location', 'northwest');
ylim([-1, 7]);
%%
plot(t, opt_u1, t, opt_u2);
xlabel('t');
legend('u_1', 'u_2', 'Location', 'northwest');
xlim([0, T]);
ylim([-1.1, 3.1])



%% Без переключений
k1 = -5;
k2 = 2;
L = -35;
S = 0.3;
T = 7;
eps = 0.5;
t_h = 0.1;
[J, opt_x1, opt_x2, opt_u1, opt_u2, opt_psi1, opt_psi2, opt_t1, opt_t2] = task1(k1, k2, T, L, S, eps, t_h);
%%
t = linspace(0, T, T/t_h + 1);
plot(opt_x1, opt_x2, 0, 0, 'og', opt_x1(end), opt_x2(end), 'or');
xlabel('x_1');
ylabel('x_2');
legend('trajectory', 'begin', 'end');
%%
plot(t, opt_x1, t, opt_x2);
xlabel('t');
legend('x_1', 'x_2');
ylim([-35, 2]);
%%
plot(t, opt_u1, t, opt_u2);
xlabel('t');
legend('u_1', 'u_2', 'Location', 'northwest');
xlim([0, T]);
ylim([-6, 1]);

%% Задача Б
%%Одно переключение
k1 = -2;
k2 = 1;
L = 6;
S = 2;
T = 2;
eps = 2;
t_h = 0.1;
%J = 16.
%%
[J2, opt_x12, opt_x22, opt_u12, opt_u22, opt_psi12, opt_psi22, opt_t12, opt_t22] = task2(k1, k2, T, L, S, eps, t_h);
%%
t_n = T/t_h + 1;
ind1 = round(t_n - (T - opt_t12)/t_h);
ind2 = round(t_n - (T - opt_t22)/t_h);
%%
t = linspace(0, T, T/t_h + 1);
plot(opt_x12, opt_x22, 0, 0, 'og', opt_x12(end), opt_x22(end), 'or', ...
    opt_x12(ind1), opt_x22(ind1), '*r', opt_x12(ind2), opt_x22(ind2), '*r');
xlabel('x_1');
ylabel('x_2');
legend('trajectory', 'begin', 'end');
%%
plot(t, opt_x12, t, opt_x22, opt_t12, opt_x12(ind1), '*r', opt_t12, opt_x22(ind1), '*r', ...
opt_t22, opt_x12(ind2), '*r', opt_t22, opt_x22(ind2), '*r');
xlabel('t');
legend('x_1', 'x_2', 'Location', 'west');
%%
plot(t, opt_psi12, t, opt_psi22, opt_t12, opt_psi12(ind1), '*r', opt_t12, opt_psi22(ind1), '*r', ...
opt_t22, opt_psi12(ind2), '*r', opt_t22, opt_psi22(ind2), '*r');
xlabel('t');
legend('\psi_1', '\psi_2', 'Location', 'west');
ylim([-0.2, 2])
%%
plot(t, opt_u12, t, opt_u22, opt_t12, opt_u12(ind1), '*r', opt_t12, opt_u22(ind1), '*r', ...
opt_t22, opt_u12(ind2), '*r', opt_t22, opt_u22(ind2), '*r');
xlabel('t');
legend('u_1', 'u_2', 'Location', 'northwest');
xlim([0, T]);
%%
plot(0, opt_psi22(1), 'og', opt_x22(end), opt_psi22(end), 'or', opt_x22, opt_psi22, opt_x22(ind1), opt_psi22(ind1), '*r', opt_x22(ind2),...
    opt_psi22(ind2), '*r');
ylabel('\psi_2');
xlabel('x_2');
legend('begin', 'end');
%% Два переключения
k1 = -1;
k2 = 3;
L = -10;
S = -4;
T = 2;
eps = 1;
t_h = 0.1;
%J = -0.5
%%
[J, opt_x1, opt_x2, opt_u1, opt_u2, opt_psi1, opt_psi2, opt_t1, opt_t2] = task2(k1, k2, T, L, S, eps, t_h);
%%
t_n = T/t_h + 1;
ind1 = round(t_n - (T - opt_t1)/t_h);
ind2 = round(t_n - (T - opt_t2)/t_h);
%%
t = linspace(0, T, T/t_h + 1);
plot(opt_x1, opt_x2, 0, 0, 'og', opt_x1(end), opt_x2(end), 'or', ...
    opt_x1(ind1), opt_x2(ind1), '*r', opt_x1(ind2), opt_x2(ind2), '*r');
xlabel('x_1');
ylabel('x_2');
legend('trajectory', 'begin', 'end');
%%
plot(t, opt_x1, t, opt_x2, opt_t1, opt_x1(ind1), '*r', opt_t1, opt_x2(ind1), '*r', ...
opt_t2, opt_x1(ind2), '*r', opt_t2, opt_x2(ind2), '*r');
xlabel('t');
legend('x_1', 'x_2', 'Location', 'west');
%%
plot(t, opt_psi1, t, opt_psi2, opt_t1, opt_psi1(ind1), '*r', opt_t1, opt_psi2(ind1), '*r', ...
opt_t2, opt_psi1(ind2), '*r', opt_t2, opt_psi2(ind2), '*r');
xlabel('t');
legend('\psi_1', '\psi_2', 'Location', 'west');
%%
plot(t, opt_u1, t, opt_u2, opt_t1, opt_u1(ind1), '*r', opt_t1, opt_u2(ind1), '*r', ...
opt_t2, opt_u1(ind2), '*r', opt_t2, opt_u2(ind2), '*r');
xlabel('t');
legend('u_1', 'u_2', 'Location', 'northwest');
ylim([-2, 4]);
%%
plot(0, opt_psi2(1), 'og', opt_x2(end), opt_psi2(end), 'or', opt_x2, opt_psi2, opt_x2(ind1), opt_psi2(ind1), '*r', opt_x2(ind2),...
    opt_psi2(ind2), '*r');
ylabel('\psi_2');
xlabel('x_2');
legend('begin', 'end');
%%
k1 = -1;
k2 = 1;
L = 6;
S = 5;
T = 2;
eps = 2;
t_h = 0.1;
    