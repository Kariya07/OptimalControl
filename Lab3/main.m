%% Анализ неподвижных точек
alpha = 1.5;
f1 = @(x) 2*x + x.*sin(x.^2) - 2*x.^2.*cos(x) - alpha;
f2 = @(x) 2*x + x.*sin(x.^2) - 2*x.^2.*cos(x) + alpha;
point1 = fzero(f1, 0);
point2 = fzero(f2, 0);
%%
k = @(x1) 2 + sin(x1^2) + 2*x1^2*cos(x1^2) - 4*x1*cos(x1) + 2*x1^2*sin(x1);
k1 = k(point1);
k2 = k(point2);
lam = linspace(-10, 10, 1000);
plot(lam, lam.^2 + lam + k1, 'r', lam, lam.^2 + lam + k2, 'b');
title('\alpha = 5');
xlabel('\lambda');
ylabel('det(J - I \lambda)');

%%
x = linspace(-10, 10, 10000);
stab1 = f1(x);
stab2 = f2(x);
eps = 0.1;
%points = find(abs(stab1) < eps);
plot(x, stab1,'r', x, stab2, 'b');
%% 
t_size = 500;
T = 5;
alpha = 1.5;
[X, Y, switchx, switchy] = reachset(alpha, T, t_size);
plot(X, Y, 'black', 'LineWidth', 1.5);
hold on
plot([0 switchx(1,1:t_size)], [0 switchy(1,1:t_size)], 'r', [0 switchx(1,t_size + 1:2*t_size)], [0 switchy(1,t_size + 1:2*t_size)], 'b', 'LineWidth', 1);
hold on
plot(switchx(2,1:t_size), switchy(2,1:t_size), 'b', switchx(2,t_size + 1:2*t_size), switchy(2,t_size + 1:2*t_size), 'r', 'LineWidth', 1);
hold on;
%plot(switchx(3,1:t_size), switchy(3,1:t_size), 'r', switchx(3,t_size + 1:2*t_size), switchy(3,t_size + 1:2*t_size), 'b', 'LineWidth', 1);
%hold on;
plot(point1, 0, 'r*', point2, 0, 'b*');
hold on;
xlabel('x_1');
ylabel('x_2');
title(['\alpha = ', num2str(alpha), '  T = ', num2str(T)]);
%%
mov = reachset_dyn(0.5, 2, 5, 20, 'res2.avi');
%%
movie(mov, 1, 10);