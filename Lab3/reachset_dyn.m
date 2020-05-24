function res = reachset_dyn(alpha, t1, t2, n, filename)
    if n == 0
        t_size = 1000;
        [X, Y, switchx, switchy] = reachset(alpha, t2, t_size);
        plot(X, Y, 'black');
        hold on;
        plot(switchx(1,1:t_size), switchy(1,1:t_size), 'r', switchx(1,t_size + 1:2*t_size), switchy(1,t_size + 1:2*t_size), 'b');
        hold on;
        plot(switchx(2,1:t_size), switchy(2,1:t_size), 'b', switchx(2,t_size + 1:2*t_size), switchy(2,t_size + 1:2*t_size), 'r');
        hold on;
        plot(switchx(3,1:t_size), switchy(3,1:t_size), 'r', switchx(3,t_size + 1:2*t_size), switchy(3,t_size + 1:2*t_size), 'b');
        hold on;
        axis([-1.1 1.1 -1.2 1.2])
        res(1) = getframe;
        cla;
        if nargin == 5
            v = VideoWriter(filename);
            open(v);
            writeVideo(v, res);
            close(v);
        end
        return
    end
    res(1:n) = struct('cdata', [], 'colormap', []);
    if ~t1
        T = linspace(t1, t2, n+1);
        T(1) = [];
    else
        T = linspace(t1, t2, n);
    end
    for j = 1:n
        t_size = 100;
        [X, Y, switchx, switchy] = reachset(alpha, T(j), t_size);
        plot(X, Y, 'black');
        hold on;
        plot([0 switchx(1,1:t_size)], [0 switchy(1,1:t_size)], 'r', [0 switchx(1,t_size + 1:2*t_size)], [0 switchy(1,t_size + 1:2*t_size)], 'b');
        hold on;
        plot(switchx(2,1:t_size), switchy(2,1:t_size), 'b', switchx(2,t_size + 1:2*t_size), switchy(2,t_size + 1:2*t_size), 'r');
        hold on;
        plot(switchx(3,1:t_size), switchy(3,1:t_size), 'r', switchx(3,t_size + 1:2*t_size), switchy(3,t_size + 1:2*t_size), 'b');
        hold on;
        axis([-1.1 1.1 -1.2 1.2])
        res(j) = getframe;
        cla;
    end
    if nargin == 5
        v = VideoWriter(filename);
        open(v);
        writeVideo(v, res);
        close(v);
    end
end
