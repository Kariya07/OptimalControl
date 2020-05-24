function [end_of_traj, sws] = find_trajectory_minus(alpha, tspan, x0, sw_cur, num_of_sw)
    if num_of_sw <= 4
        sw_cur(num_of_sw, 1) = x0(1);
        sw_cur(num_of_sw, 2) = x0(2);
    end
    opt = odeset('Events', @psi2_zero_minus);
    [~, x, te, ye, ~] = ode45(@(t, x) ode_minus(t, x, alpha), tspan, x0, opt);
    plot(x(:, 1), x(:,2), 'b');
    hold on;
    if isempty(te)
        end_of_traj = [x(end, 1), x(end, 2)];
        sws = sw_cur;
    else
        %disp(ye);
        %disp('--------------------------------------------------------');
        
        [end_of_traj, sws] = find_trajectory_plus(alpha, [te(1), tspan(end)], [ye(1, 1) ye(1, 2) ye(1, 3) ye(1, 4)], sw_cur, num_of_sw + 1);
    end
end

