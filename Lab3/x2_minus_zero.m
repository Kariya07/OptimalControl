function [value, isterminal, direction] = x2_minus_event(~, x)
    value = x(2);
    isterminal = 1;
    direction = 1;
end

