function [val, point] = rhoP(l, p1, p2, alpha, beta, gamma, delta)
    if abs(l(2)) <= 2*sqrt(gamma)*abs(l(1)) && l(1) > 0
        val = l(1)*p1 + l(2)*p2 + l(1) * (gamma - l(2)^2/(4*l(1)^2)) / alpha + l(2)^2 / (2*sqrt(delta)*l(1));
        point = [(gamma - l(2)^2 / (4*l(1)^2)) / alpha, l(2) / (2*l(1)*sqrt(delta))]';
    else
        if abs(l(2)) <= 2*sqrt(gamma)*abs(l(1)) && l(1) < 0
            val = l(1)*p1 + l(2)*p2 + l(1) * (-gamma + l(2)^2/(4*l(1)^2)) / beta - l(2)^2 / (2*sqrt(delta)*l(1));
            point = [(-gamma + l(2)^2 / (4*l(1)^2)) / beta, -l(2) / (2*l(1)*sqrt(delta))]';
        else
            if l(2) >= 0
                val = l(1)*p1 + l(2)*p2 + sqrt(gamma/delta) * l(2);
                point = [0, sqrt(gamma/delta)]';
            else
                val = l(1)*p1 + l(2)*p2 - sqrt(gamma/delta) * l(2);
                point = [0, -sqrt(gamma/delta)]';
            end
        end
    end
    point = point + [p1, p2]';
end

