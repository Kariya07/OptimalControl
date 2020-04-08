function [val, point] = rhoX0_1(l, k, x0, eps)
    a1 = 1 / (k/2 + eps)^2;
    b1 = 1 / (k*(k/2 + eps)/(2*sqrt(eps^2 + k*eps)))^2;
    if a1 < 1 && b1 < 1
        a2 = b1;
        b2 = a1;
        if abs(l(1)) > abs(l(2)) && (sqrt(a1/b1)*abs(l(2))/sqrt(b1*l(1)^2 + a1*l(2)^2)) < k/2
            val = (b1*l(1)^2 + a1*l(2)^2)/((sqrt(a1*b1)*sqrt(b1*l(1)^2 + a1*l(2)^2)));
            point = [sqrt(b1/a1)*l(1), sqrt(a1/b1)*l(2)]' / sqrt(b1*l(1)^2 + a1*l(2)^2);
        else
            if abs(l(1)) < abs(l(2)) && (sqrt(b2/a2)*abs(l(1))/sqrt(b2*l(1)^2 + a2*l(2)^2)) < k/2
                val = (b2*l(1)^2 + a2*l(2)^2)/((sqrt(a2*b2)*sqrt(b2*l(1)^2 + a2*l(2)^2)));
                point = [sqrt(b2/a2)*l(1), sqrt(a2/b2)*l(2)]' / sqrt(b2*l(1)^2 + a2*l(2)^2);
            else
                val = (k/2)*(abs(l(1)) + abs(l(2)));
                point = (k/2) * [sign(l(1)), sign(l(2))]';
            end
        end
        val = val + l(1) * x0(1) + l(2) * x0(2);
        point = point + x0;
    else
        print('a1 > 1')
    end
end
            
