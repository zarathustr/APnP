function [R_, t_, X, xs, min_val] = symbolic_pnp(image_pt, world_pt, K)
    syms q0 q1 q2 q3 t0 t1 t2 lambda real
    x = [q0; q1; q2; q3; t0; t1; t2; lambda];
    t = [t0; t1; t2];
    R = [
        q0^2 + q1^2 - q2^2 - q3^2,         2*q0*q3 + 2*q1*q2,         2*q1*q3 - 2*q0*q2;
                2*q1*q2 - 2*q0*q3, q0^2 - q1^2 + q2^2 - q3^2,         2*q0*q1 + 2*q2*q3;
                2*q0*q2 + 2*q1*q3,         2*q2*q3 - 2*q0*q1, q0^2 - q1^2 - q2^2 + q3^2];
            
    J = J_pnp_loss(image_pt, world_pt, K, R, t) * 1e-3;
    Jacob = expand(vpa(jacobian(J, x), 32));
    Jacob(1 : 4) = Jacob(1 : 4) + lambda * x(1 : 4)';
    ts = vpasolve(Jacob(5 : 7) == 0, t);
    t0 = ts.t0;
    t1 = ts.t1;
    t2 = ts.t2;
    t0_func = matlabFunction(ts.t0);
    t1_func = matlabFunction(ts.t1);
    t2_func = matlabFunction(ts.t2);
    new_eq_Jacob = vpa(expand(eval(Jacob(1 : 4))), 32);
    xs = [];
    x_count = 1;
    eqs = vpa(expand([new_eq_Jacob.'; x(1 : 4).' * x(1 : 4) - 1]), 32);
    
%     bound = 1e-5;
%     precision = 64;
%     
%     for k = 1 : 1
%         str = sprintf('s := numeric::fsolve([%s = 0, %s = 0, %s = 0, %s = 0, %s = 0], [q0 = -1..1, q1 = -1..1, q2 = -1..1, q3 = -1..1, lambda = %f..%f], MultiSolutions); Pref::outputDigits(%d): float(s)', ...
%                   char(eqs(1)), ...
%                   char(eqs(2)), ...
%                   char(eqs(3)), ...
%                   char(eqs(4)), ...
%                   char(eqs(5)), ...
%                   - bound, ...
%                   + bound, ...
%                   precision);
%         s1 = evalin(symengine, str)
%     
%         for i = 1 : 1
%             str = sprintf('s = s%d;', i);
%             eval(str);
%         
%             str_s = char(s);
%             if(strncmp(str_s, 'matrix([[[', 10))
%                 len = size(s, 2);
%                 for j = 1 : len
%                     x = s(j);
%                     q0 = sscanf(char(x(1)), 'q0 == %f');
%                     q1 = sscanf(char(x(2)), 'q1 == %f');
%                     q2 = sscanf(char(x(3)), 'q2 == %f');
%                     q3 = sscanf(char(x(4)), 'q3 == %f');
%                     lambda = sscanf(char(x(5)), 'lambda == %f');
%                     t0 = t0_func(q0, q1, q2, q3);
%                     t1 = t1_func(q0, q1, q2, q3);
%                     t2 = t2_func(q0, q1, q2, q3);
%                     xx = [q0, q1, q2, q3, t0, t1, t2, lambda];
%                     xx(1 : 4) = xx(1 : 4) ./ norm(xx(1 : 4));
%                     xs = [xs; xx];
%                     x_count = x_count + 1;
%                 end
%             elseif(~strncmp(str_s, 'FAIL', 4))
%                 q0 = sscanf(char(s(1)), 'q0 == %f');
%                 q1 = sscanf(char(s(2)), 'q1 == %f');
%                 q2 = sscanf(char(s(3)), 'q2 == %f');
%                 q3 = sscanf(char(s(4)), 'q3 == %f');
%                 lambda = sscanf(char(s(5)), 'lambda == %f');
%                 t0 = t0_func(q0, q1, q2, q3);
%                 t1 = t1_func(q0, q1, q2, q3);
%                 t2 = t2_func(q0, q1, q2, q3);
%                 xx = [q0, q1, q2, q3, t0, t1, t2, lambda];
%                 xx(1 : 4) = xx(1 : 4) ./ norm(xx(1 : 4));
%                 xs = [xs; xx];
%                 x_count = x_count + 1;
%             end
%         end
%     end
%     xs_ = xs;
    
%     lambda0 = 1e-3;
%     assumeAlso(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3 == 1);
%     eqs_ = vpa(expand(simplify(eval(subs(eqs, lambda, lambda0)))), 32);
%     str = '';
%     for i = 1 : length(eqs) - 1
%         str = strcat(str, sprintf(' PP{%d} = char(vpa(%%s, 32));', i));
%     end
%     
%     str_ = sprintf(str, char(eqs_(1)), ...
%                         char(eqs_(2)), ...
%                         char(eqs_(3)), ...
%                         char(eqs_(4)));
%     eval(str_);                    
%     [S, vars] = psolve(PP);
%     S = S.';
%     SS = S;
%     for i = 1 : length(vars)
%         if(strcmp(vars{i}, 'q0'))
%             SS(:, 1) = S(:, i);
%         elseif(strcmp(vars{i}, 'q1'))
%             SS(:, 2) = S(:, i);
%         elseif(strcmp(vars{i}, 'q2'))
%             SS(:, 3) = S(:, i);
%         elseif(strcmp(vars{i}, 'q3'))
%             SS(:, 4) = S(:, i);
%         end
%     end
%     S = real(SS);
% 
%     
%     
    str = '';
    for i = 1 : length(eqs)
        str = strcat(str, sprintf(' PP{%d} = char(vpa(%%s, 32));', i));
    end
    
    str_ = sprintf(str, char(eqs(1)), ...
                        char(eqs(2)), ...
                        char(eqs(3)), ...
                        char(eqs(4)), ...
                        char(eqs(5)));
    eval(str_);                    
    [S, vars] = psolve(PP);
    S = S.';
    SS = S;
    for i = 1 : length(vars)
        if(strcmp(vars{i}, 'q0'))
            SS(:, 1) = S(:, i);
        elseif(strcmp(vars{i}, 'q1'))
            SS(:, 2) = S(:, i);
        elseif(strcmp(vars{i}, 'q2'))
            SS(:, 3) = S(:, i);
        elseif(strcmp(vars{i}, 'q3'))
            SS(:, 4) = S(:, i);
        elseif(strcmp(vars{i}, 'lambda'))
            SS(:, 5) = S(:, i);
        end
    end
    S = real(SS);
    xs_ = S;
    x_count = length(xs_(:, 1)) + 1;
    xs = zeros(x_count - 1, 8);
    
    for i = 1 : x_count - 1
        q0 = xs_(i, 1);
        q1 = xs_(i, 2);
        q2 = xs_(i, 3);
        q3 = xs_(i, 4);
        N = sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
        q0 = q0 ./ N;
        q1 = q1 ./ N;
        q2 = q2 ./ N;
        q3 = q3 ./ N;
        lambda = xs_(i, 5);
        t0 = t0_func(q0, q1, q2, q3);
        t1 = t1_func(q0, q1, q2, q3);
        t2 = t2_func(q0, q1, q2, q3);
        xs(i, :) = [q0, q1, q2, q3, t0, t1, t2, lambda];
    end
    
    Ls = zeros(x_count - 1, 1);
    for i = 1 : x_count - 1
        RR = quat2dcm(xs(i, 1 : 4));
        tt = xs(i, 5 : 7)';
        Ls(i) = J_pnp_loss(image_pt, world_pt, K, RR, tt);
    end
    [minimum, idx] = sort(Ls);
    R_ = quat2dcm(xs(idx(1), 1 : 4))';
    t_ = xs(idx(1), 5 : 7)';
%     t_ = - R_ * t_;
%     R_ = R_';
    X = [R_', t_;
          zeros(1, 3), 1];
    min_val = minimum;
end