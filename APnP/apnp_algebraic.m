function [R_, t_, s_, xs, min_val] = apnp_algebraic(bbb, rrr, Xw, U, K)
    N = size(bbb, 1);
    b = zeros(3, N);
    r = zeros(3, N);
    for i = 1 : N
        b(:, i) = bbb(i, 1 : 3)';
        r(:, i) = rrr(i, 1 : 3)';
    end
    b_bar = zeros(3, 1);
    r_bar = zeros(3, 1);
    for i = 1 : N
        b_bar = b_bar + 1 / N * b(:, i); 
        r_bar = r_bar + 1 / N * r(:, i); 
    end
    
    P = zeros(4, 4);
    for i = 1 : N
        for j = 1 : 3
            str = sprintf('M = M%d_matrix(b(:, i) - b_bar);', j);
            eval(str);
            P = P + 1 / N * r(j, i) * M;
        end
    end
    [V, D] = eig(P);

    q11 = V(:, 1); q_1 = q11 ./ norm(q11);
    q22 = V(:, 2); q_2 = q22 ./ norm(q22);
    q33 = V(:, 3); q_3 = q33 ./ norm(q33);
    q44 = V(:, 4); q_4 = q44 ./ norm(q44);
    x_count = 4 + 1;
    xs_ = [
        q_1';
        q_2';
        q_3';
        q_4';
        ];
    
    
    for i = 1 : x_count - 1
        q0 = xs_(i, 1);
        q1 = xs_(i, 2);
        q2 = xs_(i, 3);
        q3 = xs_(i, 4);
        R = quat2dcm_([q0, q1, q2, q3]);
        aa = 0;
        bb = 0;
        for j = 1 : N
            aa = aa + 1 / N * (b(:, j)' * (b(:, j) - b_bar));
            bb = bb + 1 / N * (b(:, j)' * R * (r(:, j) - r_bar));
        end
        ss = bb / aa;
        s = ss;
        t = s * b_bar - R * r_bar;
        xs(i, :) = [q0, q1, q2, q3, t', s, ss];
    end
    
    Ls = zeros(x_count - 1, 1);
    for i = 1 : x_count - 1
        RR = quat2dcm_(xs(i, 1 : 4));
        tt = xs(i, 5 : 7)';
        Ls(i) = J_func(RR, tt, bbb, rrr, K);
    end
    [minimum, idx] = sort(Ls);
    R_ = quat2dcm_(xs(idx(1), 1 : 4));
    t_ = xs(idx(1), 5 : 7)';
    s_ = xs(idx(1), 9);
    min_val = minimum;
end

function J = J_func(R, t, b, r, K)
J = 0;
len = size(b, 2);
for j = 1 : len
    xb = b(j, :);
    bb = xb';
    xr = r(j, :)';
    xr = (R * xr + t);
    res = bb - xr;
    J = J + 1 / len * trace(res' * res);
end
end
