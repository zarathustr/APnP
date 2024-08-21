function [R, t] = apnp_opt(image_pt_, world_pt_, K_, is_opt)
    global image_pt world_pt K
    image_pt = image_pt_;
    world_pt = world_pt_;
    K = K_;
    N = size(image_pt, 1);
    [R_init, t_init] = efficient_pnp_gauss_(...
        [world_pt, ones(N, 1)], [image_pt, ones(N, 1)], K', true);

    iter = 1;
    for i = 1 : iter
        if(~is_opt)
            R = R_init;
            t = t_init;
            break;
        end
        x0 = [
            vex(logm(R_init'));
            t_init;
            ];
        opt = optimoptions('fminunc', 'Display', 'off', ...
            'Algorithm', 'quasi-newton', ...
            'SpecifyObjectiveGradient', false, ...
            'FiniteDifferenceType', 'forward', ...
            'OptimalityTolerance', 1e-15, ...
            'FiniteDifferenceStepSize', 1e-7, ...
            'UseParallel', false);
        problem.options = opt;
        problem.x0 = x0;
        problem.objective = @J_pnp;
        problem.solver = 'fminunc';
        [x1, f] = fminunc(problem);
        R = expm(skew(x1(1 : 3)))';
        t = x1(4 : 6);
    end
end