clear all
close all
clc

format long g

warning('off');
addpath('APnP');
addpath('APnP/EPnP');

load zed-M-calib-left
cam_params_left = cam_params_zed_M_left;
corners = load('zed-M-left-tracked-chessboard-img-merged-corners.txt');
num_poses = size(corners, 1);
image_points = zeros(9, 2, num_poses);
valid_images = logical(zeros(num_poses, 1));
for i = 1 : size(corners, 1)
    if(corners(i, 2) ~= 0)
        image_points(:, :, i) = reshape(corners(i, 2 : 19), [2, 9])';
        valid_images(i) = true;
    else
        valid_images(i) = false;
    end
end
square_size = 10.5; % 10.5 mm
world_points = generateCheckerboardPoints([4, 4], square_size);

img_pos = load('zed-M-left-chessboard-pos.txt');

global image_pt world_pt K
K = cam_params_left.IntrinsicMatrix;

world_points_3d = [world_points, ones(length(world_points(:, 1)), 1)];

counter = 1;
x0 = 0;
is_draw = true;
is_draw_outlier = false;

image_points_ = image_points;
for i = 1 : num_poses
    if(valid_images(i))
        for j = 1 : size(image_points, 1)
            image_points_(j, :, counter) = image_points_(j, :, counter) + [img_pos(i, 2), img_pos(i, 3)];
        end
        counter = counter + 1;
    end
end


fileFolder = fullfile('./zed-M-left');
dirOutput = dir(fullfile(fileFolder, '*.jpg'));
fileNames = {dirOutput.name};
num_poses = length(fileNames);
for i = 1 : num_poses
    image_files_shown{i} = sprintf('./zed-M-left/%s', fileNames{i}); % store image files
end

poses = [];
counter = 1;
names = {};
loss = [];
use_p3p = true;
use_apnp = true;
use_epnp = true;
markersize = 20;
for i = 1 : num_poses
    if i == 210
        pause
    end
    if(valid_images(i))
        method_counter = 0;
        image_pt = [image_points_(:, 1, counter), image_points_(:, 2, counter)];
        world_pt = world_points_3d;
        world_pt(:, 1) = world_pt(:, 1) / 10.5 * 10;
        
        if(use_p3p)
            [worldOrientation, worldLocation, ~, ~, err_p3p] = ...
                estimateWorldCameraPose(image_pt, world_points_3d, cam_params_left);
            [e1, e2, e3] = dcm2angle_(worldOrientation);
            euler = [e1, e2, e3] * 180 / pi;
            J_pnp_p3p = J_pnp_loss_(image_pt, world_pt, K, worldOrientation, worldLocation');
            prj_p3p = generateProjectedPoints(world_pt, K, worldOrientation, worldLocation');
            method_counter = method_counter + 1;
            names{method_counter} = 'P3P';
            loss(method_counter) = J_pnp_p3p;
            Js_p3p(counter) = J_pnp_p3p;
        end
        
        if(use_apnp)
            [R_apnp, t_apnp] = apnp_opt(image_pt, world_pt, K);
            J_pnp_apnp = J_pnp_loss(image_pt, world_pt, K, R_apnp, t_apnp);
            [e1, e2, e3] = dcm2angle_(R_apnp');
            euler_apnp = [e1, e2, e3] * 180 / pi;
            prj_apnp = generateProjectedPoints_(world_pt, K, R_apnp, t_apnp);
            method_counter = method_counter + 1;
            names{method_counter} = 'APnP (Ours)';
            loss(method_counter) = J_pnp_apnp;
            Js_apnp(counter) = J_pnp_apnp;
        end
        
        if(use_epnp)
            [R_epnp, t_epnp] = epnp_opt(image_pt, world_pt, K);
            J_pnp_epnp = J_pnp_loss(image_pt, world_pt, K, R_epnp, t_epnp);
            [e1, e2, e3] = dcm2angle_(R_epnp');
            euler_epnp = [e1, e2, e3] * 180 / pi;
            prj_epnp = generateProjectedPoints_(world_pt, K, R_epnp, t_epnp);
            method_counter = method_counter + 1;
            names{method_counter} = 'EPnP';
            loss(method_counter) = J_pnp_epnp;
            Js_epnp(counter) = J_pnp_epnp;
        end
        
        
        
        table(names.', loss.')
        
        if(is_draw)
            imshow(imread(image_files_shown{i})); hold on
            %plot(image_points_(:, 1, counter), image_points_(:, 2, counter), 'Marker', '+', 'MarkerSize', markersize, 'LineStyle', 'none', 'Color', 'g', 'LineWidth', 8); hold on
            if(use_p3p)
                plot(prj_p3p(:, 1), prj_p3p(:, 2), 'Marker', 'd', 'MarkerSize', markersize, 'LineStyle', 'none', 'Color', 'r', 'LineWidth', 4); hold on
            end
            if(use_apnp)
                plot(prj_apnp(:, 1), prj_apnp(:, 2), 'Marker', 'o', 'MarkerSize', markersize, 'LineStyle', 'none', 'Color', 'g', 'LineWidth', 4); hold on
            end
            if(use_epnp)
                plot(prj_epnp(:, 1), prj_epnp(:, 2), 'Marker', 'x', 'MarkerSize', markersize, 'LineStyle', 'none', 'Color', 'y', 'LineWidth', 4); hold on
            end
            
            hold off
            legend(names, 'FontSize', 24, 'Interpreter', 'LaTeX', 'Color', [1, 1, 1] / 2);
            drawnow
        end

        counter = counter + 1;
    end
end

image_points_ = zeros(9, 2, 10);
counter = 1;
for i = 1 : size(image_points, 3)
    if(mod(i, 3) == 0)
        image_points_(:, :, counter) = image_points(:, :, i);
        counter = counter + 1;
    end
end

figure(3);
if(use_p3p)
    plot(Js_p3p, 'LineStyle', '--', 'LineWidth', 1.5); hold on
end
if(use_apnp)
    plot(Js_apnp, 'LineStyle', '-.', 'LineWidth', 1.5); hold on
end
hold off
legend({'P3P - Ransac', 'APnP (Ours)'}, 'FontSize', 15);
ylim([0 50]);
xlabel('Sequence Index');
ylabel('Loss Function Value');

figure(4);
if(use_epnp)
    plot(Js_epnp, 'LineStyle', '--', 'LineWidth', 1.5); hold on
end
if(use_apnp)
    plot(Js_apnp, 'LineStyle', '-.', 'LineWidth', 1.5); hold on
end
hold off
legend({'EPnP', 'APnP (Ours)'}, 'FontSize', 15);
ylim([0 50]);
xlabel('Sequence Index');
ylabel('Loss Function Value');



function J = J_pnp(x)
global image_pt world_pt K
R = quat2dcm_([x(1), x(2), x(3), x(4)]);
t = [x(5); x(6); x(7)];

J = J_pnp_loss(image_pt, world_pt, K, R, t);
end


function ProjectedPoints = generateProjectedPoints(world_pt, K, R, t)
numPoints = size(world_pt, 1);
ProjectedPoints = zeros(numPoints, 2);
R_ = R.';
t_ = - t.' * R';
for i = 1 : numPoints
    world_point = [world_pt(i, :), 1]; % homogeneous coordinates

    cameraMatrix = [R_; t_] * K;
    projectedPoint = world_point * cameraMatrix;
    projectedPoint = projectedPoint(1 : 2) ./ projectedPoint(3);
    ProjectedPoints(i, :) = projectedPoint;
end
end


function ProjectedPoints = generateProjectedPoints_(world_pt, K, R, t)
numPoints = size(world_pt, 1);
ProjectedPoints = zeros(numPoints, 2);
for i = 1 : numPoints
    world_point = [world_pt(i, :), 1]; % homogeneous coordinates

    cameraMatrix = [R; t'] * K;
    projectedPoint = world_point * cameraMatrix;
    projectedPoint = projectedPoint(1 : 2) ./ projectedPoint(3);
    ProjectedPoints(i, :) = projectedPoint;
end
end

function J = J_pnp_loss_(image_pt, world_pt, K, R, t)
numPoints = size(image_pt, 1);
R_ = R.';
t_ = - t.' * R';
J = 0;
for i = 1 : numPoints
    world_point = [world_pt(i, :), 1]; % homogeneous coordinates
    image_point = image_pt(i, :);

    cameraMatrix = [R_; t_] * K;
    projectedPoint = world_point * cameraMatrix;
    projectedPoint = projectedPoint(1 : 2) / projectedPoint(3);
    d = image_point - projectedPoint;
    J = J + d * d.';
end
J = J ./ numPoints;
end

function J = J_pnp_loss(image_pt, world_pt, K, R, t)
numPoints = size(image_pt, 1);
J = 0;
for i = 1 : numPoints
    world_point = [world_pt(i, :), 1]; % homogeneous coordinates
    image_point = image_pt(i, :);

    cameraMatrix = [R; t.'] * K;
    projectedPoint = world_point * cameraMatrix;
    
    if(~isnumeric(R))
        image_point = image_point .* projectedPoint(3);
        projectedPoint = projectedPoint(1 : 2);
    else
        projectedPoint = projectedPoint(1 : 2) / projectedPoint(3);
    end
    d = image_point - projectedPoint;
    J = J + d * d';
end
J = J ./ numPoints;
end

function [R, t] = apnp_opt(image_pt, world_pt, K)
    N = size(image_pt, 1);
    [R_init, t_init] = efficient_pnp_(...
        [world_pt, ones(N, 1)], [image_pt, ones(N, 1)], K', true);
    
    iter = 1;
    for i = 1 : iter
        x0 = [
            dcm2quat_(R_init').';
            t_init;
            1;
            ];
        opt = optimoptions('fminunc', 'Display', 'off', ...
            'Algorithm', 'quasi-newton', ...
            'UseParallel', false);
        problem.options = opt;
        problem.x0 = x0;
        problem.objective = @J_pnp;
        problem.solver = 'fminunc';
        [x1, f] = fminunc(problem);
        R = quat2dcm_(x1(1 : 4)');
        t = x1(5 : 7);
    end
end


function [R, t] = epnp_opt(image_pt, world_pt, K)
    N = size(image_pt, 1);
    [R_init, t_init] = efficient_pnp_(...
        [world_pt, ones(N, 1)], [image_pt, ones(N, 1)], K', false);
    
    iter = 1;
    for i = 1 : iter
        x0 = [
            dcm2quat_(R_init').';
            t_init;
            1;
            ];
        opt = optimoptions('fminunc', 'Display', 'off', ...
            'Algorithm', 'quasi-newton', ...
            'UseParallel', false);
        problem.options = opt;
        problem.x0 = x0;
        problem.objective = @J_pnp;
        problem.solver = 'fminunc';
        [x1, f] = fminunc(problem);
        R = quat2dcm_(x1(1 : 4)');
        t = x1(5 : 7);
    end
end
