% MAIN Illustrates how to use the EPnP algorithm described in:
%
%       Francesc Moreno-Noguer, Vincent Lepetit, Pascal Fua.
%       Accurate Non-Iterative O(n) Solution to the PnP Problem. 
%       In Proceedings of ICCV, 2007. 
%
% Copyright (C) <2007>  <Francesc Moreno-Noguer, Vincent Lepetit, Pascal Fua>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the version 3 of the GNU General Public License
% as published by the Free Software Foundation.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
% Francesc Moreno-Noguer, CVLab-EPFL, September 2007.
% fmorenoguer@gmail.com, http://cvlab.epfl.ch/~fmoreno/ 

clear all; close all;

addpath data;
addpath error;
addpath EPnP;
format long g

%1.-Generate simulated input data------------------------------------------
load_points = 0;
iter = 1000;
time_epnp = 0;
time_apnp = 0;
time_apnp_opt = 0;
time_epnp_gn = 0;
time_apnp_gn = 0;
time_spnp = 0;
err_epnp = 0;
err_apnp = 0;
err_apnp_opt = 0;
err_epnp_gn = 0;
err_apnp_gn = 0;
err_spnp = 0;
rot_err_epnp = 0;
rot_err_apnp = 0;
rot_err_apnp_opt = 0;
rot_err_epnp_gn = 0;
rot_err_apnp_gn = 0;
rot_err_spnp = 0;
trans_err_epnp = 0;
trans_err_apnp = 0;
trans_err_apnp_opt = 0;
trans_err_epnp_gn = 0;
trans_err_apnp_gn = 0;
trans_err_spnp = 0;

for kkk = 1 : iter
if ~load_points
    n = 10; %number of points
    std_noise = 5; %noise in the measurements (in pixels)
    [A, point, Rt] = generate_noisy_input_data(n, std_noise, 'donotplot');
    save('data\input_data_noise.mat', 'A','point', 'Rt');
else
    load('data\input_data_noise.mat', 'A', 'point', 'Rt');
    n = size(point, 2);
    draw_noisy_input_data(point);
end
Rt = inv(Rt);
R_true = Rt(1 : 3, 1 : 3);
t_true = Rt(1 : 3, 4);

%2.-Inputs format--------------------------------
x3d = zeros(n, 4);
x2d = zeros(n, 3); 
A = A(:, 1 : 3);
for i = 1 : n
    x3d_h(i, :) = [point(i).Xworld', 1]; 
    x2d_h(i, :) = [point(i).Ximg(1 : 2)', 1];

    %world and camera coordinates
    X3d_world(i, :) = point(i).Xworld';
    X3d_cam(i, :) = point(i).Xcam';
end


%3.-EPnP----------------------------------------------------
Xw = x3d_h(:, 1 : 3);
U = x2d_h(:, 1 : 2);

tic;
[Rp, Tp, Xc, sol] = efficient_pnp_(x3d_h, x2d_h, A, false);
time_epnp = time_epnp + 1 / iter * toc;
J_EPnP = J_func(Rp, Tp, Xc, Xw, A)

tic;
X3d_cam_ = X3d_cam;
X3d_cam_(:, 3) = ones(n, 1);
bb = (inv(A') * x2d_h');
[Rp_, Tp_, Xc_, sol_] = apnp_algebraic(x2d_h, X3d_world, 0, 0, A, zeros(3, 4));
time_apnp = time_apnp + 1 / iter * toc;
J_APnP = J_func(Rp_, Tp_, Xc, Xw, A)


[R_RPnP, t_RPnP] = RPnP(x3d_h(:, 1 : 3)', bb(1 : 2, :));
J_RPnP = J_func(R_RPnP, t_RPnP, Xc, Xw, A)

global image_pt world_pt K
image_pt = U;
world_pt = Xw;
K = A';

tic;
[Rp_apnp_opt, Tp_apnp_opt] = apnp_opt(image_pt, world_pt, A');
time_apnp_opt = time_apnp_opt + 1 / iter * toc;

% tic;
% [R_sym, t_sym] = symbolic_pnp(image_pt, world_pt, K);
% time_spnp = time_spnp + 1 / iter * toc;

%compute error
err_epnp = err_epnp + 1 / iter * reprojection_error_usingRT(0, Xw, U, Rp, Tp, A);
rot_err_epnp = rot_err_epnp + 1 / iter * acosd((trace(Rp * R_true') - 1) / 2);
trans_err_epnp = trans_err_epnp + 1 / iter * norm(Tp - t_true)^2;

err_apnp = err_apnp + 1 / iter * reprojection_error_usingRT(0, Xw, U, Rp_, Tp_, A);
rot_err_apnp = rot_err_apnp + 1 / iter * acosd((trace(Rp_ * R_true') - 1) / 2);
trans_err_apnp = trans_err_apnp + 1 / iter * norm(Tp_ - t_true)^2;

err_apnp_opt = err_apnp_opt + 1 / iter * reprojection_error_usingRT(0, Xw, U, Rp_apnp_opt, Tp_apnp_opt, A);
rot_err_apnp_opt = rot_err_apnp_opt + 1 / iter * acosd((trace(Rp_apnp_opt * R_true') - 1) / 2);
trans_err_apnp_opt = trans_err_apnp_opt + 1 / iter * norm(Tp_apnp_opt - t_true)^2;

% err_spnp = err_spnp + 1 / iter * reprojection_error_usingRT(0, Xw, U, R_sym, t_sym, A);
% rot_err_spnp = rot_err_spnp + 1 / iter * acosd((trace(R_sym * R_true') - 1) / 2);
% trans_err_spnp = trans_err_spnp + 1 / iter * norm(t_sym - t_true)^2;

%3.-EPnP_GAUSS_NEWTON----------------------------------------------------
Xw=x3d_h(:,1:3);
U=x2d_h(:,1:2);

tic;
[Rp, Tp, Xc, sol] = efficient_pnp_gauss(x3d_h, x2d_h, A, false);
time_epnp_gn = time_epnp_gn + 1 / iter * toc;
% J_EPnP_GN = J_func(Rp, Tp, Xc, Xw, A)

tic;
[Rp__, Tp__, xxs__] = efficient_pnp_gauss(x3d_h, x2d_h, A, true);
time_apnp_gn = time_apnp_gn + 1 / iter * toc;
% J_APnP_GN = J_func(Rp__, Tp__, Xc, Xw, A)

%compute error
err_epnp_gn = err_epnp_gn + 1 / iter * reprojection_error_usingRT(0, Xw, U, Rp, Tp, A);
rot_err_epnp_gn = rot_err_epnp_gn + 1 / iter * acosd((trace(Rp * R_true') - 1) / 2);
trans_err_epnp_gn = trans_err_epnp_gn + 1 / iter * norm(Tp - t_true)^2;

err_apnp_gn = err_apnp_gn + 1 / iter * reprojection_error_usingRT(0, Xw, U, Rp__, Tp__, A);
rot_err_apnp_gn = rot_err_apnp_gn + 1 / iter * acosd((trace(Rp__ * R_true') - 1) / 2);
trans_err_apnp_gn = trans_err_apnp_gn + 1 / iter * norm(Tp__ - t_true)^2;
end

trans_err_epnp = sqrt(trans_err_epnp);
trans_err_apnp = sqrt(trans_err_apnp);
trans_err_apnp_opt = sqrt(trans_err_apnp_opt);
trans_err_epnp_gn = sqrt(trans_err_epnp_gn);
trans_err_apnp_gn = sqrt(trans_err_apnp_gn);
% trans_err_spnp = sqrt(trans_err_spnp);

LastName = {'EPnP'; 'APnP'; 'APnP-Opt'; 'EPnP-GN'; 'APnP-GN'};
Time = [time_epnp; time_apnp; time_apnp_opt; time_epnp_gn; time_apnp_gn];
ReprojectionError = [err_epnp; err_apnp; err_apnp_opt; err_epnp_gn; err_apnp_gn];
RotationError = [rot_err_epnp; rot_err_apnp; rot_err_apnp_opt; rot_err_epnp_gn; rot_err_apnp_gn];
TranslationError = [trans_err_epnp; trans_err_apnp; trans_err_apnp_opt; trans_err_epnp_gn; trans_err_apnp_gn];
table(LastName, Time, ReprojectionError, RotationError, TranslationError)

% LastName = {'EPnP'; 'APnP'; 'APnP-Opt'; 'EPnP-GN'; 'APnP-GN'; 'SPnP'};
% Time = [time_epnp; time_apnp; time_apnp_opt; time_epnp_gn; time_apnp_gn; time_spnp];
% ReprojectionError = [err_epnp; err_apnp; err_apnp_opt; err_epnp_gn; err_apnp_gn; err_spnp];
% RotationError = [rot_err_epnp; rot_err_apnp; rot_err_apnp_opt; rot_err_epnp_gn; rot_err_apnp_gn; rot_err_spnp];
% TranslationError = [trans_err_epnp; trans_err_apnp; trans_err_apnp_opt; trans_err_epnp_gn; trans_err_apnp_gn; trans_err_spnp];
% table(LastName, Time, ReprojectionError, RotationError, TranslationError)

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








