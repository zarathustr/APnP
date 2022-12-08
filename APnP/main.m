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

format long g

addpath data;
addpath error;
addpath EPnP;


fprintf('\n---------EPnP--------------\n');
%1.-Generate simulated input data------------------------------------------
load_points=0;
if ~load_points
    n = 8; %number of points
    std_noise = 1.5; %noise in the measurements (in pixels)
    [A,point,Rt]=generate_noisy_input_data(n, std_noise);
    save('data\input_data_noise.mat','A','point','Rt');
else
    load('data\input_data_noise.mat','A','point','Rt');
    n=size(point,2);
    draw_noisy_input_data(point);
end

%2.-Inputs format--------------------------------
x3d=zeros(n,3);
x2d=zeros(n,2); 
A=A(:,1:3);
for i=1:n
    x3d_h(i,:) = [point(i).Xworld',1]; 
    x2d_h(i,:) = [point(i).Ximg(1:2)',1];
    x3d(i,:) = point(i).Xworld';
    x2d(i,:) = point(i).Ximg(1:2)';

    %world and camera coordinates
    X3d_world(i,:)=point(i).Xworld';
    X3d_cam(i,:)=point(i).Xcam';
end


%3.-EPnP----------------------------------------------------
Xw=x3d_h(:,1:3);
U=x2d_h(:,1:2);

[Rp, Tp, Xc_epnp, sol] = efficient_pnp_(x3d_h, x2d_h, A, false);
J_EPnP = J_func(Rp, Tp, Xc_epnp, Xw, A)


% [Rp_, Tp_, Xc, sol] = efficient_pnp_guass(x3d_h, x2d_h, A, true);
[Rp_, Tp_, Xc_apnp, sol] = efficient_pnp_(x3d_h, x2d_h, A, true);

% Xc_apnp = apnp_campt(x3d_h, x2d_h, A);
% [Rp_, Tp_, s_, xxs, min_val] = apnp_algebraic(Xc_epnp, Xw, Xw, U, A);
% s_
J_APnP = J_func(Rp_, Tp_, Xc_apnp, Xw, A)

[R_, t_, X, xs, min_val] = symbolic_pnp(x2d, x3d, A.');
J_symbolic = J_func(R_, t_, Xc_apnp, Xw, A)




%draw Results
for i=1:n
    point(i).Xcam_est=Xc_epnp(i,:)';
end
figure; h=gcf;
plot_3d_reconstruction(point,'EPnP (Old)',h);
xlim([-2 2]); ylim([-2 2]);

%compute error
error=reprojection_error_usingRT_(0, Xw, U, Rp, Tp, A);
fprintf('error EPnP: %.3f\n', error);

error=reprojection_error_usingRT_(0, Xw, U, Rp_, Tp_, A);
fprintf('error APnP: %.3f\n', error);


%3.-EPnP_GAUSS_NEWTON----------------------------------------------------
Xw=x3d_h(:,1:3);
U=x2d_h(:,1:2);

[Rp,Tp,Xc,sol]=efficient_pnp_gauss_(x3d_h,x2d_h,A,false);
J_EPnP_GN = J_func(Rp, Tp, Xc, Xw, A)

[Rp__, Tp__, xxs__] = efficient_pnp_gauss_(x3d_h,x2d_h,A,true);
J_APnP_GN = J_func(Rp__, Tp__, Xc, Xw, A)

%draw Results
for i=1:n
    point(i).Xcam_est=Xc(i,:)';
end
figure; h=gcf;
plot_3d_reconstruction(point,'EPnP Gauss Newton',h);

%compute error
error=reprojection_error_usingRT_(0, Xw,U,Rp,Tp,A);
fprintf('error EPnP_Gauss_Newton: %.3f\n',error);
xlim([-2 2]); ylim([-2 2]);

error=reprojection_error_usingRT_(0, Xw,U,Rp__,Tp__,A);
fprintf('error APnP_Gauss_Newton: %.3f\n',error);
xlim([-2 2]); ylim([-2 2]);


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


