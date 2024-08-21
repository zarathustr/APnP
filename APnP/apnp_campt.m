function Xc = apnp_campt(x3d_h, x2d_h, A)
Xw = x3d_h(:, 1 : 3);
U = x2d_h(:, 1 : 2);

Cw = [1, 0, 0;
      0, 1, 0;
      0, 0, 1;
      0, 0, 0];

n = size(Xw, 1);
C = [Cw, ones(4, 1)];
X = [Xw, ones(n, 1)];
Alph = X * inv(C);
M = compute_M_ver2(U, Alph, A);
Km = kernel_noise(M, 4);

X1 = Km(:, end);
[~, Xc] = compute_norm_sign_scaling_factor(X1, Cw, Alph, Xw);
end