function J = J_pnp(x)
global image_pt world_pt K

R = expm(skew(x(1 : 3)));
t = x(4 : 6);

J = J_pnp_loss(image_pt, world_pt, K, R, t);
end
