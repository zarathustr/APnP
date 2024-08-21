function [R, t] = apnp_analytic(world, image)
[R, t] = apnp_opt(image.', world.', eye(3), false);
end