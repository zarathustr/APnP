function [R, t] = apnp_fminunc(world, image)
[R, t] = apnp_opt(image.', world.', eye(3), true);
end