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
%     image_point = image_point;
    projectedPoint = projectedPoint(1 : 2) / projectedPoint(3);
%     projectedPoint = projectedPoint(1 : 2) ./ projectedPoint(3);
    d = image_point - projectedPoint;
    J = J + d * d.';
end
J = J ./ numPoints;
end