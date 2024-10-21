function [theta_opt, t_opt] = polynomial_opt(pts1, pts2)
% Inputs: pts1, pts2: 3*n normalized matrices
% Outputs: theta_opt, t_opt
    
    pts1 = pts1 ./ pts1(3, :);
    pts2 = pts2 ./ pts2(3, :);

    vec = optimal_gravity_opt(pts1(1:2, :), pts2(1:2, :));
    y = vec(1);
    theta_opt = 2 * atan(y);
    t_opt = vec(2:4);
    t_opt = t_opt / norm(t_opt);
end