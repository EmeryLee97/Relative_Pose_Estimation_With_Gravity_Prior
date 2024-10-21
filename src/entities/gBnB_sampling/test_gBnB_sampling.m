clc; clear; close all;

format short

num_inlier = 50;
num_outlier = 200;
noise_level = 0.2;
epsilon = 3e-3;

[theta_gt, t_gt, pts_1, pts_2, R_v] = gen_data(num_inlier, num_outlier, noise_level);

%%

% dbstop in get_bound_theta_without_loop at 127 if upper_bound<50
tic
[t_opt, theta_opt] = solve_BnB_sampling(pts_1, pts_2, R_v, epsilon);
time1 = toc;

if t_gt(3) * t_opt(3) < 0
    t_error = pi - abs(acos(t_opt' * t_gt));
else
    t_error = abs(acos(t_opt' * t_gt));
end
theta_error = abs(theta_opt - theta_gt);

t_error = t_error * 180 / pi;
theta_error = theta_error * 180 /pi;

fprintf("Translation error: %f\n", t_error);
fprintf("Rotation error: %f\n", theta_error);
fprintf("Runtime: %f\n", time1);

