clc; clear; close all;

num_all = 100;
outlier_rate = 0.9;
num_inlier = int16(num_all * (1-outlier_rate));
num_outlier = int16(num_all * outlier_rate);
noise_level = 2;
epsilon = 1e-3;

[theta_gt, t_gt, pts1, pts2, R_v] = gen_data(num_inlier, num_outlier, noise_level);

[theta_opt, t_opt] = optimal_gravity_opt(pts1, pts2);

if t_opt(1) * t_gt(1) >= 0
    t_error = acos(t_opt' * t_gt)*180/pi;
else
    t_error = 180 - acos(t_opt' * t_gt)*180/pi;
end
theta_error = abs(theta_opt - theta_gt)*180/pi;

disp(t_error);
disp(theta_error);