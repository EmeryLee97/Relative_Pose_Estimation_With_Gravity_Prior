clc; clear; close all;

addpath("../gBnB_opt/");
addpath("../gBnB_sampling");
addpath("../polynomial_opt");

num_inlier = 50;
num_outlier = 50;
noise_level = 1;
epsilon = 1e-3;

iter_time = 1;

[theta_gt, t_gt, pts_1, pts_2, R_v] = gen_data(num_inlier, num_outlier, noise_level);

%%
tic
[theta_opt, t_opt, ~] = GBnB(pts_1, pts_2, R_v, epsilon);
toc;

theta_error_opt = abs(abs(theta_gt) - abs(theta_opt)) * 180 / pi;
t_error_opt = real(acos(abs(t_opt' * t_gt))) * 180 / pi;
fprintf("Translation error: %f degree\n", t_error_opt);
fprintf("Rotation error: %f degree\n", theta_error_opt);
% 
% 
% tic
% [t_sampling, theta_sampling] = solve_BnB_sampling(360, pts_1, pts_2, R_v, epsilon);
% toc
% 
% theta_error_sampling = abs(abs(theta_gt) - abs(theta_sampling)) * 180 / pi;
% t_error_sampling = real(acos(abs(t_sampling' * t_gt))) * 180 / pi;
% fprintf("Translation error: %f degree\n", t_error_sampling);
% fprintf("Rotation error: %f degree\n", theta_error_sampling);
% 
% 
tic
[t_stabbing, theta_stabbing] = solve_BnB_stabbing(pts_1, pts_2, R_v, epsilon);
toc;

theta_error_stabbing = abs(abs(theta_gt) - abs(theta_stabbing)) * 180 / pi;
t_error_stabbing = real(acos(abs(t_stabbing' * t_gt))) * 180 / pi;
fprintf("Translation error: %f degree\n", t_error_stabbing);
fprintf("Rotation error: %f degree\n", theta_error_stabbing);


