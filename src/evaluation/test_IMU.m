clc; clear; close all;

addpath("./gBnB_opt/");
addpath("./gBnB_sampling");
addpath("./gBnB_stabbing");
addpath("./RANSAC+3pt/");

num_pts = 100; % Number of inliers & outliers
outlier_rate = 0.5; % outlier rate
noise_level = 1; % image noise   
bias_level_list = 0.1:0.1:0.5; % biased axis angle
iter_times = 200; % repeate 200 times
epsilon = 1e-3; % threshold
rho = 0.99; % requested probability of success

% initialization
len_imu_noise = length(bias_level_list);

theta_err_ransac_3pts = zeros(len_imu_noise, iter_times);
theta_err_gbnb_opt = zeros(len_imu_noise, iter_times);
theta_err_gbnb_sampling = zeros(len_imu_noise, iter_times);
theta_err_gbnb_stabbing = zeros(len_imu_noise, iter_times);

t_err_ransac_3pts = zeros(len_imu_noise, iter_times);
t_err_gbnb_opt = zeros(len_imu_noise, iter_times);
t_err_gbnb_sampling = zeros(len_imu_noise, iter_times);
t_err_gbnb_stabbing = zeros(len_imu_noise, iter_times);


for i = 1:length(bias_level_list)

    num_outlier = int16(num_pts * outlier_rate);
    num_inlier = num_pts - num_outlier;

    ransac_3pt_iter_num = max(10, ceil(log(1-rho)./log(1-(1-outlier_rate).^3))+1);

    bias_level = bias_level_list(i);
    
    for j = 1:iter_times

        [theta_gt, t_gt, x, y, R_vb] = gen_data_imu(num_inlier, num_outlier, noise_level, bias_level);

        %--------------------------RANSAC+3pts----------------------------%
        [theta_ransac_3pts, t_ransac_3pts, ~] = ransac_3pt(x, y, R_vb, epsilon, ransac_3pt_iter_num);
        theta_err_ransac_3pts(i, j) = abs(abs(theta_gt) - abs(theta_ransac_3pts)) * 180 / pi;
        t_err_ransac_3pts(i, j) = real(acos(abs(t_ransac_3pts' * t_gt))) * 180 / pi;

        %--------------------------gBnB+opt-------------------------------%
        [theta_gbnb_opt, t_gbnb_opt, ~] = GBnB(x, y, R_vb, epsilon);
        theta_err_gbnb_opt(i, j) = abs(abs(theta_gt) - abs(theta_gbnb_opt)) * 180 / pi;
        t_err_gbnb_opt(i, j) = real(acos(abs(t_gbnb_opt' * t_gt))) * 180 / pi;

        %-------------------------gBnB+sampling---------------------------%
        [t_gbnb_sampling, theta_gbnb_sampling] = solve_BnB_sampling(360, x, y, R_vb, epsilon);
        theta_err_gbnb_sampling(i, j) = abs(abs(theta_gt) - abs(theta_gbnb_sampling)) * 180 / pi;
        t_err_gbnb_sampling(i, j) = real(acos(abs(t_gbnb_sampling' * t_gt))) * 180 / pi;   

        %-------------------------gBnB+stabbing---------------------------%
        [t_gbnb_stabbing, theta_gbnb_stabbing] = solve_BnB_stabbing(x, y, R_vb, epsilon);
        theta_err_gbnb_stabbing(i, j) = abs(abs(theta_gt) - abs(theta_gbnb_stabbing)) * 180 / pi;
        t_err_gbnb_stabbing(i, j) = real(acos(abs(t_gbnb_stabbing' * t_gt))) * 180 / pi; 

    end
end
