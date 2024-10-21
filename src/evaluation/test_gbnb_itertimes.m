clc; clear; close all

addpath("./gBnB_opt/");
addpath("./gBnB_sampling");
addpath("./gBnB_stabbing");


num_pts = 100; % Number of inliers & outliers
outlier_rate_list = 0.1:0.1:0.8; % different outlier rates
noise_level= 1; % noise level
iter_times = 100; % repeate 100 times
epsilon = 1e-3; % threshold
method_list = ["gbnb_opt", "gbnb_sampling", "gbnb_stabbing"];

gbnb_opt_iter = zeros(length(outlier_rate_list), iter_times);
gbnb_sampling_iter = zeros(length(outlier_rate_list), iter_times);
gbnb_stabbing_iter = zeros(length(outlier_rate_list), iter_times);

for i = 1:length(outlier_rate_list)

    outlier_rate = outlier_rate_list(i);
    num_outlier = int16(num_pts * outlier_rate);
    num_inlier = num_pts - num_outlier;

    for j = 1:iter_times

        [theta_gt, t_gt, x, y, R_v] = gen_data(num_inlier, num_outlier, noise_level);

        R_gt = rotationVectorToMatrix(R_v * theta_gt);

        %--------------------------gBnB+opt-------------------------------%
        [~, ~, gbnb_opt_iter(i, j)] = GBnB_itertimes(x, y, R_v, epsilon);

        %-------------------------gBnB+sampling---------------------------%
        [~, ~, gbnb_sampling_iter(i, j)] = solve_BnB_sampling_itertimes(360, x, y, R_v, epsilon);

        %-------------------------gBnB+stabbing---------------------------%
        [~, ~, gbnb_stabbing_iter(i, j)] = solve_BnB_stabbing_itertimes(x, y, R_v, epsilon);

    end
end

%% plots

gbnb_opt_iter_mean = mean(gbnb_opt_iter, 2);
gbnb_sampling_iter_mean = mean(gbnb_sampling_iter, 2);
gbnb_stabbing_iter_mean = mean(gbnb_stabbing_iter, 2);

gbnb_opt_iter_median = median(gbnb_opt_iter, 2);
gbnb_sampling_iter_median = median(gbnb_sampling_iter, 2);
gbnb_stabbing_iter_median = median(gbnb_stabbing_iter, 2);

color1 = [41, 155, 207] / 255;
color2 = [115, 180, 77] / 255;
color3 = [216, 33, 28] / 255;

% semilogy(outlier_rate_list, runtime_ransac_5pts, ...
%          outlier_rate_list, runtime_ransac_3pts, ...
%          outlier_rate_list, runtime_polynomial_opt, ...
%          outlier_rate_list, runtime_gbnb_opt, ...
%          outlier_rate_list, runtime_gbnb_sampling);

figure;

l = plot(outlier_rate_list, gbnb_opt_iter_mean, '-diamond', ...
    outlier_rate_list, gbnb_opt_iter_median, '--diamond', ...
    outlier_rate_list, gbnb_sampling_iter_mean, '-o', ...
    outlier_rate_list, gbnb_sampling_iter_median, '--o', ...
    outlier_rate_list, gbnb_stabbing_iter_mean, '-pentagram', ...
    outlier_rate_list, gbnb_stabbing_iter_median, '--pentagram');

l(1).Color = color1;
l(2).Color = color1;
l(3).Color = color2;
l(4).Color = color2;
l(5).Color = color3;
l(6).Color = color3;


legend('gbnb\_opt\_mean', 'gbnb\_opt\_median', ...
       'gbnb\_sampling\_mean', 'gbnb\_sampling\_median', ...
       'gbnb\_stabbing\_mean', 'gbnb\_stabbing\_median');

xlabel("Outlier Rate");
ylabel("iter times");



        


