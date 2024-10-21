clc; clear; close all

addpath("./gBnB_opt/");
addpath("./gBnB_sampling/");
addpath("./gBnB_stabbing/")
addpath("./RANSAC+3pt/");
addpath("./RANSAC+5pt/");
addpath("./polynomial_opt/");


num_pts_list = 100:100:1000; % Number of inliers & outliers
outlier_rate = 0.5; % different outlier rates
noise_level= 1; % noise level
iter_times = 100; % repeate 100 times
epsilon = 1e-3; % threshold
rho = 0.99; % requested probability of success
method_list = ["ransac_5pts", ...
               "ransac_3pts", ...
               "polynomial_opt", ...
               "gbnb_opt", ...
               "gbnb_sampling", ...
               "gbnb_stabbing"];

% initialization
len_pts = length(num_pts_list);

runtime_ransac_3pts = zeros(len_pts, iter_times);
runtime_ransac_5pts = zeros(len_pts, iter_times);
runtime_polynomial_opt = zeros(len_pts, iter_times);
runtime_gbnb_opt = zeros(len_pts, iter_times);
runtime_gbnb_sampling = zeros(len_pts, iter_times);
runtime_gbnb_stabbing = zeros(len_pts, iter_times);

theta_err_ransac_3pts = zeros(len_pts, iter_times);
theta_err_ransac_5pts = zeros(len_pts, iter_times);
theta_err_polynomial_opt = zeros(len_pts, iter_times);
theta_err_gbnb_opt = zeros(len_pts, iter_times);
theta_err_gbnb_sampling = zeros(len_pts, iter_times);
theta_err_gbnb_stabbing = zeros(len_pts, iter_times);

t_err_ransac_3pts = zeros(len_pts, iter_times);
t_err_ransac_5pts = zeros(len_pts, iter_times);
t_err_polynomial_opt = zeros(len_pts, iter_times);
t_err_gbnb_opt = zeros(len_pts, iter_times);
t_err_gbnb_sampling = zeros(len_pts, iter_times);
t_err_gbnb_stabbing = zeros(len_pts, iter_times);

for i = 1:length(num_pts_list)

    num_pts = num_pts_list(i);
    num_outlier = int16(num_pts * outlier_rate);
    num_inlier = num_pts - num_outlier;

    ransac_3pt_iter_num = max(10, ceil(log(1-rho)./log(1-(1-outlier_rate).^3))+1);
    ransac_5pt_iter_num = max(10, ceil(log(1-rho)./log(1-(1-outlier_rate).^5))+1);

    for j = 1:iter_times

        [theta_gt, t_gt, x, y, R_v] = gen_data(num_inlier, num_outlier, noise_level);

        R_gt = rotationVectorToMatrix(R_v * theta_gt);

        %--------------------------RANSAC+3pts----------------------------%
        tic
        [theta_ransac_3pts, t_ransac_3pts, ~] = ransac_3pt(x, y, R_v, epsilon, ransac_3pt_iter_num);
        runtime_ransac_3pts(i, j) = toc;
        theta_err_ransac_3pts(i, j) = abs(abs(theta_gt) - abs(theta_ransac_3pts)) * 180 / pi;
        t_err_ransac_3pts(i, j) = real(acos(abs(t_ransac_3pts' * t_gt))) * 180 / pi;

        %--------------------------RANSAC+5pts----------------------------%
        tic
        [R_ransac_5pts, t_ransac_5pts, ~] = ransac_5pt(x, y, R_v, epsilon, ransac_5pt_iter_num);
        runtime_ransac_5pts(i, j) = toc;
        theta_ransac_5pts = acos((trace(R_ransac_5pts)-1) / 2);
        % theta_err_ransac_5pts(i, j) = norm(rotationMatrixToVector(R_ransac_5pts' * R_gt)) * 180 / pi;
        theta_err_ransac_5pts(i, j) = abs(abs(theta_gt) - abs(theta_ransac_5pts)) * 180 / pi;
        t_err_ransac_5pts(i, j) = real(acos(abs(t_ransac_5pts' * t_gt))) * 180 / pi;

        %-------------------------polynomial+opt--------------------------%
        tic
        [theta_polynomial_opt, t_polynomial_opt] = polynomial_opt(x, y);
        runtime_polynomial_opt(i, j) = toc;
        theta_err_polynomial_opt(i, j) = abs(abs(theta_gt) - abs(theta_polynomial_opt)) * 180 / pi;
        t_err_polynomial_opt(i, j) = real(acos(abs(t_polynomial_opt' * t_gt))) * 180 / pi;

        %--------------------------gBnB+opt-------------------------------%
        tic
        [theta_gbnb_opt, t_gbnb_opt, ~] = GBnB(x, y, R_v, epsilon);
        runtime_gbnb_opt(i, j) = toc;
        theta_err_gbnb_opt(i, j) = abs(abs(theta_gt) - abs(theta_gbnb_opt)) * 180 / pi;
        t_err_gbnb_opt(i, j) = real(acos(abs(t_gbnb_opt' * t_gt))) * 180 / pi;

        %-------------------------gBnB+sampling---------------------------%
        tic
        [t_gbnb_sampling, theta_gbnb_sampling] = solve_BnB_sampling(360, x, y, R_v, epsilon);
        runtime_gbnb_sampling(i, j) = toc;
        theta_err_gbnb_sampling(i, j) = abs(abs(theta_gt) - abs(theta_gbnb_sampling)) * 180 / pi;
        t_err_gbnb_sampling(i, j) = real(acos(abs(t_gbnb_sampling' * t_gt))) * 180 / pi;  

        %-------------------------gBnB+stabbing---------------------------%
        tic
        [t_gbnb_stabbing, theta_gbnb_stabbing] = solve_BnB_stabbing(x, y, R_v, epsilon);
        runtime_gbnb_stabbing(i, j) = toc;
        theta_err_gbnb_stabbing(i, j) = abs(abs(theta_gt) - abs(theta_gbnb_stabbing)) * 180 / pi;
        t_err_gbnb_stabbing(i, j) = real(acos(abs(t_gbnb_stabbing' * t_gt))) * 180 / pi;   

    end
end

%% plots

runtime_ransac_5pts = median(runtime_ransac_5pts, 2);
runtime_ransac_3pts = median(runtime_ransac_3pts, 2);
runtime_polynomial_opt = median(runtime_polynomial_opt, 2);
runtime_gbnb_opt = median(runtime_gbnb_opt, 2);
runtime_gbnb_sampling = median(runtime_gbnb_sampling, 2);
runtime_gbnb_stabbing = median(runtime_gbnb_stabbing, 2);

color1 = [143, 185, 67] / 255; % 8fb943
color2 = [120, 185, 210] / 255; % 78b9d2
color3 = [245, 207, 54] / 255; % f5cf36
color4 = [131, 134, 168] / 255; % 8386a8
color5 = [209, 92, 107] / 255; % d15c6b
color6 = [152, 23, 213] / 255;

color = {color1; color2; color3; color4; color5; color6};

slg = semilogy(num_pts_list, runtime_gbnb_opt, '-diamond', ...
         num_pts_list, runtime_gbnb_sampling, '-o', ...
         num_pts_list, runtime_gbnb_stabbing, '-pentagram', ...
         num_pts_list, runtime_polynomial_opt, '-v', ...
         num_pts_list, runtime_ransac_3pts, '-^', ...
         num_pts_list, runtime_ransac_5pts, '-s');

slg(1).LineWidth = 2; slg(1).Color = color1; slg(1).MarkerSize = 10;
slg(2).LineWidth = 2; slg(2).Color = color2; slg(2).MarkerSize = 10;
slg(3).LineWidth = 2; slg(3).Color = color3; slg(3).MarkerSize = 10;
slg(4).LineWidth = 2; slg(4).Color = color4; slg(4).MarkerSize = 10;
slg(5).LineWidth = 2; slg(5).Color = color5; slg(5).MarkerSize = 10;
slg(6).LineWidth = 2; slg(6).Color = color6; slg(6).MarkerSize = 10;

legend('gbnb\_opt', 'gbnb\_sampling', 'gbnb\_stabbing', ...
    'polynomial\_opt', 'ransac\_3pts', 'ransac\_5pts');

xlabel("Points");
ylabel("Runtime (s)");

grid on;

% Sort according to this order:
% 1st: ransac_5pts
% 2nd: ransac_3pts
% 3rd: ransac_opt
% 4th: gbnb_opt
% 5th: gbnb_sampling
% 6th: gbnb_stabbing

% table of size 250 * [outlier_list, trans_error, rot_error, method]
test_points_result = table;

factor = [100*ones(iter_times, 1); ...
    200*ones(iter_times, 1); ...
    300*ones(iter_times, 1); ...
    400*ones(iter_times, 1); ...
    500*ones(iter_times, 1); ...
    600*ones(iter_times, 1); ...
    700*ones(iter_times, 1); ...
    800*ones(iter_times, 1); ...
    900*ones(iter_times, 1); ...
    1000*ones(iter_times, 1)];

factor_vec = unique(factor);

test_points_result.points_list = [factor; factor; factor; factor; factor; factor];

test_points_result.trans_error = [reshape(t_err_ransac_5pts', [], 1); ...
                             reshape(t_err_ransac_3pts', [], 1); ...
                             reshape(t_err_polynomial_opt', [], 1); ...
                             reshape(t_err_gbnb_opt', [], 1); ...
                             reshape(t_err_gbnb_sampling', [], 1); ...
                             reshape(t_err_gbnb_stabbing', [], 1)];

test_points_result.rot_error = [reshape(theta_err_ransac_5pts', [], 1); ...
                               reshape(theta_err_ransac_3pts', [], 1); ...
                               reshape(theta_err_polynomial_opt', [], 1); ...
                               reshape(theta_err_gbnb_opt', [], 1); ...
                               reshape(theta_err_gbnb_sampling', [], 1); ...
                               reshape(theta_err_gbnb_stabbing', [], 1)];

test_points_result.method = [string(repmat("Ransac\_5pts", iter_times*length(num_pts_list), 1)); ...
                            string(repmat("Ransac\_3pts", iter_times*length(num_pts_list), 1)); ...
                            string(repmat("Polynomial\_opt", iter_times*length(num_pts_list), 1)); ...
                            string(repmat("GBnB\_opt", iter_times*length(num_pts_list), 1)); ...
                            string(repmat("GBnB\_sampling", iter_times*length(num_pts_list), 1)); ...
                            string(repmat("GBnB\_stabbing", iter_times*length(num_pts_list), 1))];

offset = 15 * ones(size(factor));
offset_list = [-2.5*offset; -1.5*offset; -0.5*offset; 0.5*offset; 1.5*offset; 2*offset];
test_points_result.points_list = test_points_result.points_list + offset_list;

% draw translation error

figure
axesbox = boxchart(test_points_result.points_list, ...
    test_points_result.trans_error, ...
    'GroupByColor', test_points_result.method,...
    'WhiskerLineStyle', '--', ...
    'BoxFaceAlpha', 0.6, ...
    'MarkerStyle', 'x', ...
    'MarkerSize', 8, ...
    'BoxWidth', 50);

axesbox(1).WhiskerLineColor = color1;
axesbox(2).WhiskerLineColor = color2;
axesbox(3).WhiskerLineColor = color3;
axesbox(4).WhiskerLineColor = color4;
axesbox(5).WhiskerLineColor = color5;
axesbox(6).WhiskerLineColor = color6;

set(axesbox, {'BoxFaceColor'}, color, {'MarkerColor'}, color, 'LineWidth', 1.5)
set(gca, 'XTick', factor_vec, 'XTickLabel', ...
    ["100", "200", "300", "400", "500", "600", "700", "800", "900", "1000"]);
set(gca,'XGrid','on','YGrid','off', 'LineWidth', 1, 'Fontsize', 11)
set(gca,'yscale','log')
set(gca,'box','on')

xlabel("Points");
ylabel("Translation error (\circ)");
legend();

% draw rotation error

figure
axesbox = boxchart(test_points_result.points_list, ...
    test_points_result.rot_error, ...
    'GroupByColor', test_points_result.method,...
    'WhiskerLineStyle', '--', ...
    'BoxFaceAlpha', 0.6, ...
    'MarkerStyle', 'x', ...
    'MarkerSize', 8, ...
    'BoxWidth', 50);

axesbox(1).WhiskerLineColor = color1;
axesbox(2).WhiskerLineColor = color2;
axesbox(3).WhiskerLineColor = color3;
axesbox(4).WhiskerLineColor = color4;
axesbox(5).WhiskerLineColor = color5;
axesbox(6).WhiskerLineColor = color6;

set(axesbox, {'BoxFaceColor'}, color, {'MarkerColor'}, color, 'LineWidth', 1.5)
set(gca, 'XTick', factor_vec, 'XTickLabel', ...
    ["100", "200", "300", "400", "500", "600", "700", "800", "900", "1000"]);
set(gca,'XGrid','on','YGrid','off', 'LineWidth', 1, 'Fontsize', 11)
set(gca,'yscale','log')
set(gca,'box','on')

xlabel("Points");
ylabel("Rotation error (\circ)");
legend();
