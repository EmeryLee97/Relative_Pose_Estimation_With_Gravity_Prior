clc; clear; close all;
format short;

addpath('RANSAC+5pt/');
addpath('RANSAC+3pt/');
addpath('gBnB_opt/');
addpath('gBnB_sampling/');
addpath('gBnB_stabbing/');

% settings
img_pairs = 2; % number of image pairs to be tested
img_gap = 100; % e.g. test between image 0 and image 0+gap
epsilon = 1e-3; % tolerance
%rho = 0.99; % requested probability of success
ransac_iter_num = 10000; % RANSAC methods iter times
ransac_repeat = 50;

% Intrinsic matrix of camera fr3, images are undistored
K3 = [535.4, 0,     320.1; ...
      0,     539.2, 247.6; ...
      0,     0,     1];

% read data
dataset_path = "../datasets/tum_rgbd/rgbd_dataset_freiburg3_long_office_household/";
%dataset_path = "../datasets/tum_rgbd/rgbd_dataset_freiburg3_teddy/";
associated_path = strcat(dataset_path, "associate_with_groundtruth.txt");
fileid_pose = fopen(associated_path);
associated_data = textscan(fileid_pose, '%s %s %s %f %f %f %f %f %f %f');

img_path = associated_data{2};
tx = associated_data{4};
ty = associated_data{5};
tz = associated_data{6};
qx = associated_data{7};
qy = associated_data{8};
qz = associated_data{9};
qw = associated_data{10};

% initialization
runtime_gbnb_opt = zeros(1, img_pairs);
runtime_gbnb_sampling = zeros(1, img_pairs);
runtime_gbnb_stabbing = zeros(1, img_pairs);
runtime_ransac_3pts = zeros(ransac_repeat, img_pairs);
runtime_ransac_5pts = zeros(ransac_repeat, img_pairs);

theta_err_gbnb_opt = zeros(1, img_pairs);
theta_err_gbnb_sampling = zeros(1, img_pairs);
theta_err_gbnb_stabbing = zeros(1, img_pairs);
theta_err_ransac_3pts = zeros(ransac_repeat, img_pairs);
theta_err_ransac_5pts = zeros(ransac_repeat, img_pairs);

t_err_gbnb_opt = zeros(1, img_pairs);
t_err_gbnb_sampling = zeros(1, img_pairs);
t_err_gbnb_stabbing = zeros(1, img_pairs);
t_err_ransac_3pts = zeros(ransac_repeat, img_pairs);
t_err_ransac_5pts = zeros(ransac_repeat, img_pairs);

outlier_rates = zeros(1, img_pairs);
pts_nums = zeros(1, img_pairs);

% 2d cell, first dim contains index, second dim contains key points
pts_cell_1 = cell(100, 1);
pts_cell_2 = cell(100, 1);

for i = 1:img_pairs
    
    fprintf("Processing image pair %d\n", i)

    idx1 = i;
    idx2 = i + img_gap;
    % load images
    img1 = imread(strcat(dataset_path, img_path{idx1})); img1 = rgb2gray(img1);
    img2 = imread(strcat(dataset_path, img_path{idx2})); img2 = rgb2gray(img2);

    % Grond truth from KITTI poses file
    % t = [t1, t2, t3], q = [w, x, y, z]
    t1 = [tx(idx1); ty(idx1); tz(idx1)];
    t2 = [tx(idx2); ty(idx2); tz(idx2)];
    
    R1 = quat2rotm([qw(idx1) qx(idx1) qy(idx1) qz(idx1)]);
    R2 = quat2rotm([qw(idx2) qx(idx2) qy(idx2) qz(idx2)]);

    R_gt = R2' * R1;
    t_gt = R2' * (t1 - t2);
    t_gt = t_gt / norm(t_gt);

    axang = rotationMatrixToVector(R_gt);
    R_v = axang / norm(axang);
    theta_gt = norm(axang);

    Rv_mat = [0 -R_v(3) R_v(2); ...
              R_v(3) 0 -R_v(1); ...
             -R_v(2) R_v(1) 0]';

    % Detect features
    points_1 = detectSIFTFeatures(img1);
    points_2 = detectSIFTFeatures(img2);

    [f1, vpts1] = extractFeatures(img1, points_1);
    [f2, vpts2] = extractFeatures(img2, points_2);

    %indexPairs = matchFeatures(f1, f2);
    indexPairs = matchFeatures(f1, f2, 'MatchThreshold', 5, 'MaxRatio', 0.6);
    matched_pts1 = vpts1(indexPairs(:, 1));
    matched_pts2 = vpts2(indexPairs(:, 2));

    pts1 = matched_pts1.Location;
    pts2 = matched_pts2.Location;

    pts1 = [pts1'; ones(1, size(pts1, 1))];
    pts2 = [pts2'; ones(1, size(pts2, 1))];

    pts1 = K3 \ pts1; pts1 = pts1 ./ vecnorm(pts1);
    pts2 = K3 \ pts2; pts2 = pts2 ./ vecnorm(pts2);

    %-----------------------------
    if size(pts1, 2) < 5
        continue
    end

    pts_cell_1{i} = pts1;
    pts_cell_2{i} = pts2;

    % calculate the inlier rate
    b = cross(pts2, Rv_mat * pts1);
    c = -cross(pts2, Rv_mat^2 * pts1);
    a = cross(pts2, pts1) + cross(pts2, Rv_mat^2 * pts1);

    res = abs(t_gt' * (a + b*sin(theta_gt) + c*cos(theta_gt)));
    outlier_rate = sum(res > epsilon) / size(pts1, 2);
    outlier_rates(i) = outlier_rate;
    pts_nums(i) = size(pts1, 2);
    fprintf("Number of points: %d\n", size(pts1, 2));
    fprintf("Outlier rate: %d\n", outlier_rate);

    %--------------------------gBnB+opt-------------------------------%
    tic
    [theta_gbnb_opt, t_gbnb_opt, ~] = GBnB(pts1, pts2, R_v, epsilon);
    runtime_gbnb_opt(i) = toc;
    theta_err_gbnb_opt(i) = abs(theta_gt - theta_gbnb_opt) * 180 / pi;
    t_err_gbnb_opt(i) = real(acos(abs(t_gbnb_opt' * t_gt))) * 180 / pi;

    %-------------------------gBnB+sampling---------------------------%
    tic
    [t_gbnb_sampling, theta_gbnb_sampling] = solve_BnB_sampling(360, pts1, pts2, R_v, epsilon);
    runtime_gbnb_sampling(i) = toc;
    theta_err_gbnb_sampling(i) = abs(theta_gt - theta_gbnb_sampling) * 180 / pi;
    t_err_gbnb_sampling(i) = real(acos(abs(t_gbnb_sampling' * t_gt))) * 180 / pi;

    %-------------------------gBnB+stabbing---------------------------%
    tic
    [t_gbnb_stabbing, theta_gbnb_stabbing] = solve_BnB_stabbing(pts1, pts2, R_v, epsilon);
    runtime_gbnb_stabbing(i) = toc;
    theta_err_gbnb_stabbing(i) = abs(theta_gt - theta_gbnb_stabbing) * 180 / pi;
    t_err_gbnb_stabbing(i) = real(acos(abs(t_gbnb_stabbing' * t_gt))) * 180 / pi; 

%     for j = 1:ransac_repeat
%         %--------------------------RANSAC+3pts----------------------------%
%         tic
%         [theta_ransac_3pts, t_ransac_3pts, ~] = ransac_3pt(pts1, pts2, R_v, epsilon, ransac_iter_num);
%         runtime_ransac_3pts(j, i) = toc;
%         theta_err_ransac_3pts(j, i) = abs(theta_gt - theta_ransac_3pts) * 180 / pi;
%         t_err_ransac_3pts(j, i) = real(acos(abs(t_ransac_3pts' * t_gt))) * 180 / pi;
% 
%         %--------------------------RANSAC+5pts----------------------------%
%         tic
%         [R_ransac_5pts, t_ransac_5pts, ~] = ransac_5pt(pts1, pts2, R_v, epsilon, ransac_iter_num);
%         runtime_ransac_5pts(j, i) = toc;
%         %theta_ransac_5pts = acos((trace(R_ransac_5pts)-1) / 2);
%         %theta_err_ransac_5pts(i) = abs(theta_gt - theta_ransac_5pts) * 180 / pi;
%         theta_err_ransac_5pts(j, i) = norm(rotationMatrixToVector(R_ransac_5pts' * R_gt)) * 180 / pi;
%         t_err_ransac_5pts(j, i) = real(acos(abs(t_ransac_5pts' * t_gt))) * 180 / pi;
%     end
end

%%
theta_err_gbnb_opt_mean = mean(theta_err_gbnb_opt);
t_err_gbnb_opt_mean = mean(t_err_gbnb_opt);

theta_err_gbnb_sampling_mean = mean(theta_err_gbnb_sampling);
t_err_gbnb_sampling_mean = mean(t_err_gbnb_sampling);

theta_err_gbnb_stabbing_mean = mean(theta_err_gbnb_stabbing);
t_err_gbnb_stabbing_mean = mean(t_err_gbnb_stabbing);

theta_err_ransac_3pts_mean = mean(max(theta_err_ransac_3pts, [], 1));
t_err_ransac_3pts_mean = mean(max(t_err_ransac_3pts, [], 1));

theta_err_ransac_5pts_mean = mean(max(theta_err_ransac_5pts, [], 1));
t_err_ransac_5pts_mean = mean(max(t_err_ransac_5pts, [], 1));


fprintf("median_runtime_gbnb_opt: %d\n", median(runtime_gbnb_opt));
fprintf("median_runtime_gbnb_sampling: %d\n", median(runtime_gbnb_sampling));
fprintf("median_runtime_gbnb_stabbing: %d\n", median(runtime_gbnb_stabbing));
fprintf("median_runtime_ransac3pt: %d\n", median(runtime_ransac_3pts(:)));
fprintf("median_runtime_ransac5pt: %d\n", median(runtime_ransac_5pts(:)));
fprintf("\n");
fprintf("mean_runtime_gbnb_opt: %d\n", mean(runtime_gbnb_opt));
fprintf("mean_runtime_gbnb_sampling: %d\n", mean(runtime_gbnb_sampling));
fprintf("mean_runtime_gbnb_stabbing: %d\n", mean(runtime_gbnb_stabbing));
fprintf("mean_runtime_ransac3pt: %d\n", mean(runtime_ransac_3pts(:)));
fprintf("mean_runtime_ransac5pt: %d\n", mean(runtime_ransac_5pts(:)));
fprintf("\n");
fprintf("theta_err_gbnb_opt_mean: %d\n", theta_err_gbnb_opt_mean);
fprintf("theta_err_gbnb_sampling_mean: %d\n", theta_err_gbnb_sampling_mean);
fprintf("theta_err_gbnb_stabbing_mean: %d\n", theta_err_gbnb_stabbing_mean);
fprintf("theta_err_ransac_3pts_mean: %d\n", theta_err_ransac_3pts_mean);
fprintf("theta_err_ransac_5pts_mean: %d\n", theta_err_ransac_5pts_mean);
fprintf("\n");
fprintf("t_err_gbnb_opt_mean: %d\n", t_err_gbnb_opt_mean);
fprintf("t_err_gbnb_sampling_mean: %d\n", t_err_gbnb_sampling_mean);
fprintf("t_err_gbnb_stabbing_mean: %d\n", t_err_gbnb_stabbing_mean);
fprintf("t_err_ransac_3pts_mean: %d\n", t_err_ransac_3pts_mean);
fprintf("t_err_ransac_5pts_mean: %d\n", t_err_ransac_5pts_mean);

%% write key points to txt file
fileID_1 = fopen("实验结果2/realworld/tum/pts1_2e-3_120gap_sift_1-50.txt", 'w');
fileID_2 = fopen("实验结果2/realworld/tum/pts2_2e-3_120gap_sift_1-50.txt", 'w');
if fileID_1 == -1 || fileID_2 == -1
    error("无法打开文件")
end

for i = 1:size(pts_cell_1, 1)
    % each line: pts1(1, 0), pts1(2, 0), pts1(3, 0), ..., pts1(1, n), pts1(2, n), pts1(3, n)
    fprintf(fileID_1, '%f\t', reshape(pts_cell_1{i}, 1, []));
    fprintf(fileID_1, '\n');
    fprintf(fileID_2, '%f\t', reshape(pts_cell_2{i}, 1, []));
    fprintf(fileID_2, '\n');
end



runtime_box = [runtime_gbnb_opt', runtime_gbnb_sampling', runtime_gbnb_stabbing'];

r_err_box = [theta_err_gbnb_opt', theta_err_gbnb_sampling', theta_err_gbnb_stabbing', ...
             theta_err_ransac_3pts(1, :)', theta_err_ransac_5pts(1, :)'];

t_err_box = [t_err_gbnb_opt', t_err_gbnb_sampling', t_err_gbnb_stabbing', ...
             t_err_ransac_3pts(1, :)', t_err_ransac_5pts(1, :)'];

catagories_1 = {'GBnB', 'SaBnB', 'ISBnB'}; 
catagories_2 = {'GBnB', 'SaBnB', 'ISBnB', 'ransac3pt', 'ransac5pt'}; 

figure 
boxplot(runtime_box, 'Labels', catagories_1);
%ylim([1, 15]); 
ylabel("runtime (s)");

figure
boxplot(r_err_box, 'Labels', catagories_2);
%ylim([0, 5]); 
ylabel("Rotational Error (\circ)");

figure
boxplot(t_err_box, 'Labels', catagories_2);
%ylim([0, 15]); 
ylabel("Translational Error (\circ)");

