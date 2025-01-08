clc; clear; close all;
format short

addpath('gBnB_sampling/');
addpath('gBnB_stabbing/')
addpath('gBnB_opt/');
addpath('RANSAC+5pt/');
addpath('RANSAC+3pt/');

% Read image
img_path = '../datasets/kitti/sequences/00/image_0/';
images = dir(fullfile(img_path, "*.png"));

% Read ground truth pose
pose_path = '../datasets/kitti/poses/00.txt'; 
fileid_pose = fopen(pose_path);
poses = fscanf(fileid_pose, '%f');

% Reshape poses into a 3D matrix (3, 4, n)
poses_mat = zeros(3, 4, length(poses)/12);
for i = 1:length(poses) / 12
    poses_mat(:, :, i) = [
        poses((i-1)*12+1 : (i-1)*12+4)'; ...
        poses((i-1)*12+5 : (i-1)*12+8)'; ...
        poses((i-1)*12+9 : i*12)' ...
        ];
end

% Intrinsic matrix of camera 0 (gray)
K = [718.856, 0,           607.1928; ... 
     0,       718.856,     185.2157; ...
     0,       0,           1];

% Settings
img_pairs = 3;
img_gap = 11;
epsilon = 1e-3;
rho = 0.99; % requested probability of success
ransac_iter_num = 10000; % RANSAC methods iter times
ransac_repeat = 1;

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
pts_cell_1 = cell(img_pairs, 1);
pts_cell_2 = cell(img_pairs, 1);

% run experiments
for i = 1 : img_pairs

    fprintf("Processing image pair %d\n", i)

    idx_1 = i;
    idx_2 = i + img_gap;

    % Grond truth from KITTI poses file
    P1 = poses_mat(:, :, idx_1);
    R1 = P1(:, 1:3); t1 = P1(:, 4);

    P2 = poses_mat(:, :, idx_2);
    R2 = P2(:, 1:3); t2 = P2(:, 4);

    R_gt = R2' * R1;
    t_gt = R2' * (t1 - t2);

    axang = rotationMatrixToVector(R_gt);
    R_v = axang / norm(axang);
    theta_gt = norm(axang);
    t_gt = t_gt / norm(t_gt);

    Rv_mat = [0 -R_v(3) R_v(2); ...
              R_v(3) 0 -R_v(1); ...
             -R_v(2) R_v(1) 0]';

    % load images from sequence 00 in KITTI dataset
    img_1 = imread(fullfile(img_path, images(idx_1).name));
    img_2 = imread(fullfile(img_path, images(idx_2).name));

    % Detect features
    points_1 = detectSIFTFeatures(img_1);
    points_2 = detectSIFTFeatures(img_2);

    [f1, vpts1] = extractFeatures(img_1, points_1);
    [f2, vpts2] = extractFeatures(img_2, points_2);

    %indexPairs = matchFeatures(f1, f2) 
    indexPairs = matchFeatures(f1, f2, 'MatchThreshold', 20, 'MaxRatio', 0.7);
    matched_pts1 = vpts1(indexPairs(:, 1));
    matched_pts2 = vpts2(indexPairs(:, 2));

    pts1 = matched_pts1.Location;
    pts2 = matched_pts2.Location;

    pts1 = [pts1'; ones(1, size(pts1, 1))];
    pts2 = [pts2'; ones(1, size(pts2, 1))];

    pts1 = K \ pts1; pts1 = pts1 ./ vecnorm(pts1);
    pts2 = K \ pts2; pts2 = pts2 ./ vecnorm(pts2);

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

%     %-------------------------gBnB+sampling---------------------------%
%     tic
%     [t_gbnb_sampling, theta_gbnb_sampling] = solve_BnB_sampling(360, pts1, pts2, R_v, epsilon);
%     runtime_gbnb_sampling(i) = toc;
%     theta_err_gbnb_sampling(i) = abs(theta_gt - theta_gbnb_sampling) * 180 / pi;
%     t_err_gbnb_sampling(i) = real(acos(abs(t_gbnb_sampling' * t_gt))) * 180 / pi;
% 
%     %-------------------------gBnB+stabbing---------------------------%
%     tic
%     [t_gbnb_stabbing, theta_gbnb_stabbing] = solve_BnB_stabbing(pts1, pts2, R_v, epsilon);
%     runtime_gbnb_stabbing(i) = toc;
%     theta_err_gbnb_stabbing(i) = abs(theta_gt - theta_gbnb_stabbing) * 180 / pi;
%     t_err_gbnb_stabbing(i) = real(acos(abs(t_gbnb_stabbing' * t_gt))) * 180 / pi;

%     for j = 1:ransac_repeat
% 
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
%         theta_err_ransac_5pts(j, i) = norm(rotationMatrixToVector(R_ransac_5pts' * R_gt)) * 180 / pi;
%         t_err_ransac_5pts(j, i) = real(acos(abs(t_ransac_5pts' * t_gt))) * 180 / pi;
%     end
end

%% show results
theta_err_gbnb_opt_mean = mean(theta_err_gbnb_opt);
t_err_gbnb_opt_mean = mean(t_err_gbnb_opt);

theta_err_gbnb_sampling_mean = mean(theta_err_gbnb_sampling);
t_err_gbnb_sampling_mean = mean(t_err_gbnb_sampling);

theta_err_gbnb_stabbing_mean = mean(theta_err_gbnb_stabbing);
t_err_gbnb_stabbing_mean = mean(t_err_gbnb_stabbing);

theta_err_ransac_3pts_mean = mean(theta_err_ransac_3pts(1, :));
t_err_ransac_3pts_mean = mean(max(t_err_ransac_3pts, [], 1));

theta_err_ransac_5pts_mean = mean(theta_err_ransac_5pts(1, :));
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
fileID_1 = fopen("实验结果2/realworld/kitti/pts1_1e-3_40_SIFT_70_0.6_1-50.txt", 'w');
fileID_2 = fopen("实验结果2/realworld/kitti/pts2_1e-3_40_SIFT_70_0.6_1-50.txt", 'w');
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

   