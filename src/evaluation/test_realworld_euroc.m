clc; clear; close all;
format long

addpath('gBnB_sampling/');
addpath('gBnB_stabbing/')
addpath('gBnB_opt/');
addpath('RANSAC+5pt/');
addpath('RANSAC+3pt/');

% settings
img_pairs = 1; % number of image pairs to be tested
img_gap = 30; % e.g. test between image 0 and image 0+gap
epsilon = 2e-3; % tolerance
%rho = 0.99; % requested probability of success
ransac_iter_num = 10000; % RANSAC methods iter times
ransac_repeat = 50;

% Read image
dataset_path = '/Users/xuhuili/datasets/euroc/V1_01_easy/mav0/';
img_paths = dir(fullfile(strcat(dataset_path, 'cam0/data/'), "*.png"));
img_names = arrayfun(@(x) x.name, img_paths, 'UniformOutput', false);
img_timestamps = cellfun((@(str) str2double(strrep(str, '.png', ''))), {img_paths.name});
img_timestamps = img_timestamps / 1e9; % from ns to seconds

%img_path = strcat(img_paths(1).folder, '/', img_paths(1).name);

% Read ground truth pose (timestamp, x, y, z, qw, qx, qz)
pose_path = strcat(dataset_path, 'state_groundtruth_estimate0/V1_01_easy.txt');
fileid_pose = fopen(pose_path, 'r');
poses_raw = textscan(fileid_pose, '%f %f %f %f %f %f %f %f', 'CommentStyle', '#');
pose_timestamps = poses_raw{1} / 1e9;

% associate img_timestamps and pose_timestamps
associations = [];
max_dt = 0.0001; % seconds
for i = 1:length(img_paths)
    img_tstamp = img_timestamps(i);
    [min_dist, idx] = min(abs(pose_timestamps - img_tstamp));
    if min_dist < max_dt
        associations = [associations; [i, idx]];
    end
end

% camera config
K = [458.654, 0,       367.215;
     0        457.296, 248.375;
     0,       0,       1];

radial_distortion = [-0.28340811, 0.07395907];
tangential_distortion = [0.00019359, 1.76187114e-05];

% K need to be transposes
camera_params = cameraParameters(...
    'IntrinsicMatrix', K', ...
    'RadialDistortion', radial_distortion, ...
    'TangentialDistortion', tangential_distortion);


for i = 200:200
    fprintf("Processing image pair %d and %d\n", i, i+img_gap)

    % read images
    image_1_distorted = imread(strcat(dataset_path, 'cam0/data/', img_names{associations(i, 1)}));
    image_1 = undistortImage(image_1_distorted, camera_params);
    image_2_distorted = imread(strcat(dataset_path, 'cam0/data/', img_names{associations(i+img_gap, 1)}));
    image_2 = undistortImage(image_2_distorted, camera_params);

    figure; 
    subplot(1, 2, 1); imshow(image_1_distorted); title('distorted');
    subplot(1, 2, 2); imshow(image_1); title('undistorted')

    figure; 
    subplot(1, 2, 1); imshow(image_2_distorted); title('distorted');
    subplot(1, 2, 2); imshow(image_2); title('undistorted')

    % read ground truth poses and calculate relative pose
    pose_idx_1 = associations(i, 2);
    pose_idx_2 = associations(i+img_gap, 2);

    trans_gt_1 = [poses_raw{2}(pose_idx_1); poses_raw{3}(pose_idx_1); poses_raw{4}(pose_idx_1)];
    trans_gt_2 = [poses_raw{2}(pose_idx_2); poses_raw{3}(pose_idx_2); poses_raw{4}(pose_idx_2)];

    rot_gt_1 = quat2rotm([poses_raw{5}(pose_idx_1), poses_raw{6}(pose_idx_1), ...
                        poses_raw{7}(pose_idx_1), poses_raw{8}(pose_idx_1)]);
    rot_gt_2 = quat2rotm([poses_raw{5}(pose_idx_2), poses_raw{6}(pose_idx_2), ...
                        poses_raw{7}(pose_idx_2), poses_raw{8}(pose_idx_2)]);

    R_gt = rot_gt_2' * rot_gt_1;
    t_gt = rot_gt_2' * (trans_gt_1 - trans_gt_2);
    t_gt = t_gt / norm(t_gt);

    axang = rotationMatrixToVector(R_gt);
    R_v = axang / norm(axang);
    theta_gt = norm(axang);

    Rv_mat = [0,      -R_v(3), R_v(2); ...
              R_v(3),  0,      -R_v(1); ...
             -R_v(2),  R_v(1),  0]';

    % Detect features
    points_1 = detectSIFTFeatures(image_1);
    points_2 = detectSIFTFeatures(image_2);

    [f1, vpts1] = extractFeatures(image_1, points_1);
    [f2, vpts2] = extractFeatures(image_2, points_2);

    indexPairs = matchFeatures(f1, f2);
    %indexPairs = matchFeatures(f1, f2, 'MatchThreshold', 10, 'MaxRatio', 0.6);
    matched_pts1 = vpts1(indexPairs(:, 1));
    matched_pts2 = vpts2(indexPairs(:, 2));

    pts1 = matched_pts1.Location;
    pts2 = matched_pts2.Location;

    pts1 = [pts1'; ones(1, size(pts1, 1))];
    pts2 = [pts2'; ones(1, size(pts2, 1))];

    pts1 = K \ pts1; pts1 = pts1 ./ vecnorm(pts1);
    pts2 = K \ pts2; pts2 = pts2 ./ vecnorm(pts2);

    %--------------------debug-------------------
    t_mat = [0, -t_gt(3), t_gt(2); ...
            t_gt(3), 0, -t_gt(1); ...
            -t_gt(2), t_gt(1), 0];
    E_gt = t_mat * R_gt;
    for idx = 1:size(pts1, 2)
        disp(pts2(:, idx)' * E_gt * pts1(:, idx))
    end

    %-----------------------------
    if size(pts1, 2) < 5
        continue
    end

    % pts_cell_1{i} = pts1;
    % pts_cell_2{i} = pts2;

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
%     tic
%     [theta_gbnb_opt, t_gbnb_opt, ~] = GBnB(pts1, pts2, R_v, epsilon);
%     runtime_gbnb_opt(i) = toc;
%     theta_err_gbnb_opt = abs(theta_gt - theta_gbnb_opt) * 180 / pi;
%     disp(["theta_err_gbnb_opt:", theta_err_gbnb_opt])
%     t_err_gbnb_opt = real(acos(abs(t_gbnb_opt' * t_gt))) * 180 / pi;
%     disp(["t_err_gbnb_opt:", t_err_gbnb_opt])
% 
%     %-------------------------gBnB+stabbing---------------------------%
%     tic
%     [t_gbnb_stabbing, theta_gbnb_stabbing] = solve_BnB_stabbing(pts1, pts2, R_v, epsilon);
%     runtime_gbnb_stabbing(i) = toc;
%     theta_err_gbnb_stabbing = abs(theta_gt - theta_gbnb_stabbing) * 180 / pi;
%     disp(["theta_err_gbnb_stabbing", theta_err_gbnb_stabbing])
%     t_err_gbnb_stabbing = real(acos(abs(t_gbnb_stabbing' * t_gt))) * 180 / pi; 
%     disp(["t_err_gbnb_stabbing", t_err_gbnb_stabbing])

end

