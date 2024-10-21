function test_realworld(cfg_path, dataset_name)

    addpath('../entities/datasets')
    addpath('../entities/gBnB_opt/')
    addpath('../entities/gBnB_sampling/')
    addpath('../entities/gBnB_stabbing/')
    addpath('../entities/RANSAC+5pt')
    addpath('../entities/RANSAC+3pt')

    if strcmp(dataset_name, 'TUM_RGBD')
        dataset = TUM_RGBD(cfg_path);
    elseif strcmp(dataset_name, 'KITTI')
        dataset = KITTI(cfg_path);
    elseif strcmp(dataset_name, 'ScanNet')
        dataset = ScanNet(cfg_path);
    else
        fprintf('Unsupported Dataset!')
        return;
    end
    
    % initialization
    runtime_gBnB = zeros(1, dataset.frame_pairs);
    runtime_SaBnB = zeros(1, dataset.frame_pairs);
    runtime_ISBnB = zeros(1, dataset.frame_pairs);
    runtime_RANSAC_3pt = zeros(dataset.ransac_repeat, dataset.frame_pairs);
    runtime_RANSAC_5pt = zeros(dataset.ransac_repeat, dataset.frame_pairs);
    
    rot_err_gBnB = zeros(1, dataset.frame_pairs);
    rot_err_SaBnB = zeros(1, dataset.frame_pairs);
    rot_err_ISBnB = zeros(1, dataset.frame_pairs);
    rot_err_RANSAC_3pt = zeros(dataset.ransac_repeat, dataset.frame_pairs);
    rot_err_RANSAC_5pt = zeros(dataset.ransac_repeat, dataset.frame_pairs);
    
    trans_err_gBnB = zeros(1, dataset.frame_pairs);
    trans_err_SaBnB = zeros(1, dataset.frame_pairs);
    trans_err_ISBnB = zeros(1, dataset.frame_pairs);
    trans_err_RANSAC_3pt = zeros(dataset.ransac_repeat, dataset.frame_pairs);
    trans_err_RANSAC_5pt = zeros(dataset.ransac_repeat, dataset.frame_pairs);

    outlier_rates = zeros(1, dataset.frame_pairs);
    pts_nums = zeros(1, dataset.frame_pairs);

    % Run experiment
    for i = dataset.frame_start : dataset.frame_start + dataset.frame_pairs -1

        fprintf("Processing image pair %d\n", i)

        % Get rgb images and absolute camera poses
        [rgb_1, ~, pose_1] = dataset.get_rgb_depth_pose(i);
        [rgb_2, ~, pose_2] = dataset.get_rgb_depth_pose(i + dataset.frame_stride);
        gray_1 = rgb2gray(rgb_1); 
        gray_2 = rgb2gray(rgb_2);
        
        % Calculate the relative pose between them as ground truth
        rot_1 = pose_1(1:3, 1:3); trans_1 = pose_1(1:3, 4);
        rot_2 = pose_2(1:3, 1:3); trans_2 = pose_2(1:3, 4);
        rot_gt = rot_2' * rot_1; trans_gt = rot_2' * (trans_1 - trans_2);

        rot_axang = rotationMatrixToVector(rot_gt);
        rot_angle_gt = norm(rot_axang);
        rot_axis = rot_axang ./ rot_angle_gt;
        rot_axis_mat = [0          -rot_axis(3) rot_axis(2); ...
                        rot_axis(3) 0          -rot_axis(1); ...
                       -rot_axis(2) rot_axis(1) 0]';
        trans_gt = trans_gt / norm(trans_gt);

        % Detect features
        point_set_1 = detectSIFTFeatures(gray_1);
        point_set_2 = detectSIFTFeatures(gray_2);
    
        [f1, vpts1] = extractFeatures(gray_1, point_set_1);
        [f2, vpts2] = extractFeatures(gray_2, point_set_2);
    
        %indexPairs = matchFeatures(f1, f2) 
        indexPairs = matchFeatures( ...
            f1, f2, ...
            'MatchThreshold', dataset.match_threshold, ...
            'MaxRatio', dataset.max_ratio ...
            );
        matched_pts_1 = vpts1(indexPairs(:, 1));
        matched_pts_2 = vpts2(indexPairs(:, 2));
    
        pts_1 = matched_pts_1.Location;
        pts_2 = matched_pts_2.Location;
    
        pts_1 = [pts_1'; ones(1, size(pts_1, 1))];
        pts_2 = [pts_2'; ones(1, size(pts_2, 1))];
    
        pts_1 = dataset.intrinsics \ pts_1; pts_1 = pts_1 ./ vecnorm(pts_1);
        pts_2 = dataset.intrinsics \ pts_2; pts_2 = pts_2 ./ vecnorm(pts_2);

        % Calculate the inlier rate
        b = cross(pts_2, rot_axis_mat * pts_1);
        c = -cross(pts_2, rot_axis_mat^2 * pts_1);
        a = cross(pts_2, pts_1) + cross(pts_2, rot_axis_mat^2 * pts_1);
    
        residual = abs(trans_gt' * (a + b*sin(rot_angle_gt) + c*cos(rot_angle_gt)));
        outlier_rate = sum(residual > dataset.epsilon) / size(pts_1, 2);
        outlier_rates(i) = outlier_rate;
        pts_nums(i) = size(pts_1, 2);
        fprintf("Number of points: %d\n", size(pts_1, 2));
        fprintf("Outlier rate: %d\n", outlier_rate);

    end

end




