
classdef RealworldEvaluator
    properties
        dataset {mustBeA(dataset, 'BaseDataset')}
        logger {mustBeA(logger, 'Logger')}
        method_list {mustBeA(method_list, 'string')}
        feature_type {mustBeText}
        ransac_rho {mustBeA(ransac_rho, {'single', 'double'})}
        ransac_iter {mustBeInteger}
        ransac_repeat {mustBeInteger}
        frame_start {mustBeInteger}
        frame_stride {mustBeInteger}
        frame_pairs {mustBeInteger}
        epsilon {mustBeA(epsilon, {'single', 'double'})}
        match_threshold {mustBeInteger}
        max_ratio {mustBeA(max_ratio, {'single', 'double'})}
      
        outlier_rates (1, :) {mustBeA(outlier_rates, {'int', 'single', 'double'})}
        pts_nums (1, :) {mustBeInteger}
    end

    methods (Access=public)
        function obj = RealworldEvaluator(cfg_path, dataset_name)
            % Constructor

            addpath('../entities/datasets')
            addpath('../entities/gBnB_opt/')
            addpath('../entities/gBnB_sampling/')
            addpath('../entities/gBnB_stabbing/')
            addpath('../entities/RANSAC+3pt')
            addpath('../entities/RANSAC+5pt')           
            
            if strcmpi(dataset_name, 'TUM_RGBD')
                obj.dataset = TUM_RGBD(cfg_path);
            elseif strcmpi(dataset_name, 'KITTI')
                obj.dataset = KittiOdometry(cfg_path);
            elseif strcmpi(dataset_name, 'ScanNet')
                obj.dataset = ScanNet(cfg_path);
            else
                fprintf('Unsupported Dataset!')
            end
            
            obj.method_list = ["gBnB", "SaBnB", "ISBnB", "RANSAC_3pt", "RANSAC_5pt"];
            frame_pair_list = ones(1, length(obj.method_list)) * obj.frame_pairs;
            repeat_list = [1, 1, 1, obj.ransac_repeat, obj.ransac_repeat];
            obj.logger = Logger(obj.method_list, repeat_list, frame_pair_list);

            % TODO: add feature type in config file
            obj.feature_type = obj.dataset.dataset_config.feature_type;
            obj.ransac_rho = obj.dataset.dataset_config.ransac.rho;
            obj.ransac_iter = obj.dataset.dataset_config.ransac.iter;
            obj.ransac_repeat = obj.dataset.dataset_config.ransac.repeat;
            obj.frame_start = obj.dataset.dataset_config.frame_start;
            obj.frame_stride = obj.dataset.dataset_config.frame_stride;
            obj.frame_pairs = obj.dataset.dataset_config.frame_pairs;
            obj.epsilon = obj.dataset.dataset_config.epsilon;
            obj.match_threshold = obj.dataset.dataset_config.match_threshold;
            obj.max_ratio = obj.dataset.dataset_config.max_ratio;
            
            obj.outlier_rates = zeros(1, obj.frame_pairs);
            obj.pts_nums = zeros(1, obj.frame_pairs);
        end

        function [pts_1, pts_2] = detect_and_match(obj, gray_1, gray_2)
            % gray_1, gray_2: gray scale images
            % feature_name: SIFT, SURF, FAST, ORB, etc.
            if strcmpi(obj.feature_type, 'SIFT')
                point_set_1 = detectSIFTFeatures(gray_1);
                point_set_2 = detectSIFTFeatures(gray_2);
            elseif strcmpi(obj.feature_type, 'SURF')
                point_set_1 = detectSURFFeatures(gray_1);
                point_set_2 = detectSURFFeatures(gray_2);
            elseif strcmpi(obj.feature_type, 'ORB')
                point_set_1 = detectORBFeatures(gray_1);
                point_set_2 = detectORBFeatures(gray_2);
            else
                error("Unsupported Feature Type!")
            end

            [f1, vpts_1] = extractFeatures(gray, point_set_1);
            [f2, vpts_2] = extractFeatures(gray, point_set_2);

            indexPairs = matchFeatures(f1, f2, 'MatchThreshold', obj.match_threshold, 'MaxRatio', obj.max_ratio);

            pts_1 = vpts_1(indexPairs(:, 1)).Location;
            pts_2 = vpts_2(indexPairs(:, 2)).Location;

            pts_1 = [pts_1'; ones(1, size(pts_1, 1))];
            pts_2 = [pts_2'; ones(1, size(pts_2, 1))];

            pts_1 = obj.dataset.intrinsics \ pts_1; 
            pts_1 = pts_1 ./ vecnorm(pts_1);
            pts_2 = obj.dataset.intrinsics \ pts_2; 
            pts_2 = pts_2 ./ vecnorm(pts_2);
        end

        function obj = run_gBnB(obj, pts_1, pts_2, rot_axis, rot_angle_gt, trans_gt, frame_idx)
            tic
            [rot_angle, trans, ~] = GBnB(pts_1, pts_2, rot_axis, obj.epsilon);
            obj.logger.dynamic_vars.runtime_gBnB(frame_idx) = toc;
            obj.logger.dynamic_vars.rot_err_gBnB(frame_idx) = RealworldEvaluator.get_rot_angle_err(rot_angle, rot_angle_gt);     
            obj.logger.dynamic_vars.trans_err_gBnB(frame_idx) = RealworldEvaluator.get_trans_err(trans, trans_gt);
        end

        function run_SaBnB(obj, num_samples, pts_1, pts_2, rot_axis, rot_angle_gt, trans_gt, frame_idx)
            tic
            [rot_angle, trans] = solve_BnB_sampling(num_samples, pts_1, pts_2, rot_axis, obj.epsilon);
            obj.logger.dynamic_vars.runtime_SaBnB(frame_idx) = toc;
            obj.logger.dynamic_vars.rot_err_SaBnB(frame_idx) = RealworldEvaluator.get_rot_angle_err(rot_angle, rot_angle_gt);    
            obj.logger.dynamic_vars.trans_err_SaBnB(frame_idx) = RealworldEvaluator.get_trans_err(trans, trans_gt);
        end

        function run_ISBnB(obj, pts_1, pts_2, rot_axis, rot_angle_gt, trans_gt, frame_idx)
            tic
            [rot_angle, trans] = solve_BnB_stabbing(pts_1, pts_2, rot_axis, obj.epsilon);
            obj.logger.dynamic_vars.runtime_ISBnB(frame_idx) = toc;
            obj.logger.dynamic_vars.rot_err_ISBnB(frame_idx) = RealworldEvaluator.get_rot_angle_err(rot_angle, rot_angle_gt);    
            obj.logger.dynamic_vars.trans_err_ISBnB(frame_idx) = RealworldEvaluator.get_trans_err(trans, trans_gt);
        end

        function run_RANSAC_3pt(obj, pts_1, pts_2, rot_axis, rot_angle_gt, trans_gt, frame_idx, repeat_idx)
            tic
            [rot_angle, trans] = ransac_3pt(pts_1, pts_2, rot_axis, obj.epsilon, obj.ransac_iter);
            obj.logger.dynamic_vars.runtime_RANSAC_3pt(repeat_idx, frame_idx) = toc;
            obj.logger.dynamic_vars.rot_err_RANSAC_3pt(repeat_idx, frame_idx) = RealworldEvaluator.get_rot_angle_err(rot_angle, rot_angle_gt);    
            obj.logger.dynamic_vars.trans_err_RANSAC_3pt(repeat_idx, frame_idx) = RealworldEvaluator.get_trans_err(trans, trans_gt);
        end

        function run_RANSAC_5pt(obj, pts_1, pts_2, rot_axis, rot_gt, trans_gt, frame_idx, repeat_idx)
            tic
            [rot, trans] = ransac_5pt(pts_1, pts_2, rot_axis, obj.epsilon, obj.ransac_iter);
            obj.logger.dynamic_vars.runtime_RANSAC_5pt(repeat_idx, frame_idx) = toc;
            obj.logger.dynamic_vars.rot_err_RANSAC_5pt(repeat_idx, frame_idx) = RealworldEvaluator.get_rot_err(rot, rot_gt);    
            obj.logger.dynamic_vars.trans_err_RANSAC_5pt(repeat_idx, frame_idx) = RealworldEvaluator.get_trans_err(trans, trans_gt);
        end

        function obj = run(obj)

            for i = obj.frame_start : obj.frame_start+obj.frame_pairs-1
                fprintf("Processing image pair %d\n", i)

                % Get gray scale images and absolute camera poses
                [gray_1, pose_1] = obj.dataset.get_gray_pose(i);
                [gray_2, pose_2] = obj.dataset.get_gray_pose(i+obj.frame_stride);

                % Get relative pose
                [rot_axis_gt, rot_angle_gt, rot_gt, trans_gt] = RealworldEvaluator.get_relative_pose(pose_1, pose_2);

                % Detect & compute & match features
                [pts_1, pts_2] = obj.detect_and_match(gray_1, gray_2);

                % Compute number of matchings and outlier rate
                [obj.outlier_rates(i), obj.pts_nums(i)] = RealworldEvaluator.get_outlier_rate( ...
                    pts_1, pts_2, rot_axis_gt, obj.epsilon);

                % Run algorithms
                obj.run_gBnB(pts_1, pts_2, rot_axis, rot_axis_gt, rot_angle_gt, trans_gt, i);
                obj.run_SaBnB(360, pts_1, pts_2, rot_axis, rot_axis_gt, rot_angle_gt, trans_gt, i);
                obj.run_ISBnB(pts_1, pts_2, rot_axis, rot_axis_gt, rot_angle_gt, trans_gt, i);

                for j = 1:obj.ransac_repeat
                    obj.run_RANSAC_3pt(pts_1, pts_2, rot_axis_gt, rot_angle_gt, trans_gt, i, j);
                    obj.run_RANSAC_5pt(pts_1, pts_2, rot_axis_gt, rot_gt, trans_gt, i, j);
                end                
            end

            obj.logger = obj.logger.update_dict();
            fprintf("ALl DONE!")
        end
    end

    methods (Static)
        function [rot_axis_gt, rot_angle_gt, rot_gt, trans_gt] = get_relative_pose(pose_1, pose_2)
            % pose_1, pose_2: absolute poses (cam -> world)
            rot_1 = pose_1(1:3, 1:3); 
            trans_1 = pose_1(1:3, 4);

            rot_2 = pose_2(1:3, 1:3); 
            trans_2 = pose_2(1:3, 4);

            rot_gt = rot_2' * rot_1;  
            trans_gt = rot_2' * (trans_1 - trans_2);

            rot_axang = rotationMatrixToVector(rot_gt);
            rot_angle_gt = norm(rot_axang);
            rot_axis = rot_axang ./ rot_angle_gt;
            rot_axis_gt = [0          -rot_axis(3) rot_axis(2); ...
                            rot_axis(3) 0          -rot_axis(1); ...
                           -rot_axis(2) rot_axis(1) 0]';
            trans_gt = trans_gt / norm(trans_gt);
        end
    
        function [point_num, outlier_rate] = get_outlier_rate(pts_1, pts_2, rot_axis, epsilon)
            b = cross(pts_2, rot_axis * pts_1);
            c = -cross(pts_2, rot_axis^2 * pts_1);
            a = cross(pts_2, pts_1) + cross(pts_2, rot_axis^2 * pts_1);
            residual = abs(trans_gt' * (a + b*sin(rot_angle_gt) + c*cos(rot_angle_gt)));
            outlier_rate = sum(residual > epsilon) / size(pts_1, 2);
            point_num = size(pts_1, 2);
            fprintf("Number of points: %d\n", point_num);
            fprintf("Outlier rate: %d\n", outlier_rate);
        end

        function trans_err = get_trans_err(trans_est, trans_gt)
            trans_err = real(acos(abs(trans_est' * trans_gt))) * 180 / pi;
        end

        function rot_angle_err = get_rot_angle_err(rot_angle_est, rot_angle_gt)
            rot_angle_err = abs(rot_angle_est - rot_angle_gt) * 180 / pi;
        end

        function rot_err = get_rot_err(rot_est, rot_gt)
            rot_err = norm(rotationMatrixToVector(rot_est' * rot_gt)) * 180 / pi;
        end
        
    end
  
end





