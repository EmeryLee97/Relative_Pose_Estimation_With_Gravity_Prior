
classdef RealworldEvaluator
    properties
        dataset BaseDataset
        logger Logger
        method_list string
        
        % gcransac_module
        % magsacpp_module
        usac_module

        sample_num {mustBeInteger}
        feature_type string
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

            addpath('src/entities/')
            addpath('src/entities/datasets')
            addpath('src/utils')
            addpath('src/entities/gBnB_opt/')
            addpath('src/entities/gBnB_sampling/')
            addpath('src/entities/gBnB_stabbing/')
            addpath('src/entities/RANSAC+3pt')
            addpath('src/entities/RANSAC+5pt') 
            addpath('src/entities/USAC/')
            % addpath('src/entities/GCRANSAC')
            % addpath('src/entities/MAGSACPP')
            
            sys_path = string(py.sys.path);
            if ~any(strcmp(sys_path, fullfile(pwd, 'src/entities/USAC')))
                py.sys.path().append(fullfile(pwd, 'src/entities/USAC'))
            end
            obj.usac_module = py.importlib.import_module('usac');
            % obj.gcransac_module = py.importlib.import_module('gcransac');
            % obj.magsacpp_module = py.importlib.import_module('magsac');
            
            if strcmpi(dataset_name, 'TUM_RGBD')
                obj.dataset = TUM_RGBD(cfg_path);
            elseif strcmpi(dataset_name, 'KITTI')
                obj.dataset = KittiOdometry(cfg_path);
            elseif strcmpi(dataset_name, 'ScanNet')
                obj.dataset = ScanNet(cfg_path);
            else
                fprintf('Unsupported Dataset!')
            end
            
            obj.method_list = ["gBnB", "SaBnB", "ISBnB", "RANSAC_3pt", "RANSAC_5pt", "GCRANSAC", "MAGSACPP"];
            obj.sample_num = obj.dataset.dataset_config.sample_num;
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
            frame_pair_list = ones(1, length(obj.method_list)) .* obj.frame_pairs;
            repeat_list = [1, 1, 1, obj.ransac_repeat, obj.ransac_repeat, 1, 1];
            obj.logger = Logger(obj.method_list, repeat_list, frame_pair_list);
        end

        function [pts_1, pts_2] = detect_and_match(obj, gray_1, gray_2)
            % gray_1, gray_2: gray scale images
            % feature_name: SIFT, SURF, ORB, MSER, etc.
            % Binary features: FREAK, ORB, BRISK
            if strcmpi(obj.feature_type, 'SIFT')
                point_set_1 = detectSIFTFeatures(gray_1);
                point_set_2 = detectSIFTFeatures(gray_2);
            elseif strcmpi(obj.feature_type, 'SURF')
                point_set_1 = detectSURFFeatures(gray_1);
                point_set_2 = detectSURFFeatures(gray_2);
            elseif strcmpi(obj.feature_type, 'ORB')
                point_set_1 = detectORBFeatures(gray_1);
                point_set_2 = detectORBFeatures(gray_2);
            elseif strcmpi(obj.feature_type, 'MSER')
                point_set_1 = detectMSERFeatures(gray_1);
                point_set_2 = detectMSERFeatures(gray_2);
            else
                error("Unsupported Feature Type!")
            end

            [f1, vpts_1] = extractFeatures(gray_1, point_set_1);
            [f2, vpts_2] = extractFeatures(gray_2, point_set_2);

            indexPairs = matchFeatures(f1, f2, 'MatchThreshold', obj.match_threshold, 'MaxRatio', obj.max_ratio);

            pts_1 = vpts_1(indexPairs(:, 1)).Location; % (N, 2)
            pts_2 = vpts_2(indexPairs(:, 2)).Location; % (N, 2)

            pts_1 = [pts_1'; ones(1, size(pts_1, 1))]; % (3, N)
            pts_2 = [pts_2'; ones(1, size(pts_2, 1))]; % (3, N)

            if isprop(obj.dataset, 'intrinsics')
                pts_1 = obj.dataset.intrinsics \ pts_1;
                pts_2 = obj.dataset.intrinsics \ pts_2;
            elseif isprop(obj.dataset, 'intrinsics_cam0')
                pts_1 = obj.dataset.intrinsics_cam0 \ pts_1;
                pts_2 = obj.dataset.intrinsics_cam0 \ pts_2;
            else
                fprintf("No camera intrinsics!")
            end

            pts_1 = pts_1 ./ vecnorm(pts_1); % (3, N)
            pts_2 = pts_2 ./ vecnorm(pts_2); % (3, N)
        end

        function obj = run_gBnB(obj, pts_1, pts_2, rot_axis, rot_angle_gt, trans_gt, frame_idx)
            tic
            [rot_angle, trans, ~] = GBnB(pts_1, pts_2, rot_axis, obj.epsilon);
            %------------------------------------
            fprintf("rot_est = %d, rot_gt= %d \n", rot_angle, rot_angle_gt)
            obj.logger.dynamic_vars.runtime_gBnB(frame_idx) = toc;
            obj.logger.dynamic_vars.rot_err_gBnB(frame_idx) = RealworldEvaluator.get_rot_angle_err(rot_angle, rot_angle_gt);     
            obj.logger.dynamic_vars.trans_err_gBnB(frame_idx) = RealworldEvaluator.get_trans_err(trans, trans_gt);
        end

        function obj = run_SaBnB(obj, num_samples, pts_1, pts_2, rot_axis, rot_angle_gt, trans_gt, frame_idx)
            tic
            [trans, rot_angle] = solve_BnB_sampling(num_samples, pts_1, pts_2, rot_axis, obj.epsilon);
            fprintf("rot_est = %d, rot_gt= %d \n", rot_angle, rot_angle_gt)
            obj.logger.dynamic_vars.runtime_SaBnB(frame_idx) = toc;
            obj.logger.dynamic_vars.rot_err_SaBnB(frame_idx) = RealworldEvaluator.get_rot_angle_err(rot_angle, rot_angle_gt);    
            obj.logger.dynamic_vars.trans_err_SaBnB(frame_idx) = RealworldEvaluator.get_trans_err(trans, trans_gt);
        end

        function obj = run_ISBnB(obj, pts_1, pts_2, rot_axis, rot_angle_gt, trans_gt, frame_idx)
            tic
            [trans, rot_angle] = solve_BnB_stabbing(pts_1, pts_2, rot_axis, obj.epsilon);
            fprintf("rot_est = %d, rot_gt= %d \n", rot_angle, rot_angle_gt)
            obj.logger.dynamic_vars.runtime_ISBnB(frame_idx) = toc;
            obj.logger.dynamic_vars.rot_err_ISBnB(frame_idx) = RealworldEvaluator.get_rot_angle_err(rot_angle, rot_angle_gt);    
            obj.logger.dynamic_vars.trans_err_ISBnB(frame_idx) = RealworldEvaluator.get_trans_err(trans, trans_gt);
        end

        function obj = run_RANSAC_3pt(obj, pts_1, pts_2, rot_axis, rot_angle_gt, trans_gt, frame_idx, repeat_idx)
            tic
            [rot_angle, trans] = ransac_3pt(pts_1, pts_2, rot_axis, obj.epsilon, obj.ransac_iter);
            obj.logger.dynamic_vars.runtime_RANSAC_3pt(repeat_idx, frame_idx) = toc;
            obj.logger.dynamic_vars.rot_err_RANSAC_3pt(repeat_idx, frame_idx) = RealworldEvaluator.get_rot_angle_err(rot_angle, rot_angle_gt);    
            obj.logger.dynamic_vars.trans_err_RANSAC_3pt(repeat_idx, frame_idx) = RealworldEvaluator.get_trans_err(trans, trans_gt);
        end

        function obj = run_RANSAC_5pt(obj, pts_1, pts_2, rot_axis, rot_gt, trans_gt, frame_idx, repeat_idx)
            tic
            [rot, trans] = ransac_5pt(pts_1, pts_2, rot_axis, obj.epsilon, obj.ransac_iter);
            obj.logger.dynamic_vars.runtime_RANSAC_5pt(repeat_idx, frame_idx) = toc;
            obj.logger.dynamic_vars.rot_err_RANSAC_5pt(repeat_idx, frame_idx) = RealworldEvaluator.get_rot_err(rot, rot_gt);    
            obj.logger.dynamic_vars.trans_err_RANSAC_5pt(repeat_idx, frame_idx) = RealworldEvaluator.get_trans_err(trans, trans_gt);
        end

        function obj = run_GCRANSAC(obj, pts_1, pts_2, rot_gt, trans_gt, frame_idx)            
            % reshape 3D point sets from (N, 3) to 2D point sets (2, N)
            pts1_2d = pts_1 ./ pts_1(3, :); 
            pts1_2d = pts1_2d(1:2, :)';

            pts2_2d = pts_2 ./ pts_2(3, :); 
            pts2_2d = pts2_2d(1:2, :)';

            tic
            % pose_tuple = obj.gcransac_module.run_gcransac(py.numpy.array(pts1_2d), py.numpy.array(pts2_2d));
            pose_tuple = obj.usac_module.run_gcransac(py.numpy.array(pts1_2d), py.numpy.array(pts2_2d));
            obj.logger.dynamic_vars.runtime_GCRANSAC(frame_idx) = toc;
            rot = double(pose_tuple{1});
            trans = double(pose_tuple{2});
            obj.logger.dynamic_vars.rot_err_GCRANSAC(frame_idx) = RealworldEvaluator.get_rot_err(rot, rot_gt);
            obj.logger.dynamic_vars.trans_err_GCRANSAC(frame_idx) = RealworldEvaluator.get_trans_err(trans, trans_gt);
        end

        function obj = run_MAGSACPP(obj, pts_1, pts_2, rot_gt, trans_gt, frame_idx)
            % reshape 3D point sets from (N, 3) to 2D point sets (2, N)
            pts1_2d = pts_1 ./ pts_1(3, :); 
            pts1_2d = pts1_2d(1:2, :)';

            pts2_2d = pts_2 ./ pts_2(3, :); 
            pts2_2d = pts2_2d(1:2, :)';

            tic
            % pose_tuple = obj.magsacpp_module.run_magsacpp(py.numpy.array(pts1_2d), py.numpy.array(pts2_2d));
            pose_tuple = obj.usac_module.run_magsacpp(py.numpy.array(pts1_2d), py.numpy.array(pts2_2d));
            obj.logger.dynamic_vars.runtime_MAGSACPP(frame_idx) = toc;
            rot = double(pose_tuple{1});
            trans = double(pose_tuple{2});
            obj.logger.dynamic_vars.rot_err_MAGSACPP(frame_idx) = RealworldEvaluator.get_rot_err(rot, rot_gt);
            obj.logger.dynamic_vars.trans_err_MAGSACPP(frame_idx) = RealworldEvaluator.get_trans_err(trans, trans_gt);
        end

        function print_result(obj)
            obj.logger.print_result();
        end

        function obj = evaluate(obj, test_outlier)

            for i = obj.frame_start : obj.frame_start+obj.frame_pairs-1
                fprintf('\n')
                fprintf("Processing image pair %d-%d\n", i, i+obj.frame_stride)

                % Get gray scale images and absolute camera poses
                [gray_1, pose_1, ~] = obj.dataset.get_gray_pose(i);
                [gray_2, pose_2, ~] = obj.dataset.get_gray_pose(i+obj.frame_stride);

                % Get relative pose
                [rot_axis_gt, rot_angle_gt, rot_gt, trans_gt] = RealworldEvaluator.get_relative_pose(pose_1, pose_2);
                rot_axis_mat = cross_product_matrix(rot_axis_gt);

                % Detect & compute & match features
                [pts_1, pts_2] = obj.detect_and_match(gray_1, gray_2);
                % RealworldEvaluator.check_epipolar_constraint(pts_1, pts_2, rot_gt, trans_gt)

                % Compute number of matchings and outlier rate
                [obj.outlier_rates(i), obj.pts_nums(i)] = RealworldEvaluator.get_outlier_rate( ...
                    pts_1, pts_2, rot_axis_mat, rot_angle_gt, trans_gt, obj.epsilon);

                if ~test_outlier
                    % Run algorithms
                    % fprintf('Running gBnB \n')
                    % obj = obj.run_gBnB(pts_1, pts_2, rot_axis_gt, rot_angle_gt, trans_gt, i-obj.frame_start+1);

                    fprintf('Running SaBnB \n')
                    disp(obj.sample_num)
                    obj = obj.run_SaBnB(obj.sample_num, pts_1, pts_2, rot_axis_gt, rot_angle_gt, trans_gt, i-obj.frame_start+1);

                    fprintf('Running ISBnB \n')
                    obj = obj.run_ISBnB(pts_1, pts_2, rot_axis_gt, rot_angle_gt, trans_gt, i-obj.frame_start+1);

                    fprintf('Running GCRANSAC \n')
                    obj = obj.run_GCRANSAC(pts_1, pts_2, rot_gt, trans_gt, i-obj.frame_start+1);

                    fprintf('Running MAGSAC++ \n')
                    obj = obj.run_MAGSACPP(pts_1, pts_2, rot_gt, trans_gt, i-obj.frame_start+1);

                    fprintf('Running RANSAC&3pt \n')
                    fprintf('Running RANSAC&5pt \n')
                    for j = 1:obj.ransac_repeat             
                        obj = obj.run_RANSAC_3pt(pts_1, pts_2, rot_axis_gt, rot_angle_gt, trans_gt, i-obj.frame_start+1, j);   
                        obj = obj.run_RANSAC_5pt(pts_1, pts_2, rot_axis_gt, rot_gt, trans_gt, i-obj.frame_start+1, j);
                    end
                end

            end

            obj.logger = obj.logger.update_dict();
            fprintf("ALl DONE! \n")
        end
    end

    methods (Static)
        function [rot_axis_gt, rot_angle_gt, rot_gt, trans_gt] = get_relative_pose(T_c1w, T_c2w)

            T_c1c2 = T_c2w \ T_c1w;
            rot_gt = T_c1c2(1:3, 1:3);
            trans_gt = T_c1c2(1:3, 4);
            trans_gt = trans_gt / norm(trans_gt);

            rot_axang = rotationMatrixToVector(rot_gt); % postmultiply convention
            rot_angle_gt = norm(rot_axang);
            rot_axis_gt = -rot_axang ./ rot_angle_gt;   
            % rot_axang = rotmat2vec3d(rot_gt); % premultiply convention, need MATLAB R2022b
              
        end
    
        function [outlier_rate, point_num] = get_outlier_rate(pts_1, pts_2, rot_axis, rot_angle_gt, trans_gt, epsilon)
            b = cross(pts_2, rot_axis * pts_1);
            c = -cross(pts_2, rot_axis^2 * pts_1);
            a = cross(pts_2, pts_1) + cross(pts_2, rot_axis^2 * pts_1);
            residual = abs(trans_gt' * (a + b*sin(rot_angle_gt) + c*cos(rot_angle_gt)));
            outlier_rate = sum(residual > epsilon) / size(pts_1, 2);
            point_num = size(pts_1, 2);
            fprintf("Number of point pairs: %d, with outlier rate: %d\n", point_num, outlier_rate);
        end

        function trans_err = get_trans_err(trans_est, trans_gt)
            % Translation (unit vector) error in degree
            trans_err = real(acos(trans_est' * trans_gt / norm(trans_est) /norm(trans_gt))) * 180/pi;
            if trans_err > 90
                trans_err = abs(180 - trans_err);
            end
        end

        function rot_angle_err = get_rot_angle_err(rot_angle_est, rot_angle_gt)
            % Rotation angle error in degree
            rot_angle_err = abs(abs(rot_angle_est) - abs(rot_angle_gt)) * 180/pi;
        end

        function rot_err = get_rot_err(rot_est, rot_gt)
            rot_err = norm(rotationMatrixToVector(rot_est' * rot_gt)) * 180/pi;
            % rot_err = norm(rotmat2vec3d(rot_est' * rot_gt)) * 180/pi; % need MATLAB R2022b
        end

        function res = check_epipolar_constraint(pts_1, pts_2, rel_rot, rel_trans)
            % pts_1, pts_2: (3, N)
            rel_trans_mat = cross_product_matrix(rel_trans);
            E = rel_trans_mat * rel_rot;

            pts_num = size(pts_1, 2);
            res = zeros(1, pts_num);
            for i = 1:pts_num
                res(i) = abs(pts_2(:, i)' * E * pts_1(:, i));
            end
        end   
    end
  
end





