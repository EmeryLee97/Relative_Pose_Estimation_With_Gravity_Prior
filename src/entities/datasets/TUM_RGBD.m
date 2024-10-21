classdef TUM_RGBD < MonocularDataset
    properties
        poses (:, 4, 4) {mustBeA(poses, {'int', 'float'})}
    end

    methods (Access=public)
        function obj = TUM_RGBD(scene_config_path)
            % Constructor
            obj@MonocularDataset(scene_config_path);
            frame_rate = 32;
            obj = obj.load_tum(frame_rate);
        end

        function obj = load_tum(obj, frame_rate)
            % load and asssociate rgb images, depth images and poses in TUM-RGBD format
            if nargin < 2
                frame_rate = -1;
            end

            if exist(fullfile(obj.dataset_path, 'groundtruth.txt'), 'file') == 2
                pose_path = fullfile(obj.dataset_path, 'groundtruth.txt');
            elseif exist(fullfile(obj.dataset_path, 'pose.txt'), 'file') == 2
                pose_path = fullfile(obj.dataset_path, 'pose.txt');
            end
            rgb_path = fullfile(obj.dataset_path, 'rgb.txt');
            depth_path = fullfile(obj.dataset_path, 'depth.txt');

            rgb_data = TUM_RGBD.parse_list(rgb_path, 'rgb', 3);
            depth_data = TUM_RGBD.parse_list(depth_path, 'depth', 3);
            pose_data = cell2mat(TUM_RGBD.parse_list(pose_path, 'pose', 3));
            pose_vecs = pose_data(:, 2:end); % mat

            tstamp_rgb = rgb_data{1};
            tstamp_depth = depth_data{1};
            tstamp_pose = pose_data(:, 1);
            associations = TUM_RGBD.associate_frames(tstamp_rgb, tstamp_depth, tstamp_pose, 0.08);

            indices = 1;
            for idx = 2:size(associations, 1)
                t0 = tstamp_rgb(associations(indices(end), 1));
                t1 = tstamp_rgb(associations(idx, 1));
                if t1 - t0 > 1 / frame_rate
                    indices(end + 1) = idx;
                end
            end

            obj.poses = zeros(length(indices), 4, 4);
            for idx = 1:length(indices)
                i = associations(indices(idx), 1);
                j = associations(indices(idx), 2);
                k = associations(indices(idx), 3);

                obj.color_paths(end+1) = fullfile(obj.dataset_path, rgb_data{2}{i});
                obj.depth_paths(end+1) = fullfile(obj.dataset_path, depth_data{2}{j});

                c2w = TUM_RGBD.transquat2pose(pose_vecs(k, :));
                obj.poses(idx, :, :) = c2w;
            end
        end
        
        function [rgb, depth, pose] = get_rgb_depth_pose(obj, idx)
            % Return rgb image, depth image, and camera pose with given index
            rgb = imread(obj.color_paths(idx));
            depth = imread(obj.depth_paths(idx)) / obj.depth_scale;
            pose = squeeze(obj.poses(idx, :, :));

            if isfield(obj.dataset_config, 'distortion')
                rgb = undistortImage(rgb, obj.camera_params);
            end

            if obj.crop_edge > 0
                rgb = rgb(obj.crop_edge+1:end-obj.crop_edge, obj.crop_edge+1:end-obj.crop_edge, :);
                depth = depth(obj.crop_edge+1:end-obj.crop_edge, obj.crop_edge+1:end-obj.crop_edge);
            end
        end
    end

    methods (Static)
        function associations = associate_frames(tstamp_rgb, tstamp_depth, tstamp_pose, max_dt)
            % associate rgb, depth, pose according to their timestamps
            if nargin < 4
                max_dt = 0.08;
            end
            associations = [];
            for i = 1:length(tstamp_rgb)
                t = tstamp_rgb(i);
                [~, j] = min(abs(tstamp_depth - t));
                [~, k] = min(abs(tstamp_pose - t));
                if abs(tstamp_depth(j) - t) < max_dt && abs(tstamp_pose(k) - t) < max_dt
                    associations = [associations; i, j, k];
                end
            end
        end

        function data = parse_list(file_path, data_name, skip_rows)
            % skip the comment rows, 3 in our case
            if nargin < 3
                skip_rows = 3;
            end
            file_id = fopen(file_path, 'r');
            if strcmp(data_name, 'rgb') || strcmp(data_name, 'depth')
                data = textscan(file_id, '%f %s', 'Delimiter', ' ', 'HeaderLines', skip_rows);
            elseif strcmp(data_name, 'pose')
                data = textscan(file_id, '%f %f %f %f %f %f %f %f', 'Delimiter', ' ', 'HeaderLines', skip_rows);
            end
            fclose(file_id);
        end

        function pose = transquat2pose(trans_quat)
            % trans_quat: [tx, ty, tz, qx, qy, qz, qw]
            trans_quat = [trans_quat(end), trans_quat(1:end-1)];
            pose = eye(4);
            pose(1:3, 1:3) = quat2rotm(trans_quat(4:end));
            pose(1:3, 4) = trans_quat(1:3);
        end
    end
end