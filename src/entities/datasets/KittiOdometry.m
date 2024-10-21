% https://s3.eu-central-1.amazonaws.com/avg-kitti/data_odometry_gray.zip
% calib.txt: P_i = K_i[R_i|t_i] ([R_i|t_i]: cam_0 -> cam_i)
% ground truth poses: T_w->cam0 (cam0 = w when at t_0)

classdef KittiOdometry < StereoDataset
    properties
        % TODO: 
        sequence string
        poses (:, 4, 4) {mustBeA(poses, {'int', 'float'})}
        gray_cam0_paths string
        gray_cam1_paths string      
    end

    methods (Access=public)
        function obj = KittiOdometry(scene_config_path)
            % Constructor
            obj@StereoDataset(scene_config_path)

            obj.sequence = num2str(obj.dataset_config.data.sequence);

            obj = obj.load_kitti;        
        end

        function obj = load_kitti(obj)
            % load gray_cam0, gray_cam1 and poses

            cam0_files = dir(fullfile( ...
                obj.dataset_path, 'sequences', obj.sequence, 'image_0', '*.png'));
            cam1_files = dir(fullfile( ...
                obj.dataset_path, 'sequences', obj.sequence, 'image_1', '*.png'));

            obj.gray_cam0_paths = strings(length(cam0_files), 1);
            obj.gray_cam1_paths = strings(length(cam1_files), 1);
            
            for i = 1:length(cam0_files)
                obj.gray_cam0_paths(i) = fullfile(cam0_files(i).folder, cam0_files(i).name);
                obj.gray_cam1_paths(i) = fullfile(cam1_files(i).folder, cam1_files(i).name);
            end

            pose_path = fullfile(obj.dataset_path, 'poses', strcat(obj.sequence, '.txt'));
            pose_id = fopen(pose_path, 'r');
            pose = fscanf(fileid_pose, '%f');
            fclose(pose_id);

            obj.poses = zeros(length(pose)/12, 3, 4);
            for i = 1:length(pose) / 12
                obj.poses(i, :, :) = [
                pose((i-1)*12+1 : (i-1)*12+4)'; ...
                pose((i-1)*12+5 : (i-1)*12+8)'; ...
                pose((i-1)*12+9 : i*12)'];
            end

        end

        function [gray_0, gray_1, pose] = get_gray_pose(obj, idx)
            gray_0 = imread(obj.gray_cam0_paths(idx));
            gray_1 = imread(obj.gray_cam1_paths(idx));
            pose = squeeze(obj.poses(idx, :, :));
        end
    end
end