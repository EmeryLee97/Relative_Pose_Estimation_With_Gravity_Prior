% https://s3.eu-central-1.amazonaws.com/avg-kitti/data_odometry_gray.zip
% calib.txt: P_i = K_i[R_i|t_i] ([R_i|t_i]: cam_0 -> cam_i)
% ground truth poses: T_w->cam0 (cam0 = w when at t_0)

classdef KittiOdometry < StereoDataset
    properties
        sequence string
        poses_cam0 (:, 4, 4) {mustBeA(poses_cam0, {'int', 'float'})}
        poses_cam1 (:, 4, 4) {mustBeA(poses_cam1, {'int', 'float'})}
    end

    methods (Access=public)
        function obj = KittiOdometry(scene_config_path)
            % Constructor
            obj@StereoDataset(scene_config_path)
            if obj.dataset_config.data.sequence < 10
                obj.sequence = strcat('0', num2str(obj.dataset_config.data.sequence));
            else
                obj.sequence = num2str(obj.dataset_config.data.sequence);
            end
            obj = obj.load_kitti();        
        end

        function obj = load_kitti(obj)
            % load gray scale images and absolute poses

            cam0_files = dir(fullfile( ...
                obj.dataset_path, 'sequences', obj.sequence, 'image_0', '*.png'));
            cam1_files = dir(fullfile( ...
                obj.dataset_path, 'sequences', obj.sequence, 'image_1', '*.png'));
            
            for i = 1:length(cam0_files)
                obj.img_paths_cam0(i) = fullfile(cam0_files(i).folder, cam0_files(i).name);
                obj.img_paths_cam1(i) = fullfile(cam1_files(i).folder, cam1_files(i).name);
            end

            pose_path = fullfile(obj.dataset_path, 'poses', strcat(obj.sequence, '.txt'));
            pose_id = fopen(pose_path, 'r');
            pose = fscanf(pose_id, '%f');
            fclose(pose_id);

            obj.poses_cam0 = zeros(length(pose)/12, 4, 4);
            obj.poses_cam1 = zeros(length(pose)/12, 4, 4);
            for i = 1:length(pose) / 12
                pose_world = [
                    pose((i-1)*12+1 : (i-1)*12+4)'; ...
                    pose((i-1)*12+5 : (i-1)*12+8)'; ...
                    pose((i-1)*12+9 : i*12)'; ...
                    0, 0, 0, 1
                    ];
                obj.poses_cam0(i, :, :) = inv(obj.extrinsics_cam0 * pose_world);
                obj.poses_cam1(i, :, :) = inv(obj.extrinsics_cam1 * pose_world);          
            end

        end

        function [gray_cam0, pose_cam0, path] = get_gray_pose(obj, idx)
            gray_cam0 = imread(obj.img_paths_cam0(idx));
            pose_cam0 = squeeze(obj.poses_cam0(idx, :, :));
            path = obj.img_paths_cam0(idx);
        end

        function [gray_cam0, gray_cam1, pose_cam0, pose_cam1] = get_gray_pose_stereo(obj, idx)
            gray_cam0 = imread(obj.img_paths_cam0(idx));
            gray_cam1 = imread(obj.img_paths_cam1(idx));
            pose_cam0 = squeeze(obj.poses_cam0(idx, :, :));
            pose_cam1 = squeeze(obj.poses_cam1(idx, :, :));
        end
    end
end