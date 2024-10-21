classdef ScanNet < MonocularDataset
    properties
        poses (:, 4, 4) {mustBeA(poses, {'int', 'float'})}
    end

    methods (Access=public)
        function obj = ScanNet(scene_config_path)
            % Constructor
            obj@MonocularDataset(scene_config_path);

            color_files = ScanNet.sort_files_by_number( ...
                {dir(fullfile(obj.dataset_path, 'color', '*.jpg')).name});
            obj.color_paths = fullfile(obj.dataset_path, 'color', color_files);

            depth_files = ScanNet.sort_files_by_number( ...
                {dir(fullfile(obj.dataset_path, 'depth', '*.png')).name});
            obj.depth_paths = fullfile(obj.dataset_path, 'depth', depth_files);
            
            obj = obj.load_poses();
        end
        
        function obj = load_poses(obj)
            % Load camera poses with given path
            obj.poses = zeros(numel(obj.color_paths), 4, 4);
            pose_files = ScanNet.sort_files_by_number( ...
                {dir(fullfile(obj.dataset_path, 'pose', '*.txt')).name}); % cell
            pose_paths = fullfile(obj.dataset_path, 'pose', pose_files); % cell
            for i = 1:length(pose_paths)
                pose_id = fopen(pose_paths{i}, 'r');
                pose = textscan(pose_id, '%f');
                fclose(pose_id);
                c2w = reshape(pose{1}, [4, 4])';
                obj.poses(i, :, :) = c2w;
            end
        end

        function [rgb, depth, pose] = get_rgb_depth_pose(obj, idx)
            % Return rgb image, depth image, and camera pose with given index
            rgb = imread(obj.color_paths(idx));
            % rgb = cv2.resize(rgb, (obj.dataset_config["W"], obj.dataset_config["H"]))
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
        function sorted_files = sort_files_by_number(file_names)
            numbers = zeros(1, numel(file_names));
            for i = 1:numel(file_names)
                numbers(i) = str2double(regexp(file_names{i}, '\d+', 'match'));
            end
            [~, sorted_idx] = sort(numbers);
            sorted_files = file_names(sorted_idx);
        end
    end

end