classdef StereoDataset < BaseDataset
    properties (Access=public)
        
        fx_cam0 (1, 1) {mustBeA(fx_cam0, {'int', 'single', 'double'})}
        fy_cam0 (1, 1) {mustBeA(fy_cam0, {'int', 'single', 'double'})}
        cx_cam0 (1, 1) {mustBeA(cx_cam0, {'int', 'single', 'double'})}
        cy_cam0 (1, 1) {mustBeA(cy_cam0, {'int', 'single', 'double'})} 
        fovx_cam0 (1, 1) {mustBeA(fovx_cam0, {'int', 'single', 'double'})}
        fovy_cam0 (1, 1){mustBeA(fovy_cam0, {'int', 'single', 'double'})}
        intrinsics_cam0 (3, 3) {mustBeA(intrinsics_cam0, {'int', 'single', 'double'})}
        extrinsics_cam0 (4, 4) {mustBeA(extrinsics_cam0, {'int', 'single', 'double'})}
        distortion_cam0 
        params_cam0
        img_paths_cam0 string
        
        fx_cam1 (1, 1) {mustBeA(fx_cam1, {'int', 'single', 'double'})}
        fy_cam1 (1, 1) {mustBeA(fy_cam1, {'int', 'single', 'double'})}
        cx_cam1 (1, 1) {mustBeA(cx_cam1, {'int', 'single', 'double'})}
        cy_cam1 (1, 1) {mustBeA(cy_cam1, {'int', 'single', 'double'})} 
        fovx_cam1 (1, 1) {mustBeA(fovx_cam1, {'int', 'single', 'double'})}
        fovy_cam1 (1, 1){mustBeA(fovy_cam1, {'int', 'single', 'double'})}
        intrinsics_cam1 (3, 3) {mustBeA(intrinsics_cam1, {'int', 'single', 'double'})}
        extrinsics_cam1 (4, 4) {mustBeA(extrinsics_cam1, {'int', 'single', 'double'})}
        distortion_cam1 
        params_cam1
        img_paths_cam1 string
    end
    
    methods (Access=public)
        function obj = StereoDataset(scene_config_path)
            addpath('src/utils/')
            obj@BaseDataset(scene_config_path)

            % camera_0 parameters
            obj.fx_cam0 = obj.dataset_config.cam0.fx;
            obj.fy_cam0 = obj.dataset_config.cam0.fy;
            obj.cx_cam0 = obj.dataset_config.cam0.cx;
            obj.cy_cam0 = obj.dataset_config.cam0.cy;
            
            if obj.crop_edge ~= 0
                obj.cx_cam0 = obj.cx_cam0 - obj.crop_edge;
                obj.cy_cam0 = obj.cy_cam0 - obj.crop_edge;
            end

            obj.fovx_cam0 = 2 * atan(obj.width / (2 * obj.fx_cam0));
            obj.fovy_cam0 = 2 * atan(obj.height / (2 * obj.fy_cam0));
            obj.intrinsics_cam0 = [
                obj.fx_cam0, 0,           obj.cx_cam0; ...
                0,           obj.fy_cam0, obj.cy_cam0; ...
                0,           0,           1
                ];
            xyz_quat_cam0 = obj.dataset_config.cam0.T_w2c;
            obj.extrinsics_cam0 = eye(4);
            obj.extrinsics_cam0(1:3, 1:3) = quat2rotm(xyz_quat_cam0(4:end));
            obj.extrinsics_cam0(1:3, 4) = xyz_quat_cam0(1:3);

            if isfield(obj.dataset_config.cam0, 'distortion')
                obj.distortion_cam0 = obj.dataset_config.cam.distortion;
                obj.params_cam0 = cameraParameters(...
                    'IntrinsicMatrix', obj.intrinsics_cam0', ...
                    'RadialDistortion', [obj.distortion_cam0(1:2), obj.distortion_cam0(5)], ...
                    'TangentialDistortion', obj.distortion_cam0(3:4));
            end

            % camera_1 parameters
            obj.fx_cam1 = obj.dataset_config.cam1.fx;
            obj.fy_cam1 = obj.dataset_config.cam1.fy;
            obj.cx_cam1 = obj.dataset_config.cam1.cx;
            obj.cy_cam1 = obj.dataset_config.cam1.cy;
            
            if obj.crop_edge ~= 0
                obj.cx_cam1 = obj.cx_cam1 - obj.crop_edge;
                obj.cy_cam1 = obj.cy_cam1 - obj.crop_edge;
            end

            obj.fovx_cam1 = 2 * atan(obj.width / (2 * obj.fx_cam1));
            obj.fovy_cam1 = 2 * atan(obj.height / (2 * obj.fy_cam1));
            obj.intrinsics_cam1 = [
                obj.fx_cam1, 0,           obj.cx_cam1; ...
                0,           obj.fy_cam1, obj.cy_cam1; ...
                0,           0,           1
                ];
            xyz_quat_cam1 = obj.dataset_config.cam1.T_w2c;
            obj.extrinsics_cam1 = eye(4);
            obj.extrinsics_cam1(1:3, 1:3) = quat2rotm(xyz_quat_cam1(4:end));
            obj.extrinsics_cam1(1:3, 4) = xyz_quat_cam1(1:3);

            if isfield(obj.dataset_config.cam1, 'distortion')
                obj.distortion_cam1 = obj.dataset_config.cam.distortion;
                obj.params_cam1 = cameraParameters(...
                    'IntrinsicMatrix', obj.intrinsics_cam1', ...
                    'RadialDistortion', [obj.distortion_cam1(1:2), obj.distortion_cam1(5)], ...
                    'TangentialDistortion', obj.distortion_cam1(3:4));
            end
        end
    
    end
end


