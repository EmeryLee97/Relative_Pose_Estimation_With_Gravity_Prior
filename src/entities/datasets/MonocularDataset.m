classdef MonocularDataset < BaseDataset
    properties (Access=public)

        fx (1, 1) {mustBeA(fx, {'int', 'single', 'double'})}
        fy (1, 1) {mustBeA(fy, {'int', 'single', 'double'})}
        cx (1, 1) {mustBeA(cx, {'int', 'single', 'double'})}
        cy (1, 1) {mustBeA(cy, {'int', 'single', 'double'})} 
        fovx (1, 1) {mustBeA(fovx, {'int', 'single', 'double'})}
        fovy (1, 1){mustBeA(fovy, {'int', 'single', 'double'})}
        intrinsics (3, 3) {mustBeA(intrinsics, {'int', 'single', 'double'})}      
        depth_scale (1, 1) {mustBeA(depth_scale, {'int', 'single', 'double'})}
        distortion 
        camera_params

        color_paths string
        depth_paths string
        
    end
    
    methods (Access=public)
        function obj = MonocularDataset(scene_config_path)
            addpath('src/utils/')

            obj@BaseDataset(scene_config_path)
            obj.fx = obj.dataset_config.cam.fx;
            obj.fy = obj.dataset_config.cam.fy;
            obj.cx = obj.dataset_config.cam.cx;
            obj.cy = obj.dataset_config.cam.cy;

            if obj.crop_edge ~= 0
                obj.cx = obj.cx - obj.crop_edge;
                obj.cy = obj.cy - obj.crop_edge;
            end

            obj.fovx = 2 * atan(obj.width / (2 * obj.fx));
            obj.fovy = 2 * atan(obj.height / (2 * obj.fy));
            obj.intrinsics = [obj.fx, 0,       obj.cx; ...
                               0,       obj.fy, obj.cy; ...
                               0,       0,       1];
            
            if isfield(obj.dataset_config.cam, 'depth_scale')
                obj.depth_scale = obj.dataset_config.cam.depth_scale;
            end

            if isfield(obj.dataset_config.cam, 'distortion')
                obj.distortion = obj.dataset_config.cam.distortion;
                obj.camera_params = cameraParameters(...
                    'IntrinsicMatrix', obj.intrinsics', ...
                    'RadialDistortion', [obj.distortion(1:2), obj.distortion(5)], ...
                    'TangentialDistortion', obj.distortion(3:4));
            end

        end
    
    end
end


