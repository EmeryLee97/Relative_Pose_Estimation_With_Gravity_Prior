classdef BaseDataset
    properties (Access=public)
        dataset_path (1, 1) string
        dataset_config struct
        height (1, 1) {mustBeInteger}
        width (1, 1) {mustBeInteger}
        crop_edge (1, 1) {mustBeInteger}
    end
    
    methods (Access=public)
        function obj = BaseDataset(scene_config_path)
            addpath('src/utils/')

            obj.dataset_config = load_config(scene_config_path);
            obj.dataset_path = obj.dataset_config.data.input_path;
            obj.height = obj.dataset_config.cam.H;
            obj.width = obj.dataset_config.cam.W;
            
            if isfield(obj.dataset_config.cam, 'crop_edge')
                obj.crop_edge = obj.dataset_config.cam.crop_edge;
            else
                obj.crop_edge = 0;
            end
            
            if obj.crop_edge ~= 0
                obj.height = obj.height - 2 * obj.crop_edge;
                obj.width = obj.width - 2 * obj.crop_edge;
            end
        end
    
    end
end


