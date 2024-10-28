% This class logs metrics of different input methods

classdef Logger
    properties (Access=public)
        method_list
        runtime_map
        rot_err_map
        trans_err_map
        dynamic_vars
    end

    methods (Access=public)
        function obj = Logger(method_list, row_list, col_list)
            % method_list: array of strings
            % row_list, col_list: 2d int array
            obj.method_list = method_list;
            obj.runtime_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.rot_err_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.trans_err_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.dynamic_vars = struct();

            for i = 1:length(method_list)
                
                runtime_method = sprintf('runtime_%s', method_list(i));
                rot_err_method = sprintf('rot_err_%s', method_list(i));
                trans_err_method = sprintf('trans_err_%s', method_list(i));

                obj.dynamic_vars.(runtime_method) = zeros(row_list(i), col_list(i));
                obj.dynamic_vars.(rot_err_method) = zeros(row_list(i), col_list(i));
                obj.dynamic_vars.(trans_err_method) = zeros(row_list(i), col_list(i));
            end
        end

        function obj = update_dict(obj)
            for i = 1:length(obj.method_list)
                runtime_method = sprintf('runtime_%s', obj.method_list(i));
                obj.runtime_map(obj.method_list(i)) = obj.dynamic_vars.(runtime_method);

                rot_err_method = sprintf('rot_err_%s', obj.method_list(i));
                obj.rot_err_map(obj.method_list(i)) = obj.dynamic_vars.(rot_err_method);

                trans_err_method = sprintf('trans_err_%s', obj.method_list(i));
                obj.trans_err_map(obj.method_list(i)) = obj.dynamic_vars.(trans_err_method);
            end
        end

    end


end