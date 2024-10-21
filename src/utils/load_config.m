% https://github.com/MaartenJongeneel/readyaml

function cfg = load_config(yaml_path)
    cfg_scene = readyaml(yaml_path); 
    if isfield(cfg_scene, 'inherit_from')
        inherit_path = cfg_scene.inherit_from;
        cfg_dataset = readyaml(inherit_path);
    elseif isfield(cfg_scene, 'default_path')
        default_path = cfg_scene.default_path;
        cfg_dataset = readyaml(default_path);
    else
        cfg_dataset = struct();
    end
    cfg = update_config_recursive(cfg_dataset, cfg_scene);
 
end


function updated_struct = update_config_recursive(struct_1, struct_2)

    updated_struct = struct_1;
    fields_2 = fieldnames(struct_2);
    for i = 1:numel(fields_2)
        field = fields_2{i};
        value = struct_2.(field);
        
        if isfield(updated_struct, field)
            if isstruct(value)
                updated_struct.(field) = update_config_recursive(updated_struct.(field), value);
            else
                updated_struct.(field) = value;
            end
        else
            updated_struct.(field) = value;
        end
    end

end