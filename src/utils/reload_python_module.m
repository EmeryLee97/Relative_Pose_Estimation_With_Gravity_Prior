function module = reload_python_module(module_name)
    try
        module = py.importlib.import_module(module_name);
        py.importlib.reload(module);
    catch
        disp(['Module ', module_name, ' not found.']);
    end
end