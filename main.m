clear; clc; close all;
format short;

pe = pyenv;
if pe.Status == "Loaded" && pe.Version ~= "3.10"
    disp('To change the Python version, restart MATLAB, then call pyenv(Version="3.10").')
else
    pyenv(Version="3.10");
end


addpath('src/evaluation/')
addpath('src/entities/datasets/')

dataset = 'scannet';

if strcmpi(dataset, 'ScanNet')
    cfg_path = '/Users/xuhuili/Documents/SA/Relative_Pose_Estimation_With_Gravity_Prior/config/ScanNet/scene0059_00.yaml';
elseif strcmpi(dataset, 'TUM_RGBD')
    cfg_path = '/Users/xuhuili/Documents/SA/Relative_Pose_Estimation_With_Gravity_Prior/config/TUM_RGBD/rgbd_dataset_freiburg3_long_office_household.yaml';
elseif strcmpi(dataset, 'KITTI')
    cfg_path = '/Users/xuhuili/Documents/SA/Relative_Pose_Estimation_With_Gravity_Prior/config/KITTI/sequence_00.yaml';
end


evaluator = RealworldEvaluator(cfg_path, dataset);
% [rgb, ~, ~] = evaluator.dataset.get_rgb_depth_pose(50);
% imshow(rgb)
evaluator = evaluator.evaluate();

evaluator.print_result();

