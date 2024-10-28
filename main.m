clear; clc; close all;

addpath('src/evaluation/')


cfg_path = '/Users/xuhuili/Documents/SA/Relative_Pose_Estimation_With_Gravity_Prior/config/ScanNet/scene0059_00.yaml';

evaluator = RealworldEvaluator(cfg_path);