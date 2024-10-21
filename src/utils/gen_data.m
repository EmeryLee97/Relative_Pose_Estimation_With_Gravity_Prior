function [theta, t, x, y, R_v] = gen_data(num_inlier, num_outlier, noise_level)
% vecnorm(x) = 1, vecnomr(y) = 1;

    theta = (rand() * 2 - 1) * pi/2;
    y = tan(theta / 2);
    y_square = y^2;

    % Normally, R_v can be any 3*1 vector, whoes norm is 1,
    % This special case is for Ding
    R_v = [0, 1, 0];
    R = [1 - y_square, 0, 2 * y; ...
         0, 1 + y_square, 0; ...
         -2 * y, 0, 1 - y_square] / (1 + y_square);
    
    num = num_inlier + num_outlier;
    
    % Generate points in the 1st image
    X_3d = rand(3, num) * 2 - 1 + [0; 0; 2];
    x = X_3d ./ vecnorm(X_3d);
    
    X_3d_c = mean(X_3d, 2);
    T = -R * X_3d_c + X_3d_c + rand(3, 1);
    t = T ./ vecnorm(T);
    
    % Generate 2 kinds of noises in the image plane
    noise_1 = normrnd(0, noise_level, [2, num_inlier]);
    noise_2 = normrnd(0, 40, [2, num_outlier]);
    
    focal_length = 1000;
    Y_3d_0 = R * X_3d + T;
    
    % Add these 2 noises to the 2nd image
    Y_3d_1 = focal_length .* Y_3d_0(1:2, 1:num_inlier) ./ Y_3d_0(3, 1:num_inlier) + noise_1;
    Y_3d_2 = focal_length .* Y_3d_0(1:2, end-num_outlier+1:end) ./ Y_3d_0(3, end-num_outlier+1:end) + noise_2;
    
    Y_3d = [[Y_3d_1, Y_3d_2]; focal_length*ones(1, num)];
    y = Y_3d ./ vecnorm(Y_3d);

end