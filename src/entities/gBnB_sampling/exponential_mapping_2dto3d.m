function point_3d = exponential_mapping_2dto3d(point_2d)
% This function converts a 2d point to a 3d point
% point_2d = alpha * [cos(beta); sin(beta)]
% point_3d = [sin(alpha)cos(beta); sin(alpha)sin(beta); cos(alpha)]


    % alpha -> [0, pi/2]
    % beta -> [-pi, pi]
    alpha = vecnorm(point_2d);
    if alpha == 0
        point_3d = [0; 0; 1];
    else
        vec = point_2d / alpha;
        point_3d = [sin(alpha) * vec; cos(alpha)];
    end
    
end