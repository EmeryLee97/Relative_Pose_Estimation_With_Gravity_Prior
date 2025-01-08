function [upper_bound, lower_bound, theta] = get_theta_bound_sampling(branch_2d, a, b, c, num, epsilon)
    
    % Sample theta from -pi to pi, num times
    theta_set = linspace(-pi, pi, num); % (1, num)

    t_c = exponential_mapping_2dto3d(branch_2d(1:2));
    delta = sqrt(2) * branch_2d(end);
    
    % range of t'*a, t'*b, t'*c
    A = t_c' * a; 
    norm_a = vecnorm(a);

    A_angle_min = min(acos(A ./ norm_a) + delta, pi);
    % A_angle_min = acos(A ./ norm_a) + delta;
    % A_angle_min(A_angle_min > pi) = pi;
    A_min = norm_a .* cos(A_angle_min);

    A_angle_max = max(acos(A ./ norm_a) - delta, 0);
    % A_angle_max = acos(A ./ norm_a) - delta;
    % A_angle_max(A_angle_max < 0) = 0;
    A_max = norm_a .* cos(A_angle_max);
    

    B = t_c' * b; 
    norm_b = vecnorm(b);

    B_angle_min = min(acos(B ./ norm_b) + delta, pi);
    % B_angle_min = acos(B ./ norm_b) + delta;
    % B_angle_min(B_angle_min > pi) = pi;
    B_min = norm_b .* cos(B_angle_min);

    B_angle_max = max(acos(B ./ norm_b) - delta, 0);
    % B_angle_max = acos(B ./ norm_b) - delta;
    % B_angle_max(B_angle_max < 0) = 0;
    B_max = norm_b .* cos(B_angle_max);


    C = t_c' * c; 
    norm_c = vecnorm(c);

    C_angle_min = min(acos(C ./ norm_c) + delta, pi);
    % C_angle_min = acos(C ./ norm_c) + delta;
    % C_angle_min(C_angle_min > pi) = pi;
    C_min = norm_c .* cos(C_angle_min);

    C_angle_max = max(acos(C ./ norm_c) - delta, 0);
    % C_angle_max = acos(C ./ norm_c) - delta;
    % C_angle_max(C_angle_max < 0) = 0;
    C_max = norm_c .* cos(C_angle_max);

    
    sin_theta_set = sin(theta_set); % (1, num)
    cos_theta_set = cos(theta_set); % (1, num)

    % All of these are (num_theta, num_pts) matrices
    sum_1_set = sin_theta_set' * B_max + cos_theta_set' * C_max; % (num, n)
    sum_2_set = sin_theta_set' * B_max + cos_theta_set' * C_min; % (num, n)
    sum_3_set = sin_theta_set' * B_min + cos_theta_set' * C_max; % (num, n)
    sum_4_set = sin_theta_set' * B_min + cos_theta_set' * C_min; % (num, n)

    % (num_theta, num_pts, 4) 3d array
    sum_3d_matrix(:, :, 1) = sum_1_set;
    sum_3d_matrix(:, :, 2) = sum_2_set;
    sum_3d_matrix(:, :, 3) = sum_3_set;
    sum_3d_matrix(:, :, 4) = sum_4_set;

    % (num_theta, num_pts) matries
    sum_max_set = A_max + max(sum_3d_matrix, [], 3);
    sum_min_set = A_min + min(sum_3d_matrix, [], 3);    

    % Calculate upper- and lower bound
    sum1 = sum(sum_min_set < 0 & sum_max_set > 0, 2);
    sum2 = sum(-sum_max_set <= epsilon & sum_max_set <= 0, 2);
    sum3 = sum(sum_min_set <= epsilon & sum_min_set >= 0, 2);

    upper_bound_set = sum1 + sum2 + sum3;
    
    [upper_bound, idx_theta] = max(upper_bound_set);
    theta = theta_set(idx_theta);
    lower_bound = sum(abs(A + B * sin(theta) + C * cos(theta)) <= epsilon);
end

