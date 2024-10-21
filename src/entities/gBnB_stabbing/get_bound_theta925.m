function [upper_bound, lower_bound, theta] = get_bound_theta925( ...
    a, b, c, branch_2d, epsilon)

    % Branch center & branch half length
    t_c = exponential_mapping_2dto3d(branch_2d(1:2));
    delta = branch_2d(3);

    % A = t_c' * a, (1, N) vector
    angle_a_tc = acos(t_c' * a ./ vecnorm(a));
    angle_a_t_min = max(0, angle_a_tc - sqrt(2)*delta);
    angle_a_t_max = min(pi, angle_a_tc + sqrt(2)*delta);
    A_min = vecnorm(a) .* cos(angle_a_t_max);
    A_max = vecnorm(a) .* cos(angle_a_t_min);

    % B = t_c' * b, (1, N) vector
    angle_b_tc = acos(t_c' * b ./ vecnorm(b));
    angle_b_t_min = max(0, angle_b_tc - sqrt(2)*delta);
    angle_b_t_max = min(pi, angle_b_tc + sqrt(2)*delta);
    B_min = vecnorm(b) .* cos(angle_b_t_max);
    B_max = vecnorm(b) .* cos(angle_b_t_min);

    % C = t_c' * c, (1, N) vector
    angle_c_tc = acos(t_c' * c ./ vecnorm(c));
    angle_c_t_min = max(0, angle_c_tc - sqrt(2)*delta);
    angle_c_t_max = min(pi, angle_c_tc + sqrt(2)*delta);
    C_min = vecnorm(c) .* cos(angle_c_t_max);
    C_max = vecnorm(c) .* cos(angle_c_t_min);

    % B^2, (1, N) vector
    B_square_min = min(B_min.^2, B_max.^2);
    B_square_min(B_min .* B_max < 0) = 0;
    B_square_max = max(B_min.^2, B_max.^2);

    % C^2, (1, N) vector
    C_square_min = min(C_min.^2, C_max.^2);
    C_square_min(C_min .* C_max < 0) = 0;
    C_square_max = max(C_min.^2, C_max.^2);
    
    % sqrt((B^2 + C^2)), (1, N) vector
    B_C_min = sqrt(B_square_min + C_square_min);
    B_C_max = sqrt(B_square_max + C_square_max);

    % - epsilon - A, (1, N) vector
    A_ep_minus_min = -epsilon - A_max;
    A_ep_minus_max = -epsilon - A_min;

    % epsilon - A, (1, N) vector
    A_ep_plus_min = epsilon - A_max;
    A_ep_plus_max = epsilon - A_min;

    % (-epsilon - A) / sqrt(B^2 + C^2), (1, N) vector
    quotient_minus_min = min([A_ep_minus_min ./ B_C_min; ...
        A_ep_minus_min ./ B_C_max; ...
        A_ep_minus_max ./ B_C_min; ...
        A_ep_minus_max ./ B_C_max], [], 1);

    % (epsilon - A) / sqrt(B^2 + C^2), (1, N) vector
    quotient_plus_max = max([A_ep_plus_min ./ B_C_min; ...
        A_ep_plus_min ./ B_C_max; ...
        A_ep_plus_max ./ B_C_min; ...
        A_ep_plus_max ./ B_C_max], [], 1);

    % phi = arctan(C / B), (1, N-k) vector
    atan_min = atan(min([C_min./B_min; C_min./B_max; C_max./B_min; C_max./B_max], [], 1));
    atan_max = atan(max([C_min./B_min; C_min./B_max; C_max./B_min; C_max./B_max], [], 1));

    




    % Find the num and position that max theta overlap as upper bound
    [upper_bound, theta] = solve_1d_stabbing(theta_min, theta_max, "middle");

    % Find the lower bound of theta
    lower_bound = sum(abs(t_c' * (a + sin(theta)*b + cos(theta)*c)) <= epsilon);

end