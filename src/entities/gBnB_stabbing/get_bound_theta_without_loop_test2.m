function [upper_bound, lower_bound, theta] = get_bound_theta_without_loop_test2( ...
    N, a, b, c, branch_2d, epsilon)

    % Branch center & branch half length
    t_c = exponential_mapping_2dto3d(branch_2d(1:2));
    delta = branch_2d(3);

    % A = t_c' * a, (1, N) vector
    A_min = vecnorm(a) .* cos(acos(t_c' * a ./ vecnorm(a)) + sqrt(2) * delta);
    A_max = vecnorm(a) .* cos(acos(t_c' * a ./ vecnorm(a)) - sqrt(2) * delta);

    % B = t_c' * b, (1, N) vector
    B_min = vecnorm(b) .* cos(acos(t_c' * b ./ vecnorm(b)) + sqrt(2) * delta);
    B_max = vecnorm(b) .* cos(acos(t_c' * b ./ vecnorm(b)) - sqrt(2) * delta);

    % C = t_c' * c, (1, N) vector
    C_min = vecnorm(c) .* cos(acos(t_c' * c ./ vecnorm(c)) + sqrt(2) * delta);
    C_max = vecnorm(c) .* cos(acos(t_c' * c ./ vecnorm(c)) - sqrt(2) * delta);

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
    quotient_minus_max = max([A_ep_minus_min ./ B_C_min; ...
        A_ep_minus_min ./ B_C_max; ...
        A_ep_minus_max ./ B_C_min; ...
        A_ep_minus_max ./ B_C_max], [], 1);

    % (epsilon - A) / sqrt(B^2 + C^2), (1, N) vector
    quotient_plus_min = min([A_ep_plus_min ./ B_C_min; ...
        A_ep_plus_min ./ B_C_max; ...
        A_ep_plus_max ./ B_C_min; ...
        A_ep_plus_max ./ B_C_max], [], 1);
    quotient_plus_max = max([A_ep_plus_min ./ B_C_min; ...
        A_ep_plus_min ./ B_C_max; ...
        A_ep_plus_max ./ B_C_min; ...
        A_ep_plus_max ./ B_C_max], [], 1);

    % C / B, (1, N-k) vector
    atan_min = atan(min([C_min./B_min; C_min./B_max; C_max./B_min; C_max./B_max], [], 1));
    atan_max = atan(max([C_min./B_min; C_min./B_max; C_max./B_min; C_max./B_max], [], 1));

    ind_1 = quotient_minus_min >= -1 & quotient_plus_max <= 1;
    ind_3 = quotient_minus_min < -1 & quotient_plus_max <= 1 & quotient_plus_max >= -1;
    ind_2 = quotient_minus_min >= -1 & quotient_minus_min <= 1 & quotient_plus_max > 1;
    ind_4 = quotient_minus_min < -1 & quotient_plus_max > 1;

    theta_min_1 = zeros(1, N);
    theta_max_1 = zeros(1, N);

    theta_min_1(ind_1) = asin(quotient_minus_min(ind_1)) - atan_max(ind_1);
    theta_max_1(ind_1) = asin(quotient_plus_max(ind_1)) - atan_min(ind_1);

    theta_min_1(ind_2) = asin(quotient_minus_min(ind_2)) - 2*pi - atan_max(ind_2);
    theta_max_1(ind_2) = -pi - asin(quotient_minus_min(ind_2)) - atan_min(ind_2);

    % 这两个存疑-----------------------------------
    theta_min_1(ind_3) = -pi - asin(quotient_plus_max(ind_3)) - atan_max(ind_3);
    theta_max_1(ind_3) = asin(quotient_plus_max(ind_3)) - atan_min(ind_3);
    % --------------------------------------------
    
    theta_min_1(ind_4) = -pi;
    theta_max_1(ind_4) = pi;


    theta_min_2 = zeros(1, N);
    theta_max_2 = zeros(1, N);

    theta_min_2(ind_1) = pi - asin(quotient_plus_max(ind_1)) - atan_max(ind_1);
    theta_max_2(ind_1) = pi - asin(quotient_minus_min(ind_1)) - atan_min(ind_1);

    theta_min_2(ind_2) = asin(quotient_minus_min(ind_2)) - atan_max(ind_2);
    theta_max_2(ind_2) = pi - asin(quotient_minus_min(ind_2)) - atan_min(ind_2);

    theta_min_2(ind_3) = pi - asin(quotient_plus_max(ind_3)) - atan_max(ind_3);
    theta_max_2(ind_3) = asin(quotient_plus_max(ind_3)) + 2*pi - atan_min(ind_3);
    
    theta_min_2(ind_4) = -pi;
    theta_max_2(ind_4) = pi;


    % Find the num and position that max theta overlap as upper bound
    [upper_bound_1, theta_1] = Max_Stabbing(theta_min_1, theta_max_1);
    lower_bound_1 = sum(abs(t_c' * (a + sin(theta_1)*b + cos(theta_1)*c)) <= epsilon);

    [upper_bound_2, theta_2] = Max_Stabbing(theta_min_2, theta_max_2);
    lower_bound_2 = sum(abs(t_c' * (a + sin(theta_2)*b + cos(theta_2)*c)) <= epsilon);

    if upper_bound_1 > upper_bound_2
        upper_bound = upper_bound_1;
        lower_bound = lower_bound_1;
        theta = theta_1;
    elseif upper_bound_1 < upper_bound_2
        upper_bound = upper_bound_2;
        lower_bound = lower_bound_2;
        theta = theta_2;
    else
        if lower_bound_1 <= lower_bound_2
            upper_bound = upper_bound_1;
            lower_bound = lower_bound_1;
            theta = theta_1;
        else
            upper_bound = upper_bound_2;
            lower_bound = lower_bound_2;
            theta = theta_2;
        end
    end
end
