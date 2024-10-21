function [upper_bound, lower_bound, theta] = get_bound_theta929( ...
    a, b, c, branch_2d, epsilon)
    % Interval merging tool box needed (Version 1.2.0.0 by Bruno Luong)
    %dbstop in get_bound_theta917 at 148 if abs(branch_2d(3)-0.0491)<1e-3

    % Branch center & branch half length
    t_c = exponential_mapping_2dto3d(branch_2d(1:2));
    delta = branch_2d(3);

    % sample theta between [-sqrt(2)*delta, sqrt(2)*delta]
    % sample step = 0.02 (0.02rad = 1.14deg)
    gamma = -sqrt(2)*delta : 0.02 : sqrt(2)*delta; % (1, m)
    if gamma(end) ~= sqrt(2)*delta
        gamma = [gamma, sqrt(2)*delta];
    end
    gamma_set{1} = [gamma; gamma; gamma];
    gamma_set{2} = [gamma; gamma; -gamma];
    gamma_set{3} = [gamma; -gamma; gamma];
    gamma_set{4} = [gamma; -gamma; -gamma];

    % norm of vectors a, b, c
    a_norm = vecnorm(a); % (1, n)
    b_norm = vecnorm(b); % (1, n)
    c_norm = vecnorm(c); % (1, n)

    % angles between tc and a, b, c
    angle_a_tc = acos(t_c' * a ./ a_norm); % (1, n), [0, pi]
    angle_b_tc = acos(t_c' * b ./ b_norm); % (1, n), [0, pi]
    angle_c_tc = acos(t_c' * c ./ c_norm); % (1, n), [0, pi]

    max_upper_bound = 0;
    theta_opt = 0;

    for k = 1:4
    % angles between t and a, b, c              
    % angle_t_tc <= sqrt(2) * delta      
    % m = length(gamma), n = size(a, 2);      
    angle_a_t = angle_a_tc + gamma_set{k}(1, :)'; % (m, n) 
    angle_b_t = angle_b_tc + gamma_set{k}(2, :)'; % (m, n) 
    angle_c_t = angle_c_tc + gamma_set{k}(3, :)'; % (m, n) 

    % A, B, C are dot products of a, b, c and t
    A = a_norm .* cos(angle_a_t); % (m, n)
    B = b_norm .* cos(angle_b_t); % (m, n)
    C = c_norm .* cos(angle_c_t); % (m, n)

    denominator = sqrt(B.^2 + C.^2); % (m, n)

    % phi ~ [-2*pi, 2*pi]
    sin_phi = B ./ denominator; % (m, n)
    cos_phi = C ./ denominator; % (m, n)
    phi = atan(B ./ C); % (m, n), [-pi/2, pi/2]

    % determine phi with sin(phi) and cos(phi)
    % sin(phi) > 0, cos(phi) > 0: phi ~ (0, pi/2), Quadrant I
    % sin(phi) > 0, cos(phi) < 0: phi ~ (pi/2, pi), Quadrant II
    % sin(phi) < 0, cos(phi) < 0: phi ~ (-pi, -pi/2), Quadrant III
    % sin(phi) < 0, cos(phi) > 0: phi ~ (-pi/2, 0), Quadrant IV
    
    %idx_sin_pos = sin_phi > 0;
    idx_sin_neg = sin_phi < 0;
    %idx_cos_pos = cos_phi > 0;
    idx_cos_neg = cos_phi < 0;

    %idx_phi_pos = phi > 0;
    %idx_phi_neg = ~idx_phi_pos;

    idx_phi_1 = idx_sin_neg & idx_cos_neg;
    phi(idx_phi_1) = phi(idx_phi_1) - pi;
    idx_phi_2 = ~(idx_sin_neg) & idx_cos_neg;
    phi(idx_phi_2) = phi(idx_phi_2) + pi;
    
%     phi(idx_phi_pos & idx_sin_neg & idx_cos_neg) = phi(idx_phi_pos & idx_sin_neg & idx_cos_neg) - pi;
%     phi(idx_phi_neg & idx_sin_pos & idx_cos_neg) = phi(idx_phi_neg & idx_sin_pos & idx_cos_neg) + pi;

    
    A_tilde = A ./ denominator; % (m, n)
    epsilon_tilde = epsilon ./ denominator; % (m, n)

    % k1 <= cos(theta - phi) <= k2
    % note that k1 might > pi and k2 might < -pi
    k1 = max(-1, -epsilon_tilde - A_tilde); % (m, n)
    k2 = min(1, epsilon_tilde - A_tilde); % (m, n)

%     idx_k1 = k1 > 1;
%     idx_k2 = k2 < -1;
% 
%     k1(idx_k1) = nan;
%     k1(idx_k2) = nan;
% 
%     k2(idx_k1) = nan;
%     k2(idx_k2) = nan;

    k2(k2 < -1) = -1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k1(k1 > 1) = 1;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    acos_k1 = acos(k1); % (m, n)
    acos_k2 = acos(k2); % (m, n)

    % case 1: |theta - phi| <= pi
    % theta ~ [phi - arccos(k1), phi - arccos(k2)] ||
    %         [phi + arccos(k2), phi + arccos(k1)] && [-pi, pi]

    theta_min_1 = min(phi - acos_k1, [], 1); % (1, n)  [-2*pi, pi]
    theta_max_1 = max(phi - acos_k2, [], 1); % (1, n)  [-2*pi, pi]

    theta_min_2 = min(phi + acos_k2, [], 1); % (1, n)  [-pi, 2*pi]
    theta_max_2 = max(phi + acos_k1, [], 1); % (1, n)  [-pi, 2*pi]

    % case 2: pi <= |theta - phi| <= 2*pi
    % theta ~ [phi + arccos(k2) - 2*pi, phi + arccos(k1) - 2*pi] ||
    %         [phi - arccos(k1) + 2*pi, phi - arccos(k2) + 2*pi] && [-pi, pi]
    theta_min_3 = theta_min_2 - 2*pi;
    theta_max_3 = theta_max_2 - 2*pi;

    theta_min_4 = theta_min_1 + 2*pi;
    theta_max_4 = theta_max_1 + 2*pi;

    % merge intervals
    idx_l1_neg = theta_min_1 < -pi;
    idx_r1_neg = theta_max_1 < -pi;
    theta_min_1(idx_l1_neg) = -pi;
    theta_min_1(idx_r1_neg) = nan;
    theta_max_1(idx_r1_neg) = nan;

    idx_r2_pos = theta_max_2 > pi;
    idx_l2_pos = theta_min_2 > pi;
    theta_max_2(idx_r2_pos) = pi;
    theta_min_2(idx_l2_pos) = nan;
    theta_max_2(idx_l2_pos) = nan;

    idx_l3_neg = theta_min_3 < -pi;
    idx_r3_neg = theta_max_3 < -pi;
    theta_min_3(idx_l3_neg) = -pi;
    theta_min_3(idx_r3_neg) = nan;
    theta_max_3(idx_r3_neg) = nan;

    idx_l4_pos = theta_max_4 > pi;
    idx_r4_pos = theta_min_4 > pi;
    theta_max_4(idx_l4_pos) = pi;
    theta_max_4(idx_r4_pos) = nan;
    theta_min_4(idx_r4_pos) = nan;

    lower = [];
    upper = [];
    for i = 1:size(theta_min_1, 2)
        [left, right] = interval_union( ...
            [theta_min_1(i), theta_max_1(i); theta_min_2(i), theta_max_2(i); ...
             theta_min_3(i), theta_max_3(i); theta_min_4(i), theta_max_4(i)]);
        lower = [lower, left];
        upper = [upper, right];
    end

    % get upper bound
    [upper_bound_temp, theta_temp] = solve_1d_stabbing(lower, upper);
    if upper_bound_temp > max_upper_bound
        max_upper_bound = upper_bound_temp;
        theta_opt = theta_temp;
    end
    end

    upper_bound = max_upper_bound;
    theta = theta_opt;
    % get lower bound
    d = t_c' * (a + b .* sin(theta) + c .* cos(theta));
    lower_bound = sum(abs(d) <= epsilon);
    
end

