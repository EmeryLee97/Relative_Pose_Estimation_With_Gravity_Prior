function [upper_bound, lower_bound, theta] = get_theta_bound_stabbing( ...
    a, b, c, branch_2d, epsilon)
    
    % Branch center & branch half length
    t_c = exponential_mapping_2dto3d(branch_2d(1:2));
    delta = sqrt(2) * branch_2d(3); % to speed up

    % sample theta between [delta, delta]
    % sample step = 0.02 (0.02rad = 1.14deg)
    gamma = -delta : 0.002 : delta; % (1, m)

%     if gamma(end) ~= delta
%         gamma = [gamma, delta];
%     end

    % norm of vectors a, b, c
    a_norm = vecnorm(a); % (1, n)
    b_norm = vecnorm(b); % (1, n)
    c_norm = vecnorm(c); % (1, n)

    % angles between tc and a, b, c
    angle_a_tc = acos(t_c' * a ./ a_norm); % (1, n), [0, pi]
    angle_b_tc = acos(t_c' * b ./ b_norm); % (1, n), [0, pi]
    angle_c_tc = acos(t_c' * c ./ c_norm); % (1, n), [0, pi]

    % angles between t and a, b, c               % a_gamma:  a(1) a(2) a(3) ... a(n)
    % angle_t_tc <= sqrt(2) * delta              % gamma(1): 
    % m = length(gamma), n = size(a, 2);         % gamma(2):
    angle_a_t = angle_a_tc + gamma'; % (m, n)    % gamma(3):     a(i)+gamma(j)
    angle_b_t = angle_b_tc + gamma'; % (m, n)    %  ...... 
    angle_c_t = angle_c_tc + gamma'; % (m, n)    % gamma(m):

    %-----------------------------------
%     angle_a_t(angle_a_t > pi) = pi;
%     angle_a_t(angle_a_t < 0) = 0;
%     angle_b_t(angle_b_t > pi) = pi;
%     angle_b_t(angle_b_t < 0) = 0;
%     angle_c_t(angle_c_t > pi) = pi;
%     angle_c_t(angle_c_t < 0) = 0;
    %-----------------------------------

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

%     lower = [];
%     upper = [];
%     for i = 1:size(theta_min_1, 2)
%         [left, right] = MergeBrackets( ...
%             [theta_min_1(i), theta_min_2(i), theta_min_3(i), theta_min_4(i)], ...
%             [theta_max_1(i), theta_max_2(i), theta_max_3(i), theta_max_4(i)]);
%         lower = [lower, left];
%         upper = [upper, right];
%     end

    % speed up with fixed size
    num = 4 * size(theta_min_1, 2);
    lower = zeros(1, num);
    upper = zeros(1, num);
    k = 0;
    for i = 1:size(theta_min_1, 2)
        [left, right] = interval_union( ...
            [theta_min_1(i), theta_min_2(i), theta_min_3(i), theta_min_4(i)], ...
            [theta_max_1(i), theta_max_2(i), theta_max_3(i), theta_max_4(i)]);
        ii = size(left, 2);
        lower(k+1:k+ii) = left;
        upper(k+1:k+ii) = right;
        k = k + ii;
    end
    lower(k+1:end) = [];
    upper(k+1:end) = [];

    % get upper bound
    [upper_bound, theta] = solve_1d_stabbing(lower, upper);

    % get lower bound
    d = t_c' * (a + b .* sin(theta) + c .* cos(theta));
    lower_bound = sum(abs(d) <= epsilon);
    
end

