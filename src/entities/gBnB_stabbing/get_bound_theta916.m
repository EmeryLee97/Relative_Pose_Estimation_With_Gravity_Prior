function [upper_bound, lower_bound, theta] = get_bound_theta916( ...
    a, b, c, branch_2d, epsilon)
    % Interval merging tool box needed (Version 1.2.0.0 by Bruno Luong)

    % Branch center & branch half length
    t_c = exponential_mapping_2dto3d(branch_2d(1:2));
    delta = branch_2d(3);

    % sample theta between [-sqrt(2)*delta, sqrt(2)*delta]
    % sample step = 0.02 (0.02rad = 1.14deg)
    gamma = -sqrt(2)*delta : 0.001 : sqrt(2)*delta; % (1, m)

    % norm of vectors a, b, c
    a_norm = vecnorm(a); % (1, n)
    b_norm = vecnorm(b); % (1, n)
    c_norm = vecnorm(c); % (1, n)

    % angles between tc and a, b, c
    angle_a_tc = acos(t_c' * a ./ a_norm); % (1, n)
    angle_b_tc = acos(t_c' * b ./ b_norm); % (1, n)
    angle_c_tc = acos(t_c' * c ./ c_norm); % (1, n)

    % angles between t and a, b, c               % a_gamma:  a(1) a(2) a(3) ... a(n)
    % angle_t_tc <= sqrt(2) * delta              % gamma(1): 
    % m = length(gamma), n = size(a, 2);         % gamma(2):
    angle_a_t = angle_a_tc + gamma'; % (m, n)    % gamma(3):     a+gamma
    angle_b_t = angle_b_tc + gamma'; % (m, n)    %  ...... 
    angle_c_t = angle_c_tc + gamma'; % (m, n)    % gamma(m):

    % A, B, C are dot products of a, b, c and t
    A = a_norm .* cos(angle_a_t); % (m, n)
    B = b_norm .* cos(angle_b_t); % (m, n)
    C = c_norm .* cos(angle_c_t); % (m, n)

    denominator = sqrt(B.^2 + C.^2); % (m, n)

    % phi ~ [-2*pi, 2*pi]
    sin_phi = B ./ denominator; % (m, n)
    cos_phi = C ./ denominator; % (m, n)
    phi = atan(B ./ C); % (m, n)

    % determine phi with sin(phi) and cos(phi)
    % sin(phi) > 0, cos(phi) > 0: phi ~ (0, pi/2), Quadrant I
    % sin(phi) > 0, cos(phi) < 0: phi ~ (pi/2, pi), Quadrant II
    % sin(phi) < 0, cos(phi) < 0: phi ~ (-pi, -pi/2), Quadrant III
    % sin(phi) < 0, cos(phi) > 0: phi ~ (-pi/2, 0), Quadrant IV
    
    %idx_sin_pos = sin_phi > 0;
    idx_sin_neg = sin_phi < 0;
    %idx_cos_pos = cos_phi > 0;
    idx_cos_neg = cos_phi < 0;

    idx_phi_pos = phi > 0;
    %idx_phi_neg = ~idx_phi_pos;

    idx_phi_1 = idx_phi_pos & idx_sin_neg & idx_cos_neg;
    phi(idx_phi_1) = phi(idx_phi_1) - pi;
    idx_phi_2 = ~(idx_phi_pos) & ~(idx_sin_neg) & idx_cos_neg;
    phi(idx_phi_2) = phi(idx_phi_2) + pi;
    
%     phi(idx_phi_pos & idx_sin_neg & idx_cos_neg) = phi(idx_phi_pos & idx_sin_neg & idx_cos_neg) - pi;
%     phi(idx_phi_neg & idx_sin_pos & idx_cos_neg) = phi(idx_phi_neg & idx_sin_pos & idx_cos_neg) + pi;

    
    A_tilde = A ./ denominator; % (m, n)
    epsilon_tilde = epsilon ./ denominator; % (m, n)

    % k1 <= cos(theta - phi) <= k2
    % note that k1 might > pi and k2 might < -pi
    k1 = max(-1, -epsilon_tilde - A_tilde); % (m, n)
    k2 = min(1, epsilon_tilde - A_tilde); % (m, n)

    k2(k2 < -1) = nan; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k1(k1 > 1) = nan;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    lower = [];
    upper = [];
    for i = 1:size(theta_min_1, 2) 
        l1 = theta_min_1(i); r1 = theta_max_1(i);
        if (theta_min_1(i) < -pi) 
            l1 = -pi;
        end
        if (theta_max_1(i) < -pi)
            l1 = nan;
            r1 = nan;
        end

        l2 = theta_min_2(i); r2 = theta_max_2(i);
        if (theta_max_2(i) > pi)
            r2 = pi;
        end
        if (theta_min_2(i) > pi)
            l2 = nan;
            r2 = nan;
        end

        l3 = theta_min_3(i); r3 = theta_max_3(i);
        if (theta_min_3(i) < -pi) 
            l3 = -pi;
        end
        if (theta_max_3(i) < -pi)
            l3 = nan;
            r3 = nan;
        end

        l4 = theta_min_4(i); r4 = theta_max_4(i);
        if (theta_max_4(i) > pi)
            r4 = pi;
        end
        if (theta_min_4(i) > pi)
            l4 = nan;
            r4 = nan;
        end

        [l, r] = IntervalUnion([l1, l2, l3, l4], [r1, r2, r3, r4]);
        lower = [lower, l];
        upper = [upper, r];
    end

    % get upper bound
    [upper_bound, theta] = solve_1d_stabbing(lower, upper, "middle");

    % get lower bound
    d = t_c' * (a + b .* sin(theta) + c .* cos(theta));
    lower_bound = sum(abs(d) <= epsilon);
end
